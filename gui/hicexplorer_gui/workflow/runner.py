"""Workflow execution: thread budget, resume, per-step records, cancellation.

State lives in ``<workdir>/.hicexplorer-workflow/``:

* ``state.json``: per step the last key, status and the content hashes of its
  outputs. A step is skipped when its key is unchanged, its last run succeeded
  and every output still has the recorded hash.
* ``runs/<step id>/stdout.txt``, ``stderr.txt`` and ``run.json`` (argv, shell
  command, cwd, exit code, peak RSS in kB, user and system CPU seconds, wall
  seconds, start and end timestamps, key, skipped). When a step is skipped its
  ``run.json`` keeps the figures of the execution that produced the outputs
  and gains ``"skipped": true`` and ``"checked"``.

The key of a step is the sha256 of the canonical JSON of: tool, tool version,
subcommand, resolved arguments (the injected thread count excluded, since it
does not change outputs), the sha256 of every file the step reads, and the keys
of the steps it depends on.
"""

import datetime
import glob
import hashlib
import json
import os
import queue
import signal
import subprocess
import sys
import threading
import time

from .plan import make_plans
from .validate import validate_workflow

EXIT_OK = 0
EXIT_FAILED = 1
EXIT_USAGE = 2
EXIT_CANCELLED = 130

KEY_VERSION = 1
STATE_FORMAT = "hicexplorer-workflow-state"

SUCCEEDED, SKIPPED, FAILED, CANCELLED, NOT_RUN = "succeeded", "skipped", "failed", "cancelled", "not_run"


def _now():
    return datetime.datetime.now(datetime.timezone.utc).astimezone().isoformat(timespec="milliseconds")


def _write_json(path, data):
    tmp = path + ".tmp"
    with open(tmp, "w") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp, path)


def canonical_json(data):
    return json.dumps(data, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


class Hasher:
    """Content hashes of files, directories and prefixes, memoised per run."""

    def __init__(self):
        self._memo = {}

    def file(self, path):
        try:
            st = os.stat(path)
        except OSError:
            return None
        sig = (st.st_size, st.st_mtime_ns, st.st_ino)
        cached = self._memo.get(path)
        if cached and cached[0] == sig:
            return cached[1]
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
        value = digest.hexdigest()
        self._memo[path] = (sig, value)
        return value

    def path(self, path, kind="file"):
        if kind == "prefix":
            entries = sorted(p for p in glob.glob(glob.escape(path) + "*") if os.path.isfile(p))
            if not entries:
                return None
            digest = hashlib.sha256()
            for entry in entries:
                digest.update(canonical_json([os.path.basename(entry), self.file(entry)]).encode())
            return digest.hexdigest()
        if os.path.isdir(path):
            digest = hashlib.sha256()
            for root, dirs, files in os.walk(path):
                dirs.sort()
                rel_root = os.path.relpath(root, path)
                digest.update(canonical_json(["d", rel_root]).encode())
                for name in sorted(files):
                    full = os.path.join(root, name)
                    digest.update(canonical_json(["f", os.path.join(rel_root, name), self.file(full)]).encode())
            return digest.hexdigest()
        if kind == "directory":
            return None
        return self.file(path)


class _Running:
    __slots__ = ("step", "proc", "start", "start_ts", "terminated", "stdout", "stderr", "key", "argv")


class Runner:
    """Runs a validated workflow. ``run()`` returns the exit status."""

    def __init__(self, workflow, loader, threads=None, force=False, grace_seconds=5.0,
                 log=None, poll_interval=0.05):
        self.workflow = workflow
        self.loader = loader
        self.budget = workflow.budget(threads)
        self.force = force
        self.grace_seconds = grace_seconds
        self.poll_interval = poll_interval
        self.log = log if log is not None else (lambda msg: print(msg, file=sys.stderr, flush=True))
        self.status = {}
        self.keys = {}
        self.messages = []
        self._cancel = threading.Event()
        self._hasher = Hasher()

    # -- public ---------------------------------------------------------
    def cancel(self):
        """Request cancellation; safe to call from a signal handler or another thread."""
        self._cancel.set()

    @property
    def cancelled(self):
        return self._cancel.is_set()

    def run(self):
        wf = self.workflow
        self.messages = validate_workflow(wf, self.loader)
        errors = [m for m in self.messages if m.level == "error"]
        for msg in self.messages:
            self.log(str(msg))
        if errors:
            return EXIT_FAILED
        self.plans = make_plans(wf, self.loader)
        os.makedirs(os.path.join(wf.state_dir, "runs"), exist_ok=True)
        self.state_path = os.path.join(wf.state_dir, "state.json")
        self.state = self._load_state()
        return self._schedule()

    # -- state ----------------------------------------------------------
    def _load_state(self):
        try:
            with open(self.state_path) as handle:
                data = json.load(handle)
            if data.get("format") == STATE_FORMAT and isinstance(data.get("steps"), dict):
                return data
        except (OSError, ValueError):
            pass
        return {"format": STATE_FORMAT, "version": 1, "steps": {}}

    def _save_state(self):
        self.state["workflow"] = self.workflow.name
        _write_json(self.state_path, self.state)

    def _run_dir(self, step):
        path = os.path.join(self.workflow.state_dir, "runs", step.id)
        os.makedirs(path, exist_ok=True)
        return path

    # -- keys -----------------------------------------------------------
    def _key(self, step):
        plan = self.plans[step.id]
        wf = self.workflow
        inputs = {path: self._hasher.path(wf.workdir_path(path), kind)
                  for path, kind in sorted(plan.inputs.items())}
        payload = {
            "engine_key_version": KEY_VERSION,
            "tool": step.tool,
            "tool_version": plan.spec.version,
            "subcommand": step.subcommand,
            "args": plan.args,
            "inputs": inputs,
            "deps": {dep: self.keys[dep] for dep in sorted(step.deps)},
        }
        return hashlib.sha256(canonical_json(payload).encode()).hexdigest()

    def _output_hashes(self, step):
        plan = self.plans[step.id]
        return {name: self._hasher.path(os.path.join(self.workflow.workdir, rel), plan.output_kind(name))
                for name, rel in step.outputs.items()}

    def _can_skip(self, step, key):
        if self.force:
            return False
        stored = self.state["steps"].get(step.id)
        if not stored or stored.get("status") != SUCCEEDED or stored.get("key") != key:
            return False
        recorded = stored.get("outputs") or {}
        if set(recorded) != set(step.outputs):
            return False
        current = self._output_hashes(step)
        return all(current[name] is not None and current[name] == recorded[name].get("sha256")
                   for name in step.outputs)

    # -- scheduling -----------------------------------------------------
    def _schedule(self):
        wf = self.workflow
        pending = sorted(wf.steps, key=lambda s: s.index)
        running = {}
        events = queue.Queue()
        failed = False
        kill_at = None
        done_ok = (SUCCEEDED, SKIPPED)

        while True:
            if self._cancel.is_set() and kill_at is None:
                kill_at = time.monotonic() + self.grace_seconds
                if running:
                    self.log("cancelling: terminating {} running step(s)".format(len(running)))
                for item in running.values():
                    item.terminated = True
                    self._signal(item, signal.SIGTERM)
            if kill_at is not None and running and time.monotonic() >= kill_at:
                for item in running.values():
                    self._signal(item, signal.SIGKILL)
                kill_at = float("inf")

            if not failed and not self._cancel.is_set():
                progressed = True
                while progressed and not self._cancel.is_set():
                    progressed = False
                    for step in list(pending):
                        dep_status = [self.status.get(d) for d in step.deps]
                        if any(s in (FAILED, CANCELLED, NOT_RUN) for s in dep_status):
                            pending.remove(step)
                            self.status[step.id] = NOT_RUN
                            progressed = True
                            continue
                        if not all(s in done_ok for s in dep_status):
                            continue
                        if step.id not in self.keys:
                            self.keys[step.id] = self._key(step)
                        if self._can_skip(step, self.keys[step.id]):
                            pending.remove(step)
                            self._record_skip(step)
                            progressed = True
                            continue
                        in_use = sum(r.step.threads for r in running.values())
                        if running and in_use + step.threads > self.budget:
                            break  # first ready step that does not fit blocks later ones
                        pending.remove(step)
                        if self._start(step, running, events):
                            progressed = True
                        else:
                            failed = True
                            break

            if not running:
                if not pending or failed or self._cancel.is_set():
                    break
            try:
                sid, wait_status, rusage, end = events.get(timeout=self.poll_interval)
            except queue.Empty:
                continue
            item = running.pop(sid)
            if not self._finish(item, wait_status, rusage, end):
                if not item.terminated:
                    failed = True

        for step in pending:
            self.status[step.id] = CANCELLED if self._cancel.is_set() else NOT_RUN
        self._summary()
        if self._cancel.is_set():
            return EXIT_CANCELLED
        if failed or any(s == FAILED for s in self.status.values()):
            return EXIT_FAILED
        return EXIT_OK

    def _signal(self, item, sig):
        if item.proc.returncode is not None:
            return
        try:
            os.killpg(item.proc.pid, sig)
        except (ProcessLookupError, PermissionError):
            pass

    # -- one step -------------------------------------------------------
    def _start(self, step, running, events):
        wf = self.workflow
        plan = self.plans[step.id]
        run_dir = self._run_dir(step)
        for rel in plan.mkdirs:
            os.makedirs(os.path.join(wf.workdir, rel), exist_ok=True)
        key = self.keys[step.id]
        self.state["steps"][step.id] = {"status": "running", "key": key, "started": _now()}
        self._save_state()

        item = _Running()
        item.step, item.key, item.argv, item.terminated = step, key, plan.argv, False
        item.stdout = open(os.path.join(run_dir, "stdout.txt"), "wb")
        item.stderr = open(os.path.join(run_dir, "stderr.txt"), "wb")
        item.start_ts, item.start = _now(), time.monotonic()
        self.log("[start] {}: {}".format(step.id, plan.command))
        try:
            item.proc = subprocess.Popen(plan.argv, cwd=wf.workdir, stdin=subprocess.DEVNULL,
                                         stdout=item.stdout, stderr=item.stderr,
                                         start_new_session=True)
        except OSError as exc:
            item.stdout.close()
            item.stderr.close()
            self._record(item, None, None, _now(), time.monotonic(), FAILED,
                         error="cannot start: {}".format(exc))
            return False
        running[step.id] = item

        def wait():
            _pid, wait_status, rusage = os.wait4(item.proc.pid, 0)
            end = (_now(), time.monotonic())
            item.proc.returncode = os.waitstatus_to_exitcode(wait_status)
            events.put((step.id, wait_status, rusage, end))

        threading.Thread(target=wait, name="wait-" + step.id, daemon=True).start()
        return True

    def _finish(self, item, wait_status, rusage, end):
        item.stdout.close()
        item.stderr.close()
        step = item.step
        exit_code = os.waitstatus_to_exitcode(wait_status)
        error = None
        if item.terminated:
            status = CANCELLED
        elif exit_code != 0:
            status = FAILED
            error = "exit code {}".format(exit_code)
        else:
            missing = [n for n, h in self._output_hashes(step).items() if h is None]
            if missing:
                status = FAILED
                error = "declared outputs missing after the run: " + ", ".join(
                    "{} ({})".format(n, step.outputs[n]) for n in missing)
            else:
                status = SUCCEEDED
        self._record(item, exit_code, rusage, end[0], end[1], status, error)
        return status == SUCCEEDED

    def _record(self, item, exit_code, rusage, end_ts, end_mono, status, error=None):
        step = item.step
        self.status[step.id] = status
        peak_rss = user = system = None
        if rusage is not None:
            peak_rss = rusage.ru_maxrss // 1024 if sys.platform == "darwin" else rusage.ru_maxrss
            user, system = rusage.ru_utime, rusage.ru_stime
        record = {
            "step": step.id,
            "tool": step.tool,
            "tool_version": self.plans[step.id].spec.version,
            "subcommand": step.subcommand,
            "argv": item.argv,
            "command": self.plans[step.id].command,
            "cwd": self.workflow.workdir,
            "threads": step.threads,
            "exit_code": exit_code,
            "status": status,
            "error": error,
            "peak_rss_kb": peak_rss,
            "user_seconds": user,
            "sys_seconds": system,
            "cpu_seconds": None if user is None else user + system,
            "wall_seconds": round(end_mono - item.start, 6),
            "start": item.start_ts,
            "end": end_ts,
            "key": item.key,
            "skipped": False,
        }
        _write_json(os.path.join(self._run_dir(step), "run.json"), record)
        entry = {"status": status, "key": item.key, "finished": end_ts}
        if status == SUCCEEDED:
            entry["outputs"] = {name: {"path": step.outputs[name], "sha256": digest}
                                for name, digest in self._output_hashes(step).items()}
        self.state["steps"][step.id] = entry
        self._save_state()
        self.log("[{}] {} ({}{:.2f} s)".format(
            status, step.id, "" if exit_code is None else "exit {}, ".format(exit_code),
            record["wall_seconds"]) + (": " + error if error and status == FAILED else ""))

    def _record_skip(self, step):
        self.status[step.id] = SKIPPED
        path = os.path.join(self._run_dir(step), "run.json")
        try:
            with open(path) as handle:
                record = json.load(handle)
        except (OSError, ValueError):
            plan = self.plans[step.id]
            record = {"step": step.id, "tool": step.tool, "subcommand": step.subcommand,
                      "argv": plan.argv, "command": plan.command, "cwd": self.workflow.workdir,
                      "exit_code": None}
        record.update({"skipped": True, "checked": _now(), "key": self.keys[step.id]})
        _write_json(path, record)
        self.log("[skipped] {} (up to date)".format(step.id))

    def _summary(self):
        counts = {}
        for step in self.workflow.steps:
            state = self.status.get(step.id, NOT_RUN)
            counts[state] = counts.get(state, 0) + 1
        self.log("workflow {}: ".format(self.workflow.name) + ", ".join(
            "{} {}".format(n, s) for s, n in sorted(counts.items())))
