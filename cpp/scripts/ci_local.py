#!/usr/bin/env python3
"""Run the GitHub Actions jobs of .github/workflows on this machine.

    cpp/scripts/ci_local.py [--workflow cpp|python]... [--python 3.12]...
                            [--work-dir DIR] [--worktree] [--reuse-envs]
                            [--keep-going] [--list]

The steps are read from the workflow files themselves, so the pins and
commands cannot drift from what CI runs. Every job whose runner matches this
machine (ubuntu-24.04 on x86-64, ubuntu-24.04-arm on aarch64, macos-14 on
arm64 macOS) runs in a fresh conda environment of its own, on a copy of the
committed tree (git archive HEAD) or, with --worktree, of the working tree.
The two `uses:` actions the jobs need are emulated: actions/checkout is that
copy, conda-incubator/setup-miniconda is `mamba create` with the channels and
Python version the step names. Jobs for other runners are listed as skipped.

Each step's output goes to <work-dir>/logs/<job>/<nn>-<step>.log; the summary
at the end gives the result and the time of every step, and the exit code is
non-zero when any step failed.

Needs mamba (or conda) on PATH and PyYAML in the Python that runs this script.
"""
import argparse
import hashlib
import os
import platform
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import time
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
WORKFLOWS = ROOT / ".github" / "workflows"


def host_runners():
    machine = platform.machine().lower()
    if sys.platform == "darwin":
        return {"macos-14", "macos-latest"} if machine == "arm64" else {"macos-13"}
    if machine in ("aarch64", "arm64"):
        return {"ubuntu-24.04-arm"}
    return {"ubuntu-24.04", "ubuntu-latest"}


def expand_matrix(job):
    strategy = job.get("strategy", {}).get("matrix")
    if not strategy:
        return [{}]
    include = strategy.get("include")
    axes = {k: v for k, v in strategy.items() if k not in ("include", "exclude")}
    combos = [{}]
    for key, values in axes.items():
        combos = [dict(c, **{key: v}) for c in combos for v in values]
    if include:
        combos = include if combos == [{}] else combos + include
    return combos


def substitute(text, context):
    def value(match):
        expr = match.group(1).strip()
        node = context
        for part in expr.split("."):
            if not isinstance(node, dict) or part not in node:
                raise KeyError("unsupported expression ${{ %s }}" % expr)
            node = node[part]
        return str(node)
    return re.sub(r"\$\{\{\s*([^}]+?)\s*\}\}", value, str(text))


def find_conda():
    """mamba (preferred) or conda, and the base holding etc/profile.d/conda.sh."""
    exe = shutil.which("mamba") or shutil.which("conda")
    if not exe:
        sys.exit("ci_local: neither mamba nor conda is on PATH")
    candidates = [Path(exe).resolve().parents[1]]
    if os.environ.get("CONDA_EXE"):
        candidates.append(Path(os.environ["CONDA_EXE"]).resolve().parents[1])
    for tool in ("conda", "mamba"):
        if shutil.which(tool):
            out = subprocess.run([tool, "info", "--base"], capture_output=True, text=True).stdout.strip()
            # mamba 2 prints "base environment : <path>"
            candidates.append(Path(out.rsplit(" : ", 1)[-1].strip()))
    for base in candidates:
        if (base / "etc" / "profile.d" / "conda.sh").is_file():
            return exe, base
    sys.exit("ci_local: no etc/profile.d/conda.sh found for %s" % exe)


def copy_tree(dest, worktree):
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    if worktree:
        files = subprocess.run(["git", "ls-files", "-co", "--exclude-standard", "-z"], cwd=ROOT,
                               capture_output=True, check=True).stdout.split(b"\0")
        for name in filter(None, files):
            src = ROOT / name.decode()
            if src.is_file():
                target = dest / name.decode()
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, target)
    else:
        archive = dest.parent / (dest.name + ".tar")
        subprocess.run(["git", "archive", "--format=tar", "-o", str(archive), "HEAD"], cwd=ROOT, check=True)
        with tarfile.open(archive) as tar:
            tar.extractall(dest, filter="data")
        archive.unlink()


def rewrite_install(command, env_name, prefix):
    """Point mamba/conda install commands at the job's own prefix."""
    out = []
    for line in command.splitlines():
        words = line.split()
        if len(words) >= 2 and words[0] in ("mamba", "conda") and words[1] in ("install", "create", "update"):
            line = re.sub(r"(^|\s)-n\s+%s(\s|$)" % re.escape(env_name), r"\1\2", line)
            line = line.replace(words[1], "%s -p %s --override-channels -c conda-forge -c bioconda"
                                % (words[1], shlex.quote(str(prefix))), 1)
        out.append(line)
    return "\n".join(out)


def run_step(script, cwd, env, log, conda_base, prefix):
    activate = "source %s && conda activate %s || { echo 'ci_local: activation failed'; exit 97; }" % (
        shlex.quote(str(conda_base / "etc" / "profile.d" / "conda.sh")), shlex.quote(str(prefix)))
    full = "set -eo pipefail\n%s\n%s\n" % (activate, script)
    start = time.time()
    with open(log, "w") as handle:
        handle.write("$ " + script.replace("\n", "\n$ ") + "\n\n")
        handle.flush()
        code = subprocess.run(["bash", "--noprofile", "--norc", "-c", full], cwd=cwd, env=env,
                              stdout=handle, stderr=subprocess.STDOUT).returncode
    return code, time.time() - start


def slug(text):
    return re.sub(r"[^A-Za-z0-9]+", "-", text).strip("-")[:60]


def run_job(wf_name, job_id, job, matrix, args, conda, results):
    exe, conda_base = conda
    context = {"matrix": matrix}
    title = substitute(job.get("name", job_id), context)
    label = "%s/%s" % (wf_name, title)
    tag = slug("%s-%s" % (wf_name, title))
    job_dir = args.work_dir / "jobs" / tag
    logs = args.work_dir / "logs" / tag
    logs.mkdir(parents=True, exist_ok=True)
    workspace = job_dir / "src"
    temp = job_dir / "tmp"
    prefix = args.work_dir / "envs" / tag
    copy_tree(workspace, args.worktree)
    temp.mkdir(parents=True, exist_ok=True)
    context["github"] = {"workspace": str(workspace)}
    context["runner"] = {"temp": str(temp)}

    env = dict(os.environ)
    for var in ("CONDA_PREFIX", "CONDA_DEFAULT_ENV", "CONDA_SHLVL", "PYTHONPATH", "VIRTUAL_ENV",
                "HICX_BIN_DIR", "HICX_PLOT_PYTHON", "CMAKE_PREFIX_PATH", "CC", "CXX"):
        env.pop(var, None)
    env.update({"GITHUB_WORKSPACE": str(workspace), "RUNNER_TEMP": str(temp), "CI": "true",
                "GITHUB_ACTIONS": "false"})

    install_hash = hashlib.sha256()
    env_name = None
    steps = []
    for step in job["steps"]:
        uses = step.get("uses", "")
        if uses.startswith("actions/checkout"):
            continue
        if uses.startswith("conda-incubator/setup-miniconda"):
            with_ = step.get("with", {})
            env_name = with_.get("activate-environment", "test")
            python = substitute(with_.get("python-version", ""), context)
            channels = with_.get("channels", "conda-forge").split(",")
            create = [exe, "create", "--yes", "-p", str(prefix), "--override-channels"]
            for channel in channels:
                create += ["-c", channel.strip()]
            create.append("python=%s" % python if python else "python")
            steps.append(("Set up conda", create, True))
            install_hash.update(" ".join(create[1:]).encode())
            continue
        if uses:
            results.append((label, "uses " + uses, "unsupported", 0.0))
            return
        script = substitute(step["run"], context)
        if env_name:
            script = rewrite_install(script, env_name, prefix)
        step_env = {k: substitute(v, context) for k, v in step.get("env", {}).items()}
        name = step.get("name", script.splitlines()[0])
        is_install = bool(re.search(r"^\s*(mamba|conda) (install|create)", script, re.M))
        if is_install:
            install_hash.update(script.encode())
        steps.append((name, (script, step_env), is_install))

    marker = prefix / ".ci_local_install"
    reuse = args.reuse_envs and marker.exists() and marker.read_text() == install_hash.hexdigest()
    if not reuse and prefix.exists():
        shutil.rmtree(prefix)

    print("== %s" % label, flush=True)
    for number, (name, payload, is_install) in enumerate(steps, 1):
        log = logs / ("%02d-%s.log" % (number, slug(name)))
        if reuse and is_install:
            results.append((label, name, "reused", 0.0))
            print("   %-50s reused" % name, flush=True)
            continue
        if isinstance(payload, list):
            start = time.time()
            with open(log, "w") as handle:
                code = subprocess.run(payload, env=env, stdout=handle, stderr=subprocess.STDOUT).returncode
            seconds = time.time() - start
        else:
            script, step_env = payload
            code, seconds = run_step(script, workspace, dict(env, **step_env), log, conda_base, prefix)
        status = "pass" if code == 0 else "FAIL (exit %d, %s)" % (code, log)
        results.append((label, name, status, seconds))
        print("   %-50s %-6s %7.1f s" % (name, "pass" if code == 0 else "FAIL", seconds), flush=True)
        if code != 0:
            if not args.keep_going:
                return
            if is_install:
                return
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(install_hash.hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--workflow", action="append",
                        help="workflow file stem (cpp, python); default all")
    parser.add_argument("--python", action="append",
                        help="run only these matrix.python values")
    parser.add_argument("--work-dir", type=Path,
                        default=Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "hicx-ci-local",
                        help="environments, sources and logs (default %(default)s)")
    parser.add_argument("--worktree", action="store_true",
                        help="test the working tree instead of HEAD")
    parser.add_argument("--reuse-envs", action="store_true",
                        help="keep an environment whose install steps are unchanged")
    parser.add_argument("--keep-going", action="store_true",
                        help="run the remaining steps of a job after a failed one")
    parser.add_argument("--list", action="store_true", help="list the jobs and exit")
    args = parser.parse_args()
    args.work_dir = args.work_dir.resolve()

    if not args.worktree:
        dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT,
                               capture_output=True, text=True).stdout
        if dirty:
            print("ci_local: uncommitted changes are not tested (HEAD is); use --worktree to include them",
                  file=sys.stderr)

    runners = host_runners()
    planned, skipped = [], []
    for path in sorted(WORKFLOWS.glob("*.yml")):
        if args.workflow and path.stem not in args.workflow:
            continue
        workflow = yaml.safe_load(path.read_text())
        for job_id, job in workflow["jobs"].items():
            for matrix in expand_matrix(job):
                runs_on = substitute(job["runs-on"], {"matrix": matrix})
                title = substitute(job.get("name", job_id), {"matrix": matrix})
                if args.python and "python" in matrix and str(matrix["python"]) not in args.python:
                    continue
                entry = (path.stem, job_id, job, matrix)
                (planned if runs_on in runners else skipped).append((entry, title, runs_on))

    for (wf, _, _, _), title, runs_on in planned:
        print("run   %s/%s (%s)" % (wf, title, runs_on))
    for (wf, _, _, _), title, runs_on in skipped:
        print("skip  %s/%s (%s, not this machine)" % (wf, title, runs_on))
    if args.list or not planned:
        return 0

    conda = find_conda()
    results = []
    for (wf, job_id, job, matrix), _, _ in planned:
        run_job(wf, job_id, job, matrix, args, conda, results)

    print("\nSummary")
    failed = 0
    for label, step, status, seconds in results:
        if status not in ("pass", "reused"):
            failed += 1
        print("  %-40s %-48s %s%s" % (label, step, status, "" if not seconds else " (%.0f s)" % seconds))
    print("\n%d step(s) failed" % failed if failed else "\nall steps passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
