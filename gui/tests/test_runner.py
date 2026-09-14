import json
import os
import shlex
import signal
import subprocess
import sys
import threading
import time

import pytest
import yaml

from conftest import GUI_DIR, read_log, write_workflow
from hicexplorer_gui.workflow import (EXIT_CANCELLED, EXIT_FAILED, EXIT_OK, Runner, SpecLoader,
                                      load_workflow)


def run(path, tools_dir, **kwargs):
    wf = load_workflow(path, kwargs.pop("workdir", None))
    runner = Runner(wf, SpecLoader(tools_dir), log=lambda msg: None, **kwargs)
    return runner.run(), runner


def started(log_path, since=0):
    return [e["tool"] + ":" + str(e["value"]) for e in read_log(log_path)[since:] if e["event"] == "start"]


def run_json(path, step):
    with open(os.path.join(os.path.dirname(path), ".hicexplorer-workflow", "runs", step, "run.json")) as fh:
        return json.load(fh)


def state(path):
    with open(os.path.join(os.path.dirname(path), ".hicexplorer-workflow", "state.json")) as fh:
        return json.load(fh)


def edit(path, change):
    with open(path) as fh:
        data = yaml.safe_load(fh)
    change(data)
    with open(path, "w") as fh:
        yaml.safe_dump(data, fh, sort_keys=False)


def statuses(runner):
    return {k: v for k, v in sorted(runner.status.items())}


def test_run_outputs_and_records(chain, tools_dir, fake_log):
    code, runner = run(chain, tools_dir)
    assert code == EXIT_OK
    assert statuses(runner) == {"a": "succeeded", "b": "succeeded", "c": "succeeded", "d": "succeeded"}
    work = os.path.dirname(chain)
    for step in "abcd":
        assert os.path.getsize(os.path.join(work, step, "out.txt")) == 65
    rec = run_json(chain, "b")
    assert rec["argv"][0] == os.path.join(tools_dir, "combine")
    assert rec["argv"][1:] == ["--inputs", "a/out.txt", "data/x.txt", "--outFileName", "b/out.txt",
                               "--mode", "max", "--threads", "2"]
    assert shlex.split(rec["command"]) == rec["argv"]
    assert rec["cwd"] == work
    assert rec["exit_code"] == 0 and rec["status"] == "succeeded" and rec["skipped"] is False
    assert rec["peak_rss_kb"] > 1000
    assert rec["cpu_seconds"] > 0 and rec["cpu_seconds"] == pytest.approx(rec["user_seconds"] + rec["sys_seconds"])
    assert rec["wall_seconds"] > 0 and rec["start"] < rec["end"]
    assert len(rec["key"]) == 64 and rec["key"] == state(chain)["steps"]["b"]["key"]
    runs = os.path.join(work, ".hicexplorer-workflow", "runs", "b")
    assert open(os.path.join(runs, "stdout.txt")).read() == "processing combine\n"
    assert open(os.path.join(runs, "stderr.txt")).read() == "log line for combine\n"
    # the logged argv reproduces the output
    other = os.path.join(work, "rerun")
    os.makedirs(os.path.join(other, "b"))
    os.makedirs(os.path.join(other, "a"))
    os.makedirs(os.path.join(other, "data"))
    for rel in ("a/out.txt", "data/x.txt"):
        with open(os.path.join(work, rel), "rb") as src, open(os.path.join(other, rel), "wb") as dst:
            dst.write(src.read())
    subprocess.check_call(rec["argv"], cwd=other)
    assert open(os.path.join(other, "b/out.txt"), "rb").read() == open(os.path.join(work, "b/out.txt"), "rb").read()


def test_resume_skips_everything(chain, tools_dir, fake_log):
    run(chain, tools_dir)
    before = len(read_log(fake_log))
    code, runner = run(chain, tools_dir)
    assert code == EXIT_OK
    assert set(runner.status.values()) == {"skipped"}
    assert len(read_log(fake_log)) == before
    rec = run_json(chain, "a")
    assert rec["skipped"] is True and rec["exit_code"] == 0 and "checked" in rec
    code, runner = run(chain, tools_dir, force=True)
    assert set(runner.status.values()) == {"succeeded"}


def test_changed_parameter_reruns_step_and_dependents(chain, tools_dir, fake_log):
    run(chain, tools_dir)
    edit(chain, lambda d: d["steps"][1]["args"].update(mode="sum"))
    code, runner = run(chain, tools_dir)
    assert code == EXIT_OK
    assert statuses(runner) == {"a": "skipped", "b": "succeeded", "c": "skipped", "d": "succeeded"}


def test_changed_input_content_reruns_dependents(chain, tools_dir, fake_log):
    run(chain, tools_dir)
    work = os.path.dirname(chain)
    with open(os.path.join(work, "data", "y.txt"), "w") as fh:
        fh.write("changed\n")
    code, runner = run(chain, tools_dir)
    assert statuses(runner) == {"a": "skipped", "b": "skipped", "c": "succeeded", "d": "skipped"}
    with open(os.path.join(work, "data", "x.txt"), "w") as fh:
        fh.write("changed\n")
    code, runner = run(chain, tools_dir)
    assert statuses(runner) == {"a": "skipped", "b": "succeeded", "c": "skipped", "d": "succeeded"}


def test_thread_count_is_not_part_of_the_key(chain, tools_dir, fake_log):
    run(chain, tools_dir)
    edit(chain, lambda d: d["steps"][1].update(threads=3))
    code, runner = run(chain, tools_dir, threads=1)
    assert set(runner.status.values()) == {"skipped"}


@pytest.mark.parametrize("action", ["delete", "modify"])
def test_deleted_or_modified_output_reruns_that_step(chain, tools_dir, fake_log, action):
    run(chain, tools_dir)
    target = os.path.join(os.path.dirname(chain), "b", "out.txt")
    if action == "delete":
        os.remove(target)
    else:
        with open(target, "w") as fh:
            fh.write("tampered\n")
    code, runner = run(chain, tools_dir)
    assert code == EXIT_OK
    # b rewrites the same content, so d's inputs are unchanged and d stays skipped
    assert statuses(runner) == {"a": "skipped", "b": "succeeded", "c": "skipped", "d": "skipped"}


def test_prefix_and_directory_outputs_tracked(tmp_path, tools_dir, fake_log):
    work = tmp_path / "w"
    (work / "data").mkdir(parents=True)
    (work / "data" / "m.txt").write_text("m")
    data = {"version": 1, "inputs": {"m": "data/m.txt"}, "steps": [
        {"id": "p", "tool": "prefixtool", "threads": 2, "args": {"matrix": "${inputs.m}", "outPrefix": "${outputs.t}"},
         "outputs": {"t": "tads/sub/t"}},
        {"id": "q", "tool": "dirtool", "args": {"inputs": ["tads/sub/t_a.txt"], "QCfolder": "${outputs.qc}",
                                                 "summary": "${outputs.s}"},
         "after": ["p"], "outputs": {"qc": "qc/folder", "s": "qc/summary.txt"}}]}
    path = write_workflow(work, data)
    code, runner = run(path, tools_dir)
    assert code == EXIT_OK, runner.messages
    assert os.path.isfile(str(work / "tads/sub/t_b.txt")) and os.path.isfile(str(work / "qc/folder/report.txt"))
    assert [e["nproc"] for e in read_log(fake_log) if e["event"] == "start" and e["tool"] == "prefixtool"] == [2]
    code, runner = run(path, tools_dir)
    assert set(runner.status.values()) == {"skipped"}
    (work / "tads/sub/t_b.txt").write_text("x")
    code, runner = run(path, tools_dir)
    assert statuses(runner) == {"p": "succeeded", "q": "skipped"}
    (work / "qc/folder/report.txt").unlink()
    code, runner = run(path, tools_dir)
    assert statuses(runner) == {"p": "skipped", "q": "succeeded"}


def budget_workflow(threads_list, sleep=0.4, budget=4):
    return {"version": 1, "threads": budget, "steps": [
        {"id": "s{}".format(i), "tool": "makeA", "threads": t,
         "args": {"outFileName": "${outputs.o}", "value": i, "sleep": sleep},
         "outputs": {"o": "out/{}.txt".format(i)}} for i, t in enumerate(threads_list)]}


def intervals(log_path):
    events = read_log(log_path)
    starts = {e["value"]: e for e in events if e["event"] == "start"}
    ends = {e["value"]: e for e in events if e["event"] == "end"}
    return [(starts[v]["t"], ends[v]["t"], starts[v]["threads"]) for v in starts]


def max_concurrent(spans):
    points = sorted([(s, t) for s, e, t in spans] + [(e, -t) for s, e, t in spans], key=lambda p: (p[0], p[1]))
    current = peak = 0
    for _, delta in points:
        current += delta
        peak = max(peak, current)
    return peak


def test_thread_budget_never_exceeded(tmp_path, tools_dir, fake_log):
    path = write_workflow(tmp_path / "w", budget_workflow([2, 1, 2, 1, 3, 1, 2]))
    code, runner = run(path, tools_dir)
    assert code == EXIT_OK
    spans = intervals(fake_log)
    assert len(spans) == 7
    assert max_concurrent(spans) <= 4
    assert max_concurrent([(s, e, 1) for s, e, _ in spans]) >= 2  # it did run steps concurrently
    assert sorted(t for _, _, t in spans) == [1, 1, 1, 2, 2, 2, 3]


def test_oversized_step_runs_alone(tmp_path, tools_dir, fake_log):
    path = write_workflow(tmp_path / "w", budget_workflow([1, 6, 1], sleep=0.3, budget=4))
    code, runner = run(path, tools_dir)
    assert code == EXIT_OK
    spans = {t: (s, e) for s, e, t in intervals(fake_log) if t == 6}
    big_start, big_end = spans[6]
    for s, e, t in intervals(fake_log):
        if t != 6:
            assert e <= big_start or s >= big_end


def test_threads_override(tmp_path, tools_dir, fake_log):
    path = write_workflow(tmp_path / "w", budget_workflow([1, 1, 1, 1], sleep=0.3, budget=4))
    run(path, tools_dir, threads=1)
    assert max_concurrent([(s, e, 1) for s, e, _ in intervals(fake_log)]) == 1


def test_failure_stops_scheduling(tmp_path, tools_dir, fake_log):
    data = {"version": 1, "threads": 2, "steps": [
        {"id": "bad", "tool": "makeA", "args": {"outFileName": "${outputs.o}", "value": 1, "sleep": 0.2, "fail": True},
         "outputs": {"o": "bad.txt"}},
        {"id": "slow", "tool": "makeA", "args": {"outFileName": "${outputs.o}", "value": 2, "sleep": 1.0},
         "outputs": {"o": "slow.txt"}},
        {"id": "later", "tool": "makeA", "args": {"outFileName": "${outputs.o}", "value": 3},
         "outputs": {"o": "later.txt"}},
        {"id": "child", "tool": "combine", "args": {"inputs": ["${steps.bad.outputs.o}"], "outFileName": "${outputs.o}"},
         "outputs": {"o": "child.txt"}}]}
    path = write_workflow(tmp_path / "w", data)
    code, runner = run(path, tools_dir)
    assert code == EXIT_FAILED
    assert statuses(runner) == {"bad": "failed", "child": "not_run", "later": "not_run", "slow": "succeeded"}
    # bad and slow fit the budget together and start concurrently, so the
    # order in which the two processes log their start is not defined.
    assert sorted(started(fake_log)) == ["makeA:1", "makeA:2"]
    rec = run_json(path, "bad")
    assert rec["exit_code"] == 3 and rec["status"] == "failed"
    assert "failing on purpose" in open(os.path.join(tmp_path, "w", ".hicexplorer-workflow", "runs", "bad",
                                                     "stderr.txt")).read()
    assert state(path)["steps"]["bad"]["status"] == "failed"


def test_missing_declared_output_fails(tmp_path, tools_dir, fake_log):
    data = {"version": 1, "steps": [
        {"id": "s", "tool": "makeA", "args": {"outFileName": "${outputs.o}"},
         "outputs": {"o": "o.txt", "ghost": "ghost.txt"}}]}
    path = write_workflow(tmp_path / "w", data)
    code, runner = run(path, tools_dir)
    assert code == EXIT_FAILED
    assert "ghost" in run_json(path, "s")["error"]


def _wait_for_start(log_path, count=1, timeout=20):
    deadline = time.time() + timeout
    while time.time() < deadline:
        events = [e for e in read_log(log_path) if e["event"] == "start"]
        if len(events) >= count:
            return events
        time.sleep(0.05)
    raise AssertionError("tool did not start")


def _alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


@pytest.mark.parametrize("ignore_term", [False, True])
def test_cancel_terminates_children(tmp_path, tools_dir, fake_log, ignore_term):
    data = {"version": 1, "threads": 4, "steps": [
        {"id": "sleeper", "tool": "makeA",
         "args": {"outFileName": "${outputs.o}", "sleep": 60, "ignoreTerm": ignore_term},
         "outputs": {"o": "sleep.txt"}},
        {"id": "after", "tool": "combine", "args": {"inputs": ["${steps.sleeper.outputs.o}"],
                                                    "outFileName": "${outputs.o}"},
         "outputs": {"o": "after.txt"}}]}
    path = write_workflow(tmp_path / "w", data)
    wf = load_workflow(path)
    runner = Runner(wf, SpecLoader(tools_dir), log=lambda m: None, grace_seconds=0.5)
    result = {}
    thread = threading.Thread(target=lambda: result.update(code=runner.run()))
    t0 = time.time()
    thread.start()
    pid = _wait_for_start(fake_log)[0]["pid"]
    runner.cancel()
    thread.join(20)
    assert not thread.is_alive()
    assert result["code"] == EXIT_CANCELLED
    assert time.time() - t0 < 15
    assert not _alive(pid)
    assert statuses(runner) == {"after": "cancelled", "sleeper": "cancelled"}
    assert state(path)["steps"]["sleeper"]["status"] == "cancelled"
    assert run_json(path, "sleeper")["exit_code"] == -(signal.SIGKILL if ignore_term else signal.SIGTERM)
    # the next run does not treat it as done
    edit(path, lambda d: d["steps"][0]["args"].update(sleep=0))
    code, runner = run(path, tools_dir)
    assert code == EXIT_OK and statuses(runner) == {"after": "succeeded", "sleeper": "succeeded"}


def test_cli_sigint_exits_130(tmp_path, tools_dir, fake_log):
    data = {"version": 1, "steps": [
        {"id": "sleeper", "tool": "makeA", "args": {"outFileName": "${outputs.o}", "sleep": 60},
         "outputs": {"o": "sleep.txt"}}]}
    path = write_workflow(tmp_path / "w", data)
    env = dict(os.environ, PYTHONPATH=GUI_DIR)
    proc = subprocess.Popen([sys.executable, "-m", "hicexplorer_gui.workflow.cli", "run", path,
                             "--tools-dir", tools_dir], env=env,
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        pid = _wait_for_start(fake_log)[0]["pid"]
        proc.send_signal(signal.SIGINT)
        assert proc.wait(20) == 130
    finally:
        if proc.poll() is None:
            proc.kill()
    deadline = time.time() + 5
    while _alive(pid) and time.time() < deadline:
        time.sleep(0.05)
    assert not _alive(pid)
    assert state(path)["steps"]["sleeper"]["status"] == "cancelled"
