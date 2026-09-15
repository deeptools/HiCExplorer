"""pytest over the reference cache of equiv.py (reference_cache.py).

    python -m pytest cpp/scripts/test_equiv_cache.py

Runs the Python side of two synthetic cases in a throwaway repository, data
directory and cache: case A runs hicexplorer.toolA, which imports
hicexplorer.helper, on a.txt; case B runs hicexplorer.toolB, which imports
hicexplorer.other_helper, on b.txt. Every test first shows both cases run
fresh and then come from the cache, then changes one thing that can change a
Python result and checks that the affected case runs fresh again while an
unaffected one stays cached. The interpreter fingerprint is replaced by a
dictionary the tests control, so a package version can be faked.
"""

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import equiv  # noqa: E402
import reference_cache  # noqa: E402

TOOL = '''import sys
from hicexplorer.{helper} import transform


def main():
    source, target = sys.argv[1], sys.argv[2]
    with open(source) as handle:
        text = handle.read()
    with open(target, "w") as handle:
        handle.write(transform(text) + " ".join(sys.argv[3:]))
'''

LAUNCHER = '''#!/usr/bin/env python
from hicexplorer.{tool} import main

main()
'''


class World:
    def __init__(self, root, monkeypatch):
        self.repo = root / "repo"
        package = self.repo / "hicexplorer"
        package.mkdir(parents=True)
        (self.repo / "bin").mkdir()
        (package / "__init__.py").write_text("")
        (package / "helper.py").write_text("def transform(text):\n    return text.upper()\n")
        (package / "other_helper.py").write_text("def transform(text):\n    return text[::-1]\n")
        for tool, helper in (("toolA", "helper"), ("toolB", "other_helper")):
            (package / f"{tool}.py").write_text(TOOL.format(helper=helper))
            (self.repo / "bin" / tool).write_text(LAUNCHER.format(tool=tool))
        self.data = root / "data"
        self.data.mkdir()
        (self.data / "a.txt").write_text("alpha")
        (self.data / "b.txt").write_text("beta")
        self.tmp = root / "tmp"
        self.tmp.mkdir()
        self.fingerprint = {"python": "3.12.test", "packages": {"numpy": "1.26.4"},
                            "modules": {"fit_nbinom": None}, "path": []}
        monkeypatch.setattr(equiv, "REPO_ROOT", self.repo)
        monkeypatch.setattr(equiv, "PY_BIN", self.repo / "bin")
        monkeypatch.setattr(reference_cache, "interpreter_fingerprint",
                            lambda python, env, cache_dir, repo_root:
                            json.loads(json.dumps(self.fingerprint)))
        self.options = argparse.Namespace(py_python=sys.executable, noise_runs=5, cache="use",
                                          cache_dir=str(root / "cache"), tmpdir=str(self.tmp))
        self.cases = {
            "A": self._case("synthetic.A", "toolA", ["{data}/a.txt", "{out}/out.txt"]),
            "B": self._case("synthetic.B", "toolB", ["{data}/b.txt", "{out}/out.txt"]),
        }

    @staticmethod
    def _case(identifier, tool, args):
        return {"id": identifier, "tool": tool, "tier": 0, "args": args, "expect_exit": 0,
                "outputs": [{"path": "{out}/out.txt", "format": "text", "class": "E0"}],
                "large": False, "validators": [], "py_script": None, "path_prepend": []}

    def run(self, name, tmpdir=None):
        """True when the case's Python side came from the cache."""
        reference_cache.clear_hash_memo()
        case = self.cases[name]
        _, workdir = equiv.make_workdir(case["id"], tmpdir or self.tmp)
        (workdir / "out_cpp").mkdir()
        options = argparse.Namespace(**vars(self.options))
        side = equiv.run_python_side(case, options, workdir, str(self.data),
                                     equiv.python_environment())
        assert side["measurement"]["exit_code"] == 0, (workdir / "out_py" / "stderr.txt").read_text()
        output = workdir / "out_py" / "out.txt"
        assert output.is_file()
        self.last_output = output.read_text()
        self.last_side = side
        return side["cached"]

    def both_cached_after_first_run(self):
        assert self.run("A") is False
        assert self.run("B") is False
        assert self.run("A") is True
        assert self.run("B") is True


@pytest.fixture
def world(tmp_path, monkeypatch):
    return World(tmp_path, monkeypatch)


def test_a_cached_side_restores_the_outputs_and_the_measurement(world):
    assert world.run("A") is False
    fresh_output, fresh_measurement = world.last_output, dict(world.last_side["measurement"])
    assert world.run("A") is True
    assert world.last_output == fresh_output == "ALPHA"
    for name in ("cpu_seconds", "peak_rss_kb", "seconds", "exit_code"):
        assert world.last_side["measurement"][name] == fresh_measurement[name]
    assert world.last_side["processes"] == 0


def test_one_changed_input_byte(world):
    world.both_cached_after_first_run()
    (world.data / "a.txt").write_text("alphb")
    assert world.run("A") is False
    assert world.last_output == "ALPHB"
    assert world.run("B") is True


def test_one_changed_line_in_an_imported_source_file(world):
    world.both_cached_after_first_run()
    helper = world.repo / "hicexplorer" / "helper.py"
    # A different size and a later mtime, so the interpreter does not reuse
    # the bytecode it cached for the old line (it checks size and mtime).
    helper.write_text("def transform(text):\n    return text.lower()  # lower case\n")
    stat = helper.stat()
    os.utime(helper, ns=(stat.st_atime_ns, stat.st_mtime_ns + 5_000_000_000))
    assert world.run("A") is False
    assert world.last_output == "alpha"
    assert world.run("B") is True


def test_a_changed_argument(world):
    world.both_cached_after_first_run()
    world.cases["A"]["args"] = world.cases["A"]["args"] + ["extra"]
    assert world.run("A") is False
    assert world.last_output == "ALPHAextra"
    assert world.run("B") is True


def test_a_faked_package_version(world):
    world.both_cached_after_first_run()
    world.fingerprint["packages"]["numpy"] = "9.9.9"
    assert world.run("A") is False
    assert world.run("A") is True
    world.fingerprint["packages"]["numpy"] = "1.26.4"
    assert world.run("B") is True


def test_a_changed_harness_environment_variable(world, monkeypatch):
    world.both_cached_after_first_run()
    monkeypatch.setattr(equiv, "PYTHON_SIDE_ENVIRONMENT", {"COLUMNS": "81"})
    assert world.run("A") is False
    assert world.run("A") is True
    monkeypatch.setattr(equiv, "PYTHON_SIDE_ENVIRONMENT", {"COLUMNS": "80"})
    assert world.run("B") is True


def test_refresh_reruns_and_off_neither_reads_nor_writes(world):
    world.both_cached_after_first_run()
    world.options.cache = "refresh"
    assert world.run("A") is False
    world.options.cache = "off"
    assert world.run("B") is False
    assert world.last_side["key"] is None
    world.options.cache = "use"
    assert world.run("A") is True


def test_an_unexpected_exit_status_is_not_cached(world):
    world.cases["A"]["expect_exit"] = 1
    reference_cache.clear_hash_memo()
    case = world.cases["A"]
    for _ in range(2):
        workdir = Path(tempfile.mkdtemp(prefix=f"equiv-{case['id']}-", dir=world.tmp))
        side = equiv.run_python_side(case, argparse.Namespace(**vars(world.options)), workdir,
                                     str(world.data), equiv.python_environment())
        assert side["cached"] is False and side["stored"] is False


def test_an_embedded_working_directory_is_rewritten_into_another_tmpdir(world):
    world.cases["A"]["args"] = ["{data}/a.txt", "{out}/out.txt", "{out}"]
    assert world.run("A") is False
    first = world.last_output
    longer = world.tmp / "a_much_longer_temporary_directory_name"
    longer.mkdir()
    assert world.run("A", tmpdir=longer) is True
    # the embedded out_py path is the new working directory's, not the old one
    assert world.last_output != first
    assert world.last_output.startswith("ALPHA" + str(longer))
    assert world.last_output.endswith("/out_py")


def test_working_directories_have_one_path_length(tmp_path):
    lengths = set()
    for name in ("t", "tmp_cold", "tmp_warm2", "x" * 150):
        directory = tmp_path / name
        directory.mkdir()
        for case_id in ("short", "hicPCA.bedgraph.lieberman.geneTrack.mm9_reduced_chr1"):
            _, workdir = equiv.make_workdir(case_id, directory)
            assert workdir.is_dir()
            lengths.add(len(str(workdir)))
    assert lengths == {equiv.WORKDIR_PATH_LENGTH}


def test_scheduler_runs_thin_margin_cases_alone_and_auto_uses_half_the_physical_cores():
    options = argparse.Namespace(noise_runs=5, determinism=False)
    case = {"id": "x", "outputs": [], "large": False, "memory": {}}
    measured = {"py_peak_rss_kb": 100_000, "py_seconds": 10.0, "py_cpu_seconds": 10.0,
                "cpp_peak_rss_kb": 50_000, "cpp_seconds": 2.0, "cpp_cpu_seconds": 2.0}
    assert equiv._demand(case, options, dict(measured, time_ratio=0.2), 8, False)[1] == 1
    assert equiv._demand(case, options, dict(measured, time_ratio=0.6), 8, False)[1] == 8
    assert equiv._demand(case, options, dict(measured, time_ratio=0.6), 8, True)[1] == 8
    assert equiv.resolve_jobs("auto") == max(1, equiv.physical_cores() // 2)
    assert equiv.resolve_jobs("3") == 3


def test_no_history_drawing_and_thin_margin_tools_run_alone():
    options = argparse.Namespace(noise_runs=5, determinism=False)
    def case(tool, args=("x",), cpp_args=()):
        return {"id": tool, "tool": tool, "args": list(args), "cpp_args": list(cpp_args),
                "outputs": [], "large": False, "memory": {}}
    measured = {"py_peak_rss_kb": 100_000, "py_seconds": 10.0, "py_cpu_seconds": 10.0,
                "cpp_peak_rss_kb": 50_000, "cpp_seconds": 2.0, "cpp_cpu_seconds": 2.0,
                "time_ratio": 0.2}
    assert equiv._demand(case("hicPlotMatrix"), options, None, 8, False)[1] == 8
    assert equiv._demand(case("hicPlotMatrix"), options, measured, 8, False)[1] == 1
    assert equiv._demand(case("hicBuildMatrix"), options, None, 8, False)[1] == 8
    assert equiv._demand(case("hicTransform"), options, None, 8, False)[1] == 8
    assert equiv._demand(case("hicInfo"), options, None, 8, False)[1] == 1
    assert equiv._demand(case("hicCorrectMatrix", ["correct"]), options, None, 8, False)[1] == 1
    assert equiv._demand(case("hicCorrectMatrix", ["diagnostic_plot"]), options, None, 8, False)[1] == 8
    assert equiv._demand(case("hicCompartmentalization", cpp_args=["--noPlot"]), options, None, 8,
                         False)[1] == 1


def test_a_time_gate_failure_beside_other_cases_is_rerun_alone(tmp_path):
    """A fake runner: case "flaky" fails the time gate only while another case
    runs, "slow" fails it alone too, "memory" fails a different gate beside
    others. Every other case passes and keeps the machine busy."""
    import threading
    import time as time_module

    active = set()
    lock = threading.Lock()

    def runner(case, options):
        with lock:
            active.add(case["id"])
        time_module.sleep(0.4)
        with lock:
            others = len(active - {case["id"]})
            active.discard(case["id"])
        crowded = others > 0
        gates = []
        ratio = 0.6
        if case["id"] == "flaky" and crowded:
            gates, ratio = ["time"], 1.3
        elif case["id"] == "slow":
            gates, ratio = ["time"], 1.2
        elif case["id"] == "memory" and crowded:
            gates = ["memory"]
        return {"id": case["id"], "tier": 0, "passed": not gates, "failed_gates": gates,
                "time_gate": {"ratio": ratio}, "py_cpu_seconds": 1.0,
                "cpp_cpu_seconds": ratio, "cache_mode_seen": options.cache}
    runner.reruns_time_gate = True

    ids = ["flaky", "slow", "memory"] + [f"filler{i}" for i in range(6)]
    cases = [{"id": i, "tool": "hicInfo", "tier": 0, "args": ["x"], "outputs": [],
              "large": False, "memory": {}} for i in ids]
    options = argparse.Namespace(jobs="4", stop_after=None, cache="use", noise_runs=5,
                                 determinism=False, data=str(tmp_path), keep_workdirs=False,
                                 cache_dir=str(tmp_path / "cache"), expect_from=None)
    persisted = []
    results = {r["id"]: r for r in equiv._execute(cases, options, runner, on_result=persisted.append)}
    assert len(persisted) == len(ids)

    flaky = results["flaky"]
    assert flaky["passed"] is True
    record = flaky["rerun_alone"]
    assert record["counted"] == "rerun"
    assert record["first_attempt"]["passed"] is False
    assert record["first_attempt"]["time_gate_ratio"] == 1.3
    assert record["first_attempt"]["neighbours"] > 0
    assert record["rerun"]["passed"] is True and record["rerun"]["neighbours"] == 0
    assert flaky["scheduled"]["rerun"] is True and flaky["scheduled"]["alone"] is True
    assert flaky["cache_mode_seen"] == "refresh"

    slow = results["slow"]
    assert slow["passed"] is False
    assert slow["rerun_alone"]["rerun"]["passed"] is False
    assert slow["rerun_alone"]["rerun"]["neighbours"] == 0

    assert results["memory"]["passed"] is False and "rerun_alone" not in results["memory"]
    for i in ids[3:]:
        assert results[i]["passed"] is True and "rerun_alone" not in results[i]
    report = equiv._report_skeleton(argparse.Namespace(cpp_bin="b", py_python="p", cache="off",
                                                       jobs_resolved=4, deferred=[]),
                                    list(results.values()), "run")
    assert {entry["id"] for entry in report["reruns_alone"]} == {"flaky", "slow"}


def test_stop_after_starts_unmeasured_cases_while_time_remains(tmp_path):
    import time as time_module

    def runner(case, options):
        time_module.sleep(0.05)
        return {"id": case["id"], "tier": 0, "passed": True, "failed_gates": []}
    runner.reruns_time_gate = True
    cases = [{"id": f"large{i}", "tool": "hicInfo", "tier": 0, "args": ["x"], "outputs": [],
              "large": True, "memory": {}} for i in range(6)]
    options = argparse.Namespace(jobs="1", stop_after=30.0, cache="off", noise_runs=5,
                                 determinism=False, data=str(tmp_path), keep_workdirs=False,
                                 cache_dir=str(tmp_path / "cache"), expect_from=None)
    results = equiv._execute(cases, options, runner)
    assert len(results) == 6 and options.deferred == []
    options.stop_after = 0.0
    results = equiv._execute(cases, options, runner)
    assert len(results) == 0 and len(options.deferred) == 6


def test_reports_seed_durations_but_a_new_cache_stays_conservative():
    options = argparse.Namespace(noise_runs=5, determinism=False)
    drawing = {"id": "d", "tool": "hicPlotMatrix", "args": ["x"], "cpp_args": [], "outputs": [],
               "large": False, "memory": {}}
    large = {"id": "l", "tool": "hicInfo", "args": ["x"], "cpp_args": [], "outputs": [],
             "large": True, "memory": {}}
    reported = {"py_peak_rss_kb": 100_000, "py_seconds": 10.0, "py_cpu_seconds": 10.0,
                "cpp_peak_rss_kb": 50_000, "cpp_seconds": 2.0, "cpp_cpu_seconds": 2.0,
                "time_ratio": 0.2}
    # from a report only: measured durations, but no history in this cache
    peak, slots, seconds, measured = equiv._demand(drawing, options, reported, 8, False,
                                                   has_history=False)
    assert slots == 8 and measured and seconds == 10.0 * 1 + 2.0
    assert equiv._demand(drawing, options, reported, 8, False, has_history=True)[1] == 1
    # a large case without any measurement runs alone
    assert equiv._demand(large, options, None, 8, False, has_history=False)[1] == 8
    assert equiv._demand(large, options, reported, 8, False, has_history=False)[1] == 1
