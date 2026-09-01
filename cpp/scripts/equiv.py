#!/usr/bin/env python3
"""Python versus C++ equivalence harness for HiCExplorer v4.

Runs the same case against the 3.7.6 Python tool and against the C++ port,
compares every declared output with the comparator for its format at the
equivalence class the case declares, and records runtime and peak RSS for both.

    cpp/scripts/equiv.py run     [--tool NAME]... [--tier N]... [--case ID]...
                                 [--cpp-bin DIR] [--py-python PATH] [--jobs N]
                                 [--out DIR] [--keep-workdirs]
    cpp/scripts/equiv.py compare --format FMT --class E0 A B
    cpp/scripts/equiv.py report  [--out DIR] [--format {md,json}]
    cpp/scripts/equiv.py list    [--tool NAME]

Cases live in cpp/scripts/cases/<tool>.json. See cpp/PLAN.md section 9 for the
specification this implements. A case normally names a `tool` and the runner
executes bin/<tool> against <cpp-bin>/<tool>; a tier 0 case that exercises the
file layer before any tool exists gives `py_script` (a path relative to the
repository root) and `cpp_binary` (a path relative to --cpp-bin) instead.

Run it with the reference venv interpreter named in cpp/AGENTS_CONTRACT.md, so
that the cool and h5 comparators find h5py, PyTables and numpy:

    $VENV/bin/python cpp/scripts/equiv.py run --tool hicInfo

Exit code 0 only when every selected case passed.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import comparators  # noqa: E402  pylint: disable=C0413

HARNESS_VERSION = "1.0"

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
CASES_DIR = SCRIPT_DIR / "cases"
DEFAULT_DATA = REPO_ROOT / "hicexplorer" / "test" / "test_data"
DEFAULT_CPP_BIN = REPO_ROOT / "cpp" / "build" / "tools"
DEFAULT_OUT = REPO_ROOT / "cpp" / "build" / "equivalence"



def _default_py_python():
    """Interpreter that runs the Python reference implementation.

    HICX_PY_PYTHON wins; otherwise the throwaway venv named in
    cpp/AGENTS_CONTRACT.md is used if it still exists, and finally the
    interpreter running the harness.
    """
    from_env = os.environ.get("HICX_PY_PYTHON")
    if from_env:
        return Path(from_env)
    contract_venv = Path(
        "/tmp/claude-1000167729/-home-mh-hannover-local-wolffjoa/"
        "3de46a4d-41bb-4fb9-aea9-26cceef9ad01/scratchpad/hicx-venv/bin/python")
    if contract_venv.exists():
        return contract_venv
    return Path(sys.executable)


DEFAULT_PY_PYTHON = _default_py_python()

# The Python entry points live in bin/, one script per tool.
PY_BIN = REPO_ROOT / "bin"


# --------------------------------------------------------------------------
# case loading


def load_cases(tools=None, tiers=None, ids=None):
    cases = []
    for path in sorted(CASES_DIR.glob("*.json")):
        with open(path, encoding="utf-8") as handle:
            payload = json.load(handle)
        for case in payload:
            case.setdefault("tool", path.stem)
            case.setdefault("tier", 0)
            case.setdefault("expect_exit", 0)
            case.setdefault("outputs", [])
            case.setdefault("large", False)
            case.setdefault("notes", "")
            case.setdefault("threads_arg", None)
            # A tier 0 case names its two programs instead of a tool.
            case.setdefault("py_script", None)
            case.setdefault("cpp_binary", None)
            cases.append(case)
    if tools:
        cases = [case for case in cases if case["tool"] in tools]
    if tiers:
        cases = [case for case in cases if case["tier"] in tiers]
    if ids:
        cases = [case for case in cases if case["id"] in ids]
    return cases


def expand(value, mapping):
    for key, replacement in mapping.items():
        value = value.replace("{" + key + "}", str(replacement))
    return value


# --------------------------------------------------------------------------
# measured execution


class Measurement(dict):
    """seconds, peak_rss_kb, user_seconds, sys_seconds, exit_code."""


def run_measured(argv, cwd, stdout_path, stderr_path, env=None):
    """Run a command under /usr/bin/time and capture stdout and stderr."""
    time_file = stdout_path.parent / (stdout_path.name + ".time")
    wrapper = ["/usr/bin/time", "-f", "%e %M %U %S", "-o", str(time_file)]
    if not Path(wrapper[0]).exists():
        wrapper = []
    start = time.perf_counter()
    with open(stdout_path, "wb") as out, open(stderr_path, "wb") as err:
        completed = subprocess.run(wrapper + list(argv), cwd=str(cwd), stdout=out,
                                   stderr=err, env=env, check=False)
    wall = time.perf_counter() - start

    seconds, peak_rss, user, system = wall, 0, 0.0, 0.0
    if wrapper and time_file.exists():
        text = time_file.read_text().strip().splitlines()
        if text:
            parts = text[-1].split()
            if len(parts) == 4:
                try:
                    seconds, peak_rss, user, system = (float(parts[0]),
                                                       int(parts[1]),
                                                       float(parts[2]),
                                                       float(parts[3]))
                except ValueError:
                    pass
        time_file.unlink()
    return Measurement(seconds=seconds, peak_rss_kb=peak_rss,
                       user_seconds=user, sys_seconds=system,
                       exit_code=completed.returncode,
                       command=" ".join(shlex.quote(part) for part in argv))


# --------------------------------------------------------------------------
# a single case


def run_case(case, options):
    workdir = Path(tempfile.mkdtemp(prefix=f"equiv-{case['id']}-",
                                    dir=options.tmpdir))
    out_py = workdir / "out_py"
    out_cpp = workdir / "out_cpp"
    out_py.mkdir()
    out_cpp.mkdir()

    data = str(options.data)
    args_py = [expand(arg, {"data": data, "out": out_py}) for arg in case["args"]]
    args_cpp = [expand(arg, {"data": data, "out": out_cpp}) for arg in case["args"]]

    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    env["COLUMNS"] = "80"

    # A tool case runs bin/<tool> against <cpp-bin>/<tool>. A tier 0 case has
    # no tool and names the two programs itself; that is how the file layer is
    # exercised before any tool exists (PLAN.md 8.2, tier 0).
    python_tool = (REPO_ROOT / case["py_script"]) if case.get("py_script") \
        else (PY_BIN / case["tool"])
    cpp_tool = (Path(options.cpp_bin) / case["cpp_binary"]).resolve() \
        if case.get("cpp_binary") else (Path(options.cpp_bin) / case["tool"])

    load_average = os.getloadavg()[0]

    result = {
        "id": case["id"],
        "tool": case["tool"],
        "tier": case["tier"],
        "class_declared": _declared_class(case),
        "notes": case["notes"],
        "workdir": str(workdir),
        "timing_unreliable": load_average > 2.0,
        "load_average": load_average,
    }

    if not cpp_tool.exists():
        result.update(passed=False, class_met=None, outputs=[],
                      error=f"missing C++ binary {cpp_tool}")
        return result

    measure_py = run_measured([str(options.py_python), str(python_tool)] + args_py,
                              workdir, out_py / "stdout.txt", out_py / "stderr.txt",
                              env)
    measure_cpp = run_measured([str(cpp_tool)] + args_cpp, workdir,
                               out_cpp / "stdout.txt", out_cpp / "stderr.txt")

    result.update({
        "py_seconds": measure_py["seconds"],
        "cpp_seconds": measure_cpp["seconds"],
        "py_peak_rss_kb": measure_py["peak_rss_kb"],
        "cpp_peak_rss_kb": measure_cpp["peak_rss_kb"],
        "py_user_seconds": measure_py["user_seconds"],
        "cpp_user_seconds": measure_cpp["user_seconds"],
        "py_sys_seconds": measure_py["sys_seconds"],
        "cpp_sys_seconds": measure_cpp["sys_seconds"],
        "py_exit": measure_py["exit_code"],
        "cpp_exit": measure_cpp["exit_code"],
        "py_command": measure_py["command"],
        "cpp_command": measure_cpp["command"],
        "stderr_py": _tail(out_py / "stderr.txt"),
        "stderr_cpp": _tail(out_cpp / "stderr.txt"),
    })

    passed = True
    if measure_py["exit_code"] != case["expect_exit"]:
        passed = False
        result["error"] = (f"the Python tool exited {measure_py['exit_code']}, "
                           f"expected {case['expect_exit']}")
    if measure_cpp["exit_code"] != case["expect_exit"]:
        passed = False
        result["error"] = (f"the C++ tool exited {measure_cpp['exit_code']}, "
                           f"expected {case['expect_exit']}")

    outputs = []
    for declared in case["outputs"]:
        path_py = Path(expand(declared["path"], {"data": data, "out": out_py}))
        path_cpp = Path(expand(declared["path"], {"data": data, "out": out_cpp}))
        entry = {
            "path": declared["path"],
            "format": declared["format"],
            "class": declared.get("class", "E0"),
        }
        if not path_py.exists() or not path_cpp.exists():
            entry.update(passed=False, class_met=None, metrics={},
                         diffs=[f"missing output: python={path_py.exists()} "
                                f"cpp={path_cpp.exists()}"])
            outputs.append(entry)
            passed = False
            continue
        comparison = comparators.compare(declared["format"], str(path_py),
                                         str(path_cpp), entry["class"],
                                         declared.get("options"))
        entry.update(comparison.to_json())
        outputs.append(entry)
        passed = passed and comparison.passed

    result["outputs"] = outputs
    result["passed"] = passed
    result["class_met"] = result["class_declared"] if passed else None

    if not options.keep_workdirs and passed:
        shutil.rmtree(workdir, ignore_errors=True)
        result["workdir"] = None
    return result


def _declared_class(case):
    classes = {output.get("class", "E0") for output in case["outputs"]}
    return "/".join(sorted(classes)) if classes else "E0"


def _tail(path, limit=4000):
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return ""
    return text[-limit:]


# --------------------------------------------------------------------------
# reporting


def write_report(report, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "report.json"
    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=False)
        handle.write("\n")
    markdown_path = out_dir / "report.md"
    markdown_path.write_text(render_markdown(report), encoding="utf-8")
    return json_path, markdown_path


def render_markdown(report):
    lines = []
    lines.append("# Equivalence report")
    lines.append("")
    lines.append(f"- harness version: {report['harness_version']}")
    lines.append(f"- generated: {report['timestamp']}")
    lines.append(f"- git commit: {report['git_commit']}")
    lines.append(f"- host: {report['host']}")
    cases = report["cases"]
    passed = sum(1 for case in cases if case.get("passed"))
    lines.append(f"- cases: {len(cases)}, passed: {passed}, "
                 f"failed: {len(cases) - passed}")
    lines.append("")

    tiers = sorted({case["tier"] for case in cases})
    for tier in tiers:
        lines.append(f"## Tier {tier}")
        lines.append("")
        lines.append("| case | class | result | py s | cpp s | speedup | "
                     "py RSS MB | cpp RSS MB | mem ratio |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        for case in [case for case in cases if case["tier"] == tier]:
            speedup = _ratio(case.get("py_seconds"), case.get("cpp_seconds"))
            memory = _ratio(case.get("cpp_peak_rss_kb"), case.get("py_peak_rss_kb"))
            lines.append(
                "| {id} | {cls} | {result} | {py:.3f} | {cpp:.3f} | {speedup} | "
                "{pyrss:.1f} | {cpprss:.1f} | {memory} |".format(
                    id=case["id"],
                    cls=case.get("class_declared", "?"),
                    result="pass" if case.get("passed") else "FAIL",
                    py=case.get("py_seconds", 0.0),
                    cpp=case.get("cpp_seconds", 0.0),
                    speedup=speedup,
                    pyrss=case.get("py_peak_rss_kb", 0) / 1024.0,
                    cpprss=case.get("cpp_peak_rss_kb", 0) / 1024.0,
                    memory=memory))
        lines.append("")

    failures = [case for case in cases if not case.get("passed")]
    if failures:
        lines.append("## Failures")
        lines.append("")
        for case in failures:
            lines.append(f"### {case['id']}")
            lines.append("")
            if case.get("error"):
                lines.append(f"- {case['error']}")
            lines.append(f"- python: `{case.get('py_command', '')}`")
            lines.append(f"- c++:    `{case.get('cpp_command', '')}`")
            for output in case.get("outputs", []):
                if output.get("passed"):
                    continue
                lines.append(f"- output `{output['path']}` "
                             f"({output['format']}, {output['class']}):")
                for diff in output.get("diffs", [])[:20]:
                    lines.append(f"      {diff}")
            if case.get("stderr_cpp"):
                lines.append("- C++ stderr:")
                for line in case["stderr_cpp"].strip().splitlines()[-10:]:
                    lines.append(f"      {line}")
            lines.append("")
    return "\n".join(lines) + "\n"


def _ratio(numerator, denominator):
    if not numerator or not denominator:
        return "n/a"
    return f"{numerator / denominator:.2f}"


def git_commit():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(REPO_ROOT),
                              capture_output=True, text=True,
                              check=False).stdout.strip() or "unknown"
    except OSError:
        return "unknown"


# --------------------------------------------------------------------------
# entry points


def command_run(options):
    cases = load_cases(options.tool, options.tier, options.case)
    if not cases:
        print("no cases selected", file=sys.stderr)
        return 1
    results = []
    if options.jobs > 1:
        with concurrent.futures.ThreadPoolExecutor(max_workers=options.jobs) as pool:
            futures = {pool.submit(run_case, case, options): case for case in cases}
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())
    else:
        for case in cases:
            results.append(run_case(case, options))
    results.sort(key=lambda item: (item["tier"], item["id"]))

    report = {
        "harness_version": HARNESS_VERSION,
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": git_commit(),
        "host": platform.node(),
        "cpp_bin": str(options.cpp_bin),
        "py_python": str(options.py_python),
        "cases": results,
    }
    json_path, markdown_path = write_report(report, options.out)

    for case in results:
        status = "pass" if case.get("passed") else "FAIL"
        speedup = _ratio(case.get("py_seconds"), case.get("cpp_seconds"))
        print(f"{status:4} {case['id']:<48} "
              f"py {case.get('py_seconds', 0):7.3f}s  "
              f"cpp {case.get('cpp_seconds', 0):7.3f}s  "
              f"speedup {speedup:>6}  "
              f"rss {case.get('py_peak_rss_kb', 0) / 1024:8.1f} / "
              f"{case.get('cpp_peak_rss_kb', 0) / 1024:8.1f} MB")
        if not case.get("passed"):
            if case.get("error"):
                print(f"       {case['error']}")
            for output in case.get("outputs", []):
                for diff in output.get("diffs", [])[:5]:
                    print(f"       {diff}")
    print(f"\nreport: {json_path}\n        {markdown_path}")
    return 0 if all(case.get("passed") for case in results) else 1


def command_compare(options):
    result = comparators.compare(options.format, options.a, options.b,
                                 getattr(options, "class"), None)
    print(json.dumps(result.to_json(), indent=2))
    return 0 if result.passed else 1


def command_report(options):
    json_path = Path(options.out) / "report.json"
    if not json_path.exists():
        print(f"no report at {json_path}", file=sys.stderr)
        return 1
    with open(json_path, encoding="utf-8") as handle:
        report = json.load(handle)
    if options.format == "json":
        print(json.dumps(report, indent=2))
    else:
        text = render_markdown(report)
        (Path(options.out) / "report.md").write_text(text, encoding="utf-8")
        print(text)
    return 0


def command_list(options):
    for case in load_cases(options.tool, None, None):
        print(f"{case['tier']}  {case['id']:<48} {case['tool']:<16} "
              f"{_declared_class(case)}  {case['notes']}")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="run cases and compare")
    run_parser.add_argument("--tool", action="append")
    run_parser.add_argument("--tier", action="append", type=int)
    run_parser.add_argument("--case", action="append")
    run_parser.add_argument("--cpp-bin", default=str(DEFAULT_CPP_BIN))
    run_parser.add_argument("--py-python", default=str(DEFAULT_PY_PYTHON))
    run_parser.add_argument("--data", default=str(DEFAULT_DATA))
    run_parser.add_argument("--jobs", type=int, default=1)
    run_parser.add_argument("--out", default=str(DEFAULT_OUT))
    run_parser.add_argument("--tmpdir", default=tempfile.gettempdir())
    run_parser.add_argument("--keep-workdirs", action="store_true")
    run_parser.set_defaults(handler=command_run)

    compare_parser = subparsers.add_parser("compare", help="compare two files")
    compare_parser.add_argument("--format", required=True,
                                choices=comparators.formats())
    compare_parser.add_argument("--class", required=True,
                                choices=list(comparators.CLASSES))
    compare_parser.add_argument("a")
    compare_parser.add_argument("b")
    compare_parser.set_defaults(handler=command_compare)

    report_parser = subparsers.add_parser("report", help="rerender a report")
    report_parser.add_argument("--out", default=str(DEFAULT_OUT))
    report_parser.add_argument("--format", choices=("md", "json"), default="md")
    report_parser.set_defaults(handler=command_report)

    list_parser = subparsers.add_parser("list", help="list the known cases")
    list_parser.add_argument("--tool", action="append")
    list_parser.set_defaults(handler=command_list)

    options = parser.parse_args(argv)
    return options.handler(options)


if __name__ == "__main__":
    sys.exit(main())
