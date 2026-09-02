#!/usr/bin/env python3
"""Numeric validation for the hicCorrectMatrix cases the harness cannot compare.

Why this exists. cpp/scripts/comparators/cool.py compares every column of a
cooler's /bins group with an exact array equality, at every equivalence class:

    for column in columns_a:
        if not _arrays_identical(group_a["bins/" + column][:], ...):
            diffs.append(...)

That is right for chrom, start and end, which are integers, and it was right
for every tool ported so far, none of which writes a float bin column.
hicCorrectMatrix is the first that does: the correction factors leave the tool
through /bins/weight. A balancing result is never bit identical to another
implementation's, so under the current comparator no cool output of this tool
can pass at any class, including at ED. Making _compare_bins class aware is a
one function change, but comparators/ is owned by the harness and not by this
tool, so the cases keep their gates and the numeric comparison is done here.

This script also implements class EN of cpp/PLAN.md 5.7, which the harness does
not have yet: the tolerance is measured from N runs of the reference rather than
chosen.

    hicCorrectMatrix_validate.py [--runs N] [--case ID]... [--keep]

For every configuration it prints the class the C++ output met, the measured
oracle envelope where the reference is nondeterministic, and the peak RSS and
CPU time of both implementations.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
DATA = REPO_ROOT / "hicexplorer" / "test" / "test_data"
CPP_BIN = REPO_ROOT / "cpp" / "build" / "tools" / "hicCorrectMatrix"
PY_BIN = REPO_ROOT / "bin" / "hicCorrectMatrix"
VENV = Path("/tmp/claude-1000167729/-home-mh-hannover-local-wolffjoa/"
            "3de46a4d-41bb-4fb9-aea9-26cceef9ad01/scratchpad/hicx-venv/bin/python")

# The acceptance gate and the stricter classes, cpp/PLAN.md 5.1.
CLASS_TOLERANCE = {"E2": (0.0, 1.0), "E3": (1e-12, 1.0), "E4": (1e-6, 1e-6),
                   "ED": (1e-3, 0.0)}
CLASS_ORDER = ["E2", "E3", "E4", "ED"]

CASES = [
    {"id": "KR.gm12878_raw_values.cool", "class": "EN", "output": "out.cool",
     "args": ["correct", "-m", str(DATA / "hicCorrectMatrix/gm12878_raw_values.cool"),
              "--correctionMethod", "KR", "-o", "{out}"]},
    {"id": "KR.gm12878_raw_values.cool.v3", "class": "EN", "output": "out.cool",
     "cpp_extra": ["--compatMode", "v3"],
     "args": ["correct", "-m", str(DATA / "hicCorrectMatrix/gm12878_raw_values.cool"),
              "--correctionMethod", "KR", "-o", "{out}"]},
    {"id": "KR.Li_et_al_2015.h5.v3", "class": "EN", "output": "out.h5",
     "cpp_extra": ["--compatMode", "v3"], "reader": "h5",
     "args": ["correct", "-m", str(DATA / "Li_et_al_2015.h5"),
              "--correctionMethod", "KR", "-o", "{out}"]},
    # The configuration of the existing test_correct_matrix_KR_partial_cool,
    # whose assert_allclose(rtol=1.0) cannot fail. A single --chromosomes on a
    # cool input takes hicCorrectMatrix.py:603's fast path and loads only that
    # chromosome's block, which gives a different NaN bin list from a whole file
    # load followed by a selection.
    {"id": "KR.one_chromosome.gm12878_raw_values.cool", "class": "EN", "output": "out.cool",
     "args": ["correct", "-m", str(DATA / "hicCorrectMatrix/gm12878_raw_values.cool"),
              "--correctionMethod", "KR", "--chromosomes", "3", "-o", "{out}"]},
    {"id": "KR.perchr.gm12878_raw_values.cool", "class": "EN", "output": "out.cool",
     "args": ["correct", "-m", str(DATA / "hicCorrectMatrix/gm12878_raw_values.cool"),
              "--correctionMethod", "KR", "--perchr", "-o", "{out}"]},
    {"id": "KR.perchr.h5name.gm12878_raw_values.cool", "class": "EN", "output": "out.h5",
     "args": ["correct", "-m", str(DATA / "hicCorrectMatrix/gm12878_raw_values.cool"),
              "--correctionMethod", "KR", "--perchr", "-o", "{out}"]},
    {"id": "ICE.gm12878_raw_values.cool", "class": "E3", "output": "out.cool",
     "args": ["correct", "-m", str(DATA / "hicCorrectMatrix/gm12878_raw_values.cool"),
              "--correctionMethod", "ICE", "--filterThreshold", "-1.5", "5.0",
              "-o", "{out}"]},
    {"id": "ICE.perchr.gm12878_raw_values.cool", "class": "E3", "output": "out.cool",
     "args": ["correct", "-m", str(DATA / "hicCorrectMatrix/gm12878_raw_values.cool"),
              "--correctionMethod", "ICE", "--filterThreshold", "-1.5", "5.0",
              "--perchr", "-o", "{out}"]},
    {"id": "ICE.gm12878_chr1.cool", "class": "E3", "output": "out.cool", "large": True,
     "args": ["correct", "-m", str(DATA / "hicTADClassifier/gm12878_chr1.cool"),
              "--correctionMethod", "ICE", "--filterThreshold", "-1.5", "5",
              "-o", "{out}"]},
    # The float64 default cannot lie inside the oracle's envelope on this input,
    # because at 61.8 M terms the reference's float32 accumulator is biased by
    # more than its own spread. The EN comparison is therefore run in v3.
    {"id": "KR.gm12878_chr1.cool.v3", "class": "EN", "output": "out.cool", "large": True,
     "cpp_extra": ["--compatMode", "v3"],
     "args": ["correct", "-m", str(DATA / "hicTADClassifier/gm12878_chr1.cool"),
              "--correctionMethod", "KR", "-o", "{out}"]},
    {"id": "KR.gm12878_chr1.cool.v4", "class": "EN", "output": "out.cool", "large": True,
     "args": ["correct", "-m", str(DATA / "hicTADClassifier/gm12878_chr1.cool"),
              "--correctionMethod", "KR", "-o", "{out}"]},
]


def read_h5(path):
    """The numeric arrays of a HiCExplorer h5 matrix."""
    import tables

    arrays = {}
    with tables.open_file(path) as handle:
        arrays["matrix/data"] = handle.root.matrix.data.read()
        arrays["matrix/indices"] = handle.root.matrix.indices.read()
        arrays["matrix/indptr"] = handle.root.matrix.indptr.read()
        for name in ("correction_factors", "distance_counts"):
            if name in handle.root:
                arrays[name] = np.asarray(getattr(handle.root, name).read()).ravel()
    return arrays


def read_matrix(path, reader):
    return read_h5(path) if reader == "h5" else read_cool(path)


def read_cool(path):
    """Every numeric array of a cooler, as a flat mapping."""
    import h5py

    arrays = {}
    with h5py.File(path, "r") as handle:
        for name in ("pixels/bin1_id", "pixels/bin2_id", "pixels/count",
                     "bins/start", "bins/end", "chroms/length"):
            arrays[name] = handle[name][:]
        if "weight" in handle["bins"]:
            arrays["bins/weight"] = handle["bins/weight"][:]
        arrays["_attrs"] = {key: handle.attrs[key] for key in handle.attrs
                            if key not in ("creation-date", "generated-by",
                                           "generated-by-cooler-lib", "tool-url")}
    return arrays


def relative_difference(a, b):
    """max |a - b| / |b| over the entries, with an exact zero matched exactly."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return float("inf")
    finite = np.isfinite(a) & np.isfinite(b)
    if not np.array_equal(np.isfinite(a), np.isfinite(b)):
        return float("inf")
    zero = finite & (b == 0.0)
    if np.any(a[zero] != 0.0):
        return float("inf")
    live = finite & (b != 0.0)
    if not live.any():
        return 0.0
    return float(np.max(np.abs(a[live] - b[live]) / np.abs(b[live])))


def compare(reference, candidate):
    """Per array relative difference, and the coordinates that must match."""
    report = {}
    for name in reference:
        if name == "_attrs":
            continue
        if name in ("pixels/bin1_id", "pixels/bin2_id",
                    "matrix/indices", "matrix/indptr"):
            same = (reference[name].shape == candidate[name].shape
                    and np.array_equal(reference[name], candidate[name]))
            report[name] = 0.0 if same else float("inf")
            continue
        report[name] = relative_difference(candidate[name], reference[name])
    return report


def class_met(report):
    worst = max(report.values()) if report else 0.0
    for name in CLASS_ORDER:
        tolerance, _ = CLASS_TOLERANCE[name]
        if worst <= tolerance:
            return name, worst
    return None, worst


def run(argv, out_path, label):
    """Runs a command under /usr/bin/time and returns (seconds cpu, peak rss kb)."""
    marker = tempfile.NamedTemporaryFile(delete=False)
    marker.close()
    command = ["/usr/bin/time", "-f", "%U %S %M", "-o", marker.name] + argv
    environment = dict(os.environ, PYTHONPATH=str(REPO_ROOT))
    started = time.perf_counter()
    completed = subprocess.run(command, capture_output=True, text=True,
                               env=environment, cwd=str(REPO_ROOT), check=False)
    wall = time.perf_counter() - started
    with open(marker.name, encoding="utf-8") as handle:
        fields = handle.read().split()
    os.unlink(marker.name)
    user, system, rss = (float(fields[0]), float(fields[1]), float(fields[2]))
    if completed.returncode != 0:
        print(f"    {label} exited {completed.returncode}: "
              f"{completed.stderr.strip().splitlines()[-1:]}")
    return {"cpu": user + system, "wall": wall, "rss_kb": rss,
            "exit": completed.returncode, "stderr": completed.stderr}


def median_run(outputs):
    """The run whose weight column is closest to the elementwise median."""
    key = "bins/weight" if "bins/weight" in outputs[0] else "correction_factors"
    if key not in outputs[0]:
        return 0
    stack = np.vstack([entry[key] for entry in outputs])
    reference = np.median(stack, axis=0)
    distances = [float(np.nanmax(np.abs(entry - reference))) for entry in stack]
    return int(np.argmin(distances))


def envelope(outputs):
    """S: the maximum pairwise relative difference among the reference runs."""
    worst = {}
    for i in range(len(outputs)):
        for j in range(i + 1, len(outputs)):
            for name, value in compare(outputs[i], outputs[j]).items():
                worst[name] = max(worst.get(name, 0.0), value)
    return worst


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=int, default=5)
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--keep", action="store_true")
    parser.add_argument("--tmpdir", default="/tmp")
    options = parser.parse_args()

    report = []
    failures = 0
    for case in CASES:
        if options.case and case["id"] not in options.case:
            continue
        runs = 1 if case["class"] != "EN" else options.runs
        directory = Path(tempfile.mkdtemp(prefix="hcm-validate-", dir=options.tmpdir))
        print(f"=== {case['id']}  class {case['class']}  reference runs {runs}",
              flush=True)
        try:
            reference_outputs = []
            reference_cost = []
            for index in range(runs):
                path = directory / f"py{index}" / case["output"]
                path.parent.mkdir(parents=True, exist_ok=True)
                argv = [str(VENV), str(PY_BIN)] + [
                    argument.replace("{out}", str(path)) for argument in case["args"]]
                cost = run(argv, path, "python")
                reference_cost.append(cost)
                reference_outputs.append(read_matrix(path, case.get("reader")))

            cpp_path = directory / "cpp" / case["output"]
            cpp_path.parent.mkdir(parents=True, exist_ok=True)
            argv = [str(CPP_BIN)] + [argument.replace("{out}", str(cpp_path))
                                     for argument in case["args"]]
            argv += case.get("cpp_extra", [])
            cpp_cost = run(argv, cpp_path, "c++")
            cpp_output = read_matrix(cpp_path, case.get("reader"))

            index = median_run(reference_outputs)
            report_entry = {
                "id": case["id"], "class_declared": case["class"],
                "reference_runs": runs, "median_run": index,
                "py_cpu": statistics.median(c["cpu"] for c in reference_cost),
                "py_rss_mb": max(c["rss_kb"] for c in reference_cost) * 1024 / 1e6,
                "cpp_cpu": cpp_cost["cpu"],
                "cpp_rss_mb": cpp_cost["rss_kb"] * 1024 / 1e6,
            }
            difference = compare(reference_outputs[index], cpp_output)
            report_entry["difference"] = difference

            if case["class"] == "EN":
                spread = envelope(reference_outputs)
                report_entry["oracle_envelope"] = spread
                passed = True
                for name, value in difference.items():
                    limit = max(2.0 * spread.get(name, 0.0), 1e-12)
                    if value > limit:
                        passed = False
                        print(f"    {name}: {value:.3e} outside "
                              f"max(2S, 1e-12) = {limit:.3e}")
                report_entry["passed"] = passed
                report_entry["class_met"] = "EN" if passed else None
                print(f"    oracle envelope S: "
                      + ", ".join(f"{k} {v:.3e}" for k, v in sorted(spread.items())
                                  if v > 0) or "    oracle envelope S: 0 (reproducible)")
            else:
                met, worst = class_met(difference)
                report_entry["class_met"] = met
                report_entry["worst_relative"] = worst
                report_entry["passed"] = met is not None
                print(f"    worst relative difference {worst:.3e}, class met {met}")

            print("    per array: " + ", ".join(
                f"{k} {v:.3e}" for k, v in sorted(difference.items())))
            print(f"    cpu  python {report_entry['py_cpu']:.2f} s   "
                  f"c++ {report_entry['cpp_cpu']:.2f} s   "
                  f"speedup {report_entry['py_cpu'] / max(report_entry['cpp_cpu'], 1e-9):.1f}x")
            print(f"    rss  python {report_entry['py_rss_mb']:.0f} MB  "
                  f"c++ {report_entry['cpp_rss_mb']:.0f} MB  "
                  f"reduction {report_entry['py_rss_mb'] / max(report_entry['cpp_rss_mb'], 1e-9):.1f}x")
            if not report_entry["passed"]:
                failures += 1
            print(f"    {'PASS' if report_entry['passed'] else 'FAIL'}", flush=True)
            report.append(report_entry)
        finally:
            if not options.keep:
                shutil.rmtree(directory, ignore_errors=True)
            else:
                print(f"    kept {directory}")

    output = REPO_ROOT / "cpp" / "build" / "equivalence" / "hicCorrectMatrix_validate.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=1, default=float)
    print(f"\nreport: {output}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
