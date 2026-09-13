#!/usr/bin/env python3
"""Python versus C++ equivalence harness for HiCExplorer v4.

Runs the same case against the 3.7.6 Python tool and against the C++ port,
compares every declared output with the comparator for its format at the
equivalence class the case declares, and then applies three resource gates:
memory, CPU time and determinism.

    cpp/scripts/equiv.py run          [--tool NAME]... [--tier N]... [--case ID]...
                                      [--cpp-bin DIR] [--py-python PATH] [--jobs N]
                                      [--out DIR] [--keep-workdirs]
                                      [--skip-memory-gate] [--skip-time-gate]
                                      [--determinism] [--noise-runs N]
    cpp/scripts/equiv.py determinism  [--tool NAME]... [--case ID]... [--noise-runs N]
                                      [--threads-high N] [--out DIR]
    cpp/scripts/equiv.py compare      --format FMT --class E0 A B
    cpp/scripts/equiv.py report       [--out DIR] [--format {md,json}]
    cpp/scripts/equiv.py list         [--tool NAME]

Cases live in cpp/scripts/cases/<tool>.json. See cpp/PLAN.md section 9 for the
specification this implements. A case normally names a `tool` and the runner
executes bin/<tool> against <cpp-bin>/<tool>; a tier 0 case that exercises the
file layer before any tool exists gives `py_script` (a path relative to the
repository root) and `cpp_binary` (a path relative to --cpp-bin) instead.
A case may also give `cpp_args`, appended to the C++ command line only, for an
option the port has and the Python does not (hicCompartmentalization --noPlot);
its notes must say why the two command lines differ.

Run it with the reference venv interpreter named in cpp/AGENTS_CONTRACT.md, so
that the cool and h5 comparators find h5py, PyTables and numpy:

    $VENV/bin/python cpp/scripts/equiv.py run --tool hicInfo

Exit code 0 only when every selected case passed every comparator and every
gate.

Units
-----
All memory figures are SI: 1 MB is 1e6 bytes. `/usr/bin/time -f %M` reports KiB
(it prints `getrusage`'s `ru_maxrss`), so the harness multiplies it by 1024 to
get bytes and divides by 1e6 to print MB. The budget formula of PLAN.md 4.5 is
in SI bytes as well (`W = 741.9 MB` for 741,856,800 bytes), so the two sides of
the comparison use one convention. Mixing MiB against an SI budget would
overstate the headroom by 4.9 %, which on the tightest case in the corpus is
the whole margin.

The memory gate (PLAN.md 4.5 and 8.3 criterion 4)
-------------------------------------------------
    W      = nnz_stored * (8 + 4) + (nbins + 1) * 8   stored upper-triangle CSR
    D      = max_chromosome_bins^2 * 8                largest dense chrom block
    C      = 64 MB
    budget = alpha * W + beta * D + C

A case declares `alpha` and `beta`; `W` and `D` are measured from the actual
input, never taken from the case, so a budget cannot drift away from the data
it is supposed to describe. Inputs are the arguments that resolve to an
existing cool, mcool or h5 matrix; `W` and `D` are the maxima over them, since
the budget's alpha already accounts for how many of them are live at once. When
a case has no matrix input (text or hicpro input) the matrix the Python run
*produced* is measured instead, and when there is neither, the gate is recorded
as not applicable rather than silently passing on `budget = C`. A case that
also declares `nnz_stored`/`nbins` gets those cross-checked against the
measurement and a mismatch is reported, but the measurement wins.

The C++ peak RSS exceeding the budget fails the case exactly as a comparator
mismatch does. `alpha`/`beta` come from the case where it declares them and
otherwise from ALPHA_BETA_BY_TOOL, which is PLAN.md 4.5's table; the report
records which. Raising a budget is a reviewed edit to the case file and
STATUS.md, never an inline exemption (PLAN.md 9.5).

The time gate
-------------
Design decision, 2026-09-01, in answer to the project owner's request that time
be measured alongside memory:

1. **The gate is on CPU time (user + sys), not on wall clock.** This machine
   runs unrelated training jobs; the one-minute load average has been between
   13 and 57 throughout the project and was 8 to 12 during the run that
   calibrated this code. Wall clock under that load is not reproducible and
   would make the gate a lottery. CPU time is charged to the process whatever
   else the machine is doing.
2. **The gate is `cpp_cpu <= py_cpu`.** A rewrite that burns more CPU than the
   interpreted original is a defect, whatever the wall clock says. The ratio
   `cpp_cpu / py_cpu` is recorded; so is the speedup `py_cpu / cpp_cpu`.
3. **Wall clock and the load average are recorded but not gated**, marked
   indicative. `timing_unreliable` stays and is driven by the load average
   exceeding LOAD_AVERAGE_LIMIT, as PLAN.md 10 specifies.
4. **A floor below which the gate is skipped.** `/usr/bin/time` resolves CPU
   time to a 10 ms clock tick. Measured on this machine, 15 repeats each: the
   C++ `hicInfo --version` and the cool metadata path report 0.00 s every time,
   the h5 load path scatters over 0.02, 0.03, 0.04 and 0.05 s (a 150 % spread
   on a 0.02 s reading), and the 0.07 s cool load path is stable to the tick in
   all 15. Below roughly five ticks the number is quantisation, not work.
   TIME_FLOOR_SECONDS is therefore 0.05 s, and:
     - `py_cpu < floor`  -> the gate is skipped and the case is marked
       `time_gate.applied = false`, because the reference itself is inside the
       noise and the comparison would be meaningless in either direction. In
       practice this never fires: importing hicexplorer costs 1.9 to 2.1 s of
       CPU (0.4 s of wall time across several threads) before a tool does any
       work, so every Python case is far above the floor.
     - `cpp_cpu < floor` -> the verdict still stands, because `cpp < py` then
       holds by a margin far larger than the measurement error, but the ratio
       is meaningless and is reported as a lower bound with
       `startup_dominated = true` instead of a spurious speedup figure.

Determinism (PLAN.md 4.1 and 8.3 criterion 3)
---------------------------------------------
`equiv.py determinism`, and `run --determinism`, run the C++ tool
`--noise-runs` times (default 5) and require every declared output to be
byte-identical across the repeats, and, for a case that names `threads_arg`,
require the `--threads 1` and `--threads <--threads-high>` outputs to be
byte-identical too. Byte-identical means byte-identical: this check does not go
through the comparators and has no tolerance. No tool ported so far takes
`--threads`, so the thread half is a no-op today; hicCorrectMatrix and
hicFindTADs need it and it exists before they arrive.

Two things make a raw byte comparison alone a poor check for HDF5 output, and
both are handled explicitly rather than by loosening it:

1. An HDF5 object header carries an optional modification time (message 0x12),
   written unless the creator clears `H5Pset_obj_track_times`. Two runs a
   second apart therefore differ, and two runs in the same second do not, so a
   byte comparison of a fast tool passes or fails by luck. The check does not
   rely on it: it scans the output for embedded wall-clock stamps directly and
   fails the case if it finds any, which is deterministic. The PyTables and
   cooler reference files carry the same messages, so this is the HDF5 default
   rather than a regression, but it is still a reproducibility defect in the
   output and the fix is one line.
2. A cool file carries a `creation-date` attribute that cannot be identical
   between two writes, whoever writes them. PLAN.md 5.1 already rules that
   field and three other provenance fields non-significant, so a byte
   difference confined to them is recorded as a qualification rather than a
   failure. Every differing byte is classified, so an object timestamp is never
   mistaken for a numeric difference and, just as important, a numeric
   difference is never excused as a timestamp: when the content differs at the
   strictest class the format admits, the case fails.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import filecmp
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import comparators  # noqa: E402  pylint: disable=C0413

HARNESS_VERSION = "1.2"

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
CASES_DIR = SCRIPT_DIR / "cases"
DEFAULT_DATA = REPO_ROOT / "hicexplorer" / "test" / "test_data"
DEFAULT_CPP_BIN = REPO_ROOT / "cpp" / "build" / "tools"
DEFAULT_OUT = REPO_ROOT / "cpp" / "build" / "equivalence"

# --- memory budget (PLAN.md 4.5) ------------------------------------------
MB = 1_000_000.0                      # SI, as in PLAN.md and STATUS.md
BUDGET_CONSTANT_BYTES = 64 * MB       # C: process, HDF5, buffers, output staging
VALUE_BYTES = 8                       # float64 values
INDEX_BYTES = 4                       # int32 while nbins <= INT32_MAX, rule 3
INDPTR_BYTES = 8

# alpha and beta per tool, transcribed from the table in PLAN.md 4.5. This is
# the fallback for a case that does not declare its own; the case wins when it
# does. `roundtrip` is the tier 0 load-and-write exerciser, whose exit criterion
# in PLAN.md 8.2 is 1.3 * W + C.
ALPHA_BETA_BY_TOOL = {
    "hicInfo": (1.05, 0.0),
    "hicQuickQC": (1.05, 0.0),
    "hicConvertFormat": (1.3, 0.0),
    "hicNormalize": (1.3, 0.0),
    "hicAdjustMatrix": (1.3, 0.0),
    "hicMergeMatrixBins": (1.3, 0.0),
    "hicMergeTADbins": (1.3, 0.0),
    "hicAverageRegions": (1.3, 0.0),
    # Two matrices live plus the result. scipy's binary operations drop
    # exactly-zero results, so the output is neither a superset nor a subset of
    # either operand and there is no safe in-place merge: the peak is
    # nnz_a + nnz_b + nnz_result, about 2.9 W. Corrected from 2.2 on
    # 2026-09-01; at 2.2 the largest case in the corpus,
    # hicSumMatrices.cool.GSM_pair_chr1_chr2, measures 100.4 % of its budget.
    "hicSumMatrices": (3.0, 0.0),
    "hicCompareMatrices": (3.0, 0.0),
    "hicCorrectMatrix": (1.2, 0.0),
    "hicTransform": (1.2, 0.0),
    "hicFindTADs": (2.2, 0.0),
    "hicDetectLoops": (2.2, 0.0),
    "hicAggregateContacts": (2.2, 0.0),
    "roundtrip": (1.3, 0.0),
}
DEFAULT_ALPHA_BETA = (1.3, 0.0)       # "everything else" in PLAN.md 4.5

MATRIX_SUFFIXES = (".cool", ".mcool", ".scool", ".h5")

# --- time gate -------------------------------------------------------------
# Five 10 ms clock ticks. See the module docstring for the measurement.
TIME_FLOOR_SECONDS = 0.05
LOAD_AVERAGE_LIMIT = 2.0              # PLAN.md 10

# --- determinism -----------------------------------------------------------
DEFAULT_NOISE_RUNS = 5                # PLAN.md 9.1
DEFAULT_THREADS_HIGH = 16             # PLAN.md 4.1: --threads 1 versus 16

RSS_TRACE_INTERVAL = 0.1              # PLAN.md 10, for large: true cases


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
            case.setdefault("memory", {})
            # A tier 0 case names its two programs instead of a tool.
            case.setdefault("py_script", None)
            case.setdefault("cpp_binary", None)
            # Arguments appended to the C++ command line only, for an option
            # the port has and the Python does not (hicCompartmentalization
            # --noPlot). Never a way to make the two runs do different work on
            # the outputs being compared: the case notes must say why.
            case.setdefault("cpp_args", [])
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
# the working set of a matrix, measured from the file


_MATRIX_INFO_CACHE = {}


def matrix_info(path):
    """(nbins, nnz_stored, max_chromosome_bins) of a matrix file, or None.

    Accepts a cooler URI (`file.mcool::/resolutions/10000`). Without a `::`
    suffix every cooler group in the file is inspected and the largest is
    taken, which is the conservative choice for the budget. Returns None for
    anything that is not a readable matrix, so that a BED, BAM or text argument
    is simply not an input for budget purposes.
    """
    key = str(path)
    if key in _MATRIX_INFO_CACHE:
        return _MATRIX_INFO_CACHE[key]
    info = _matrix_info_uncached(key)
    _MATRIX_INFO_CACHE[key] = info
    return info


def _matrix_info_uncached(uri):
    filename, _, group = uri.partition("::")
    if not os.path.isfile(filename):
        return None
    lowered = filename.lower()
    if not lowered.endswith(MATRIX_SUFFIXES):
        return None
    try:
        if lowered.endswith(".h5"):
            info = _h5_matrix_info(filename)
        else:
            info = _cool_matrix_info(filename, group or None)
    except Exception as error:  # pylint: disable=W0718
        return {"error": f"{type(error).__name__}: {error}", "path": filename}
    if info is not None:
        info["path"] = filename
    return info


def _cool_matrix_info(filename, group):
    """Largest cooler group in a cool / mcool / scool file.

    Reads only dataset shapes and the per-bin chromosome column, so the cost is
    independent of the number of pixels.
    """
    import h5py
    import numpy as np

    required = {"chroms", "bins", "pixels", "indexes"}
    with h5py.File(filename, "r") as handle:
        if group:
            groups = [group if group.startswith("/") else "/" + group]
        else:
            groups = []
            if required <= set(handle.keys()):
                groups.append("/")

            def visit(name, obj):
                if isinstance(obj, h5py.Group) and required <= set(obj.keys()):
                    groups.append("/" + name)

            handle.visititems(visit)
        best = None
        for name in groups:
            node = handle[name]
            if not required <= set(node.keys()):
                continue
            chrom = node["bins/chrom"][:]
            counts = np.bincount(np.asarray(chrom, dtype=np.int64))
            candidate = {
                "nbins": int(chrom.shape[0]),
                "nnz_stored": int(node["pixels/bin1_id"].shape[0]),
                "max_chromosome_bins": int(counts.max()) if counts.size else 0,
                "group": name,
            }
            if best is None or candidate["nnz_stored"] > best["nnz_stored"]:
                best = candidate
        return best


def _h5_matrix_info(filename):
    """The HiCExplorer PyTables layout.

    Read with PyTables rather than h5py, because the blosc filter these files
    use is a PyTables filter that plain h5py cannot decode here. The stored
    arrays are the upper triangle already (`/matrix/data` on Li_et_al_2015.h5
    holds 1,661,678 entries, not the 3,313,107 of the symmetrised matrix), so
    the length of `/matrix/data` is nnz_stored directly.
    """
    import collections

    import tables

    with tables.open_file(filename, "r") as handle:
        shape = handle.get_node("/matrix/shape").read()
        chromosomes = handle.get_node("/intervals/chr_list").read()
        counts = collections.Counter(chromosomes.tolist())
        return {
            "nbins": int(shape[0]),
            "nnz_stored": int(handle.get_node("/matrix/data").shape[0]),
            "max_chromosome_bins": max(counts.values()) if counts else 0,
            "group": None,
        }


def working_set_bytes(info):
    """W of PLAN.md 4.5, the stored upper-triangle CSR working set."""
    return (info["nnz_stored"] * (VALUE_BYTES + INDEX_BYTES)
            + (info["nbins"] + 1) * INDPTR_BYTES)


def dense_block_bytes(info):
    """D of PLAN.md 4.5, the largest dense per-chromosome block."""
    return info["max_chromosome_bins"] ** 2 * VALUE_BYTES


# --------------------------------------------------------------------------
# the memory gate


def alpha_beta(case):
    """(alpha, beta, source). The case wins; PLAN.md 4.5's table is the fallback."""
    memory = case.get("memory") or {}
    if "alpha" in memory or "beta" in memory:
        return float(memory.get("alpha", 0.0)), float(memory.get("beta", 0.0)), "case"
    alpha, beta = ALPHA_BETA_BY_TOOL.get(case["tool"], DEFAULT_ALPHA_BETA)
    return alpha, beta, "tool-default"


def case_matrix_inputs(case, data):
    """Every argument of the case that resolves to a readable matrix."""
    found = []
    for argument in case["args"]:
        if "{out}" in argument:
            continue
        info = matrix_info(expand(argument, {"data": data}))
        if info is not None:
            found.append((argument, info))
    return found


def evaluate_memory_gate(case, options, peak_rss_kb, produced_dir, data):
    """The hard gate of PLAN.md 8.3 criterion 4.

    Returns a record that always says what it did, so that a case escaping the
    gate is visible in the report rather than silently passing.
    """
    alpha, beta, source = alpha_beta(case)
    memory = case.get("memory") or {}
    gate = {
        "applied": False,
        "passed": True,
        "alpha": alpha,
        "beta": beta,
        "alpha_source": source,
        "constant_mb": BUDGET_CONSTANT_BYTES / MB,
        "cpp_peak_rss_mb": (peak_rss_kb or 0) * 1024 / MB,
        "warnings": [],
    }
    if options.skip_memory_gate:
        gate["reason"] = "skipped by --skip-memory-gate"
        return gate

    if "budget_mb" in memory:
        if "alpha" in memory or "beta" in memory:
            gate["passed"] = False
            gate["applied"] = True
            gate["reason"] = ("the case declares budget_mb together with "
                              "alpha/beta; PLAN.md 9.2 allows exactly one "
                              "source of truth for the budget")
            return gate
        budget_bytes = float(memory["budget_mb"]) * MB
        gate.update(applied=True, budget_bytes=budget_bytes,
                    w_bytes=None, d_bytes=None, w_source="declared budget_mb",
                    inputs=[])
        return _finish_memory_gate(gate, peak_rss_kb, budget_bytes)

    inputs = case_matrix_inputs(case, data)
    broken = [entry for entry in inputs if "error" in entry[1]]
    inputs = [entry for entry in inputs if "error" not in entry[1]]
    for argument, info in broken:
        gate["warnings"].append(f"could not measure {argument}: {info['error']}")
    w_source = "input"
    if not inputs and produced_dir is not None:
        # A tool whose input is text or hicpro still produces a matrix, and the
        # matrix it produced is the working set it had to hold.
        inputs = []
        for candidate in sorted(Path(produced_dir).rglob("*")):
            if candidate.is_file() and candidate.name.lower().endswith(MATRIX_SUFFIXES):
                info = matrix_info(str(candidate))
                if info is not None and "error" not in info:
                    inputs.append((candidate.name, info))
        if inputs:
            w_source = "python output"

    if not inputs:
        if alpha == 0.0 and beta == 0.0:
            # The budget is the constant alone, so no matrix is needed. This is
            # the hicInfo metadata path and the argument-parsing cases.
            gate.update(applied=True, w_bytes=0, d_bytes=0, inputs=[],
                        w_source="not needed, alpha = beta = 0")
            return _finish_memory_gate(gate, peak_rss_kb, BUDGET_CONSTANT_BYTES)
        gate["reason"] = ("no matrix input or output could be measured, so W "
                          "and D are unknown; declare memory.budget_mb in the "
                          "case to gate it")
        return gate

    w_bytes = max(working_set_bytes(info) for _, info in inputs)
    d_bytes = max(dense_block_bytes(info) for _, info in inputs)
    budget_bytes = alpha * w_bytes + beta * d_bytes + BUDGET_CONSTANT_BYTES

    measured_nnz = max(info["nnz_stored"] for _, info in inputs)
    measured_nbins = max(info["nbins"] for _, info in inputs)
    for field, measured in (("nnz_stored", measured_nnz), ("nbins", measured_nbins)):
        if field in memory and int(memory[field]) != measured:
            gate["warnings"].append(
                f"the case declares {field} = {memory[field]} but the input "
                f"measures {measured}; the measurement is used")

    gate.update(applied=True, w_bytes=w_bytes, d_bytes=d_bytes,
                w_mb=w_bytes / MB, d_mb=d_bytes / MB, w_source=w_source,
                inputs=[{"argument": argument, **{key: value
                                                  for key, value in info.items()
                                                  if key != "path"}}
                        for argument, info in inputs])
    return _finish_memory_gate(gate, peak_rss_kb, budget_bytes)


def _finish_memory_gate(gate, peak_rss_kb, budget_bytes):
    peak_bytes = (peak_rss_kb or 0) * 1024
    gate["budget_bytes"] = budget_bytes
    gate["budget_mb"] = budget_bytes / MB
    gate["budget_kb"] = budget_bytes / 1024.0
    gate["cpp_peak_rss_bytes"] = peak_bytes
    gate["ratio_to_budget"] = peak_bytes / budget_bytes if budget_bytes else None
    gate["headroom"] = 1.0 - gate["ratio_to_budget"] if budget_bytes else None
    gate["passed"] = peak_bytes <= budget_bytes
    if not gate["passed"]:
        gate["reason"] = (
            f"peak RSS {peak_bytes / MB:.1f} MB exceeds the budget "
            f"{budget_bytes / MB:.1f} MB "
            f"({gate['ratio_to_budget'] * 100:.1f} % of it)")
    return gate


# --------------------------------------------------------------------------
# the time gate


def evaluate_time_gate(case, options, py_cpu, cpp_cpu):
    """CPU-time gate: the C++ must not burn more CPU than the Python.

    See the module docstring for why this is CPU time rather than wall clock
    and where the floor comes from.
    """
    gate = {
        "applied": False,
        "passed": True,
        "py_cpu_seconds": py_cpu,
        "cpp_cpu_seconds": cpp_cpu,
        "floor_seconds": TIME_FLOOR_SECONDS,
        "startup_dominated": cpp_cpu < TIME_FLOOR_SECONDS,
    }
    if options.skip_time_gate:
        gate["reason"] = "skipped by --skip-time-gate"
        return gate
    if py_cpu < TIME_FLOOR_SECONDS:
        gate["reason"] = (
            f"the Python reference used {py_cpu:.2f} s of CPU, at or below the "
            f"{TIME_FLOOR_SECONDS:.2f} s measurement floor, so the comparison "
            f"would be quantisation noise")
        return gate
    gate["applied"] = True
    gate["ratio"] = cpp_cpu / py_cpu if py_cpu else None
    if gate["startup_dominated"]:
        # cpp_cpu is one or two clock ticks, so its magnitude means nothing.
        # The verdict is still safe: the gap to py_cpu is far wider than the
        # error. Report the speedup as a lower bound instead of a number that
        # would swing by 150 % between runs.
        gate["speedup"] = None
        gate["speedup_lower_bound"] = py_cpu / TIME_FLOOR_SECONDS
        gate["note"] = (f"the C++ CPU time {cpp_cpu:.2f} s is below the "
                        f"{TIME_FLOOR_SECONDS:.2f} s floor, so the ratio is "
                        f"dominated by process start-up and clock-tick "
                        f"quantisation; the pass verdict still holds")
    else:
        gate["speedup"] = py_cpu / cpp_cpu if cpp_cpu else None
    gate["passed"] = cpp_cpu <= py_cpu
    if not gate["passed"]:
        gate["reason"] = (f"the C++ used {cpp_cpu:.2f} s of CPU against the "
                          f"Python's {py_cpu:.2f} s")
    return gate


# --------------------------------------------------------------------------
# measured execution


class Measurement(dict):
    """seconds, peak_rss_kb, user_seconds, sys_seconds, exit_code."""


def _descendant_pids(root):
    pids = []
    pending = [root]
    while pending:
        pid = pending.pop()
        pids.append(pid)
        try:
            with open(f"/proc/{pid}/task/{pid}/children", encoding="ascii") as handle:
                pending.extend(int(part) for part in handle.read().split())
        except OSError:
            continue
    return pids


def _vm_hwm_kb(pid):
    try:
        with open(f"/proc/{pid}/status", encoding="ascii") as handle:
            for line in handle:
                if line.startswith("VmHWM:"):
                    return int(line.split()[1])
    except (OSError, ValueError, IndexError):
        pass
    return 0


def _trace_rss(popen, trace_path, stop):
    """Sample VmHWM of the process tree every 100 ms (PLAN.md 10).

    A single peak says a tool went over budget; the trace says where.
    """
    start = time.perf_counter()
    try:
        with open(trace_path, "w", encoding="ascii") as handle:
            handle.write("seconds\tvm_hwm_kb\n")
            while not stop.wait(RSS_TRACE_INTERVAL):
                if popen.poll() is not None:
                    break
                total = max((_vm_hwm_kb(pid) for pid in _descendant_pids(popen.pid)),
                            default=0)
                if total:
                    handle.write(f"{time.perf_counter() - start:.3f}\t{total}\n")
                    handle.flush()
    except OSError:
        return


def run_measured(argv, cwd, stdout_path, stderr_path, env=None, trace_path=None):
    """Run a command under /usr/bin/time and capture stdout and stderr."""
    time_file = stdout_path.parent / (stdout_path.name + ".time")
    wrapper = ["/usr/bin/time", "-f", "%e %M %U %S", "-o", str(time_file)]
    if not Path(wrapper[0]).exists():
        wrapper = []
    start = time.perf_counter()
    with open(stdout_path, "wb") as out, open(stderr_path, "wb") as err:
        popen = subprocess.Popen(wrapper + list(argv), cwd=str(cwd), stdout=out,
                                 stderr=err, env=env)
        tracer, stop = None, None
        if trace_path is not None:
            stop = threading.Event()
            tracer = threading.Thread(target=_trace_rss,
                                      args=(popen, trace_path, stop), daemon=True)
            tracer.start()
        returncode = popen.wait()
        if tracer is not None:
            stop.set()
            tracer.join(timeout=1.0)
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
                       cpu_seconds=user + system,
                       exit_code=returncode,
                       command=" ".join(shlex.quote(part) for part in argv))


# --------------------------------------------------------------------------
# determinism (PLAN.md 8.3 criterion 3)


def _declared_outputs(case):
    """The declared outputs as (name relative to the work directory, format)."""
    found = []
    for declared in case["outputs"]:
        path = declared["path"]
        if "{out}" not in path:
            continue
        found.append((path.replace("{out}/", "").replace("{out}", ""),
                      declared.get("format", "plain")))
    return found


def _declared_options(case, name):
    """The comparator options a case declares for the output named `name`."""
    for declared in case["outputs"]:
        path = declared["path"].replace("{out}/", "").replace("{out}", "")
        if path == name:
            return declared.get("options") or {}
    return {}


# The strictest class each format admits, used to say *what* differs when two
# runs are not byte-identical. A cool file carries a `creation-date` attribute,
# so two writes of the same matrix can never be byte-identical, whoever writes
# them; PLAN.md 5.1 already rules that field and three other provenance fields
# non-significant, and class E1 still requires every dataset to decode to
# identical bytes, which is what a reduction-order difference would break.
STRICTEST_CLASS_BY_FORMAT = {"cool": "E1", "mcool": "E1", "h5": "E2"}

# An HDF5 object header carries an optional modification time (message type
# 0x12, H5O_MTIME_NEW): four reserved bytes, then a 4-byte Unix time. It is
# written unless the creator clears `H5Pset_obj_track_times`, and it makes
# every output differ between two runs a second apart. That is a defect with a
# one-line fix, not a provenance field anyone agreed to, so the determinism
# check names it and fails on it rather than normalising it away.
HDF5_MTIME_PREFIX = b"\x12\x00\x08\x00"
HDF5_MTIME_VERSION = b"\x01\x00\x00\x00"


def _byte_differences(left, right, max_offsets=64):
    """Offsets at which two files differ, capped. (offsets, truncated)."""
    offsets = []
    chunk = 1 << 20
    with open(left, "rb") as handle_a, open(right, "rb") as handle_b:
        position = 0
        while True:
            block_a = handle_a.read(chunk)
            block_b = handle_b.read(chunk)
            if not block_a and not block_b:
                return offsets, False
            if block_a != block_b:
                shared = min(len(block_a), len(block_b))
                for start in range(0, shared, 4096):
                    end = min(start + 4096, shared)
                    if block_a[start:end] == block_b[start:end]:
                        continue
                    for index in range(start, end):
                        if block_a[index] != block_b[index]:
                            offsets.append(position + index)
                            if len(offsets) >= max_offsets:
                                return offsets, True
                if len(block_a) != len(block_b):
                    offsets.append(position + shared)
                    return offsets, True
            position += len(block_a)


def _is_hdf5_mtime_byte(handle, offset):
    """Does this offset fall inside an HDF5 object modification time?"""
    start = max(0, offset - 14)
    handle.seek(start)
    window = handle.read(offset - start + 1)
    base = offset - start
    for back in range(0, 4):
        head = base - back - 12
        if head < 0:
            continue
        if window[head:head + 4] == HDF5_MTIME_PREFIX \
                and window[head + 8:head + 12] == HDF5_MTIME_VERSION:
            return True
    return False


def _classify_byte_differences(path, offsets):
    """Split the differing offsets into HDF5 object times and everything else."""
    mtime, other = 0, 0
    try:
        with open(path, "rb") as handle:
            for offset in offsets:
                if _is_hdf5_mtime_byte(handle, offset):
                    mtime += 1
                else:
                    other += 1
    except OSError:
        return 0, len(offsets)
    return mtime, other


def _compare_run_outputs(case, reference_dir, other_dir):
    """Compare every declared output of two C++ runs.

    Byte-identity is the check. When it fails, the difference is diagnosed
    rather than merely reported: the strictest class of the format says whether
    the *content* is identical, and the differing bytes are classified, so that
    an HDF5 object timestamp is never mistaken for a numeric difference and,
    just as important, a numeric difference is never excused as a timestamp.
    Only a difference confined to the provenance fields PLAN.md 5.1 already
    normalises is downgraded to a qualification.
    """
    diffs, qualified = [], []
    for name, fmt in _declared_outputs(case):
        left = Path(reference_dir) / name
        right = Path(other_dir) / name
        if left.exists() != right.exists():
            diffs.append(f"{name}: produced by only one of the two runs")
            continue
        if not left.exists():
            continue
        if filecmp.cmp(str(left), str(right), shallow=False):
            continue

        offsets, truncated = _byte_differences(str(left), str(right))
        mtime_bytes, other_bytes = _classify_byte_differences(str(left), offsets)
        where = (f"{len(offsets)}{'+' if truncated else ''} differing bytes, "
                 f"first at offset {offsets[0]}" if offsets else "sizes differ")

        # A text output whose case declares a named normalisation, such as the
        # random temporary matrix name hicQuickQC prints, is compared through
        # that normalisation and nothing else: a difference it absorbs is a
        # qualification, any other byte is a failure.
        normalise = _declared_options(case, name).get("normalise")
        if normalise and fmt in ("plain", "text"):
            comparison = comparators.compare(fmt, str(left), str(right), "E0",
                                             {"normalise": normalise})
            if comparison.passed:
                qualified.append(f"{name}: identical after the named normalisation "
                                 f"{normalise} but not byte for byte ({where})")
            else:
                diffs.append(f"{name}: differs after the named normalisation "
                             f"{normalise}: " + "; ".join(comparison.diffs[:3]))
            continue

        strictest = STRICTEST_CLASS_BY_FORMAT.get(fmt)
        content_identical = None
        if strictest is not None:
            comparison = comparators.compare(fmt, str(left), str(right),
                                             strictest, None)
            content_identical = comparison.passed
            if not content_identical:
                diffs.append(f"{name}: differs at {strictest}, which is a real "
                             f"difference in the data: "
                             + "; ".join(comparison.diffs[:3]))
                continue

        if content_identical:
            cause = (f"{mtime_bytes} of them are HDF5 object modification times "
                     f"(message 0x12), {other_bytes} are provenance fields"
                     if mtime_bytes else
                     "the difference is confined to the provenance fields "
                     "PLAN.md 5.1 normalises")
            qualified.append(f"{name}: identical at {strictest} but not byte for "
                             f"byte ({where}); {cause}")
            continue
        diffs.append(f"{name}: not byte-identical ({where}), and the format has "
                     f"no structural comparator to diagnose it")
    return diffs, qualified


def _embedded_wall_clock(path):
    """Count HDF5 object modification times that hold a plausible wall clock.

    A file that embeds the time it was written cannot be byte-identical to the
    same file written a second later, so this is checked directly rather than
    left to whether two repeats happened to straddle a second boundary. HDF5
    writes the message unless the creator clears `H5Pset_obj_track_times`.
    """
    import struct

    pattern = HDF5_MTIME_PREFIX + b"\x00\x00\x00\x00" + HDF5_MTIME_VERSION
    count, example = 0, None
    # One byte short of a whole match, so a match straddling a block boundary
    # is found once and a match at the end of a block is not counted twice.
    overlap = len(pattern) + 3
    carry = b""
    try:
        with open(path, "rb") as handle:
            while True:
                block = handle.read(8 << 20)
                if not block:
                    break
                data = carry + block
                position = 0
                while True:
                    index = data.find(pattern, position)
                    if index < 0 or index + len(pattern) + 4 > len(data):
                        break
                    stamp, = struct.unpack_from("<I", data, index + len(pattern))
                    # A plausible recent Unix time, so that compressed payload
                    # bytes matching the pattern by chance are not counted.
                    if 1_400_000_000 < stamp < 2_500_000_000:
                        count += 1
                        example = example or stamp
                    position = index + 1
                carry = data[-overlap:] if len(data) > overlap else data
    except OSError:
        return 0, None
    return count, example


def check_determinism(case, options, workdir, cpp_tool, data, reference_dir):
    """Repeat runs byte-identical, and --threads 1 identical to --threads N.

    `reference_dir` is the output directory of the run that has already
    happened, so the first repeat is free.
    """
    record = {
        "repeats": options.noise_runs,
        "repeats_passed": True,
        "byte_identical": True,
        "reproducible_output": True,
        "embedded_wall_clock": [],
        "threads_arg": case.get("threads_arg"),
        "threads_1_vs_high": None,
        "threads_high": options.threads_high,
        "passed": True,
        "diffs": [],
        "qualifications": [],
    }

    # Deterministic, and independent of whether two repeats happened to fall in
    # the same second: does the output embed the wall clock at all?
    for name, fmt in _declared_outputs(case):
        if fmt not in STRICTEST_CLASS_BY_FORMAT:
            continue
        produced = Path(reference_dir) / name
        if not produced.exists():
            continue
        count, example = _embedded_wall_clock(str(produced))
        if count:
            record["reproducible_output"] = False
            record["embedded_wall_clock"].append(
                {"output": name, "count": count, "example": example})
            record["diffs"].append(
                f"{name} embeds {count} HDF5 object modification times "
                f"(example {example}), so two runs a second apart cannot be "
                f"byte-identical. Clear H5Pset_obj_track_times on the file, "
                f"group and dataset creation property lists. The PyTables and "
                f"cooler reference files carry the same messages, so this is "
                f"the HDF5 default rather than a regression")

    def run_into(directory, extra_args):
        directory.mkdir(parents=True, exist_ok=True)
        argv = [expand(argument, {"data": data, "out": directory})
                for argument in case["args"] + case.get("cpp_args", [])] \
            + list(extra_args)
        return run_measured([str(cpp_tool)] + argv, workdir,
                            directory / "stdout.txt", directory / "stderr.txt")

    for index in range(1, max(1, options.noise_runs)):
        repeat_dir = workdir / f"out_cpp_repeat{index}"
        run_into(repeat_dir, [])
        diffs, qualified = _compare_run_outputs(case, reference_dir, repeat_dir)
        if qualified:
            record["byte_identical"] = False
            for note in qualified:
                if note not in record["qualifications"]:
                    record["qualifications"].append(note)
        if diffs:
            record["repeats_passed"] = False
            record["diffs"].extend(f"repeat {index}: {diff}" for diff in diffs)
        else:
            shutil.rmtree(repeat_dir, ignore_errors=True)

    threads_arg = case.get("threads_arg")
    if threads_arg:
        low_dir = workdir / "out_cpp_threads1"
        high_dir = workdir / f"out_cpp_threads{options.threads_high}"
        run_into(low_dir, [threads_arg, "1"])
        run_into(high_dir, [threads_arg, str(options.threads_high)])
        diffs, qualified = _compare_run_outputs(case, low_dir, high_dir)
        if qualified:
            record["byte_identical"] = False
            for note in qualified:
                entry = f"{threads_arg} 1 against {options.threads_high}: {note}"
                if entry not in record["qualifications"]:
                    record["qualifications"].append(entry)
        record["threads_1_vs_high"] = not diffs
        if diffs:
            record["diffs"].extend(
                f"{threads_arg} 1 against {options.threads_high}: {diff}"
                for diff in diffs)
        else:
            shutil.rmtree(low_dir, ignore_errors=True)
            shutil.rmtree(high_dir, ignore_errors=True)

    record["passed"] = (record["repeats_passed"]
                        and record["reproducible_output"]
                        and record["threads_1_vs_high"] is not False)
    return record


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
    args_cpp = [expand(arg, {"data": data, "out": out_cpp})
                for arg in case["args"] + case["cpp_args"]]

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
        "timing_unreliable": load_average > LOAD_AVERAGE_LIMIT,
        "load_average": load_average,
        "large": bool(case.get("large")),
    }

    if not cpp_tool.exists():
        result.update(passed=False, class_met=None, outputs=[],
                      error=f"missing C++ binary {cpp_tool}")
        return result

    trace_py = (workdir / "rss_py.tsv") if case.get("large") else None
    trace_cpp = (workdir / "rss_cpp.tsv") if case.get("large") else None
    measure_py = run_measured([str(options.py_python), str(python_tool)] + args_py,
                              workdir, out_py / "stdout.txt", out_py / "stderr.txt",
                              env, trace_path=trace_py)
    measure_cpp = run_measured([str(cpp_tool)] + args_cpp, workdir,
                               out_cpp / "stdout.txt", out_cpp / "stderr.txt",
                               trace_path=trace_cpp)

    result.update({
        "py_seconds": measure_py["seconds"],
        "cpp_seconds": measure_cpp["seconds"],
        "py_peak_rss_kb": measure_py["peak_rss_kb"],
        "cpp_peak_rss_kb": measure_cpp["peak_rss_kb"],
        "py_user_seconds": measure_py["user_seconds"],
        "cpp_user_seconds": measure_cpp["user_seconds"],
        "py_sys_seconds": measure_py["sys_seconds"],
        "cpp_sys_seconds": measure_cpp["sys_seconds"],
        "py_cpu_seconds": measure_py["cpu_seconds"],
        "cpp_cpu_seconds": measure_cpp["cpu_seconds"],
        "py_exit": measure_py["exit_code"],
        "cpp_exit": measure_cpp["exit_code"],
        "py_command": measure_py["command"],
        "cpp_command": measure_cpp["command"],
        "stderr_py": _tail(out_py / "stderr.txt"),
        "stderr_cpp": _tail(out_cpp / "stderr.txt"),
    })

    passed = True
    errors = []
    if measure_py["exit_code"] != case["expect_exit"]:
        passed = False
        errors.append(f"the Python tool exited {measure_py['exit_code']}, "
                      f"expected {case['expect_exit']}")
    if measure_cpp["exit_code"] != case["expect_exit"]:
        passed = False
        errors.append(f"the C++ tool exited {measure_cpp['exit_code']}, "
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
    result["class_met"] = result["class_declared"] if passed else None

    # --- the gates ---------------------------------------------------------
    memory_gate = evaluate_memory_gate(case, options, measure_cpp["peak_rss_kb"],
                                       out_py, data)
    memory_gate["ratio_to_python"] = (
        measure_cpp["peak_rss_kb"] / measure_py["peak_rss_kb"]
        if measure_py["peak_rss_kb"] else None)
    result["memory_gate"] = memory_gate
    result["budget_kb"] = memory_gate.get("budget_kb")
    if not memory_gate["passed"]:
        passed = False
        errors.append("memory gate: " + memory_gate.get("reason", "over budget"))

    time_gate = evaluate_time_gate(case, options, measure_py["cpu_seconds"],
                                   measure_cpp["cpu_seconds"])
    result["time_gate"] = time_gate
    if not time_gate["passed"]:
        passed = False
        errors.append("time gate: " + time_gate.get("reason", "too slow"))

    if options.determinism:
        determinism = check_determinism(case, options, workdir, cpp_tool, data,
                                        out_cpp)
        result["determinism"] = determinism
        if not determinism["passed"]:
            passed = False
            errors.append("determinism: " + "; ".join(determinism["diffs"][:3]))

    result["passed"] = passed
    if errors:
        result["error"] = "; ".join(errors)

    if not options.keep_workdirs and passed:
        shutil.rmtree(workdir, ignore_errors=True)
        result["workdir"] = None
    return result


def run_determinism_case(case, options):
    """The determinism mode: the C++ against itself, no Python, no comparators."""
    workdir = Path(tempfile.mkdtemp(prefix=f"determinism-{case['id']}-",
                                    dir=options.tmpdir))
    out_cpp = workdir / "out_cpp"
    out_cpp.mkdir()
    data = str(options.data)
    cpp_tool = (Path(options.cpp_bin) / case["cpp_binary"]).resolve() \
        if case.get("cpp_binary") else (Path(options.cpp_bin) / case["tool"])

    result = {
        "id": case["id"],
        "tool": case["tool"],
        "tier": case["tier"],
        "class_declared": _declared_class(case),
        "notes": case["notes"],
        "workdir": str(workdir),
        "load_average": os.getloadavg()[0],
        "large": bool(case.get("large")),
    }
    if not cpp_tool.exists():
        result.update(passed=False, error=f"missing C++ binary {cpp_tool}")
        return result

    args_cpp = [expand(arg, {"data": data, "out": out_cpp})
                for arg in case["args"] + case.get("cpp_args", [])]
    measure = run_measured([str(cpp_tool)] + args_cpp, workdir,
                           out_cpp / "stdout.txt", out_cpp / "stderr.txt")
    result["cpp_seconds"] = measure["seconds"]
    result["cpp_cpu_seconds"] = measure["cpu_seconds"]
    result["cpp_peak_rss_kb"] = measure["peak_rss_kb"]
    result["cpp_command"] = measure["command"]
    if measure["exit_code"] != case["expect_exit"]:
        result.update(passed=False,
                      error=f"the C++ tool exited {measure['exit_code']}, "
                            f"expected {case['expect_exit']}")
        return result

    determinism = check_determinism(case, options, workdir, cpp_tool, data, out_cpp)
    result["determinism"] = determinism
    result["passed"] = determinism["passed"]
    if not determinism["passed"]:
        result["error"] = "; ".join(determinism["diffs"][:5])
    if not options.keep_workdirs and result["passed"]:
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


def _percent(value):
    return "n/a" if value is None else f"{value * 100:.1f} %"


def _memory_summary(lines, cases):
    """PLAN.md 9.4: every case sorted by cpp_peak_rss / budget, tightest first."""
    gated = [case for case in cases
             if (case.get("memory_gate") or {}).get("ratio_to_budget") is not None]
    ungated = [case for case in cases
               if (case.get("memory_gate") or {}).get("ratio_to_budget") is None]
    lines.append("## Memory gate")
    lines.append("")
    lines.append("`budget = alpha*W + beta*D + C`, `C = 64 MB`, SI MB throughout "
                 "(PLAN.md 4.5). `W` and `D` are measured from the actual input, "
                 "not taken from the case. Sorted by how close the tool is to "
                 "its gate.")
    lines.append("")
    lines.append("| case | alpha | beta | source | W MB | D MB | budget MB | "
                 "cpp RSS MB | of budget | headroom | py RSS MB | result |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for case in sorted(gated, key=lambda item: -item["memory_gate"]["ratio_to_budget"]):
        gate = case["memory_gate"]
        lines.append(
            "| {id} | {alpha} | {beta} | {source} | {w} | {d} | {budget:.1f} | "
            "{rss:.1f} | {ratio} | {headroom} | {pyrss:.1f} | {result} |".format(
                id=case["id"], alpha=gate["alpha"], beta=gate["beta"],
                source=gate["alpha_source"],
                w="-" if gate.get("w_mb") is None else f"{gate['w_mb']:.1f}",
                d="-" if gate.get("d_mb") is None else f"{gate['d_mb']:.1f}",
                budget=gate["budget_mb"], rss=gate["cpp_peak_rss_mb"],
                ratio=_percent(gate["ratio_to_budget"]),
                headroom=_percent(gate["headroom"]),
                pyrss=case.get("py_peak_rss_kb", 0) * 1024 / MB,
                result="pass" if gate["passed"] else "**FAIL**"))
    lines.append("")
    if ungated:
        lines.append("Not gated:")
        lines.append("")
        for case in ungated:
            gate = case.get("memory_gate") or {}
            lines.append(f"- `{case['id']}`: {gate.get('reason', 'no gate record')}")
        lines.append("")
    warned = [(case["id"], warning) for case in cases
              for warning in (case.get("memory_gate") or {}).get("warnings", [])]
    if warned:
        lines.append("Warnings:")
        lines.append("")
        for case_id, warning in warned:
            lines.append(f"- `{case_id}`: {warning}")
        lines.append("")


def _time_summary(lines, cases):
    lines.append("## Time gate")
    lines.append("")
    lines.append("Gated on CPU time (user + sys): the C++ must not exceed the "
                 "Python. Wall clock and the load average are recorded but not "
                 "gated, because this machine is never idle. A C++ CPU time "
                 f"below the {TIME_FLOOR_SECONDS:.2f} s measurement floor is "
                 "marked `startup` and its speedup is a lower bound.")
    lines.append("")
    lines.append("| case | py CPU s | cpp CPU s | cpp/py | CPU speedup | "
                 "py wall s | cpp wall s | load | result |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    ordered = sorted(cases, key=lambda item: -((item.get("time_gate") or {}).get("ratio") or 0.0))
    for case in ordered:
        gate = case.get("time_gate") or {}
        if gate.get("speedup") is not None:
            speedup = f"{gate['speedup']:.1f}x"
        elif gate.get("speedup_lower_bound") is not None:
            speedup = f"> {gate['speedup_lower_bound']:.0f}x (startup)"
        else:
            speedup = "n/a"
        result = "pass" if gate.get("passed") else "**FAIL**"
        if not gate.get("applied"):
            result = "not gated"
        lines.append(
            "| {id} | {py:.2f} | {cpp:.2f} | {ratio} | {speedup} | {pyw:.2f} | "
            "{cppw:.2f} | {load:.1f}{flag} | {result} |".format(
                id=case["id"], py=gate.get("py_cpu_seconds", 0.0),
                cpp=gate.get("cpp_cpu_seconds", 0.0),
                ratio="n/a" if gate.get("ratio") is None else f"{gate['ratio']:.3f}",
                speedup=speedup, pyw=case.get("py_seconds", 0.0),
                cppw=case.get("cpp_seconds", 0.0),
                load=case.get("load_average", 0.0),
                flag=" (unreliable)" if case.get("timing_unreliable") else "",
                result=result))
    lines.append("")


def _determinism_summary(lines, cases):
    checked = [case for case in cases if case.get("determinism")]
    if not checked:
        return
    lines.append("## Determinism")
    lines.append("")
    lines.append("PLAN.md 8.3 criterion 3: repeat C++ runs byte-identical, and "
                 "for a tool that takes `--threads`, `--threads 1` byte-identical "
                 "to the higher thread count. A byte difference is diagnosed, not "
                 "excused: the strictest class of the format says whether the "
                 "content differs, and only a difference confined to the "
                 "provenance fields PLAN.md 5.1 normalises is a qualification "
                 "rather than a failure. An output that embeds the wall clock "
                 "fails on its own, whether or not two repeats happened to "
                 "straddle a second boundary.")
    lines.append("")
    lines.append("| case | repeats | repeats identical | embeds wall clock | "
                 "threads 1 vs N | result |")
    lines.append("|---|---|---|---|---|---|")
    for case in checked:
        record = case["determinism"]
        threads = record["threads_1_vs_high"]
        if not record["repeats_passed"]:
            repeats = "**no**"
        elif record.get("byte_identical", True):
            repeats = "yes, byte for byte"
        else:
            repeats = "content only"
        stamps = sum(entry["count"] for entry in record.get("embedded_wall_clock", []))
        lines.append("| {id} | {repeats} | {rep} | {stamps} | {threads} | "
                     "{result} |".format(
                         id=case["id"], repeats=record["repeats"], rep=repeats,
                         stamps="no" if not stamps else f"**{stamps} times**",
                         threads="not threaded" if threads is None
                         else ("yes" if threads else "**no**"),
                         result="pass" if record["passed"] else "**FAIL**"))
    lines.append("")
    qualified = [(case["id"], note) for case in checked
                 for note in case["determinism"].get("qualifications", [])]
    if qualified:
        lines.append("Qualified, not byte for byte:")
        lines.append("")
        for case_id, note in qualified:
            lines.append(f"- `{case_id}`: {note}")
        lines.append("")


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
    if report.get("memory_gate") == "skipped":
        lines.append("- **memory gate: skipped. This report cannot record a "
                     "pass (PLAN.md 9.1).**")
    if report.get("time_gate") == "skipped":
        lines.append("- **time gate: skipped. This report cannot record a pass.**")
    lines.append("")

    if any(case.get("memory_gate") for case in cases):
        _memory_summary(lines, cases)
    if any(case.get("time_gate") for case in cases):
        _time_summary(lines, cases)
    _determinism_summary(lines, cases)

    tiers = sorted({case["tier"] for case in cases})
    for tier in tiers:
        lines.append(f"## Tier {tier}")
        lines.append("")
        lines.append("Wall clock and the speedup derived from it are indicative "
                     "only; the gates are the CPU time and the budget columns.")
        lines.append("")
        lines.append("| case | class | result | py wall s | cpp wall s | "
                     "wall speedup | cpu speedup | py RSS MB | cpp RSS MB | "
                     "budget MB | of budget |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for case in [case for case in cases if case["tier"] == tier]:
            gate = case.get("memory_gate") or {}
            time_gate = case.get("time_gate") or {}
            speedup = _ratio(case.get("py_seconds"), case.get("cpp_seconds"))
            if time_gate.get("speedup") is not None:
                cpu_speedup = f"{time_gate['speedup']:.1f}x"
            elif time_gate.get("speedup_lower_bound") is not None:
                cpu_speedup = f"> {time_gate['speedup_lower_bound']:.0f}x"
            else:
                cpu_speedup = "n/a"
            lines.append(
                "| {id} | {cls} | {result} | {py:.3f} | {cpp:.3f} | {speedup} | "
                "{cpu} | {pyrss:.1f} | {cpprss:.1f} | {budget} | {ratio} |".format(
                    id=case["id"],
                    cls=case.get("class_declared", "?"),
                    result="pass" if case.get("passed") else "FAIL",
                    py=case.get("py_seconds", 0.0),
                    cpp=case.get("cpp_seconds", 0.0),
                    speedup=speedup,
                    cpu=cpu_speedup,
                    pyrss=case.get("py_peak_rss_kb", 0) * 1024 / MB,
                    cpprss=case.get("cpp_peak_rss_kb", 0) * 1024 / MB,
                    budget="-" if gate.get("budget_mb") is None
                    else f"{gate['budget_mb']:.1f}",
                    ratio=_percent(gate.get("ratio_to_budget"))))
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
            for diff in (case.get("determinism") or {}).get("diffs", [])[:20]:
                lines.append(f"- determinism: {diff}")
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


def _is_memory_gated(case, options):
    """Is this case's peak RSS a gate rather than a report line?

    Every case is, unless the gate is switched off for development: the budget
    formula has a value for every tool (PLAN.md 4.5, "everything else"), so
    there is no case whose memory is merely informational. A gated case must
    not share the machine with another (PLAN.md 9.1): a peak-RSS measurement
    taken while other cases compete for memory is not a measurement.
    """
    return not getattr(options, "skip_memory_gate", False)


def _execute(cases, options, runner):
    """Gated cases serially, the rest in the pool."""
    serial = [case for case in cases if _is_memory_gated(case, options)]
    serial_ids = {case["id"] for case in serial}
    parallel = [case for case in cases if case["id"] not in serial_ids]
    results = [runner(case, options) for case in serial]
    if parallel:
        if options.jobs > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=options.jobs) as pool:
                futures = [pool.submit(runner, case, options) for case in parallel]
                results.extend(future.result() for future in futures)
        else:
            results.extend(runner(case, options) for case in parallel)
    results.sort(key=lambda item: (item["tier"], item["id"]))
    return results


def _report_skeleton(options, results, mode):
    report = {
        "harness_version": HARNESS_VERSION,
        "mode": mode,
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": git_commit(),
        "host": platform.node(),
        "cpp_bin": str(options.cpp_bin),
        "py_python": str(getattr(options, "py_python", "")),
        "budget_constant_mb": BUDGET_CONSTANT_BYTES / MB,
        "time_floor_seconds": TIME_FLOOR_SECONDS,
        "cases": results,
    }
    if getattr(options, "skip_memory_gate", False):
        report["memory_gate"] = "skipped"
    if getattr(options, "skip_time_gate", False):
        report["time_gate"] = "skipped"
    return report


def command_run(options):
    cases = load_cases(options.tool, options.tier, options.case)
    if not cases:
        print("no cases selected", file=sys.stderr)
        return 1
    results = _execute(cases, options, run_case)
    report = _report_skeleton(options, results, "run")
    json_path, markdown_path = write_report(report, options.out)

    for case in results:
        status = "pass" if case.get("passed") else "FAIL"
        gate = case.get("memory_gate") or {}
        time_gate = case.get("time_gate") or {}
        ratio = gate.get("ratio_to_budget")
        budget = ("      n/a" if gate.get("budget_mb") is None
                  else f"{gate['budget_mb']:9.1f}")
        share = "    n/a" if ratio is None else f"{ratio * 100:6.1f}%"
        print(f"{status:4} {case['id']:<48} "
              f"cpu {case.get('py_cpu_seconds', 0):6.2f} / "
              f"{case.get('cpp_cpu_seconds', 0):6.2f} s  "
              f"rss {case.get('py_peak_rss_kb', 0) * 1024 / MB:7.1f} / "
              f"{case.get('cpp_peak_rss_kb', 0) * 1024 / MB:7.1f} MB  "
              f"budget {budget} MB {share}"
              f"{'  time gate not applied' if not time_gate.get('applied') else ''}")
        if not case.get("passed"):
            if case.get("error"):
                print(f"       {case['error']}")
            for output in case.get("outputs", []):
                for diff in output.get("diffs", [])[:5]:
                    print(f"       {diff}")
        for warning in gate.get("warnings", []):
            print(f"       warning: {warning}")
    if report.get("memory_gate") == "skipped":
        print("\nthe memory gate was skipped: this run cannot record a pass")
    print(f"\nreport: {json_path}\n        {markdown_path}")
    return 0 if all(case.get("passed") for case in results) else 1


def command_determinism(options):
    cases = load_cases(options.tool, options.tier, options.case)
    if not cases:
        print("no cases selected", file=sys.stderr)
        return 1
    options.determinism = True
    results = _execute(cases, options, run_determinism_case)
    report = _report_skeleton(options, results, "determinism")
    json_path, markdown_path = write_report(report, options.out)
    threaded = [case for case in results
                if (case.get("determinism") or {}).get("threads_1_vs_high") is not None]
    for case in results:
        record = case.get("determinism") or {}
        status = "pass" if case.get("passed") else "FAIL"
        if not record.get("repeats_passed"):
            repeats = "NO"
        elif record.get("byte_identical", True):
            repeats = "yes, byte for byte"
        else:
            repeats = "content only"
        if not record.get("reproducible_output", True):
            repeats += ", output embeds the wall clock"
        print(f"{status:4} {case['id']:<48} "
              f"{record.get('repeats', 0)} repeats identical: {repeats}  "
              f"threads: {'not threaded' if record.get('threads_1_vs_high') is None else ('identical' if record['threads_1_vs_high'] else 'DIFFER')}")
        if case.get("error"):
            print(f"       {case['error']}")
    print(f"\n{len(results)} cases, {len(threaded)} of them threaded "
          f"(the thread half of criterion 3 is a no-op until a threaded tool "
          f"lands and its case declares threads_arg)")
    print(f"report: {json_path}\n        {markdown_path}")
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
        alpha, beta, source = alpha_beta(case)
        print(f"{case['tier']}  {case['id']:<48} {case['tool']:<16} "
              f"{_declared_class(case)}  alpha={alpha} beta={beta} ({source})  "
              f"{case['notes']}")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_selection(target):
        target.add_argument("--tool", action="append")
        target.add_argument("--tier", action="append", type=int)
        target.add_argument("--case", action="append")
        target.add_argument("--cpp-bin", default=str(DEFAULT_CPP_BIN))
        target.add_argument("--data", default=str(DEFAULT_DATA))
        target.add_argument("--out", default=str(DEFAULT_OUT))
        target.add_argument("--tmpdir", default=tempfile.gettempdir())
        target.add_argument("--jobs", type=int, default=1)
        target.add_argument("--keep-workdirs", action="store_true")
        target.add_argument("--noise-runs", type=int, default=DEFAULT_NOISE_RUNS,
                            help="repeat runs for the determinism check and for "
                                 "class EN (PLAN.md 9.1)")
        target.add_argument("--threads-high", type=int, default=DEFAULT_THREADS_HIGH,
                            help="the high thread count in the determinism check")

    run_parser = subparsers.add_parser("run", help="run cases and compare")
    add_selection(run_parser)
    run_parser.add_argument("--py-python", default=str(DEFAULT_PY_PYTHON))
    run_parser.add_argument("--skip-memory-gate", action="store_true",
                            help="development only; the report is marked and "
                                 "cannot record a pass")
    run_parser.add_argument("--skip-time-gate", action="store_true",
                            help="development only")
    run_parser.add_argument("--determinism", action="store_true",
                            help="also run the determinism check of PLAN.md "
                                 "8.3 criterion 3 inside the run")
    run_parser.set_defaults(handler=command_run)

    determinism_parser = subparsers.add_parser(
        "determinism", help="C++ against itself: repeat runs and thread counts")
    add_selection(determinism_parser)
    determinism_parser.set_defaults(handler=command_determinism, determinism=True,
                                    skip_memory_gate=False, skip_time_gate=False,
                                    py_python=str(DEFAULT_PY_PYTHON))

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
    if not hasattr(options, "determinism"):
        options.determinism = False
    return options.handler(options)


if __name__ == "__main__":
    sys.exit(main())
