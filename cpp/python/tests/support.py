"""Shared helpers of the hicx_matrix tests.

Environment variables:
    HICX_PYTHON_MODULE_DIR       directory holding the built hicx_matrix module
                                 (<build>/python); PYTHONPATH works as well
    HICX_TEST_DATA               test data directory, default
                                 hicexplorer/test/test_data of this repository
    HICX_ORACLE_PYTHON_COOLER    interpreter with cooler and hicmatrix
    HICX_ORACLE_PYTHON_HICSTRAW  interpreter with hicstraw 1.3.1
    HICX_ORACLE_PYTHON           fallback for both; default: this interpreter
    HICX_LARGE_HIC               a large .hic file (tests skip without it)
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[3]
TEST_DATA = Path(os.environ.get("HICX_TEST_DATA", REPO / "hicexplorer" / "test" / "test_data"))
ORACLE = Path(__file__).with_name("oracle.py")
RTOL = 1e-6

# oracle name -> [compared arrays, compared values, non-zero values, max relative difference]
SUMMARY = {}


def import_hicx_matrix():
    module_dir = os.environ.get("HICX_PYTHON_MODULE_DIR")
    if module_dir and module_dir not in sys.path:
        sys.path.insert(0, module_dir)
    return pytest.importorskip("hicx_matrix")


def data(relative):
    return str(TEST_DATA / relative)


def oracle_python(kind):
    specific = {"cooler": "HICX_ORACLE_PYTHON_COOLER",
                "hicmatrix": "HICX_ORACLE_PYTHON_COOLER",
                "hicstraw": "HICX_ORACLE_PYTHON_HICSTRAW"}[kind]
    return os.environ.get(specific) or os.environ.get("HICX_ORACLE_PYTHON") or sys.executable


def run_oracle(kind, queries, directory):
    """Runs oracle.py once for all queries; returns {id: {key: array}}."""
    query_file = Path(directory) / f"{kind}.json"
    out_file = Path(directory) / f"{kind}.npz"
    query_file.write_text(json.dumps(queries))
    proc = subprocess.run([oracle_python(kind), str(ORACLE), kind, str(query_file), str(out_file)],
                          capture_output=True, text=True)
    if proc.returncode == 3:
        pytest.skip(f"{kind} oracle unavailable: {proc.stderr.strip()}")
    if proc.returncode != 0:
        raise RuntimeError(f"{kind} oracle failed ({proc.returncode}):\n{proc.stderr[-4000:]}")
    results = {}
    with np.load(out_file) as archive:
        for name in archive.files:
            query_id, key = name.rsplit("__", 1)
            results.setdefault(query_id, {})[key] = archive[name]
    return results


def assert_e2(actual, expected, oracle):
    """E2: equal shape and NaN positions, relative tolerance 1e-6 elsewhere."""
    actual = np.asarray(actual)
    expected = np.asarray(expected, dtype=np.float64)
    assert actual.dtype == np.float64
    assert actual.shape == expected.shape, f"shape {actual.shape} != reference {expected.shape}"
    nan_actual = np.isnan(actual)
    nan_expected = np.isnan(expected)
    assert np.array_equal(nan_actual, nan_expected), \
        f"NaN positions differ, e.g. at {np.argwhere(nan_actual != nan_expected)[:5].tolist()}"
    infinite = np.isinf(expected)
    assert np.array_equal(actual[infinite], expected[infinite]), "infinite values differ"
    finite = np.isfinite(expected)
    a = actual[finite]
    e = expected[finite]
    assert np.all(np.isfinite(a)), "finite reference values are not finite here"
    diff = np.abs(a - e)
    nonzero = e != 0
    max_rel = float((diff[nonzero] / np.abs(e[nonzero])).max()) if nonzero.any() else 0.0
    assert np.all(diff[~nonzero] == 0), "values differ where the reference is 0"
    assert max_rel <= RTOL, f"max relative difference {max_rel:.3e} > {RTOL}"
    entry = SUMMARY.setdefault(oracle, [0, 0, 0, 0.0])
    entry[0] += 1
    entry[1] += expected.size
    entry[2] += int(np.count_nonzero(expected[finite]))
    entry[3] = max(entry[3], max_rel)
