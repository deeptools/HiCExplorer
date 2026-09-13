#!/usr/bin/env python3
"""Class E5 for the cHi-C background model: do the significance calls change?

The fitted size and prob of chicViewpointBackgroundModel are compared at class
EN, because the reference is not reproducible against itself (PLAN.md 5.7).
What a user sees downstream of them is a set of significant interactions, and
that set is class E5 (Jaccard >= 0.99). This script measures it:

    1. the background model is built three times: by the Python, by the C++,
       and by the Python again with its distributions presented to the fit in
       a different order (noise_runner.py, shim fit_nbinom_order), which is a
       second draw of the reference's own noise;
    2. each model is handed to the *same* Python chicViewpoint and the *same*
       Python chicSignificantInteractions, so any difference in the calls comes
       from the model alone;
    3. the significant interactions are read back as a set of
       (matrix, reference point group, start, end) and compared.

It prints the Jaccard index of C++ against Python and, as the control, of the
second Python run against the first, together with how many p-values of the
interaction files differ at all.

    chic_background_jaccard.py [--cpp-bin DIR] [--keep DIR] [--pvalue P]
                               [--xfold X] [--range UP DOWN]
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
DATA = REPO_ROOT / "hicexplorer" / "test" / "test_data" / "cHi-C"
VENV = Path("/tmp/claude-1000167729/-home-mh-hannover-local-wolffjoa/"
            "3de46a4d-41bb-4fb9-aea9-26cceef9ad01/scratchpad/hicx-venv/bin/python")
PYTHON = Path(os.environ.get("HICX_PY_PYTHON", VENV if VENV.exists() else sys.executable))
MATRICES = [str(DATA / "FL-E13-5_chr1.cool"), str(DATA / "MB-E10-5_chr1.cool")]
REFERENCE_POINTS = str(DATA / "referencePoints.bed")


def run(argv, cwd, env=None):
    result = subprocess.run([str(part) for part in argv], cwd=cwd, env=env,
                            capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stdout + result.stderr)
        raise SystemExit(f"failed: {' '.join(str(part) for part in argv)}")


def python_env(extra=None):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    env.update(extra or {})
    return env


def significant_set(path):
    calls = set()
    with h5py.File(path, "r") as handle:
        def visit(name, obj):
            if isinstance(obj, h5py.Group) and "start_list" in obj and "end_list" in obj:
                if name.split("/")[1] == "genes":
                    return
                for start, end in zip(obj["start_list"][()], obj["end_list"][()]):
                    calls.add((name, int(start), int(end)))
        handle.visititems(visit)
    return calls


def pvalues(path):
    values = {}
    with h5py.File(path, "r") as handle:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset) and name.endswith("/pvalue") and "/genes/" not in name:
                values[name] = obj[()]
        handle.visititems(visit)
    return values


def jaccard(left, right):
    union = left | right
    return 1.0 if not union else len(left & right) / len(union)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cpp-bin", default=str(REPO_ROOT / "cpp" / "build" / "tools"))
    parser.add_argument("--keep", default=None)
    parser.add_argument("--pvalue", default="0.2")
    parser.add_argument("--xfold", default="1.5")
    parser.add_argument("--range", nargs=2, default=["200000", "200000"])
    options = parser.parse_args()

    work = Path(options.keep) if options.keep else Path(tempfile.mkdtemp(prefix="chic_e5_"))
    work.mkdir(parents=True, exist_ok=True)
    models = {
        "python": work / "background_python.txt",
        "cpp": work / "background_cpp.txt",
        "python_reordered": work / "background_python_reordered.txt",
    }
    common = ["--matrices", *MATRICES, "--referencePoints", REFERENCE_POINTS, "-t", "1", "-o"]
    run([PYTHON, REPO_ROOT / "bin" / "chicViewpointBackgroundModel", *common, models["python"]],
        work, python_env())
    run([Path(options.cpp_bin) / "chicViewpointBackgroundModel", *common, models["cpp"]], work)
    run([PYTHON, SCRIPT_DIR / "noise_runner.py", REPO_ROOT / "bin" / "chicViewpointBackgroundModel",
         *common, models["python_reordered"]], work,
        python_env({"HICX_NOISE_SHIM": "fit_nbinom_order", "HICX_NOISE_RUN": "1"}))

    calls = {}
    interaction_pvalues = {}
    for name, model in models.items():
        interactions = work / f"interactions_{name}.hdf5"
        run([PYTHON, REPO_ROOT / "bin" / "chicViewpoint", "--matrices", *MATRICES,
             "--referencePoints", REFERENCE_POINTS, "--backgroundModelFile", model,
             "--range", *options.range, "-o", interactions, "-t", "1"], work, python_env())
        significant = work / f"significant_{name}.hdf5"
        target = work / f"target_{name}.hdf5"
        run([PYTHON, REPO_ROOT / "bin" / "chicSignificantInteractions",
             "--interactionFile", interactions, "--backgroundModelFile", model,
             "--range", *options.range, "--pValue", options.pvalue,
             "--xFoldBackground", options.xfold,
             "--outFileNameSignificant", significant, "--outFileNameTarget", target,
             "-t", "1", "--combinationMode", "dual"], work, python_env())
        calls[name] = significant_set(significant)
        interaction_pvalues[name] = pvalues(interactions)

    def pvalue_difference(left, right):
        differing = 0
        total = 0
        for key, values in interaction_pvalues[left].items():
            other = interaction_pvalues[right][key]
            total += values.size
            differing += int(np.count_nonzero(values != other))
        return differing, total

    for candidate, label in (("cpp", "C++ model"), ("python_reordered", "reordered Python model")):
        differing, total = pvalue_difference("python", candidate)
        left = calls["python"]
        right = calls[candidate]
        print(f"{label}: {len(right)} significant interactions against {len(left)} from the "
              f"Python model; Jaccard {jaccard(left, right):.6f} "
              f"(only Python {len(left - right)}, only {label} {len(right - left)}); "
              f"interaction file p-values differing: {differing} of {total}")
    if not options.keep:
        shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
