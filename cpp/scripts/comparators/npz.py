"""Comparator for the scipy sparse .npz that hicAverageRegions writes.

cpp/PLAN.md 9.3 specifies it as `scipy.sparse.load_npz` on both files followed
by the section 5.5 rules, and that is what this does.

A byte comparison would be the wrong tool even though numpy's npz is otherwise
reproducible (it stamps every ZIP entry with the fixed 1980-01-01 timestamp
rather than the current time). The deflate streams come from two different zlib
builds, the interpreter's and the one in the conda prefix, and two zlib builds
do not emit the same bytes for the same input. So the archive is opened and the
arrays are compared, which is the same decision PLAN.md 2.5 takes for h5.

The container is still compared, not only the numbers: the entry names and
their order, the format tag, the shape and the dtype of every array. A port
that wrote a coo instead of a csr, or float32 instead of float64, fails here
even when the dense matrix it represents is the same.
"""
from __future__ import annotations

from .base import Result, fail, values_agree


def _load(path):
    import numpy as np
    import zipfile

    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        arrays = {}
        for name in names:
            with archive.open(name) as handle:
                arrays[name] = np.load(handle, allow_pickle=False)
    return names, arrays


def compare(path_a, path_b, cls, opts=None):
    import numpy as np

    with open(path_a, "rb") as handle:
        raw_a = handle.read()
    with open(path_b, "rb") as handle:
        raw_b = handle.read()

    try:
        names_a, arrays_a = _load(path_a)
        names_b, arrays_b = _load(path_b)
    except Exception as error:  # noqa: BLE001 - report, do not crash the run
        return fail(f"cannot read npz: {error}")

    metrics = {"entries_py": names_a, "entries_cpp": names_b}
    if names_a != names_b:
        return fail("npz entry names or order differ", f"{names_a} vs {names_b}",
                    **metrics)

    diffs = []
    # E0 means byte identical, and an npz written by two different zlib builds
    # is not, so the best a matching archive can claim is E2: same objects,
    # same dtypes, every stored value bit identical.
    class_met = "E0" if raw_a == raw_b else "E2"
    for name in names_a:
        left = arrays_a[name]
        right = arrays_b[name]
        if left.dtype != right.dtype:
            diffs.append(f"{name}: dtype {left.dtype} vs {right.dtype}")
            continue
        if left.shape != right.shape:
            diffs.append(f"{name}: shape {left.shape} vs {right.shape}")
            continue
        if left.dtype.kind in "SU" or left.dtype.kind in "iub":
            if not np.array_equal(left, right):
                diffs.append(f"{name}: values differ")
            continue
        if np.array_equal(left, right):
            continue
        # Float array: fall back to the case's class.
        class_met = cls if cls not in ("E0", "E1", "E2") else None
        if cls in ("E0", "E1", "E2"):
            differing = int(np.count_nonzero(left != right))
            diffs.append(f"{name}: {differing} of {left.size} values differ, "
                         f"max abs {float(np.max(np.abs(left - right)))}")
            continue
        bad = 0
        for index in range(left.size):
            if not values_agree(float(left.flat[index]), float(right.flat[index]), cls):
                if bad < 10:
                    diffs.append(f"{name}[{index}]: {left.flat[index]!r} vs "
                                 f"{right.flat[index]!r}")
                bad += 1
        if bad:
            diffs.append(f"{name}: {bad} of {left.size} values outside {cls}")

    if "data.npy" in arrays_a:
        metrics["nnz_py"] = int(arrays_a["data.npy"].size)
        metrics["nnz_cpp"] = int(arrays_b["data.npy"].size)
        if arrays_a["data.npy"].size and arrays_a["data.npy"].dtype.kind == "f":
            metrics["max_abs"] = float(
                np.max(np.abs(arrays_a["data.npy"] - arrays_b["data.npy"])))

    if diffs:
        return Result(False, None, metrics, diffs)
    return Result(True, class_met or cls, metrics)
