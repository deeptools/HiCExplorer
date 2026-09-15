"""Comparator for the native HiCExplorer h5 format.

The node set is compared, then the CSR arrays and the interval lists are
decoded and compared at the requested class. Chunk shapes and blosc parameters
are deliberately not compared (cpp/PLAN.md section 2.5): a C++ writer will pick
different chunking than PyTables. The filter identity is compared, so a file
that PyTables could not read still fails.
"""
from __future__ import annotations

try:
    import numpy as np
except ImportError as error:  # pragma: no cover
    raise SystemExit(
        "the h5 comparator needs numpy and PyTables; run equiv.py with the "
        "reference venv interpreter named in cpp/AGENTS_CONTRACT.md") from error

from .base import Result, fail, values_agree

NODES = (
    "/matrix/data",
    "/matrix/indices",
    "/matrix/indptr",
    "/matrix/shape",
    "/intervals/chr_list",
    "/intervals/start_list",
    "/intervals/end_list",
    "/intervals/extra_list",
)
OPTIONAL_NODES = ("/nan_bins", "/correction_factors", "/distance_counts")


def _open(path):
    import tables

    return tables.open_file(path, "r")


def _read(handle, node):
    return handle.get_node(node).read()


def _has(handle, node):
    try:
        handle.get_node(node)
        return True
    except Exception:  # pylint: disable=W0718
        return False


def compare(path_a, path_b, cls, opts=None):
    opts = opts or {}
    handle_a = _open(path_a)
    handle_b = _open(path_b)
    try:
        diffs = []
        metrics = {}
        for node in NODES:
            if not _has(handle_a, node) or not _has(handle_b, node):
                diffs.append(f"{node}: present in "
                             f"{'py' if _has(handle_a, node) else 'cpp'} only")
        if diffs:
            return Result(False, None, metrics, diffs)
        for node in OPTIONAL_NODES:
            if _has(handle_a, node) != _has(handle_b, node):
                diffs.append(f"{node}: present in "
                             f"{'py' if _has(handle_a, node) else 'cpp'} only")

        for node in ("/matrix/shape", "/matrix/indptr", "/matrix/indices"):
            left = _read(handle_a, node)
            right = _read(handle_b, node)
            if left.shape != right.shape or not np.array_equal(left, right):
                diffs.append(f"{node}: differs "
                             f"({left.shape} vs {right.shape})")
        data_a = _read(handle_a, "/matrix/data")
        data_b = _read(handle_b, "/matrix/data")
        metrics["nnz_py"] = int(data_a.size)
        metrics["nnz_cpp"] = int(data_b.size)
        if data_a.shape != data_b.shape:
            diffs.append(f"/matrix/data: {data_a.shape} vs {data_b.shape}")
        elif data_a.dtype != data_b.dtype:
            diffs.append(f"/matrix/data: dtype {data_a.dtype} vs {data_b.dtype}")
        else:
            diffs += _compare_values("/matrix/data", data_a, data_b, cls, metrics)

        # "ignore_nodes": interval lists whose values a case declares out of
        # scope, with the reason in its notes. The node must still exist in
        # both files. Used by hicBuildMatrix --pairsFile against the BAM
        # route: a .pairs file carries no read length, so the per bin coverage
        # maximum of /intervals/extra_list cannot be computed.
        ignored = set(opts.get("ignore_nodes", []))
        for node in ("/intervals/chr_list", "/intervals/start_list",
                     "/intervals/end_list", "/intervals/extra_list"):
            if node in ignored:
                continue
            left = _read(handle_a, node)
            right = _read(handle_b, node)
            if left.shape != right.shape:
                diffs.append(f"{node}: {left.shape} vs {right.shape}")
            elif left.dtype.kind == "f":
                diffs += _compare_values(node, left, right, cls, metrics)
            elif not np.array_equal(left, right):
                diffs.append(f"{node}: differs")

        for node in OPTIONAL_NODES:
            if _has(handle_a, node) and _has(handle_b, node):
                left = _read(handle_a, node)
                right = _read(handle_b, node)
                if left.shape != right.shape:
                    diffs.append(f"{node}: {left.shape} vs {right.shape}")
                elif left.dtype.kind == "f":
                    diffs += _compare_values(node, left, right, cls, metrics)
                elif not np.array_equal(left, right):
                    diffs.append(f"{node}: differs")
    finally:
        handle_a.close()
        handle_b.close()

    if diffs:
        return Result(False, None, metrics, diffs[:20])
    return Result(True, cls, metrics)


def _compare_values(node, left, right, cls, metrics):
    diffs = []
    if cls in ("E0", "E1", "E2"):
        equal = np.array_equal(left, right, equal_nan=left.dtype.kind == "f")
        if left.dtype.kind == "f" and equal:
            equal = np.array_equal(np.signbit(left), np.signbit(right))
        if not equal:
            differing = np.flatnonzero(_differing_mask(left, right))
            metrics[f"differing{node}"] = int(differing.size)
            for position in differing[:20]:
                diffs.append(f"{node}[{position}]: {left[position]!r} vs "
                             f"{right[position]!r}")
        return diffs
    bad = []
    for position in np.flatnonzero(_differing_mask(left, right)):
        if not values_agree(float(left[position]), float(right[position]), cls):
            bad.append(position)
    metrics[f"differing{node}"] = len(bad)
    if bad:
        with np.errstate(invalid="ignore"):
            metrics[f"max_abs{node}"] = float(
                np.nanmax(np.abs(left.astype(float) - right.astype(float))))
        for position in bad[:20]:
            diffs.append(f"{node}[{position}]: {left[position]!r} vs "
                         f"{right[position]!r}")
    return diffs


def _differing_mask(left, right):
    if left.dtype.kind == "f":
        both_nan = np.isnan(left) & np.isnan(right)
        return (left != right) & ~both_nan
    return left != right
