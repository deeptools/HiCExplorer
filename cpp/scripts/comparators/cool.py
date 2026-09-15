"""Comparator for cool and mcool files.

E1 compares the HDF5 structure: the object tree, dtypes, shapes, chunking, the
filter pipeline and the decoded bytes of every dataset, plus the attributes
after normalising the four fields that legitimately differ between two
implementations writing the same matrix.

E2 and above ignore the layout and compare the pixel table as a COO set, so a
file written with a different chunk size still passes as long as the matrix is
the same.
"""
from __future__ import annotations

try:
    import numpy as np
except ImportError as error:  # pragma: no cover
    raise SystemExit(
        "the cool comparator needs numpy and h5py; run equiv.py with the "
        "reference venv interpreter named in cpp/AGENTS_CONTRACT.md") from error

from .base import Result, fail

# Attributes that carry the identity of the writing implementation, not the
# content of the matrix.
NORMALISED_ATTRS = {
    "creation-date",
    "generated-by",
    "generated-by-cooler-lib",
    "tool-url",
}


def _coolers(handle):
    """Every cooler group in the file, '/' for a plain .cool."""
    import h5py

    groups = []

    def visit(name, obj):
        if isinstance(obj, h5py.Group) and {"chroms", "bins", "pixels", "indexes"} <= set(obj.keys()):
            groups.append("/" + name)

    if {"chroms", "bins", "pixels", "indexes"} <= set(handle.keys()):
        groups.append("/")
    handle.visititems(visit)
    return sorted(set(groups))


def _datasets(group):
    import h5py

    found = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            found[name] = obj

    group.visititems(visit)
    return found


def compare(path_a, path_b, cls, opts=None):
    import h5py

    opts = opts or {}
    # "ignore_attributes": attributes a case declares out of scope, with the
    # reason in its notes. Used where the reference is another writer (cooler
    # cload pairs against hicBuildMatrix --pairsFile), whose provenance and
    # metadata attributes differ by design; the pixels, bins and chromosomes
    # are still compared.
    ignored = frozenset(opts.get("ignore_attributes", []))
    with h5py.File(path_a, "r") as file_a, h5py.File(path_b, "r") as file_b:
        groups_a = _coolers(file_a)
        groups_b = _coolers(file_b)
        if groups_a != groups_b:
            return fail(f"different cooler groups: {groups_a} vs {groups_b}")
        diffs = []
        metrics = {}
        if "/" not in groups_a:
            # An mcool: the coolers live in groups and the file root carries
            # the provenance hicmatrix writes there for the first resolution
            # only (hicmatrix/lib/cool.py:422-426). Nothing else compares it.
            diffs += _compare_attrs(file_a, file_b, "/", ignored)
            if sorted(file_a.keys()) != sorted(file_b.keys()):
                diffs.append(f"root children {sorted(file_a.keys())} vs "
                             f"{sorted(file_b.keys())}")
        for group_path in groups_a:
            group_a = file_a[group_path]
            group_b = file_b[group_path]
            if cls == "E1":
                diffs += _compare_structure(group_a, group_b, group_path)
            diffs += _compare_attrs(group_a, group_b, group_path, ignored)
            pixel_diffs, pixel_metrics = _compare_pixels(group_a, group_b, cls,
                                                         group_path)
            diffs += pixel_diffs
            metrics.update(pixel_metrics)
            diffs += _compare_bins(group_a, group_b, cls, group_path)
            if len(diffs) >= 20:
                break
    if diffs:
        return Result(False, None, metrics, diffs)
    return Result(True, cls, metrics)


def _compare_structure(group_a, group_b, prefix):
    diffs = []
    datasets_a = _datasets(group_a)
    datasets_b = _datasets(group_b)
    only_a = sorted(set(datasets_a) - set(datasets_b))
    only_b = sorted(set(datasets_b) - set(datasets_a))
    if only_a:
        diffs.append(f"{prefix}: datasets only in the Python output: {only_a}")
    if only_b:
        diffs.append(f"{prefix}: datasets only in the C++ output: {only_b}")
    for name in sorted(set(datasets_a) & set(datasets_b)):
        left = datasets_a[name]
        right = datasets_b[name]
        if left.dtype != right.dtype:
            diffs.append(f"{prefix}/{name}: dtype {left.dtype} vs {right.dtype}")
        if left.shape != right.shape:
            diffs.append(f"{prefix}/{name}: shape {left.shape} vs {right.shape}")
            continue
        if left.chunks != right.chunks:
            diffs.append(f"{prefix}/{name}: chunks {left.chunks} vs {right.chunks}")
        filters_left = _filters(left)
        filters_right = _filters(right)
        if filters_left != filters_right:
            diffs.append(f"{prefix}/{name}: filters {filters_left} vs {filters_right}")
        if not _arrays_identical(left[:], right[:]):
            diffs.append(f"{prefix}/{name}: dataset content differs")
    return diffs


def _filters(dataset):
    plist = dataset.id.get_create_plist()
    return [plist.get_filter(i)[:3] for i in range(plist.get_nfilters())]


def _arrays_identical(a, b):
    if a.dtype != b.dtype or a.shape != b.shape:
        return False
    if a.dtype.kind in "fc":
        return np.array_equal(a, b, equal_nan=True) and np.array_equal(
            np.signbit(a), np.signbit(b))
    return np.array_equal(a, b)


def _compare_attrs(group_a, group_b, prefix, ignored=frozenset()):
    diffs = []
    keys_a = set(group_a.attrs) - NORMALISED_ATTRS - set(ignored)
    keys_b = set(group_b.attrs) - NORMALISED_ATTRS - set(ignored)
    for key in sorted(keys_a - keys_b):
        diffs.append(f"{prefix}: attribute {key} only in the Python output")
    for key in sorted(keys_b - keys_a):
        diffs.append(f"{prefix}: attribute {key} only in the C++ output")
    for key in sorted(keys_a & keys_b):
        left = group_a.attrs[key]
        right = group_b.attrs[key]
        if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
            same = np.array_equal(left, right)
        else:
            same = left == right
        if not same:
            diffs.append(f"{prefix}: attribute {key}: {left!r} vs {right!r}")
    return diffs


def _compare_pixels(group_a, group_b, cls, prefix):
    from .base import values_agree

    metrics = {}
    diffs = []
    bin1_a = group_a["pixels/bin1_id"][:]
    bin2_a = group_a["pixels/bin2_id"][:]
    count_a = group_a["pixels/count"][:]
    bin1_b = group_b["pixels/bin1_id"][:]
    bin2_b = group_b["pixels/bin2_id"][:]
    count_b = group_b["pixels/count"][:]
    metrics[f"nnz_py{prefix}"] = int(len(count_a))
    metrics[f"nnz_cpp{prefix}"] = int(len(count_b))

    if len(bin1_a) != len(bin1_b) or not (
            np.array_equal(bin1_a, bin1_b) and np.array_equal(bin2_a, bin2_b)):
        diffs.append(f"{prefix}: the pixel coordinates differ "
                     f"({len(bin1_a)} vs {len(bin1_b)} entries)")
        return diffs, metrics

    if count_a.dtype.kind in "fc" or count_b.dtype.kind in "fc":
        differing = ~np.isclose(count_a, count_b, rtol=0, atol=0, equal_nan=True)
    else:
        differing = count_a != count_b
    if differing.any():
        if cls in ("E0", "E1", "E2"):
            positions = np.flatnonzero(differing)[:20]
            for position in positions:
                diffs.append(f"{prefix}: pixel ({bin1_a[position]}, "
                             f"{bin2_a[position]}): {count_a[position]!r} vs "
                             f"{count_b[position]!r}")
            metrics[f"differing_pixels{prefix}"] = int(differing.sum())
        else:
            bad = [position for position in np.flatnonzero(differing)
                   if not values_agree(float(count_a[position]),
                                       float(count_b[position]), cls)]
            metrics[f"differing_pixels{prefix}"] = len(bad)
            for position in bad[:20]:
                diffs.append(f"{prefix}: pixel ({bin1_a[position]}, "
                             f"{bin2_a[position]}): {count_a[position]!r} vs "
                             f"{count_b[position]!r}")
    if len(count_a) and len(count_b) and len(count_a) == len(count_b):
        with np.errstate(invalid="ignore"):
            difference = np.abs(count_a.astype(float) - count_b.astype(float))
        metrics[f"max_abs{prefix}"] = float(np.nanmax(difference)) if len(difference) else 0.0
    return diffs, metrics


def _compare_bins(group_a, group_b, cls, prefix):
    diffs = []
    columns_a = sorted(group_a["bins"].keys())
    columns_b = sorted(group_b["bins"].keys())
    if columns_a != columns_b:
        diffs.append(f"{prefix}: bin columns {columns_a} vs {columns_b}")
        return diffs
    for column in columns_a:
        left = group_a["bins/" + column][:]
        right = group_b["bins/" + column][:]
        if not _arrays_identical(left, right):
            diffs.append(f"{prefix}: bins/{column} differs")
    for name in ("chroms/name", "chroms/length"):
        if not _arrays_identical(group_a[name][:], group_b[name][:]):
            diffs.append(f"{prefix}: {name} differs")
    return diffs
