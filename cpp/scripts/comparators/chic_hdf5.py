"""Comparator for the cHi-C HDF5 files (PLAN.md 9.3, chic_hdf5).

The layout is written with h5py directly by lib/viewpoint.py and the cHi-C
tools: a group per matrix, below it a group per chromosome and a `genes` group,
a group per reference point with scalar and one dimensional datasets, and
attributes on the root. Nothing here knows those names; the comparator walks
whatever tree the reference file has, so the same module serves every cHi-C
tool.

Compared, in this order:

  1. the link tree: the set of paths, and for each whether it is a group or a
     dataset;
  2. hard links: paths that are the same object in the reference (the `genes`
     group links to the chromosome groups) must be the same object in the
     candidate, and paths that are different objects must stay different;
  3. attributes of every object: names, dtype, shape, and value;
  4. every dataset: dtype (including a string's encoding and whether it is of
     variable length), shape, chunk shape, compression filter and level,
     shuffle and fletcher32, then the values.

Integer, boolean and string values are compared exactly at every class.
Floating point values are compared at the requested class through
base.values_agree, with NaN equal to NaN and a stored zero required to stay
zero under ED; at class EN every float dataset and attribute is compared
against the envelope of the reference runs in opts["noise_paths"].
"""
from __future__ import annotations

try:
    import h5py
    import numpy as np
except ImportError as error:  # pragma: no cover
    raise SystemExit("the chic_hdf5 comparator needs h5py and numpy; run equiv.py with "
                     "the reference venv interpreter") from error

from . import noise
from .base import Result, values_agree


def _identity(obj):
    info = h5py.h5o.get_info(obj.id)
    return (info.fileno, info.addr)


def _walk(handle):
    """{path: (kind, identity)} over every link, following each object once."""
    tree = {}
    seen = set()

    def visit(group, prefix):
        for name in sorted(group.keys()):
            path = f"{prefix}/{name}" if prefix else name
            obj = group[name]
            identity = _identity(obj)
            kind = "group" if isinstance(obj, h5py.Group) else "dataset"
            tree[path] = (kind, identity)
            if kind == "group" and identity not in seen:
                seen.add(identity)
                visit(obj, path)

    seen.add(_identity(handle))
    visit(handle, "")
    return tree


def _canonical(tree):
    """{path: the first path, in sorted order, naming the same object}"""
    first = {}
    for path in sorted(tree):
        first.setdefault(tree[path][1], path)
    return {path: first[identity] for path, (_, identity) in tree.items()}


def _dtype_signature(dtype):
    string = h5py.check_string_dtype(dtype)
    if string is not None:
        return ("string", string.encoding, string.length)
    return ("numeric", dtype.str)


def _values(node):
    value = node[()]
    return np.atleast_1d(np.asarray(value))


def _compare_arrays(label, left, right, cls, diffs, stats):
    if left.shape != right.shape:
        diffs.append(f"{label}: shape {left.shape} vs {right.shape}")
        return
    if left.dtype.kind == "f" and right.dtype.kind == "f":
        flat_left = left.ravel()
        flat_right = right.ravel()
        stats["float_values"] += flat_left.size
        bad = 0
        for index in range(flat_left.size):
            a = float(flat_left[index])
            b = float(flat_right[index])
            if not values_agree(b, a, cls):
                if bad < 5:
                    diffs.append(f"{label}[{index}]: {a!r} vs {b!r}")
                bad += 1
            elif a != b and not (np.isnan(a) and np.isnan(b)) and a != 0.0:
                stats["max_relative"] = max(stats["max_relative"], abs(a - b) / abs(a))
        stats["float_values_differing"] += bad
        return
    if not np.array_equal(left, right):
        differing = np.flatnonzero(left.ravel() != right.ravel()) if left.dtype == right.dtype \
            else np.arange(min(left.size, right.size))
        position = int(differing[0]) if differing.size else 0
        diffs.append(f"{label}: values differ, first at [{position}]: "
                     f"{left.ravel()[position]!r} vs {right.ravel()[position]!r}")


def compare(path_a, path_b, cls, opts=None):
    opts = opts or {}
    diffs = []
    stats = {"groups": 0, "datasets": 0, "attributes": 0, "float_values": 0,
             "float_values_differing": 0, "max_relative": 0.0}
    noise_paths = list(opts.get("noise_paths") or [])
    if cls == "EN" and len(noise_paths) < 2:
        return Result(False, None, {}, ["class EN needs at least two reference runs "
                                        "(opts['noise_paths'])"])
    compare_cls = "E0" if cls == "EN" else cls
    en_values = {}

    with h5py.File(path_a, "r") as reference, h5py.File(path_b, "r") as candidate:
        tree_a = _walk(reference)
        tree_b = _walk(candidate)
        only_a = sorted(set(tree_a) - set(tree_b))
        only_b = sorted(set(tree_b) - set(tree_a))
        for path in only_a[:10]:
            diffs.append(f"{path}: only in the python output")
        for path in only_b[:10]:
            diffs.append(f"{path}: only in the cpp output")
        common = sorted(set(tree_a) & set(tree_b))
        for path in common:
            if tree_a[path][0] != tree_b[path][0]:
                diffs.append(f"{path}: a {tree_a[path][0]} in python, a {tree_b[path][0]} in cpp")
        canonical_a = _canonical({path: tree_a[path] for path in common})
        canonical_b = _canonical({path: tree_b[path] for path in common})
        for path in common:
            if canonical_a[path] != canonical_b[path]:
                diffs.append(f"{path}: is the object at {canonical_a[path]} in python but "
                             f"{canonical_b[path]} in cpp (hard link structure differs)")

        objects = [("", reference, candidate)] + [
            (path, reference[path], candidate[path]) for path in common
            if canonical_a[path] == path and tree_a[path][0] == tree_b[path][0]]
        for path, left, right in objects:
            label = path or "/"
            if isinstance(left, h5py.Group):
                stats["groups"] += 1
            names_a = sorted(left.attrs.keys())
            names_b = sorted(right.attrs.keys())
            if names_a != names_b:
                diffs.append(f"{label}: attributes {names_a} vs {names_b}")
            for name in sorted(set(names_a) & set(names_b)):
                stats["attributes"] += 1
                attribute_a = left.attrs.get_id(name)
                attribute_b = right.attrs.get_id(name)
                if _dtype_signature(attribute_a.dtype) != _dtype_signature(attribute_b.dtype):
                    diffs.append(f"{label}@{name}: dtype {_dtype_signature(attribute_a.dtype)} vs "
                                 f"{_dtype_signature(attribute_b.dtype)}")
                    continue
                if attribute_a.shape != attribute_b.shape:
                    diffs.append(f"{label}@{name}: shape {attribute_a.shape} vs {attribute_b.shape}")
                    continue
                _compare_arrays(f"{label}@{name}", np.atleast_1d(np.asarray(left.attrs[name])),
                                np.atleast_1d(np.asarray(right.attrs[name])), compare_cls,
                                diffs if cls != "EN" else [], stats)
                if cls == "EN" and np.asarray(left.attrs[name]).dtype.kind == "f":
                    en_values[f"{label}@{name}"] = None
                elif cls == "EN":
                    local = []
                    _compare_arrays(f"{label}@{name}", np.atleast_1d(np.asarray(left.attrs[name])),
                                    np.atleast_1d(np.asarray(right.attrs[name])), "E0", local, stats)
                    diffs.extend(local)
            if not isinstance(left, h5py.Dataset):
                continue
            stats["datasets"] += 1
            if _dtype_signature(left.dtype) != _dtype_signature(right.dtype):
                diffs.append(f"{label}: dtype {_dtype_signature(left.dtype)} vs "
                             f"{_dtype_signature(right.dtype)}")
                continue
            for property_name in ("shape", "chunks", "compression", "compression_opts",
                                  "shuffle", "fletcher32"):
                value_a = getattr(left, property_name)
                value_b = getattr(right, property_name)
                if value_a != value_b:
                    diffs.append(f"{label}: {property_name} {value_a!r} vs {value_b!r}")
            if left.shape != right.shape:
                continue
            values_a = _values(left)
            values_b = _values(right)
            if cls == "EN" and values_a.dtype.kind == "f":
                en_values[label] = None
                continue
            _compare_arrays(label, values_a, values_b, compare_cls, diffs, stats)

    metrics = dict(stats)
    if cls == "EN" and en_values and not diffs:
        runs = []
        for reference_path in noise_paths:
            with h5py.File(reference_path, "r") as handle:
                runs.append(np.concatenate([_read_float(handle, key) for key in en_values]))
        with h5py.File(path_b, "r") as handle:
            candidate_values = np.concatenate([_read_float(handle, key) for key in en_values])
        labels = []
        with h5py.File(path_a, "r") as handle:
            for key in en_values:
                labels.extend(f"{key}[{i}]" for i in range(_read_float(handle, key).size))
        passed, noise_metrics, noise_diffs = noise.compare(runs, candidate_values, labels)
        metrics.update(noise_metrics)
        diffs.extend(noise_diffs)
        if not passed and not noise_diffs:
            diffs.append("outside the reference envelope")
    if diffs:
        return Result(False, None, metrics, diffs[:20])
    return Result(True, cls, metrics)


def _read_float(handle, key):
    if "@" in key:
        path, name = key.split("@", 1)
        node = handle if path == "/" else handle[path]
        return np.atleast_1d(np.asarray(node.attrs[name], dtype=np.float64)).ravel()
    return np.atleast_1d(np.asarray(handle[key][()], dtype=np.float64)).ravel()
