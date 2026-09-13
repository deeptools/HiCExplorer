"""Comparator registry.

Adding a format means adding a module here that exposes
compare(path_a, path_b, cls, opts) -> Result and registering it below.
"""
from __future__ import annotations

from . import chic_background, chic_hdf5, cool, h5, interval, npz, text
from .base import CLASSES, Result, fail

_REGISTRY = {
    "text": text,
    "plain": text,
    "bed": text,
    "bedgraph": text,
    "bedpe": text,
    "tsv": text,
    "cool": cool,
    "mcool": cool,
    "h5": h5,
    # Set agreement over called regions, for a tool whose reference breaks a
    # tie non-reproducibly; see comparators/interval.py.
    "interval": interval,
    # scipy.sparse .npz, for hicAverageRegions; see comparators/npz.py.
    "npz": npz,
    # The cHi-C HDF5 layout written with h5py by lib/viewpoint.py and the
    # cHi-C tools; see comparators/chic_hdf5.py.
    "chic_hdf5": chic_hdf5,
    "hdf5-chic": chic_hdf5,
    # chicViewpointBackgroundModel; the fit columns are left to
    # validators/chic_background_model.py.
    "chic_background_model": chic_background,
}

# Default text schema per format name, so a case can just say format: bedgraph.
_DEFAULT_SCHEMA = {
    "text": "plain",
    "plain": "plain",
    "bed": "bed",
    "bedgraph": "bedgraph",
    "bedpe": "bedpe",
    "tsv": "tsv-header",
}


def formats():
    return sorted(_REGISTRY)


def compare(fmt, path_a, path_b, cls, opts=None):
    if fmt not in _REGISTRY:
        return fail(f"no comparator for format {fmt!r}; known: {formats()}")
    if cls not in CLASSES:
        return fail(f"unknown equivalence class {cls!r}")
    if cls == "E7":
        return Result(True, "E7", {"note": "deliberate deviation, not checked"})
    options = dict(opts or {})
    if fmt in _DEFAULT_SCHEMA:
        options.setdefault("schema", _DEFAULT_SCHEMA[fmt])
    return _REGISTRY[fmt].compare(path_a, path_b, cls, options)
