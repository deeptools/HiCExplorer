"""Comparator registry.

Adding a format means adding a module here that exposes
compare(path_a, path_b, cls, opts) -> Result and registering it below.
"""
from __future__ import annotations

from . import cool, h5, interval, npz, rendered, text
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
    # Files an external renderer draws from a source the tool writes, for
    # hicMergeDomains' graphviz trees; see comparators/rendered.py.
    "rendered": rendered,
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
    options = dict(opts or {})
    if cls == "E7":
        # The content is not compared, but a format that can check the
        # structure of what it is given (comparators/rendered.py) still does.
        structure = getattr(_REGISTRY[fmt], "check_structure", None)
        if structure is not None:
            checked = structure(path_a, path_b, options)
            if not checked.passed:
                return checked
        return Result(True, "E7", {"note": "deliberate deviation, not checked"})
    if fmt in _DEFAULT_SCHEMA:
        options.setdefault("schema", _DEFAULT_SCHEMA[fmt])
    return _REGISTRY[fmt].compare(path_a, path_b, cls, options)
