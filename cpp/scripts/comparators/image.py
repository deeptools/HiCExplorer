"""Figures, class E6: the same pixel dimensions and an RMS difference of at
most 5 on the 0 to 255 scale.

The RMS is matplotlib.testing.compare's, the measure HiCExplorer's own figure
tests use through compare_images: both images are loaded as 8 bit RGB, and the
root mean square of the per channel differences is taken. Images of different
dimensions fail outright. The tolerance is the class definition (base.py)
and a case cannot raise it; a case may lower it with opts["tolerance"].

The plotting tools draw through the same matplotlib calls from the same data
(cpp/PLAN.md tier 7, option (a)), so most figures reach RMS 0; the metric is
recorded either way. E0 is a byte comparison. Only raster formats matplotlib
can load without an external converter are accepted (png).
"""
from __future__ import annotations

from pathlib import Path

from .base import Result, fail

E6_TOLERANCE = 5.0


def _rms(path_a, path_b):
    import numpy as np
    from matplotlib.testing.compare import _load_image, calculate_rms

    a = _load_image(str(path_a))
    b = _load_image(str(path_b))
    if a.shape != b.shape:
        return None, a.shape, b.shape
    return float(calculate_rms(a.astype(np.int16), b.astype(np.int16))), a.shape, b.shape


def check_structure(path_a, path_b, opts=None):
    for path in (path_a, path_b):
        with open(path, "rb") as handle:
            if not handle.read(8).startswith(b"\x89PNG\r\n\x1a\n"):
                return fail(f"{path}: not a PNG file")
    return Result(True, "E7", {})


def compare(path_a, path_b, cls, opts=None):
    opts = dict(opts or {})
    if Path(path_a).suffix.lower() != ".png" or Path(path_b).suffix.lower() != ".png":
        return fail("the image comparator reads png files only")
    structure = check_structure(path_a, path_b)
    if not structure.passed:
        return structure
    with open(path_a, "rb") as handle:
        a = handle.read()
    with open(path_b, "rb") as handle:
        b = handle.read()
    if cls == "E0":
        if a == b:
            return Result(True, "E0", {"bytes": len(a)})
        return fail("figures differ byte for byte", bytes_py=len(a), bytes_cpp=len(b))
    if cls != "E6":
        return fail(f"class {cls} is not defined for figures; use E6 or E0")
    tolerance = min(float(opts.get("tolerance", E6_TOLERANCE)), E6_TOLERANCE)
    rms, shape_a, shape_b = _rms(path_a, path_b)
    if rms is None:
        return fail(f"figure dimensions differ: python {shape_a[1]}x{shape_a[0]}, "
                    f"cpp {shape_b[1]}x{shape_b[0]}")
    metrics = {"rms": rms, "width": shape_a[1], "height": shape_a[0],
               "tolerance": tolerance, "byte_identical": a == b}
    if rms > tolerance:
        return fail(f"figure RMS {rms:.3f} exceeds {tolerance}", **metrics)
    return Result(True, "E0" if a == b else "E6", metrics)
