"""Comparator interface.

A comparator answers one question: do these two output files agree at the
requested equivalence class? Classes are defined in cpp/PLAN.md section 5.1:

    E0  byte identical
    E1  structurally identical HDF5, every dataset bit identical
    E2  same sparsity pattern, every stored value bit identical
    E3  same pattern, max |a-b| / max(1, |b|) <= 1e-12
    E4  pattern equal to within 0.1 % of the nonzeros, relative 1e-6, r >= 1-1e-9
    E5  Jaccard of called intervals >= 0.99
    E6  image RMS <= 5
    E7  deliberate deviation, not checked automatically

Every comparator module exposes

    compare(path_a, path_b, cls, opts) -> Result

where path_a is the Python output and path_b is the C++ output.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List

CLASSES = ("E0", "E1", "E2", "E3", "E4", "E5", "E6", "E7")


@dataclass
class Result:
    passed: bool
    class_met: str | None = None
    metrics: Dict[str, Any] = field(default_factory=dict)
    diffs: List[str] = field(default_factory=list)

    def to_json(self):
        return {
            "passed": self.passed,
            "class_met": self.class_met,
            "metrics": self.metrics,
            "diffs": self.diffs[:20],
        }


def fail(message, **metrics):
    return Result(passed=False, class_met=None, metrics=metrics, diffs=[message])


def float_tolerance(cls):
    """Relative tolerance and its denominator floor for the float classes."""
    if cls == "E3":
        return 1e-12, 1.0
    if cls == "E4":
        return 1e-6, 1e-6
    return 0.0, 1.0


def values_agree(a, b, cls):
    """Compare two floats at the given class, treating NaN as equal to NaN."""
    if cls in ("E0", "E1", "E2"):
        if math.isnan(a) and math.isnan(b):
            return True
        # Bit identical, so the sign of zero matters as well.
        return math.copysign(1.0, a) == math.copysign(1.0, b) and a == b
    if math.isnan(a) and math.isnan(b):
        return True
    if math.isinf(a) or math.isinf(b):
        return a == b
    tolerance, floor = float_tolerance(cls)
    return abs(a - b) / max(floor, abs(b)) <= tolerance
