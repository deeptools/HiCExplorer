"""Comparator interface.

A comparator answers one question: do these two output files agree at the
requested equivalence class? Classes are defined in cpp/PLAN.md section 5.1:

    ED  the project acceptance gate: every item agrees to three significant
        digits, that is |a-b| / |b| <= 1e-3, with a stored zero required to
        stay a stored zero
    E0  byte identical
    E1  structurally identical HDF5, every dataset bit identical
    E2  same sparsity pattern, every stored value bit identical
    E3  same pattern, max |a-b| / max(1, |b|) <= 1e-12
    E4  pattern equal to within 0.1 % of the nonzeros, relative 1e-6, r >= 1-1e-9
    E5  Jaccard of called intervals >= 0.99
    E6  image RMS <= 5
    E7  deliberate deviation, not checked automatically

ED is what a port must meet. E0 to E4 are stricter and are recorded when they
are met, because a stricter result is a better regression signal and several
tools reach E0 for free, but they are not required. Declaring a case stricter
than ED is a statement that the tool is expected to achieve it.

Every comparator module exposes

    compare(path_a, path_b, cls, opts) -> Result

where path_a is the Python output and path_b is the C++ output.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List

CLASSES = ("ED", "E0", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "EN")

# The project acceptance gate, set 2026-09-01: three significant digits per
# item. Expressed as a relative tolerance because the corpus spans many orders
# of magnitude (corrected matrices hold values from 7.9e-06 to 0.09, raw counts
# reach 10^8), so a fixed number of decimal places would be meaningless at one
# end and unreachable at the other.
ACCEPTANCE_SIGNIFICANT_DIGITS = 3
ACCEPTANCE_RTOL = 10.0 ** (-ACCEPTANCE_SIGNIFICANT_DIGITS)


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
    if cls == "ED":
        # Pure relative, no floor: three significant digits means the same
        # thing at 1e-06 as at 1e+08. The zero case is handled separately in
        # values_agree, since a floor would silently accept a value appearing
        # where Python stored an exact zero.
        return ACCEPTANCE_RTOL, 0.0
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
    if floor == 0.0:
        # Pure relative comparison. An exact zero in the reference has no
        # significant digits to agree with, so it has to be matched exactly;
        # otherwise a spurious nonzero would pass by dividing by nothing.
        if b == 0.0:
            return a == 0.0
        return abs(a - b) / abs(b) <= tolerance
    return abs(a - b) / max(floor, abs(b)) <= tolerance
