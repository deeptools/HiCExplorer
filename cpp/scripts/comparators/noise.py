"""Class EN: the reference is not reproducible against itself (PLAN.md 5.7).

Given N runs of the Python reference and one C++ output, per element:

    S       = the largest relative difference between any two reference runs,
              |a - b| / |b| over ordered pairs
    median  = the elementwise median of the reference runs
    passes  = |cpp - median| / |median| <= max(2 * S, 1e-12)

A reference value that is exactly zero has no relative scale. When the runs
agree on it (S = 0) the C++ value must be zero too; when some runs are zero
and others are not, S is infinite and the element cannot fail.

This module also measures the same rule for the reference against itself: the
last reference run is taken out, the envelope is rebuilt from the remaining
runs, and the held-out run is tested exactly as the C++ output is. That
number is reported next to the C++ one, because an elementwise envelope over
a handful of runs is a sample, and on a heavy tailed reference a fresh run of
the reference itself can land outside it. A failure rate of the C++ output
that matches the held-out reference run's is the signature of an oracle
whose noise the envelope under-samples, not of a port that is wrong; the gate
itself is not relaxed for it.
"""
from __future__ import annotations

import numpy as np

FLOOR = 1e-12


def _pairwise_spread(runs):
    """S per element; runs has shape (N, k)."""
    count = runs.shape[0]
    spread = np.zeros(runs.shape[1])
    with np.errstate(divide="ignore", invalid="ignore"):
        for a in range(count):
            for b in range(count):
                if a == b:
                    continue
                left = runs[a]
                right = runs[b]
                both_nan = np.isnan(left) & np.isnan(right)
                one_nan = np.isnan(left) ^ np.isnan(right)
                equal = (left == right) | both_nan
                relative = np.abs(left - right) / np.abs(right)
                relative = np.where(equal, 0.0, relative)
                relative = np.where(one_nan, np.inf, relative)
                relative = np.where((right == 0.0) & ~equal, np.inf, relative)
                spread = np.maximum(spread, relative)
    return spread


def _evaluate(runs, candidate):
    spread = _pairwise_spread(runs)
    median = np.median(runs, axis=0)
    tolerance = np.maximum(2.0 * spread, FLOOR)
    with np.errstate(divide="ignore", invalid="ignore"):
        deviation = np.abs(candidate - median) / np.abs(median)
    median_nan = np.isnan(median)
    deviation = np.where(median_nan, np.where(np.isnan(candidate), 0.0, np.inf), deviation)
    deviation = np.where(~median_nan & (median == 0.0),
                         np.where(candidate == 0.0, 0.0, np.inf), deviation)
    deviation = np.where(~median_nan & (median != 0.0) & (candidate == median), 0.0, deviation)
    passes = (deviation <= tolerance) | np.isinf(spread)
    return passes, spread, deviation, tolerance, median


def compare(reference_runs, candidate, labels=None):
    """reference_runs: sequence of 1-D float arrays, one per Python run.

    Returns (passed, metrics, diffs).
    """
    runs = np.vstack([np.asarray(run, dtype=np.float64) for run in reference_runs])
    candidate = np.asarray(candidate, dtype=np.float64)
    if runs.shape[0] < 2:
        return False, {}, ["class EN needs at least two reference runs"]
    if runs.shape[1] != candidate.shape[0]:
        return False, {}, [f"{runs.shape[1]} reference values against "
                           f"{candidate.shape[0]} in the candidate"]
    passes, spread, deviation, tolerance, median = _evaluate(runs, candidate)
    finite_spread = spread[np.isfinite(spread)]
    metrics = {
        "reference_runs": int(runs.shape[0]),
        "values": int(candidate.shape[0]),
        "within_envelope": int(np.count_nonzero(passes)),
        "outside_envelope": int(np.count_nonzero(~passes)),
        "envelope_S_median": float(np.median(finite_spread)) if finite_spread.size else 0.0,
        "envelope_S_max": float(np.max(finite_spread)) if finite_spread.size else 0.0,
        "deviation_from_median_max": float(np.max(np.where(np.isfinite(deviation), deviation, 0.0)))
        if deviation.size else 0.0,
    }
    if runs.shape[0] >= 3:
        held_passes, _, _, _, _ = _evaluate(runs[:-1], runs[-1])
        metrics["held_out_reference_run_outside_envelope"] = int(np.count_nonzero(~held_passes))
        metrics["held_out_reference_run_envelope_runs"] = int(runs.shape[0] - 1)
    diffs = []
    for index in np.flatnonzero(~passes)[:20]:
        label = labels[index] if labels is not None else str(index)
        diffs.append(f"{label}: cpp {candidate[index]!r}, reference median {median[index]!r}, "
                     f"deviation {deviation[index]:.3g} > tolerance {tolerance[index]:.3g} "
                     f"(S {spread[index]:.3g})")
    return bool(np.all(passes)), metrics, diffs
