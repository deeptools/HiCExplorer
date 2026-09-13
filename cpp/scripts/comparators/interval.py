"""Set agreement comparator for tools whose output is a set of called regions.

Class E5 of cpp/PLAN.md 5.1: the Jaccard index of the called intervals must be
at least 0.99. The lines are compared as whole records, so a call that differs
in any column counts as a disagreement.

This exists for hicMergeLoops. Its reference breaks a tie between two equally
wide overlapping loop anchors by CPython's set iteration order, which depends
on the insertion order into a set built by intervaltree's tree traversal and is
therefore not reproducible without reimplementing both. The port reproduces the
part that is reproducible, the hash and the table size, and the residue is a
handful of lines out of about 16,000. Comparing that at E0 would either fail a
correct port or force the tolerance into the case file, and PLAN.md 9.5 forbids
a per-case tolerance; a class is the right place for it.

E0 through E4 fall back to the byte comparison, so declaring a case at E0 still
means byte identical.
"""
from __future__ import annotations

from .base import Result, fail

JACCARD_FLOOR = 0.99


def _read_lines(path):
    with open(path, "rb") as handle:
        text = handle.read().decode("utf-8", errors="replace")
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    return lines


def compare(path_a, path_b, cls, opts=None):
    with open(path_a, "rb") as handle:
        raw_a = handle.read()
    with open(path_b, "rb") as handle:
        raw_b = handle.read()
    lines_a = _read_lines(path_a)
    lines_b = _read_lines(path_b)

    if lines_a == lines_b:
        # E0 is reserved for a byte comparison; equal lines with a different
        # final newline are E2, the value exact class.
        return Result(True, "E0" if raw_a == raw_b else "E2",
                      {"lines": len(lines_a), "jaccard": 1.0})
    if cls != "E5":
        diffs = []
        for index, (left, right) in enumerate(zip(lines_a, lines_b), start=1):
            if left != right:
                diffs.append(f"line {index}: {left!r} vs {right!r}")
            if len(diffs) >= 20:
                break
        if len(lines_a) != len(lines_b):
            diffs.insert(0, f"line count differs: {len(lines_a)} vs {len(lines_b)}")
        return fail("files differ", *diffs, lines_py=len(lines_a),
                    lines_cpp=len(lines_b))

    set_a = set(lines_a)
    set_b = set(lines_b)
    union = len(set_a | set_b)
    intersection = len(set_a & set_b)
    jaccard = 1.0 if union == 0 else intersection / union
    metrics = {
        "lines_py": len(lines_a),
        "lines_cpp": len(lines_b),
        "only_py": len(set_a - set_b),
        "only_cpp": len(set_b - set_a),
        "jaccard": jaccard,
    }
    if jaccard < JACCARD_FLOOR:
        diffs = [f"only in python: {line!r}" for line in sorted(set_a - set_b)[:10]]
        diffs += [f"only in cpp: {line!r}" for line in sorted(set_b - set_a)[:10]]
        return Result(False, None, metrics, diffs)
    diffs = [f"only in python: {line!r}" for line in sorted(set_a - set_b)[:10]]
    diffs += [f"only in cpp: {line!r}" for line in sorted(set_b - set_a)[:10]]
    return Result(True, "E5", metrics, diffs)
