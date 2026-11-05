"""Set agreement comparator for tools whose output is a set of called regions.

Class E5 of cpp/PLAN.md 5.1: the Jaccard index of the called intervals must be
at least 0.99. The lines are compared as whole records, so a call that differs
in any column counts as a disagreement.

It was written for hicMergeLoops, whose reference breaks a tie between two
equally wide overlapping loop anchors by CPython's set iteration order. Since
2026-10-22 the port reimplements intervaltree's tree and CPython's set layout
and is byte identical, so its cases are E0; hicDifferentialTAD still declares
E5 outputs. A per-case tolerance is forbidden by PLAN.md 9.5, so a class is the
place for set agreement.

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
    columns = (opts or {}).get("columns")
    if columns is not None:
        # The records are the first `columns` tab separated fields of a line
        # (for hicDifferentialTAD the region and its name), so that the
        # last digits of the p-values and statistics that follow do not make
        # a called region a different record. Those columns are compared by
        # the text comparator of the same file. Lines of a comment keep their
        # full text.
        def key(line):
            return line if line.startswith("#") else "\t".join(line.split("\t")[:columns])
        lines_a = [key(line) for line in lines_a]
        lines_b = [key(line) for line in lines_b]

    if lines_a == lines_b:
        # E0 is reserved for a byte comparison; equal lines with a different
        # final newline are E2, the value exact class.
        met = "E0" if raw_a == raw_b else ("E5" if columns is not None else "E2")
        return Result(True, met,
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
