"""Comparator for chicViewpointBackgroundModel's text output.

    Relative position  size nbinom  prob nbinom  max value  mean value

Compared here, at the declared class: the header line exactly, the line
count and order exactly, the position column exactly, and the max value and
mean value columns as floats at the class (E0 compares the printed text).

Not compared here: the size nbinom and prob nbinom columns. They are the
fitted parameters of a negative binomial whose reference is not reproducible
against itself (see validators/chic_background_model.py for the measurement),
and they are judged by the validators chic_background_likelihood (fit
likelihood) and chic_background_downstream (class E5 on the calls a
downstream tool makes). equiv.py refuses a case that uses this format without
declaring the likelihood validator on the same output, so the two columns can
never go unchecked. The only check applied to them here is that both parse as
floats in both files.
"""
from __future__ import annotations

from .base import Result, fail, values_agree

FIT_COLUMNS = (1, 2)
REQUIRED_VALIDATOR = "chic_background_likelihood"


def _lines(path):
    with open(path, "rb") as handle:
        text = handle.read().decode("utf-8", errors="replace")
    lines = text.split("\n")
    trailing_newline = bool(lines) and lines[-1] == ""
    if trailing_newline:
        lines.pop()
    return lines, trailing_newline


def compare(path_a, path_b, cls, opts=None):
    lines_a, newline_a = _lines(path_a)
    lines_b, newline_b = _lines(path_b)
    if newline_a != newline_b:
        return fail("the final newline differs")
    if len(lines_a) != len(lines_b):
        return fail(f"line count differs: {len(lines_a)} vs {len(lines_b)}")
    if not lines_a or lines_a[0] != lines_b[0]:
        return fail(f"header differs: {lines_a[:1]!r} vs {lines_b[:1]!r}")
    diffs = []
    exact_fields = 0
    for number, (line_a, line_b) in enumerate(zip(lines_a[1:], lines_b[1:]), start=2):
        fields_a = line_a.split("\t")
        fields_b = line_b.split("\t")
        if len(fields_a) != 5 or len(fields_b) != 5:
            diffs.append(f"line {number}: {len(fields_a)} vs {len(fields_b)} columns")
            continue
        if fields_a[0] != fields_b[0]:
            diffs.append(f"line {number} position: {fields_a[0]!r} vs {fields_b[0]!r}")
        for column in FIT_COLUMNS:
            try:
                float(fields_a[column])
                float(fields_b[column])
            except ValueError:
                diffs.append(f"line {number} column {column + 1}: not a float "
                             f"({fields_a[column]!r} vs {fields_b[column]!r})")
        for column in (3, 4):
            exact_fields += 1
            if fields_a[column] == fields_b[column]:
                continue
            if cls != "E0":
                try:
                    if values_agree(float(fields_b[column]), float(fields_a[column]), cls):
                        continue
                except ValueError:
                    pass
            diffs.append(f"line {number} column {column + 1}: {fields_a[column]!r} vs "
                         f"{fields_b[column]!r}")
        if len(diffs) >= 20:
            break
    metrics = {"lines": len(lines_a), "exact_fields_compared": exact_fields,
               "fit_columns": "judged by the chic_background_likelihood and "
                              "chic_background_downstream validators"}
    if diffs:
        return Result(False, None, metrics, diffs)
    return Result(True, cls, metrics)
