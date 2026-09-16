"""The result tables of hicDifferentialAnalysis (class EX, cpp/PLAN.md 9.7).

The tool has no Python counterpart, so a case cannot compare its output with
a reference. This validator checks what a correct table must satisfy on its
own:

- the provenance header, and the exploratory label exactly when the case
  expects one (option "exploratory");
- one field per column on every row;
- p-values and FDRs in [0, 1] or nan, and an FDR never below its p-value
  (Benjamini-Hochberg only raises a p-value);
- the differential flag set exactly where the FDR is at most the threshold
  the header records;
- option "min_rows": the least number of rows;
- option "called": "chrom:start-end" rows (a TAD, or a window for boundaries)
  that must be called differential, for a case with planted differences;
- option "max_called_fraction": the largest fraction of rows called.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

from comparators.base import Result, fail


def _float(text):
    try:
        return float(text)
    except ValueError:
        return None


def differential_analysis_table(context):
    path = Path(context["out_cpp"]) / context["output"]
    options = context.get("options") or {}
    if not path.is_file():
        return fail(f"{context['output']} was not written")
    lines = path.read_text().splitlines()
    comments = [line for line in lines if line.startswith("#")]
    if not comments or not comments[0].startswith("# hicDifferentialAnalysis "):
        return fail("the provenance header '# hicDifferentialAnalysis ...' is missing")
    headers = [line for line in comments if line.startswith("#chrom")]
    if not headers:
        return fail("the column header '#chrom ...' is missing")
    columns = headers[-1][1:].split("\t")
    diffs = []
    exploratory = any(line.startswith("# EXPLORATORY") for line in comments)
    if "exploratory" in options and exploratory != bool(options["exploratory"]):
        diffs.append(f"exploratory label is {exploratory}, the case expects {options['exploratory']}")
    threshold = None
    for line in comments:
        match = re.search(r"FDR ([0-9.eE+-]+)", line)
        if match and line.startswith("# minimum fold change"):
            threshold = float(match.group(1))
    if threshold is None:
        diffs.append("the FDR threshold is not recorded in the header")
        threshold = 0.05
    rows = [line.split("\t") for line in lines if line and not line.startswith("#")]
    for name in ("pvalue", "fdr", "differential"):
        if name not in columns:
            return fail(f"column {name} is missing", columns=columns)
    ip, iq, iflag = columns.index("pvalue"), columns.index("fdr"), columns.index("differential")
    tested = called = low = 0
    called_keys = set()
    for number, row in enumerate(rows, start=1):
        if len(row) != len(columns):
            diffs.append(f"row {number}: {len(row)} fields for {len(columns)} columns")
            continue
        p, q = _float(row[ip]), _float(row[iq])
        if p is None or q is None:
            diffs.append(f"row {number}: p-value {row[ip]!r} or FDR {row[iq]!r} is not a number")
            continue
        for label, value in (("p-value", p), ("FDR", q)):
            if not math.isnan(value) and not 0.0 <= value <= 1.0:
                diffs.append(f"row {number}: {label} {value} outside [0, 1]")
        if not math.isnan(p):
            tested += 1
            low += p <= 0.05
            if math.isnan(q) or q < p * (1 - 1e-9):
                diffs.append(f"row {number}: FDR {q} below its p-value {p}")
        expected_flag = "1" if (not math.isnan(q) and q <= threshold) else "0"
        if row[iflag] != expected_flag:
            diffs.append(f"row {number}: differential {row[iflag]} but FDR {q} and threshold {threshold}")
        if row[iflag] == "1":
            called += 1
            called_keys.add(f"{row[0]}:{row[1]}-{row[2]}")
            if "windowStart" in columns:
                called_keys.add(f"{row[0]}:{row[columns.index('windowStart')]}-{row[columns.index('windowEnd')]}")
    if len(rows) < options.get("min_rows", 1):
        diffs.append(f"{len(rows)} rows, the case expects at least {options.get('min_rows', 1)}")
    for key in options.get("called", []):
        if key not in called_keys:
            diffs.append(f"{key} is not called differential")
    if "max_called_fraction" in options and rows and called / len(rows) > options["max_called_fraction"]:
        diffs.append(f"{called} of {len(rows)} rows called, above {options['max_called_fraction']}")
    metrics = {"rows": len(rows), "tested": tested, "called": called,
               "fraction_p_at_most_0.05": (low / tested) if tested else None,
               "exploratory": exploratory, "fdr_threshold": threshold}
    return Result(passed=not diffs, class_met="EX" if not diffs else None, metrics=metrics,
                  diffs=diffs)
