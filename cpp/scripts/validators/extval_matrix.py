"""Validators of the external-data matrix cases (cpp/scripts/cases_external).

plotdata_ranking
    hicCorrelate --plotData JSON of a C++-only method (scc, genomedisco,
    hicspector, quasar). There is no 3.7.6 reference; the check is that the
    biological replicate pair scores higher than every pair across cell types.
    options: {"higher": [[i, j], ...], "lower": [[i, j], ...]} (indices into
    the --matrices list).

cool_pixels_match
    A C++ cool or mcool resolution against an independently written cool file
    (for example the hic2cool output of the 3.7.6 tool): same chromosomes,
    same bins and bit-identical pixels. options: {"expected": "{data}/...",
    "uri": "::/resolutions/1000000" (appended to the output path, optional)}.

qc_pretty_report
    The C++-only hicQC_report.html of hicPrepareQCreport/hicQC --prettyReport:
    the file is written, and every PASS/WARN/FAIL badge agrees with the status
    recomputed from QC_table.txt with the cutoffs and the WARN margin the tool
    documents. options: {"folder": "qc", "thresholds": {...}} (defaults are the
    tool's defaults).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from comparators.base import Result, fail


def _out(context, relative):
    return Path(context["out_cpp"]) / relative


def plotdata_ranking(context):
    path = _out(context, context["output"])
    if not path.is_file():
        return fail(f"{context['output']} was not written")
    payload = json.loads(path.read_text())
    results = np.array(payload["results"], dtype=float)
    results = np.maximum(results, results.T)  # the JSON stores the upper triangle
    options = context["options"]
    higher = [results[i, j] for i, j in options["higher"]]
    lower = [results[i, j] for i, j in options["lower"]]
    metrics = {"method": payload.get("method"), "labels": payload.get("labels"),
               "higher": higher, "lower": lower}
    if not all(np.isfinite(higher + lower)):
        return Result(False, None, metrics, ["a score is not finite"])
    if min(higher) <= max(lower):
        return Result(False, None, metrics,
                      [f"replicate pair score {min(higher):.4f} is not above the "
                       f"cross-cell-type score {max(lower):.4f}"])
    return Result(True, "EX", metrics, [])


def cool_pixels_match(context):
    import cooler  # pylint: disable=C0415

    options = context["options"]
    produced = str(_out(context, context["output"])) + options.get("uri", "")
    expected = options["expected"].replace("{data}", context["data"])
    try:
        a = cooler.Cooler(expected)
        b = cooler.Cooler(produced)
    except Exception as error:  # pylint: disable=W0718
        return fail(f"cannot open: {error}")
    diffs = []
    if list(a.chromnames) != list(b.chromnames):
        diffs.append(f"chromosomes differ: {a.chromnames[:5]} vs {b.chromnames[:5]}")
    elif not np.array_equal(a.chromsizes.values, b.chromsizes.values):
        diffs.append("chromosome lengths differ")
    bins_a, bins_b = a.bins()[["chrom", "start", "end"]][:], b.bins()[["chrom", "start", "end"]][:]
    if len(bins_a) != len(bins_b) or not (bins_a.start.values == bins_b.start.values).all() \
            or not (bins_a.end.values == bins_b.end.values).all():
        diffs.append(f"bins differ: {len(bins_a)} vs {len(bins_b)}")
    metrics = {"nnz_expected": int(a.info["nnz"]), "nnz_produced": int(b.info["nnz"]),
               "binsize": None if a.binsize is None else int(a.binsize)}
    if not diffs:
        pa, pb = a.pixels()[:], b.pixels()[:]
        if len(pa) != len(pb):
            diffs.append(f"pixel count {len(pa)} vs {len(pb)}")
        else:
            for column in ("bin1_id", "bin2_id"):
                if not np.array_equal(pa[column].values, pb[column].values):
                    diffs.append(f"{column} differs")
            va, vb = pa["count"].values.astype(float), pb["count"].values.astype(float)
            if not np.array_equal(va, vb):
                where = np.flatnonzero(va != vb)
                diffs.append(f"{len(where)} counts differ, first at pixel {where[0]}: "
                             f"{va[where[0]]} vs {vb[where[0]]}")
            metrics["sum"] = float(va.sum())
    return Result(not diffs, "EX" if not diffs else None, metrics, diffs)


_QC_LABELS = ("PASS", "WARN", "FAIL")


def _status_at_least(value, minimum, margin):
    if value >= minimum:
        return "PASS"
    return "WARN" if value >= minimum * (1.0 - margin) else "FAIL"


def _status_at_most(value, maximum, margin):
    if value <= maximum:
        return "PASS"
    return "WARN" if value <= maximum * (1.0 + margin) else "FAIL"


def _column(header, row, name):
    for title, value in zip(header, row):
        if title == name:
            return float(value)
    for title, value in zip(header, row):
        if title.startswith(name):
            return float(value)
    return None


def _labels(header, rows, html):
    """The name each sample goes by in the report: the --labels value when
    QC_table.txt is indexed by Labels, else the File entry."""
    return [row[0] for row in rows]


def qc_pretty_report(context):
    options = context["options"]
    folder = _out(context, options.get("folder", "qc"))
    report = folder / "hicQC_report.html"
    table = folder / "QC_table.txt"
    if not report.is_file() or report.stat().st_size == 0:
        return fail("hicQC_report.html was not written")
    if not table.is_file():
        return fail("QC_table.txt was not written")
    thresholds = {"valid": 25.0, "inter": 10.0, "duplicated": 20.0, "dangling": 10.0,
                  "self_ligation": 10.0, "margin": 0.2}
    thresholds.update(options.get("thresholds", {}))
    lines = [line.rstrip("\n").split("\t") for line in table.read_text().splitlines() if line]
    header, rows = lines[0], lines[1:]
    expected = []
    for row in rows:
        sequenced = _column(header, row, "Sequenced reads")
        valid = _column(header, row, "Hi-C contacts")
        mappable = _column(header, row, "Pairs mappable, unique and high quality")
        inter = _column(header, row, "inter chromosomal")
        statuses = {
            "valid": _status_at_least(100.0 * valid / sequenced, thresholds["valid"], thresholds["margin"]),
            "inter": _status_at_most(100.0 * inter / valid, thresholds["inter"], thresholds["margin"]),
            "duplicated": _status_at_most(100.0 * _column(header, row, "duplicated pairs") / mappable,
                                          thresholds["duplicated"], thresholds["margin"]),
            "dangling": _status_at_most(100.0 * _column(header, row, "dangling end") / mappable,
                                        thresholds["dangling"], thresholds["margin"]),
            "self_ligation": _status_at_most(
                100.0 * _column(header, row, "self ligation (removed)") / mappable,
                thresholds["self_ligation"], thresholds["margin"]),
        }
        expected.append(statuses)
    html = report.read_text(errors="replace")
    # Every badge is attributed to the metric and the sample named in the text
    # between it and the previous badge: the overview cards use the names
    # below, the section bars the QC log's own row names.
    names = {"valid": ("Valid pairs",), "inter": ("Inter-chromosomal contacts", "inter chromosomal"),
             "duplicated": ("Duplicated pairs", "duplicated pairs"),
             "dangling": ("Dangling ends", "dangling end"),
             "self_ligation": ("Self-ligation", "self ligation (removed)")}
    samples = [row[0] for row in rows]
    badge = re.compile(r'<span class="qc-badge[^"]*\bqc-(pass|warn|fail)\b[^"]*">\s*(PASS|WARN|FAIL)')
    previous = 0
    checked = 0
    current = 0 if len(samples) == 1 else None
    diffs = []
    labels = _labels(header, rows, html)
    for match in badge.finditer(html):
        text = re.sub(r"<[^>]+>", " ", html[previous:match.start()])
        previous = match.end()
        metric = max(((text.rfind(name), key) for key, alias in names.items() for name in alias),
                     default=(-1, None))
        sample = max(((text.rfind(label), index) for index, label in enumerate(labels)),
                     default=(-1, None))
        if metric[0] < 0:
            diffs.append(f"badge {match.group(2)} without a metric name before it")
            continue
        if sample[0] >= 0:
            current = sample[1]
        index = current
        if index is None:
            diffs.append(f"badge {match.group(2)} for {metric[1]} without a sample name")
            continue
        want = expected[index][metric[1]]
        checked += 1
        if match.group(2) != want or match.group(1) != want.lower():
            diffs.append(f"sample {labels[index]!r} {metric[1]}: badge {match.group(2)}, "
                         f"recomputed {want}")
    metrics = {"expected_status": expected, "badges_checked": checked}
    if checked < 5 * len(rows):
        diffs.append(f"only {checked} badges found for {len(rows)} samples x 5 metrics")
    return Result(not diffs, "EX" if not diffs else None, metrics, diffs)
