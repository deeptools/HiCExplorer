"""bigWig files (chicExportData), read through pyBigWig.

Compared, in this order:

  1. the chromosome list with its lengths, in file order;
  2. the header pyBigWig reports: version, number of zoom levels, bases
     covered, minimum, maximum, sum and sum of squares;
  3. every interval of every chromosome: start and end exactly, the value at
     the requested class (the file stores float32; E0 to E2 require the same
     bits, looser classes base.values_agree);
  4. the zoom levels, through pyBigWig's summary statistics (mean, min, max,
     coverage and std over 1 and 100 bins per chromosome, which pyBigWig
     answers from the zoom levels), at the same class as the values.

Run with the reference venv interpreter, which has pyBigWig.
"""
from __future__ import annotations

import math

try:
    import pyBigWig
except ImportError as error:  # pragma: no cover
    raise SystemExit("the bigwig comparator needs pyBigWig; run equiv.py with the reference venv "
                     "interpreter") from error

from .base import Result, fail, values_agree

_HEADER = ("version", "nLevels", "nBasesCovered", "minVal", "maxVal", "sumData", "sumSquared")


def _agree(candidate, reference, cls):
    if candidate is None or reference is None:
        return candidate is None and reference is None
    return values_agree(float(candidate), float(reference), cls)


def compare(path_a, path_b, cls, opts=None):
    try:
        reference = pyBigWig.open(str(path_a))
        candidate = pyBigWig.open(str(path_b))
    except RuntimeError as error:
        return fail(f"not a readable bigWig file: {error}")
    diffs = []
    metrics = {"chromosomes": 0, "intervals": 0, "values_differing": 0}
    try:
        if not reference.isBigWig() or not candidate.isBigWig():
            return fail("not a bigWig file")
        chroms_a = list(reference.chroms().items())
        chroms_b = list(candidate.chroms().items())
        if chroms_a != chroms_b:
            diffs.append(f"chromosomes {chroms_a[:5]} vs {chroms_b[:5]}")
        header_a = reference.header()
        header_b = candidate.header()
        for key in _HEADER:
            if header_a.get(key) != header_b.get(key):
                diffs.append(f"header {key}: {header_a.get(key)!r} vs {header_b.get(key)!r}")
        for chrom, _length in chroms_a:
            if chrom not in dict(chroms_b):
                continue
            metrics["chromosomes"] += 1
            intervals_a = reference.intervals(chrom) or ()
            intervals_b = candidate.intervals(chrom) or ()
            if len(intervals_a) != len(intervals_b):
                diffs.append(f"{chrom}: {len(intervals_a)} vs {len(intervals_b)} intervals")
                continue
            for index, ((start_a, end_a, value_a), (start_b, end_b, value_b)) in enumerate(
                    zip(intervals_a, intervals_b)):
                metrics["intervals"] += 1
                if (start_a, end_a) != (start_b, end_b):
                    if len(diffs) < 20:
                        diffs.append(f"{chrom}[{index}]: interval {start_a}-{end_a} vs {start_b}-{end_b}")
                elif not _agree(value_b, value_a, cls):
                    metrics["values_differing"] += 1
                    if len(diffs) < 20:
                        diffs.append(f"{chrom}[{index}] {start_a}-{end_a}: value {value_a!r} vs {value_b!r}")
            for statistic in ("mean", "min", "max", "coverage", "std"):
                for bins in (1, 100):
                    stats_a = reference.stats(chrom, type=statistic, nBins=bins)
                    stats_b = candidate.stats(chrom, type=statistic, nBins=bins)
                    for position, (value_a, value_b) in enumerate(zip(stats_a, stats_b)):
                        if not _agree(value_b, value_a, cls) and not (
                                value_a is not None and value_b is not None
                                and math.isnan(value_a) and math.isnan(value_b)):
                            if len(diffs) < 20:
                                diffs.append(f"{chrom} zoom {statistic} over {bins} bins [{position}]: "
                                             f"{value_a!r} vs {value_b!r}")
    finally:
        reference.close()
        candidate.close()
    if diffs:
        return Result(False, None, metrics, diffs[:20])
    return Result(True, cls, metrics)
