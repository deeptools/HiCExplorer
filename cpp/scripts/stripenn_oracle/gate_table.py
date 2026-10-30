#!/usr/bin/env python
"""Per-chromosome gate table: hicDetectStripes (one run per chromosome, with
/usr/bin/time -v logs) against Stripenn's result_unfiltered.tsv of the same
chromosome. Usage: gate_table.py STRIPENN_RUNS_DIR HICX_DIR   (stripenn_<chr>/ and hx_<chr>_unfiltered.tsv, hx_<chr>.log)"""
import re
import sys

import numpy as np
import pandas as pd

sd, hd = sys.argv[1], sys.argv[2]
print("chr\tstripenn\thicx\tidentical_set\tmax_dp\tstripiness_bit_equal\tcalled_stripenn\tcalled_hicx\tcalled_identical\twall_s\tpeak_MB\tload_before")
tot = [0, 0, 0]
for c in range(1, 23):
    o = pd.read_csv("%s/stripenn_%d/result_unfiltered.tsv" % (sd, c), sep="\t", dtype={"chr": str})
    h = pd.read_csv("%s/hx_%d_unfiltered.tsv" % (hd, c), sep="\t", dtype={"chr": str})
    key = ["chr", "pos1", "pos2", "pos3", "pos4", "maxpixel"]
    m = o.merge(h, on=key, how="outer", suffixes=("_s", "_h"), indicator=True)
    b = m[m["_merge"] == "both"]
    ident = (m["_merge"] == "both").all()
    dp = float(np.abs(b["pvalue_s"] - b["pvalue_h"]).max()) if len(b) else 0.0
    sb = int((b["Stripiness_s"].fillna(-9e9) == b["Stripiness_h"].fillna(-9e9)).sum())
    cs = set(map(tuple, o[o.pvalue < 0.1][key].astype(str).values))
    ch = set(map(tuple, h[h.pvalue < 0.1][key].astype(str).values))
    log = open("%s/hx_%d.log" % (hd, c)).read()
    wall = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", log).group(1)
    parts = [float(x) for x in wall.split(":")]
    secs = parts[-1] + 60 * parts[-2] + (3600 * parts[0] if len(parts) == 3 else 0)
    rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", log).group(1)) / 1024
    load = re.search(r"load before: (\S+)", log).group(1)
    print("%d\t%d\t%d\t%s\t%g\t%d/%d\t%d\t%d\t%s\t%.0f\t%.0f\t%s" % (c, len(o), len(h), ident, dp, sb, len(b), len(cs), len(ch), cs == ch, secs, rss, load))
    tot[0] += len(o)
    tot[1] += len(h)
    tot[2] += len(cs)
print("total\t%d\t%d" % (tot[0], tot[1]))
