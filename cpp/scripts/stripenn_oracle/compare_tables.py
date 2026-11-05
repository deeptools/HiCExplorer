#!/usr/bin/env python
"""Gate comparison of two Stripenn result tables (result_unfiltered.tsv
layout): the candidate sets (chr, pos1..pos4), the p-values and stripiness of
the candidates present in both, and the called sets (pvalue < cutoff).
Usage: compare_tables.py STRIPENN_TSV HICX_TSV [--chrom 21] [--cutoff 0.1]
Prints one summary line; exit status 1 when anything differs."""
import argparse
import sys

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stripenn")
    ap.add_argument("hicx")
    ap.add_argument("--chrom", default=None)
    ap.add_argument("--cutoff", type=float, default=0.1)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    o = pd.read_csv(a.stripenn, sep="\t", dtype={"chr": str})
    c = pd.read_csv(a.hicx, sep="\t", dtype={"chr": str})
    if a.chrom:
        o = o[o["chr"] == a.chrom]
        c = c[c["chr"] == a.chrom]
    key = ["chr", "pos1", "pos2", "pos3", "pos4", "maxpixel"]
    mo = o.merge(c, on=key, how="outer", suffixes=("_s", "_h"), indicator=True)
    both = mo[mo["_merge"] == "both"]
    only_s = int((mo["_merge"] == "left_only").sum())
    only_h = int((mo["_merge"] == "right_only").sum())
    dp = np.abs(both["pvalue_s"] - both["pvalue_h"])
    ds = np.abs(both["Stripiness_s"] - both["Stripiness_h"])
    dm = np.abs(both["Mean_s"] - both["Mean_h"])
    same_str = int((both["Stripiness_s"].fillna(-999).astype(float) == both["Stripiness_h"].fillna(-999).astype(float)).sum())
    cs = set(map(tuple, o[o["pvalue"] < a.cutoff][key].astype(str).values))
    ch = set(map(tuple, c[c["pvalue"] < a.cutoff][key].astype(str).values))
    ident = only_s == 0 and only_h == 0 and (dp.max() if len(dp) else 0) == 0 and cs == ch
    print("%s stripenn=%d hicx=%d matched=%d only_stripenn=%d only_hicx=%d  max|dp|=%.3g max|dStripiness|=%.3g max|dMean|=%.3g stripiness_bit_equal=%d/%d  called stripenn=%d hicx=%d identical=%s  => %s" % (
        a.chrom or "all", len(o), len(c), len(both), only_s, only_h, dp.max() if len(dp) else 0,
        np.nanmax(ds) if len(ds) else 0, dm.max() if len(dm) else 0, same_str, len(both),
        len(cs), len(ch), cs == ch, "IDENTICAL" if ident else "DIFFERENT"))
    sys.exit(0 if ident else 1)


if __name__ == "__main__":
    main()
