#!/usr/bin/env python
"""Planted-stripe recovery of Stripenn and hicDetectStripes on one planted
matrix, with the rule of cpp/scripts/stripe_calibration.py: a plant is
recovered when a call of the same orientation has its anchor within 2 bins and
covers at least half of the planted length. The calls of both tools are the
rows of their unfiltered tables with pvalue < cutoff, read with one function,
so that the two sets are compared on the same footing.
Usage: planted_recovery.py PLANTS_JSON STRIPENN_UNFILTERED HICX_UNFILTERED [--bin 10000] [--cutoff 0.1]"""
import argparse
import collections
import json

import pandas as pd


def calls_from(path, cutoff):
    t = pd.read_csv(path, sep="\t", dtype={"chr": str})
    t = t[t["pvalue"] < cutoff]
    calls = []
    for r in t.itertuples():
        a0, a1, e0, e1 = r.pos1 - 1, r.pos2, r.pos3 - 1, r.pos4
        orientation = "vertical" if e1 == a1 else ("horizontal" if e0 == a0 else "unknown")
        calls.append(dict(chrom=r.chr, anchor_start=a0, orientation=orientation, extent_start=e0, extent_end=e1))
    return calls


def recovered(plant, calls, bin_size):
    for c in calls:
        if c["chrom"] != plant["chrom"] or c["orientation"] != plant["orientation"]:
            continue
        if abs(c["anchor_start"] - plant["anchor_bp"]) > 2 * bin_size:
            continue
        overlap = max(0, min(c["extent_end"], plant["extent_end"]) - max(c["extent_start"], plant["extent_start"]))
        if overlap >= 0.5 * plant["length_bp"]:
            return c
    return None


def candidates_for(plant, table, bin_size):
    """Plants with a candidate (any p-value) of the right orientation and anchor and length."""
    t = pd.read_csv(table, sep="\t", dtype={"chr": str})
    allc = []
    for r in t.itertuples():
        a0, a1, e0, e1 = r.pos1 - 1, r.pos2, r.pos3 - 1, r.pos4
        orientation = "vertical" if e1 == a1 else ("horizontal" if e0 == a0 else "unknown")
        allc.append(dict(chrom=r.chr, anchor_start=a0, orientation=orientation, extent_start=e0, extent_end=e1))
    return allc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("plants")
    ap.add_argument("stripenn")
    ap.add_argument("hicx")
    ap.add_argument("--bin", type=int, default=10000)
    ap.add_argument("--cutoff", type=float, default=0.1)
    a = ap.parse_args()
    plants = json.load(open(a.plants))
    out = {"plants": len(plants)}
    for name, path in [("Stripenn", a.stripenn), ("hicDetectStripes", a.hicx)]:
        calls = calls_from(path, a.cutoff)
        cands = candidates_for(None, path, a.bin)
        hit = [recovered(p, calls, a.bin) is not None for p in plants]
        has_candidate = [recovered(p, cands, a.bin) is not None for p in plants]
        strata = collections.OrderedDict()
        for p, h, hc in zip(plants, hit, has_candidate):
            key = "%dkb_x%s" % (p["length_bp"] // 1000, p["enrichment"])
            s = strata.setdefault(key, [0, 0, 0])
            s[0] += int(h); s[1] += int(hc); s[2] += 1
        out[name] = {"calls": len(calls), "candidates": len(cands), "recovered_as_call": int(sum(hit)),
                     "plants_with_candidate": int(sum(has_candidate)),
                     "by_length_enrichment [recovered, with candidate, plants]": strata}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
