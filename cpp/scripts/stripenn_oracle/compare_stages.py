#!/usr/bin/env python
"""Compares hicDetectStripes --dumpStages output with the Stripenn oracle's
dump (oracle_dump.py), stage by stage, in the order of Stripenn's compute():
  a maxpixel quantiles, b expected values, c null distribution,
  d per frame: image, edges, vertical lines, StripeSearch boxes,
  e extract() tables (RemoveRedundant by size, p-values),
  f the final table (RemoveRedundant by p-value, stripiness).
Usage: compare_stages.py ORACLE_DIR CPP_DIR [--chroms 21,22] [--stage a,b,c,d,e,f]
Run with any Python that has numpy and pandas."""
import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd


def same(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        return False, "shape %s vs %s" % (a.shape, b.shape)
    if a.dtype.kind == "f" or b.dtype.kind == "f":
        eq = (a == b) | (np.isnan(a.astype(float)) & np.isnan(b.astype(float)))
    else:
        eq = a == b
    n = int((~eq).sum())
    return n == 0, "%d of %d differ" % (n, a.size)


def tables_equal(o, c, cols):
    """Row by row equality of two result tables on `cols`."""
    msg = []
    if len(o) != len(c):
        msg.append("rows %d vs %d" % (len(o), len(c)))
    n = min(len(o), len(c))
    bad = {}
    for col in cols:
        if col not in o.columns or col not in c.columns:
            continue
        a = o[col].iloc[:n].to_numpy()
        b = c[col].iloc[:n].to_numpy()
        if a.dtype.kind in "fiu" and b.dtype.kind in "fiu":
            eq = (a == b) | (pd.isna(a) & pd.isna(b))
        else:
            eq = a.astype(str) == b.astype(str)
        if (~eq).sum():
            bad[col] = int((~eq).sum())
    if bad:
        msg.append("differing rows by column: %s" % bad)
    return not msg, "; ".join(msg) if msg else "%d rows identical" % n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("oracle")
    ap.add_argument("cpp")
    ap.add_argument("--chroms", default="")
    ap.add_argument("--stage", default="a,b,c,d,e,f")
    ap.add_argument("--max-frames", type=int, default=0)
    a = ap.parse_args()
    stages = a.stage.split(",")
    rc = 0

    def say(name, ok, detail=""):
        nonlocal rc
        print(("ok     " if ok else "DIFFER ") + name + ("  " + detail if detail else ""))
        if not ok:
            rc = 1

    O, C = a.oracle, a.cpp
    if "a" in stages:
        jo = json.load(open(O + "/maxpixel.json"))
        jc = json.load(open(C + "/maxpixel.json"))
        for k in jo:
            if a.chroms and k not in a.chroms.split(","):
                continue
            say("a maxpixel chr%s" % k, jo[k] == jc.get(k), "oracle %s cpp %s" % (jo[k], jc.get(k)))
    if "b" in stages:
        for f in sorted(glob.glob(O + "/expected_*.npy")):
            k = os.path.basename(f)[9:-4]
            if a.chroms and k not in a.chroms.split(","):
                continue
            ok, d = same(np.load(f), np.load(C + "/expected_%s.npy" % k))
            say("b expected chr%s" % k, ok, d)
    if "c" in stages:
        for n in ["leftup", "rightup", "leftdown", "rightdown"]:
            o = np.load(O + "/null_%s.npy" % n)
            c = np.load(C + "/null_%s.npy" % n)
            ok, d = same(o, c)
            if not ok and o.shape == c.shape:
                cols = np.where(~(((o == c) | (np.isnan(o) & np.isnan(c))).all(axis=0)))[0]
                d += "; first differing sample columns %s" % cols[:8]
            say("c null_%s" % n, ok, d)
    if "d" in stages:
        dirs = sorted(os.listdir(O + "/frames"))
        if a.max_frames:
            dirs = dirs[:a.max_frames]
        tot = {"frames": 0, "bad": 0}
        for d in dirs:
            for f in sorted(glob.glob(O + "/frames/%s/*_result.tsv" % d)):
                chrom = os.path.basename(f).split("_")[0]
                if a.chroms and chrom not in a.chroms.split(","):
                    continue
                tot["frames"] += 1
                pre = "/frames/%s/%s_" % (d, chrom)
                problems = []
                for k, src in [("submat", "submat"), ("start_array", "start_array"), ("end_array", "end_array"),
                               ("gray", "gray"), ("edges", "edges"), ("vert", "vert")]:
                    fo = O + pre + src + ".npy"
                    fc = C + pre + src + ".npy"
                    if not os.path.exists(fc):
                        problems.append("%s missing in cpp" % k)
                        continue
                    ok, det = same(np.load(fo), np.load(fc))
                    if not ok:
                        problems.append("%s: %s" % (k, det))
                to = pd.read_csv(f, sep="\t")
                fc = C + pre + "result.tsv"
                if os.path.exists(fc):
                    tc = pd.read_csv(fc, sep="\t")
                    ok, det = tables_equal(to, tc, ["pos1", "pos2", "pos3", "pos4", "length", "width", "total", "Mean", "num", "start", "end", "x", "y", "h", "w"])
                    if not ok:
                        problems.append("boxes: " + det)
                else:
                    problems.append("result.tsv missing in cpp")
                if problems:
                    tot["bad"] += 1
                    print("DIFFER d frame %s chr%s: %s" % (d, chrom, " | ".join(problems)))
        say("d frames compared", tot["bad"] == 0, "%d frames, %d with differences" % (tot["frames"], tot["bad"]))
    if "e" in stages:
        for f in sorted(glob.glob(O + "/extract_[0-9].tsv")):
            i = os.path.basename(f)[8:-4]
            to = pd.read_csv(f, sep="\t")
            tc = pd.read_csv(C + "/extract_%s.tsv" % i, sep="\t")
            if a.chroms:
                to = to[to["chr"].astype(str).isin(a.chroms.split(","))]
            ok, det = tables_equal(to.reset_index(drop=True), tc.reset_index(drop=True),
                                   ["chr", "pos1", "pos2", "pos3", "pos4", "length", "width", "total", "Mean", "num", "start", "end", "x", "y", "h", "w", "pvalue"])
            say("e extract %s" % i, ok, det)
    if "f" in stages and os.path.exists(O + "/result_unfiltered.tsv"):
        to = pd.read_csv(O + "/result_unfiltered.tsv", sep="\t")
        tc = pd.read_csv(C + "/result_unfiltered.tsv", sep="\t")
        ok, det = tables_equal(to, tc, ["chr", "pos1", "pos2", "pos3", "pos4", "length", "width", "Mean", "maxpixel", "pvalue", "Stripiness"])
        say("f final table", ok, det)
    sys.exit(rc)


if __name__ == "__main__":
    main()
