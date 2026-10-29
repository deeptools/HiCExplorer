#!/usr/bin/env python
"""Stage-by-stage oracle for hicDetectStripes: runs Stripenn's own Python
(1.1.65.22, an external oracle, never a dependency of HiCExplorer) and dumps
the intermediates of its compute() as .npy files, so that the C++ port can be
compared stage by stage (see compare_stages.py).

Run it with the Stripenn environment's interpreter, e.g.
  ~/miniconda3/envs/stripenn-bench/bin/python oracle_dump.py \
      --cool X.cool --out DIR --chrom 21,22 --numcores 8 stage1
Subcommands:
  stage1  maxpixel quantiles, expected values, null distribution
  frames  per frame / percentile: StripeSearch inputs and outputs, canny
          edges, vertical-line image (runs the search in-process)
  final   RemoveRedundant, p-values, stripiness, the final tables
Stripenn is never modified; the instrumentation wraps its functions.
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd


def build(a):
    """The set-up part of stripenn.stripenn.compute, line for line."""
    import cooler
    from stripenn import getStripe

    np.seterr(divide="ignore", invalid="ignore")
    maxpixel = list(map(float, a.maxpixel.split(",")))
    chroms = a.chrom.split(",")
    Lib = cooler.Cooler(a.cool)
    norm = a.norm
    PossibleNorm = Lib.bins().columns
    if norm == "None":
        norm = False
    elif norm == "weight":
        norm = True
    elif norm not in PossibleNorm:
        norm = False
    all_chromnames = Lib.chromnames
    all_chromsizes = Lib.chromsizes
    keep = [i for i in range(len(all_chromnames))
            if "JH5" not in all_chromnames[i] and "GL4" not in all_chromnames[i]
            and "RANDOM" not in all_chromnames[i] and all_chromnames[i] != "M"
            and all_chromnames[i] != "chrM" and all_chromnames[i] != "Y"
            and all_chromnames[i] != "chrY"]
    all_chromnames = [all_chromnames[i] for i in keep]
    all_chromsizes = all_chromsizes[keep]
    chromnames = all_chromnames
    chromsizes = all_chromsizes
    if chroms[0] != "all":
        idx = [all_chromnames.index(c) for c in chroms if c in all_chromnames]
        chromnames = chroms
        chromsizes = all_chromsizes[idx]
    unbalLib = Lib.matrix(balance=norm)
    obj = getStripe.getStripe(unbalLib, Lib.binsize, a.minL, a.maxW, a.canny,
                              all_chromnames, chromnames, all_chromsizes,
                              chromsizes, a.numcores, a.bfilter, a.seed)
    return Lib, obj, maxpixel, chromnames


def stage1(a):
    from stripenn import getStripe
    Lib, obj, maxpixel, chromnames = build(a)
    os.makedirs(a.out, exist_ok=True)
    t = time.time()
    MP = getStripe.getStripe.getQuantile_original(obj, Lib, chromnames, maxpixel)
    json.dump({k: [float(x) for x in v] for k, v in MP.items()},
              open(os.path.join(a.out, "maxpixel.json"), "w"))
    print("maxpixel", time.time() - t, flush=True)
    t = time.time()
    EV = getStripe.getStripe.mpmean(obj)
    for k, v in EV.items():
        np.save(os.path.join(a.out, "expected_%s.npy" % k), np.array(v, dtype=float))
    print("expected", time.time() - t, flush=True)
    if a.skip_null:
        return
    t = time.time()
    bg = getStripe.getStripe.nulldist(obj)
    for name, arr in zip(["leftup", "rightup", "leftdown", "rightdown"], bg):
        np.save(os.path.join(a.out, "null_%s.npy" % name), arr)
    print("null", time.time() - t, flush=True)


def frames(a):
    """Run extract() in-process (numcores forced to 1 so the wrappers see
    every call) and record each StripeSearch call."""
    from stripenn import getStripe, ImageProcessing
    import skimage.feature
    Lib, obj, maxpixel, chromnames = build(a)
    obj.core = 1
    os.makedirs(a.out, exist_ok=True)
    MP = {k: v for k, v in json.load(open(os.path.join(a.out, "maxpixel.json"))).items()}
    nul = [np.load(os.path.join(a.out, "null_%s.npy" % n)) for n in
           ["leftup", "rightup", "leftdown", "rightdown"]]

    cur = {}
    orig_canny = skimage.feature.canny
    orig_vert = ImageProcessing.verticalLine
    orig_ss = getStripe.getStripe.StripeSearch

    def canny(gray, *args, **kw):
        e = orig_canny(gray, *args, **kw)
        cur.setdefault("gray", []).append(np.array(gray))
        cur.setdefault("edges", []).append(np.array(e))
        return e

    def vert(M, *args, **kw):
        v = orig_vert(M, *args, **kw)
        cur.setdefault("vert", []).append(np.array(v))
        return v

    records = []

    def ss(self, submat, num, start, end, M, perc, chr, framesize, start_array, end_array):
        cur.clear()
        res = orig_ss(self, submat, num, start, end, M, perc, chr, framesize,
                      start_array, end_array)
        records.append(dict(num=num, perc=perc, chr=chr, M=M, submat=np.array(submat),
                            start_array=np.array(start_array), end_array=np.array(end_array),
                            gray=np.array(cur.get("gray", [])), edges=np.array(cur.get("edges", [])),
                            vert=np.array(cur.get("vert", [])), result=res))
        return res

    skimage.feature.canny = canny
    getStripe.feature.canny = canny
    ImageProcessing.verticalLine = vert
    getStripe.getStripe.StripeSearch = ss
    try:
        out = []
        for i, perc in enumerate(maxpixel):
            res = obj.extract(MP, i, perc, *nul)
            res.to_csv(os.path.join(a.out, "extract_%d.tsv" % i), sep="\t", index=False)
    finally:
        skimage.feature.canny = orig_canny
    for r in records:
        d = os.path.join(a.out, "frames", "p%s_f%d" % (r["perc"], r["num"]))
        os.makedirs(d, exist_ok=True)
        for k in ["submat", "start_array", "end_array", "gray", "edges", "vert"]:
            arr = r[k]
            if k in ("edges", "vert"):
                arr = arr.astype(np.uint8)
            np.save(os.path.join(d, "%s_%s.npy" % (r["chr"], k)), arr)
        r["result"].to_csv(os.path.join(d, "%s_result.tsv" % r["chr"]), sep="\t", index=False)
    print("frames recorded", len(records))


def final(a):
    from stripenn import getStripe
    Lib, obj, maxpixel, chromnames = build(a)
    obj.core = a.numcores
    MP = json.load(open(os.path.join(a.out, "maxpixel.json")))
    nul = [np.load(os.path.join(a.out, "null_%s.npy" % n)) for n in
           ["leftup", "rightup", "leftdown", "rightdown"]]
    EV = {os.path.basename(f)[9:-4]: np.load(os.path.join(a.out, f))
          for f in os.listdir(a.out) if f.startswith("expected_")}
    tabs = []
    for i, perc in enumerate(maxpixel):
        r = obj.extract(MP, i, perc, *nul)
        tabs.append(r)
    t = pd.concat(tabs)
    t.to_csv(os.path.join(a.out, "extract_concat.tsv"), sep="\t", index=False)
    t = getStripe.getStripe.RemoveRedundant(obj, df=t, by="pvalue")
    t.to_csv(os.path.join(a.out, "after_pvalue_redundant.tsv"), sep="\t", index=False)
    s = obj.scoringstripes(t, EV, a.mask)[0]
    t = t.drop(columns=["total", "num", "start", "end", "x", "y", "h", "w", "medpixel"])
    t.insert(t.shape[1], "Stripiness", s, True)
    t.to_csv(os.path.join(a.out, "result_unfiltered.tsv"), sep="\t", index=False)
    f = t[t["pvalue"] < a.pvalue].sort_values(by=["Stripiness"], ascending=False)
    f.to_csv(os.path.join(a.out, "result_filtered.tsv"), sep="\t", index=False)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("cmd", choices=["stage1", "frames", "final"])
    p.add_argument("--cool", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--norm", default="weight")
    p.add_argument("--chrom", default="all")
    p.add_argument("--canny", type=float, default=2.0)
    p.add_argument("--minL", type=int, default=10)
    p.add_argument("--maxW", type=int, default=8)
    p.add_argument("--maxpixel", default="0.95,0.96,0.97,0.98,0.99")
    p.add_argument("--numcores", type=int, default=8)
    p.add_argument("--pvalue", type=float, default=0.1)
    p.add_argument("--mask", default="0")
    p.add_argument("--bfilter", type=int, default=3)
    p.add_argument("--seed", type=int, default=123456789)
    p.add_argument("--skip-null", action="store_true")
    a = p.parse_args()
    {"stage1": stage1, "frames": frames, "final": final}[a.cmd](a)


if __name__ == "__main__":
    main()
