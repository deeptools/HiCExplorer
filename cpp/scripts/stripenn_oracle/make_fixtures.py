#!/usr/bin/env python
"""Writes the small fixtures of cpp/tests/test_detect_stripes.cpp and
cpp/tests/cli/hicDetectStripes_oracle.sh: inputs and the outputs that
Stripenn 1.1.65.22, numpy, scipy, scikit-image, OpenCV and Python's random
module produce for them. Run with the Stripenn environment's interpreter:
  ~/miniconda3/envs/stripenn-bench/bin/python make_fixtures.py OUT_DIR [FRAMES_DIR]
FRAMES_DIR (oracle_dump.py frames output) supplies a real gray image for the
Canny fixture."""
import json
import os
import random
import sys
import tempfile
from types import SimpleNamespace

import cv2 as cv
import numpy as np
import pandas as pd
import scipy.ndimage as ndi
import skimage.feature as feat
from skimage.feature import _canny_cy
from skimage.feature._canny import _preprocess

out = sys.argv[1]
frames = sys.argv[2] if len(sys.argv) > 2 else None
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(20261029)


def save(name, arr):
    np.save(os.path.join(out, name), np.ascontiguousarray(arr))


# --- Python random: random() and choices()
r = random.Random(123456789)
save("random.npy", np.array([r.random() for _ in range(3000)]))
r = random.Random(123456789)
picks = []
for n in [1, 2, 7, 100, 1234]:
    picks += r.choices(list(range(n)), k=50)
save("choices.npy", np.array(picks, dtype=np.int64))

# --- numpy sums and means of views
big = rng.random((60, 60)) * 3
q = []
for _ in range(400):
    R = int(rng.integers(1, 30)); C = int(rng.integers(1, 30))
    q.append((int(rng.integers(0, 60 - R + 1)), int(rng.integers(0, 60 - C + 1)), R, C))
q += [(0, 0, 60, 60), (0, 0, 1, 60), (0, 0, 60, 1)]
q = np.array(q, dtype=np.int64)
save("sums_big.npy", big)
save("sums_queries.npy", q)
save("sums_expected.npy", np.array([[np.sum(big[a:a + R, b:b + C]), np.mean(big[a:a + R, b:b + C])] for a, b, R, C in q]))
save("rowmeans_expected.npy", np.concatenate([np.mean(big[a:a + R, b:b + C], axis=1) for a, b, R, C in q]))

# --- numpy median and quantile
vals = np.round(rng.random(5000) * 20, 1) + 0.1
upper = []
n = 90
mat = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        if rng.random() < 0.6:
            mat[i, j] = mat[j, i] = round(float(rng.random() * 30), 2) + 0.01
upper = [(i, j, mat[i, j]) for i in range(n) for j in range(i, n) if mat[i, j] > 0]
save("quantile_upper.npy", np.array(upper))
qs = [0.0, 0.25, 0.5, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0]
save("quantile_q.npy", np.array(qs))
save("quantile_expected.npy", np.array([np.quantile(mat[mat > 0], x) for x in qs]))
save("median_in.npy", vals[:1001]); save("median_expected.npy", np.array([np.median(vals[:1001])]))
save("median_in_even.npy", vals[:1000]); save("median_expected_even.npy", np.array([np.median(vals[:1000])]))

# --- scipy gaussian filter, sobel and skimage canny, stage by stage
def canny_stages(gray, sigma):
    smoothed, eroded = _preprocess(gray, None, sigma, "constant", 0.0)
    j = ndi.sobel(smoothed, axis=1)
    i = ndi.sobel(smoothed, axis=0)
    mag = i * i
    mag += j * j
    np.sqrt(mag, out=mag)
    low = _canny_cy._nonmaximum_suppression_bilinear(i, j, mag, eroded, 0.1)
    return smoothed, i, j, mag, low

grays = []
if frames:
    import glob
    f = sorted(glob.glob(frames + "/frames/*/*_gray.npy"))[3]
    g = np.load(f)[2]
    grays.append(("real", np.ascontiguousarray(g[100:190, 60:150]), 2.0))
g = ndi.uniform_filter(rng.random((70, 85)).astype(np.float32), 3).astype(np.float32)
grays.append(("random", g, 1.0))
for name, g, sigma in grays:
    sm, i, j, mag, low = canny_stages(g, sigma)
    save("canny_%s_in.npy" % name, g)
    save("canny_%s_sigma.npy" % name, np.array([sigma]))
    save("canny_%s_smoothed.npy" % name, sm)
    save("canny_%s_isobel.npy" % name, i)
    save("canny_%s_magnitude.npy" % name, mag)
    save("canny_%s_suppressed.npy" % name, low)
    save("canny_%s_edges.npy" % name, feat.canny(g, sigma=sigma).astype(np.uint8))

from scipy.ndimage._filters import _gaussian_kernel1d
save("gauss_kernel_sigma2.npy", _gaussian_kernel1d(2.0, 0, 8))

# --- OpenCV box filter and RGB2GRAY
plane = rng.random((23, 31))
save("box_plane.npy", plane)
save("box_out3.npy", cv.filter2D(plane, -1, np.ones((3, 3)) / 9))
save("box_out5.npy", cv.filter2D(plane, -1, np.ones((5, 5)) / 25))
px = rng.random((101, 3)).astype(np.float32)
px[:, 1] = px[:, 2]
save("gray_px.npy", px)
save("gray_out.npy", cv.cvtColor(px.reshape(1, -1, 3), cv.COLOR_RGB2GRAY).ravel())

# --- Python repr of floats
reprs = [0.1 + 0.2, 1e-05, 1e16, 123456789012345680.0, 1.5e-7, -0.0, 0.0, 5e-324, 1.7976931348623157e308, 0.001, 0.0001, 12345.678,
         175.22304601696706, 0.734, -4.56370604659623, 1e22, 100.0, 99.99999999999999, 2.5e-05, 95.0]
json.dump({"values": [float.hex(v) for v in reprs], "reprs": [repr(v) for v in reprs]}, open(os.path.join(out, "repr.json"), "w"))

# --- nulldist's allocation of the samples over the chromosomes
def allocate(sizes, available):
    sizes = np.array(sizes)
    samplesize = (sizes / np.sum(sizes)) * 1000
    samplesize = np.uint64(samplesize)
    dif = 1000 - np.sum(samplesize)
    notzero = np.where(samplesize != 0)
    names = list(notzero[0])
    avail = [available[i] for i in names]
    ss = (avail / np.sum(avail)) * 1000
    ss = np.uint64(ss)
    dif = 1000 - np.sum(ss)
    notzero = np.where(ss != 0)
    names = [names[i] for i in notzero[0]]
    ss[0] += dif
    res = [-1] * len(sizes)
    for c, name in enumerate(names):
        res[name] = int(ss[c])
    return res

cases = []
for sizes, available in [
    ([248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309],
     [20000, 19000, 16000, 15000, 14500, 13900, 13000, 11800, 11000, 10900, 11100, 10800]),
    ([100000000, 90000000, 80000000, 70000000, 3000000, 2000000],
     [9000, 8000, 7000, 6000, 3, 1]),
    ([100000000, 90000000, 80000000, 70000000, 50000000],
     [9000, 0, 7000, 6000, 5000]),
]:
    cases.append({"sizes": sizes, "available": available, "expected": allocate(sizes, available)})
json.dump(cases, open(os.path.join(out, "allocate.json"), "w"))

# --- RemoveRedundant by size and by p-value
from stripenn import getStripe
obj = getStripe.getStripe.__new__(getStripe.getStripe)
obj.core = 1
rows = []
for k in range(70):
    num = int(rng.integers(0, 5))
    x0 = int(rng.integers(0, 40)) * 10000 + 1
    w = int(rng.integers(2, 9)) * 10000
    y0 = x0 - int(rng.integers(5, 30)) * 10000
    ln = int(rng.integers(10, 40)) * 10000
    rows.append(dict(chr="1", pos1=x0, pos2=x0 + w - 1, pos3=y0, pos4=y0 + ln - 1, h=int(ln // 10000), w=int(w // 10000),
                     num=num, pvalue=float(np.round(rng.random(), 3)), row=k))
df = pd.DataFrame(rows).sort_values("num", kind="stable").reset_index(drop=True)
df["row"] = range(len(df))
df.to_csv(os.path.join(out, "redundant_in.tsv"), sep="\t", index=False)
for by in ["size", "pvalue"]:
    kept = getStripe.getStripe.RemoveRedundant(obj, df.copy(), by)
    json.dump([int(v) for v in kept["row"]], open(os.path.join(out, "redundant_%s.json" % by), "w"))

# --- StripeSearch on a synthetic frame
S = 100
res = 10000
sub = rng.gamma(2.0, 1.0, size=(S, S))
sub = np.triu(sub) + np.triu(sub, 1).T
d = np.abs(np.arange(S)[:, None] - np.arange(S)[None, :])
sub = sub * np.exp(-d / 40.0) * 6
for col, top, wid in [(30, 5, 3), (62, 20, 4), (80, 40, 2)]:
    for c in range(col, col + wid):
        for r in range(top, c):
            sub[r, c] *= 5
            sub[c, r] = sub[r, c]
sub = np.round(sub, 3)
np.fill_diagonal(sub, np.round(sub.diagonal(), 3))
M = float(np.quantile(sub[sub > 0], 0.95))
save("search_submat.npy", sub)
save("search_M.npy", np.array([M]))
obj2 = getStripe.getStripe.__new__(getStripe.getStripe)
obj2.minH = 10; obj2.maxW = 8; obj2.canny = 2.0; obj2.bfilter = 3; obj2.core = 1; obj2.resol = res
sa = [j * res + 1 for j in range(S)]
ea = [s + res - 1 for s in sa]
ea[-1] = S * res
result = getStripe.getStripe.StripeSearch(obj2, sub.copy(), 0, 0, S - 1, M, 0.95, "c1", S, sa, ea)
result.to_csv(os.path.join(out, "search_expected.tsv"), sep="\t", index=False)
print("search boxes:", len(result))

# --- a small cooler and Stripenn's own result on it
import cooler
tmp = tempfile.mkdtemp(dir=os.environ.get("TMPDIR", "/tmp"))
chrom_bins = {"c1": 520, "c2": 530, "c3": 560}
bins = []
for name, nb in chrom_bins.items():
    for k in range(nb):
        bins.append((name, k * res, min((k + 1) * res, nb * res - 4000)))
bins = pd.DataFrame(bins, columns=["chrom", "start", "end"])
band = 130
pix = []
offset = 0
stripes = {"c1": [(100, 60, 3), (300, 90, 2)], "c2": [(150, 70, 3)], "c3": [(120, 80, 4), (400, 60, 3)]}
for name, nb in chrom_bins.items():
    mean = np.zeros((nb, nb))
    for i in range(nb):
        for j in range(i, min(nb, i + band)):
            mean[i, j] = 3.0 * np.exp(-(j - i) / 25.0) + 0.05
    for (col, length, wid) in stripes[name]:
        for c in range(col, col + wid):
            for r in range(max(0, c - length), c):
                mean[r, c] *= 6
    counts = rng.poisson(mean)
    for i in range(nb):
        for j in range(i, min(nb, i + band)):
            if counts[i, j] > 0 and rng.random() < 0.8:
                pix.append((offset + i, offset + j, int(counts[i, j])))
    offset += nb
pix = pd.DataFrame(pix, columns=["bin1_id", "bin2_id", "count"])
weight = rng.normal(1.0, 0.12, size=len(bins))
weight[[7, 8, 560, 900]] = np.nan
bins["weight"] = weight
path = os.path.join(out, "tiny_stripes.cool")
if os.path.exists(path):
    os.remove(path)
cooler.create_cooler(path, bins[["chrom", "start", "end"]], pix, columns=["count"], dtypes={"count": np.int32}, symmetric_upper=True,
                     ordered=True, assembly="synthetic")
import h5py
with h5py.File(path, "a") as h:
    h["bins"].create_dataset("weight", data=weight, compression="gzip")
print("cool bytes", os.path.getsize(path))

# --- Stripenn's own tables on that cooler, unmasked and masked (two cores)
import shutil
import subprocess
stripenn_bin = os.path.join(os.path.dirname(sys.executable), "stripenn")
for tag, extra in [("tiny", []), ("tiny_mask", ["--mask", "c1:600000-630000"])]:
    run_dir = os.path.join(tmp, tag + "_out")
    subprocess.run([stripenn_bin, "compute", "--cool", path, "--out", run_dir, "--norm", "weight", "--numcores", "2"] + extra,
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    shutil.copy(os.path.join(run_dir, "result_unfiltered.tsv"), os.path.join(out, tag + "_unfiltered.tsv"))
    shutil.copy(os.path.join(run_dir, "result_filtered.tsv"), os.path.join(out, tag + "_filtered.tsv"))
