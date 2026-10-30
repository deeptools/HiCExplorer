#!/usr/bin/env python
"""Randomised fixtures for tests/test_stripes_primitives.cpp: the expected
results of numpy, Python's random module, scikit-image, OpenCV and SciPy on
random inputs (ties, extremes, degenerate shapes). Run with the Stripenn
environment's interpreter (numpy 1.26.4, scikit-image 0.24.0, OpenCV 4.9.0):
  ~/miniconda3/envs/stripenn-bench/bin/python make_primitive_fixtures.py OUTDIR
The library is used as a black box; nothing of it is copied. The large
randomised comparison against a freshly built probe is check_primitives.py.
"""
import os
import random
import sys
import warnings

import cv2 as cv
import numpy as np
import scipy.signal as sig
import skimage.feature as feat
from skimage.feature import _canny_cy
from skimage.feature._canny import _preprocess
import scipy.ndimage as ndi

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(20261030)


def save(name, array):
    np.save(os.path.join(out, "prim_" + name + ".npy"), array)


# random.Random: 12 seeds (below 2^63), 1300 draws each, so the state is regenerated twice
seeds = np.array([0, 1, 2, 2**31, 2**32 - 1, 2**32, 2**32 + 1, 2**62 + 12345, 123456789, 987654321987, 7, 42], dtype=np.int64)
save("seeds", seeds)
draws = np.zeros((len(seeds), 1300))
for i, s in enumerate(seeds):
    r = random.Random(int(s))
    draws[i] = [r.random() for _ in range(1300)]
save("draws", draws)
picks = np.zeros((len(seeds), 3, 200), dtype=np.int64)
sizes = [1, 17, 99991]
for i, s in enumerate(seeds):
    for j, size in enumerate(sizes):
        r = random.Random(int(s))
        picks[i, j] = [r.choices(range(size), k=1)[0] for _ in range(200)]
save("picks", picks)
cum = np.cumsum(np.round(rng.random(23) * 4) + 0.5)
save("cumweights", cum)
wpicks = np.zeros((len(seeds), 200), dtype=np.int64)
for i, s in enumerate(seeds):
    r = random.Random(int(s))
    wpicks[i] = [r.choices(range(len(cum)), cum_weights=list(cum), k=1)[0] for _ in range(200)]
save("wpicks", wpicks)

# median / nanmedian / quantile / round over ragged arrays: values concatenated, lengths apart
vals, lens, med, nanmed = [], [], [], []
for t in range(300):
    n = int(rng.integers(1, 40))
    x = rng.integers(-4, 5, n).astype(float) / 2 if t % 3 == 0 else rng.normal(size=n) * 100
    if t % 7 == 0:
        x[0] = -0.0
    if t % 11 == 0:
        x[-1] = np.inf
    if t % 5 == 0:
        x[rng.random(n) < 0.3] = np.nan
    vals.append(x)
    lens.append(n)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        med.append(np.median(x))
        nanmed.append(np.nanmedian(x))
save("ragged_values", np.concatenate(vals))
save("ragged_lengths", np.array(lens, dtype=np.int64))
save("ragged_median", np.array(med))
save("ragged_nanmedian", np.array(nanmed))
q_in = np.round(rng.normal(size=777) * 4, 1)
qs = np.concatenate([rng.random(60), [0.0, 1.0, 0.5, 0.25, 0.75, 0.999]])
save("quantile_in", q_in)
save("quantile_q", qs)
save("quantile_out", np.quantile(q_in, qs))
q_small = np.array([3.0])
save("quantile_one_in", q_small)
save("quantile_one_out", np.quantile(q_small, qs))
r_in = np.concatenate([rng.normal(size=500) * 50, [0.5, 1.5, 2.5, -0.5, -1.5, -2.5, 0.49999999999999994, 2.5000000000000004, 1e16, -0.0]])
save("round_in", r_in)
save("round_out", np.round(r_in))

# canny: 36 images 37 x 41 of five kinds, and 12 tiny ones (4 x 5), with sigmas and thresholds
def make(kind, h, w):
    if kind == 0:
        g = rng.random((h, w))
    elif kind == 1:
        g = rng.integers(0, 5, (h, w)) / 4.0
    elif kind == 2:
        g = np.zeros((h, w))
        for _ in range(4):
            a, b = rng.integers(0, h), rng.integers(0, w)
            g[a:a + int(rng.integers(2, 20)), b:b + int(rng.integers(2, 20))] = rng.random()
    elif kind == 3:
        g = np.full((h, w), float(rng.random()))
    else:
        yy, xx = np.mgrid[0:h, 0:w]
        g = ((yy + xx) % 7) / 12.0 + (yy > xx) * 0.3
    return g.astype(np.float32)


for tag, (h, w, count) in {"big": (37, 41, 36), "tiny": (4, 5, 12)}.items():
    imgs = np.stack([make(i % 5, h, w) for i in range(count)])
    sig_ = np.array([[0.5, 1.0, 2.0, 2.0, 2.5, 3.0][i % 6] for i in range(count)])
    low = np.array([0.1 if i % 3 == 0 else float(np.round(rng.random() * 0.3 + 0.001, 4)) for i in range(count)])
    high = np.array([0.2 if i % 3 == 0 else float(low[i] + rng.random() * 0.4) for i in range(count)])
    edges = np.zeros((count, h, w), dtype=np.uint8)
    supp = np.zeros((count, h, w), dtype=np.float32)
    smooth = np.zeros((count, h, w), dtype=np.float32)
    for i in range(count):
        s, e = _preprocess(imgs[i], None, float(sig_[i]), "constant", 0.0)
        isb, jsb = ndi.sobel(s, axis=0), ndi.sobel(s, axis=1)
        mag = isb * isb
        mag += jsb * jsb
        np.sqrt(mag, out=mag)
        supp[i] = _canny_cy._nonmaximum_suppression_bilinear(isb, jsb, mag, e, float(low[i]))
        smooth[i] = s
        edges[i] = feat.canny(imgs[i], sigma=float(sig_[i]), low_threshold=float(low[i]), high_threshold=float(high[i]))
    save("canny_%s_in" % tag, imgs)
    save("canny_%s_params" % tag, np.stack([sig_, low, high], axis=1))
    save("canny_%s_edges" % tag, edges)
    save("canny_%s_suppressed" % tag, supp)
    save("canny_%s_smoothed" % tag, smooth)

# box filter: 24 planes 17 x 23 (and 3 x 2, 1 x 9) with k = 1, 3, 5, 7
planes, ks, outs = [], [], []
for i in range(24):
    k = [1, 3, 5, 7][i % 4]
    p = rng.integers(0, 6, (17, 23)).astype(float) / 3 if i % 2 else rng.random((17, 23)) * 100
    planes.append(p)
    ks.append(k)
    outs.append(cv.filter2D(p, -1, np.ones((k, k)) / (k * k)))
save("box_in", np.stack(planes))
save("box_k", np.array(ks, dtype=np.int64))
save("box_out", np.stack(outs))
for name, shape in {"small": (3, 2), "row": (1, 9)}.items():
    p = rng.random(shape)
    save("box_%s_in" % name, p)
    save("box_%s_out" % name, np.stack([cv.filter2D(p, -1, np.ones((k, k)) / (k * k)) for k in (1, 3, 5)]))

# RGB to grey: float32 and uint8
px = rng.random((20000, 3)).astype(np.float32)
px[:5000] = np.round(px[:5000] * 8) / 8
save("gray_px", px)
save("gray_out", cv.cvtColor(px.reshape(-1, 1, 3), cv.COLOR_RGB2GRAY).ravel())
px8 = rng.integers(0, 256, (60000, 3)).astype(np.uint8)
corner = np.array([[a, b, c] for a in (0, 1, 254, 255) for b in (0, 1, 254, 255) for c in (0, 1, 254, 255)], dtype=np.uint8)
px8 = np.concatenate([px8, corner])
save("gray8_px", px8)
save("gray8_out", cv.cvtColor(px8.reshape(-1, 1, 3), cv.COLOR_RGB2GRAY).ravel())

# convolve2d same: 30 binary matrices 12 x 15 with random 3 x 3 integer kernels
ms = (rng.random((30, 12, 15)) < 0.4).astype(np.uint8)
kern = rng.integers(-3, 4, (30, 3, 3)).astype(np.int32)
save("conv_in", ms)
save("conv_kernel", kern)
save("conv_out", np.stack([sig.convolve2d(ms[i], kern[i], mode="same") for i in range(30)]).astype(np.int32))
