#!/usr/bin/env python
"""Bit-for-bit comparison of the C++ primitives of hicDetectStripes
(tests/stripes_probe.cpp over tools/stripes_numpy.cpp) with numpy, scipy,
scikit-image, OpenCV and Python's random module. Run with the Stripenn
environment's interpreter:
  ~/miniconda3/envs/stripenn-bench/bin/python check_primitives.py PROBE_BINARY [FRAMES_DIR]
With FRAMES_DIR (oracle_dump.py frames output) the canny stages are also
checked on the gray images Stripenn itself produced."""
import os
import random
import subprocess
import sys
import tempfile

import cv2 as cv
import numpy as np
import scipy.ndimage as ndi
import skimage.feature as feat
from skimage.feature import _canny_cy
from skimage.filters import gaussian

probe = sys.argv[1]
frames = sys.argv[2] if len(sys.argv) > 2 else None
tmp = tempfile.mkdtemp(prefix="probe_", dir=os.environ.get("TMPDIR", "/tmp"))
fails = 0


def run(*args):
    subprocess.run([probe] + [str(a) for a in args], check=True)


def report(name, ok, detail=""):
    global fails
    print(("ok     " if ok else "FAILED ") + name + (" " + detail if detail else ""))
    if not ok:
        fails += 1


rng = np.random.default_rng(11)

# random.Random
for seed in [123456789, 0, 1, 2**32 + 5]:
    run("random", seed, 5000, tmp + "/r.npy")
    r = random.Random(seed)
    ref = np.array([r.random() for _ in range(5000)])
    report("random seed %d" % seed, np.array_equal(np.load(tmp + "/r.npy"), ref))

# gaussian kernels
bad = []
for sigma in [0.5, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.3, 4.0, 5.0, 6.0, 7.5]:
    run("kernel", sigma, tmp + "/k.npy")
    r = int(4.0 * sigma + 0.5)
    from scipy.ndimage._filters import _gaussian_kernel1d
    ref = _gaussian_kernel1d(sigma, 0, r)
    if not np.array_equal(np.load(tmp + "/k.npy"), ref):
        bad.append(sigma)
report("gaussian kernel1d", not bad, "mismatching sigmas %s" % bad)

# block sums and row means on views of a big array
big = rng.random((300, 300)) * 3
np.save(tmp + "/big.npy", big)
q = []
for _ in range(20000):
    R = rng.integers(1, 40); C = rng.integers(1, 40)
    r0 = rng.integers(0, 300 - R + 1); c0 = rng.integers(0, 300 - C + 1)
    q.append((r0, c0, R, C))
for (R, C) in [(300, 300), (1, 300), (300, 1), (200, 150), (60, 300), (9, 300)]:
    q.append((0, 0, R, C))
q = np.array(q, dtype=np.int64)
np.save(tmp + "/q.npy", q)
run("blocksums", tmp + "/big.npy", tmp + "/q.npy", tmp + "/o.npy")
o = np.load(tmp + "/o.npy")
ref = np.array([[np.sum(big[a:a + R, b:b + C]), np.mean(big[a:a + R, b:b + C])] for a, b, R, C in q])
mism = np.where((o != ref).any(axis=1))[0]
report("np.sum / np.mean of 2-D views", len(mism) == 0, "%d of %d differ" % (len(mism), len(q)))
run("rowmeans", tmp + "/big.npy", tmp + "/q.npy", tmp + "/o.npy")
o = np.load(tmp + "/o.npy")
ref = np.concatenate([np.mean(big[a:a + R, b:b + C], axis=1) for a, b, R, C in q])
report("np.mean(view, axis=1)", np.array_equal(o, ref), "%d of %d differ" % ((o != ref).sum(), len(ref)))

# gaussian_filter, float32, mode constant, as skimage.filters.gaussian
img = rng.random((137, 211)).astype(np.float32)
np.save(tmp + "/img.npy", img)
bad = []
for sigma in [0.7, 1.0, 2.0, 2.5, 3.0]:
    run("gauss", tmp + "/img.npy", sigma, tmp + "/g.npy")
    ref = ndi.gaussian_filter(img, sigma, mode="constant", cval=0.0)
    if not np.array_equal(np.load(tmp + "/g.npy"), ref):
        bad.append(sigma)
report("ndimage.gaussian_filter float32 constant", not bad, "sigmas %s" % bad)


def canny_stages(gray, sigma):
    # skimage.feature._canny.canny, stage by stage
    from skimage.feature._canny import _preprocess
    smoothed, eroded = _preprocess(gray, None, sigma, "constant", 0.0)
    jsobel = ndi.sobel(smoothed, axis=1)
    isobel = ndi.sobel(smoothed, axis=0)
    mag = isobel * isobel
    mag += jsobel * jsobel
    np.sqrt(mag, out=mag)
    low = _canny_cy._nonmaximum_suppression_bilinear(isobel, jsobel, mag, eroded, 0.1)
    return smoothed, isobel, jsobel, mag, low


def compare_canny(gray, sigma, name):
    np.save(tmp + "/gray.npy", gray)
    os.makedirs(tmp + "/c", exist_ok=True)
    run("canny", tmp + "/gray.npy", sigma, tmp + "/c")
    sm, i, j, m, low = canny_stages(gray, sigma)
    ref_edges = feat.canny(gray, sigma=sigma)
    res = {}
    for key, refv in [("smoothed", sm), ("isobel", i), ("jsobel", j), ("magnitude", m), ("suppressed", low)]:
        res[key] = int((np.load(tmp + "/c/%s.npy" % key) != refv).sum())
    res["edges"] = int((np.load(tmp + "/c/edges.npy").astype(bool) != ref_edges).sum())
    report("canny %s sigma %s" % (name, sigma), all(v == 0 for v in res.values()), str(res))


for sigma in [1.0, 2.0, 2.5]:
    g = rng.random((90, 120)).astype(np.float32)
    g = ndi.uniform_filter(g, 3).astype(np.float32)
    compare_canny(g, sigma, "random smooth")

if frames:
    import glob
    files = sorted(glob.glob(frames + "/frames/*/*_gray.npy"))[:40]
    for f in files:
        grays = np.load(f)
        for k in range(grays.shape[0]):
            compare_canny(np.ascontiguousarray(grays[k]), 2.0, "%s b%d" % (os.path.relpath(f, frames), k))

# filter2D box
for ksz in [3, 5]:
    plane = rng.random((57, 83))
    np.save(tmp + "/plane.npy", plane)
    run("box", tmp + "/plane.npy", ksz, tmp + "/b.npy")
    ref = cv.filter2D(plane, -1, np.ones((ksz, ksz)) / (ksz * ksz))
    o = np.load(tmp + "/b.npy")
    report("cv.filter2D box %d" % ksz, np.array_equal(o, ref), "%d of %d differ" % ((o != ref).sum(), o.size))

# cvtColor
for S in [37, 41, 400, 301, 1, 3, 9]:
    x = rng.random((S, S, 3)).astype(np.float32)
    x[:, :, 1] = x[:, :, 2]
    np.save(tmp + "/px.npy", x.reshape(-1, 3))
    run("gray", tmp + "/px.npy", tmp + "/gr.npy")
    ref = cv.cvtColor(x, cv.COLOR_RGB2GRAY).ravel()
    o = np.load(tmp + "/gr.npy")
    report("cv.cvtColor gray S=%d" % S, np.array_equal(o, ref), "%d of %d differ" % ((o != ref).sum(), o.size))

print("failures:", fails)
sys.exit(1 if fails else 0)
