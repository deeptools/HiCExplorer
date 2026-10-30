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

# gaussian kernels: the default sigma 2 must be numpy's exactly; other sigmas
# use the C library's exp and may differ in the last place (numpy's exp on
# AVX512 CPUs is Intel's SVML), which is reported but not a failure
from scipy.ndimage._filters import _gaussian_kernel1d
other = []
for sigma in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.3, 4.0, 5.0, 6.0, 7.5]:
    run("kernel", sigma, tmp + "/k.npy")
    r = int(4.0 * sigma + 0.5)
    same = np.array_equal(np.load(tmp + "/k.npy"), _gaussian_kernel1d(sigma, 0, r))
    if sigma == 2.0:
        report("gaussian kernel1d sigma 2 (default)", same)
    elif not same:
        other.append(sigma)
print("info   gaussian kernels of other sigmas differing in the last bit from numpy's: %s" % other)

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


def canny_stages(gray, sigma, low=0.1):
    # skimage.feature._canny.canny, stage by stage
    from skimage.feature._canny import _preprocess
    smoothed, eroded = _preprocess(gray, None, sigma, "constant", 0.0)
    jsobel = ndi.sobel(smoothed, axis=1)
    isobel = ndi.sobel(smoothed, axis=0)
    mag = isobel * isobel
    mag += jsobel * jsobel
    np.sqrt(mag, out=mag)
    low = _canny_cy._nonmaximum_suppression_bilinear(isobel, jsobel, mag, eroded, low)
    return smoothed, isobel, jsobel, mag, low


def compare_canny(gray, sigma, name, low_t=None, high_t=None, quiet_ok=False):
    np.save(tmp + "/gray.npy", gray)
    os.makedirs(tmp + "/c", exist_ok=True)
    if low_t is None:
        run("canny", tmp + "/gray.npy", sigma, tmp + "/c")
        sm, i, j, m, low = canny_stages(gray, sigma)
        ref_edges = feat.canny(gray, sigma=sigma)
    else:
        run("canny", tmp + "/gray.npy", sigma, tmp + "/c", repr(low_t), repr(high_t))
        sm, i, j, m, low = canny_stages(gray, sigma, low_t)
        ref_edges = feat.canny(gray, sigma=sigma, low_threshold=low_t, high_threshold=high_t)
    res = {}
    for key, refv in [("smoothed", sm), ("isobel", i), ("jsobel", j), ("magnitude", m), ("suppressed", low)]:
        res[key] = int((np.load(tmp + "/c/%s.npy" % key) != refv).sum())
    res["edges"] = int((np.load(tmp + "/c/edges.npy").astype(bool) != ref_edges).sum())
    if quiet_ok and all(v == 0 for v in res.values()):
        return True
    report("canny %s sigma %s thresholds %s" % (name, sigma, (low_t, high_t)), all(v == 0 for v in res.values()), str(res))
    return all(v == 0 for v in res.values())


for sigma in [1.0, 2.0, 2.5]:
    g = rng.random((90, 120)).astype(np.float32)
    g = ndi.uniform_filter(g, 3).astype(np.float32)
    compare_canny(g, sigma, "random smooth")

if frames:
    import glob
    files = sorted(glob.glob(frames + "/frames/*/*_gray.npy"))
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


# ---------------------------------------------------------------------------
# Randomised comparisons added with the own-code rewrite of stripes_numpy.cpp

def check(name, ok, detail=""):
    report(name, bool(ok), detail)


# random.Random: many seeds (including two-word seeds), random() and choices
seeds = [0, 1, 2, 3, 2**31, 2**32 - 1, 2**32, 2**32 + 1, 2**63, 2**64 - 1] + [int(x) for x in rng.integers(0, 2**63, 90, dtype=np.int64)] \
    + [int(x) for x in rng.integers(0, 2**32, 100, dtype=np.int64)]
bad = 0
for seed in seeds:
    n = 1500
    run("random", seed, n, tmp + "/r.npy")
    r = random.Random(seed)
    if not np.array_equal(np.load(tmp + "/r.npy"), np.array([r.random() for _ in range(n)])):
        bad += 1
check("random.Random.random() for %d seeds x 1500 draws (the state is regenerated twice)" % len(seeds), bad == 0, "%d seeds differ" % bad)
bad = 0
for seed in seeds[:100]:
    for size in [1, 2, 3, 10, 1000, 123457]:
        run("choices", seed, size, 400, tmp + "/c.npy")
        r = random.Random(seed)
        ref = np.array([r.choices(range(size), k=1)[0] for _ in range(400)])
        bad += int(not np.array_equal(np.load(tmp + "/c.npy"), ref))
check("random.choices(population) for 100 seeds x 6 sizes", bad == 0, "%d differ" % bad)
bad = 0
for seed in seeds[:100]:
    n = int(rng.integers(1, 40))
    w = rng.random(n) * (rng.random() * 10)
    if rng.random() < 0.3:
        w = np.round(w)  # ties in the cumulative weights
    cum = np.cumsum(w)
    if cum[-1] <= 0:
        continue
    np.save(tmp + "/cw.npy", cum)
    run("wchoices", seed, tmp + "/cw.npy", 300, tmp + "/c.npy")
    r = random.Random(seed)
    ref = np.array([r.choices(range(n), cum_weights=list(cum), k=1)[0] for _ in range(300)])
    bad += int(not np.array_equal(np.load(tmp + "/c.npy"), ref))
check("random.choices(cum_weights) for 100 seeds", bad == 0, "%d differ" % bad)

# median / nanmedian / quantile / round
bad = 0
for trial in range(1500):
    n = int(rng.integers(1, 60))
    x = rng.integers(-5, 6, n).astype(float) / (1 + int(rng.integers(0, 3))) if trial % 3 == 0 else rng.normal(size=n) * 10.0 ** rng.integers(-3, 6)
    if trial % 7 == 0:
        x[rng.integers(0, n)] = -0.0
    if trial % 11 == 0:
        x[rng.integers(0, n)] = np.inf
    np.save(tmp + "/x.npy", x)
    run("median", tmp + "/x.npy", tmp + "/o.npy")
    bad += int(np.load(tmp + "/o.npy").tobytes() != np.array([np.median(x)]).tobytes() and not np.isnan(np.median(x)))
    y = x.copy()
    y[rng.random(n) < 0.3] = np.nan
    np.save(tmp + "/x.npy", y)
    run("nanmedian", tmp + "/x.npy", tmp + "/o.npy")
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = np.nanmedian(y)
    bad += int(np.load(tmp + "/o.npy").tobytes() != np.array([ref]).tobytes() and not np.isnan(ref))
check("np.median / np.nanmedian of 1500 random arrays (ties, -0.0, inf, NaN)", bad == 0, "%d differ" % bad)
np.save(tmp + "/x.npy", np.array([]))
run("median", tmp + "/x.npy", tmp + "/o.npy")
check("np.median of an empty array is NaN", np.isnan(np.load(tmp + "/o.npy")[0]))
xr = np.concatenate([rng.normal(size=100000) * 100, [0.5, -0.5, 1.5, 2.5, -1.5, 2.5000000000000004, 1e16, -1e16, 0.49999999999999994, -0.0]])
np.save(tmp + "/x.npy", xr)
run("round", tmp + "/x.npy", tmp + "/o.npy")
check("np.round (half to even)", np.array_equal(np.load(tmp + "/o.npy"), np.round(xr)))
bad = 0
for trial in range(600):
    n = int(rng.integers(1, 300))
    x = np.round(rng.normal(size=n) * 5, int(rng.integers(0, 3))) if trial % 2 else rng.normal(size=n)
    qs = np.concatenate([rng.random(40), [0.0, 1.0, 0.5, 0.25, 0.75, 0.9, 0.99, 0.999, 0.05, 1 / 3]])
    np.save(tmp + "/x.npy", x)
    np.save(tmp + "/q.npy", qs)
    run("quantile", tmp + "/x.npy", tmp + "/q.npy", tmp + "/o.npy")
    bad += int(not np.array_equal(np.load(tmp + "/o.npy"), np.quantile(x, qs)))
check("np.quantile (linear) at 50 q for 600 arrays (ties, n = 1..300)", bad == 0, "%d differ" % bad)

# block sums at the buffer boundaries
q = []
for R, C in [(91, 91), (92, 90), (90, 92), (8192, 1), (1, 8192), (128, 65), (65, 128), (129, 64), (64, 129), (200, 41), (41, 200)]:
    q.append((0, 0, R, C))
big2 = rng.random((300, 8300)) * 3
np.save(tmp + "/big2.npy", big2)
q2 = []
for _ in range(3000):
    R = int(rng.integers(1, 120)); C = int(rng.integers(1, 120))
    q2.append((int(rng.integers(0, 300 - R + 1)), int(rng.integers(0, 8300 - C + 1)), R, C))
q2 = np.array(q2 + [(0, 0, 100, 82), (0, 0, 82, 100), (3, 5, 128, 64)], dtype=np.int64)
q2 = q2[(q2[:, 2] <= 300) & (q2[:, 0] + q2[:, 2] <= 300)]
np.save(tmp + "/q2.npy", q2)
run("blocksums", tmp + "/big2.npy", tmp + "/q2.npy", tmp + "/o.npy")
o = np.load(tmp + "/o.npy")
ref = np.array([[np.sum(big2[a:a + R, b:b + C]), np.mean(big2[a:a + R, b:b + C])] for a, b, R, C in q2])
check("np.sum / np.mean of 2-D views with wide strides and more than 8192 elements", np.array_equal(o, ref), "%d of %d differ" % ((o != ref).any(axis=1).sum(), len(q2)))

# the sign of a zero sum: blocks made of -0.0 and 0.0 (numpy 1.26 sums from +0.0)
zz = np.where(rng.random((120, 9000)) < 0.5, -0.0, 0.0)
zz[:, ::3] = -0.0
np.save(tmp + "/zz.npy", zz)
qz = np.array([(0, 0, 1, 1), (0, 0, 1, 5), (0, 3, 2, 3), (0, 0, 7, 1), (0, 0, 12, 700), (5, 7, 100, 90), (0, 0, 1, 9000), (0, 0, 120, 9000)], dtype=np.int64)
qz = np.concatenate([qz, np.array([(int(rng.integers(0, 100)), int(rng.integers(0, 8000)), int(rng.integers(1, 20)), int(rng.integers(1, 900))) for _ in range(300)], dtype=np.int64)])
np.save(tmp + "/qz.npy", qz)
run("blocksums", tmp + "/zz.npy", tmp + "/qz.npy", tmp + "/o.npy")
o = np.load(tmp + "/o.npy")
ref = np.array([[np.sum(zz[a:a + R, b:b + C]), np.mean(zz[a:a + R, b:b + C])] for a, b, R, C in qz])
check("np.sum / np.mean of blocks of signed zeros: same bits incl. the sign", o.tobytes() == ref.tobytes(), "%d differ" % (o.view(np.int64) != ref.view(np.int64)).sum())

# canny: random images with ties, several sigmas and thresholds, small and degenerate shapes
def random_image(kind):
    h = int(rng.integers(1, 90)); w = int(rng.integers(1, 90))
    if kind == 0:
        g = rng.random((h, w))
    elif kind == 1:  # values on a coarse grid: exact ties in every stage
        g = rng.integers(0, 5, (h, w)) / 4.0
    elif kind == 2:  # blocks and lines
        g = np.zeros((h, w))
        for _ in range(int(rng.integers(1, 6))):
            a, b = rng.integers(0, h), rng.integers(0, w)
            g[a:a + int(rng.integers(1, 30)), b:b + int(rng.integers(1, 30))] = rng.random()
        g += 0.01 * rng.random((h, w)) * (rng.random() < 0.5)
    elif kind == 3:  # constant
        g = np.full((h, w), float(rng.random()))
    else:  # smooth gradient plus diagonal structure (exact |gi| == |gj| ties)
        yy, xx = np.mgrid[0:h, 0:w]
        g = ((yy + xx) % int(rng.integers(3, 12))) / 12.0 + (yy > xx) * 0.3
    return np.ascontiguousarray(g.astype(np.float32))


bad = 0
total = 0
details = []
for trial in range(900):
    g = random_image(trial % 5)
    sigma = float(rng.choice([0.5, 0.8, 1.0, 1.5, 2.0, 2.0, 2.5, 3.0, 4.0]))
    if trial % 3 == 0:
        lo, hi = None, None
    else:
        lo = float(np.round(rng.random() * 0.3, int(rng.integers(1, 7))))
        hi = float(lo + rng.random() * 0.4)
        if rng.random() < 0.2:
            hi = lo  # equal thresholds
        lo, hi = max(lo, 1e-6), max(hi, 1e-6)
    total += 1
    if not compare_canny(g, sigma, "random kind %d %s" % (trial % 5, g.shape), lo, hi, quiet_ok=True):
        bad += 1
        details.append((trial % 5, g.shape, sigma, lo, hi))
check("canny stages and edges for %d random images (ties, several sigmas and thresholds, shapes 1x1 to 90x90)" % total, bad == 0, "%d differ %s" % (bad, details[:3]))

# box filter: kernel sizes, shapes smaller than the kernel, ties
bad = 0
for trial in range(300):
    h = int(rng.integers(1, 40)); w = int(rng.integers(1, 40)); k = int(rng.choice([1, 3, 5, 7]))
    plane = rng.integers(0, 6, (h, w)).astype(float) / 3 if trial % 2 else rng.random((h, w)) * 100
    np.save(tmp + "/plane.npy", plane)
    run("box", tmp + "/plane.npy", k, tmp + "/b.npy")
    ref = cv.filter2D(plane, -1, np.ones((k, k)) / (k * k))
    bad += int(not np.array_equal(np.load(tmp + "/b.npy"), ref))
check("cv.filter2D normalised box (k = 1, 3, 5, 7, shapes from 1x1, ties)", bad == 0, "%d of 300 differ" % bad)

# cvtColor: float32 random pixels with ties and extremes; all 2^24 uint8 colours
px = rng.random((200000, 3)).astype(np.float32)
px[:20000] = np.round(px[:20000] * 8) / 8
px[20000:20010] = [[1, 1, 1], [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [3e38, 3e38, 3e38], [-1, 2, 5], [1e-40, 1e-40, 1e-40], [np.nan, 1, 1], [0.5, 0.5, 0.5]]
np.save(tmp + "/px.npy", px)
run("gray", tmp + "/px.npy", tmp + "/gr.npy")
ref = cv.cvtColor(px.reshape(-1, 1, 3), cv.COLOR_RGB2GRAY).ravel()
o = np.load(tmp + "/gr.npy")
check("cv.cvtColor RGB2GRAY float32, 200000 pixels incl. extremes", np.array_equal(o, ref, equal_nan=True), "%d differ" % (~((o == ref) | (np.isnan(o) & np.isnan(ref)))).sum())
allc = np.arange(2**24, dtype=np.uint32)
rgb = np.stack([(allc >> 16) & 255, (allc >> 8) & 255, allc & 255], axis=1).astype(np.uint8)
np.save(tmp + "/px8.npy", rgb)
run("gray8", tmp + "/px8.npy", tmp + "/gr8.npy")
ref8 = cv.cvtColor(rgb.reshape(-1, 4096, 3), cv.COLOR_RGB2GRAY).ravel()
check("cv.cvtColor RGB2GRAY uint8, all 16777216 colours", np.array_equal(np.load(tmp + "/gr8.npy"), ref8), "%d differ" % (np.load(tmp + "/gr8.npy") != ref8).sum())

# signal.convolve2d same
import scipy.signal as sig
bad = 0
for trial in range(300):
    h = int(rng.integers(1, 30)); w = int(rng.integers(1, 30))
    m = (rng.random((h, w)) < rng.random()).astype(np.uint8)
    kern = rng.integers(-3, 4, (3, 3)).astype(np.int32)
    np.save(tmp + "/m.npy", m)
    np.save(tmp + "/k.npy", kern)
    run("convolve", tmp + "/m.npy", tmp + "/k.npy", tmp + "/o.npy")
    bad += int(not np.array_equal(np.load(tmp + "/o.npy"), sig.convolve2d(m, kern, mode="same")))
check("scipy.signal.convolve2d mode same, 3x3 kernels, shapes from 1x1", bad == 0, "%d differ" % bad)

print("failures:", fails)
sys.exit(1 if fails else 0)
