#!/usr/bin/env python
"""Write the oracle fixtures for the k-means, agglomerative and RandomState tests.

Run with the HiCExplorer 3.7.6 reference environment (numpy 1.26.4,
scikit-learn 1.3.2):

    ~/miniconda3/envs/hicx-ref-3.7.6/bin/python cpp/scripts/sklearn_oracle_fixture.py \
        cpp/tests/fixtures/sklearn

kmeans_oracle.bin (little endian):
    int32 cases; per case: int32 n, d, k, seed; float64 X[n*d];
    int64 kmeans_labels[n]; float64 kmeans_centers[k*d]; int64 ward_labels[n]
random_state_oracle.bin:
    int32 seeds; per seed: uint32 seed; float64 random_sample[1500];
    float64 uniform[7]; int64 choice[40] from p = the fixed vector below,
    each block drawn from a fresh RandomState(seed).
einsum_oracle.bin:
    int32 rows; per row: int32 length; float64 x[length]; float64 einsum('ij,ij->i')
"""
import os
import struct
import sys
import warnings

# Two OpenMP threads: scikit-learn then forms the centre sums of the 256-sample
# blocks of the assignment step per thread and adds the threads' sums, which for
# up to 512 samples (two blocks) is the same as adding the blocks in order. With
# more samples the bits depend on the thread count, so the C++ test compares
# centres of larger datasets only to rounding error.
os.environ["OMP_NUM_THREADS"] = "2"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

import numpy as np
from sklearn.cluster import AgglomerativeClustering, KMeans

CHOICE_P = [0.05, 0.2, 0.1, 0.25, 0.15, 0.05, 0.2]


def datasets(rng):
    cases = []
    sizes = [5, 8, 12, 20, 33, 50, 64, 100, 150, 255, 256, 257, 300, 511, 512, 513, 700]
    for index in range(96):
        n = sizes[index % len(sizes)]
        d = int(rng.choice([1, 2, 3, 4, 5, 7, 8, 9, 13, 16, 31, 40]))
        kind = index % 4
        if kind == 0:
            X = rng.normal(size=(n, d))
        elif kind == 1:
            # clustered blobs
            centres = rng.normal(scale=5.0, size=(4, d))
            X = centres[rng.integers(0, 4, size=n)] + rng.normal(size=(n, d))
        elif kind == 2:
            # integer valued with many duplicate points
            X = rng.integers(0, 3, size=(n, d)).astype(float)
        else:
            # contact-profile like: non-negative, skewed
            X = rng.gamma(0.5, 2.0, size=(n, d))
        k = int(rng.integers(2, 11))
        k = min(k, n)
        seed = 0 if index % 8 else int(rng.integers(1, 2**31))
        cases.append((X, k, seed))
    # exact duplicates of a few rows, and a constant matrix
    base = rng.normal(size=(6, 4))
    cases.append((np.repeat(base, 7, axis=0), 5, 0))
    cases.append((np.repeat(base, 40, axis=0), 6, 0))
    cases.append((np.ones((20, 3)), 3, 0))
    cases.append((rng.normal(size=(30, 3)), 1, 0))
    cases.append((rng.normal(size=(10, 2)), 10, 0))
    # many clusters, for the label order of the dendrogram cut (heap layout)
    cases.append((rng.normal(size=(300, 5)), 25, 0))
    cases.append((rng.integers(0, 4, size=(200, 3)).astype(float), 30, 0))
    cases.append((rng.normal(size=(120, 2)), 60, 0))
    return cases


def write_kmeans(path):
    rng = np.random.default_rng(20261030)
    cases = datasets(rng)
    with open(path, "wb") as out, warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out.write(struct.pack("<i", len(cases)))
        for X, k, seed in cases:
            X = np.ascontiguousarray(X, dtype=np.float64)
            n, d = X.shape
            fit = KMeans(n_clusters=k, random_state=seed).fit(X)
            ward = AgglomerativeClustering(n_clusters=k).fit(X).labels_
            out.write(struct.pack("<iiii", n, d, k, seed))
            out.write(X.astype("<f8").tobytes())
            out.write(fit.labels_.astype("<i8").tobytes())
            out.write(fit.cluster_centers_.astype("<f8").tobytes())
            out.write(ward.astype("<i8").tobytes())
    print(path, len(cases), "cases")


def write_random(path):
    seeds = [0, 1, 2, 3, 5, 42, 123, 1234, 65535, 65536, 2**31 - 1, 2**31, 2**32 - 1]
    gen = np.random.default_rng(7)
    seeds += [int(s) for s in gen.integers(0, 2**32, size=40)]
    with open(path, "wb") as out:
        out.write(struct.pack("<i", len(seeds)))
        for seed in seeds:
            state = np.random.RandomState(seed)
            out.write(struct.pack("<I", seed))
            out.write(np.array([state.random_sample() for _ in range(1500)], "<f8").tobytes())
            out.write(np.asarray(state.uniform(size=7), "<f8").tobytes())
            state = np.random.RandomState(seed)
            choices = [int(state.choice(len(CHOICE_P), p=CHOICE_P)) for _ in range(40)]
            out.write(np.array(choices, "<i8").tobytes())
    print(path, len(seeds), "seeds")


def write_einsum(path):
    rng = np.random.default_rng(99)
    rows = []
    for length in list(range(1, 131)) + [257, 1000, 4097]:
        for scale in (1.0, 1e-3):
            rows.append(rng.normal(scale=scale, size=length) * rng.gamma(1.0, 3.0, size=length))
    with open(path, "wb") as out:
        out.write(struct.pack("<i", len(rows)))
        for x in rows:
            x = np.ascontiguousarray(x.reshape(1, -1))
            out.write(struct.pack("<i", x.shape[1]))
            out.write(x.astype("<f8").tobytes())
            out.write(np.einsum("ij,ij->i", x, x).astype("<f8").tobytes())
    print(path, len(rows), "rows")


if __name__ == "__main__":
    target = sys.argv[1]
    os.makedirs(target, exist_ok=True)
    write_kmeans(os.path.join(target, "kmeans_oracle.bin"))
    write_random(os.path.join(target, "random_state_oracle.bin"))
    write_einsum(os.path.join(target, "einsum_oracle.bin"))
