"""Equivalence cases for hicexplorer.lib.cnb (apieq.py).

The arguments are the shapes the callers use: hicDetectLoops passes the
obs/exp values of one genomic distance as a float64 array with the fitted
size and prob as scalars, chicSignificantInteractions and lib.viewpoint
pass Python floats of a background model. The background model values are
the ones of the repository's cHi-C test data.

Deviation: pdf evaluates np.exp and np.log with libm, numpy 1.26 with its
own AVX512F routines, which differ by one ulp on about a quarter of the
arguments. exp() turns an absolute error of its argument into a relative
error of the result, and the argument is a sum of log-gamma terms of a few
hundred, so the pdf cases allow rtol 1e-12. cdf is exact.
"""
import os

import numpy as np

from apieq import Case

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
C = 'hicexplorer.lib.cnb.'

PDF_NOTE = 'np.exp/np.log of numpy 1.26 (AVX512F) differ from libm by one ulp'


def background_model():
    """(relative distance, size, prob) rows of the cHi-C background model."""
    rows = []
    with open(os.path.join(TEST_DATA, 'cHi-C', 'background.txt')) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t')
            try:
                rows.append((int(fields[0]), float(fields[1]), float(fields[2])))
            except ValueError:
                continue
    return rows


def model_arrays(limit=200):
    rows = background_model()[:limit]
    size = np.array([r[1] for r in rows])
    prob = np.array([r[2] for r in rows])
    return size, prob


def obs_exp_values(seed=3, n=500):
    rng = np.random.default_rng(seed)
    return np.round(rng.gamma(2.0, 3.0, n), 3)


def per_element():
    size, prob = model_arrays()
    x = np.linspace(0, 40, len(size))
    return (x, size, prob)


CASES = []
for name in ('pdf', 'cdf'):
    tol = dict(rtol=1e-12, note=PDF_NOTE) if name == 'pdf' else {}
    CASES += [
        Case(name + '/array_scalar_params', C + name, args=lambda: (obs_exp_values(), 2.5, 0.3), **tol),
        Case(name + '/per_element_model', C + name, args=per_element, **tol),
        Case(name + '/python_floats', C + name, args=(3.0, 4.2, 0.55), **tol),
        Case(name + '/python_ints', C + name, args=(3, 4, 0.5), **tol),
        Case(name + '/numpy_scalars', C + name, args=lambda: (np.float64(7.5), np.float64(1.3), np.float64(0.2)), **tol),
        Case(name + '/zero_d_array', C + name, args=lambda: (np.array(2.0), 1.5, 0.4), **tol),
        Case(name + '/int_array', C + name, args=lambda: (np.arange(20, dtype=np.int64), 3.0, 0.25), **tol),
        Case(name + '/int32_array', C + name, args=lambda: (np.arange(20, dtype=np.int32), 3.0, 0.25), **tol),
        Case(name + '/list', C + name, args=([0.0, 1.0, 2.5, 10.0], 2.0, 0.5), **tol),
        Case(name + '/tuple', C + name, args=((0.0, 1.0, 2.5), 2.0, 0.5), **tol),
        Case(name + '/list_params', C + name, args=lambda: (np.array([0.0, 1.0, 2.5]), [2.0, 3.0, 4.0], [0.5, 0.25, 0.1]), **tol),
        Case(name + '/none', C + name, args=(None, 2.0, 0.5)),
        Case(name + '/two_d_broadcast', C + name,
             args=lambda: (np.arange(12, dtype=np.float64).reshape(3, 4), np.array([[1.0], [2.0], [3.0]]),
                           np.array([0.1, 0.3, 0.5, 0.7])), **tol),
        Case(name + '/edge_probabilities', C + name,
             args=lambda: (np.array([0.0, 1.0, 5.0, 50.0]), 2.0, np.array([0.0, 1.0, 1e-300, 1 - 1e-16])), **tol),
        Case(name + '/nan_inf', C + name,
             args=lambda: (np.array([np.nan, np.inf, -1.0, 2.0]), np.array([1.0, 1.0, 1.0, np.nan]), 0.5), **tol),
        Case(name + '/invalid_params', C + name,
             args=lambda: (np.array([1.0, 2.0, 3.0]), np.array([-1.0, 0.0, 2.0]), np.array([0.5, 1.5, -0.1])), **tol),
        Case(name + '/large_values', C + name,
             args=lambda: (np.array([1e3, 1e5, 1e8]), np.array([50.0, 1e4, 3.0]), 0.999), **tol),
        Case(name + '/empty', C + name, args=lambda: (np.array([], dtype=np.float64), 2.0, 0.5)),
        Case(name + '/broadcast_error', C + name, args=lambda: (np.ones(3), np.ones(2), 0.5)),
    ]
CASES += [
    Case('cdf/broadcast_error_three', C + 'cdf', args=lambda: (np.ones(3), np.ones(2), np.ones(2))),
    Case('pdf/broadcast_error_prob', C + 'pdf', args=lambda: (np.ones(3), 1.0, np.ones(2))),
]
