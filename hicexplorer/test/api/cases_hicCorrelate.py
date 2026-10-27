"""Equivalence cases for hicexplorer.hicCorrelate (apieq.py).

get_vectors is compared on csr, csc and coo matrices (sorted and unsorted
indices, duplicates, explicit zeros, cancelling entries, every dtype of the
upcast table), including the order of the returned entries.

plot_correlation draws; the driver compares what the drawing is given
(the linkage handed to scipy's dendrogram, the colour range handed to
matshow) and the size of the written figure, the exceptions, and the file
that is written. Pixel content depends on the matplotlib version of the
interpreter and is not compared.

Deviations from 3.7.6 (hicexplorer/hicCorrelate.py, cpp/python/src/api/hicCorrelate.cpp):
  * get_vectors takes csr/csc/coo matrices of dtype int32, int64, float32,
    float64; other formats and dtypes raise TypeError.
  * plot_correlation takes a two dimensional ndarray of a real dtype.
"""
import os
import struct

import numpy as np
import scipy.sparse as sp

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
H = 'hicexplorer.hicCorrelate.'


# ---------------------------------------------------------------------------
# inputs

def real_block(window=0, scale=1.0, dtype=np.float64):
    """A 260 x 260 window of the Dekker et al. Hi-C counts (test_data/dekker.txt.gz,
    a dense text matrix, read with the standard library and numpy) as csr, its
    values multiplied by `scale`."""
    import gzip
    import itertools
    first = 900 + 400 * window
    rows = []
    with gzip.open(os.path.join(TEST_DATA, 'dekker.txt.gz'), 'rt') as handle:
        for line in itertools.islice(handle, first + 1, first + 261):
            fields = line.rstrip('\n').split('\t')
            rows.append([float(x) for x in fields[first + 1:first + 261]])
    block = sp.csr_matrix(np.array(rows))
    if dtype in (np.int32, np.int64):
        block.data = np.rint(block.data * scale).astype(dtype)
    else:
        block.data = (block.data * scale).astype(dtype)
    return block


def random_matrix(shape, density, seed, dtype, fmt='csr'):
    rng = np.random.default_rng(seed)
    m = sp.random(shape[0], shape[1], density, random_state=int(seed), format=fmt)
    if dtype in (np.int32, np.int64):
        m.data = rng.integers(-5, 20, size=m.data.shape).astype(dtype)
    else:
        m.data = (rng.normal(3.0, 2.0, size=m.data.shape)).astype(dtype)
    return m


def raw(fmt, data, indices, indptr, shape, dtype=np.float64, index_dtype=np.int32):
    data = np.array(data, dtype=dtype)
    indices = np.array(indices, dtype=index_dtype)
    indptr = np.array(indptr, dtype=index_dtype)
    if fmt == 'csr':
        return sp.csr_matrix((data, indices, indptr), shape=shape)
    return sp.csc_matrix((data, indices, indptr), shape=shape)


def unsorted_csr(dtype=np.float64):
    # row 0 stores columns 2, 0, 1; row 1 stores 3, 1; row 2 is empty
    return raw('csr', [3, 1, 2, 7, 5], [2, 0, 1, 3, 1], [0, 3, 5, 5], (3, 4), dtype)


def duplicate_csr(dtype=np.float64):
    # row 0 stores column 1 twice and column 0; row 1 column 3 twice
    return raw('csr', [2, 3, 1, 4, 6], [1, 1, 0, 3, 3], [0, 3, 5, 5], (3, 4), dtype)


def sorted_csr(dtype=np.float64):
    return raw('csr', [1, 5, 2, 8, 6, 3], [0, 3, 1, 1, 3, 2], [0, 2, 4, 6], (3, 4), dtype)


def zeros_csr(dtype=np.float64):
    # explicit zeros and cancelling partners of sorted_csr
    return raw('csr', [-1, 0, -2, 4, 0, 3], [0, 3, 1, 2, 3, 0], [0, 2, 4, 6], (3, 4), dtype)


def coo(rows, cols, data, shape, dtype=np.float64):
    return sp.coo_matrix((np.array(data, dtype=dtype), (np.array(rows), np.array(cols))),
                         shape=shape)


def coo_duplicates(dtype=np.float64):
    return coo([0, 0, 1, 0, 2, 2], [1, 1, 3, 0, 2, 2], [2, 3, 1, 4, 5, 1], (3, 4), dtype)


# ---------------------------------------------------------------------------
# drivers

def PNG_size(path):
    with open(path, 'rb') as handle:
        head = handle.read(24)
    return struct.unpack('>II', head[16:24])


def plot(target, matrix, labels, name, **kwargs):
    """Calls plot_correlation with the figure file in the case directory;
    returns what the drawing was given and the figure's size."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scipy.cluster.hierarchy as sch
    from matplotlib.axes import Axes

    seen = {}
    dendrogram = sch.dendrogram
    matshow = Axes.matshow

    def spy_dendrogram(z, *args, **kw):
        seen['linkage'] = np.array(z)
        return dendrogram(z, *args, **kw)

    def spy_matshow(self, a, *args, **kw):
        seen['vmax'] = kw.get('vmax')
        seen['vmin'] = kw.get('vmin')
        seen['matrix'] = np.array(a)
        return matshow(self, a, *args, **kw)

    sch.dendrogram = spy_dendrogram
    Axes.matshow = spy_matshow
    try:
        result = target(matrix, labels, name, **kwargs)
    finally:
        sch.dendrogram = dendrogram
        Axes.matshow = matshow
        plt.close('all')
    return {'result': result, 'seen': seen, 'size': PNG_size(name),
            'argument_matrix': np.array(matrix)}


def correlation(n, seed, negative=False):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(n, 12))
    m = np.corrcoef(a)
    if not negative:
        m = (m + 1) / 2
    return m


CASES = []

# --- get_vectors --------------------------------------------------------------
for dtype in (np.float64, np.float32, np.int64, np.int32):
    tag = np.dtype(dtype).name
    CASES += [
        Case('get_vectors/random_csr_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (random_matrix((30, 30), 0.3, 1, dtype),
                                       random_matrix((30, 30), 0.3, 2, dtype))),
        Case('get_vectors/sorted_zeros_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (sorted_csr(dtype), zeros_csr(dtype))),
        Case('get_vectors/unsorted_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (unsorted_csr(dtype), sorted_csr(dtype)[:3, :])),
        Case('get_vectors/unsorted_both_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (unsorted_csr(dtype), unsorted_csr(dtype))),
        Case('get_vectors/duplicates_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (duplicate_csr(dtype), sorted_csr(dtype))),
        Case('get_vectors/csc_csc_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (random_matrix((17, 23), 0.4, 3, dtype, 'csc'),
                                       random_matrix((17, 23), 0.4, 4, dtype, 'csc'))),
        Case('get_vectors/csr_csc_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (random_matrix((17, 23), 0.4, 5, dtype, 'csr'),
                                       random_matrix((17, 23), 0.4, 6, dtype, 'csc'))),
        Case('get_vectors/csc_csr_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (random_matrix((17, 23), 0.4, 7, dtype, 'csc'),
                                       random_matrix((17, 23), 0.4, 8, dtype, 'csr'))),
        Case('get_vectors/coo_csr_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (coo_duplicates(dtype), sorted_csr(dtype))),
        Case('get_vectors/csr_coo_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (sorted_csr(dtype), coo_duplicates(dtype))),
        Case('get_vectors/coo_coo_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (coo_duplicates(dtype), coo_duplicates(dtype))),
        Case('get_vectors/same_matrix_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (sorted_csr(dtype), sorted_csr(dtype))),
        Case('get_vectors/negated_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (sorted_csr(dtype), -sorted_csr(dtype))),
        Case('get_vectors/empty_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (sp.csr_matrix((4, 5), dtype=dtype),
                                       sp.csr_matrix((4, 5), dtype=dtype))),
        Case('get_vectors/one_empty_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (sp.csr_matrix((3, 4), dtype=dtype), sorted_csr(dtype))),
        Case('get_vectors/zero_size_' + tag, H + 'get_vectors',
             args=lambda dtype=dtype: (sp.csr_matrix((0, 0), dtype=dtype),
                                       sp.csr_matrix((0, 0), dtype=dtype))),
    ]

for name, a, b in [('f64_f32', np.float64, np.float32), ('f32_f64', np.float32, np.float64),
                   ('i32_i64', np.int32, np.int64), ('i64_i32', np.int64, np.int32),
                   ('i32_f32', np.int32, np.float32), ('f32_i64', np.float32, np.int64),
                   ('i64_f64', np.int64, np.float64), ('f64_i32', np.float64, np.int32)]:
    CASES.append(Case('get_vectors/mixed_' + name, H + 'get_vectors',
                      args=lambda a=a, b=b: (random_matrix((20, 20), 0.3, 11, a),
                                             random_matrix((20, 20), 0.3, 12, b))))

CASES += [
    Case('get_vectors/real_chr_same', H + 'get_vectors',
         args=lambda: (real_block(0), real_block(0, 0.5))),
    Case('get_vectors/real_chr_other_scale', H + 'get_vectors',
         args=lambda: (real_block(1), real_block(1, 3.0))),
    Case('get_vectors/real_int64', H + 'get_vectors',
         args=lambda: (real_block(0, 1.0, np.int64), real_block(0, 0.25, np.int64))),
    Case('get_vectors/real_float32', H + 'get_vectors',
         args=lambda: (real_block(0, 1.0, np.float32), real_block(0, 0.37, np.float32))),
    Case('get_vectors/int64_index', H + 'get_vectors',
         args=lambda: (raw('csr', [1, 5, 2, 8, 6, 3], [0, 3, 1, 1, 3, 2], [0, 2, 4, 6], (3, 4),
                           np.float64, np.int64), sorted_csr())),
    Case('get_vectors/nan_inf', H + 'get_vectors',
         args=lambda: (raw('csr', [np.nan, np.inf, 2.0], [0, 1, 2], [0, 2, 3], (2, 3)),
                       raw('csr', [1.0, -np.inf, 3.0], [0, 1, 1], [0, 2, 3], (2, 3)))),
    Case('get_vectors/rounding', H + 'get_vectors',
         args=lambda: (raw('csr', [0.1, 1e-17, 1e16], [0, 1, 2], [0, 3], (1, 3)),
                       raw('csr', [0.2, 3.0, 1.0], [0, 1, 2], [0, 3], (1, 3)))),
    Case('get_vectors/int_overflow', H + 'get_vectors',
         args=lambda: (raw('csr', [2 ** 31 - 1, 5], [0, 1], [0, 2], (1, 2), np.int32),
                       raw('csr', [3, 7], [0, 1], [0, 2], (1, 2), np.int32))),
    Case('get_vectors/args_unchanged', H + 'get_vectors',
         args=lambda: (sorted_csr(), unsorted_csr()[:3, :])),
    # exceptions
    Case('get_vectors/shape_mismatch', H + 'get_vectors',
         args=lambda: (sp.csr_matrix((3, 4)), sp.csr_matrix((4, 3)))),
    Case('get_vectors/shape_mismatch_one_dim', H + 'get_vectors',
         args=lambda: (sp.csr_matrix((3, 4)), sp.csr_matrix((3, 5)))),
    Case('get_vectors/no_shape', H + 'get_vectors', args=lambda: ([1, 2], sp.csr_matrix((3, 4)))),
    Case('get_vectors/no_shape_second', H + 'get_vectors',
         args=lambda: (sp.csr_matrix((3, 4)), None)),
]

# --- plot_correlation ---------------------------------------------------------
CASES += [
    Case('plot_correlation/positive', H + 'plot_correlation', driver='plot',
         args=lambda: (correlation(5, 1), list('abcde'), File('p.png')), compare_args=True),
    Case('plot_correlation/negative', H + 'plot_correlation', driver='plot',
         args=lambda: (correlation(6, 2, True), list('abcdef'), File('p.png'))),
    Case('plot_correlation/numbers', H + 'plot_correlation', driver='plot',
         args=lambda: (correlation(4, 3), ['w', 'x', 'y', 'z'], File('p.png')),
         kwargs={'pPlotNumbers': True}),
    Case('plot_correlation/range_given', H + 'plot_correlation', driver='plot',
         args=lambda: (correlation(5, 4), list('abcde'), File('p.png')),
         kwargs={'vmax': 0.9, 'vmin': 0.1, 'colormap': 'viridis'}),
    Case('plot_correlation/vmax_int', H + 'plot_correlation', driver='plot',
         args=lambda: (correlation(4, 5), list('abcd'), File('p.png')),
         kwargs={'vmax': 2, 'vmin': 0}),
    Case('plot_correlation/two_rows', H + 'plot_correlation', driver='plot',
         args=lambda: (np.array([[1.0, 0.3], [0.3, 1.0]]), ['a', 'b'], File('p.png'))),
    Case('plot_correlation/integer_matrix', H + 'plot_correlation', driver='plot',
         args=lambda: (np.array([[3, 1, 2], [1, 3, 0], [2, 0, 3]]), ['a', 'b', 'c'], File('p.png'))),
    Case('plot_correlation/float32', H + 'plot_correlation', driver='plot',
         args=lambda: (correlation(5, 6).astype(np.float32), list('abcde'), File('p.png'))),
    Case('plot_correlation/hollow', H + 'plot_correlation', driver='plot',
         args=lambda: (np.array([[0.0, 0.2, 0.7], [0.2, 0.0, 0.5], [0.7, 0.5, 0.0]]),
                       ['a', 'b', 'c'], File('p.png'))),
    Case('plot_correlation/nan_min', H + 'plot_correlation', driver='plot',
         args=lambda: (np.array([[1.0, np.nan], [np.nan, 1.0]]), ['a', 'b'], File('p.png'))),
    Case('plot_correlation/nan_range_given', H + 'plot_correlation', driver='plot',
         args=lambda: (np.array([[1.0, np.nan], [np.nan, 1.0]]), ['a', 'b'], File('p.png')),
         kwargs={'vmin': 0, 'vmax': 1}),
    Case('plot_correlation/inf', H + 'plot_correlation', driver='plot',
         args=lambda: (np.array([[1.0, np.inf], [np.inf, 1.0]]), ['a', 'b'], File('p.png'))),
    Case('plot_correlation/overflow_distance', H + 'plot_correlation', driver='plot',
         args=lambda: (np.array([[1e308, 0.0], [-1e308, 0.0]]), ['a', 'b'], File('p.png'))),
    Case('plot_correlation/one_row', H + 'plot_correlation', driver='plot',
         args=lambda: (np.array([[1.0]]), ['a'], File('p.png'))),
    Case('plot_correlation/empty', H + 'plot_correlation', driver='plot',
         args=lambda: (np.zeros((0, 0)), [], File('p.png'))),
    Case('plot_correlation/empty_range_given', H + 'plot_correlation', driver='plot',
         args=lambda: (np.zeros((0, 0)), [], File('p.png')), kwargs={'vmin': 0}),
    Case('plot_correlation/no_columns', H + 'plot_correlation', driver='plot',
         args=lambda: (np.zeros((3, 0)), ['a', 'b', 'c'], File('p.png'))),
    Case('plot_correlation/rectangular', H + 'plot_correlation', driver='plot',
         args=lambda: (np.arange(12.0).reshape(3, 4) / 12, ['a', 'b', 'c'], File('p.png'))),
    Case('plot_correlation/no_shape', H + 'plot_correlation', driver='plot',
         args=lambda: ([[1.0, 0.5], [0.5, 1.0]], ['a', 'b'], File('p.png'))),
    Case('plot_correlation/scalar', H + 'plot_correlation', driver='plot',
         args=lambda: (np.float64(1.0), ['a'], File('p.png'))),
    Case('plot_correlation/three_dimensional', H + 'plot_correlation', driver='plot',
         args=lambda: (np.ones((2, 2, 2)), ['a', 'b'], File('p.png'))),
]
