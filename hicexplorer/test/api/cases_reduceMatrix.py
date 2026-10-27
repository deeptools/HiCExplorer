"""Equivalence cases for hicexplorer.reduceMatrix (apieq.py).

The matrices are real chromosome blocks of small_test_matrix_50kb_res.h5
(read with PyTables), symmetrised as hiCMatrix does, merged with the groups
hicMergeMatrixBins.merge_bins and hicMergeTADbins build (consecutive runs,
dropped bins, non consecutive groups), in every dtype and format the
function is called with, plus the docstring matrices and constructed ones
for NaN values, explicit zeros, the early return and the index errors.

The two debug records carry elapsed times and are compared with the time
masked (driver 'masked_times').
"""
import logging
import os

import numpy as np
import scipy.sparse as sp

from apieq import Case

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
T = 'hicexplorer.reduceMatrix.reduce_matrix'


def real_block(chromosome_index=0, dtype=np.float64, fmt='csr', corrected=False):
    import tables
    with tables.open_file(os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5'), 'r') as f:
        m = sp.csr_matrix((f.root.matrix.data.read(), f.root.matrix.indices.read(), f.root.matrix.indptr.read()),
                          shape=tuple(f.root.matrix.shape.read()))
        chroms = [c.decode() if isinstance(c, bytes) else c for c in f.root.intervals.chr_list.read()]
    names = []
    for c in chroms:
        if c not in names:
            names.append(c)
    first = chroms.index(names[chromosome_index])
    last = len(chroms) - chroms[::-1].index(names[chromosome_index])
    block = m[first:last, first:last]
    block = sp.csr_matrix(block + sp.triu(block, 1).T)
    if corrected:
        rng = np.random.default_rng(11)
        bias = rng.uniform(0.3, 3.0, block.shape[0])
        coo = block.tocoo()
        coo.data = coo.data * bias[coo.row] * bias[coo.col] / 7.0
        block = coo.tocsr()
    if dtype in (np.int32, np.int64):
        block.data = np.rint(block.data).astype(dtype)
    else:
        block.data = block.data.astype(dtype)
    return block.asformat(fmt)


def consecutive(n, size):
    return [list(range(i, min(i + size, n))) for i in range(0, n, size)]


def args_consecutive(size, **block):
    def make():
        m = real_block(**block)
        return (m, consecutive(m.shape[0], size))
    return make


def args_dropped(**block):
    def make():
        m = real_block(**block)
        groups = consecutive(m.shape[0], 3)
        return (m, [g for i, g in enumerate(groups) if i % 4 != 1])
    return make


def args_interleaved(**block):
    def make():
        m = real_block(**block)
        n = m.shape[0]
        return (m, [list(range(i, n, 7)) for i in range(7)])
    return make


def args_reversed(**block):
    def make():
        m = real_block(**block)
        return (m, consecutive(m.shape[0], 4)[::-1])
    return make


def docstring_upper():
    A = sp.csr_matrix(np.array([[5, 5, 2, 2, 0], [0, 5, 2, 2, 1], [0, 0, 1, 1, 0], [0, 0, 0, 1, 0],
                                [0, 0, 0, 0, 0]]), dtype=np.int32)
    return (A, [(0, 1), (2, 3), (4,)])


def float_nan():
    A = sp.csr_matrix(np.array([[0.1, 0.1, 0.2, 0.2, np.nan], [0, 0.1, 0.2, 0.2, 1.1], [0, 0, 0.2, 0.2, 0],
                                [0, 0, 0, 0.1, 0], [0, 0, 0, 0, 0]]))
    dia = sp.dia_matrix(([A.diagonal()], [0]), shape=A.shape)
    return (sp.csr_matrix(A + A.T - dia), [(0, 1), (2, 3), (4,)])


def explicit_zeros():
    m = real_block(1)
    m.data[::4] = 0.0
    return (m, consecutive(m.shape[0], 2))


def cancelling():
    # merged cells that sum to exactly zero are eliminated
    A = sp.csr_matrix(np.array([[1.0, -1.0, 2.0], [-1.0, 1.0, 0.0], [2.0, 0.0, 3.0]]))
    return (A, [(0, 1), (2,)])


class _MaskTimes(logging.Filter):
    def filter(self, record):
        message = record.getMessage()
        if message.startswith('time '):
            record.msg = message.split(':')[0] + ': <masked>'
            record.args = ()
        return True


def masked_times(target, *args, **kwargs):
    logger = logging.getLogger('hicexplorer.reduceMatrix')
    mask = _MaskTimes()
    logger.addFilter(mask)
    try:
        return target(*args, **kwargs)
    finally:
        logger.removeFilter(mask)


def identity(target, *args, **kwargs):
    result = target(*args, **kwargs)
    return result is args[0], result


D = dict(driver='masked_times')
CASES = [
    Case('real_merge2', T, args=args_consecutive(2), **D),
    Case('real_merge3_chr2', T, args=args_consecutive(3, chromosome_index=1), **D),
    Case('real_merge5_keep_diagonal', T, args=args_consecutive(5), kwargs={'diagonal': True}, **D),
    Case('real_no_triu', T, args=args_consecutive(2), kwargs={'use_triu': False}, **D),
    Case('real_no_triu_keep_diagonal', T, args=args_consecutive(4, chromosome_index=2),
         kwargs={'use_triu': False, 'diagonal': True}, **D),
    Case('real_int32', T, args=args_consecutive(2, dtype=np.int32), kwargs={'diagonal': True}, **D),
    Case('real_int64', T, args=args_consecutive(3, dtype=np.int64), **D),
    Case('real_float32_corrected', T, args=args_consecutive(3, dtype=np.float32, corrected=True),
         kwargs={'diagonal': True}, **D),
    Case('real_float64_corrected', T, args=args_consecutive(4, corrected=True), kwargs={'diagonal': True}, **D),
    Case('real_float64_corrected_no_triu', T, args=args_consecutive(4, corrected=True),
         kwargs={'use_triu': False}, **D),
    Case('real_csc', T, args=args_consecutive(3, fmt='csc', corrected=True), kwargs={'diagonal': True}, **D),
    Case('real_coo', T, args=args_consecutive(3, fmt='coo', corrected=True), kwargs={'diagonal': True}, **D),
    Case('real_coo_no_triu', T, args=args_consecutive(3, fmt='coo', corrected=True), kwargs={'use_triu': False}, **D),
    Case('real_dropped_groups', T, args=args_dropped(corrected=True), kwargs={'diagonal': True}, **D),
    Case('real_interleaved', T, args=args_interleaved(corrected=True), kwargs={'diagonal': True}, **D),
    Case('real_reversed_groups', T, args=args_reversed(corrected=True), kwargs={'diagonal': True}, **D),
    Case('real_reversed_no_triu', T, args=args_reversed(corrected=True), kwargs={'use_triu': False}, **D),
    Case('docstring_upper', T, args=docstring_upper, kwargs={'diagonal': True}, **D),
    Case('docstring_no_triu', T, args=docstring_upper, kwargs={'diagonal': True, 'use_triu': False}, **D),
    Case('docstring_diagonal_false', T, args=docstring_upper, **D),
    Case('diagonal_zero_is_not_false', T, args=docstring_upper, kwargs={'diagonal': 0}, **D),
    Case('use_triu_truthy', T, args=docstring_upper, kwargs={'use_triu': 1, 'diagonal': True}, **D),
    Case('float_nan', T, args=float_nan, kwargs={'diagonal': True, 'use_triu': False}, **D),
    Case('float_nan_triu', T, args=float_nan, kwargs={'diagonal': True}, **D),
    Case('explicit_zeros', T, args=explicit_zeros, kwargs={'diagonal': True}, **D),
    Case('cancelling', T, args=cancelling, kwargs={'diagonal': True}, **D),
    Case('cancelling_no_triu', T, args=cancelling, kwargs={'use_triu': False}, **D),
    Case('same_size_returns_argument', T, args=lambda: (real_block(0), [[i] for i in range(real_block(0).shape[0])]),
         driver='identity'),
    Case('tuples_and_numpy_ints', T,
         args=lambda: (real_block(0), [tuple(np.arange(i, i + 2)) for i in range(0, 40, 2)]), **D),
    Case('negative_index', T, args=lambda: (real_block(0), [[0, 1], [-1, -2]]), kwargs={'diagonal': True}, **D),
    Case('index_out_of_range', T, args=lambda: (real_block(0), [[0, 1], [10000]]), **D),
    Case('float_index', T, args=lambda: (real_block(0), [[0, 1.0]]), **D),
    Case('empty_groups', T, args=lambda: (real_block(0), []), **D),
    Case('sparse_array', T, args=lambda: (sp.csr_array(real_block(0)), consecutive(real_block(0).shape[0], 2)), **D),
]
