"""Equivalence cases for hicexplorer.iterativeCorrection (apieq.py).

The matrices are real chromosome blocks of small_test_matrix_50kb_res.h5
(read with PyTables so that both interpreters build the same input),
symmetrised as hiCMatrix does, in the dtypes and sparse formats the function
is called with, plus small constructed ones for the NaN replacement, the
symmetry check and the two exit(1) guards.

The verbose log lines that print an estimated time are compared with the
time masked (driver 'masked_times'); exit(1) is compared as the SystemExit
code (driver 'catch_exit').
"""
import logging
import os

import numpy as np
import scipy.sparse as sp

from apieq import Case

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
T = 'hicexplorer.iterativeCorrection.iterativeCorrection'


def real_block(chromosome_index=0, dtype=np.float64, fmt='csr'):
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
    if dtype in (np.int32, np.int64):
        block.data = np.rint(block.data).astype(dtype)
    else:
        block.data = block.data.astype(dtype)
    return block.asformat(fmt)


def with_nan(dtype=np.float64):
    m = real_block(1, dtype)
    m.data[::7] = np.nan
    return m


def asymmetric():
    m = real_block(1).tolil()
    m[0, 5] = m[0, 5] + 1000.0
    return m.tocsr()


def tiny_bin():
    return sp.csr_matrix(np.array([[1e-200, 0.0], [0.0, 1.0]]))


def large_at_end():
    # converges in one pass but the final rescale overflows 1e10
    return sp.csr_matrix(np.array([[1e11, 0.0], [0.0, 1e11]]))


def with_explicit_zeros():
    m = real_block(2)
    m.data[::5] = 0.0
    return m


class _MaskTimes(logging.Filter):
    def filter(self, record):
        message = record.getMessage()
        if 'Estimated time' in message:
            record.msg = message.split('Estimated time')[0] + 'Estimated time <masked>'
            record.args = ()
        return True


def masked_times(target, *args, **kwargs):
    logger = logging.getLogger('hicexplorer.iterativeCorrection')
    mask = _MaskTimes()
    logger.addFilter(mask)
    try:
        return target(*args, **kwargs)
    finally:
        logger.removeFilter(mask)
        logger.setLevel(logging.NOTSET)


def catch_exit(target, *args, **kwargs):
    try:
        return target(*args, **kwargs)
    except SystemExit as error:
        return ('SystemExit', error.code)


CASES = [
    Case('real_default', T, args=lambda: (real_block(0),)),
    Case('real_chr2', T, args=lambda: (real_block(1),)),
    Case('real_converges', T, args=lambda: (real_block(2),), kwargs={'M': 500}),
    Case('real_tolerance', T, args=lambda: (real_block(0),), kwargs={'M': 200, 'tolerance': 1e-3}),
    Case('real_loose_tolerance_first_pass', T, args=lambda: (real_block(0),), kwargs={'tolerance': 10.0}),
    Case('real_int32', T, args=lambda: (real_block(0, np.int32),)),
    Case('real_int64', T, args=lambda: (real_block(1, np.int64),)),
    Case('real_float32', T, args=lambda: (real_block(1, np.float32),)),
    Case('real_coo', T, args=lambda: (real_block(0, fmt='coo'),)),
    Case('real_csc', T, args=lambda: (real_block(0, fmt='csc'),)),
    Case('verbose', T, args=lambda: (real_block(2),), kwargs={'M': 30, 'verbose': True}, driver='masked_times'),
    Case('M_zero', T, args=lambda: (real_block(0),), kwargs={'M': 0}),
    Case('M_positional', T, args=lambda: (real_block(0), None, 7)),
    Case('M_float', T, args=lambda: (real_block(0),), kwargs={'M': 3.0}),
    Case('nan_float64_in_place', T, args=lambda: (with_nan(),)),
    Case('nan_float32_in_place', T, args=lambda: (with_nan(np.float32),)),
    Case('explicit_zeros', T, args=lambda: (with_explicit_zeros(),)),
    Case('asymmetric', T, args=lambda: (asymmetric(),)),
    Case('large_in_loop', T, args=lambda: (tiny_bin(),), driver='catch_exit'),
    Case('large_at_end', T, args=lambda: (large_at_end(),), driver='catch_exit'),
    Case('all_zero', T, args=lambda: (sp.csr_matrix((4, 4)),), kwargs={'M': 3}),
]
