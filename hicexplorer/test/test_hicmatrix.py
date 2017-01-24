from hicexplorer import HiCMatrix as hm
#import os.path
from os import unlink
import numpy as np
import numpy.testing as nt
from scipy.sparse import csr_matrix


# ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"

def test_save_load():
    outfile = '/tmp/matrix.h5'
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
    ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]
    hic = hm.hiCMatrix()
    hic.nan_bins = []
    matrix = np.array([[ 1,  8,  5, 3, 0],
                       [ 0,  4, 15, 5, 1],
                       [ 0,  0,  0, np.nan, 2],
                       [ 0,  0,  0, 0, 1],
                       [ 0,  0,  0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    # make matrix symmetric
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.matrix = hm.hiCMatrix.fillLowerTriangle(hic.matrix)
    hic.correction_factors = np.array([0.5, 1, 2, 3, 4])
    hic.nan_bins = np.array([4])

    hic.save(outfile)

    h5 = hm.hiCMatrix(outfile)

    nt.assert_equal(hic.correction_factors, h5.correction_factors)
    nt.assert_equal(hic.matrix.data, h5.matrix.data)
    nt.assert_equal(hic.matrix.indices, h5.matrix.indices)
    nt.assert_equal(hic.matrix.indptr, h5.matrix.indptr)
    nt.assert_equal(hic.nan_bins, h5.nan_bins)
    assert hic.cut_intervals == h5.cut_intervals
    unlink(outfile)


def test_convert_to_zscore_matrix():

    # make test matrix
    m_size = 100
    mat = np.triu(np.random.random_integers(0, 100, (m_size, m_size)))
    mu = dict([(idx, mat.diagonal(idx).mean()) for idx in range(mat.shape[0])])
    std = dict([(idx, np.std(mat.diagonal(idx))) for idx in range(mat.shape[0])])

    # compute z-score for test matrix
    zscore_mat = np.zeros((m_size, m_size))
    for _i in range(mat.shape[0]):
        for _j in range(mat.shape[0]):
            if _j >= _i:
                diag = _j- _i
                zscore = (mat[_i, _j] - mu[diag]) / std[diag]
                zscore_mat[_i, _j] = zscore

    # make Hi-C matrix based on test matrix
    hic = hm.hiCMatrix()
    hic.matrix = csr_matrix(mat)
    cut_intervals = [('chr', idx, idx+10, 0) for idx in range(0, mat.shape[0] * 10, 10)]
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.convert_to_zscore_matrix()

    from numpy.testing import assert_almost_equal
    assert_almost_equal(hic.matrix.todense(), zscore_mat)