from hicexplorer import HiCMatrix as hm
import os.path
import sys
from os import unlink
import numpy as np
import numpy.testing as nt
from scipy.sparse import csr_matrix
from past.builtins import zip
from six import iteritems


ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"


def test_save_load():
    outfile = '/tmp/matrix.h5'
    cut_intervals = [(b'a', 0, 10, 1), (b'a', 10, 20, 1),
                     (b'a', 20, 30, 1), (b'a', 30, 40, 1), (b'b', 40, 50, 1)]
    hic = hm.hiCMatrix()
    hic.nan_bins = []
    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, np.nan, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

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

    nt.assert_equal(hic.cut_intervals, h5.cut_intervals)
    unlink(outfile)


def test_convert_to_zscore_matrix():

    # make test matrix
    m_size = 100
    mat = np.triu(np.random.random_integers(0, 100, (m_size, m_size)))
    # add a number of zeros
    mat[mat < 90] = 0
    # import ipdb;ipdb.set_trace()
    mu = dict([(idx, mat.diagonal(idx).mean()) for idx in range(mat.shape[0])])
    std = dict([(idx, np.std(mat.diagonal(idx))) for idx in range(mat.shape[0])])

    # compute z-score for test matrix
    zscore_mat = np.zeros((m_size, m_size))
    for _i in range(mat.shape[0]):
        for _j in range(mat.shape[0]):
            if _j >= _i:
                diag = _j - _i
                if std[diag] == 0:
                    zscore = np.nan
                else:
                    zscore = (mat[_i, _j] - mu[diag]) / std[diag]
                zscore_mat[_i, _j] = zscore

    # make Hi-C matrix based on test matrix
    hic = hm.hiCMatrix()
    hic.matrix = csr_matrix(mat)
    cut_intervals = [('chr', idx, idx + 10, 0) for idx in range(0, mat.shape[0] * 10, 10)]
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.convert_to_zscore_matrix()

    from numpy.testing import assert_almost_equal
    assert_almost_equal(hic.matrix.todense(), zscore_mat)


def test_convert_to_zscore_matrix_2():

    # load test matrix
    hic = hm.hiCMatrix(ROOT + '/Li_et_al_2015.h5')
    hic.maskBins(hic.nan_bins)

    mat = hic.matrix.todense()
    max_depth = 10000
    bin_size = hic.getBinSize()
    max_depth_in_bins = int(float(max_depth) / bin_size)

    m_size = mat.shape[0]
    # compute matrix values per distance
    chrom, start, end, extra = zip(*hm.hiCMatrix.fit_cut_intervals(hic.cut_intervals))
#    chrom, start, end, extra = zip(*hic.cut_intervals)
    dist_values = {}
    sys.stderr.write("Computing values per distance for each matrix entry\n")

    for _i in range(mat.shape[0]):
        for _j in range(mat.shape[0]):
            if _j >= _i:
                # dist is translated to bins
                dist = int(float(start[_j] - start[_i]) / bin_size)
                if dist <= max_depth_in_bins:
                    if dist not in dist_values:
                        dist_values[dist] = []
                    dist_values[dist].append(mat[_i, _j])

    mu = {}
    std = {}
    for dist, values in iteritems(dist_values):
        mu[dist] = np.mean(values)
        std[dist] = np.std(values)

    # compute z-score for test matrix
    sys.stderr.write("Computing zscore for each matrix entry\n")
    zscore_mat = np.full((m_size, m_size), np.nan)
    for _i in range(mat.shape[0]):
        for _j in range(mat.shape[0]):
            if _j >= _i:
                dist = int(float(start[_j] - start[_i]) / bin_size)
                if dist <= max_depth_in_bins:
                    zscore = (mat[_i, _j] - mu[dist]) / std[dist]
                    zscore_mat[_i, _j] = zscore

    # compare with zscore from class
    hic.convert_to_zscore_matrix(maxdepth=max_depth)

    from numpy.testing import assert_almost_equal
    # only the main diagonal is check. Other diagonals show minimal differences
    assert_almost_equal(hic.matrix.todense().diagonal(0).A1, zscore_mat.diagonal(0))
