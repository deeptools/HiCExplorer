from hicexplorer import HiCMatrix as hm
import os.path
import sys
from os import unlink
import numpy as np
import numpy.testing as nt
from scipy.sparse import csr_matrix
import warnings
from past.builtins import zip
from six import iteritems
import pytest
import logging

log = logging.getLogger(__name__)


warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"


def test_save_load():
    outfile = '/tmp/matrix.cool'
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]
    hic = hm.hiCMatrix(pCooler_only_init=True)
    hic.nan_bins = []
    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    # hic.matrix = csr_matrix(matrix)
    # make matrix symmetric
    hic.setMatrix(csr_matrix(matrix), cut_intervals)
    log.debug('hic.matrix {}'.format(hic.matrix))
    
    hic.fillLowerTriangle()
    log.debug('hic.matrix {}'.format(hic.matrix))
    log.debug('hic.matrix.shape {}'.format(hic.matrix.shape))


    hic.correction_factors = np.array([1, 1, 1, 1, 1])
    hic.nan_bins = np.array([])

    hic.save(outfile)

    cool_file = hm.hiCMatrix(outfile)
    # nt.assert_equal(hic.correction_factors, cool_file.correction_factors)
    log.debug('muh: {}'.format(hic.matrix.data))
    nt.assert_equal(hic.matrix.data, cool_file.matrix.data)
    nt.assert_equal(hic.matrix.indices, cool_file.matrix.indices)
    nt.assert_equal(hic.matrix.indptr, cool_file.matrix.indptr)
    nt.assert_equal(hic.nan_bins, cool_file.nan_bins)

    nt.assert_equal(hic.cut_intervals, cool_file.cut_intervals)
    unlink(outfile)


def test_convert_to_zscore_matrix():

    # make test matrix
    m_size = 100
    mat = np.triu(np.random.randint(0, 101, (m_size, m_size)))
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
    hic = hm.hiCMatrix(ROOT + '/Li_et_al_2015.cool')
    hic.maskBins(hic.nan_bins)

    mat = hic.matrix.todense()
    max_depth = 10000
    bin_size = hic.getBinSize()
    max_depth_in_bins = int(float(max_depth) / bin_size)

    m_size = mat.shape[0]
    # compute matrix values per distance
    chrom, start, end, extra = zip(*hm.hiCMatrix.fit_cut_intervals(hic.cut_intervals))
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


def test_save_load_cooler_format():
    outfile = '/tmp/matrix2.cool'
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]
    hic = hm.hiCMatrix()
    hic.nan_bins = []
    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    # make matrix symmetric
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.fillLowerTriangle()

    hic.save(outfile)

    matrix_cool = hm.hiCMatrix(outfile)

    log.info('original data: {}'.format(hic.matrix))
    log.info('cool data: {}'.format(matrix_cool.matrix))
    nt.assert_equal(hic.matrix.data, matrix_cool.matrix.data)
    nt.assert_equal(hic.matrix.indices, matrix_cool.matrix.indices)
    nt.assert_equal(hic.matrix.indptr, matrix_cool.matrix.indptr)

    # nan_bins and correction_factor are not supported by cool-format

    nt.assert_equal(hic.cut_intervals, matrix_cool.cut_intervals)
    unlink(outfile)


@pytest.mark.xfail
def test_load_mcooler_format_fail():
    matrix = hm.hiCMatrix(ROOT + 'matrix.mcool')  # noqa: F841


def test_load_mcooler_format_success():
    matrix = hm.hiCMatrix(ROOT + "matrix.mcool::/1")  # noqa: F841

def test_load_homer_format():
    pass

def test_load_hicpro_format():
    pass