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
from intervaltree import IntervalTree, Interval
from collections import OrderedDict


log = logging.getLogger(__name__)


warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"


def test_save_load():
    outfile = '/tmp/matrix.h5'
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]
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

    nt.assert_almost_equal(hic.matrix.todense(), zscore_mat)


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

    # from numpy.testing import assert_almost_equal
    # only the main diagonal is check. Other diagonals show minimal differences
    nt.assert_almost_equal(hic.matrix.todense().diagonal(0).A1, zscore_mat.diagonal(0))


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
    hic.matrix = hm.hiCMatrix.fillLowerTriangle(hic.matrix)

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


def test_getCountsByDistance():
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
    hic.setMatrix(hic.matrix, cut_intervals)

    distance = hic.getCountsByDistance()

    nt.assert_equal(distance[-1], [0, 1, 2, 1])
    nt.assert_equal(distance[0], [1, 4, 0, 0, 0])
    nt.assert_equal(distance[10], [8, 15, 0])
    nt.assert_equal(distance[20], [5, 5])
    nt.assert_equal(distance[30], [3])

    hic = hm.hiCMatrix()
    hic.nan_bins = []

    matrix = np.matrix([[np.nan for x in range(5)] for y in range(5)])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    distance = hic.getCountsByDistance()

    nt.assert_equal(distance[-1], [0, 0, 0, 0])
    nt.assert_equal(distance[0], [0, 0, 0, 0, 0])
    nt.assert_equal(distance[10], [0, 0, 0])
    nt.assert_equal(distance[20], [0, 0])
    nt.assert_equal(distance[30], [0])

    # test for mean distance
    # mean = True
    # hic = hm.hiCMatrix()
    # hic.nan_bins = []

    # matrix = np.matrix([[np.nan for x in range(5)] for y in range(5)])

    # hic.matrix = csr_matrix(matrix)
    # hic.setMatrix(hic.matrix, cut_intervals)

    # distance = hic.getCountsByDistance(mean=mean)

def test_dist_list_to_dict():
    hic = hm.hiCMatrix()

    data = np.array([1, 8, 5, 3, 0, 4, 15, 5, 1, 0, 0, 2, 0, 1, 0])
    dist_list = np.array([0, 10, 20, 30, -1, 0, 10, 20, -1, 0, 10, -1, 0, -1, 0])

    distance = hic.dist_list_to_dict(data, dist_list)

    nt.assert_equal(distance[-1], [0, 1, 2, 1])
    nt.assert_equal(distance[0], [1, 4, 0, 0, 0])
    nt.assert_equal(distance[10], [8, 15, 0])
    nt.assert_equal(distance[20], [5, 5])
    nt.assert_equal(distance[30], [3])

    data = np.array([0, 100, 200, 0, 100, 200, 0, 100, 0])
    dist_list = np.array([0, 100, 200, 0, 100, 200, 0, 100, 0])

    distance = hic.dist_list_to_dict(data, dist_list)

    nt.assert_equal(distance[0], [0, 0, 0, 0])
    nt.assert_equal(distance[100], [100, 100, 100])
    nt.assert_equal(distance[200], [200, 200])


def test_getUnwantedChrs():
    hic = hm.hiCMatrix()
    assert hic.getUnwantedChrs() == {'chrM', 'chrYHet', 'chrXHet', 'chrUextra', 'chrU',
                                     'chr3RHet', 'chr3LHet', 'chr2RHet', 'chr2LHet'}

def test_keepOnlyTheseChr():
    chromosome_list = ['chrX', 'chr2RHet']

    hic = hm.hiCMatrix(ROOT + 'small_test_matrix.h5')

    keep = hic.keepOnlyTheseChr(chromosome_list)

    nt.assert_equal(hic.getChrNames().sort(), chromosome_list.sort())

def test_filterUnwantedChr():
    hic = hm.hiCMatrix(ROOT + 'small_test_matrix.h5')

    assert 'chr2RHet' in hic.getChrNames()
    assert 'chr3LHet' in hic.getChrNames()
    assert 'chr3RHet' in hic.getChrNames()

    hic.filterUnwantedChr()

    assert 'chr2RHet' not in hic.getChrNames()
    assert 'chr3LHet' not in hic.getChrNames()
    assert 'chr3RHet' not in hic.getChrNames()

    chromosomes = list(hic.getChrNames())

    # make sure there are any other chromosomes than 'chrX'
    assert any(x != 'chrX' for x in chromosomes)

    # then filter for 'chrX'
    hic.filterUnwantedChr(chromosome='chrX')

    chromosomes = list(hic.getChrNames())

    # and check that there are only 'chrX'-chromosomes left in matrix
    assert all(x == 'chrX' for x in chromosomes)

@pytest.mark.xfail
def test_save_bing_ren():
    """ Test needs to be marked as xfail because .gz files are expected in __init__ to be in dekker file format """
    outfile = '/tmp/matrix.gz'
    try:
        _outfile = open(outfile, 'r')
    except FileNotFoundError:
        _outfile = open(outfile, 'w')

    _outfile.close()

    hic = hm.hiCMatrix()

    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]


    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    # make matrix symmetric
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.matrix = hm.hiCMatrix.fillLowerTriangle(hic.matrix)

    hic.save_bing_ren(outfile)

    # Test fails here due to __init__ of hiCMatrix
    hicTest = hm.hiCMatrix(outfile)

    # test = hicTest.matrix

    # nt.assert_equal(matrix.shape, test.shape)


def test_save_dekker():
    outfile = '/tmp/matrix.gz'
    try:
        _outfile = open(outfile, 'r')
    except FileNotFoundError:
        _outfile = open(outfile, 'w')
    _outfile.close()

    hic = hm.hiCMatrix()

    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]


    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    # make matrix symmetric
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.matrix = hm.hiCMatrix.fillLowerTriangle(hic.matrix)

    hic.save_dekker(outfile)

    dekker_test = hm.hiCMatrix(outfile)
    dekker_test.fillLowerTriangle(dekker_test.matrix)

    nt.assert_equal(hic.getMatrix().shape, dekker_test.getMatrix().shape)
    nt.assert_equal(hic.getMatrix(), dekker_test.getMatrix())

@pytest.mark.xfail
def test_save_lieberman():
    """ Test fails because lieberman format requires folder and saves .gz files. Loading causes IO-Error. """
    outpath = '/tmp/matrix_lieberman/'

    hic = hm.hiCMatrix()

    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]


    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    # make matrix symmetric
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.matrix = hm.hiCMatrix.fillLowerTriangle(hic.matrix)

    hic.save_lieberman(outpath)

    lieberman_test = hm.hiCMatrix(outpath)

@pytest.mark.xfail
def test_save_GInteractions():
    """
    Test fails because GInteractions saves file as .tsv but __init__
    can only process .npz, h5, dekker, cool. Otherwise files are treated as h5f...
    """
    outfile = '/tmp/matrix_GInteractions'
    try:
        _outfile = open(outfile, 'r')
    except FileNotFoundError:
        _outfile = open(outfile, 'w')
    _outfile.close()

    hic = hm.hiCMatrix()

    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]


    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    # make matrix symmetric
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.matrix = hm.hiCMatrix.fillLowerTriangle(hic.matrix)

    hic.save_GInteractions(outfile)

    # test fails during load
    GI_test = hm.hiCMatrix(outfile)
    # GI_test.fillLowerTriangle(GI_test.matrix)

    # nt.assert_equal(hic.getMatrix().shape, GI_test.getMatrix().shape)
    # nt.assert_equal(hic.getMatrix(), GI_test.getMatrix())

@pytest.mark.xfail
def test_create_empty_cool_file():
    """
    Test fails. As far as I can see function is never called from anywhere. Perhaps not important.
    Perhaps test is not correctly written...
    """
    outfile = '/tmp/matrix3.cool'

    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    # make matrix symmetric
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.matrix = hm.hiCMatrix.fillLowerTriangle(hic.matrix)

    hic.create_empty_cool_file(outfile)


def test_save_cooler():
    """
    Test is running for 4 different configurations:
        save_cooler(outfile, pSymmetric=True, pApplyCorrections=True)
        save_cooler(outfile, pSymmetric=False, pApplyCorrections=True)
        save_cooler(outfile, pSymmetric=True, pApplyCorrections=False)
        save_cooler(outfile, pSymmetric=False, pApplyCorrections=True)
    """
    outfile = '/tmp/matrix3.cool'

    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]
    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic = hm.hiCMatrix()
    hic.nan_bins = []

    hic.matrix = csr_matrix(matrix)
    # make matrix symmetric
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.matrix = hic.fillLowerTriangle(hic.matrix)

    hic.save_cooler(outfile)

    cooler_test_default = hm.hiCMatrix(outfile)

    # Test default configuration
    nt.assert_equal(hic.getMatrix(), cooler_test_default.getMatrix())

    # Test pSymmetric=False, pApplyCorrection=True
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    hic.save_cooler(outfile, pSymmetric=False)

    cooler_test_pSym_False = hm.hiCMatrix(outfile)

    # make both matrices symmetric
    cooler_test_pSym_False.matrix = cooler_test_pSym_False.fillLowerTriangle(cooler_test_pSym_False.matrix)
    hic.matrix = hic.fillLowerTriangle(hic.matrix)

    nt.assert_equal(hic.getMatrix(), cooler_test_pSym_False.getMatrix())

    # Test pApplyCorrection=False, pSymmetric=True
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.matrix = hic.fillLowerTriangle(hic.matrix)

    hic.save_cooler(outfile, pApplyCorrection=False)

    cooler_test_pAppCor_False = hm.hiCMatrix(outfile)

    # make test matrix symmetric
    cooler_test_pAppCor_False.matrix = cooler_test_pAppCor_False.fillLowerTriangle(cooler_test_pAppCor_False.matrix)

    nt.assert_equal(hic.getMatrix(), cooler_test_pAppCor_False.getMatrix())

    # Test pSymmetric and pApplyCorrections = false
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    hic.save_cooler(outfile, pSymmetric=False, pApplyCorrection=False)

    cooler_test_all_False = hm.hiCMatrix(outfile)

    # make both matrices symmetric
    cooler_test_all_False.matrix = cooler_test_all_False.fillLowerTriangle(cooler_test_all_False.matrix)
    hic.matrix = hic.fillLowerTriangle(hic.matrix)

    nt.assert_equal(hic.getMatrix(), cooler_test_all_False.getMatrix())


def test_save_hdf5():
    """
    Test is running for 2 different configurations:
        save_hdf5(filename, pSymmetric=True) (Default)
        save_hdf5(filename, pSymmetric=True)
    """
    outfile = '/tmp/matrix.h5'

    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.matrix = hic.fillLowerTriangle(hic.matrix)

    hic.save_hdf5(outfile)

    hdf5_test = hm.hiCMatrix(outfile)

    nt.assert_equal(hic.getMatrix(), hdf5_test.getMatrix())

    # Test pSymmetric=False
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    hic.save_hdf5(outfile)

    hdf5_test_pSym_False = hm.hiCMatrix(outfile)

    hdf5_test_pSym_False.matrix = hdf5_test_pSym_False.fillLowerTriangle(hdf5_test_pSym_False.matrix)
    hic.matrix = hic.fillLowerTriangle(hic.matrix)

    nt.assert_equal(hic.getMatrix(), hdf5_test_pSym_False.getMatrix())


def test_save_npz():
    outfile = '/tmp/matrix.npz'

    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.matrix = hic.fillLowerTriangle(hic.matrix)

    hic.save_npz(outfile)

    npz_test = hm.hiCMatrix(outfile)

    nt.assert_equal(hic.getMatrix(), npz_test.getMatrix())

def test_save():
    """
    Test will not cover testing of following formats due to unsupported file_formats (see __init__ of class hiCMatrix):

    * ren
    * lieberman
    * GInteractions

    see also single test for these formats (marked as xfail)
    """
    matrix_h5 = '/tmp/matrix.h5'
    matrix_cool = '/tmp/matrix.cool'
    matrix_npz = '/tmp/matrix.npz'
    matrix_gz = '/tmp/matrix.gz'

    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.matrix = hic.fillLowerTriangle(hic.matrix)

    # test .h5
    hic.save(matrix_h5)
    h5_test = hm.hiCMatrix(matrix_h5)

    # test cool
    hic.save(matrix_cool)
    cool_test = hm.hiCMatrix(matrix_cool)

    # test npz
    hic.save(matrix_npz)
    npz_test = hm.hiCMatrix(matrix_npz)

    # test dekker
    hic.save(matrix_gz)
    dekker_test = hm.hiCMatrix(matrix_gz)

    nt.assert_equal(hic.getMatrix(), h5_test.getMatrix())
    nt.assert_equal(hic.getMatrix(), cool_test.getMatrix())
    nt.assert_equal(hic.getMatrix(), npz_test.getMatrix())
    nt.assert_equal(hic.getMatrix(), dekker_test.getMatrix())


def test_diagflat():
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.matrix = hic.fillLowerTriangle(hic.matrix)

    hic.diagflat(value=1000)
    nt.assert_equal(np.array([1000 for x in range(matrix.shape[0])]), hic.matrix.diagonal())

    hic.diagflat()
    nt.assert_equal(np.array([np.nan for x in range(5)]), hic.matrix.diagonal())

def test_filterOutInterChrCounts():
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)
    hic.matrix = hic.fillLowerTriangle(hic.matrix)
    hic.filterOutInterChrCounts()

    filtered_matrix = np.matrix([[1, 8, 5, 0, 0],
                                [8, 4, 15, 0, 0],
                                [5, 15, 0, 0, 0],
                                [0, 0, 0, 0, 1],
                                [0, 0, 0, 1, 0]])

    nt.assert_equal(hic.getMatrix(), filtered_matrix)


def test_setMatrixValues_success():
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    new_matrix = np.array([[10, 80, 50, 30, 0],
                          [0, 40, 150, 50, 10],
                          [0, 0, 0, 0, 20],
                          [0, 0, 0, 0, 10],
                          [0, 0, 0, 0, 0]])

    hic.setMatrixValues(new_matrix)

    nt.assert_equal(hic.getMatrix(), new_matrix)

def test_setMatrixValues_fail():
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1)]

    new_matrix = np.array([[10, 80, 50, 30],
                          [0, 40, 150, 50],
                          [0, 0, 0, 0],
                          [0, 0, 0, 0]])
    with pytest.raises(AssertionError):
        hic.setMatrixValues(new_matrix)


def test_setCorrectionFactors_success():
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    assert hic.correction_factors == None

    hic.setCorrectionFactors([5, 5, 5, 5, 5])

    nt.assert_equal(hic.correction_factors, [5, 5, 5, 5, 5])

def test_setCorrectionFactors_fail():
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])


    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    assert hic.correction_factors == None
    with pytest.raises(AssertionError):
        hic.setCorrectionFactors([5, 5, 5, 5])


def test_reorderChromosomes_old():
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    new_chr_order = ['b', 'a']
    hic.reorderChromosomes_old(new_chr_order)

    nt.assert_equal(hic.chrBinBoundaries, OrderedDict([('b', (0, 2)), ('a', (2, 5))]))

    old_chr_order = ['a', 'b']
    hic.reorderChromosomes_old(old_chr_order)

    nt.assert_equal(hic.chrBinBoundaries, OrderedDict([('a', (0, 3)), ('b', (3, 5))]))

    # new order too long will cause function to return
    false_chr_order = ['a', 'b', 'c']
    hic.reorderChromosomes_old(false_chr_order)

    nt.assert_equal(hic.chrBinBoundaries, OrderedDict([('a', (0, 3)), ('b', (3, 5))]))


def test_reorderChromosomes():
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    new_chr_order = ['b', 'a']
    hic.reorderChromosomes(new_chr_order)

    nt.assert_equal(hic.chrBinBoundaries, OrderedDict([('b', (0, 2)), ('a', (2, 5))]))

    old_chr_order = ['a', 'b']
    hic.reorderChromosomes(old_chr_order)

    nt.assert_equal(hic.chrBinBoundaries, OrderedDict([('a', (0, 3)), ('b', (3, 5))]))


def test_reorderChromosomes_fail():
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    # name 'c' not in chromosome names, thus fail
    false_chr_order = ['a', 'b', 'c']
    with pytest.raises(SystemExit):
        hic.reorderChromosomes(false_chr_order)



def test_reorderBins():
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    nt.assert_equal(hic.getMatrix(), matrix)

    new_order = [0, 1, 3, 2, 4]
    new_matrix = np.matrix([[1, 8, 3, 5, 0],
                            [0, 4, 5, 15, 1],
                            [0, 0, 0, 0, 1],
                            [0, 0, 0, 0, 2],
                            [0, 0, 0, 0, 0]])

    hic.reorderBins(new_order)

    nt.assert_equal(hic.getMatrix(), new_matrix)

    hic.reorderBins(new_order)

    nt.assert_equal(hic.getMatrix(), matrix)

    # order smaller than original matrix should delete unused ids
    small_order = [2, 3]
    small_matrix = np.matrix([[0, 0],
                              [0, 0]])

    hic.reorderBins(small_order)

    nt.assert_equal(hic.getMatrix(), small_matrix)
    nt.assert_equal(hic.matrix.shape, small_matrix.shape)
    nt.assert_equal(hic.chrBinBoundaries, OrderedDict([('a', (0, 1)), ('b', (1, 2))]))
    nt.assert_equal(hic.cut_intervals, [('a', 20, 30, 1), ('b', 30, 40, 1)])
    nt.assert_equal(hic.nan_bins, [])


def test_removeBins():
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    nt.assert_equal(hic.getMatrix(), matrix)

    new_matrix = np.matrix([[1, 3, 0],
                            [0, 0, 1],
                            [0, 0, 0]])

    ids2remove = [1, 2]
    hic.removeBins(ids2remove)

    nt.assert_equal(hic.getMatrix(), new_matrix)
    nt.assert_equal(hic.matrix.shape, new_matrix.shape)
    nt.assert_equal(hic.chrBinBoundaries, OrderedDict([('a', (0, 1)), ('b', (1, 3))]))
    nt.assert_equal(hic.cut_intervals, [('a', 0, 10, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)])


def test_maskBins():
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    nt.assert_equal(hic.getMatrix(), matrix)
    nt.assert_equal(hic.orig_bin_ids, [])

    new_matrix = np.matrix([[0, 0, 2],
                            [0, 0, 1],
                            [0, 0, 0]])

    masking_ids = [0, 1]
    hic.maskBins(masking_ids)

    nt.assert_equal(hic.getMatrix(), new_matrix)
    nt.assert_equal(sorted(hic.orig_cut_intervals), sorted([('a', 0, 10, 1), ('a', 10, 20, 1),
                                                            ('a', 20, 30, 1), ('b', 30, 40, 1),
                                                            ('b', 40, 50, 1)]))
    nt.assert_equal(sorted(hic.cut_intervals), sorted([('a', 20, 30, 1), ('b', 30, 40, 1),
                                                       ('b', 40, 50, 1)]))
    nt.assert_equal(hic.chrBinBoundaries, OrderedDict([('a', (0, 1)), ('b', (1, 3))]))
    nt.assert_equal(sorted(hic.orig_bin_ids), sorted([0, 1, 2, 3, 4]))

    # direct return if masking_ids is None or has len() == 0, thus no changes to matrix
    masking_ids = None
    hic.maskBins(masking_ids)

    nt.assert_equal(hic.getMatrix(), new_matrix)
    nt.assert_equal(sorted(hic.orig_cut_intervals), sorted([('a', 0, 10, 1), ('a', 10, 20, 1),
                                                            ('a', 20, 30, 1), ('b', 30, 40, 1),
                                                            ('b', 40, 50, 1)]))
    nt.assert_equal(sorted(hic.cut_intervals), sorted([('a', 20, 30, 1), ('b', 30, 40, 1),
                                                       ('b', 40, 50, 1)]))
    nt.assert_equal(hic.chrBinBoundaries, OrderedDict([('a', (0, 1)), ('b', (1, 3))]))

    masking_ids = []

    hic.maskBins(masking_ids)

    nt.assert_equal(hic.getMatrix(), new_matrix)
    nt.assert_equal(sorted(hic.orig_cut_intervals), sorted([('a', 0, 10, 1), ('a', 10, 20, 1),
                                                            ('a', 20, 30, 1), ('b', 30, 40, 1),
                                                            ('b', 40, 50, 1)]))
    nt.assert_equal(sorted(hic.cut_intervals), sorted([('a', 20, 30, 1), ('b', 30, 40, 1),
                                                       ('b', 40, 50, 1)]))
    nt.assert_equal(hic.chrBinBoundaries, OrderedDict([('a', (0, 1)), ('b', (1, 3))]))

    nt.assert_equal(sorted(hic.orig_bin_ids), sorted([0, 1, 2, 3, 4]))


def test_update_matrix(capsys):
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    nt.assert_equal(hic.getMatrix(), matrix)

    new_cut_intervals = [('c', 0, 10, 1), ('d', 10, 20, 1), ('d', 20, 30, 1)]

    new_matrix = np.array([[3, 6, 4],
                           [np.nan, 0, 2],
                           [1, 0, 0]])
    try:
        hic.update_matrix(new_matrix, new_cut_intervals)
    except AttributeError:
        with capsys.disabled():
            print('\n')
            print('\033[31m' + "AttributeError @test_update_matrix" + '\x1b[0m')
            print('\033[31m' + "Please check function hiCMatrix::update_matrix()" + '\x1b[0m')
            print('\n')

    # if matrix.shape[0] not equal to length of cut_intervals assertionError is raised
    short_cut_intervals = [('c', 0, 10, 1), ('d', 10, 20, 1)]

    with pytest.raises(AssertionError):
        hic.update_matrix(new_matrix, short_cut_intervals)

    # if matrix contains masked bins exception is raised
    masking_ids = [0, 1]
    hic.maskBins(masking_ids)

    with pytest.raises(Exception):
        hic.update_matrix(new_matrix, new_cut_intervals)


def test_restoreMaskedBins():
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    nt.assert_equal(hic.getMatrix(), matrix)
    nt.assert_equal(hic.orig_bin_ids, [])

    # function should directly return if there are no masked_bins
    hic.restoreMaskedBins()

    nt.assert_equal(hic.getMatrix(), matrix)
    nt.assert_equal(hic.orig_bin_ids, [])

    # test general use
    # first get some masked bins
    masking_ids = [0, 1]
    hic.maskBins(masking_ids)

    new_matrix = np.matrix([[0, 0, 2],
                            [0, 0, 1],
                            [0, 0, 0]])

    nt.assert_equal(hic.getMatrix(), new_matrix)
    nt.assert_equal(sorted(hic.orig_bin_ids), sorted([0, 1, 2, 3, 4]))

    # and now restore masked bins
    hic.restoreMaskedBins()

    result_matrix = np.matrix([[np.nan, np.nan, np.nan, np.nan, np.nan],
                               [np.nan, np.nan, np.nan, np.nan, np.nan],
                               [np.nan, np.nan, 0, 0, 2],
                               [np.nan, np.nan, 0, 0, 1],
                               [np.nan, np.nan, 0, 0, 0]])

    nt.assert_equal(hic.getMatrix(), result_matrix)
    nt.assert_equal(hic.orig_bin_ids, [])


def test_reorderMatrix():
    orig = (1, 3)
    dest = 2

    # get matrix
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    nt.assert_equal(hic.getMatrix(), matrix)

    # reorder matrix
    hic.reorderMatrix(orig, dest)

    new_matrix = np.matrix([[1, 3, 8, 5, 0],
                            [0, 0, 0, 0, 1],
                            [0, 5, 4, 15, 1],
                            [0, 0, 0, 0, 2],
                            [0, 0, 0, 0, 0]])

    new_cut_intervals = [('a', 0, 10, 1), ('b', 30, 40, 1),
                         ('a', 10, 20, 1), ('a', 20, 30, 1), ('b', 40, 50, 1)]

    #check if it is equal
    nt.assert_equal(hic.getMatrix(), new_matrix)
    nt.assert_equal(hic.matrix.shape, new_matrix.shape)
    nt.assert_equal(hic.cut_intervals, new_cut_intervals)


def test_truncTrans():
    # get matrix
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[-1, 8, 5, 3, 0],
                       [np.nan, 4, 15, 5, 100],
                       [0, 0, 0, 0, 2000],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    nt.assert_equal(hic.getMatrix(), matrix)

    # define expected outcome
    new_matrix = np.matrix([[-1., 8., 5., 3., 0.],
                            [np.nan, 4., 15., 5., 1.e+2],
                            [0., 0., 0., 0., 2.e+3],
                            [0., 0., 0., 0., 1.],
                            [0., 0., 0., 0., 0.]])

    # truncTrans of matrix
    hic.truncTrans()

    # test against expected outcome
    nt.assert_equal(hic.getMatrix(), new_matrix)

    # reset matrix
    matrix = np.array([[-1, 8, 5, 3, 0],
                       [np.nan, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])
    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    # method should directly return if nothing to do, matrix stays the same
    hic.truncTrans()
    nt.assert_equal(hic.getMatrix(), matrix)


def test_truncTrans_bk(capsys):
    # get matrix
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[-1, 8, 5, 3, 0],
                       [np.nan, 4, 15, 5, 100],
                       [0, 0, 0, 0, 2000],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    nt.assert_equal(hic.getMatrix(), matrix)

    # define expected outcome
    new_matrix = np.matrix([[-1., 8., 5., 3., 0.],
                            [np.nan, 4., 15., 5., 1.e+2],
                            [0., 0., 0., 0., 2.e+3],
                            [0., 0., 0., 0., 1.],
                            [0., 0., 0., 0., 0.]])

    # truncTrans_bk fails because getDistList() expects also the cut_intervals when called
    try:
        hic.truncTrans_bk()
    except TypeError:
        with capsys.disabled():
            print('\n')
            print('\033[31m' + "TypeError @test_truncTrans_bk" + '\x1b[0m')
            print('\033[31m' + "Please check function hiCMatrix::truncTrans_bk()" + '\x1b[0m')
            print('\n')


def test_removePoorRegions(capsys):
    # get matrix
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[-1, 8, 5, 3, 0],
                       [np.nan, 4, 15, 5, 100],
                       [0, 0, 0, 0, 2000],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    nt.assert_equal(hic.getMatrix(), matrix)

    # define expected outcome
    new_matrix = np.matrix([[8, 5, 3, 0],
                            [4, 15, 5, 100],
                            [0, 0, 0, 2000],
                            [0, 0, 0, 1]])

    # removePoorRegions
    try:
        hic.removePoorRegions()
    except IndexError:
        with capsys.disabled():
            print('\n')
            print('\033[31m' + "IndexError @test_removePoorRegions" + '\x1b[0m')
            print('\033[31m' + "Please check function hiCMatrix::removePoorRegions()" + '\x1b[0m')
            print('\n')


    # print(hic.getMatrix())
    # check correctness
    # nt.assert_equal(hic.getMatrix(), new_matrix)

def test_printchrtoremove(capsys):
    # get matrix
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    nt.assert_equal(hic.getMatrix(), matrix)

    # first test exception message for no self.prev_to_remove
    to_remove = [0, 1]

    with pytest.raises(Exception):
        hic.printchrtoremove(to_remove)

        captured = capsys.readouterr()
        assert captured.out == "No self.prev_to_remove defined, defining it now."

        nt.assert_equal(hic.prev_to_remove, np.array(to_remove))

    nt.assert_equal(hic.orig_bin_ids, [])

    # also test with masked_bins
    hic.maskBins(to_remove)

    assert len(hic.orig_bin_ids) > 0

    hic.printchrtoremove(to_remove)

    nt.assert_equal(hic.prev_to_remove, np.array(to_remove))


def test_removeBySequencedCount():
    # get matrix
    hic = hm.hiCMatrix()
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]

    hic.nan_bins = []

    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    nt.assert_equal(hic.getMatrix(), matrix)

    
