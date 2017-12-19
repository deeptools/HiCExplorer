from __future__ import division
from hicexplorer import hicTransform
from hicexplorer import HiCMatrix as hm
import numpy.testing as nt

from tempfile import NamedTemporaryFile
import os
from os.path import basename, dirname


ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
original_matrix = ROOT + "small_test_matrix_50kb_res.h5"


def test_hic_transfer_all():
    outfile = NamedTemporaryFile(suffix='all.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method all".format(original_matrix, outfile.name).split()
    hicTransform.main(args)

    dirname_new = dirname(outfile.name)
    basename_new = basename(outfile.name)
    # obs_exp
    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp_small_50kb.h5")
    new = hm.hiCMatrix(dirname_new + "/obs_exp_" + basename_new)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data)
    os.unlink(dirname_new + "/obs_exp_" + basename_new)

    # pearson
    test = hm.hiCMatrix(ROOT + "hicTransform/pearson_small_50kb.h5")
    new = hm.hiCMatrix(dirname_new + "/pearson_" + basename_new)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data)
    os.unlink(dirname_new + "/pearson_" + basename_new)

    # covariance
    test = hm.hiCMatrix(ROOT + "hicTransform/covariance_small_50kb.h5")
    new = hm.hiCMatrix(dirname_new + "/covariance_" + basename_new)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data)
    os.unlink(dirname_new + "/covariance_" + basename_new)
    os.unlink(outfile.name)


def test_hic_transfer_obs_exp():
    outfile = NamedTemporaryFile(suffix='obs_exp_.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method obs_exp".format(original_matrix, outfile.name).split()
    hicTransform.main(args)

    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp_small_50kb.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data)
    os.unlink(outfile.name)


def test_hic_transfer_pearson():
    outfile = NamedTemporaryFile(suffix='pearson_.h5', delete=False)
    outfile.close()

    matrix = ROOT + "/hicTransform/obs_exp_small_50kb.h5"
    args = "--matrix {} --outFileName {} --method pearson".format(matrix, outfile.name).split()
    hicTransform.main(args)

    test = hm.hiCMatrix(ROOT + "hicTransform/pearson_small_50kb.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data)
    os.unlink(outfile.name)


def test_hic_transfer_covariance():
    outfile = NamedTemporaryFile(suffix='covariance_.h5', delete=False)
    outfile.close()
    matrix = ROOT + "/hicTransform/pearson_small_50kb.h5"

    args = "--matrix {} --outFileName {} --method covariance".format(matrix, outfile.name).split()
    hicTransform.main(args)
    test = hm.hiCMatrix(ROOT + "hicTransform/covariance_small_50kb.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data)
    os.unlink(outfile.name)
