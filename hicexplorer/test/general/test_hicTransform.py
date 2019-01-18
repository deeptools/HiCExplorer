from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicTransform
from hicmatrix import HiCMatrix as hm
import numpy.testing as nt

from tempfile import NamedTemporaryFile
import os


ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")
original_matrix = ROOT + "small_test_matrix_50kb_res.h5"
original_matrix_cool = ROOT + "small_test_matrix.cool"

DELTA_DECIMAL = 0


def test_hic_transfer_obs_exp():

    outfile = NamedTemporaryFile(suffix='obs_exp_.cool', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method obs_exp".format(original_matrix_cool, outfile.name).split()
    hicTransform.main(args)

    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp.cool")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    os.unlink(outfile.name)


def test_hic_transfer_obs_exp_perChromosome():

    outfile = NamedTemporaryFile(suffix='obs_exp_.cool', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method obs_exp --perChromosome".format(original_matrix_cool, outfile.name).split()
    hicTransform.main(args)

    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp_per_chromosome.cool")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    os.unlink(outfile.name)


def test_hic_transfer_obs_exp_non_zero():

    outfile = NamedTemporaryFile(suffix='obs_exp_.cool', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method obs_exp_non_zero".format(original_matrix_cool, outfile.name).split()
    hicTransform.main(args)

    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp_non_zero.cool")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    os.unlink(outfile.name)


def test_hic_transfer_obs_exp_non_zero_perChromosome():

    outfile = NamedTemporaryFile(suffix='obs_exp_.cool', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method obs_exp_non_zero --perChromosome".format(original_matrix_cool, outfile.name).split()
    hicTransform.main(args)

    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp_non_zero_per_chromosome.cool")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    os.unlink(outfile.name)


def test_hic_transfer_obs_exp_lieberman():
    outfile = NamedTemporaryFile(suffix='obs_exp_lieberman_.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method obs_exp_lieberman".format(original_matrix, outfile.name).split()
    hicTransform.main(args)

    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp_lieberman.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    os.unlink(outfile.name)


def test_hic_transfer_obs_exp_norm():
    outfile = NamedTemporaryFile(suffix='obs_exp_norm_.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method obs_exp_norm".format(original_matrix, outfile.name).split()
    hicTransform.main(args)

    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp_norm.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    os.unlink(outfile.name)


def test_hic_transfer_obs_exp_norm_perChromosome():
    outfile = NamedTemporaryFile(suffix='obs_exp_norm_.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method obs_exp_norm --perChromosome".format(original_matrix, outfile.name).split()
    hicTransform.main(args)

    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp_norm_perChromosome.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    os.unlink(outfile.name)


def test_hic_transfer_pearson():
    outfile = NamedTemporaryFile(suffix='pearson_.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method pearson".format(original_matrix, outfile.name).split()
    hicTransform.main(args)

    test = hm.hiCMatrix(ROOT + "hicTransform/pearson.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    os.unlink(outfile.name)


def test_hic_transfer_pearson_perChromosome():
    outfile = NamedTemporaryFile(suffix='pearson_.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method pearson --perChromosome".format(original_matrix, outfile.name).split()
    hicTransform.main(args)

    test = hm.hiCMatrix(ROOT + "hicTransform/pearson_perChromosome.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    os.unlink(outfile.name)


def test_hic_transfer_covariance():
    outfile = NamedTemporaryFile(suffix='covariance_.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method covariance".format(original_matrix, outfile.name).split()
    hicTransform.main(args)
    test = hm.hiCMatrix(ROOT + "hicTransform/covariance.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    os.unlink(outfile.name)


def test_hic_transfer_covariance_perChromosome():
    outfile = NamedTemporaryFile(suffix='covariance_.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method covariance --perChromosome".format(original_matrix, outfile.name).split()
    hicTransform.main(args)
    test = hm.hiCMatrix(ROOT + "hicTransform/covariance_perChromosome.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    os.unlink(outfile.name)
