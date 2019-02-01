import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicmatrix import HiCMatrix as hm
from hicexplorer import hicNormalize

from tempfile import NamedTemporaryFile

import os
import numpy.testing as nt


ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/hicNormalize")

matrix_one_h5 = ROOT + '/small_test_matrix.h5'
matrix_two_h5 = ROOT + '/small_test_matrix_scaled_up.h5'

matrix_one_cool = ROOT + '/small_test_matrix.cool'
matrix_two_cool = ROOT + '/small_test_matrix_scaled_up.cool'


def test_normalize_smallest(capsys):
    outfile_one = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_one.close()

    outfile_two = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_two.close()

    args = "--matrices {} {} --normalize smallest -o {} {}".format(matrix_one_h5, matrix_two_h5,
                                                                   outfile_one.name, outfile_two.name).split()
    hicNormalize.main(args)

    test_one = hm.hiCMatrix(ROOT + "/smallest_one.h5")
    test_two = hm.hiCMatrix(ROOT + "/smallest_two.h5")

    new_one = hm.hiCMatrix(outfile_one.name)
    new_two = hm.hiCMatrix(outfile_two.name)

    nt.assert_equal(test_one.matrix.data, new_one.matrix.data)
    nt.assert_equal(test_one.cut_intervals, new_one.cut_intervals)

    nt.assert_equal(test_two.matrix.data, new_two.matrix.data)
    nt.assert_equal(test_two.cut_intervals, new_two.cut_intervals)

    os.unlink(outfile_one.name)
    os.unlink(outfile_two.name)


def test_normalize_smallest_h5(capsys):
    outfile_one = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_one.close()

    outfile_two = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_two.close()

    args = "--matrices {} {} --normalize smallest -o {} {}".format(matrix_one_h5, matrix_two_h5,
                                                                   outfile_one.name, outfile_two.name).split()
    hicNormalize.main(args)

    test_one = hm.hiCMatrix(ROOT + "/smallest_one.h5")
    test_two = hm.hiCMatrix(ROOT + "/smallest_two.h5")

    new_one = hm.hiCMatrix(outfile_one.name)
    new_two = hm.hiCMatrix(outfile_two.name)

    nt.assert_equal(test_one.matrix.data, new_one.matrix.data)
    nt.assert_equal(test_one.cut_intervals, new_one.cut_intervals)

    nt.assert_equal(test_two.matrix.data, new_two.matrix.data)
    nt.assert_equal(test_two.cut_intervals, new_two.cut_intervals)

    os.unlink(outfile_one.name)
    os.unlink(outfile_two.name)


def test_normalize_smallest_cool(capsys):
    outfile_one = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile_one.close()

    outfile_two = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile_two.close()

    args = "--matrices {} {} --normalize smallest -o {} {}".format(matrix_one_cool, matrix_two_cool,
                                                                   outfile_one.name, outfile_two.name).split()
    hicNormalize.main(args)

    test_one = hm.hiCMatrix(ROOT + "/smallest_one.cool")
    test_two = hm.hiCMatrix(ROOT + "/smallest_two.cool")

    new_one = hm.hiCMatrix(outfile_one.name)
    new_two = hm.hiCMatrix(outfile_two.name)

    nt.assert_equal(test_one.matrix.data, new_one.matrix.data)
    nt.assert_equal(test_one.cut_intervals, new_one.cut_intervals)

    nt.assert_equal(test_two.matrix.data, new_two.matrix.data)
    nt.assert_equal(test_two.cut_intervals, new_two.cut_intervals)

    os.unlink(outfile_one.name)
    os.unlink(outfile_two.name)


def test_normalize_norm_range(capsys):
    outfile_one = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_one.close()

    outfile_two = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_two.close()

    args = "--matrices {} {} --normalize norm_range -o {} {}".format(matrix_one_h5, matrix_two_h5,
                                                                     outfile_one.name, outfile_two.name).split()
    hicNormalize.main(args)

    test_one = hm.hiCMatrix(ROOT + "/norm_range_one.h5")
    test_two = hm.hiCMatrix(ROOT + "/norm_range_two.h5")

    new_one = hm.hiCMatrix(outfile_one.name)
    new_two = hm.hiCMatrix(outfile_two.name)

    nt.assert_equal(test_one.matrix.data, new_one.matrix.data)
    nt.assert_equal(test_one.cut_intervals, new_one.cut_intervals)

    nt.assert_equal(test_two.matrix.data, new_two.matrix.data)
    nt.assert_equal(test_two.cut_intervals, new_two.cut_intervals)

    os.unlink(outfile_one.name)
    os.unlink(outfile_two.name)


def test_normalize_norm_range_cool(capsys):
    outfile_one = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile_one.close()

    outfile_two = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile_two.close()

    args = "--matrices {} {} --normalize norm_range -o {} {}".format(matrix_one_cool, matrix_two_cool,
                                                                     outfile_one.name, outfile_two.name).split()
    hicNormalize.main(args)

    test_one = hm.hiCMatrix(ROOT + "/norm_range_one.cool")
    test_two = hm.hiCMatrix(ROOT + "/norm_range_two.cool")

    new_one = hm.hiCMatrix(outfile_one.name)
    new_two = hm.hiCMatrix(outfile_two.name)

    nt.assert_equal(test_one.matrix.data, new_one.matrix.data)
    nt.assert_equal(test_one.cut_intervals, new_one.cut_intervals)

    nt.assert_equal(test_two.matrix.data, new_two.matrix.data)
    nt.assert_equal(test_two.cut_intervals, new_two.cut_intervals)

    os.unlink(outfile_one.name)
    os.unlink(outfile_two.name)


def test_normalize_norm_range_h5_cool_equal(capsys):
    outfile_one = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile_one.close()

    outfile_two = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_two.close()

    args = "--matrices {} --normalize norm_range -o {}".format(matrix_one_cool,
                                                               outfile_one.name).split()
    hicNormalize.main(args)

    args = "--matrices {} --normalize norm_range -o {}".format(matrix_one_h5,
                                                               outfile_two.name).split()
    hicNormalize.main(args)

    test_one = hm.hiCMatrix(ROOT + "/norm_range_one.cool")
    test_two = hm.hiCMatrix(ROOT + "/norm_range_one.h5")

    new_one = hm.hiCMatrix(outfile_one.name)
    new_two = hm.hiCMatrix(outfile_two.name)

    nt.assert_equal(test_one.matrix.data, new_one.matrix.data)
    nt.assert_equal(test_one.cut_intervals, new_one.cut_intervals)

    nt.assert_equal(test_two.matrix.data, new_two.matrix.data)
    nt.assert_equal(test_two.cut_intervals, new_two.cut_intervals)

    nt.assert_equal(new_one.matrix.data, new_two.matrix.data)
    nt.assert_equal(len(new_one.cut_intervals), len(new_two.cut_intervals))

    os.unlink(outfile_one.name)
    os.unlink(outfile_two.name)


def test_normalize_smallest_h5_cool_equal(capsys):
    outfile_one = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile_one.close()

    outfile_one_cool = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile_one.close()

    outfile_two = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_two.close()
    outfile_two_h5 = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_two.close()

    args = "--matrices {} {} --normalize smallest -o {} {}".format(matrix_one_cool, matrix_two_cool,
                                                                   outfile_one.name, outfile_one_cool.name).split()
    hicNormalize.main(args)

    args = "--matrices {} {} --normalize smallest -o {} {}".format(matrix_one_h5, matrix_two_h5,
                                                                   outfile_two.name, outfile_two_h5.name).split()
    hicNormalize.main(args)

    test_one = hm.hiCMatrix(ROOT + "/smallest_one.cool")
    test_two = hm.hiCMatrix(ROOT + "/smallest_one.h5")

    new_one = hm.hiCMatrix(outfile_one_cool.name)
    new_two = hm.hiCMatrix(outfile_two_h5.name)

    nt.assert_equal(test_one.matrix.data, new_one.matrix.data)
    nt.assert_equal(test_one.cut_intervals, new_one.cut_intervals)

    nt.assert_equal(test_two.matrix.data, new_two.matrix.data)
    nt.assert_equal(test_two.cut_intervals, new_two.cut_intervals)

    nt.assert_equal(new_one.matrix.data, new_two.matrix.data)
    nt.assert_equal(len(new_one.cut_intervals), len(new_two.cut_intervals))

    os.unlink(outfile_one.name)
    os.unlink(outfile_two.name)
