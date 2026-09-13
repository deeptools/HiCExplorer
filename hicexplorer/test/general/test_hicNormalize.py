import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicmatrix import HiCMatrix as hm
from hicexplorer import hicNormalize

from tempfile import NamedTemporaryFile

import os
import numpy as np
import numpy.testing as nt

from hicexplorer.test.test_compute_function import compute

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/hicNormalize")

matrix_one_h5 = ROOT + '/small_test_matrix.h5'
matrix_two_h5 = ROOT + '/small_test_matrix_scaled_up.h5'

matrix_one_cool = ROOT + '/small_test_matrix.cool'
matrix_two_cool = ROOT + '/small_test_matrix_scaled_up.cool'


def assert_matrix_identical(pExpected, pObserved):
    """Full CSR identity, not only the value array.

    The original version of this file compared `matrix.data` and
    `cut_intervals` and nothing else. That leaves the sparsity pattern, the
    matrix shape, the stored dtype and the NaN bin list unconstrained, so a
    port that scaled the right numbers into the wrong cells, dropped the
    float32 cast of hicNormalize.py:76,105,131, or wrote a different bin table
    would still pass. All four are observable in the output file, so all four
    are asserted here.

    `data` is compared with assert_equal, i.e. bit exactly, because
    hicNormalize does nothing but elementwise scaling: there is no reduction
    whose order could legitimately differ.
    """
    nt.assert_equal(pExpected.matrix.shape, pObserved.matrix.shape)
    nt.assert_equal(pExpected.matrix.dtype, pObserved.matrix.dtype)
    nt.assert_equal(pExpected.matrix.nnz, pObserved.matrix.nnz)
    nt.assert_equal(pExpected.matrix.data, pObserved.matrix.data)
    nt.assert_equal(pExpected.matrix.indices, pObserved.matrix.indices)
    nt.assert_equal(pExpected.matrix.indptr, pObserved.matrix.indptr)
    nt.assert_equal(pExpected.cut_intervals, pObserved.cut_intervals)
    nt.assert_equal(pExpected.nan_bins, pObserved.nan_bins)


def assert_thresholded(pReference, pObserved, pThreshold):
    """--setToZeroThreshold t is exactly `drop every value below t`.

    hicNormalize.py:96-98, :123-125 and :149-151 apply the threshold after the
    normalization is complete and after eliminate_zeros, so the result is the
    untresholded output with the entries below the threshold removed and
    nothing else changed. Expressing the expectation that way pins the
    threshold path against an existing master instead of against a new binary
    file, and it also pins the comparison operator: the mask is `< t`, so a
    value exactly equal to t survives.
    """
    keep = pReference.matrix.data >= pThreshold
    nt.assert_equal(pReference.matrix.data[keep], pObserved.matrix.data)
    nt.assert_equal(pReference.matrix.indices[keep], pObserved.matrix.indices)
    nt.assert_equal(pObserved.matrix.dtype, np.float32)
    nt.assert_equal(pReference.matrix.shape, pObserved.matrix.shape)
    nt.assert_equal(pReference.cut_intervals, pObserved.cut_intervals)


def test_normalize_smallest(capsys):
    outfile_one = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_one.close()

    outfile_two = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_two.close()

    args = "--matrices {} {} --normalize smallest -o {} {}".format(matrix_one_h5, matrix_two_h5,
                                                                   outfile_one.name, outfile_two.name).split()
    # hicNormalize.main(args)
    compute(hicNormalize.main, args, 5)
    test_one = hm.hiCMatrix(ROOT + "/smallest_one.h5")
    test_two = hm.hiCMatrix(ROOT + "/smallest_two.h5")

    new_one = hm.hiCMatrix(outfile_one.name)
    new_two = hm.hiCMatrix(outfile_two.name)

    assert_matrix_identical(test_one, new_one)
    assert_matrix_identical(test_two, new_two)

    os.unlink(outfile_one.name)
    os.unlink(outfile_two.name)


def test_normalize_smallest_h5(capsys):
    outfile_one = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_one.close()

    outfile_two = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_two.close()

    args = "--matrices {} {} --normalize smallest -o {} {}".format(matrix_one_h5, matrix_two_h5,
                                                                   outfile_one.name, outfile_two.name).split()
    # hicNormalize.main(args)
    compute(hicNormalize.main, args, 5)

    test_one = hm.hiCMatrix(ROOT + "/smallest_one.h5")
    test_two = hm.hiCMatrix(ROOT + "/smallest_two.h5")

    new_one = hm.hiCMatrix(outfile_one.name)
    new_two = hm.hiCMatrix(outfile_two.name)

    assert_matrix_identical(test_one, new_one)
    assert_matrix_identical(test_two, new_two)

    os.unlink(outfile_one.name)
    os.unlink(outfile_two.name)


def test_normalize_smallest_cool(capsys):
    outfile_one = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile_one.close()

    outfile_two = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile_two.close()

    args = "--matrices {} {} --normalize smallest -o {} {}".format(matrix_one_cool, matrix_two_cool,
                                                                   outfile_one.name, outfile_two.name).split()
    # hicNormalize.main(args)
    compute(hicNormalize.main, args, 5)

    test_one = hm.hiCMatrix(ROOT + "/smallest_one.cool")
    test_two = hm.hiCMatrix(ROOT + "/smallest_two.cool")

    new_one = hm.hiCMatrix(outfile_one.name)
    new_two = hm.hiCMatrix(outfile_two.name)

    assert_matrix_identical(test_one, new_one)
    assert_matrix_identical(test_two, new_two)

    os.unlink(outfile_one.name)
    os.unlink(outfile_two.name)


def test_normalize_norm_range(capsys):
    outfile_one = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_one.close()

    outfile_two = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_two.close()

    args = "--matrices {} {} --normalize norm_range -o {} {}".format(matrix_one_h5, matrix_two_h5,
                                                                     outfile_one.name, outfile_two.name).split()
    # hicNormalize.main(args)
    compute(hicNormalize.main, args, 5)

    test_one = hm.hiCMatrix(ROOT + "/norm_range_one.h5")
    test_two = hm.hiCMatrix(ROOT + "/norm_range_two.h5")

    new_one = hm.hiCMatrix(outfile_one.name)
    new_two = hm.hiCMatrix(outfile_two.name)

    assert_matrix_identical(test_one, new_one)
    assert_matrix_identical(test_two, new_two)

    os.unlink(outfile_one.name)
    os.unlink(outfile_two.name)


def test_normalize_norm_range_cool(capsys):
    outfile_one = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile_one.close()

    outfile_two = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile_two.close()

    args = "--matrices {} {} --normalize norm_range -o {} {}".format(matrix_one_cool, matrix_two_cool,
                                                                     outfile_one.name, outfile_two.name).split()
    # hicNormalize.main(args)
    compute(hicNormalize.main, args, 5)

    test_one = hm.hiCMatrix(ROOT + "/norm_range_one.cool")
    test_two = hm.hiCMatrix(ROOT + "/norm_range_two.cool")

    new_one = hm.hiCMatrix(outfile_one.name)
    new_two = hm.hiCMatrix(outfile_two.name)

    assert_matrix_identical(test_one, new_one)
    assert_matrix_identical(test_two, new_two)

    os.unlink(outfile_one.name)
    os.unlink(outfile_two.name)


def test_normalize_norm_range_h5_cool_equal(capsys):
    outfile_one = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile_one.close()

    outfile_two = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile_two.close()

    args = "--matrices {} --normalize norm_range -o {}".format(matrix_one_cool,
                                                               outfile_one.name).split()
    # hicNormalize.main(args)
    compute(hicNormalize.main, args, 5)

    args = "--matrices {} --normalize norm_range -o {}".format(matrix_one_h5,
                                                               outfile_two.name).split()
    # hicNormalize.main(args)
    compute(hicNormalize.main, args, 5)

    test_one = hm.hiCMatrix(ROOT + "/norm_range_one.cool")
    test_two = hm.hiCMatrix(ROOT + "/norm_range_one.h5")

    new_one = hm.hiCMatrix(outfile_one.name)
    new_two = hm.hiCMatrix(outfile_two.name)

    assert_matrix_identical(test_one, new_one)
    assert_matrix_identical(test_two, new_two)

    # The cross-format check: the same matrix in cool and in h5 must normalize
    # to the same values in the same cells. Only the value array was compared
    # before, which a permutation of the pattern would have survived.
    nt.assert_equal(new_one.matrix.data, new_two.matrix.data)
    nt.assert_equal(new_one.matrix.indices, new_two.matrix.indices)
    nt.assert_equal(new_one.matrix.indptr, new_two.matrix.indptr)
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
    # hicNormalize.main(args)
    compute(hicNormalize.main, args, 5)

    args = "--matrices {} {} --normalize smallest -o {} {}".format(matrix_one_h5, matrix_two_h5,
                                                                   outfile_two.name, outfile_two_h5.name).split()
    # hicNormalize.main(args)
    compute(hicNormalize.main, args, 5)

    test_one = hm.hiCMatrix(ROOT + "/smallest_one.cool")
    test_two = hm.hiCMatrix(ROOT + "/smallest_one.h5")

    new_one = hm.hiCMatrix(outfile_one_cool.name)
    new_two = hm.hiCMatrix(outfile_two_h5.name)

    # The second output of each run is the scaled-up matrix divided by its own
    # adjust factor, so it must come back equal to the first matrix's master.
    assert_matrix_identical(test_one, new_one)
    assert_matrix_identical(test_two, new_two)

    nt.assert_equal(new_one.matrix.data, new_two.matrix.data)
    nt.assert_equal(new_one.matrix.indices, new_two.matrix.indices)
    nt.assert_equal(new_one.matrix.indptr, new_two.matrix.indptr)
    nt.assert_equal(len(new_one.cut_intervals), len(new_two.cut_intervals))

    os.unlink(outfile_one.name)
    os.unlink(outfile_two.name)


def test_normalize_multiplicative_h5_cool(capsys):
    outfile_one = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile_one.close()

    args = "--matrices {} --normalize multiplicative --multiplicativeValue {} -o {}".format(matrix_one_cool, 2,
                                                                                            outfile_one.name).split()
    # hicNormalize.main(args)
    compute(hicNormalize.main, args, 5)

    test_one = hm.hiCMatrix(ROOT + "/small_test_matrix_scaled_by_2.cool")

    new_one = hm.hiCMatrix(outfile_one.name)

    assert_matrix_identical(test_one, new_one)

    os.unlink(outfile_one.name)


def test_normalize_multiplicative_h5(capsys):
    """multiplicative on an h5 input, which the file did not cover.

    The cool master is the only multiplicative reference in the repository, so
    the h5 path is pinned against it: the two inputs hold the same matrix, and
    scaling is elementwise, so the h5 output must carry the same values in the
    same cells. Only the bin tables differ, because the cool loader derives NaN
    bins and the h5 loader does not.
    """
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()

    args = "--matrices {} --normalize multiplicative --multiplicativeValue {} -o {}".format(matrix_one_h5, 2,
                                                                                            outfile.name).split()
    compute(hicNormalize.main, args, 5)

    reference = hm.hiCMatrix(ROOT + "/small_test_matrix_scaled_by_2.cool")
    observed = hm.hiCMatrix(outfile.name)

    nt.assert_equal(reference.matrix.shape, observed.matrix.shape)
    nt.assert_equal(observed.matrix.dtype, np.float32)
    nt.assert_equal(reference.matrix.data, observed.matrix.data)
    nt.assert_equal(reference.matrix.indices, observed.matrix.indices)
    nt.assert_equal(reference.matrix.indptr, observed.matrix.indptr)

    os.unlink(outfile.name)


def test_normalize_multiplicative_default_value(capsys):
    """--multiplicativeValue defaults to 1, so the tool only casts to float32.

    This pins the default, which no test exercised, and with it the fact that
    the multiplicative branch still rewrites the dtype and still runs
    eliminate_zeros even when it scales by one.
    """
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()

    args = "--matrices {} --normalize multiplicative -o {}".format(matrix_one_h5, outfile.name).split()
    compute(hicNormalize.main, args, 5)

    source = hm.hiCMatrix(matrix_one_h5)
    observed = hm.hiCMatrix(outfile.name)

    nt.assert_equal(observed.matrix.dtype, np.float32)
    nt.assert_equal(source.matrix.data.astype(np.float32), observed.matrix.data)
    nt.assert_equal(source.matrix.indices, observed.matrix.indices)
    nt.assert_equal(source.matrix.indptr, observed.matrix.indptr)
    nt.assert_equal(source.cut_intervals, observed.cut_intervals)

    os.unlink(outfile.name)


def test_normalize_smallest_set_to_zero_threshold(capsys):
    """--setToZeroThreshold on the smallest mode, which no test exercised.

    The smallest master holds values from 1.0 to 8.0, so a threshold of 2.0
    removes a large, well-defined part of the matrix: 69,213 stored entries
    fall to 1,970. Asserting against the master filtered by the same rule
    checks the threshold, the comparison operator and the eliminate_zeros that
    follows it, all three exactly.
    """
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()

    args = "--matrices {} --normalize smallest -o {} --setToZeroThreshold 2.0".format(
        matrix_one_h5, outfile.name).split()
    compute(hicNormalize.main, args, 5)

    reference = hm.hiCMatrix(ROOT + "/smallest_one.h5")
    observed = hm.hiCMatrix(outfile.name)

    nt.assert_equal(observed.matrix.nnz, 1970)
    assert_thresholded(reference, observed, 2.0)

    os.unlink(outfile.name)


def test_normalize_norm_range_set_to_zero_threshold(capsys):
    """--setToZeroThreshold on the norm_range mode, through the cool writer.

    norm_range maps the values into [0, 1], so a threshold of 0.5 keeps only
    the upper half of the range: 1,970 stored entries fall to 23. The short
    form -sz is used here and the long form above, so both spellings are
    covered.
    """
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()

    args = "--matrices {} --normalize norm_range -o {} -sz 0.5".format(
        matrix_one_cool, outfile.name).split()
    compute(hicNormalize.main, args, 5)

    reference = hm.hiCMatrix(ROOT + "/norm_range_one.cool")
    observed = hm.hiCMatrix(outfile.name)

    nt.assert_equal(observed.matrix.nnz, 23)
    assert_thresholded(reference, observed, 0.5)

    os.unlink(outfile.name)
