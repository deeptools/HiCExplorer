import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicCompareMatrices
from hicmatrix import HiCMatrix as hm
from tempfile import NamedTemporaryFile
import os
import numpy.testing as nt
from hicexplorer.test.test_compute_function import compute


ROOT = os.path.join(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))), "test_data/")


# I doubled the values in ./hicexplorer/test/test_data/hicConvertFormat/small_test_matrix_chr4.cool
# hicNormalize -m ./hicexplorer/test/test_data/small_test_matrix.cool --normalize multiplicative -mv 2 -o ./hicexplorer/test/test_data/hicCompareMatrices/small_test_matrix_twice.cool

def test_hicCompareMatrices_doubleMinusOneEqual0():
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()

    args = "--matrices {} {} --outFileName {} " \
           "--operation diff" \
           .format(ROOT + "hicCompareMatrices/small_test_matrix_twice.cool",
                   ROOT + "small_test_matrix.cool",
                   outfile.name).split()

    compute(hicCompareMatrices.main, args, 5)
    input = hm.hiCMatrix(
        ROOT + "small_test_matrix.cool")
    new = hm.hiCMatrix(outfile.name)
    nt.assert_equal([], new.matrix.data)
    nt.assert_equal(input.cut_intervals, new.cut_intervals)

    os.unlink(outfile.name)


def test_hicCompareMatrices_noNorm_doubleMinusOneEqualOne():
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()

    args = "--matrices {} {} --outFileName {} " \
           "--operation diff --noNorm" \
           .format(ROOT + "hicCompareMatrices/small_test_matrix_twice.cool",
                   ROOT + "small_test_matrix.cool",
                   outfile.name).split()

    compute(hicCompareMatrices.main, args, 5)
    input = hm.hiCMatrix(
        ROOT + "small_test_matrix.cool")
    new = hm.hiCMatrix(outfile.name)
    nt.assert_equal(input.matrix.data, new.matrix.data)
    nt.assert_equal(input.cut_intervals, new.cut_intervals)

    os.unlink(outfile.name)


# ---------------------------------------------------------------------------
# Characterization tests added for the v4 C++ port.
#
# The two tests above cover only --operation diff, and only on a matrix pair
# that is an exact factor of two apart. The tests below add ratio and
# log2ratio, both with and without --noNorm, and check every stored value.
#
# The expected matrices are recomputed densely with numpy, which is a different
# code path from the sparse elementwise multiply the tool uses, and the
# comparison is bit exact. The dense model deliberately mirrors the tool's
# arithmetic order (a * (1/b), not a/b) because that order is part of the
# behaviour being pinned; reproducing it is what makes exact equality the right
# assertion instead of a tolerance chosen to make the test pass.
# ---------------------------------------------------------------------------
import numpy as np  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402

DTAD = ROOT + "hicDifferentialTAD/"
UNTREATED_H5 = DTAD + "GSM2644945_Untreated-R1.100000_chr1.h5"
AUXIN_H5 = DTAD + "GSM2644947_Auxin2days-R1.100000_chr1.h5"
UNTREATED_COOL = DTAD + "GSM2644945_Untreated-R1.100000_chr1.cool"
AUXIN_COOL = DTAD + "GSM2644947_Auxin2days-R1.100000_chr1.cool"
TWICE_COOL = ROOT + "hicCompareMatrices/small_test_matrix_twice.cool"
SINGLE_COOL = ROOT + "small_test_matrix.cool"


def run_compare(matrix_a, matrix_b, operation, no_norm, suffix):
    outfile = NamedTemporaryFile(suffix=suffix, delete=False)
    outfile.close()
    args = ['--matrices', matrix_a, matrix_b, '--outFileName', outfile.name,
            '--operation', operation]
    if no_norm:
        args.append('--noNorm')
    try:
        hicCompareMatrices.main(args)
        return hm.hiCMatrix(outfile.name)
    finally:
        if os.path.exists(outfile.name):
            os.unlink(outfile.name)


def expected_compare(matrix_a, matrix_b, operation, no_norm):
    """Dense model of hicCompareMatrices.

    Returns the expected csr matrix and the union of the two nan_bins lists.
    """
    hic1 = hm.hiCMatrix(matrix_a)
    hic2 = hm.hiCMatrix(matrix_b)
    dense_a = hic1.matrix.toarray().astype(float)
    dense_b = hic2.matrix.toarray().astype(float)
    if not no_norm:
        dense_a = dense_a / hic1.matrix.data.sum()
        dense_b = dense_b / hic2.matrix.data.sum()

    if operation == 'diff':
        result = dense_a - dense_b
    else:
        # The tool inverts hic2 in place and multiplies, so entries where hic2
        # is zero drop out of the sparse product entirely.
        inverse_b = np.zeros_like(dense_b)
        nonzero = hic2.matrix.toarray() != 0
        inverse_b[nonzero] = 1.0 / dense_b[nonzero]
        result = dense_a * inverse_b
        if operation == 'log2ratio':
            with np.errstate(divide='ignore', invalid='ignore'):
                result = np.where(result != 0, np.log2(result), 0.0)

    nan_union = sorted(set(hic1.nan_bins) | set(hic2.nan_bins))
    result[nan_union, :] = 0
    result[:, nan_union] = 0
    expected = csr_matrix(result)
    expected.eliminate_zeros()
    expected.sort_indices()
    return expected, nan_union


def assert_csr_identical(got, expected):
    got = got.tocsr()
    got.sort_indices()
    assert got.shape == expected.shape
    assert got.nnz == expected.nnz
    nt.assert_array_equal(got.indptr, expected.indptr)
    nt.assert_array_equal(got.indices, expected.indices)
    nt.assert_array_equal(got.data, expected.data)


def test_ratio_noNorm_of_a_doubled_matrix_is_exactly_two():
    """small_test_matrix_twice.cool is small_test_matrix.cool times two."""
    new = run_compare(TWICE_COOL, SINGLE_COOL, 'ratio', True, '.cool')
    assert new.matrix.nnz == 69213
    nt.assert_array_equal(new.matrix.data,
                          np.full(new.matrix.nnz, 2.0))


def test_log2ratio_noNorm_of_a_doubled_matrix_is_exactly_one():
    new = run_compare(TWICE_COOL, SINGLE_COOL, 'log2ratio', True, '.cool')
    assert new.matrix.nnz == 69213
    nt.assert_array_equal(new.matrix.data,
                          np.full(new.matrix.nnz, 1.0))


def test_ratio_with_normalisation_is_one_up_to_one_ulp():
    """After normalisation the two matrices are identical, so every ratio is 1.

    It is not exactly 1 everywhere: 190 of the 69,213 entries come out as
    0.9999999999999999, one ulp below one, because the tool computes
    a * (1 / b) rather than a / b. The exact counts are pinned so that a port
    which reorders the division is detected.
    """
    new = run_compare(TWICE_COOL, SINGLE_COOL, 'ratio', False, '.cool')
    values, counts = np.unique(new.matrix.data, return_counts=True)
    nt.assert_array_equal(values, np.array([0.9999999999999999, 1.0]))
    nt.assert_array_equal(counts, np.array([190, 69023]))


def test_log2ratio_with_normalisation_keeps_only_the_rounding_residue():
    """log2 of the ratio above is zero except for the 190 one-ulp entries.

    eliminate_zeros then throws away the 69,023 exact zeros, so the written
    matrix holds 190 values of -1.6017132519074588e-16 and the file loads back
    with 33,568 nan bins instead of the 14,845 of the inputs. Pinned as it is.
    """
    new = run_compare(TWICE_COOL, SINGLE_COOL, 'log2ratio', False, '.cool')
    assert new.matrix.nnz == 190
    nt.assert_array_equal(new.matrix.data,
                          np.full(190, -1.6017132519074588e-16))
    assert len(new.nan_bins) == 33568


def test_diff_h5_is_bit_exact():
    expected, nan_union = expected_compare(UNTREATED_H5, AUXIN_H5, 'diff', False)
    new = run_compare(UNTREATED_H5, AUXIN_H5, 'diff', False, '.h5')
    assert new.matrix.nnz == 3152621
    assert_csr_identical(new.matrix, expected)
    nt.assert_array_equal(sorted(new.nan_bins), nan_union)


def test_diff_noNorm_h5_is_bit_exact():
    expected, _ = expected_compare(UNTREATED_H5, AUXIN_H5, 'diff', True)
    new = run_compare(UNTREATED_H5, AUXIN_H5, 'diff', True, '.h5')
    assert_csr_identical(new.matrix, expected)


def test_ratio_h5_is_bit_exact():
    expected, _ = expected_compare(UNTREATED_H5, AUXIN_H5, 'ratio', False)
    new = run_compare(UNTREATED_H5, AUXIN_H5, 'ratio', False, '.h5')
    # The product keeps only the entries where both matrices are non zero,
    # 2,277,549 of the 3,152,621 entries of the difference.
    assert new.matrix.nnz == 2277549
    assert_csr_identical(new.matrix, expected)


def test_ratio_noNorm_h5_is_bit_exact():
    expected, _ = expected_compare(UNTREATED_H5, AUXIN_H5, 'ratio', True)
    new = run_compare(UNTREATED_H5, AUXIN_H5, 'ratio', True, '.h5')
    assert_csr_identical(new.matrix, expected)


def test_log2ratio_h5_is_bit_exact():
    expected, _ = expected_compare(UNTREATED_H5, AUXIN_H5, 'log2ratio', False)
    new = run_compare(UNTREATED_H5, AUXIN_H5, 'log2ratio', False, '.h5')
    assert new.matrix.nnz == 2277549
    assert_csr_identical(new.matrix, expected)


def test_log2ratio_noNorm_h5_is_bit_exact():
    expected, _ = expected_compare(UNTREATED_H5, AUXIN_H5, 'log2ratio', True,)
    new = run_compare(UNTREATED_H5, AUXIN_H5, 'log2ratio', True, '.h5')
    assert_csr_identical(new.matrix, expected)


def test_log2ratio_cool_is_bit_exact():
    """cool output takes a different writer than h5."""
    expected, nan_union = expected_compare(UNTREATED_COOL, AUXIN_COOL,
                                           'log2ratio', False)
    new = run_compare(UNTREATED_COOL, AUXIN_COOL, 'log2ratio', False, '.cool')
    assert_csr_identical(new.matrix, expected)
    nt.assert_array_equal(sorted(new.nan_bins), nan_union)


def test_normalisation_divides_by_the_sum_of_the_stored_values():
    """The normaliser is the sum over the stored symmetric data, not the
    upper triangle, and it is pinned at full precision.
    """
    hic1 = hm.hiCMatrix(UNTREATED_H5)
    hic2 = hm.hiCMatrix(AUXIN_H5)
    nt.assert_allclose(hic1.matrix.data.sum(), 1514.2970371802683,
                       rtol=1e-15, atol=0)
    nt.assert_allclose(hic2.matrix.data.sum(), 1385.1134232101283,
                       rtol=1e-15, atol=0)


def test_expected_model_is_not_vacuous():
    """diff, ratio and log2ratio must produce three different matrices."""
    diff, _ = expected_compare(UNTREATED_H5, AUXIN_H5, 'diff', False)
    ratio, _ = expected_compare(UNTREATED_H5, AUXIN_H5, 'ratio', False)
    log2, _ = expected_compare(UNTREATED_H5, AUXIN_H5, 'log2ratio', False)
    assert diff.nnz != ratio.nnz
    assert ratio.nnz == log2.nnz
    assert not np.array_equal(ratio.data, log2.data)
