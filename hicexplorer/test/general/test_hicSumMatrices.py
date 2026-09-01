"""Characterization tests for hicSumMatrices.

There was no test file for this tool before. The tests below pin the exact
numbers the 3.7.6 Python implementation produces on the real matrices in
hicexplorer/test/test_data, so that the C++ port can be checked against a fixed
reference.

The expected matrices are not taken from a stored reference file. They are
recomputed in the test from the inputs with a deliberately different code path
(dense masking instead of hicmatrix's maskBins/restoreMaskedBins), and the
comparison is bit exact, so the assertions constrain every stored value.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

import os
from tempfile import NamedTemporaryFile

import numpy as np
import numpy.testing as nt
import pytest
from scipy.sparse import csr_matrix

from hicmatrix import HiCMatrix as hm
from hicexplorer import hicSumMatrices

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "test_data")
DTAD = os.path.join(ROOT, 'hicDifferentialTAD')

SMALL_H5 = os.path.join(ROOT, 'small_test_matrix_50kb_res.h5')
UNTREATED_H5 = os.path.join(DTAD, 'GSM2644945_Untreated-R1.100000_chr1.h5')
AUXIN_H5 = os.path.join(DTAD, 'GSM2644947_Auxin2days-R1.100000_chr1.h5')
UNTREATED_COOL = os.path.join(DTAD, 'GSM2644945_Untreated-R1.100000_chr1.cool')
AUXIN_COOL = os.path.join(DTAD, 'GSM2644947_Auxin2days-R1.100000_chr1.cool')
LI_H5 = os.path.join(ROOT, 'Li_et_al_2015.h5')


def run_sum(matrices, suffix):
    """Run hicSumMatrices and return the loaded result matrix."""
    outfile = NamedTemporaryFile(suffix=suffix, delete=False)
    outfile.close()
    try:
        hicSumMatrices.main(['-m'] + list(matrices) + ['-o', outfile.name])
        return hm.hiCMatrix(outfile.name)
    finally:
        if os.path.exists(outfile.name):
            os.unlink(outfile.name)


def expected_sum(path_a, path_b):
    """Independent model of what hicSumMatrices writes.

    hicSumMatrices adds the two sparse matrices, then calls maskBins on the
    union of the two nan_bins lists. maskBins deletes those rows and columns;
    save() restores them as empty rows. The net effect on the stored values is
    that every entry that lives in a nan bin of either input is dropped. This
    is modelled here by zeroing those rows and columns in a dense copy, which
    is a different code path from the one the tool takes.
    """
    a = hm.hiCMatrix(path_a)
    b = hm.hiCMatrix(path_b)
    nan_union = sorted(set(a.nan_bins) | set(b.nan_bins))
    summed = (a.matrix + b.matrix).tolil()
    summed[nan_union, :] = 0
    summed[:, nan_union] = 0
    summed = summed.tocsr()
    summed.eliminate_zeros()
    summed.sort_indices()
    return summed, nan_union, a


def assert_csr_identical(got, expected):
    got = got.tocsr()
    got.sort_indices()
    assert got.shape == expected.shape
    assert got.nnz == expected.nnz
    nt.assert_array_equal(got.indptr, expected.indptr)
    nt.assert_array_equal(got.indices, expected.indices)
    nt.assert_array_equal(got.data, expected.data)


def test_sum_matrix_with_itself_doubles_every_value():
    """Adding a matrix to itself must double every stored value exactly."""
    original = hm.hiCMatrix(SMALL_H5)
    result = run_sum([SMALL_H5, SMALL_H5], '.h5')

    assert result.matrix.shape == (2794, 2794)
    assert result.matrix.nnz == 45626
    nt.assert_array_equal(result.matrix.data, 2 * original.matrix.data)
    nt.assert_equal(result.cut_intervals, original.cut_intervals)
    nt.assert_array_equal(sorted(result.nan_bins), sorted(original.nan_bins))


def test_sum_of_an_integer_matrix_is_written_as_float64():
    """small_test_matrix_50kb_res.h5 stores int64, the sum stores float64.

    The upcast comes from maskBins/restoreMaskedBins in hicmatrix, not from the
    addition itself. Pinned because the C++ writer has to choose the same dtype.
    """
    original = hm.hiCMatrix(SMALL_H5)
    assert original.matrix.dtype == np.int64
    result = run_sum([SMALL_H5, SMALL_H5], '.h5')
    assert result.matrix.dtype == np.float64


def test_sum_two_h5_matrices_is_bit_exact():
    expected, nan_union, _ = expected_sum(UNTREATED_H5, AUXIN_H5)
    result = run_sum([UNTREATED_H5, AUXIN_H5], '.h5')

    assert_csr_identical(result.matrix, expected)
    nt.assert_array_equal(sorted(result.nan_bins), nan_union)


def test_sum_drops_every_entry_that_lives_in_a_nan_bin():
    """The masking step silently discards counts.

    GSM2644945 has 83 nan bins, GSM2644947 has 81, and the union is the 83 of
    the first matrix. The plain sparse sum has 3,157,763 stored entries; the
    file hicSumMatrices writes has 3,152,621, so 5,142 entries that sit in a
    nan row or column are lost. This is pinned, not fixed.
    """
    a = hm.hiCMatrix(UNTREATED_H5)
    b = hm.hiCMatrix(AUXIN_H5)
    assert len(a.nan_bins) == 83
    assert len(b.nan_bins) == 81
    assert len(set(a.nan_bins) | set(b.nan_bins)) == 83

    plain = (a.matrix + b.matrix).tocsr()
    plain.eliminate_zeros()
    assert plain.nnz == 3157763

    result = run_sum([UNTREATED_H5, AUXIN_H5], '.h5')
    assert result.matrix.nnz == 3152621
    assert len(result.nan_bins) == 83


def test_sum_two_cool_matrices_is_bit_exact():
    """cool goes through a different writer than h5, so it is checked too."""
    expected, nan_union, _ = expected_sum(UNTREATED_COOL, AUXIN_COOL)
    result = run_sum([UNTREATED_COOL, AUXIN_COOL], '.cool')

    assert_csr_identical(result.matrix, expected)
    nt.assert_array_equal(sorted(result.nan_bins), nan_union)


def test_sum_of_three_matrices_is_bit_exact():
    """The accumulation loop is run more than once."""
    a = hm.hiCMatrix(UNTREATED_H5)
    b = hm.hiCMatrix(AUXIN_H5)
    nan_union = sorted(set(a.nan_bins) | set(b.nan_bins))
    summed = (a.matrix + b.matrix + a.matrix).tolil()
    summed[nan_union, :] = 0
    summed[:, nan_union] = 0
    summed = summed.tocsr()
    summed.eliminate_zeros()
    summed.sort_indices()

    result = run_sum([UNTREATED_H5, AUXIN_H5, UNTREATED_H5], '.h5')
    assert_csr_identical(result.matrix, summed)


def test_matrices_with_a_different_chromosome_layout_exit_1():
    """The chrBinBoundaries check at hicSumMatrices.py:53 calls exit(1)."""
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    os.unlink(outfile.name)
    with pytest.raises(SystemExit) as excinfo:
        hicSumMatrices.main(['-m', LI_H5, UNTREATED_H5, '-o', outfile.name])
    assert excinfo.value.code == 1
    # Nothing is written when the check fires.
    assert not os.path.exists(outfile.name)


def test_the_shape_mismatch_branch_is_unreachable_for_valid_matrices():
    """hicSumMatrices.py:60-67 guards against a shape mismatch, but the
    chrBinBoundaries comparison above it already rejects every pair of matrices
    whose shapes differ, because the boundaries encode the shape. The branch is
    therefore dead for any file pair that can be loaded. Pinned so the port does
    not have to reproduce an error message that cannot be triggered.
    """
    a = hm.hiCMatrix(LI_H5)
    b = hm.hiCMatrix(UNTREATED_H5)
    assert a.matrix.shape != b.matrix.shape
    assert a.chrBinBoundaries != b.chrBinBoundaries


def test_summing_a_single_matrix_rewrites_it_through_the_mask_step():
    """One input matrix means no addition at all, but the mask step still runs.

    small_test_matrix_50kb_res.h5 has 139 nan bins that hold no entries, so the
    values survive unchanged while the dtype still becomes float64.
    """
    original = hm.hiCMatrix(SMALL_H5)
    result = run_sum([SMALL_H5], '.h5')
    nt.assert_array_equal(result.matrix.data, original.matrix.data.astype(np.float64))
    nt.assert_array_equal(result.matrix.indices, original.matrix.indices)
    nt.assert_array_equal(result.matrix.indptr, original.matrix.indptr)
    nt.assert_equal(result.cut_intervals, original.cut_intervals)


def test_expected_model_is_not_vacuous():
    """Guard for the helper used by the tests above.

    expected_sum() must produce something that differs from either input,
    otherwise a broken tool that copied its first input would pass.
    """
    expected, _, first = expected_sum(UNTREATED_H5, AUXIN_H5)
    assert expected.nnz != first.matrix.nnz
    assert isinstance(expected, csr_matrix)
