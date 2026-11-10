"""Characterization tests for hicSumMatrices.

The tests below use the real matrices in hicexplorer/test/test_data. 3.7.6
masked the bins that are empty in either input after the addition, which
discarded the counts of those bins (#814). HiCExplorer 4 keeps every pixel of
every input, so the sum of the result equals the sum of the inputs.

The expected matrices are not taken from a stored reference file. They are
recomputed in the test from the inputs with the plain sparse addition, and the
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
    """Independent model of what hicSumMatrices writes: the plain sparse sum.

    3.7.6 zeroed every row and column that is a nan bin of either input (it
    called maskBins on the union); that lost counts (#814) and is not modelled.
    """
    a = hm.hiCMatrix(path_a)
    b = hm.hiCMatrix(path_b)
    summed = (a.matrix + b.matrix).tocsr()
    summed.eliminate_zeros()
    summed.sort_indices()
    empty = np.flatnonzero(np.diff(summed.indptr) == 0).tolist()
    return summed, empty, a


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


def test_sum_of_an_integer_matrix_stays_integer():
    """small_test_matrix_50kb_res.h5 stores int64, and so does the sum.

    3.7.6 wrote float64, because restoreMaskedBins pads with a float64 block.
    The masking is gone, so the dtype of the addition is kept.
    """
    original = hm.hiCMatrix(SMALL_H5)
    assert original.matrix.dtype == np.int64
    result = run_sum([SMALL_H5, SMALL_H5], '.h5')
    assert result.matrix.dtype == np.int64


def test_sum_two_h5_matrices_is_bit_exact():
    expected, _, _ = expected_sum(UNTREATED_H5, AUXIN_H5)
    result = run_sum([UNTREATED_H5, AUXIN_H5], '.h5')

    assert_csr_identical(result.matrix, expected)


def test_sum_keeps_every_entry_that_lives_in_a_nan_bin():
    """No count is lost in a bin that is empty in one input (#814).

    GSM2644945 has 83 nan bins and GSM2644947 has 81. The plain sparse sum has
    3,157,763 stored entries. 3.7.6 masked the union of the nan bins and wrote
    3,152,621, so 5,142 entries that sit in a nan row or column of one input
    were lost. All of them are written now.
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
    assert result.matrix.nnz == 3157763
    assert result.matrix.sum() == a.matrix.sum() + b.matrix.sum()


def test_sum_of_two_halves_of_a_cool_file_is_the_original(tmp_path):
    """Split the pixels of a matrix into two cool files and sum them (#814)."""
    import cooler
    source = cooler.Cooler(UNTREATED_COOL)
    pixels = source.pixels()[:]
    bins = source.bins()[:][['chrom', 'start', 'end']]
    paths = []
    for part in (0, 1):
        path = str(tmp_path / 'part{}.cool'.format(part))
        cooler.create_cooler(path, bins, pixels.iloc[part::2].reset_index(drop=True),
                             dtypes={'count': 'float64'})
        paths.append(path)
    out = str(tmp_path / 'sum.cool')
    hicSumMatrices.main(['-m'] + paths + ['-o', out])

    result = cooler.Cooler(out).pixels()[:]
    merged = pixels.merge(result, on=['bin1_id', 'bin2_id'], how='outer',
                          suffixes=('_original', '_sum'))
    assert len(merged) == len(pixels) == len(result)
    assert (merged['count_original'] == merged['count_sum']).all()
    assert result['count'].sum() == pixels['count'].sum()


def test_sum_two_cool_matrices_is_bit_exact():
    """cool goes through a different writer than h5, so it is checked too."""
    expected, _, _ = expected_sum(UNTREATED_COOL, AUXIN_COOL)
    result = run_sum([UNTREATED_COOL, AUXIN_COOL], '.cool')

    assert_csr_identical(result.matrix, expected)


def test_sum_of_three_matrices_is_bit_exact():
    """The accumulation loop is run more than once."""
    a = hm.hiCMatrix(UNTREATED_H5)
    b = hm.hiCMatrix(AUXIN_H5)
    summed = (a.matrix + b.matrix + a.matrix).tocsr()
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


def test_summing_a_single_matrix_rewrites_it():
    """One input matrix means no addition at all; the matrix is written as it is.

    small_test_matrix_50kb_res.h5 has 139 nan bins that hold no entries.
    """
    original = hm.hiCMatrix(SMALL_H5)
    result = run_sum([SMALL_H5], '.h5')
    nt.assert_array_equal(result.matrix.data, original.matrix.data)
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
    assert expected.sum() > first.matrix.sum()
    assert isinstance(expected, csr_matrix)
