"""Characterization tests for hicMergeTADbins.

There was no test file for this tool before. hicMergeTADbins is the only caller
of reduceMatrix.reduce_matrix with diagonal=True outside hicMergeMatrixBins, and
it deliberately clears correction_factors (hicMergeTADbins.py:87), so both are
pinned here.

The expected matrix is recomputed from the input with numpy.add.at over the
upper triangle instead of the complex-number np.unique plus np.bincount trick
that reduce_matrix uses. That is a genuinely different reduction, and the
comparison below is bit exact on float64 data.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

import os
from tempfile import NamedTemporaryFile

import numpy as np
import numpy.testing as nt
import pytest
from scipy.sparse import triu

from hicmatrix import HiCMatrix as hm
from hicexplorer import hicMergeTADbins

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "test_data")

LI_H5 = os.path.join(ROOT, 'Li_et_al_2015.h5')
LI_COOL = os.path.join(ROOT, 'Li_et_al_2015.cool')
LI_CUT_H5 = os.path.join(ROOT, 'Li_cut.h5')
DOMAINS = os.path.join(ROOT, 'domains.bed')


def run_merge(matrix, suffix):
    outfile = NamedTemporaryFile(suffix=suffix, delete=False)
    outfile.close()
    try:
        hicMergeTADbins.main(['-m', matrix, '--domains', DOMAINS,
                              '-o', outfile.name])
        return hm.hiCMatrix(outfile.name)
    finally:
        if os.path.exists(outfile.name):
            os.unlink(outfile.name)


def groups_from_output(input_cut_intervals, output_cut_intervals):
    """Recover which input bins were merged into each output bin.

    The mapping is read back from the output intervals, so this helper does not
    reimplement hicMergeTADbins' boundary logic. The bin counts are asserted
    separately by the tests, which is what pins that logic.
    """
    groups = []
    for chrom, start, end, _ in output_cut_intervals:
        groups.append([i for i, (c, s, e, _) in enumerate(input_cut_intervals)
                       if c == chrom and s >= start and e <= end])
    return groups


def expected_dense(input_matrix, groups):
    """Sum the upper triangle of every block, then mirror it."""
    upper = triu(input_matrix, k=0, format='coo')
    mapping = np.full(input_matrix.shape[0], -1, dtype=int)
    for k, group in enumerate(groups):
        for bin_id in group:
            mapping[bin_id] = k
    new_row = mapping[upper.row]
    new_col = mapping[upper.col]
    keep = (new_row > -1) & (new_col > -1)
    size = len(groups)
    reduced = np.zeros((size, size), dtype=np.float64)
    np.add.at(reduced, (new_row[keep], new_col[keep]), upper.data[keep])
    return reduced + reduced.T - np.diag(np.diag(reduced))


def restored_input(path):
    matrix = hm.hiCMatrix(path)
    matrix.restoreMaskedBins()
    return matrix


def test_merge_tad_bins_h5_values_are_bit_exact():
    """Li_et_al_2015.h5 (11,104 bins, chromosome X) against domains.bed."""
    original = restored_input(LI_H5)
    result = run_merge(LI_H5, '.h5')

    assert result.matrix.shape == (500, 500)
    assert result.matrix.nnz == 13547
    assert result.matrix.dtype == np.float64

    groups = groups_from_output(original.cut_intervals, result.cut_intervals)
    # Every input bin ends up in exactly one output bin: domains.bed tiles all
    # of chromosome X, so nothing is dropped.
    assert sum(len(g) for g in groups) == 11104
    nt.assert_array_equal(result.matrix.toarray(),
                          expected_dense(original.matrix, groups))


def test_merge_tad_bins_h5_cut_intervals():
    """The first, the last and the count of the merged intervals."""
    result = run_merge(LI_H5, '.h5')
    assert len(result.cut_intervals) == 500
    assert result.cut_intervals[0] == ('X', 0, 72656, 10.32972972972973)
    assert result.cut_intervals[1] == ('X', 72656, 98011, 13.114285714285717)
    assert result.cut_intervals[-1] == ('X', 22421267, 22422827, 47.0)
    assert list(result.chrBinBoundaries.items()) == [('X', (0, 500))]


def test_merge_tad_bins_clears_correction_factors():
    """hicMergeTADbins.py:87 sets correction_factors to None on purpose.

    Note that Li_et_al_2015.h5 carries its stored correction factors in
    distance_counts, not in correction_factors, because the hicmatrix loader
    unpacks the two in the wrong order (finding F10). Both are None on the
    output, which is what the port has to reproduce.
    """
    original = hm.hiCMatrix(LI_H5)
    assert original.correction_factors is None
    assert np.shape(original.distance_counts) == (11104,)

    result = run_merge(LI_H5, '.h5')
    assert result.correction_factors is None


def test_merge_tad_bins_cool_values_are_bit_exact():
    """The cool writer is a separate code path from the h5 writer."""
    original = restored_input(LI_COOL)
    result = run_merge(LI_COOL, '.cool')

    assert result.matrix.shape == (500, 500)
    groups = groups_from_output(original.cut_intervals, result.cut_intervals)
    assert sum(len(g) for g in groups) == original.matrix.shape[0]
    nt.assert_array_equal(result.matrix.toarray(),
                          expected_dense(original.matrix, groups))


def test_merge_tad_bins_conserves_the_upper_triangle_but_not_the_full_sum():
    """The upper triangle is conserved, the symmetric sum is not.

    reduce_matrix sums the upper triangle block by block, so the upper triangle
    of the result carries the same total as the upper triangle of the input.
    The symmetric matrix is then rebuilt as R + R.T - diag(R), and with
    diagonal=True the subtracted diagonal is the whole within-TAD block sum,
    not the original main diagonal. The full symmetric sum therefore drops from
    30,482,969.637651745 to 23,695,524.863065504, a loss of 22 percent, and a
    merged matrix is not comparable to its input by total count. Pinned, not
    fixed.
    """
    original = restored_input(LI_H5)
    result = run_merge(LI_H5, '.h5')

    nt.assert_allclose(triu(result.matrix, k=0).sum(),
                       triu(original.matrix, k=0).sum(), rtol=1e-15, atol=0)
    nt.assert_allclose(triu(original.matrix, k=0).sum(),
                       17548966.536917932, rtol=1e-15, atol=0)
    nt.assert_allclose(original.matrix.sum(), 30482969.637651745,
                       rtol=1e-15, atol=0)
    nt.assert_allclose(result.matrix.sum(), 23695524.863065504,
                       rtol=1e-15, atol=0)


def test_domains_outside_the_matrix_raise_a_type_error():
    """Li_cut.h5 holds only the first 306 bins of chromosome X.

    domains.bed addresses the whole chromosome, so getRegionBinRange returns
    None for the regions past the end, and get_boundary_bin_id unpacks it
    (hicMergeTADbins.py:129). The tool dies with a TypeError rather than a
    diagnostic. Pinned as it stands.
    """
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    os.unlink(outfile.name)
    with pytest.raises(TypeError):
        hicMergeTADbins.main(['-m', LI_CUT_H5, '--domains', DOMAINS,
                              '-o', outfile.name])
    assert not os.path.exists(outfile.name)


def test_expected_model_is_not_vacuous():
    """The reduction really changes the matrix, so the comparison can fail."""
    original = restored_input(LI_H5)
    result = run_merge(LI_H5, '.h5')
    groups = groups_from_output(original.cut_intervals, result.cut_intervals)
    reduced = expected_dense(original.matrix, groups)
    assert reduced.shape == (500, 500)
    assert np.count_nonzero(reduced) == 13547
    assert max(len(g) for g in groups) > 1
