import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicMergeMatrixBins
from hicmatrix import HiCMatrix as hm
from tempfile import NamedTemporaryFile
import os
import numpy.testing as nt
from hicexplorer.test.test_compute_function import compute


ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")


def test_correct_matrix():
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()

    args = "--matrix {} --numBins 5 " \
        " --outFileName {}".format(ROOT + "small_test_matrix.h5",
                                   outfile.name).split()
    # hicMergeMatrixBins.main(args)
    compute(hicMergeMatrixBins.main, args, 5)
    test = hm.hiCMatrix(ROOT + "hicMergeMatrixBins/result.h5")
    new = hm.hiCMatrix(outfile.name)
    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)

    os.unlink(outfile.name)


# ---------------------------------------------------------------------------
# Characterization tests added for the v4 C++ port.
#
# The test above covers one merge_bins run against a stored reference file.
# --runningWindow, the cool writer, the chromosome-skipping rule and the nan
# handling were not covered at all. hicConvertFormat calls merge_bins directly
# for its mcool output, so a regression here breaks that tool too.
#
# The expected matrices are recomputed with numpy from the input. For
# merge_bins the block sums are accumulated with numpy.add.at over the upper
# triangle, which is a different reduction from the complex-number np.unique
# plus np.bincount trick in reduceMatrix.reduce_matrix. For --runningWindow the
# window is applied by dense array shifts instead of by concatenating shifted
# coo index arrays. Both comparisons are bit exact.
# ---------------------------------------------------------------------------
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from scipy.sparse import triu  # noqa: E402

SMALL_50KB_H5 = ROOT + "small_test_matrix_50kb_res.h5"
SMALL_50KB_COOL = ROOT + "small_test_matrix_50kb_res.cool"
SMALL_H5 = ROOT + "small_test_matrix.h5"


def run_merge(matrix, num_bins, running_window, suffix):
    outfile = NamedTemporaryFile(suffix=suffix, delete=False)
    outfile.close()
    args = ['--matrix', matrix, '--outFileName', outfile.name,
            '--numBins', str(num_bins)]
    if running_window:
        args.append('--runningWindow')
    try:
        hicMergeMatrixBins.main(args)
        return hm.hiCMatrix(outfile.name)
    finally:
        if os.path.exists(outfile.name):
            os.unlink(outfile.name)


def nan_free_input(path):
    """What remove_nans_if_needed leaves behind: the matrix without nan bins.

    Returned as a dense array plus the surviving cut intervals.
    """
    original = hm.hiCMatrix(path)
    dropped = set(original.nan_bins)
    keep = [i for i in range(original.matrix.shape[0]) if i not in dropped]
    sub = original.matrix.tocsr()[keep, :][:, keep]
    intervals = [original.cut_intervals[i] for i in keep]
    return sub, intervals


def dense_running_window(dense, num_bins):
    """Dense model of running_window_merge.

    Every entry of the upper triangle is added into every cell of the
    num_bins x num_bins window centred on it, out-of-range positions are
    dropped, the result is folded back to the upper triangle and mirrored.
    """
    half = (num_bins - 1) // 2
    size = dense.shape[0]
    upper = np.triu(dense)
    accumulated = np.zeros_like(upper)
    for row_shift in range(-half, half + 1):
        for col_shift in range(-half, half + 1):
            r0, r1 = max(0, row_shift), min(size, size + row_shift)
            c0, c1 = max(0, col_shift), min(size, size + col_shift)
            accumulated[r0:r1, c0:c1] += upper[r0 - row_shift:r1 - row_shift,
                                               c0 - col_shift:c1 - col_shift]
    folded = np.triu(accumulated)
    return folded + folded.T - np.diag(np.diag(folded))


def dense_block_reduce(sparse_input, groups):
    upper = triu(sparse_input, k=0, format='coo')
    mapping = np.full(sparse_input.shape[0], -1, dtype=int)
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


def groups_from_output(input_intervals, output_intervals):
    """Recover which input bins ended up in each output bin.

    The grouping rule itself is pinned separately by the bin counts and the
    cut intervals asserted in the tests; this helper only makes the value
    comparison possible.
    """
    groups = []
    for chrom, start, end, _ in output_intervals:
        groups.append([i for i, (c, s, e, _) in enumerate(input_intervals)
                       if c == chrom and s >= start and e <= end])
    return groups


def test_running_window_3_is_bit_exact():
    """--runningWindow was never exercised by the previous test file."""
    sub, _ = nan_free_input(SMALL_50KB_H5)
    expected = dense_running_window(sub.toarray().astype(np.int64), 3)

    new = run_merge(SMALL_50KB_H5, 3, True, '.h5')
    assert new.matrix.shape == (2655, 2655)
    assert new.matrix.nnz == 273674
    nt.assert_array_equal(new.matrix.toarray(), expected)


def test_running_window_5_is_bit_exact():
    sub, _ = nan_free_input(SMALL_50KB_H5)
    expected = dense_running_window(sub.toarray().astype(np.int64), 5)

    new = run_merge(SMALL_50KB_H5, 5, True, '.h5')
    assert new.matrix.shape == (2655, 2655)
    assert new.matrix.nnz == 616057
    nt.assert_array_equal(new.matrix.toarray(), expected)


def test_running_window_keeps_the_resolution_and_the_intervals():
    """The running window does not change the bins, only the values.

    It also does not stop at chromosome borders: the window is applied to raw
    bin indices, so the last bins of one chromosome pick up counts from the
    first bins of the next. Pinned, not fixed.
    """
    _, intervals = nan_free_input(SMALL_50KB_H5)
    new = run_merge(SMALL_50KB_H5, 3, True, '.h5')
    nt.assert_equal(new.cut_intervals, intervals)

    boundaries = hm.hiCMatrix(SMALL_50KB_H5).chrBinBoundaries
    assert list(boundaries) == ['chr2RHet', 'chr3RHet', 'chr2LHet', 'chr4',
                               'chr3L', 'chr2L', 'chrU', 'chrX', 'chrXHet',
                               'chr2R', 'chr3R', 'chr3LHet']
    # A cell that spans the chr2RHet / chr3RHet border of the nan free matrix
    # is non zero although the two bins are on different chromosomes.
    sub, sub_intervals = nan_free_input(SMALL_50KB_H5)
    border = next(i for i in range(1, len(sub_intervals))
                  if sub_intervals[i][0] != sub_intervals[i - 1][0])
    assert sub[border - 1, border] == 0
    assert new.matrix[border - 1, border] != 0


def test_running_window_of_one_bin_only_drops_the_nan_bins():
    """num_bins == 1 returns early, so only remove_nans_if_needed runs."""
    sub, intervals = nan_free_input(SMALL_50KB_H5)
    new = run_merge(SMALL_50KB_H5, 1, True, '.h5')
    assert new.matrix.shape == (2655, 2655)
    nt.assert_array_equal(new.matrix.toarray(), sub.toarray())
    nt.assert_equal(new.cut_intervals, intervals)


def test_running_window_rejects_an_even_number_of_bins():
    """running_window_merge asserts num_bins % 2 == 1."""
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    os.unlink(outfile.name)
    with pytest.raises(AssertionError, match="num_bins has to be an odd number"):
        hicMergeMatrixBins.main(['--matrix', SMALL_50KB_H5,
                                 '--outFileName', outfile.name,
                                 '--numBins', '4', '--runningWindow'])
    assert not os.path.exists(outfile.name)


def test_running_window_cool_output_is_bit_exact():
    """The cool writer is a different code path from the h5 writer.

    small_test_matrix_50kb_res.cool holds 3,383 bins and 664 nan bins, so the
    running window runs on a 2,719 bin matrix here.
    """
    sub, _ = nan_free_input(SMALL_50KB_COOL)
    expected = dense_running_window(sub.toarray().astype(np.int64), 3)

    new = run_merge(SMALL_50KB_COOL, 3, True, '.cool')
    assert new.matrix.shape == (2719, 2719)
    nt.assert_array_equal(new.matrix.toarray(), expected)


def test_merge_bins_4_is_bit_exact():
    sub, intervals = nan_free_input(SMALL_50KB_H5)
    new = run_merge(SMALL_50KB_H5, 4, False, '.h5')

    assert new.matrix.shape == (665, 665)
    assert new.matrix.nnz == 27109
    groups = groups_from_output(intervals, new.cut_intervals)
    # 4 of the 2,655 nan free bins are not part of any merged bin: they are the
    # trailing bins of chr2RHet, chr3RHet, chr4 and chrX, each of which is left
    # over as a group of one, below the numBins/2 threshold of
    # hicMergeMatrixBins.py:248. Their counts are dropped silently.
    assert sum(len(g) for g in groups) == 2651
    nt.assert_array_equal(new.matrix.toarray().astype(np.float64),
                          dense_block_reduce(sub, groups))


def test_merge_bins_drops_chromosomes_with_too_few_bins():
    """A chromosome with fewer than numBins/2 remaining bins is discarded.

    With --numBins 20 chr2LHet (8 bins) and chrXHet (5 bins) fall below the
    numBins/2 threshold at hicMergeMatrixBins.py:248 and disappear from the
    output entirely, without an error. Trailing part-groups on chr3RHet, chr4,
    chrX and chr2R are dropped for the same reason, 33 of the 2,655 nan free
    bins in total. Their counts go with them: the sum over the merged matrix
    drops from 59,360 at --numBins 4 to 52,061 at --numBins 20. Pinned, not
    fixed.
    """
    sub, intervals = nan_free_input(SMALL_50KB_H5)
    new = run_merge(SMALL_50KB_H5, 20, False, '.h5')

    assert list(new.chrBinBoundaries) == ['chr2RHet', 'chr3RHet', 'chr4',
                                          'chr3L', 'chr2L', 'chrU', 'chrX',
                                          'chr2R', 'chr3R', 'chr3LHet']
    assert 'chr2LHet' in hm.hiCMatrix(SMALL_50KB_H5).chrBinBoundaries
    assert new.matrix.shape == (133, 133)

    groups = groups_from_output(intervals, new.cut_intervals)
    assert sum(len(g) for g in groups) == 2622
    dense = new.matrix.toarray().astype(np.float64)
    nt.assert_array_equal(dense, dense_block_reduce(sub, groups))
    assert dense.sum() == 52061.0

    at_four = run_merge(SMALL_50KB_H5, 4, False, '.h5')
    assert at_four.matrix.toarray().sum() == 59360.0


def test_merge_bins_cut_intervals_span_the_removed_nan_bins():
    """Merged intervals are built from the nan free bin list.

    The second merged interval of chr2RHet at --numBins 4 is 250 kb wide, not
    200 kb, because a nan bin was removed from the middle of it. The coverage
    of a merged bin is the mean of the coverages of the bins it holds.
    """
    new = run_merge(SMALL_50KB_H5, 4, False, '.h5')
    assert new.cut_intervals[0] == ('chr2RHet', 0, 200000, 0.625)
    assert new.cut_intervals[1] == ('chr2RHet', 200000, 450000, 0.425)
    assert new.cut_intervals[-1] == ('chr3LHet', 2400000, 2550000,
                                     0.5666666666666667)


def test_merge_bins_cool_output_is_bit_exact():
    sub, intervals = nan_free_input(SMALL_50KB_COOL)
    new = run_merge(SMALL_50KB_COOL, 4, False, '.cool')

    groups = groups_from_output(intervals, new.cut_intervals)
    nt.assert_array_equal(new.matrix.toarray().astype(np.float64),
                          dense_block_reduce(sub, groups))


def test_merge_bins_5_reproduces_the_stored_reference_exactly():
    """The same case as test_correct_matrix, but every stored value is
    compared, including the csr index arrays, and the reference is also
    checked against an independent block reduction.
    """
    reference = hm.hiCMatrix(ROOT + "hicMergeMatrixBins/result.h5")
    new = run_merge(SMALL_H5, 5, False, '.h5')

    got = new.matrix.tocsr()
    got.sort_indices()
    expected = reference.matrix.tocsr()
    expected.sort_indices()
    assert got.shape == expected.shape
    nt.assert_array_equal(got.indptr, expected.indptr)
    nt.assert_array_equal(got.indices, expected.indices)
    nt.assert_array_equal(got.data, expected.data)
    nt.assert_equal(new.cut_intervals, reference.cut_intervals)

    sub, intervals = nan_free_input(SMALL_H5)
    groups = groups_from_output(intervals, new.cut_intervals)
    nt.assert_array_equal(new.matrix.toarray().astype(np.float64),
                          dense_block_reduce(sub, groups))


def test_dense_models_are_not_vacuous():
    """The two models must change the input, otherwise nothing is tested."""
    sub, _ = nan_free_input(SMALL_50KB_H5)
    dense = sub.toarray().astype(np.int64)
    windowed = dense_running_window(dense, 3)
    assert not np.array_equal(windowed, dense)
    assert np.count_nonzero(windowed) > np.count_nonzero(dense)
