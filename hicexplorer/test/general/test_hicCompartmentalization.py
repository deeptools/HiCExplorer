import warnings
import os
from matplotlib.testing.compare import compare_images
from matplotlib.testing.exceptions import ImageComparisonFailure
import pytest

from hicexplorer import hicCompartmentalization
from tempfile import NamedTemporaryFile

warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
tolerance = 60  # default matplotlib
ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")
from hicexplorer.test.test_compute_function import compute


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
def test_compartmentalization():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)
    outfile.close()

    args = " -m {} --pca {} -o {} --outliers 0.0 --quantile 30".format(ROOT + "hicPCA/obsexp_norm.h5",
                                                                       ROOT + "hicCompartmentalization/pca1.bedgraph",
                                                                       outfile.name).split()
    # hicCompartmentalization.main(args)
    compute(hicCompartmentalization.main, args, 5)
    test = ROOT + "hicCompartmentalization/compartmentalizationRatio.png"
    res = compare_images(test, outfile.name, tolerance)
    assert res is None, res

    os.unlink(outfile.name)


# ---------------------------------------------------------------------------
# Full precision characterization tests, added 2026-09-02 for the C++ port.
#
# The test above is the whole of the existing coverage and it asserts nothing
# that can fail: it is xfail(ImageComparisonFailure) against a checked-in PNG
# at matplotlib's default tolerance of 60, so the tool is effectively
# uncovered. Neither --outputMatrix nor --offset nor --outliers nor a
# --quantile other than 30 has ever been exercised.
#
# What is pinned below is the numeric output, which the image is only a
# rendering of. hicCompartmentalization writes two numeric artefacts:
#
#   <outputFileName>_dat   np.savetxt of the polarization ratios, written
#                          unconditionally, one row per input matrix
#   --outputMatrix         np.savez of the normalised sum per quantile
#
# and both are compared here at float64 precision, with a small table of
# literal values so that the arithmetic is pinned and not only the plumbing.
#
# Three properties of the reference are pinned deliberately rather than
# corrected, because a port has to reproduce them:
#
#  1. **The first ratio is infinite.** within_vs_between_compartments starts at
#     q = 1, where the between-compartment block is the single cell
#     [0, quantile-1] and the input has nothing there, so the ratio is a
#     division by zero. np.savetxt writes it as the literal `inf`.
#  2. **Every quantile pair is counted twice on the diagonal.**
#     count_interactions adds each submatrix sum to both [qi, qj] and [qj, qi],
#     so a block with qi == qj is added to the same cell twice. The bin count
#     is doubled in the same breath, so the normalised value is unchanged and
#     this is invisible in the output; a port that fixed it would still agree.
#  3. **The per chromosome loop does not use its loop variable.** The body of
#     `for chrom in chromosomes` in count_interactions ignores `chrom` (the
#     per chromosome slicing is commented out at :107-110), so the whole
#     genome-wide computation is repeated once per chromosome in the PCA
#     bedgraph and both accumulators are multiplied by the chromosome count.
#     The ratio is again unchanged, and the cost is linear in the number of
#     chromosomes. test_result_does_not_depend_on_the_number_of_chromosomes
#     pins that the *values* do not depend on it.

import numpy as np
import numpy.testing as nt
import pandas as pd

from hicmatrix import HiCMatrix as hm
from hicexplorer.hicCompartmentalization import count_interactions
from hicexplorer.hicCompartmentalization import within_vs_between_compartments

OBSEXP = ROOT + "hicPCA/obsexp_norm.h5"
PCA = ROOT + "hicCompartmentalization/pca1.bedgraph"

# Produced by the reference on this tree, 2026-09-02, with
#   -m hicPCA/obsexp_norm.h5 --pca hicCompartmentalization/pca1.bedgraph
#   --outliers 0.0 --quantile 30
# Enough of the vector to pin the arithmetic; the whole of it is compared
# against the recomputed pipeline in the test below.
EXPECTED_HEAD = [np.inf, 56.20433807953578, 4.517170735129827,
                 2.5919137980265394, 2.4669236379464374, 1.986223966366796]
EXPECTED_TAIL = [1.024983838203027, 1.0167862225595856, 1.0043537932115325]


def _run(args):
    hicCompartmentalization.main(args.split())


def _quantiled_pc1(pQuantiles=30, pOutliers=0.0):
    """hicCompartmentalization.main:181-199, written out.

    Note the dtype the bedgraph is read with: the pc1 column is float32, so
    the quantile boundaries are computed from float32 values widened back to
    float64 and are not the quantiles of the file's decimal text.
    """
    pc1 = pd.read_table(PCA, header=None, sep="\t",
                        dtype={0: "object", 1: "Int64", 2: "Int64", 3: "float32"})
    pc1 = pc1.rename(columns={0: "chr", 1: "start", 2: "end", 3: "pc1"})
    if pOutliers != 0:
        quantile = [pOutliers / 100, (100 - pOutliers) / 100]
        boundaries = np.nanquantile(pc1['pc1'].values.astype(float), quantile)
        quantiled_bins = np.linspace(boundaries[0], boundaries[1], pQuantiles)
    else:
        quantile = [j / (pQuantiles - 1) for j in range(0, pQuantiles)]
        quantiled_bins = np.nanquantile(pc1['pc1'].values.astype(float), quantile)
    pc1["quantile"] = np.searchsorted(quantiled_bins,
                                      pc1['pc1'].values.astype(float),
                                      side="right")
    return pc1


def _reference_ratios(pQuantiles=30, pOutliers=0.0, pOffset=None):
    """The whole pipeline, recomputed in the test at float64 precision."""
    pc1 = _quantiled_pc1(pQuantiles, pOutliers)
    obs_exp = hm.hiCMatrix(OBSEXP)
    pc1["bin_id"] = pc1.apply(
        lambda row: np.arange(
            obs_exp.getRegionBinRange(row['chr'], row['start'], row['end'] - 1)[0],
            obs_exp.getRegionBinRange(row['chr'], row['start'], row['end'] - 1)[1] + 1),
        axis=1)
    per_quantile = count_interactions(obs_exp, pc1, pQuantiles, pOffset)
    per_quantile = np.nan_to_num(per_quantile)
    return per_quantile, np.array(
        within_vs_between_compartments(per_quantile, pQuantiles))


def _read_dat(pPath):
    return np.loadtxt(pPath + '_dat')


def test_polarization_ratio_dat_matches_the_pipeline_at_full_precision():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)
    outfile.close()
    _run("-m {} --pca {} -o {} --outliers 0.0 --quantile 30".format(
        OBSEXP, PCA, outfile.name))

    actual = _read_dat(outfile.name)
    _, expected = _reference_ratios()

    assert actual.shape == (29,), \
        'quantile N produces N - 1 ratios, one per prefix length'
    nt.assert_array_equal(actual, expected)
    os.unlink(outfile.name)
    os.unlink(outfile.name + '_dat')


def test_polarization_ratio_matches_the_literal_reference_values():
    """A table of literal values, so that the test above is not comparing the
    module against itself alone."""
    outfile = NamedTemporaryFile(suffix='.png', delete=False)
    outfile.close()
    _run("-m {} --pca {} -o {} --outliers 0.0 --quantile 30".format(
        OBSEXP, PCA, outfile.name))

    actual = _read_dat(outfile.name)
    assert np.isinf(actual[0]), \
        'the first ratio divides by an empty between-compartment block'
    nt.assert_allclose(actual[1:6], EXPECTED_HEAD[1:], rtol=1e-12, atol=0)
    nt.assert_allclose(actual[-3:], EXPECTED_TAIL, rtol=1e-12, atol=0)
    # The ratio falls monotonically after the first few quantiles, which is
    # what "the two compartments separate" means here.
    assert actual[-1] < actual[5] < actual[2]
    os.unlink(outfile.name)
    os.unlink(outfile.name + '_dat')


def test_outputMatrix_holds_the_normalised_sum_per_quantile():
    """--outputMatrix has never been tested.

    np.savez is called with a *list* of arrays, so the file has one entry,
    `arr_0`, whose shape is (number of input matrices, quantiles, quantiles),
    not one entry per matrix.
    """
    outfile = NamedTemporaryFile(suffix='.png', delete=False)
    outfile.close()
    matrixfile = NamedTemporaryFile(suffix='.npz', delete=False)
    matrixfile.close()
    _run("-m {} --pca {} -o {} --outliers 0.0 --quantile 30 "
         "--outputMatrix {}".format(OBSEXP, PCA, outfile.name, matrixfile.name))

    stored = np.load(matrixfile.name)
    assert list(stored.keys()) == ['arr_0']
    actual = stored['arr_0']
    assert actual.shape == (1, 30, 30)
    assert actual.dtype == np.float64

    expected, _ = _reference_ratios()
    nt.assert_array_equal(actual[0], expected)
    # Quantile 0 is empty on this input: np.searchsorted with side='right'
    # never returns 0 for a value that equals the smallest boundary, so every
    # bin lands in 1..29 and the first row and column stay zero.
    assert np.all(actual[0][0, :] == 0.0)
    assert np.all(actual[0][:, 0] == 0.0)
    # The matrix is symmetric by construction: every sum is added to both
    # [qi, qj] and [qj, qi].
    nt.assert_array_equal(actual[0], actual[0].T)
    os.unlink(outfile.name)
    os.unlink(outfile.name + '_dat')
    os.unlink(matrixfile.name)


def test_offset_sets_the_named_diagonals_to_nan():
    """--offset has never been tested.

    It writes NaN into the given diagonals of the matrix *before* the
    quantile sums are taken, and count_interactions then drops every NaN, so
    the effect is to exclude those diagonals from both the sum and the bin
    count. Note the side effect a port has to reproduce: assigning NaN to a
    position that held no stored entry makes it a stored entry, so --offset
    changes the sparsity of the matrix as well as its values.
    """
    plain = NamedTemporaryFile(suffix='.png', delete=False)
    plain.close()
    offset0 = NamedTemporaryFile(suffix='.png', delete=False)
    offset0.close()
    offset01 = NamedTemporaryFile(suffix='.png', delete=False)
    offset01.close()

    _run("-m {} --pca {} -o {} --outliers 0.0 --quantile 30".format(
        OBSEXP, PCA, plain.name))
    _run("-m {} --pca {} -o {} --outliers 0.0 --quantile 30 --offset 0".format(
        OBSEXP, PCA, offset0.name))
    _run("-m {} --pca {} -o {} --outliers 0.0 --quantile 30 "
         "--offset 0 1".format(OBSEXP, PCA, offset01.name))

    without = _read_dat(plain.name)
    with_main = _read_dat(offset0.name)
    with_both = _read_dat(offset01.name)

    nt.assert_array_equal(with_main, _reference_ratios(pOffset=[0])[1])
    nt.assert_array_equal(with_both, _reference_ratios(pOffset=[0, 1])[1])
    # Excluding the main diagonal removes the strongest within-compartment
    # signal there is, so the ratios have to move.
    assert not np.allclose(without[1:], with_main[1:])
    assert not np.allclose(with_main[1:], with_both[1:])

    for handle in (plain, offset0, offset01):
        os.unlink(handle.name)
        os.unlink(handle.name + '_dat')


def test_quantile_count_changes_the_binning_not_only_the_length():
    """--quantile has only ever been passed as 30."""
    ten = NamedTemporaryFile(suffix='.png', delete=False)
    ten.close()
    _run("-m {} --pca {} -o {} --outliers 0.0 --quantile 10".format(
        OBSEXP, PCA, ten.name))

    actual = _read_dat(ten.name)
    assert actual.shape == (9,)
    nt.assert_array_equal(actual, _reference_ratios(pQuantiles=10)[1])
    os.unlink(ten.name)
    os.unlink(ten.name + '_dat')


def test_outliers_switches_the_boundaries_to_a_linear_range():
    """--outliers has only ever been passed as 0.0, which takes the other
    branch: a non-zero value replaces the empirical quantiles by
    np.linspace between the two trimmed boundaries, which is a different
    binning and not a trimmed version of the same one."""
    trimmed = NamedTemporaryFile(suffix='.png', delete=False)
    trimmed.close()
    _run("-m {} --pca {} -o {} --outliers 5.0 --quantile 30".format(
        OBSEXP, PCA, trimmed.name))

    actual = _read_dat(trimmed.name)
    nt.assert_array_equal(actual, _reference_ratios(pOutliers=5.0)[1])

    plain = NamedTemporaryFile(suffix='.png', delete=False)
    plain.close()
    _run("-m {} --pca {} -o {} --outliers 0.0 --quantile 30".format(
        OBSEXP, PCA, plain.name))
    assert not np.allclose(actual[1:], _read_dat(plain.name)[1:])

    for handle in (trimmed, plain):
        os.unlink(handle.name)
        os.unlink(handle.name + '_dat')


def test_two_matrices_produce_one_row_each():
    """-m takes nargs='+' and no test ever passed more than one matrix."""
    outfile = NamedTemporaryFile(suffix='.png', delete=False)
    outfile.close()
    _run("-m {} {} --pca {} -o {} --outliers 0.0 --quantile 30".format(
        OBSEXP, ROOT + "hicPCA/obsexp_norm.cool", PCA, outfile.name))

    actual = _read_dat(outfile.name)
    assert actual.shape == (2, 29), \
        'one row of ratios per input matrix'
    os.unlink(outfile.name)
    os.unlink(outfile.name + '_dat')


def _count_interactions_written_out(pQuantiles, pOffset=None):
    """count_interactions, written out in the test rather than imported.

    The tests above import count_interactions, so a mutation inside it moves
    both sides equally: measured, replacing the NaN filter at :123 by a no-op
    is invisible to every one of them. This is the independent form, and it is
    what makes the --offset assertions mean something.
    """
    pc1 = _quantiled_pc1(pQuantiles)
    obs_exp = hm.hiCMatrix(OBSEXP)
    bin_ids = [
        np.arange(
            obs_exp.getRegionBinRange(row['chr'], row['start'], row['end'] - 1)[0],
            obs_exp.getRegionBinRange(row['chr'], row['start'], row['end'] - 1)[1] + 1)
        for _, row in pc1.iterrows()]

    matrix = obs_exp.matrix.copy()
    if pOffset:
        for dist in pOffset:
            indices = np.arange(0, matrix.shape[0] - dist)
            matrix[indices, indices + dist] = np.nan
            matrix[indices + dist, indices] = np.nan

    quantiles = pc1['quantile'].values
    total = np.zeros((pQuantiles, pQuantiles))
    counts = np.zeros((pQuantiles, pQuantiles))
    # The loop variable is unused in the reference, so the whole computation
    # is repeated once per chromosome and both accumulators scale with it.
    for _ in pc1['chr'].unique():
        for qi in range(pQuantiles):
            if not (quantiles == qi).any():
                continue
            rows = np.concatenate([bin_ids[k]
                                   for k in np.flatnonzero(quantiles == qi)])
            for qj in range(pQuantiles):
                if not (quantiles == qj).any():
                    continue
                cols = np.concatenate([bin_ids[k]
                                       for k in np.flatnonzero(quantiles == qj)])
                block = np.asarray(matrix[np.ix_(rows, cols)].todense())
                finite = block[np.isfinite(block)]
                total[qi, qj] += finite.sum()
                total[qj, qi] += finite.sum()
                counts[qi, qj] += finite.size
                counts[qj, qi] += finite.size
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.nan_to_num(total / counts)


def test_normalised_sum_per_quantile_written_out_with_and_without_offset():
    """Six quantiles, the whole count independent of the module's own code.

    Also pins what --offset actually removes: with --offset 0 the diagonal is
    NaN and therefore excluded from both the sum and the bin count, so the
    within-quantile cells fall while the off-diagonal ones are untouched.
    """
    for offset, flag in ((None, ''), ([0], ' --offset 0')):
        outfile = NamedTemporaryFile(suffix='.png', delete=False)
        outfile.close()
        matrixfile = NamedTemporaryFile(suffix='.npz', delete=False)
        matrixfile.close()
        _run("-m {} --pca {} -o {} --outliers 0.0 --quantile 6 "
             "--outputMatrix {}{}".format(OBSEXP, PCA, outfile.name,
                                          matrixfile.name, flag))
        actual = np.load(matrixfile.name)['arr_0'][0]
        nt.assert_allclose(actual, _count_interactions_written_out(6, offset),
                           rtol=1e-12, atol=0)
        os.unlink(outfile.name)
        os.unlink(outfile.name + '_dat')
        os.unlink(matrixfile.name)

    plain = NamedTemporaryFile(suffix='.npz', delete=False)
    plain.close()
    dropped = NamedTemporaryFile(suffix='.npz', delete=False)
    dropped.close()
    without = _count_interactions_written_out(6, None)
    with_offset = _count_interactions_written_out(6, [0])
    diagonal_falls = all(with_offset[q, q] < without[q, q]
                         for q in range(6) if without[q, q] != 0.0)
    assert diagonal_falls, \
        '--offset 0 must remove the main diagonal from the within-quantile cells'
    os.unlink(plain.name)
    os.unlink(dropped.name)


# ---------------------------------------------------------------------------
# Further characterization, added 2026-09-13 with the C++ port.
#
# The comment at the top of the 2026-09-02 block names
# test_result_does_not_depend_on_the_number_of_chromosomes, but no such test
# was committed; it is added here. The other tests pin behaviour a port has to
# reproduce and that nothing above constrains:
#
#  4. The bin holding the largest pc1 value is dropped. The boundaries run from
#     the minimum to the maximum, and np.searchsorted(side='right') puts a value
#     equal to the last boundary at index quantile, one past the last quantile,
#     so that bin is in no block. This is the mirror image of quantile 0 being
#     empty.
#  5. A NaN pc1 is dropped the same way. :199 intends to move NaN rows out of
#     the way, but `pc1["pc1"] == np.nan` is never true and the chained
#     indexing assigns into a copy, so the line does nothing; searchsorted
#     sorts NaN last and returns quantile as well.
#  6. np.savez appends '.npz' to an --outputMatrix name that lacks it.
#  7. --quantile 1 raises ZeroDivisionError at :192 unless --outliers is given;
#     --quantile 0 succeeds and writes one empty line.

import shutil
from tempfile import mkdtemp


def _write_bedgraph(pPath, pRows):
    with open(pPath, 'w') as handle:
        for row in pRows:
            handle.write("\t".join(row) + "\n")


def _bedgraph_rows():
    with open(PCA) as handle:
        return [line.rstrip("\n").split("\t") for line in handle if line.strip()]


def _run_matrix(pPca, pQuantiles, pFolder, pName):
    plot = os.path.join(pFolder, pName + '.png')
    matrix = os.path.join(pFolder, pName + '.npz')
    _run("-m {} --pca {} -o {} --quantile {} --outputMatrix {}".format(
        OBSEXP, pPca, plot, pQuantiles, matrix))
    return np.load(matrix)['arr_0'][0]


def test_the_bin_with_the_largest_pc1_is_in_no_quantile():
    folder = mkdtemp(prefix="testCompartments_")
    actual = _run_matrix(PCA, 6, folder, 'plain')

    # The written out count with the reference binning agrees ...
    nt.assert_allclose(actual, _count_interactions_written_out(6), rtol=1e-12, atol=0)

    # ... because the maximum really is outside the range, and exactly one row is.
    pc1 = _quantiled_pc1(6)
    assert (pc1['quantile'] == 6).sum() == 1
    assert pc1.loc[pc1['quantile'] == 6, 'pc1'].iloc[0] == pc1['pc1'].max()
    assert (pc1['quantile'] == 0).sum() == 0

    # Putting that bin into the top quantile would have changed the output.
    rows = _bedgraph_rows()
    top = int(np.argmax([float(row[3]) for row in rows]))
    capped = [list(row) for row in rows]
    second = sorted(float(row[3]) for row in rows)[-2]
    capped[top][3] = repr(second)
    capped_path = os.path.join(folder, 'capped.bedgraph')
    _write_bedgraph(capped_path, capped)
    moved = _run_matrix(capped_path, 6, folder, 'capped')
    assert not np.allclose(actual, moved)
    shutil.rmtree(folder)


def test_a_nan_pc1_row_is_dropped_as_if_it_were_absent():
    folder = mkdtemp(prefix="testCompartments_")
    rows = _bedgraph_rows()
    victim = 100
    assert rows[victim][0] == 'chrX'

    with_nan = [list(row) for row in rows]
    with_nan[victim][3] = 'nan'
    nan_path = os.path.join(folder, 'nan.bedgraph')
    _write_bedgraph(nan_path, with_nan)

    without = [row for index, row in enumerate(rows) if index != victim]
    without_path = os.path.join(folder, 'without.bedgraph')
    _write_bedgraph(without_path, without)

    nan_result = _run_matrix(nan_path, 10, folder, 'nan')
    without_result = _run_matrix(without_path, 10, folder, 'without')
    nt.assert_array_equal(nan_result, without_result)
    assert not np.array_equal(nan_result, _run_matrix(PCA, 10, folder, 'plain'))
    shutil.rmtree(folder)


def test_outputMatrix_without_the_suffix_is_written_with_it():
    folder = mkdtemp(prefix="testCompartments_")
    plot = os.path.join(folder, 'ratio.png')
    matrix = os.path.join(folder, 'matrix')
    _run("-m {} --pca {} -o {} --quantile 6 --outputMatrix {}".format(
        OBSEXP, PCA, plot, matrix))
    assert not os.path.exists(matrix)
    assert np.load(matrix + '.npz')['arr_0'].shape == (1, 6, 6)
    shutil.rmtree(folder)


def test_quantile_one_fails_and_quantile_zero_writes_an_empty_line():
    folder = mkdtemp(prefix="testCompartments_")
    plot = os.path.join(folder, 'ratio.png')
    with pytest.raises(ZeroDivisionError):
        _run("-m {} --pca {} -o {} --quantile 1".format(OBSEXP, PCA, plot))
    assert not os.path.exists(plot + '_dat')

    _run("-m {} --pca {} -o {} --quantile 1 --outliers 2.5".format(OBSEXP, PCA, plot))
    with open(plot + '_dat') as handle:
        assert handle.read() == "\n"
    os.unlink(plot + '_dat')

    _run("-m {} --pca {} -o {} --quantile 0".format(OBSEXP, PCA, plot))
    with open(plot + '_dat') as handle:
        assert handle.read() == "\n"
    shutil.rmtree(folder)


def test_result_does_not_depend_on_the_number_of_chromosomes():
    """count_interactions repeats the whole count once per chromosome of the
    pca file, because its chromosome loop ignores the loop variable. Both
    accumulators scale together, so relabelling every row to a single
    chromosome, which leaves the bins and the quantiles alone, gives the same
    normalised sums to rounding.

    Only to rounding: s + s + s is not always 3 * s in floating point, so the
    last bits do depend on the chromosome count. Measured on a six bin fixture
    with non dyadic values, one against three chromosomes differ in 5 of 16
    cells by one ulp. On this input the two happen to agree bit for bit, which
    is not a property to pin, hence rtol=1e-12. A port has to repeat the
    additions to reproduce the bits (cpp/tools/compartmentalization_impl.cpp
    replays them)."""
    pc1 = _quantiled_pc1(10)
    obs_exp = hm.hiCMatrix(OBSEXP)
    pc1["bin_id"] = pc1.apply(
        lambda row: np.arange(
            obs_exp.getRegionBinRange(row['chr'], row['start'], row['end'] - 1)[0],
            obs_exp.getRegionBinRange(row['chr'], row['start'], row['end'] - 1)[1] + 1),
        axis=1)
    assert len(pc1['chr'].unique()) == 2

    two = np.nan_to_num(count_interactions(obs_exp, pc1, 10, None))
    relabelled = pc1.copy()
    relabelled['chr'] = 'chrX'
    one = np.nan_to_num(count_interactions(hm.hiCMatrix(OBSEXP), relabelled, 10, None))
    nt.assert_allclose(two, one, rtol=1e-12, atol=0)
    assert np.count_nonzero(two) > 0
