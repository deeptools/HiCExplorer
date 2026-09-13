
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import os
from tempfile import NamedTemporaryFile
from hicexplorer import hicAverageRegions
import numpy as np
import numpy.testing as nt

from scipy.sparse import load_npz
from hicexplorer.test.test_compute_function import compute

import logging
log = logging.getLogger(__name__)

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")

# The exact float64 values of every case, as (nnz, [(value, count), ...]).
#
# The old assertions in this file were assert_almost_equal(decimal=0) on the
# data array alone: agreement to the nearest integer, on values that all lie in
# [0.16, 1.17], so every one of them would have passed against an array of
# zeros. These literals pin the same runs at full float64 precision without
# adding a new binary master, and they are what makes the file able to fail.
#
# Two properties they encode that no master in the repository does:
#
#  * the output dtype is **float64**. The checked-in .npz masters are float32
#    because the scipy of the day returned a float32 matrix from
#    `lil(float32) /= ndarray(float64)`; scipy 1.14 returns float64.
#  * the values are `sum * (1 / n)` and not `sum / n`. scipy's
#    `_spbase._divide` takes the reciprocal of the dense operand and multiplies
#    (scipy/sparse/_base.py), and the two differ in the last bit: with six
#    regions and a summed count of five, the division gives
#    0.8333333333333334 and the multiplication 0.8333333333333333, which is the
#    value below.
EXPECTED = {
    'multi_start': (191, [(0.16666666666666666, 156), (0.3333333333333333, 28),
                          (0.5, 7)]),
    'multi_center': (252, [(0.16666666666666666, 204), (0.3333333333333333, 35),
                           (0.5, 8), (0.6666666666666666, 2),
                           (0.8333333333333333, 3)]),
    'multi_end': (233, [(0.16666666666666666, 192), (0.3333333333333333, 32),
                        (0.5, 7), (0.6666666666666666, 2)]),
    'multi_bins_start': (1814, [(0.16666666666666666, 1583),
                                (0.3333333333333333, 178), (0.5, 39),
                                (0.6666666666666666, 9), (0.8333333333333333, 3),
                                (1.0, 2)]),
    'multi_bins_center': (2005, [(0.16666666666666666, 1772),
                                 (0.3333333333333333, 182), (0.5, 35),
                                 (0.6666666666666666, 10), (0.8333333333333333, 4),
                                 (1.0, 2)]),
    'multi_bins_end': (2043, [(0.16666666666666666, 1818),
                              (0.3333333333333333, 172), (0.5, 40),
                              (0.6666666666666666, 11), (1.1666666666666665, 2)]),
    'single': (191, [(0.16666666666666666, 156), (0.3333333333333333, 28),
                     (0.5, 7)]),
    'single_bins': (1814, [(0.16666666666666666, 1583), (0.3333333333333333, 178),
                           (0.5, 39), (0.6666666666666666, 9),
                           (0.8333333333333333, 3), (1.0, 2)]),
}


def assert_values_exact(pKey, pPath):
    """The stored values, bit for bit, against the literals above."""
    observed = load_npz(pPath)
    nnz, histogram = EXPECTED[pKey]
    nt.assert_equal(observed.dtype, np.float64)
    nt.assert_equal(observed.nnz, nnz)
    values, counts = np.unique(observed.data, return_counts=True)
    nt.assert_equal(list(values), [value for value, _ in histogram])
    nt.assert_equal(list(counts), [count for _, count in histogram])


def assert_matches_master(pMaster, pPath):
    """The sparsity pattern exactly, and the values to float32 precision.

    The masters were written when the tool still produced float32, so they can
    only be compared to float32 accuracy; that is still seven orders of
    magnitude tighter than the decimal=0 this file used, and it is the
    assertion that catches a value landing in the wrong cell, which a
    comparison of the data array alone never could.
    """
    expected = load_npz(pMaster)
    observed = load_npz(pPath)
    nt.assert_equal(expected.shape, observed.shape)
    nt.assert_equal(expected.nnz, observed.nnz)
    nt.assert_equal(expected.indices, observed.indices)
    nt.assert_equal(expected.indptr, observed.indptr)
    nt.assert_allclose(observed.data, expected.data, rtol=1e-6, atol=0)


def run(pOutFile, pArguments, pRegions='hicAverageRegions/regions_multi.bed'):
    args = "--matrix {} --regions {} -o {} {}".format(
        ROOT + 'small_test_matrix.cool', ROOT + pRegions, pOutFile, pArguments).split()
    compute(hicAverageRegions.main, args, 5)


def test_average_regions_start():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    run(outfile.name, "--range 100000 100000 -cb start")

    assert_matches_master(ROOT + 'hicAverageRegions/regions_multi_start.npz',
                          outfile.name)
    assert_values_exact('multi_start', outfile.name)

    os.remove(outfile.name)


def test_average_regions_range_in_bins_start():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    run(outfile.name, "--rangeInBins 100 100 -cb start")

    assert_matches_master(
        ROOT + 'hicAverageRegions/regions_multi_range_in_bins_start.npz', outfile.name)
    assert_values_exact('multi_bins_start', outfile.name)

    os.remove(outfile.name)


def test_average_regions_center():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    run(outfile.name, "--range 100000 100000 -cb center")

    assert_matches_master(ROOT + 'hicAverageRegions/regions_multi_center.npz',
                          outfile.name)
    assert_values_exact('multi_center', outfile.name)

    os.remove(outfile.name)


def test_average_regions_range_in_bins_center():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    run(outfile.name, "--rangeInBins 100 100 -cb center")

    assert_matches_master(
        ROOT + 'hicAverageRegions/regions_multi_range_in_bins_center.npz', outfile.name)
    assert_values_exact('multi_bins_center', outfile.name)

    os.remove(outfile.name)


def test_average_regions_end():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    run(outfile.name, "--range 100000 100000 -cb end")

    assert_matches_master(ROOT + 'hicAverageRegions/regions_multi_end.npz',
                          outfile.name)
    assert_values_exact('multi_end', outfile.name)

    os.remove(outfile.name)


def test_average_regions_range_in_bins_end():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    run(outfile.name, "--rangeInBins 100 100 -cb end")

    assert_matches_master(
        ROOT + 'hicAverageRegions/regions_multi_range_in_bins_end.npz', outfile.name)
    assert_values_exact('multi_bins_end', outfile.name)

    os.remove(outfile.name)


def test_average_regions_single():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    run(outfile.name, "--range 100000 100000", 'hicAverageRegions/regions.bed')

    assert_matches_master(ROOT + 'hicAverageRegions/result_range_100000.npz',
                          outfile.name)
    assert_values_exact('single', outfile.name)

    os.remove(outfile.name)


def test_average_regions_range_in_bins_single():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    run(outfile.name, "--rangeInBins 100 100", 'hicAverageRegions/regions.bed')

    assert_matches_master(ROOT + 'hicAverageRegions/result_rangeInBins_100.npz',
                          outfile.name)
    assert_values_exact('single_bins', outfile.name)

    os.remove(outfile.name)


def test_two_column_region_file_is_a_point():
    """regions.bed has two columns and regions_multi.bed three.

    With two columns the viewpoint is (chrom, start, start)
    (hicAverageRegions.py:147), so start, centre and end all fall on the same
    coordinate and the three --coordinatesToBinMapping values must give the
    same file. Never asserted, and it is the only thing that makes the
    two-column branch distinguishable from the three-column one.
    """
    results = []
    for mapping in ('start', 'center', 'end'):
        outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region',
                                     delete=False)
        run(outfile.name, "--range 100000 100000 -cb " + mapping,
            'hicAverageRegions/regions.bed')
        results.append(load_npz(outfile.name))
        os.remove(outfile.name)

    for other in results[1:]:
        nt.assert_equal(results[0].indices, other.indices)
        nt.assert_equal(results[0].indptr, other.indptr)
        nt.assert_equal(results[0].data, other.data)


def test_consider_strand_direction_changes_nothing():
    """--considerStrandDirection is a no-op, and this pins that it is.

    The option was never exercised. Its help text promises that "the contacts
    of a reverse strand region are inverted e.g. [1,2,3] becomes [3,2,1]", and
    the implementation transposes the submatrix (hicAverageRegions.py:199).
    But the submatrix is `matrix[start:end, start:end]`, a diagonal block of a
    Hi-C matrix, and a Hi-C matrix is symmetric: hicmatrix fills the lower
    triangle on load. A symmetric block is its own transpose, so the reverse
    strand branch adds exactly what the forward branch would have added.

    So the option cannot change any output, and a transpose is not a reversal
    in any case. Both are worth reporting upstream. The test asserts the
    behaviour as it stands: with three of the six regions on the minus strand,
    the result is identical to the run without the option, byte for byte in
    every array.
    """
    with_option = NamedTemporaryFile(suffix='.npz', prefix='average_region',
                                     delete=False)
    without_option = NamedTemporaryFile(suffix='.npz', prefix='average_region',
                                        delete=False)
    run(with_option.name, "--range 100000 100000 -cb start --considerStrandDirection",
        'hicAverageRegions/regions_multi_strand.bed')
    run(without_option.name, "--range 100000 100000 -cb start",
        'hicAverageRegions/regions_multi_strand.bed')

    left = load_npz(with_option.name)
    right = load_npz(without_option.name)
    nt.assert_equal(left.shape, right.shape)
    nt.assert_equal(left.indices, right.indices)
    nt.assert_equal(left.indptr, right.indptr)
    nt.assert_equal(left.data, right.data)

    # And the strand file gives the same numbers as the plain three-column one,
    # since the extra columns carry no information the tool acts on.
    assert_values_exact('multi_start', with_option.name)

    os.remove(with_option.name)
    os.remove(without_option.name)


def test_consider_strand_direction_needs_six_columns():
    """--considerStrandDirection on a three-column file exits 1.

    hicAverageRegions.py:150-152 logs an error and calls exit(1) when the
    strand column is missing. Never tested, and it is the tool's only
    validation.
    """
    import pytest

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    with pytest.raises(SystemExit) as raised:
        run(outfile.name, "--range 100000 100000 --considerStrandDirection")
    nt.assert_equal(raised.value.code, 1)
    os.remove(outfile.name)


def test_output_is_a_csr_npz_with_the_expected_arrays():
    """The container, not just the numbers.

    scipy.sparse.save_npz writes five entries in a fixed order and with fixed
    dtypes, and hicPlotAverageRegions reads the file back with load_npz. A port
    that wrote a valid but differently laid out archive would still load, so
    the layout is pinned here: the entry names and order, the format tag, and
    the dtype of each array.
    """
    import zipfile

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    run(outfile.name, "--range 100000 100000 -cb start")

    with zipfile.ZipFile(outfile.name) as archive:
        nt.assert_equal(archive.namelist(),
                        ['indices.npy', 'indptr.npy', 'format.npy', 'shape.npy',
                         'data.npy'])
        arrays = {name: np.load(archive.open(name), allow_pickle=False)
                  for name in archive.namelist()}

    nt.assert_equal(arrays['format.npy'].tobytes(), b'csr')
    nt.assert_equal(arrays['indices.npy'].dtype, np.int32)
    nt.assert_equal(arrays['indptr.npy'].dtype, np.int32)
    nt.assert_equal(arrays['shape.npy'].dtype, np.int64)
    nt.assert_equal(list(arrays['shape.npy']), [40, 40])
    nt.assert_equal(arrays['data.npy'].dtype, np.float64)

    os.remove(outfile.name)


def test_out_of_range_region_is_skipped_not_clipped():
    """A window that does not fit is dropped, and does not enter the count.

    hicAverageRegions.py:190 compares the shape of the submatrix against the
    accumulator's and skips the region with a warning. Never asserted, and it
    is what makes the `start_out`/`end_out` bookkeeping above it dead code.

    regions_multi_out_of_range.bed is regions_multi.bed with one extra region
    at the very start of chr2L, whose 100 kb upstream window is clipped to the
    chromosome boundary and is therefore 20 bins wide instead of 40. The result
    must be exactly the six-region result, not a seven-region average and not a
    partially accumulated corner.
    """
    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    run(outfile.name, "--range 100000 100000 -cb start",
        'hicAverageRegions/regions_multi_out_of_range.bed')

    assert_values_exact('multi_start', outfile.name)
    assert_matches_master(ROOT + 'hicAverageRegions/regions_multi_start.npz',
                          outfile.name)

    os.remove(outfile.name)


def test_float32_accumulator_on_a_float64_matrix():
    """The running sum is rounded to float32 after every region.

    `summed_matrix` is a lil_matrix of dtype float32 (hicAverageRegions.py:171)
    and `lil[a:b, a:b] += csr` writes the promoted sum back into float32
    storage, so the accumulation loses precision once per region. On the small
    integer counts of small_test_matrix that is invisible, which is why no
    existing test constrains it: the sums are at most 48 and every float32 in
    that range is exact.

    Li_et_al_2015.h5 is the case where it is visible. Its values are float64
    between 0.17 and 1914.01, and six of them summed and rounded to float32
    differ from the float64 sum in the low bits: the first value below,
    377.5341796875, is a float32 to the last digit and could not come out of a
    float64 accumulator.

    The file is 11,104 restriction fragment bins of chromosome X, so --range is
    not usable (the tool's own warning) and --rangeInBins is used instead.
    """
    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    args = ("--matrix {} --regions {} -o {} --rangeInBins 20 20 -cb start").format(
        ROOT + 'Li_et_al_2015.h5',
        ROOT + 'hicAverageRegions/regions_Li_chrX.bed', outfile.name).split()
    compute(hicAverageRegions.main, args, 5)

    observed = load_npz(outfile.name)
    nt.assert_equal(observed.shape, (40, 40))
    nt.assert_equal(observed.dtype, np.float64)
    nt.assert_equal(observed.nnz, 1600)
    nt.assert_equal(list(observed.data[:5]),
                    [377.5341796875, 131.0463155110677, 72.43905639648438,
                     54.21789042154948, 85.37013753255208])
    nt.assert_equal(list(observed.data[-5:]),
                    [41.99386088053385, 50.94233194986979, 62.0654551188151,
                     172.66215006510416, 467.67118326822913])
    nt.assert_equal(float(observed.data.min()), 3.104411443074544)
    nt.assert_equal(float(observed.data.max()), 540.2255045572916)

    os.remove(outfile.name)
