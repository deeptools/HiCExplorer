import warnings
import pytest
import os
from tempfile import NamedTemporaryFile
from hicmatrix import HiCMatrix as hm
from hicexplorer import hicAdjustMatrix
import numpy.testing as np
from hicexplorer.test.test_compute_function import compute


warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)


ROOT = os.path.join(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))), "test_data/")
matrix = ROOT + 'small_test_matrix_50kb_res.h5'
outfile = NamedTemporaryFile(suffix='.h5', prefix='test_matrix', delete=True)
bed_file = ROOT + 'regions.bed'
bed_file_xfail = ROOT + 'regions_xfail.bed'


@pytest.mark.parametrize("matrix", [matrix])  # required
@pytest.mark.parametrize("outFileName", [outfile])  # required
@pytest.mark.parametrize("chromosomes", ['chrX', 'chr3R'])  # optional
@pytest.mark.parametrize("action", ['keep', 'remove', 'mask'])  # optional
@pytest.mark.parametrize("regions", [bed_file, None])  # optional
def test_trivial_run(matrix, outFileName, chromosomes, action, regions):
    """
        Test checks if all commandline args work in general.
    """
    args = "--matrix {} --outFileName {} --chromosomes {} --action {}".format(
        matrix,
        outFileName.name,
        chromosomes,
        action,
    ).split()

    if regions:
        args = "--matrix {} --outFileName {} --regions {} --action {}".format(
            matrix,
            outFileName.name,
            regions,
            action,
        ).split()

    # hicAdjustMatrix.main(args)
    compute(hicAdjustMatrix.main, args, 5)


@pytest.mark.parametrize("matrix", [matrix])  # required
@pytest.mark.parametrize("outFileName", [outfile])  # required
@pytest.mark.parametrize("chromosomes", ['chr10'])  # optional
@pytest.mark.parametrize("action", ['keep', 'remove', 'mask'])  # optional
@pytest.mark.xfail
def test_trivial_run_xfail(matrix, outFileName, chromosomes, action):
    """
        Test checks if all commandline args work in general.
    """
    args = "--matrix {} --outFileName {} --chromosomes {} --action {}".format(
        matrix,
        outFileName.name,
        chromosomes,
        action,
    ).split()

    # hicAdjustMatrix.main(args)
    compute(hicAdjustMatrix.main, args, 5)


@pytest.mark.parametrize("matrix", [matrix])  # required
@pytest.mark.parametrize("outFileName", [outfile])  # required
@pytest.mark.parametrize("chromosomes", ['chr10', 'chr11'])  # optional
@pytest.mark.parametrize("action", ['keep', 'remove', 'mask'])  # optional
@pytest.mark.xfail
def test_trivial_run_xfail_multichromosomes(matrix, outFileName, chromosomes, action):
    """
        Test checks if all commandline args work in general.
    """
    args = "--matrix {} --outFileName {} --chromosomes {} --action {}".format(
        matrix,
        outFileName.name,
        chromosomes,
        action,
    ).split()

    # hicAdjustMatrix.main(args)
    compute(hicAdjustMatrix.main, args, 5)


@pytest.mark.parametrize("matrix", [matrix])  # required
@pytest.mark.parametrize("outFileName", [outfile])  # required
@pytest.mark.parametrize("action", ['keep', 'remove', 'mask'])  # optional
@pytest.mark.parametrize("regions", [bed_file_xfail])  # optional
@pytest.mark.xfail
def test_trivial_run_xfail_regions(matrix, outFileName, action, regions):
    """
        Test checks if all commandline args work in general.
    """

    if regions:
        args = "--matrix {} --outFileName {} --regions {} --action {}".format(
            matrix,
            outFileName.name,
            regions,
            action,
        ).split()

    # hicAdjustMatrix.main(args)
    compute(hicAdjustMatrix.main, args, 5)


def test_keep():
    outfile = NamedTemporaryFile(
        suffix='.h5', prefix='test_matrix', delete=True)
    outfile.close()
    args = "--matrix {} --outFileName {} --regions {} --action {}".format(
        ROOT + 'small_test_matrix_50kb_res.h5',
        outfile.name,
        ROOT + 'hicAdjustMatrix/keep_region.bed',
        "keep").split()

    compute(hicAdjustMatrix.main, args, 5)
    test = hm.hiCMatrix(
        ROOT + "hicAdjustMatrix/small_test_matrix_50kb_res_keep.h5")
    new = hm.hiCMatrix(outfile.name)
    np.assert_almost_equal(test.matrix.data, new.matrix.data, decimal=5)
    np.assert_equal(test.cut_intervals, new.cut_intervals)

    os.unlink(outfile.name)


def test_remove():
    outfile = NamedTemporaryFile(
        suffix='.h5', prefix='test_matrix', delete=True)
    outfile.close()
    args = "--matrix {} --outFileName {} --regions {} --action {}".format(
        ROOT + 'small_test_matrix_50kb_res.h5',
        outfile.name,
        ROOT + 'hicAdjustMatrix/remove.bed',
        "remove").split()

    compute(hicAdjustMatrix.main, args, 5)
    # HiCExplorer 4 output because of the fix: --action remove deletes only the
    # bins of the regions. The reference of 3.7.6 had lost the empty (NaN) bins
    # as well, so it was regenerated with the fixed tool.
    test = hm.hiCMatrix(
        ROOT + "hicAdjustMatrix/small_test_matrix_50kb_res_remove.h5")
    new = hm.hiCMatrix(outfile.name)
    np.assert_almost_equal(test.matrix.data, new.matrix.data, decimal=5)
    np.assert_equal(test.cut_intervals, new.cut_intervals)

    os.unlink(outfile.name)


def test_remove_inter():
    outfile = NamedTemporaryFile(
        suffix='.cool', prefix='test_matrix', delete=True)
    outfile.close()
    args = "--matrix {} --outFileName {} --chromosomes 1 2 3 --action {} --interIntraHandling {} ".format(
        ROOT + 'hicAdjustMatrix/gm12878_1_2_3.cool',
        outfile.name,
        "keep", "inter").split()

    compute(hicAdjustMatrix.main, args, 5)
    test = hm.hiCMatrix(
        ROOT + "hicAdjustMatrix/inter-removed.cool")
    new = hm.hiCMatrix(outfile.name)
    np.assert_almost_equal(test.matrix.data, new.matrix.data, decimal=5)
    np.assert_equal(test.cut_intervals, new.cut_intervals)

    os.unlink(outfile.name)


def test_remove_intra():
    outfile = NamedTemporaryFile(
        suffix='.cool', prefix='test_matrix', delete=True)
    outfile.close()
    args = "--matrix {} --outFileName {} --chromosomes 1 2 3 --action {} --interIntraHandling {} ".format(
        ROOT + 'hicAdjustMatrix/gm12878_1_2_3.cool',
        outfile.name,
        "keep", "intra").split()

    compute(hicAdjustMatrix.main, args, 5)
    test = hm.hiCMatrix(
        ROOT + "hicAdjustMatrix/intra-removed.cool")
    new = hm.hiCMatrix(outfile.name)
    np.assert_almost_equal(test.matrix.data, new.matrix.data, decimal=5)
    np.assert_equal(test.cut_intervals, new.cut_intervals)

    os.unlink(outfile.name)


# ---------------------------------------------------------------------------
# Characterization tests added for the v4 C++ port.
#
# Of the 30 collected items above, 24 are "did not crash" runs and the four
# real ones compare only matrix.data at decimal=5. test_remove_inter and
# test_remove_intra pass the BED path positionally right after
# --chromosomes 1 2 3, and since --chromosomes takes nargs='+' the path is
# swallowed as a fourth chromosome name and only warned about, so --regions is
# not exercised there at all. --maskBadRegions was never run.
#
# The tests below assert every stored value at full precision against a dense
# model built from the input, cover .h5 and .cool output, and pin the two
# --maskBadRegions behaviours.
# ---------------------------------------------------------------------------
import numpy as numpy_lib  # noqa: E402

ADJ = ROOT + "hicAdjustMatrix/"
GM_COOL = ADJ + "gm12878_1_2_3.cool"
SMALL_50KB_H5 = ROOT + "small_test_matrix_50kb_res.h5"
SMALL_50KB_COOL = ROOT + "small_test_matrix_50kb_res.cool"
REGIONS_BED = ROOT + "regions.bed"
MASK_BED = ADJ + "mask.bed"


def run_adjust(matrix_path, extra_args, suffix):
    out = NamedTemporaryFile(suffix=suffix, prefix='test_matrix', delete=False)
    out.close()
    try:
        hicAdjustMatrix.main(['--matrix', matrix_path,
                              '--outFileName', out.name] + extra_args)
        return hm.hiCMatrix(out.name)
    finally:
        if os.path.exists(out.name):
            os.unlink(out.name)


def region_bin_ids(matrix_path, bed_path):
    """The bin ids --regions selects, resolved the way the tool resolves them."""
    hic = hm.hiCMatrix(matrix_path)
    chromosomes = list(hic.chrBinBoundaries)
    ids = []
    with open(bed_path) as bed:
        for line in bed:
            fields = line.strip().split('\t')
            if len(fields) < 3 or fields[0] not in chromosomes:
                continue
            bin_range = hic.getRegionBinRange(fields[0], int(fields[1]),
                                              int(fields[2]))
            if bin_range is not None:
                # end is inclusive
                ids.extend(range(bin_range[0], bin_range[1] + 1))
    return ids


def chromosome_bin_ids(matrix_path, chromosomes):
    hic = hm.hiCMatrix(matrix_path)
    ids = []
    for chromosome in chromosomes:
        start, end = hic.chrBinBoundaries[chromosome]
        ids.extend(range(start, end))
    return ids


def dense_of(matrix_path):
    hic = hm.hiCMatrix(matrix_path)
    return hic.matrix.toarray().astype(float), hic


def test_regions_keep_h5_is_bit_exact():
    """--regions with --action keep, on its own and not behind --chromosomes.

    regions.bed selects chrX:10,000,000-15,000,000 and chrX:18,000,000-22,000,000,
    182 bins of small_test_matrix_50kb_res.h5. reorderBins keeps them in the
    order the BED file lists them.
    """
    dense, original = dense_of(SMALL_50KB_H5)
    selected = region_bin_ids(SMALL_50KB_H5, REGIONS_BED)
    assert len(selected) == 182

    new = run_adjust(SMALL_50KB_H5, ['--regions', REGIONS_BED,
                                     '--action', 'keep'], '.h5')
    assert new.matrix.shape == (182, 182)
    assert new.matrix.nnz == 995
    np.assert_array_equal(new.matrix.toarray().astype(float),
                          dense[numpy_lib.ix_(selected, selected)])
    np.assert_equal(new.cut_intervals,
                    [original.cut_intervals[i] for i in selected])


def test_regions_mask_h5_is_bit_exact():
    """--action mask keeps the matrix size and blanks the selected bins.

    The masked bins are the union of the 182 selected bins and the 139 nan bins
    the file already carries, 318 after the overlap, and the written matrix is
    upcast to float64 by the mask/restore round trip.
    """
    dense, original = dense_of(SMALL_50KB_H5)
    union = sorted(set(region_bin_ids(SMALL_50KB_H5, REGIONS_BED))
                   | set(original.nan_bins))
    assert len(union) == 318

    expected = dense.copy()
    expected[union, :] = 0
    expected[:, union] = 0

    new = run_adjust(SMALL_50KB_H5, ['--regions', REGIONS_BED,
                                     '--action', 'mask'], '.h5')
    assert new.matrix.shape == (2794, 2794)
    assert new.matrix.nnz == 43083
    assert new.matrix.dtype == numpy_lib.float64
    np.assert_array_equal(new.matrix.toarray().astype(float), expected)
    np.assert_array_equal(sorted(new.nan_bins), union)
    np.assert_equal(new.cut_intervals, original.cut_intervals)


def test_regions_remove_h5_is_bit_exact():
    """--action remove deletes exactly the bins of the regions.

    Fixed in HiCExplorer 4. hicAdjustMatrix.py:157-159 of 3.7.6 cleared
    orig_bin_ids, orig_cut_intervals and nan_bins after masking, so every
    empty (NaN) bin of the matrix went with the regions (2,794 - 318 = 2,476
    bins). Now only the 182 bins of the regions are removed (2,612 bins); the
    empty bins that remain stay in the matrix and keep their place in nan_bins,
    as they do for --chromosomes --action remove.
    """
    dense, original = dense_of(SMALL_50KB_H5)
    removed = set(region_bin_ids(SMALL_50KB_H5, REGIONS_BED))
    survivors = [i for i in range(dense.shape[0]) if i not in removed]
    position = {old: new for new, old in enumerate(survivors)}
    expected_nan = sorted(position[i] for i in original.nan_bins if i in position)

    new = run_adjust(SMALL_50KB_H5, ['--regions', REGIONS_BED,
                                     '--action', 'remove'], '.h5')
    assert new.matrix.shape == (2794 - len(removed), 2794 - len(removed))
    assert len(removed) == 182
    np.assert_array_equal(new.matrix.toarray().astype(float),
                          dense[numpy_lib.ix_(survivors, survivors)])
    np.assert_array_equal(sorted(new.nan_bins), expected_nan)
    np.assert_equal(new.cut_intervals,
                    [original.cut_intervals[i] for i in survivors])


def test_regions_mask_cool_is_bit_exact():
    """cool output takes a different writer than h5.

    The values match the same dense model as the h5 case, but the nan bin list
    that comes back does not: a cool file carries no explicit nan bin list, so
    hicmatrix rebuilds it on load as every all-zero row. That is 848 bins here,
    5 more than the 843 the tool masked, because masking emptied 5 further
    rows. The h5 path keeps the exact list instead. Pinned for both formats.
    """
    dense, original = dense_of(SMALL_50KB_COOL)
    union = sorted(set(region_bin_ids(SMALL_50KB_COOL, REGIONS_BED))
                   | set(original.nan_bins))
    assert len(union) == 843
    expected = dense.copy()
    expected[union, :] = 0
    expected[:, union] = 0

    new = run_adjust(SMALL_50KB_COOL, ['--regions', REGIONS_BED,
                                       '--action', 'mask'], '.cool')
    np.assert_array_equal(new.matrix.toarray().astype(float), expected)
    all_zero_rows = numpy_lib.flatnonzero(expected.sum(axis=0) == 0)
    assert len(new.nan_bins) == 848
    np.assert_array_equal(sorted(new.nan_bins), all_zero_rows)
    assert set(union).issubset(set(new.nan_bins))


def test_regions_mask_with_a_second_bed_file():
    """mask.bed selects two chrX intervals, one of them 2 Mb wide."""
    dense, original = dense_of(SMALL_50KB_H5)
    union = sorted(set(region_bin_ids(SMALL_50KB_H5, MASK_BED))
                   | set(original.nan_bins))
    expected = dense.copy()
    expected[union, :] = 0
    expected[:, union] = 0

    new = run_adjust(SMALL_50KB_H5, ['--regions', MASK_BED,
                                     '--action', 'mask'], '.h5')
    np.assert_array_equal(new.matrix.toarray().astype(float), expected)
    np.assert_array_equal(sorted(new.nan_bins), union)


def test_chromosomes_keep_reorders_the_matrix():
    """--chromosomes keeps the chromosomes in the order they are given.

    gm12878_1_2_3.cool holds chromosomes 1, 2 and 3; asking for '3 1' produces
    a matrix whose first block is chromosome 3.
    """
    dense, original = dense_of(GM_COOL)
    selected = chromosome_bin_ids(GM_COOL, ['3', '1'])

    new = run_adjust(GM_COOL, ['--chromosomes', '3', '1',
                               '--action', 'keep'], '.cool')
    assert new.matrix.shape == (449, 449)
    assert list(new.chrBinBoundaries) == ['3', '1']
    np.assert_array_equal(new.matrix.toarray().astype(float),
                          dense[numpy_lib.ix_(selected, selected)])
    np.assert_equal(new.cut_intervals,
                    [original.cut_intervals[i] for i in selected])


def test_chromosomes_remove_is_bit_exact():
    dense, _ = dense_of(GM_COOL)
    survivors = chromosome_bin_ids(GM_COOL, ['1', '3'])

    new = run_adjust(GM_COOL, ['--chromosomes', '2',
                               '--action', 'remove'], '.cool')
    assert new.matrix.shape == (449, 449)
    assert list(new.chrBinBoundaries) == ['1', '3']
    np.assert_array_equal(new.matrix.toarray().astype(float),
                          dense[numpy_lib.ix_(survivors, survivors)])


def test_chromosomes_mask_is_bit_exact():
    """maskChromosomes blanks chromosome 2 but keeps the 693 bins.

    The 244 bins of chromosome 2 plus the 25 nan bins the file carries, minus
    the 2 that overlap, give the 267 nan bins of the output.
    """
    dense, _ = dense_of(GM_COOL)
    masked = chromosome_bin_ids(GM_COOL, ['2'])
    expected = dense.copy()
    expected[masked, :] = 0
    expected[:, masked] = 0

    new = run_adjust(GM_COOL, ['--chromosomes', '2',
                               '--action', 'mask'], '.cool')
    assert new.matrix.shape == (693, 693)
    assert len(new.nan_bins) == 267
    assert new.matrix.dtype == numpy_lib.float64
    np.assert_array_equal(new.matrix.toarray().astype(float), expected)


def test_single_chromosome_keep_on_cool_takes_the_cooler_fast_path():
    """One chromosome, --action keep and a cool input load only that chromosome.

    hicAdjustMatrix.py:73-76 passes pChrnameList to hicmatrix instead of
    loading the whole file, so the dtype stays int32 rather than being upcast.
    The values must still equal the corresponding block of the full matrix.
    """
    dense, _ = dense_of(GM_COOL)
    selected = chromosome_bin_ids(GM_COOL, ['1'])

    new = run_adjust(GM_COOL, ['--chromosomes', '1',
                               '--action', 'keep'], '.cool')
    assert new.matrix.shape == (250, 250)
    assert new.matrix.dtype == numpy_lib.int32
    assert len(new.nan_bins) == 20
    np.assert_array_equal(new.matrix.toarray().astype(float),
                          dense[numpy_lib.ix_(selected, selected)])


def test_maskBadRegions_matches_regions_mask():
    """--maskBadRegions masks the bins of its BED file.

    Fixed in HiCExplorer 4. hicAdjustMatrix.py:161-165 of 3.7.6 never opened
    the file: an h5 matrix was written back unchanged and a cool matrix raised
    a TypeError. The option now does what its help says, on both formats, and
    gives the matrix --regions --action mask gives for the same BED file.
    """
    for source, suffix in [(SMALL_50KB_H5, '.h5'), (SMALL_50KB_COOL, '.cool')]:
        dense, original = dense_of(source)
        masked = (set(region_bin_ids(source, MASK_BED))
                  | set(original.nan_bins))
        expected = dense.copy()
        for i in masked:
            expected[i, :] = 0
            expected[:, i] = 0
        assert len(masked) > len(original.nan_bins)

        new = run_adjust(source, ['--maskBadRegions', MASK_BED], suffix)
        reference = run_adjust(source, ['--regions', MASK_BED,
                                        '--action', 'mask'], suffix)
        np.assert_array_equal(new.matrix.toarray().astype(float), expected)
        np.assert_array_equal(new.matrix.toarray(), reference.matrix.toarray())
        np.assert_equal(new.cut_intervals, original.cut_intervals)
        np.assert_array_equal(sorted(new.nan_bins), sorted(reference.nan_bins))
        assert set(new.nan_bins) >= masked


def test_maskBadRegions_ignores_action():
    """The bins are masked whatever --action says (it describes --regions)."""
    plain = run_adjust(SMALL_50KB_H5, ['--maskBadRegions', MASK_BED], '.h5')
    removed = run_adjust(SMALL_50KB_H5, ['--maskBadRegions', MASK_BED,
                                         '--action', 'remove'], '.h5')
    np.assert_array_equal(plain.matrix.toarray(), removed.matrix.toarray())


def test_no_selection_is_an_error():
    """Without --chromosomes, --regions or --maskBadRegions the tool exits 1.

    Fixed in HiCExplorer 4. 3.7.6 logged one line, wrote nothing and exited 0,
    which a pipeline could not tell from success.
    """
    out = NamedTemporaryFile(suffix='.h5', prefix='test_matrix', delete=False)
    out.close()
    os.unlink(out.name)
    with pytest.raises(SystemExit) as excinfo:
        hicAdjustMatrix.main(['--matrix', SMALL_50KB_H5,
                              '--outFileName', out.name])
    assert excinfo.value.code == 1
    assert not os.path.exists(out.name)


def test_interIntraHandling_without_the_stray_bed_argument():
    """The same result as test_remove_inter / test_remove_intra, cleanly.

    Those two tests pass the BED path positionally after --chromosomes 1 2 3,
    where nargs='+' turns it into a fourth chromosome name that is only
    warned about. Dropping it must not change the output, and every stored
    value is compared here instead of matrix.data at decimal=5.
    """
    for handling, reference in [('inter', 'inter-removed.cool'),
                                ('intra', 'intra-removed.cool')]:
        new = run_adjust(GM_COOL, ['--chromosomes', '1', '2', '3',
                                   '--action', 'keep',
                                   '--interIntraHandling', handling], '.cool')
        expected = hm.hiCMatrix(ADJ + reference)
        got = new.matrix.tocsr()
        got.sort_indices()
        want = expected.matrix.tocsr()
        want.sort_indices()
        assert got.shape == want.shape
        np.assert_array_equal(got.indptr, want.indptr)
        np.assert_array_equal(got.indices, want.indices)
        np.assert_array_equal(got.data, want.data)
        np.assert_equal(new.cut_intervals, expected.cut_intervals)


def test_unknown_chromosome_is_an_error():
    """A chromosome name that is not in the matrix exits 1 and writes nothing.

    Fixed in HiCExplorer 4. 3.7.6 warned about the name and went on with the
    others, so a misspelt chromosome silently changed the selection (with
    --action remove it stayed in the output).
    """
    for names, action in [(['1', 'not_a_chromosome'], 'keep'),
                          (['not_a_chromosome'], 'remove'),
                          (['1', 'not_a_chromosome'], 'remove'),
                          (['2', 'also_missing', 'not_a_chromosome'], 'mask')]:
        out = NamedTemporaryFile(suffix='.cool', prefix='test_matrix',
                                 delete=False)
        out.close()
        os.unlink(out.name)
        with pytest.raises(SystemExit) as excinfo:
            hicAdjustMatrix.main(['--matrix', GM_COOL, '--outFileName', out.name,
                                  '--chromosomes'] + names + ['--action', action])
        assert excinfo.value.code == 1
        assert not os.path.exists(out.name)


def test_unknown_chromosome_in_a_bed_file_is_still_only_a_warning(tmp_path):
    """Regions of other chromosomes in a BED file are skipped, not rejected."""
    bed = tmp_path / 'regions.bed'
    bed.write_text('chrX\t10000000\t15000000\nchrNotThere\t0\t1000\n')
    new = run_adjust(SMALL_50KB_H5, ['--regions', str(bed),
                                     '--action', 'keep'], '.h5')
    selected = region_bin_ids(SMALL_50KB_H5, str(bed))
    assert new.matrix.shape == (len(selected), len(selected))


def test_remove_warning_names_the_chromosome_of_the_region(tmp_path, capfd):
    """The 'discounted regions' warning reports the region it is about.

    Fixed in HiCExplorer 4: 3.7.6 used the chromosome of the last BED line.
    """
    hic = hm.hiCMatrix(SMALL_50KB_H5)
    start, end = hic.chrBinBoundaries['chrX']
    middle = hic.cut_intervals[(start + end) // 2]
    bed = tmp_path / 'regions.bed'
    bed.write_text('chrX\t{}\t{}\nchr2L\t0\t100000\n'.format(
        middle[1], middle[2]))
    run_adjust(SMALL_50KB_H5, ['--regions', str(bed), '--action', 'remove'],
               '.h5')
    captured = capfd.readouterr().err
    assert 'chrX:' in captured and 'discounted regions' in captured
    assert 'chr2L:' not in captured
