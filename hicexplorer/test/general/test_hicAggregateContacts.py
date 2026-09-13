import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from matplotlib.testing.exceptions import ImageComparisonFailure
from hicexplorer.test.test_compute_function import compute
import os.path
from tempfile import NamedTemporaryFile
import hicexplorer.hicAggregateContacts
import pytest
from psutil import virtual_memory
mem = virtual_memory()
memory = mem.total / 2 ** 30

# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 4
HIGH_MEMORY = 120

REMOVE_OUTPUT = True
# DIFF = 60


ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")
tolerance = 50  # default matplotlib pixed difference tolerance


def are_files_equal(file1, file2, delta=2, skip=0):
    equal = True
    if delta:
        mismatches = 0
    with open(file1) as textfile1, open(file2) as textfile2:
        for i, (x, y) in enumerate(zip(textfile1, textfile2)):
            # if x.startswith('File'):
            #     continue
            if i < skip:
                continue
            if x != y:
                if delta:
                    mismatches += 1
                    if mismatches > delta:
                        equal = False
                        break
                else:
                    equal = False
                    break
    return equal


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_intra_perChr():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)
    outfile_aggregate_matrix = NamedTemporaryFile(suffix='.tab', prefix='hicaggregate_test_', delete=False)
    args = "--matrix {root}/Li_et_al_2015.h5 --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --disable_bbox_tight --dpi 100 "\
           "--mode intra-chr --perChr --outFilePrefixMatrix {out_mat}".\
        format(root=ROOT, out_agg=outfile_aggregate_plots.name, out_mat=outfile_aggregate_matrix.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate_intra_perChr.png'
    test_matrix_agg = ROOT + 'hicAggregateContacts/master_aggregate_intra_perChr.tab'

    # hicexplorer.hicAggregateContacts.main(args.split())
    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)
    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res
    assert are_files_equal(test_matrix_agg, outfile_aggregate_matrix.name)
    os.remove(outfile_aggregate_plots.name)
    os.remove(outfile_aggregate_matrix.name)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_intra():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)
    outfile_aggregate_matrix = NamedTemporaryFile(suffix='.tab', prefix='hicaggregate_test_', delete=False)
    args = "--matrix {root}/Li_et_al_2015.h5 --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --disable_bbox_tight --dpi 100 "\
           "--mode intra-chr --outFilePrefixMatrix {out_mat}".\
        format(root=ROOT, out_agg=outfile_aggregate_plots.name, out_mat=outfile_aggregate_matrix.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate_intra.png'
    test_matrix_agg = ROOT + 'hicAggregateContacts/master_aggregate_intra.tab'
    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)
    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res
    assert are_files_equal(test_matrix_agg, outfile_aggregate_matrix.name)
    os.remove(outfile_aggregate_plots.name)
    os.remove(outfile_aggregate_matrix.name)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_inter():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)
    outfile_aggregate_matrix = NamedTemporaryFile(suffix='.tab', prefix='hicaggregate_test_', delete=False)
    args = "--matrix {root}/small_test_matrix.h5 --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --disable_bbox_tight --dpi 100 "\
           "--mode inter-chr  --operationType mean --outFilePrefixMatrix {out_mat}".\
           format(root=ROOT, out_agg=outfile_aggregate_plots.name, out_mat=outfile_aggregate_matrix.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate_inter.png'
    test_matrix_agg = ROOT + 'hicAggregateContacts/master_aggregate_inter.tab'
    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)
    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res
    assert are_files_equal(test_matrix_agg, outfile_aggregate_matrix.name)
    os.remove(outfile_aggregate_plots.name)
    os.remove(outfile_aggregate_matrix.name)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_all():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)
    outfile_aggregate_matrix = NamedTemporaryFile(suffix='.tab', prefix='hicaggregate_test_', delete=False)
    args = "--matrix {root}/Li_et_al_2015.h5 --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --disable_bbox_tight --dpi 100 "\
           "--mode all --outFilePrefixMatrix {out_mat} --keep_outlier".\
           format(root=ROOT, out_agg=outfile_aggregate_plots.name, out_mat=outfile_aggregate_matrix.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate_all.png'
    test_matrix_agg = ROOT + 'hicAggregateContacts/master_aggregate_all.tab'
    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)
    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res
    assert are_files_equal(test_matrix_agg, outfile_aggregate_matrix.name)
    os.remove(outfile_aggregate_plots.name)
    os.remove(outfile_aggregate_matrix.name)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_chromosome_not_given():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)

    args = "--matrix {root}/Li_et_al_2015.h5 --BED {root}/hicAggregateContacts/test_regions_region_not_given.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --disable_bbox_tight --dpi 100 "\
           "--mode intra-chr --perChr ".format(root=ROOT, out_agg=outfile_aggregate_plots.name)

    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_cooler():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)
    outfile_aggregate_matrix = NamedTemporaryFile(suffix='.tab', prefix='hicaggregate_test_', delete=False)

    args = "--matrix {root}/Li_et_al_2015.cool --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --disable_bbox_tight --dpi 100 "\
           "--mode intra-chr --perChr --outFilePrefixMatrix {out_mat}".\
           format(root=ROOT, out_agg=outfile_aggregate_plots.name, out_mat=outfile_aggregate_matrix.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate_intra_perChr.png'  # noqa: F841
    test_matrix_agg = ROOT + 'hicAggregateContacts/master_aggregate_intra_perChr.tab'
    # hicexplorer.hicAggregateContacts.main(args.split())
    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)

    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res
    assert are_files_equal(test_matrix_agg, outfile_aggregate_matrix.name)

    os.remove(outfile_aggregate_plots.name)
    os.remove(outfile_aggregate_matrix.name)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_clustering():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)
    outfile_heatmaps = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_heatmap_', delete=False)
    outfile_aggregate_matrix = NamedTemporaryFile(suffix='.tab', prefix='hicaggregate_test_', delete=False)

    args = "--matrix {root}/Li_et_al_2015.h5 --BED {root}/hicAggregateContacts/test_regions.bed " \
        "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --hclust 4 " \
        "--diagnosticHeatmapFile {out_heat} --howToCluster diagonal  --disable_bbox_tight --dpi 100 " \
        "--BED2 {root}/hicAggregateContacts/test_regions.bed --mode intra-chr --perChr  --outFilePrefixMatrix {out_mat}".\
        format(root=ROOT, out_agg=outfile_aggregate_plots.name,
               out_heat=outfile_heatmaps.name, out_mat=outfile_aggregate_matrix.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate_hclust4.png'
    test_image_heatmap = ROOT + 'hicAggregateContacts/master_heatmap.png'
    test_matrix1_agg = ROOT + 'hicAggregateContacts/master_aggregate_clust1_hclust4.tab'
    test_matrix2_agg = ROOT + 'hicAggregateContacts/master_aggregate_clust2_hclust4.tab'
    test_matrix3_agg = ROOT + 'hicAggregateContacts/master_aggregate_clust3_hclust4.tab'
    test_matrix4_agg = ROOT + 'hicAggregateContacts/master_aggregate_clust4_hclust4.tab'

    # hicexplorer.hicAggregateContacts.main(args.split())
    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)

    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res

    res = compare_images(test_image_heatmap, outfile_heatmaps.name, tolerance)
    assert res is None, res

    assert are_files_equal(test_matrix1_agg, outfile_aggregate_matrix.name + "_X_cluster_1.tab")
    assert are_files_equal(test_matrix2_agg, outfile_aggregate_matrix.name + "_X_cluster_2.tab")
    assert are_files_equal(test_matrix3_agg, outfile_aggregate_matrix.name + "_X_cluster_3.tab")
    assert are_files_equal(test_matrix4_agg, outfile_aggregate_matrix.name + "_X_cluster_4.tab")

    os.unlink(outfile_aggregate_plots.name)
    os.unlink(outfile_heatmaps.name)
    for i in range(1, 5):
        file = outfile_aggregate_matrix.name + "_X_cluster_" + str(i) + ".tab"
        os.unlink(file)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_clustering_cool():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)
    outfile_heatmaps = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_heatmap_', delete=False)
    outfile_aggregate_matrix = NamedTemporaryFile(suffix='.tab', prefix='hicaggregate_test_', delete=False)

    args = "--matrix {root}/Li_et_al_2015.cool --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --hclust 4 " \
           "--diagnosticHeatmapFile {out_heat} --howToCluster diagonal  --disable_bbox_tight --dpi 100 " \
           "--BED2 {root}/hicAggregateContacts/test_regions.bed  --mode intra-chr --perChr --outFilePrefixMatrix {out_mat}".\
           format(root=ROOT, out_agg=outfile_aggregate_plots.name,
                  out_heat=outfile_heatmaps.name, out_mat=outfile_aggregate_matrix.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate_hclust4.png'
    test_image_heatmap = ROOT + 'hicAggregateContacts/master_heatmap.png'
    test_matrix1_agg = ROOT + 'hicAggregateContacts/master_aggregate_clust1_hclust4.tab'
    test_matrix2_agg = ROOT + 'hicAggregateContacts/master_aggregate_clust2_hclust4.tab'
    test_matrix3_agg = ROOT + 'hicAggregateContacts/master_aggregate_clust3_hclust4.tab'
    test_matrix4_agg = ROOT + 'hicAggregateContacts/master_aggregate_clust4_hclust4.tab'

    # hicexplorer.hicAggregateContacts.main(args.split())
    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)

    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res

    res = compare_images(test_image_heatmap, outfile_heatmaps.name, tolerance)
    assert res is None, res

    assert are_files_equal(test_matrix1_agg, outfile_aggregate_matrix.name + "_X_cluster_1.tab")
    assert are_files_equal(test_matrix2_agg, outfile_aggregate_matrix.name + "_X_cluster_2.tab")
    assert are_files_equal(test_matrix3_agg, outfile_aggregate_matrix.name + "_X_cluster_3.tab")
    assert are_files_equal(test_matrix4_agg, outfile_aggregate_matrix.name + "_X_cluster_4.tab")

    os.unlink(outfile_aggregate_plots.name)
    os.unlink(outfile_heatmaps.name)
    for i in range(1, 5):
        file = outfile_aggregate_matrix.name + "_X_cluster_" + str(i) + ".tab"
        os.unlink(file)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_3d():

    outfile_aggregate_3d = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_3d', delete=False)

    args = "--matrix {root}/Li_et_al_2015.h5 --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --hclust 2 --dpi 100 " \
           "--plotType 3d --disable_bbox_tight " \
           "--BED2 {root}/hicAggregateContacts/test_regions.bed  --mode intra-chr --perChr".\
           format(root=ROOT, out_agg=outfile_aggregate_3d.name)

    test_image_agg_3d = ROOT + 'hicAggregateContacts/master_aggregate_3d.png'

    # hicexplorer.hicAggregateContacts.main(args.split())
    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)

    res = compare_images(test_image_agg_3d, outfile_aggregate_3d.name, tolerance)
    assert res is None, res

    os.remove(outfile_aggregate_3d.name)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_3d_cooler():

    outfile_aggregate_3d = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_3d', delete=False)

    args = "--matrix {root}/Li_et_al_2015.cool --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --hclust 2 --dpi 100 " \
           "--plotType 3d --disable_bbox_tight  --mode intra-chr --perChr " \
           "--BED2 {root}/hicAggregateContacts/test_regions.bed".format(root=ROOT, out_agg=outfile_aggregate_3d.name)

    test_image_agg_3d = ROOT + 'hicAggregateContacts/master_aggregate_3d.png'

    # hicexplorer.hicAggregateContacts.main(args.split())
    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)

    res = compare_images(test_image_agg_3d, outfile_aggregate_3d.name, tolerance)
    assert res is None, res

    os.remove(outfile_aggregate_3d.name)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_row_wise_intra_perChr():

    outfile_aggregate_row_wise = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_row_wise', delete=False)

    args = "--matrix {root}/Li_et_al_2015.h5 --BED {root}/hicAggregateContacts/bed1_row-wise.bed " \
           "--BED2 {root}/hicAggregateContacts/bed2_row-wise.bed "\
           "--outFileName {out_agg} --numberOfBins 30 --row_wise "\
           "--range 50000:6000000 --dpi 100 --mode intra-chr --perChr --keep_outlier".\
           format(root=ROOT, out_agg=outfile_aggregate_row_wise.name)

    test_image_agg_row_wise = ROOT + 'hicAggregateContacts/master_aggregate_row_wise_intra_perChr.png'

    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)
    res = compare_images(test_image_agg_row_wise, outfile_aggregate_row_wise.name, tolerance)
    assert res is None, res

    os.remove(outfile_aggregate_row_wise.name)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_row_wise_inter():

    outfile_aggregate_row_wise = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_row_wise', delete=False)

    args = "--matrix {root}/small_test_matrix_50kb_res.h5 --BED {root}/hicAggregateContacts/bed1_row-wise.bed " \
        "--BED2 {root}/hicAggregateContacts/bed2_row-wise.bed "\
        "--outFileName {out_agg} --numberOfBins 30 --row_wise "\
        "--dpi 100 --mode inter-chr ".\
        format(root=ROOT, out_agg=outfile_aggregate_row_wise.name)

    test_image_agg_row_wise = ROOT + 'hicAggregateContacts/master_aggregate_row_wise_inter.png'

    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)
    res = compare_images(test_image_agg_row_wise, outfile_aggregate_row_wise.name, tolerance)
    assert res is None, res

    os.remove(outfile_aggregate_row_wise.name)


# ---------------------------------------------------------------------------
# Characterization tests of the numeric outputs, added 2026-09-13 for the C++
# port.
#
# Everything above asserts on images only, and every one of those tests is
# xfail(ImageComparisonFailure). The .tab comparisons they also make are
# vacuous with two exceptions: the tool writes `<prefix>_<chrom>.tab`, while
# the tests compare the empty NamedTemporaryFile `<prefix>` itself, and
# are_files_equal zips the two files, so an empty file always compares equal.
# Only the two clustering tests read the real `_X_cluster_N.tab` names. The
# trivial_runs tests assert nothing, and test_hicAggregateContacts_trivial_runs
# _three.py builds a fixed argument string, so its 72 parametrised cases are
# the same run. --considerStrandDirection, --largeRegionsOperation,
# --outFileObsExp, --spectral and --max_deviation were never exercised.
#
# What is pinned below is what the tool writes besides its figures:
#
#   --outFilePrefixMatrix   np.savetxt('%0.5f') of the aggregate submatrix, one
#                           file per chromosome (or `genome`) and cluster
#   --outFileContactPairs   the contact positions and centre values, one file
#                           per chromosome (or `genome`) and cluster
#   --outFileObsExp         the obs/exp matrix, with --transform obs/exp
#
# Reproduced defects, pinned rather than corrected, because a port has to
# reproduce them:
#
#  1. --spectral is parsed but main() never reads it, so it is silently
#     ignored and no clustering happens (test_spectral_option_is_ignored).
#  2. With more than one cluster, the contact pair file of a cluster writes
#     the centre value of its cl_idx-th member but the coordinates of the
#     cl_idx-th submatrix of *all* clusters, because coords is indexed with the
#     position inside the cluster (plot_aggregated_contacts, `coords[cl_idx]`).
#     test_contact_pairs_of_a_cluster_take_coordinates_by_position pins it.
#  3. A '-'/'-' strand pair is transposed (`mat_to_append.T`), not flipped on
#     both axes as the comment says (test_strand_direction_flips_and_transposes).
#  4. Outlier removal only ever triggers with --howToCluster center: for the
#     full and diagonal vectors more than half of all values are zero on these
#     inputs, so the median absolute value, which get_outlier_indices uses as
#     the scale where a median absolute *deviation* was meant, is zero and the
#     function returns None.
#  5. In row-wise mode the coordinates are written as given, without the
#     reordering the non row-wise mode applies when the second bin is smaller.
#  6. With --perChr, a chromosome with fewer submatrices than clusters sets k
#     to 1 for itself and for every chromosome after it
#     (test_per_chromosome_k_falls_to_one_and_stays).
#  7. --chromosomes on a matrix with NaN bins undoes maskBins and enlarge_bins,
#     because keepOnlyTheseChr starts with restoreMaskedBins
#     (test_chromosomes_restore_the_masked_bins).

import glob

import numpy as np

AGG = ROOT + "hicAggregateContacts/"
INTRA = "--matrix {root}Li_et_al_2015.h5 --mode intra-chr --range 50000:900000 "


def _aggregate(pTmpPath, pArgs, pName="run"):
    """Runs the tool into a fresh directory and returns that directory."""
    out = pTmpPath / pName
    out.mkdir()
    args = ("--outFileName {out}/aggregate.png --numberOfBins 30 "
            "--outFilePrefixMatrix {out}/m --outFileContactPairs {out}/p "
            + pArgs).format(root=ROOT, agg=AGG, out=out).split()
    hicexplorer.hicAggregateContacts.main(args)
    return out


def _tables(pDirectory):
    return sorted(os.path.basename(path) for path in glob.glob(str(pDirectory) + "/*.tab"))


def _matrix(pDirectory, pName):
    return np.loadtxt(os.path.join(str(pDirectory), pName))


def _pairs(pDirectory, pName):
    with open(os.path.join(str(pDirectory), pName)) as handle:
        return handle.read().splitlines()


def _assert_matrix(pMatrix, pSum, pCenter, pRightOfCenter, pNonZero):
    assert pMatrix.shape == (31, 31)
    assert pMatrix.sum() == pytest.approx(pSum, rel=1e-9)
    assert pMatrix[15, 15] == pytest.approx(pCenter, rel=1e-9, abs=1e-12)
    assert pMatrix[15, 16] == pytest.approx(pRightOfCenter, rel=1e-9, abs=1e-12)
    assert int((pMatrix != 0).sum()) == pNonZero


def test_aggregate_intra_per_chromosome_values(tmp_path):
    out = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed --perChr")
    assert _tables(out) == ["m_X.tab", "p_X.tab"]
    _assert_matrix(_matrix(out, "m_X.tab"), 421.56695, 1.09234, 1.07726, 447)
    pairs = _pairs(out, "p_X.tab")
    assert len(pairs) == 279
    assert pairs[0] == "X\t1268823\t1268824\tX\t1356647\t1356648\t24.85811715426542"
    assert pairs[1] == "X\t1588015\t1588016\tX\t1679121\t1679122\t18.38930207884974"
    # sorted by the centre value, highest first, one line per submatrix
    values = [float(line.split("\t")[6]) for line in pairs]
    assert values == sorted(values, reverse=True)
    assert len(set(tuple(line.split("\t")[:6]) for line in pairs)) == 279


def test_operation_types_sum_and_mean(tmp_path):
    summed = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed --operationType sum", "sum")
    mean = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed --operationType mean", "mean")
    matrix_sum = _matrix(summed, "m_genome.tab")
    matrix_mean = _matrix(mean, "m_genome.tab")
    _assert_matrix(matrix_sum, 462055.15265, 569.53375, 548.22945, 961)
    _assert_matrix(matrix_mean, 1656.1117, 2.04134, 1.96498, 961)
    # 279 submatrices, and both files are rounded to five decimals
    assert np.allclose(matrix_sum / 279, matrix_mean, rtol=0, atol=1e-5)


def test_kmeans_partition_and_values(tmp_path):
    out = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed --kmeans 4")
    assert _tables(out) == ["m_genome_cluster_{}.tab".format(i) for i in range(1, 5)] + \
        ["p_genome_cluster_{}.tab".format(i) for i in range(1, 5)]
    sizes = [len(_pairs(out, "p_genome_cluster_{}.tab".format(i))) for i in range(1, 5)]
    assert sizes == [241, 6, 1, 31]
    _assert_matrix(_matrix(out, "m_genome_cluster_1.tab"), 135.74855, 0.9317, 0.78936, 173)
    _assert_matrix(_matrix(out, "m_genome_cluster_2.tab"), 9345.84237, 5.92278, 3.40061, 957)
    _assert_matrix(_matrix(out, "m_genome_cluster_3.tab"), 13787.64353, 8.52615, 16.62724, 916)
    _assert_matrix(_matrix(out, "m_genome_cluster_4.tab"), 4349.26925, 6.55486, 5.32221, 961)


def test_contact_pairs_of_a_cluster_take_coordinates_by_position(tmp_path):
    """Pinned defect 2: coords[cl_idx] instead of coords[cluster_indices[cl_idx]].

    The coordinates written for a cluster of n members are those of the first
    n submatrices overall, so every smaller cluster repeats coordinates of the
    largest one, while the centre values are the cluster's own.
    """
    out = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed --kmeans 4")
    coordinates = []
    for i in range(1, 5):
        coordinates.append(set(tuple(line.split("\t")[:6])
                               for line in _pairs(out, "p_genome_cluster_{}.tab".format(i))))
    for smaller in coordinates[1:]:
        assert smaller <= coordinates[0]
    assert _pairs(out, "p_genome_cluster_4.tab")[0] == \
        "X\t366894\t366895\tX\t548547\t548548\t24.85811715426542"


def test_hclust_per_chromosome_diagonal(tmp_path):
    out = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed --perChr --hclust 4 "
                     "--howToCluster diagonal")
    sizes = [len(_pairs(out, "p_X_cluster_{}.tab".format(i))) for i in range(1, 5)]
    assert sizes == [25, 81, 169, 4]
    _assert_matrix(_matrix(out, "m_X_cluster_1.tab"), 4811.53988, 7.63756, 4.81992, 961)
    _assert_matrix(_matrix(out, "m_X_cluster_2.tab"), 1372.97318, 1.97454, 1.5188, 923)
    _assert_matrix(_matrix(out, "m_X_cluster_3.tab"), 2.54214, 0.0, 0.0, 4)
    _assert_matrix(_matrix(out, "m_X_cluster_4.tab"), 10910.09729, 13.45773, 17.02556, 960)


def test_spectral_option_is_ignored(tmp_path):
    """Pinned defect 1: --spectral never reaches the clustering."""
    plain = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed", "plain")
    spectral = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed --spectral 3", "spectral")
    assert _tables(spectral) == ["m_genome.tab", "p_genome.tab"]
    for name in _tables(plain):
        with open(os.path.join(str(plain), name)) as a, open(os.path.join(str(spectral), name)) as b:
            assert a.read() == b.read()


@pytest.mark.parametrize("pOption, pPairs", [("--max_deviation 1", 229),
                                              ("", 251),
                                              ("--max_deviation 5", 272),
                                              ("--keep_outlier", 279)])
def test_max_deviation_controls_outlier_removal(tmp_path, pOption, pPairs):
    """Outliers are removed from the 279 submatrices before clustering.

    Only the centre layout has a non-zero scale on this input (defect 4), so
    that is the one that shows the threshold.
    """
    out = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed --kmeans 2 "
                     "--howToCluster center " + pOption)
    total = sum(len(_pairs(out, name)) for name in _tables(out) if name.startswith("p_"))
    assert total == pPairs


def _write_bed_with_strand(pPath, pStrand):
    with open(AGG + "test_regions.bed") as source, open(str(pPath), "w") as target:
        for line in source:
            fields = line.split()
            if fields[0] == "X":
                target.write("\t".join(fields[:5] + [pStrand]) + "\n")
    return str(pPath)


def test_consider_strand_direction_values(tmp_path):
    out = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions_strand.bed "
                     "--considerStrandDirection")
    _assert_matrix(_matrix(out, "m_genome.tab"), 324.85484, 1.09234, 1.14894, 407)
    assert len(_pairs(out, "p_genome.tab")) == 279


def test_strand_direction_flips_and_transposes(tmp_path):
    """'+'/'-' flips the columns, '-'/'+' the rows, '-'/'-' transposes.

    Every submatrix is oriented the same way here, and an elementwise median
    commutes with any fixed permutation of the cells, so the aggregate of the
    oriented submatrices is the oriented aggregate. The '-'/'-' case is a
    transpose, not a rotation (defect 3).
    """
    plus = _write_bed_with_strand(tmp_path / "plus.bed", "+")
    minus = _write_bed_with_strand(tmp_path / "minus.bed", "-")
    unstranded = _matrix(_aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed", "plain"),
                         "m_genome.tab")

    def stranded(pFirst, pSecond, pName):
        return _matrix(_aggregate(tmp_path, INTRA + "--BED {} --BED2 {} --considerStrandDirection"
                                  .format(pFirst, pSecond), pName), "m_genome.tab")

    assert np.array_equal(stranded(plus, plus, "pp"), unstranded)
    assert np.array_equal(stranded(plus, minus, "pm"), np.fliplr(unstranded))
    assert np.array_equal(stranded(minus, plus, "mp"), np.flipud(unstranded))
    assert np.array_equal(stranded(minus, minus, "mm"), unstranded.T)
    assert not np.array_equal(unstranded.T, np.flipud(np.fliplr(unstranded)))


def test_consider_strand_direction_requires_six_columns(tmp_path):
    with pytest.raises(SystemExit) as error:
        _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed --considerStrandDirection")
    assert error.value.code == 1


@pytest.mark.parametrize("pOperation, pPairs, pSum, pCenter, pRight, pNonZero",
                         [("first", 279, 421.56695, 1.09234, 1.07726, 447),
                          ("last", 280, 375.81971, 0.68732, 0.0, 407),
                          ("center", 278, 416.90271, 0.0, 0.0, 442)])
def test_large_regions_operation(tmp_path, pOperation, pPairs, pSum, pCenter, pRight, pNonZero):
    """Regions of 12 kb span several bins of Li_et_al_2015 (about 2 kb each)."""
    out = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions_large.bed "
                     "--largeRegionsOperation " + pOperation)
    _assert_matrix(_matrix(out, "m_genome.tab"), pSum, pCenter, pRight, pNonZero)
    assert len(_pairs(out, "p_genome.tab")) == pPairs


def test_out_file_obs_exp(tmp_path):
    import tables
    import cooler
    out = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed --transform obs/exp "
                     "--outFileObsExp {out}/obs_exp.h5")
    _assert_matrix(_matrix(out, "m_genome.tab"), 164.25439, 0.63621, 0.56861, 447)
    # the NaN bins are masked out of the matrix before the transform, so the
    # written matrix has 10,249 of the 11,104 bins, upper triangle only
    with tables.open_file(str(out / "obs_exp.h5")) as handle:
        assert list(handle.root.matrix.shape.read()) == [10249, 10249]
        data = handle.root.matrix.data.read()
        assert len(data) == 1661678
        assert float(np.sum(data)) == pytest.approx(3418201.0000000014, rel=1e-12)
    # Li_et_al_2015 holds no contacts beyond the band that --range 50000:900000
    # keeps, so the depth passed to the transform, int(max) * 2.5, only shows
    # with a narrower range: 60000 * 2.5 * 1.5 / 1843 bp keeps 122 diagonals.
    narrow = _aggregate(tmp_path, "--matrix {root}Li_et_al_2015.h5 --mode intra-chr "
                        "--range 20000:60000 --BED {agg}test_regions.bed --transform obs/exp "
                        "--outFileObsExp {out}/obs_exp.h5", "narrow")
    _assert_matrix(_matrix(narrow, "m_genome.tab"), 761.66451, 0.9869, 0.13327, 919)
    with tables.open_file(str(narrow / "obs_exp.h5")) as handle:
        data = handle.root.matrix.data.read()
        assert len(data) == 1034113
        assert float(np.sum(data)) == pytest.approx(2971253.999999999, rel=1e-12)
    # the file type follows the output name, not the input
    out_cool = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed --transform obs/exp "
                          "--outFileObsExp {out}/obs_exp.cool", "cool")
    assert cooler.fileops.is_cooler(str(out_cool / "obs_exp.cool"))
    assert cooler.Cooler(str(out_cool / "obs_exp.cool")).info["nbins"] == 10249
    # without the obs/exp transform nothing is written
    out_none = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed "
                          "--outFileObsExp {out}/obs_exp.h5", "none")
    assert not os.path.exists(str(out_none / "obs_exp.h5"))


def test_transforms_z_score_and_total_counts(tmp_path):
    zscore = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed --transform z-score", "z")
    # z-score leaves NaN wherever a diagonal has no spread, and the median of a
    # cell with a NaN is NaN
    assert np.isnan(_matrix(zscore, "m_genome.tab")).all()
    assert len(_pairs(zscore, "p_genome.tab")) == 443
    total = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed --transform total-counts", "t")
    _assert_matrix(_matrix(total, "m_genome.tab"), 0.3018, 0.00093, 0.00068, 447)


def test_row_wise_with_strand(tmp_path):
    out = _aggregate(tmp_path, "--matrix {root}Li_et_al_2015.h5 --mode intra-chr "
                     "--range 50000:6000000 --BED {agg}bed1_row-wise_strand.bed "
                     "--BED2 {agg}bed2_row-wise_strand.bed --row_wise --perChr --keep_outlier "
                     "--considerStrandDirection")
    _assert_matrix(_matrix(out, "m_X.tab"), 594.56093, 0.98566, 0.92617, 603)
    pairs = _pairs(out, "p_X.tab")
    assert len(pairs) == 10
    # defect 5: written as given, the larger start first
    assert pairs[0] == "X\t1775025\t1775026\tX\t1375616\t1375617\t2.901700013243535"


def test_row_wise_inter_chromosomal(tmp_path):
    out = _aggregate(tmp_path, "--matrix {root}small_test_matrix_50kb_res.h5 --mode inter-chr "
                     "--BED {agg}bed1_row-wise.bed --BED2 {agg}bed2_row-wise.bed --row_wise")
    assert int((_matrix(out, "m_genome.tab") != 0).sum()) == 1
    # an integer matrix keeps its dtype for the centre value
    assert _pairs(out, "p_genome.tab") == ["chrX\t1775025\t1775026\tchr3R\t1956476\t1956477\t0"]


def test_inter_chromosomal_mean(tmp_path):
    out = _aggregate(tmp_path, "--matrix {root}small_test_matrix.h5 --mode inter-chr "
                     "--BED {agg}test_regions.bed --operationType mean")
    matrix = _matrix(out, "m_genome.tab")
    assert matrix.sum() == pytest.approx(0.99999, rel=1e-9)
    assert int((matrix != 0).sum()) == 9
    assert len(_pairs(out, "p_genome.tab")) == 9


def test_per_chromosome_k_falls_to_one_and_stays(tmp_path):
    """Pinned defect 6: `k = 1` in cluster_matrices is never reset.

    chr3L comes first among the chromosomes with submatrices and has fewer
    than four, so it and every chromosome after it are clustered with k=1,
    while the names keep the `_cluster_1` suffix of --kmeans 4.
    """
    out = _aggregate(tmp_path, "--matrix {root}small_test_matrix.cool --mode intra-chr "
                     "--BED {agg}test_regions.bed --perChr --kmeans 4")
    assert _tables(out) == ["m_chr2L_cluster_1.tab", "m_chr2R_cluster_1.tab",
                            "m_chr3L_cluster_1.tab", "m_chrX_cluster_1.tab",
                            "p_chr2L_cluster_1.tab", "p_chr2R_cluster_1.tab",
                            "p_chr3L_cluster_1.tab", "p_chrX_cluster_1.tab"]
    assert len(_pairs(out, "p_chrX_cluster_1.tab")) == 38
    assert len(_pairs(out, "p_chr2L_cluster_1.tab")) == 2


def test_chromosomes_restore_the_masked_bins(tmp_path):
    """Pinned defect 7: keepOnlyTheseChr undoes maskBins and enlarge_bins.

    Li_et_al_2015 has 855 NaN bins. Without --chromosomes they are removed
    and the gaps closed, and 279 submatrices are aggregated; with
    --chromosomes X the bins come back with their original, gapped
    coordinates, two BED positions fall between bins, and 277 remain.
    """
    plain = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed", "plain")
    only_x = _aggregate(tmp_path, INTRA + "--BED {agg}test_regions.bed --chromosomes X", "x")
    assert len(_pairs(plain, "p_genome.tab")) == 279
    assert len(_pairs(only_x, "p_genome.tab")) == 277
    # an integer matrix turns float64 on the same path
    ints = _aggregate(tmp_path, "--matrix {root}small_test_matrix_50kb_res.h5 --mode intra-chr "
                      "--range 100000:20000000 --BED {agg}test_regions.bed "
                      "--chromosomes chrX chr2L --operationType mean", "ints")
    assert _pairs(ints, "p_genome.tab")[0] == "chrX\t1268823\t1268824\tchrX\t1356647\t1356648\t2.0"


def test_mode_all(tmp_path):
    out = _aggregate(tmp_path, "--matrix {root}Li_et_al_2015.h5 --mode all "
                     "--BED {agg}test_regions.bed --keep_outlier")
    _assert_matrix(_matrix(out, "m_genome.tab"), 715.4623, 1.34548, 1.27351, 685)
    pairs = _pairs(out, "p_genome.tab")
    assert len(pairs) == 308
    assert pairs[0] == "X\t3593688\t3593689\tX\t3596515\t3596516\t224.10480373945023"
