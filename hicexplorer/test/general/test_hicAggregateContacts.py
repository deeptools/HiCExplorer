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


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_intra_perChr():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)

    args = "--matrix {root}/Li_et_al_2015.h5 --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --disable_bbox_tight --dpi 100 "\
           "--mode intra-chr --perChr".\
        format(root=ROOT, out_agg=outfile_aggregate_plots.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate_intra_perChr.png'

    # hicexplorer.hicAggregateContacts.main(args.split())
    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)
    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res

    os.remove(outfile_aggregate_plots.name)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_intra():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)

    args = "--matrix {root}/Li_et_al_2015.h5 --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --disable_bbox_tight --dpi 100 "\
           "--mode intra-chr".\
        format(root=ROOT, out_agg=outfile_aggregate_plots.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate_intra.png'

    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)
    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res

    os.remove(outfile_aggregate_plots.name)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_inter():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)

    args = "--matrix {root}/small_test_matrix.h5 --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --disable_bbox_tight --dpi 100 "\
           "--mode inter-chr".format(root=ROOT, out_agg=outfile_aggregate_plots.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate_inter.png'

    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)
    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res

    os.remove(outfile_aggregate_plots.name)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_all():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)

    args = "--matrix {root}/Li_et_al_2015.h5 --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --disable_bbox_tight --dpi 100 "\
           "--mode all".format(root=ROOT, out_agg=outfile_aggregate_plots.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate_all.png'

    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)
    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res

    os.remove(outfile_aggregate_plots.name)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_chromosome_not_given():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)

    args = "--matrix {root}/Li_et_al_2015.h5 --BED {root}/hicAggregateContacts/test_regions_region_not_given.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --disable_bbox_tight --dpi 100 "\
           "--mode intra-chr --perChr ".format(root=ROOT, out_agg=outfile_aggregate_plots.name)


    # hicexplorer.hicAggregateContacts.main(args.split())
    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_cooler():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)

    args = "--matrix {root}/Li_et_al_2015.cool --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --disable_bbox_tight --dpi 100 "\
           "--mode intra-chr --perChr ".format(root=ROOT, out_agg=outfile_aggregate_plots.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate_intra_perChr.png'  # noqa: F841

    # hicexplorer.hicAggregateContacts.main(args.split())
    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)

    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res

    os.remove(outfile_aggregate_plots.name)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_clustering():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)
    outfile_heatmaps = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_heatmap_', delete=False)

    args = "--matrix {root}/Li_et_al_2015.h5 --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --hclust 4 " \
           "--diagnosticHeatmapFile {out_heat} --howToCluster diagonal  --disable_bbox_tight --dpi 100 " \
           "--BED2 {root}/hicAggregateContacts/test_regions.bed --mode intra-chr --perChr".\
           format(root=ROOT, out_agg=outfile_aggregate_plots.name, out_heat=outfile_heatmaps.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate_hclust4.png'
    test_image_heatmap = ROOT + 'hicAggregateContacts/master_heatmap.png'

    # hicexplorer.hicAggregateContacts.main(args.split())
    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)

    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res

    res = compare_images(test_image_heatmap, outfile_heatmaps.name, tolerance)
    assert res is None, res

    os.remove(outfile_aggregate_plots.name)
    os.remove(outfile_heatmaps.name)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_clustering_cool():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)
    outfile_heatmaps = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_heatmap_', delete=False)

    args = "--matrix {root}/Li_et_al_2015.cool --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --hclust 4 " \
           "--diagnosticHeatmapFile {out_heat} --howToCluster diagonal  --disable_bbox_tight --dpi 100 " \
           "--BED2 {root}/hicAggregateContacts/test_regions.bed  --mode intra-chr --perChr".\
           format(root=ROOT, out_agg=outfile_aggregate_plots.name, out_heat=outfile_heatmaps.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate_hclust4.png'
    test_image_heatmap = ROOT + 'hicAggregateContacts/master_heatmap.png'

    # hicexplorer.hicAggregateContacts.main(args.split())
    compute(hicexplorer.hicAggregateContacts.main, args.split(), 5)

    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res

    res = compare_images(test_image_heatmap, outfile_heatmaps.name, tolerance)
    assert res is None, res

    os.remove(outfile_aggregate_plots.name)
    os.remove(outfile_heatmaps.name)


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
           "--range 50000:6000000 --dpi 100 --mode intra-chr --perChr ".\
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
