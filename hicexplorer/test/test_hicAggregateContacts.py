import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
import os.path
from tempfile import NamedTemporaryFile
import hicexplorer.hicAggregateContacts
import pytest
from psutil import virtual_memory
mem = virtual_memory()
memory = mem.total / 2**30

# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 7
HIGH_MEMORY = 200

REMOVE_OUTPUT = True
# DIFF = 60


ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
tolerance = 13  # default matplotlib pixed difference tolerance


@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)

    args = "--matrix {root}/Li_et_al_2015.cool --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000".\
        format(root=ROOT, out_agg=outfile_aggregate_plots.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate.png'

    hicexplorer.hicAggregateContacts.main(args.split())

    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res

    os.remove(outfile_aggregate_plots.name)


@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_clustering():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)
    outfile_heatmaps = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_heatmap_', delete=False)

    args = "--matrix {root}/Li_et_al_2015.cool --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --hclust 4 " \
           "--diagnosticHeatmapFile {out_heat} --howToCluster diagonal " \
           "--BED2 {root}/hicAggregateContacts/test_regions.bed".format(root=ROOT, out_agg=outfile_aggregate_plots.name,
                                                                        out_heat=outfile_heatmaps.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate_hclust4.png'
    test_image_heatmap = ROOT + 'hicAggregateContacts/master_heatmap.png'

    hicexplorer.hicAggregateContacts.main(args.split())

    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res

    res = compare_images(test_image_heatmap, outfile_heatmaps.name, tolerance)
    assert res is None, res

    os.remove(outfile_aggregate_plots.name)
    os.remove(outfile_heatmaps.name)


@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicAggregateContacts_3d():

    outfile_aggregate_3d = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_3d', delete=False)

    args = "--matrix {root}/Li_et_al_2015.cool --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --hclust 2 " \
           "--plotType 3d " \
           "--BED2 {root}/hicAggregateContacts/test_regions.bed".format(root=ROOT, out_agg=outfile_aggregate_3d.name)

    test_image_agg_3d = ROOT + 'hicAggregateContacts/master_aggregate_3d.png'

    hicexplorer.hicAggregateContacts.main(args.split())

    res = compare_images(test_image_agg_3d, outfile_aggregate_3d.name, tolerance)
    assert res is None, res

    os.remove(outfile_aggregate_3d.name)
