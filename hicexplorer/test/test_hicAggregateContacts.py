import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
import os.path
from tempfile import NamedTemporaryFile
import hicexplorer.hicAggregateContacts

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
tolerance = 13  # default matplotlib pixed difference tolerance


def test_hicAggregateContacts():

    outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)
    outfile_heatmaps = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_heatmap_', delete=False)

    args = "--matrix {root}/Li_et_al_2015.cool --BED {root}/hicAggregateContacts/test_regions.bed " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --hclust 4 " \
           "--diagnosticHeatmapFile {out_heat} --clusterOnDiagonal " \
           "--BED2 {root}/hicAggregateContacts/test_regions.bed".format(root=ROOT, out_agg=outfile_aggregate_plots.name,
                                                                        out_heat=outfile_heatmaps.name)

    test_image_agg = ROOT + 'hicAggregateContacts/master_aggregate.png'
    test_image_heatmap = ROOT + 'hicAggregateContacts/master_heatmap.png'

    hicexplorer.hicAggregateContacts.main(args.split())

    res = compare_images(test_image_agg, outfile_aggregate_plots.name, tolerance)
    assert res is None, res

    res = compare_images(test_image_heatmap, outfile_heatmaps.name, tolerance)
    assert res is None, res

    os.remove(outfile_aggregate_plots.name)
    os.remove(outfile_heatmaps.name)
