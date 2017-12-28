from hicexplorer import hicPlotDistVsCounts
from tempfile import NamedTemporaryFile
import os

import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
import os.path


ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"


def test_plot():

    outfile = NamedTemporaryFile(suffix='.png', prefix='plotFile', delete=False)
    matrix = ROOT + 'small_test_matrix_50kb_res.h5'
    args = "--matrices {} --plotFile {}".format(matrix, outfile.name).split()
    hicPlotDistVsCounts.main(args)

    res = compare_images(ROOT + 'hicPlotDistVsCounts/dist_vs_counts.png', outfile.name, tol=40)
    assert res is None, res

    # os.remove(outfile.name)
