from hicexplorer import hicPlotDistVsCounts
from tempfile import NamedTemporaryFile
import os

import matplotlib as mpl
mpl.use('agg')
import os.path


ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")


def test_plot():

    outfile = NamedTemporaryFile(suffix='.png', prefix='plotFile', delete=False)
    matrix = ROOT + 'small_test_matrix_50kb_res.h5'
    args = "--matrices {} --plotFile {} --plotsize 6 4".format(matrix, outfile.name).split()
    hicPlotDistVsCounts.main(args)

    # local computer: test passes with delta of 3000
    # travis: needs to be at least 4500 to pass
    # I love this voodoo :(
    size_new = os.path.getsize(outfile.name)
    size_reference = os.path.getsize(ROOT + 'hicPlotDistVsCounts/dist_vs_counts.png',)
    assert abs(size_new - size_reference) < 5000

    os.remove(outfile.name)
