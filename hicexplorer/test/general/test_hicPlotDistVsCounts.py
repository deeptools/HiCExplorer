import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicPlotDistVsCounts
from tempfile import NamedTemporaryFile
import os
from hicexplorer.test.test_compute_function import compute

import matplotlib as mpl
mpl.use('agg')
import os.path
import pytest


ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")
matrix = ROOT + 'small_test_matrix_50kb_res.h5'


def test_plot():

    outfile = NamedTemporaryFile(suffix='.png', prefix='plotFile', delete=False)
    matrix = ROOT + 'small_test_matrix_50kb_res.h5'
    args = "--matrices {} --plotFile {} --plotsize 8 4".format(matrix, outfile.name).split()
    hicPlotDistVsCounts.main(args)

    # local computer: test passes with delta of 3000
    # travis: needs to be at least 4500 to pass
    # I love this voodoo :(
    # After changing the plotsize, let see if 2000 works for both...
    size_new = os.path.getsize(outfile.name)
    size_reference = os.path.getsize(ROOT + 'hicPlotDistVsCounts/dist_vs_counts.png',)
    assert abs(size_new - size_reference) < 2000

    # os.remove(outfile.name)


@pytest.mark.parametrize("matrices", [matrix])  # required
@pytest.mark.parametrize("labels", ['', '--labels label_1 label_2'])  # don't know if this will work
@pytest.mark.parametrize("plotsize1, plotsize2", [[6, 4], [6, 5]])  # (6, 5 is default)
@pytest.mark.parametrize("perchr", ['--perchr', ''])  # Not sure if this works due to store_True option of argparse
@pytest.mark.parametrize("skipDiagonal", ['--skipDiagonal', ''])  # not sure if this will work, see comment on perchr
@pytest.mark.parametrize("maxdepth", [int(3e6), '1000000'])  # default is 3e6
@pytest.mark.parametrize("chromosomeExclude", ['chrX'])
def test_trivial_run(matrices, labels, plotsize1, plotsize2, perchr, skipDiagonal, maxdepth, chromosomeExclude):
    """
        Simple test for general behaviour with all commandline args.
    """

    outfile_new = NamedTemporaryFile(suffix='.png', prefix='plotFile', delete=False)
    outFileData_new = NamedTemporaryFile(suffix='.txt', prefix='dataFile', delete=False)
    args = "--matrices {} --plotFile {} {} --plotsize {} {} {} --outFileData {} --maxdepth {} --chromosomeExclude {}".format(
        matrices,
        outfile_new.name,
        labels,
        plotsize1,
        plotsize2,
        perchr,
        outFileData_new.name,
        maxdepth,
        chromosomeExclude,
    ).split()

    # with open("test.txt", 'a') as fo:
    #     fo.write(str(args))
    #     fo.write("\n")
    # hicPlotDistVsCounts.main(args)
    compute(hicPlotDistVsCounts.main, args, 5)
