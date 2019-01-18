import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicPlotDistVsCounts
from tempfile import NamedTemporaryFile
import os

import matplotlib as mpl
mpl.use('agg')
import os.path
import pytest


ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")
matrix = ROOT + 'small_test_matrix_50kb_res.h5'
outfile = NamedTemporaryFile(suffix='.png', prefix='plotFile', delete=False)
outFileData = NamedTemporaryFile(suffix='.txt', prefix='dataFile', delete=True)


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


@pytest.mark.parametrize("matrices", [matrix])  # required
@pytest.mark.parametrize("plotFile", [outfile])  # required
@pytest.mark.parametrize("labels", [None, ['label_1', 'label_2']])  # don't know if this will work
@pytest.mark.parametrize("plotsize1, plotsize2", [[6, 4], [6, 5]])  # (6, 5 is default)
@pytest.mark.parametrize("perchr", ['--perchr', ''])  # Not sure if this works due to store_True option of argparse
@pytest.mark.parametrize("outFileData", [outFileData])
@pytest.mark.parametrize("skipDiagonal", ['--skipDiagonal', ''])  # not sure if this will work, see comment on perchr
@pytest.mark.parametrize("maxdepth", [int(3e6)])  # default is 3e6
@pytest.mark.parametrize("chromosomeExclude", ['chrX'])
def test_trivial_run(matrices, plotFile, labels, plotsize1, plotsize2, perchr, outFileData, skipDiagonal, maxdepth, chromosomeExclude):
    """
        Simple test for general behaviour with all commandline args.
    """

    args = "--matrices {} --plotFile {} --labels {} --plotsize {} {} {} --outFileData {} --maxdepth {} --chromosomeExclude {}".format(
        matrices,
        plotFile.name,
        labels,
        plotsize1,
        plotsize2,
        perchr,
        outFileData.name,
        maxdepth,
        chromosomeExclude,
    ).split()

    hicPlotDistVsCounts.main(args)
