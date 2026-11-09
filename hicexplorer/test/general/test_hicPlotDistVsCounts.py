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


def test_no_curve_to_draw_exits_1_with_a_message(tmp_path, capfd):
    """Issue #836: with --maxdepth equal to the bin size (5000) only the diagonal
    is left, no curve has more than one point, and 3.7.6 ended with
    UnboundLocalError for 'lgd'. The tool stops with a message and status 1."""
    plot = tmp_path / 'dv.png'
    args = "--matrices {} --plotFile {} --maxdepth 5000".format(ROOT + 'small_test_matrix.cool', plot).split()
    with pytest.raises(SystemExit) as exit_info:
        hicPlotDistVsCounts.main(args)
    assert exit_info.value.code == 1
    captured = capfd.readouterr()
    assert 'no curve with more than one point can be drawn' in captured.err
    assert 'UnboundLocalError' not in captured.err
    assert not plot.exists() or plot.stat().st_size == 0


def _plot_package_directory():
    """The directory holding hicexplorer_plot, where the tools look for it.

    The tools put <bin>/../plot (source and build tree) or
    <bin>/../share/hicexplorer/plot (installed package) on PYTHONPATH; the
    test looks in the source tree of this checkout and in the installed
    share directory of the running Python.
    """
    import sys
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [os.path.join(here, "..", "..", "..", "plot"),
                  os.path.join(sys.prefix, "share", "hicexplorer", "plot")]
    for candidate in candidates:
        if os.path.exists(os.path.join(candidate, "hicexplorer_plot", "__main__.py")):
            return os.path.abspath(candidate)
    pytest.skip("hicexplorer_plot is neither in the checkout nor installed")


def test_plot_layer_without_series_does_not_fail(tmp_path, monkeypatch):
    """The plot code survives a data file without any series (lgd = None)."""
    monkeypatch.syspath_prepend(_plot_package_directory())
    from hicexplorer_plot import hicPlotDistVsCounts as plot_layer
    plot = tmp_path / 'empty.png'
    plot_layer.draw({'plotFile': str(plot), 'perchr': False, 'num_chroms': 0, 'maxdepth': 5000,
                     'plotsize': None, 'matrices': [{'label': 'a', 'series': []}]})
    assert plot.exists()
