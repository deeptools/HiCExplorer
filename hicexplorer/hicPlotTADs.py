import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from pygenometracks import plotTracks


def _main_python(args=None):

    plotTracks.main(args)


def main(args=None):
    """Run the C++ implementation; the Python implementation above is kept as _main_python."""
    from hicexplorer._cpp import entry_point
    entry_point('hicPlotTADs')(args)
