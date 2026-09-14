"""hicPlotTADs: pygenometracks.plotTracks.main with the command line the C++
binary checked, as hicexplorer/hicPlotTADs.py calls it.

Data:
    argv   the command line tokens, without the C++-only --plotData
"""

import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

from pygenometracks import plotTracks  # noqa: E402


def draw(data):
    plotTracks.main(data["argv"])
