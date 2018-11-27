from pygenometracks import plotTracks

import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)


def main(args=None):

    plotTracks.main(args)
