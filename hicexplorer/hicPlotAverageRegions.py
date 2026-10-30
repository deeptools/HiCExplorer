"""hicPlotAverageRegions: a thin layer over the C++ executable hicPlotAverageRegions.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import logging
log = logging.getLogger(__name__)
import matplotlib
matplotlib.use('Agg')


def parse_arguments(args=None):
    """The argparse parser of hicPlotAverageRegions, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicPlotAverageRegions')


def main(args=None):
    """Run the C++ implementation of hicPlotAverageRegions with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicPlotAverageRegions')(args)
