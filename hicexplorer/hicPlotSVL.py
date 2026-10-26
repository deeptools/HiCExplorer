"""hicPlotSVL, implemented in C++.

main() runs the C++ executable hicPlotSVL and parse_arguments() is built from
its --help-json (hicexplorer._cpp). compute_relation_short_long_range keeps
the name, the signature and the behaviour it had in HiCExplorer 3.7.6; its
body calls hicexplorer._hicx.hicPlotSVL (cpp/python/src/api/hicPlotSVL.cpp).
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import matplotlib
matplotlib.use('Agg')
import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import hicPlotSVL as _impl


def parse_arguments(args=None):
    """The argparse parser of hicPlotSVL, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicPlotSVL')


def compute_relation_short_long_range(pHiCMatrix, pChromosomes, pDistance, pIsCooler, pQueue):
    return _impl.compute_relation_short_long_range(pHiCMatrix, pChromosomes, pDistance, pIsCooler, pQueue)


def main(args=None):
    """Run the C++ implementation of hicPlotSVL with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicPlotSVL')(args)
