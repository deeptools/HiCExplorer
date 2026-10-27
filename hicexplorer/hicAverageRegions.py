"""hicAverageRegions, implemented in C++.

main() runs the C++ executable hicAverageRegions and parse_arguments() is
built from its --help-json (hicexplorer._cpp). The other public functions
keep the name, the signature and the behaviour they had in HiCExplorer 3.7.6;
their bodies call hicexplorer._hicx.hicAverageRegions
(cpp/python/src/api/hicAverageRegions.cpp).
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import hicAverageRegions as _impl


def parse_arguments(args=None):
    """The argparse parser of hicAverageRegions, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicAverageRegions')


def calculateViewpointRange(pHiCMatrix, pViewpoint, pRange, pCoordinatesToBinMapping):
    '''
    This function computes the correct start and end position of a viewpoint given the reference and the range.
    '''
    return _impl.calculateViewpointRange(pHiCMatrix, pViewpoint, pRange, pCoordinatesToBinMapping)


def getBinIndices(pHiCMatrix, pViewpoint):
    return _impl.getBinIndices(pHiCMatrix, pViewpoint)


def calculateViewpointRangeBins(pHiCMatrix, pViewpoint, pRange, pCoordinatesToBinMapping):
    return _impl.calculateViewpointRangeBins(pHiCMatrix, pViewpoint, pRange, pCoordinatesToBinMapping)


def main(args=None):
    """Run the C++ implementation of hicAverageRegions with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicAverageRegions')(args)
