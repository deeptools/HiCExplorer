"""hicexplorer.hicDifferentialTAD, implemented in C++.

main() runs the C++ executable hicDifferentialTAD and parse_arguments() is built from its
--help-json (hicexplorer._cpp). The other functions keep their 3.7.6
signatures, docstrings and behaviour, and their bodies call
hicexplorer._hicx.hicDifferentialTAD (cpp/python/src/api/hicDifferentialTAD.cpp,
cpp/core/include/hicx/tad_contacts.hpp). The deviations from 3.7.6 are listed
in hicexplorer/test/api/cases_hicDifferentialTAD.py.
"""
import logging
log = logging.getLogger(__name__)

from hicexplorer.utilities import check_cooler
from hicexplorer._version import __version__

from hicexplorer._hicx import hicDifferentialTAD as _impl


def parse_arguments(args=None):
    """The argparse parser of hicDifferentialTAD, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicDifferentialTAD')


def readDomainBoundaries(pFile):
    return _impl.readDomainBoundaries(pFile)


def computeDifferentialTADs(pMatrixTarget, pMatrixControl, pDomainList, pCoolOrH5, pPValue, pThreadId, pQueue):
    return _impl.computeDifferentialTADs(pMatrixTarget, pMatrixControl, pDomainList, pCoolOrH5, pPValue, pThreadId, pQueue)


def main(args=None):
    """Run the C++ implementation of hicDifferentialTAD with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicDifferentialTAD')(args)
