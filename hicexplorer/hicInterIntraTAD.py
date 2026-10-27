"""hicexplorer.hicInterIntraTAD, implemented in C++.

main() runs the C++ executable hicInterIntraTAD and parse_arguments() is built from its
--help-json (hicexplorer._cpp). The other functions keep their 3.7.6
signatures, docstrings and behaviour, and their bodies call
hicexplorer._hicx.hicInterIntraTAD (cpp/python/src/api/hicInterIntraTAD.cpp,
cpp/core/include/hicx/tad_contacts.hpp). The deviations from 3.7.6 are listed
in hicexplorer/test/api/cases_hicInterIntraTAD.py.
"""
import logging
log = logging.getLogger(__name__)

from hicexplorer.utilities import check_cooler
from hicexplorer._version import __version__

from hicexplorer._hicx import hicInterIntraTAD as _impl


def parse_arguments(args=None):
    """The argparse parser of hicInterIntraTAD, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicInterIntraTAD')


def readDomainBoundaries(pFile):
    return _impl.readDomainBoundaries(pFile)


def computeInterIntraTADs(pMatrix, pDomainList, pCoolOrH5, pThreadId, pQueue):
    return _impl.computeInterIntraTADs(pMatrix, pDomainList, pCoolOrH5, pThreadId, pQueue)


def main(args=None):
    """Run the C++ implementation of hicInterIntraTAD with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicInterIntraTAD')(args)
