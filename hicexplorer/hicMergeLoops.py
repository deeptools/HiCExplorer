"""hicMergeLoops: the tool is the C++ executable hicMergeLoops.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). mergeLoops and readFile keep their 3.7.6 signatures and
behaviour and call hicexplorer._hicx.hicMergeLoops
(cpp/python/src/api/hicMergeLoops.cpp, cpp/core/include/hicx/merge_loops.hpp).
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import logging
log = logging.getLogger(__name__)

from hicexplorer._version import __version__  # noqa: F401

from hicexplorer._hicx import hicMergeLoops as _impl


def parse_arguments(args=None):
    """The argparse parser of hicMergeLoops, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicMergeLoops')


def mergeLoops(pDataFrame, pLowestResolution, pTupleX, pTupleY):
    return _impl.mergeLoops(pDataFrame, pLowestResolution, pTupleX, pTupleY)


def readFile(pFile):
    return _impl.readFile(pFile)


def main(args=None):
    """Run the C++ implementation of hicMergeLoops with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicMergeLoops')(args)
