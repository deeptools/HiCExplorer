"""hicAdjustMatrix: the tool is the C++ executable hicAdjustMatrix.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). adjustMatrix keeps its 3.7.6 signature and behaviour and
calls hicexplorer._hicx.hicAdjustMatrix (cpp/python/src/api/hicAdjustMatrix.cpp),
which loads the matrix with the C++ matrix layer.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse  # noqa: F401
from hicmatrix import HiCMatrix as hm  # noqa: F401
from hicexplorer._version import __version__  # noqa: F401
from hicexplorer.utilities import check_cooler  # noqa: F401 (adjustMatrix calls it)
import numpy as np  # noqa: F401
import cooler  # noqa: F401
import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import hicAdjustMatrix as _impl


def parse_arguments(args=None):
    """The argparse parser of hicAdjustMatrix, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicAdjustMatrix')


def adjustMatrix(pArgs):
    return _impl.adjustMatrix(pArgs)


def main(args=None):
    """Run the C++ implementation of hicAdjustMatrix with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicAdjustMatrix')(args)
