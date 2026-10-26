"""hicNormalize: a thin layer over the C++ executable hicNormalize.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
import logging
log = logging.getLogger(__name__)

import numpy as np


def parse_arguments(args=None):
    """The argparse parser of hicNormalize, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicNormalize')


def main(args=None):
    """Run the C++ implementation of hicNormalize with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicNormalize')(args)
