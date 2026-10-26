"""hicQuickQC: a thin layer over the C++ executable hicQuickQC.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

from tempfile import NamedTemporaryFile
import shutil
import argparse
import os
import errno

from hicexplorer import hicBuildMatrix as hicBuildMatrix
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """The argparse parser of hicQuickQC, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicQuickQC')


def main(args=None):
    """Run the C++ implementation of hicQuickQC with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicQuickQC')(args)
