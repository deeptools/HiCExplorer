"""hicBuildMatrixMicroC: a thin layer over the C++ executable hicBuildMatrixMicroC.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import argparse
import os

import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)


# own tools
from hicexplorer.utilities import genomicRegion
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)

from hicexplorer.lib.buildMatrixMethods import *


def parse_arguments(args=None):
    """The argparse parser of hicBuildMatrixMicroC, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicBuildMatrixMicroC')


class Tester(object):
    def __init__(self):
        hic_test_data_dir = os.environ.get('HIC_TEST_DATA_DIR', False)
        if hic_test_data_dir:
            self.root = hic_test_data_dir
        else:
            self.root = os.path.dirname(
                os.path.abspath(__file__)) + "/test/test_data/"
        self.bam_file_1 = os.path.join(self.root, "hic.bam")


def main(args=None):
    """Run the C++ implementation of hicBuildMatrixMicroC with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicBuildMatrixMicroC')(args)
