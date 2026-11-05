"""hicBuildMatrixMicroC: a thin layer over the C++ executable hicBuildMatrixMicroC.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool; they are C++ too
(hicexplorer._hicx.lib_buildMatrixMethods through the star import below, and
hicexplorer._hicx.hicBuildMatrixMicroC for Tester).
"""

import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)


# own tools

import logging
log = logging.getLogger(__name__)

from hicexplorer.lib.buildMatrixMethods import *  # noqa: F401,F403 (tests and the equivalence cases use the build functions through this module)

from hicexplorer._hicx import hicBuildMatrixMicroC as _impl


def parse_arguments(args=None):
    """The argparse parser of hicBuildMatrixMicroC, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicBuildMatrixMicroC')


class Tester(object):
    def __init__(self):
        _impl.Tester___init__(self, __file__)


def main(args=None):
    """Run the C++ implementation of hicBuildMatrixMicroC with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicBuildMatrixMicroC')(args)
