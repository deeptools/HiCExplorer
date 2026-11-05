"""hicBuildMatrix: a thin layer over the C++ executable hicBuildMatrix.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""

import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)


# own tools

import logging
log = logging.getLogger(__name__)

from hicexplorer.lib.buildMatrixMethods import *  # noqa: F401,F403 (the module re-exports the build functions)


def parse_arguments(args=None):
    """The argparse parser of hicBuildMatrix, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicBuildMatrix')


def main(args=None):
    """Run the C++ implementation of hicBuildMatrix with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicBuildMatrix')(args)
