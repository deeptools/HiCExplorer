"""hicTransform: a thin layer over the C++ executable hicTransform.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The module has no other code.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)


import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """The argparse parser of hicTransform, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicTransform')


def main(args=None):
    """Run the C++ implementation of hicTransform with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicTransform')(args)
