"""chicViewpointBackgroundModel: a thin layer over the C++ executable chicViewpointBackgroundModel.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). compute_background keeps its 3.7.6 signature and behaviour and calls
hicexplorer._hicx.chicViewpointBackgroundModel (cpp/python/src/api/chicViewpointBackgroundModel.cpp);
the tool does not use it.
"""
import logging
log = logging.getLogger(__name__)


from hicexplorer._hicx import chicViewpointBackgroundModel as _impl


def parse_arguments(args=None):
    """The argparse parser of chicViewpointBackgroundModel, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicViewpointBackgroundModel')


def compute_background(pReferencePoints, pViewpointObj, pArgs, pQueue):
    return _impl.compute_background(pReferencePoints, pViewpointObj, pArgs, pQueue)


def main(args=None):
    """Run the C++ implementation of chicViewpointBackgroundModel with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicViewpointBackgroundModel')(args)
