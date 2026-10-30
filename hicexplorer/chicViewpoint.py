"""chicViewpoint: a thin layer over the C++ executable chicViewpoint.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other functions keep their 3.7.6 signatures and
behaviour and call hicexplorer._hicx.chicViewpoint
(cpp/python/src/api/chicViewpoint.cpp); the tool does not use them.
"""
import logging
log = logging.getLogger(__name__)


from hicexplorer._hicx import chicViewpoint as _impl


def parse_arguments(args=None):
    """The argparse parser of chicViewpoint, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicViewpoint')


def adjustViewpointData(pViewpointObj, pData, pBackground, pReferencePoint, pRegionStart, pRegionEnd):
    return _impl.adjustViewpointData(pViewpointObj, pData, pBackground, pReferencePoint, pRegionStart, pRegionEnd)


def compute_x_fold(pDataList, pBackgroundList):
    return _impl.compute_x_fold(pDataList, pBackgroundList)


def compute_viewpoint(pViewpointObj, pArgs, pQueue, pReferencePoints, pGeneList, pMatrix, pBackgroundModel, pBackgroundModelRelativeInteractions):
    return _impl.compute_viewpoint(pViewpointObj, pArgs, pQueue, pReferencePoints, pGeneList, pMatrix, pBackgroundModel, pBackgroundModelRelativeInteractions)


def main(args=None):
    """Run the C++ implementation of chicViewpoint with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicViewpoint')(args)
