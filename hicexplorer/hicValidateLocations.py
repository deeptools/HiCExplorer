"""hicValidateLocations: a thin layer over the C++ executable hicValidateLocations.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other functions take and return the pandas objects
they did in HiCExplorer 3.7.6 and call hicexplorer._hicx.hicValidateLocations
(C++, cpp/python/src/api/hicValidateLocations.cpp), the computations the
executable runs.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import logging
log = logging.getLogger(__name__)


from hicexplorer._hicx import hicValidateLocations as _impl


def parse_arguments(args=None):
    """The argparse parser of hicValidateLocations, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicValidateLocations')


def readThreeColumn(pFile, pChrPrefixProtein):
    return _impl.readThreeColumn(pFile, pChrPrefixProtein)


def readLoopFile(pInputFile, pChrPrefixLoops):
    return _impl.readLoopFile(pInputFile, pChrPrefixLoops)


def overlapLoop(pDataFrameLoop, pDataFrameProtein):
    return _impl.overlapLoop(pDataFrameLoop, pDataFrameProtein)


def overlapTAD(pDataFrameTAD, pDataFrameProtein):
    return _impl.overlapTAD(pDataFrameTAD, pDataFrameProtein)


def correlateCool(pCoolFile, pDataFrameLoop):
    return _impl.correlateCool(pCoolFile, pDataFrameLoop)


# def correlate2DText(pDataFrameProtein, pDataFrameLoop):
#     loop_bedtool_x = BedTool.from_dataframe(pDataFrameLoop[[0, 1, 2, 3, 4, 5]])

#     protein_bedtool = BedTool.from_dataframe(pDataFrameProtein)
#     x = loop_bedtool_x.intersect(protein_bedtool, c=True).to_dataframe()

#     mask_x = x['name'] >= 1

#     log.debug(mask_x)
#     return mask_x


def applyBinning(pDataFrame, pBinSize, pMerge=True):
    return _impl.applyBinning(pDataFrame, pBinSize, pMerge)


def writeLoopFile(pOutFileName, pLoopDataFrame):
    return _impl.writeLoopFile(pOutFileName, pLoopDataFrame)


def main(args=None):
    """Run the C++ implementation of hicValidateLocations with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicValidateLocations')(args)
