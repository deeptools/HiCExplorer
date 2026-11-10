"""chicSignificantInteractions: a thin layer over the C++ executable chicSignificantInteractions.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other functions keep their 3.7.6 signatures and
behaviour and call hicexplorer._hicx.chicSignificantInteractions
(cpp/python/src/api/chicSignificantInteractions.cpp); the tool does not use
them. call_multi_core runs compute_interaction_file in this process, once per
thread.
"""
import logging
log = logging.getLogger(__name__)


from hicexplorer._hicx import chicSignificantInteractions as _impl


def parse_arguments(args=None):
    """The argparse parser of chicSignificantInteractions, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicSignificantInteractions')


def compute_interaction_file(pInteractionFilesList, pArgs, pViewpointObj, pBackground, pFilePath, pResolution, pQueue=None):
    return _impl.compute_interaction_file(pInteractionFilesList, pArgs, pViewpointObj, pBackground, pFilePath, pResolution, pQueue)


def filter_by_pvalue(pData, pPValue, pMergedLinesDict, pPeakInteractionsThreshold):
    return _impl.filter_by_pvalue(pData, pPValue, pMergedLinesDict, pPeakInteractionsThreshold)


def compute_new_p_values(pData, pBackgroundModel, pPValue, pMergedLinesDict, pPeakInteractionsThreshold):
    return _impl.compute_new_p_values(pData, pBackgroundModel, pPValue, pMergedLinesDict, pPeakInteractionsThreshold)


def merge_neighbors_x_fold(pXfold, pData, pViewpointObj, pResolution):
    return _impl.merge_neighbors_x_fold(pXfold, pData, pViewpointObj, pResolution)


def merge_neighbors_loose_p_value(pLoosePValue, pData, pViewpointObj, pResolution, pTruncateZeroPvalues):
    return _impl.merge_neighbors_loose_p_value(pLoosePValue, pData, pViewpointObj, pResolution, pTruncateZeroPvalues)


def call_multi_core(pInteractionFilesList, pArgs, pViewpointObj, pBackground, pFilePath, pResolution):
    return _impl.call_multi_core(pInteractionFilesList, pArgs, pViewpointObj, pBackground, pFilePath, pResolution)


def read_threshold_file(pFile):
    return _impl.read_threshold_file(pFile)


def writeSignificantHDF(pOutFileName, pSignificantDataList, pSignificantKeyList, pViewpointObj, pReferencePointsList, pArgs):
    return _impl.writeSignificantHDF(pOutFileName, pSignificantDataList, pSignificantKeyList, pViewpointObj, pReferencePointsList, pArgs)


def writeTargetHDF(pOutFileName, pTargetDataList, pTargetKeyList, pViewpointObj, pResolution, pReferencePointsList, pArgs):
    return _impl.writeTargetHDF(pOutFileName, pTargetDataList, pTargetKeyList, pViewpointObj, pResolution, pReferencePointsList, pArgs)


def main(args=None):
    """Run the C++ implementation of chicSignificantInteractions with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicSignificantInteractions')(args)
