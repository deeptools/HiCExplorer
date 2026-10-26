"""chicViewpoint: a thin layer over the C++ executable chicViewpoint.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import argparse
import sys
import os
import errno
from multiprocessing import Process, Queue
import time
import math
import logging
log = logging.getLogger(__name__)

import h5py
import numpy as np

import hicmatrix.HiCMatrix as hm
from hicexplorer import utilities
from .lib import Viewpoint
from hicexplorer._version import __version__


def parse_arguments(args=None):
    """The argparse parser of chicViewpoint, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicViewpoint')


def adjustViewpointData(pViewpointObj, pData, pBackground, pReferencePoint, pRegionStart, pRegionEnd):
    data_viewpoint = {}
    data_background = {}
    view_point_start, _ = pViewpointObj.getReferencePointAsMatrixIndices(
        pReferencePoint)
    view_point_range_start, view_point_range_end = \
        pViewpointObj.getViewpointRangeAsMatrixIndices(
            pReferencePoint[0], pRegionStart, pRegionEnd)

    for i, data in zip(range(view_point_range_start, view_point_range_end, 1), pData):
        relative_position = i - view_point_start
        data_viewpoint[relative_position] = data
    for i, data in zip(range(view_point_range_start, view_point_range_end, 1), pBackground):
        relative_position = i - view_point_start
        data_background[relative_position] = data

    for i in data_background:
        if i in data_viewpoint:
            continue
        else:
            data_viewpoint[i] = 0

    data = np.fromiter(data_viewpoint.values(), dtype=np.float32)
    background = list(data_background.values())

    return data, background


def compute_x_fold(pDataList, pBackgroundList):
    return pDataList / pBackgroundList


def compute_viewpoint(pViewpointObj, pArgs, pQueue, pReferencePoints, pGeneList, pMatrix, pBackgroundModel, pBackgroundModelRelativeInteractions):
    interaction_data_list = []

    try:
        for i, referencePoint in enumerate(pReferencePoints):
            # range of viewpoint with reference point in the middle in genomic units
            # get fixateRange for relative interaction computation denominator
            region_start_fixed, region_end_fixed, range_fixed = pViewpointObj.calculateViewpointRange(
                referencePoint, (pArgs.fixateRange, pArgs.fixateRange))

            intermediate_viewpoint, _ = pViewpointObj.computeViewpoint(
                referencePoint, referencePoint[0], region_start_fixed, region_end_fixed)
            denominator_relative_interactions = np.sum(intermediate_viewpoint)

            # viewpoint data uses full range
            region_start, region_end, _range = pViewpointObj.calculateViewpointRange(
                referencePoint, pArgs.range)

            data_list, index_reference_point = pViewpointObj.computeViewpoint(
                referencePoint, referencePoint[0], region_start, region_end)

            # background uses fixed range, handles fixate range implicitly by same range used in background computation

            background_relative_interaction = pViewpointObj.interactionBackgroundData(
                pBackgroundModelRelativeInteractions, _range).flatten()
            data_list_relative = data_list
            if len(data_list) != len(background_relative_interaction):
                data_list, background_relative_interaction = adjustViewpointData(
                    pViewpointObj, data_list_relative, background_relative_interaction, referencePoint, region_start, region_end)

            if pArgs.averageContactBin > 0 and len(data_list) >= pArgs.averageContactBin:
                data_list = pViewpointObj.smoothInteractionValues(
                    data_list, pArgs.averageContactBin)

            data_list_raw = np.copy(data_list)

            data_list = pViewpointObj.computeRelativeValues(
                data_list, denominator_relative_interactions)

            x_fold_list = compute_x_fold(
                data_list, background_relative_interaction)
            p_value_list = pViewpointObj.pvalues(
                pBackgroundModel, data_list_raw, index_reference_point)

            # add values if range is larger than fixate range

            region_start_range, region_end_range, _ = pViewpointObj.calculateViewpointRange(
                referencePoint, (pArgs.range[0], pArgs.range[1]))

            interaction_data = pViewpointObj.createInteractionFileDataHDF5(referencePoint, referencePoint[0],
                                                                           region_start_range, region_end_range, data_list, data_list_raw,
                                                                           pGeneList[i], denominator_relative_interactions, p_value_list, x_fold_list)

            referencePointString = '_'.join(str(j) for j in referencePoint)

            region_start_in_units = utilities.in_units(region_start)
            region_end_in_units = utilities.in_units(region_end)
            denominator_relative_interactions_str = 'Sum of interactions in fixate range: '
            denominator_relative_interactions_str += str(
                denominator_relative_interactions)
            header_information = '# Interaction file, created with HiCExplorer\'s chicViewpoint version ' + \
                __version__ + '\n# '
            header_information += '\t'.join([pMatrix, referencePointString, str(region_start_in_units), str(
                region_end_in_units), pGeneList[i], denominator_relative_interactions_str])
            header_information += '\n# Chromosome\tStart\tEnd\tGene\tSum of interactions\tRelative position\tRelative Interactions\tp-value\tx-fold\tRaw\n#'
            matrix_name = '.'.join(pMatrix.split('/')[-1].split('.')[:-1])
            matrix_name = '_'.join(
                [matrix_name, referencePointString, pGeneList[i]])
            referencePointGenename = '_'.join([referencePointString, pGeneList[i]])
            referencePointGenename = str(referencePointGenename)
            interaction_data_list.append([matrix_name, interaction_data, header_information, referencePointGenename])
    except Exception as exp:
        log.debug('Error! {}'.format(str(exp)))
        pQueue.put('Fail: ' + str(exp))
        return
    pQueue.put(interaction_data_list)
    return


def main(args=None):
    """Run the C++ implementation of chicViewpoint with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicViewpoint')(args)
