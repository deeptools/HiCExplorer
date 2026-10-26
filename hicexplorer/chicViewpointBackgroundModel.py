"""chicViewpointBackgroundModel: a thin layer over the C++ executable chicViewpointBackgroundModel.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import argparse
import math
from multiprocessing import Process, Queue
import time
import logging
log = logging.getLogger(__name__)

import numpy as np
import fit_nbinom

from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from .lib import Viewpoint


def parse_arguments(args=None):
    """The argparse parser of chicViewpointBackgroundModel, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicViewpointBackgroundModel')


def compute_background(pReferencePoints, pViewpointObj, pArgs, pQueue):

    background_model_data = {}
    relative_positions = set()
    try:
        for i, referencePoint in enumerate(pReferencePoints):

            region_start, region_end, _ = pViewpointObj.calculateViewpointRange(
                referencePoint, (pArgs.fixateRange, pArgs.fixateRange))

            data_list, _ = pViewpointObj.computeViewpoint(
                referencePoint, referencePoint[0], region_start, region_end)

            # set data in relation to viewpoint, upstream are negative values, downstream positive, zero is viewpoint
            view_point_start, _ = pViewpointObj.getReferencePointAsMatrixIndices(
                referencePoint)
            view_point_range_start, view_point_range_end = \
                pViewpointObj.getViewpointRangeAsMatrixIndices(
                    referencePoint[0], region_start, region_end)

            for i, data in zip(range(view_point_range_start, view_point_range_end, 1), data_list):
                relative_position = i - view_point_start

            if pArgs.averageContactBin > 0:
                data_list = pViewpointObj.smoothInteractionValues(
                    data_list, pArgs.averageContactBin)

            for i, data in zip(range(view_point_range_start, view_point_range_end, 1), data_list):
                relative_position = i - view_point_start
                if relative_position in background_model_data:
                    background_model_data[relative_position].append(data)
                else:
                    background_model_data[relative_position] = [data]
                    relative_positions.add(relative_position)
    except Exception as exp:
        pQueue.put('Fail: ' + str(exp))
        return
    pQueue.put([background_model_data, relative_positions])
    return


def main(args=None):
    """Run the C++ implementation of chicViewpointBackgroundModel with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicViewpointBackgroundModel')(args)
