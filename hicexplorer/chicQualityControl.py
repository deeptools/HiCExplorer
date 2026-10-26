"""chicQualityControl: a thin layer over the C++ executable chicQualityControl.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import argparse
import math
from multiprocessing import Process, Queue
import time
import os
import logging
logging.getLogger('hicmatrix').setLevel(logging.CRITICAL)
log = logging.getLogger(__name__)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from .lib import Viewpoint


def parse_arguments(args=None):
    """The argparse parser of chicQualityControl, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicQualityControl')


def compute_sparsity(pReferencePoints, pViewpointObj, pArgs, pQueue):

    sparsity_list = []
    try:
        chromosome_names = pViewpointObj.hicMatrix.getChrNames()

        for i, referencePoint in enumerate(pReferencePoints):

            if referencePoint is not None and referencePoint[0] in chromosome_names:

                region_start, region_end, _ = pViewpointObj.calculateViewpointRange(
                    referencePoint, (pArgs.fixateRange, pArgs.fixateRange))
                try:
                    data_list, _ = pViewpointObj.computeViewpoint(
                        referencePoint, referencePoint[0], region_start, region_end)
                    sparsity = (np.count_nonzero(data_list) / len(data_list))
                except (TypeError, IndexError):
                    sparsity = -1.0
                sparsity_list.append(sparsity)
            else:
                sparsity_list.append(-1.0)
    except Exception as exp:
        pQueue.put('Fail: ' + str(exp))
        return

    pQueue.put(sparsity_list)
    return


def main(args=None):
    """Run the C++ implementation of chicQualityControl with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicQualityControl')(args)
