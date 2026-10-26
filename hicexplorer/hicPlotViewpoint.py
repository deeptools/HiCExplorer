"""hicPlotViewpoint: a thin layer over the C++ executable hicPlotViewpoint.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
import numpy as np
from hicmatrix import HiCMatrix as hm
from hicexplorer.utilities import toString

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """The argparse parser of hicPlotViewpoint, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicPlotViewpoint')


def relabelTicks(pTick):
    if pTick < 1e6:
        xlabels = "{:.2f} Kb".format(int(pTick) / 1e3)
    else:
        xlabels = "{:.2f} Mb".format(int(pTick) / 1e6)
    return xlabels


def getViewpointValues(pMatrix, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end, pInteractionList=None, pChromosome=None):

    hic = hm.hiCMatrix(pMatrix)
    if pChromosome is not None:
        hic.keepOnlyTheseChr(pChromosome)

    if len(pReferencePoint) == 2:
        view_point_start, view_point_end = hic.getRegionBinRange(pReferencePoint[0], int(pReferencePoint[1]), int(pReferencePoint[1]))
    elif len(pReferencePoint) == 3:
        view_point_start, view_point_end = hic.getRegionBinRange(pReferencePoint[0], int(pReferencePoint[1]), int(pReferencePoint[2]))
    else:
        log.error("No valid reference point given. {}".format(pReferencePoint))
        exit(1)

    view_point_range = hic.getRegionBinRange(pChromViewpoint, pRegion_start, pRegion_end)
    elements_of_viewpoint = view_point_range[1] - view_point_range[0]
    data_list = np.zeros(elements_of_viewpoint)
    view_point_start_ = view_point_start
    interactions_list = None
    if pInteractionList is not None:
        interactions_list = []
    while view_point_start_ <= view_point_end:
        chrom, start, end, _ = hic.getBinPos(view_point_start_)
        for j, idx in zip(range(elements_of_viewpoint), range(view_point_range[0], view_point_range[1], 1)):
            data_list[j] += hic.matrix[view_point_start_, idx]
            if interactions_list is not None:
                chrom_second, start_second, end_second, _ = hic.getBinPos(idx)
                interactions_list.append((chrom, start, end, chrom_second, start_second, end_second, hic.matrix[view_point_start_, idx]))
        view_point_start_ += 1

    return [view_point_start, view_point_end, view_point_range, data_list, interactions_list]


def main(args=None):
    """Run the C++ implementation of hicPlotViewpoint with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicPlotViewpoint')(args)
