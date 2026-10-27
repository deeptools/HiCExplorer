"""hicPlotViewpoint: a thin layer over the C++ executable hicPlotViewpoint.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). relabelTicks and getViewpointValues keep their 3.7.6
signatures and behaviour and call hicexplorer._hicx.hicPlotViewpoint
(cpp/python/src/api/hicPlotViewpoint.cpp); the tool does not use them.
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
from hicexplorer._hicx import hicPlotViewpoint as _impl

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """The argparse parser of hicPlotViewpoint, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicPlotViewpoint')


def relabelTicks(pTick):
    return _impl.relabelTicks(pTick)


def getViewpointValues(pMatrix, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end, pInteractionList=None, pChromosome=None):
    return _impl.getViewpointValues(pMatrix, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end, pInteractionList, pChromosome)


def main(args=None):
    """Run the C++ implementation of hicPlotViewpoint with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicPlotViewpoint')(args)
