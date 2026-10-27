"""chicQualityControl: a thin layer over the C++ executable chicQualityControl.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). compute_sparsity keeps its 3.7.6 signature and behaviour and calls
hicexplorer._hicx.chicQualityControl (cpp/python/src/api/chicQualityControl.cpp);
the tool does not use it.
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
from hicexplorer._hicx import chicQualityControl as _impl
from .lib import Viewpoint


def parse_arguments(args=None):
    """The argparse parser of chicQualityControl, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicQualityControl')


def compute_sparsity(pReferencePoints, pViewpointObj, pArgs, pQueue):
    return _impl.compute_sparsity(pReferencePoints, pViewpointObj, pArgs, pQueue)


def main(args=None):
    """Run the C++ implementation of chicQualityControl with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicQualityControl')(args)
