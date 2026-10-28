"""chicExportData: a thin layer over the C++ executable chicExportData.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). exportData keeps its 3.7.6 signature and behaviour and calls
hicexplorer._hicx.chicExportData (cpp/python/src/api/chicExportData.cpp); the tool does not use it.
"""
import argparse
import sys
import os
import errno
import math
from multiprocessing import Process, Queue
import time
import traceback

import logging
log = logging.getLogger(__name__)

import numpy as np
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

import hicmatrix.HiCMatrix as hm
from hicexplorer import utilities
from hicexplorer._version import __version__
from hicexplorer._hicx import chicExportData as _impl
from .lib import Viewpoint

import h5py
import io
import tarfile
from contextlib import closing
from collections import OrderedDict

from tempfile import mkdtemp
import shutil


def parse_arguments(args=None):
    """The argparse parser of chicExportData, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicExportData')


def exportData(pFileList, pArgs, pViewpointObject, pDecimalPlace, pChromosomeSizes, pBackgroundData, pFileType, pQueue):
    return _impl.exportData(pFileList, pArgs, pViewpointObject, pDecimalPlace, pChromosomeSizes, pBackgroundData, pFileType, pQueue)


def main(args=None):
    """Run the C++ implementation of chicExportData with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicExportData')(args)
