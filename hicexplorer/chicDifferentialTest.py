"""chicDifferentialTest: a thin layer over the C++ executable chicDifferentialTest.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other functions keep their 3.7.6 signatures and
behaviour and call hicexplorer._hicx.chicDifferentialTest
(cpp/python/src/api/chicDifferentialTest.cpp); the tool does not use them.
"""
import argparse
import sys
import errno
import os
import math
from multiprocessing import Process, Queue
import time
import logging
log = logging.getLogger(__name__)

import numpy as np
from scipy import stats
import h5py

import hicmatrix.HiCMatrix as hm
from hicexplorer import utilities
from hicexplorer._version import __version__
from hicexplorer._hicx import chicDifferentialTest as _impl
from .lib import Viewpoint


def parse_arguments(args=None):
    """The argparse parser of chicDifferentialTest, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicDifferentialTest')


def readInteractionFile(pInteractionFile):
    return _impl.readInteractionFile(pInteractionFile)


def chisquare_test(pDataFile1, pDataFile2, pAlpha):
    return _impl.chisquare_test(pDataFile1, pDataFile2, pAlpha)


def fisher_exact_test(pDataFile1, pDataFile2, pAlpha):
    return _impl.fisher_exact_test(pDataFile1, pDataFile2, pAlpha)


def writeResult(pOutFileName, pData, pHeaderOld, pHeaderNew, pAlpha, pTest):
    return _impl.writeResult(pOutFileName, pData, pHeaderOld, pHeaderNew, pAlpha, pTest)


def writeResultHDF(pOutFileName, pAcceptedData, pRejectedData, pAllResultData, pInputData, pAlpha, pTest):
    return _impl.writeResultHDF(pOutFileName, pAcceptedData, pRejectedData, pAllResultData, pInputData, pAlpha, pTest)


def run_statistical_tests(pInteractionFilesList, pArgs, pViewpointObject, pQueue=None):
    return _impl.run_statistical_tests(pInteractionFilesList, pArgs, pViewpointObject, pQueue)


def main(args=None):
    """Run the C++ implementation of chicDifferentialTest with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicDifferentialTest')(args)
