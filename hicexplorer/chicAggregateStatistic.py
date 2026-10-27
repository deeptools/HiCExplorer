"""chicAggregateStatistic: a thin layer over the C++ executable chicAggregateStatistic.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other functions keep their 3.7.6 signatures and
behaviour and call hicexplorer._hicx.chicAggregateStatistic
(cpp/python/src/api/chicAggregateStatistic.cpp); the tool does not use them.
call_multi_core runs the worker function in this process, once per thread.
"""
import argparse
import sys
import errno
import os
import math
from multiprocessing import Process, Queue
import time
import traceback
import logging
log = logging.getLogger(__name__)

import numpy as np
import h5py

from intervaltree import IntervalTree, Interval
import hicmatrix.HiCMatrix as hm

from hicexplorer import utilities
from hicexplorer._version import __version__
from hicexplorer._hicx import chicAggregateStatistic as _impl
from .lib import Viewpoint


def parse_arguments(args=None):
    """The argparse parser of chicAggregateStatistic, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicAggregateStatistic')


def filter_scores_target_list(pScoresDictionary, pTargetFType, pTargetPosDict, pTargetList=None, pTargetIntervalTree=None, pTargetFile=None):
    return _impl.filter_scores_target_list(pScoresDictionary, pTargetFType, pTargetPosDict, pTargetList, pTargetIntervalTree, pTargetFile)


def writeAggregateHDF(pOutFileName, pOutfileNamesList, pAcceptedScoresList, pArgs):
    return _impl.writeAggregateHDF(pOutFileName, pOutfileNamesList, pAcceptedScoresList, pArgs)


def run_target_list_compilation(pInteractionFilesList, pTargetList, pTargetFType, pTargetPosDict, pArgs, pViewpointObj, pQueue=None, pOneTarget=False):
    return _impl.run_target_list_compilation(pInteractionFilesList, pTargetList, pTargetFType, pTargetPosDict, pArgs, pViewpointObj, pQueue, pOneTarget)


def call_multi_core(pInteractionFilesList, pTargetFileList, pTargetFType, pTargetPosDict, pFunctionName, pArgs, pViewpointObj):
    return _impl.call_multi_core(pInteractionFilesList, pTargetFileList, pTargetFType, pTargetPosDict, pFunctionName, pArgs, pViewpointObj)


def main(args=None):
    """Run the C++ implementation of chicAggregateStatistic with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicAggregateStatistic')(args)
