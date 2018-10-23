import argparse
import sys
import numpy as np
import hicmatrix.HiCMatrix as hm
from hicexplorer import utilities

from hicexplorer._version import __version__
from .lib import Viewpoint

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import os

import math
import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     description='Plots the number of interactions around a given reference point in a region.')

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--interactionFile', '-if',
                                help='path to the interaction files which should be used for plotting',
                                required=True,
                                nargs='+')

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the result.',
                                required=True)


    parserRequired.add_argument('--range',
                           help='Defines the region upstream and downstream of a reference point which should be included. '
                           'Format is --region upstream downstream',
                           required=True,
                           type=int,
                           nargs=2)
    parserRequired.add_argument('--rbzScore', '-rbz',
                           help='Detect all rbz score above this value as significant',
                           type=float,
                           default=1.96,
                           required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    
    
    parserOpt.add_argument("--mergeBins", "-mb", action='store_true', help="Merge neighboring significant interactions to one. The value is averaged.")


    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser



def main(args=None):
    args = parse_arguments().parse_args(args)
    viewpointObj = Viewpoint()
    background_data = None

    # read all interaction files.
    for interactionFile in args.interactionFile:
        _, interaction_data, z_score, interaction_file_data = viewpointObj.readInteractionFile(interactionFile)

    