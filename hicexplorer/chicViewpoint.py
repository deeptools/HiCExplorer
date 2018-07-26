# compute viewpoint interaction data and write it to an interaction file
# a)
# - option to smooth values
# - n matrices:
# - m referencePoints via input file: should one referencePoint via bash be supported?
# - output: n * m interaction files
# - +/- around it
# - option: values as counts or percentage
#

import argparse
import sys
import numpy as np
import hicexplorer.HiCMatrix as hm
from hicexplorer.utilities import toString
from .lib import Viewpoint
from .lib import Utilities
from hicexplorer._version import __version__
from scipy.stats import zscore
# from hicexplorer.chicViewpointBackgroundModel import getViewpointValues
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

import math
import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     description='Plots the number of interactions around a given reference point in a region.')

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='path of the Hi-C matrices to plot',
                                required=True,
                                nargs='+')

    parserRequired.add_argument('--range',
                                help='Defines the region upstream and downstream of a reference point which should be included. '
                                'Format is --region upstream downstream',
                                required=True,
                                type=int,
                                nargs=2)

    parserRequired.add_argument('--referencePoint', '-rp', help='Reference point file. Needs to be in the format: \'chr 100\' for a '
                                'single reference point or \'chr 100 200\' for a reference region and per line one reference point',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--averageContactBin',
                           help='Average the contacts of n bins, written to last column.',
                           type=int,
                           default=0)

    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)

    viewpointObj = Viewpoint()
    utilitiesObj = Utilities()

    referencePoints, gene_list = viewpointObj.readReferencePointFile(args.referencePoint)

    for matrix in args.matrices:
        hic_ma = hm.hiCMatrix(matrix)
        viewpointObj.hicMatrix = hic_ma

        for i, referencePoint in enumerate(referencePoints):
            region_start, region_end = viewpointObj.calculateViewpointRange(referencePoint, args.range)

            data_list = viewpointObj.computeViewpoint(referencePoint, referencePoint[0], region_start, region_end)
            if args.averageContactBin > 0:
                data_list = viewpointObj.smoothInteractionValues(data_list, args.averageContactBin)
            data_list_raw = np.copy(data_list)
            z_score_data = zscore(data_list)

            data_list = viewpointObj.computeRelativeValues(data_list)

            interaction_data = viewpointObj.createInteractionFileData(referencePoint, referencePoint[0], region_start, region_end, data_list, data_list_raw)

            referencePointString = '_'.join(str(j) for j in referencePoint)

            region_start_in_units = utilitiesObj.in_units(region_start)
            region_end_in_units = utilitiesObj.in_units(region_end)

            header_information = '\t'.join([matrix, referencePointString, str(region_start_in_units), str(region_end_in_units), gene_list[i]])
            header_information += '\n# ChrViewpoint\tStart\tEnd\tChrInteraction\tStart\tEnd\tRelative position\tRelative Interactions\tZ-score\tRaw\n#'
            matrix_name = '.'.join(matrix.split('.')[:-1])
            viewpointObj.writeInteractionFile('_'.join([matrix_name, referencePointString, gene_list[i]]), interaction_data, header_information, z_score_data)
