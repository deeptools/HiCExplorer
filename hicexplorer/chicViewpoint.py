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

    parserRequired.add_argument('--outFileName', '-o', 
                                help='One bed file each reference point and for each matrix is created.'
                                ' Naming is: interactionOutFileName_matrixName_ReferencePoint.bed',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')


    parserOpt.add_argument('--averageContactBin',
                           help='Average the contacts of n bins, written to last column.',
                           type=int,
                           default=0)
    parserOpt.add_argument('--relativeValues',
                           help='Write relative contacts instead of absolute ones to interaction file. Relative to total number of counts of a viewpoint.',
                           required=False,
                           action='store_true')
    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)

    viewpointObj = Viewpoint()

    referencePoints = viewpointObj.readReferencePointFile(args.referencePoint)


    for matrix in args.matrices:
        hic_ma = hm.hiCMatrix(matrix)
        viewpointObj.setHiCMatrixObj(hic_ma)

        # bin_size = hic_ma.getBinSize()
        for referencePoint in referencePoints:
            log.debug('referencePoint {}'.format(referencePoint))
            region_start, region_end = viewpointObj.calculateViewpointRange(referencePoint, args.range)
            
            data_list = viewpointObj.computeViewpoint(referencePoint, referencePoint[0], region_start, region_end)

            if args.averageContactBin > 0:
                data_list = viewpointObj.smoothInteractionValues(data_list, args.averageContactBin)
            if args.relativeValues:
                viewpointObj.computeRelativeValues(data_list)
            interaction_data = viewpointObj.createInteractionFileData(referencePoint, referencePoint[0], region_start, region_end, data_list)
            referencePointString = ':'.join(str(i) for i in referencePoint)
            header_information = matrix + '\t' + referencePointString + '\t' + str(args.range[0]) + '\t' + str(args.range[1])
            viewpointObj.writeInteractionFile(args.outFileName + '_' + matrix + '_' + referencePointString, interaction_data, header_information)
            
