# as input 
# - n bedfiles with interactions
# - 1 background model file
# - plot:
#       - one image with all in one file
#       - create n bedfiles
# how to show background model? 

import argparse
import sys
import numpy as np
import hicexplorer.HiCMatrix as hm
from hicexplorer.utilities import toString
from hicexplorer.chicViewpointBackgroundModel import getViewpointValues
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

    parserRequired.add_argument('--interactionFile', '-if',
                                help='path to the interaction files which should be used for plotting',
                                required=True,
                                nargs='+')

    parserRequired.add_argument('--region',
                                help='The format is chr:start-end ',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the image.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')


    parserOpt.add_argument('--plotOneImage', '-i', 
                            help='Optional parameter:  If set all viewpoints are plotted in one image.',
                            required=False,
                            type=bool,
                            default=False)

    parserOpt.add_argument('--dpi',
                           help='Optional parameter: Resolution for the image in case the'
                           'ouput is a raster graphics image (e.g png, jpg)',
                           type=int,
                           default=300)

    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    return parser




def main(args=None):
    args = parse_arguments().parse_args(args)
    viewpointObj = Viewpoint()
    fig = plt.figure(figsize=(6.4, 4.8))
    ax = plt.subplot(111)

    for interactionFile in args.interactionFile:
        header, data = viewpointObj.readInteractionFile(interactionFile)
        # xticklabels = viewpointObj.createXlabels()


    # fig = plt.figure(figsize=(6.4, 4.8))
    # ax = plt.subplot(111)
    # matrices_plot_legend = []
    # for i, data in enumerate(data_list):
    #     matrices_plot_legend.append(ax.plot(range(len(data)), data, alpha=0.7, label=matrix_name_legend[i])[0])
    
    # if 

        ax.plot(range(len(data)), data, alpha=0.7, label=header)
        # ax.set_xticklabels(xticklabels)
        ax.set_ylabel('Number of interactions')
        # left, width = .45, .5
        # bottom, height = .25, .7
        # right = left + width
        # top = bottom + height
        # if 
        # plt.legend(header)
        if not args.plotOneImage:
            plt.savefig(args.outFileName, dpi=args.dpi)
            plt.close(fig)
    if args.plotOneImage:
        plt.savefig(args.outFileName, dpi=args.dpi)
        plt.close(fig)