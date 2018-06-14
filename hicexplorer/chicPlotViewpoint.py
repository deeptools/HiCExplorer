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
# from hicexplorer.chicViewpointBackgroundModel import getViewpointValues
from .lib import Viewpoint

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


    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the image.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--backgroundModelFile', '-bmf',
                                help='path to the background file which should be used for plotting',
                                required=False)

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
    if args.backgroundModelFile:
        background_data = viewpointObj.readBackgroundDataFile(args.backgroundModelFile)
        # log.debug('background data {}'.format(background_data))
       
        background_data_sorted = sorted(background_data)
        # log.debug('background data {}'.format(background_data))
        background_data_list= list(background_data.values())
    

    for interactionFile in args.interactionFile:
        fig = plt.figure(figsize=(6.4, 4.8))
        ax = plt.subplot(111)
        header, interaction_data = viewpointObj.readInteractionFile(interactionFile)
        matrix_name, viewpoint, upstream_range, downstream_range = header.split('\t')
        matrix_name = matrix_name[1:]
        data = []

        
        if args.backgroundModelFile:
            viewpoint_index = background_data_sorted.index(0)

            data_background = []
            # interaction_data = sorted(interaction_data)
            for key in background_data_sorted:
                if key in interaction_data:
                    data.append(interaction_data[key])
                else:
                    data.append(0)
                data_background.append(background_data[key][0])

        else:
            data = list(interaction_data.values())
            viewpoint_index = list(interaction_data.keys()).index(0)

        
        # xticklabels = viewpointObj.createXlabels()
        # log.debug('header {}'.format(header))
        # log.debug('data {}'.format(data))


    # fig = plt.figure(figsize=(6.4, 4.8))
    # ax = plt.subplot(111)
    # matrices_plot_legend = []
    # for i, data in enumerate(data_list):
    #     matrices_plot_legend.append(ax.plot(range(len(data)), data, alpha=0.7, label=matrix_name_legend[i])[0])
    
    # # if 

        legend = [matrix_name, 'background model']
        ax.plot(range(len(data)), data, alpha=0.7, label=header)
        if args.backgroundModelFile:
            ax.plot(range(len(data_background)), data_background, '--r', alpha=0.7, label=header)
        ax.set_xticks([0, viewpoint_index, len(data)])

        ax.set_xticklabels([upstream_range, viewpoint, downstream_range])
        ax.set_ylabel('Number of interactions')
    #     # left, width = .45, .5
    #     # bottom, height = .25, .7
    #     # right = left + width
    #     # top = bottom + height
    #     # if 
        plt.legend(legend, loc='upper right')
    #     if not args.plotOneImage:
    #         plt.savefig(args.outFileName, dpi=args.dpi)
    #         plt.close(fig)
    # if args.plotOneImage:
        outFileName = '.'.join(args.outFileName.split('.')[:-1])
        fileFormat = args.outFileName.split('.')[-1]

        plt.savefig(outFileName + '_' + header + '.' + fileFormat, dpi=args.dpi)
        plt.close(fig)