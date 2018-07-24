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
from hicexplorer._version import __version__
from .lib import Viewpoint
from .lib import Utilities

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
                                help='File name to save the image.',
                                required=True)
    parserRequired.add_argument('--outputViewpointFile', '-ovf',
                                help='path to data file which holds the used data of the viewpoint and the backgroundmodel per bin.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--backgroundModelFile', '-bmf',
                           help='path to the background file which should be used for plotting',
                           required=False)

    parserOpt.add_argument('--dpi',
                           help='Optional parameter: Resolution for the image in case the'
                           'output is a raster graphics image (e.g png, jpg)',
                           type=int,
                           default=300)
    parserOpt.add_argument('--range',
                           help='Defines the region upstream and downstream of a reference point which should be included. '
                           'Format is --region upstream downstream',
                           required=False,
                           type=int,
                           nargs=2)

    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)
    viewpointObj = Viewpoint()
    utilitiesObj = Utilities()
    background_data = None
    background_data_sorted = None

    if args.backgroundModelFile:
        background_data = viewpointObj.readBackgroundDataFile(args.backgroundModelFile)

        # background_data_list = list(background_data.values())

    number_of_rows_plot = len(args.interactionFile)

    fig = plt.figure(figsize=(9.4, 4.8))
    # ax = plt.subplot(111)
    z_score_heights = [0.1] * number_of_rows_plot
    # TODO add check that ratio is greater than 0

    gs = gridspec.GridSpec(1 + len(args.interactionFile), 2, height_ratios=[0.90 - (0.1 * number_of_rows_plot), *z_score_heights], width_ratios=[0.97, 0.03])
    gs.update(hspace=0.5, wspace=0.05)
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0])

    # z_score_axis = ax2.subplot(len(args.interactionFile))
    colors = ['g', 'b', 'c', 'm', 'y', 'k']
    background_plot = True
    data_plot_label = None
    for i, interactionFile in enumerate(args.interactionFile):

        header, data, background_data_plot, data_background_mean, z_score, interaction_file_data_raw, viewpoint_index = viewpointObj.getDataForPlotting(interactionFile, args.range, background_data)
        log.debug('header {}'.format(header))
        matrix_name, viewpoint, upstream_range, downstream_range, gene = header.split('\t')

        if data_plot_label:
            data_plot_label += viewpointObj.plotViewpoint(pAxis=ax1, pData=data, pColor=colors[i % len(colors)], pLabelName=gene)
        else:
            data_plot_label = viewpointObj.plotViewpoint(pAxis=ax1, pData=data, pColor=colors[i % len(colors)], pLabelName=gene)

        if background_plot:
            data_plot_label += viewpointObj.plotBackgroundModel(pAxis=ax1, pBackgroundData=background_data_plot, pBackgroundDataMean=data_background_mean)
            background_plot = False

        viewpointObj.plotZscore(pAxis=plt.subplot(gs[1 + i, 0]), pAxisLabel=plt.subplot(gs[1 + i, 1]), pZscoreData=z_score, pLabelText=gene)

        # if args.outputViewpointFile:
        viewpointObj.writePlotData(interaction_file_data_raw, args.outputViewpointFile, args.backgroundModelFile)

    ax1.set_ylabel('Number of interactions')
    ax1.set_xticks([0, viewpoint_index, len(background_data_plot)])
    ax1.set_xticklabels([upstream_range, 'Viewpoint', downstream_range])

    # multiple legends in one figure
    data_legend = [label.get_label() for label in data_plot_label]
    ax1.legend(data_plot_label, data_legend, loc=0)
    # outFileName = '.'.join(args.outFileName.split('.')[:-1])
    # fileFormat = args.outFileName.split('.')[-1]
    # plt.savefig(outFileName + '_' + header + '.' + fileFormat, dpi=args.dpi)
    plt.savefig(args.outFileName, dpi=args.dpi)
    plt.close(fig)
