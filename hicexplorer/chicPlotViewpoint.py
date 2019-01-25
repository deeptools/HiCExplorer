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
                                help='File name to save the image.',
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
    parserOpt.add_argument('--colorMapZscore',
                           help='Color map to use for the z-score. Available '
                           'values can be seen here: '
                           'http://matplotlib.org/examples/color/colormaps_reference.html',
                           default='viridis')
    parserOpt.add_argument('--maxZscore', '-maz',
                           help='Maximal value for z-score. Values above are set to this value.',
                           type=float,
                           default=None)
    parserOpt.add_argument('--minZscore', '-mz',
                           help='Minimal value for z-score. Values below are set to this value.',
                           type=float,
                           default=None)
    parserOpt.add_argument('--rbzScore', '-rbz',
                           help='Plot rbz-score as a colorbar',
                           choices=['integrated', 'heatmap', ''],
                           default=''
                           )

    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)
    viewpointObj = Viewpoint()
    background_data = None

    if args.backgroundModelFile:
        background_data = viewpointObj.readBackgroundDataFile(args.backgroundModelFile, args.range)

    number_of_rows_plot = len(args.interactionFile)
    matplotlib.rcParams.update({'font.size': 9})
    fig = plt.figure(figsize=(9.4, 4.8))

    z_score_heights = [0.07] * number_of_rows_plot
    viewpoint_height_ratio = 0.95 - (0.07 * number_of_rows_plot)
    if viewpoint_height_ratio < 0.4:
        viewpoint_height_ratio = 0.4
        _ratio = 0.6 / len(number_of_rows_plot)
        z_score_heights = [_ratio] * number_of_rows_plot

    if args.rbzScore == 'heatmap':
        gs = gridspec.GridSpec(1 + len(args.interactionFile), 2, height_ratios=[0.95 - (0.07 * number_of_rows_plot), *z_score_heights], width_ratios=[0.75, 0.25])
        gs.update(hspace=0.5, wspace=0.05)
        ax1 = plt.subplot(gs[0, 0])
        ax1.margins(x=0)
    else:
        ax1 = plt.subplot()
    colors = ['g', 'b', 'c', 'm', 'y', 'k']
    background_plot = True
    data_plot_label = None
    for i, interactionFile in enumerate(args.interactionFile):

        header, data, background_data_plot, data_background_mean, z_score, interaction_file_data_raw, viewpoint_index = viewpointObj.getDataForPlotting(interactionFile, args.range, background_data)
        if len(data) <= 1 or len(z_score) <= 1:
            log.warning('Only one data point in given range, no plot is created! Interaction file {} Range {}'.format(interactionFile, args.range))
            continue
        matrix_name, viewpoint, upstream_range, downstream_range, gene, _ = header.strip().split('\t')
        matrix_name = matrix_name[1:].split('.')[0]
        if data_plot_label:
            data_plot_label += viewpointObj.plotViewpoint(pAxis=ax1, pData=data, pColor=colors[i % len(colors)], pLabelName=gene + ': ' + matrix_name)
        else:
            data_plot_label = viewpointObj.plotViewpoint(pAxis=ax1, pData=data, pColor=colors[i % len(colors)], pLabelName=gene + ': ' + matrix_name)

        if background_plot:
            data_plot_label += viewpointObj.plotBackgroundModel(pAxis=ax1, pBackgroundData=background_data_plot,
                                                                pBackgroundDataMean=data_background_mean)
            background_plot = False
        if args.minZscore is not None or args.maxZscore is not None:
            z_score = np.array(z_score, dtype=np.float32)
            z_score.clip(args.minZscore, args.maxZscore, z_score)
        if args.rbzScore == 'heatmap':
            viewpointObj.plotZscore(pAxis=plt.subplot(gs[1 + i, 0]), pAxisLabel=plt.subplot(gs[1 + i, 1]), pZscoreData=z_score,
                                    pLabelText=gene + ': ' + matrix_name, pCmap=args.colorMapZscore,
                                    pFigure=fig,)
        elif args.rbzScore == 'integrated':
            data_plot_label += viewpointObj.plotViewpoint(pAxis=ax1, pData=z_score, pColor=colors[i % len(colors)], pLabelName=gene + ': ' + matrix_name + ' rbz-score')

        viewpointObj.writePlotData(interaction_file_data_raw, matrix_name + '_' + gene + '_raw_plot', args.backgroundModelFile)

    if data_plot_label is not None:
        ax1.set_ylabel('Number of interactions')
        ax1.set_xticks([0, viewpoint_index, len(data) - 1])

        if args.range:
            ax1.set_xticklabels([str(-args.range[0]), 'Viewpoint', str(args.range[1])])
        else:
            ax1.set_xticklabels([upstream_range, 'Viewpoint', downstream_range])

        # multiple legends in one figure
        data_legend = [label.get_label() for label in data_plot_label]
        ax1.legend(data_plot_label, data_legend, loc=0)
        plt.savefig(args.outFileName, dpi=args.dpi)
    plt.close(fig)
