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
    parserOpt.add_argument('--outputViewpointFile', '-ovf',
                           help='path to data file which holds the used data of the viewpoint and the backgroundmodel per bin.',
                           required=False)
    parserOpt.add_argument('--dpi',
                           help='Optional parameter: Resolution for the image in case the'
                           'output is a raster graphics image (e.g png, jpg)',
                           type=int,
                           default=300)

    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)
    viewpointObj = Viewpoint()
    utilitiesObj = Utilities()
    if args.backgroundModelFile:
        background_data = viewpointObj.readBackgroundDataFile(args.backgroundModelFile)

        background_data_sorted = sorted(background_data)
        background_data_list = list(background_data.values())

    for interactionFile in args.interactionFile:
        fig = plt.figure(figsize=(9.4, 4.8))
        ax = plt.subplot(111)
        plt.subplots_adjust(left=0.1, right=0.85)
        header, interaction_data, z_score_data, interaction_file_data_raw = viewpointObj.readInteractionFile(interactionFile)
        matrix_name, viewpoint, upstream_range, downstream_range = header.split('\t')
        viewpoint_ = viewpoint.split(':')
        referencePointString = viewpoint_[0] + ':' + utilitiesObj.in_units(viewpoint_[1]) + " - " + utilitiesObj.in_units(viewpoint_[2])

        matrix_name = matrix_name[1:]
        data = []
        z_score = []

        if args.backgroundModelFile:
            viewpoint_index = background_data_sorted.index(0)

            data_background = []
            data_background_mean = []

            for key in background_data_sorted:
                if key in interaction_data:
                    data.append(interaction_data[key])
                else:
                    data.append(0)

                if key in z_score_data:
                    z_score.append(z_score_data[key])
                else:
                    z_score.append(0)

                data_background.append(background_data[key][0])
                data_background_mean.append(background_data[key][1])

                if args.outputViewpointFile:

                    if key in interaction_file_data_raw:
                        line_data_raw = interaction_file_data_raw[key]
                        line_data_raw.append(background_data[key][0])
                        interaction_file_data_raw[key] = line_data_raw

                        line_data_raw = interaction_file_data_raw[key]
                        line_data_raw.append(background_data[key][1])
                        interaction_file_data_raw[key] = line_data_raw

        else:
            data = []
            interaction_key = sorted(interaction_data)
            for key in interaction_key:
                data.append(interaction_data[key])
            log.debug('data {}'.format(interaction_key))
            viewpoint_index = interaction_key.index(0)

        legend = [matrix_name, 'background model']
        data_plot_label = ax.plot(range(len(data)), data, alpha=0.9, label=matrix_name)

        # z_score
        ax_sub = ax.twinx()
        data_plot_label += ax_sub.plot(range(len(z_score)), z_score, '--y', alpha=0.7, label='z-score')
        ax_sub.set_ylabel('z-score', color='y')

        if args.backgroundModelFile:
            ax_sub_background = ax.twinx()
            ax_sub_background.spines["right"].set_position(("axes", 1.1))
            ax_sub_background.set_ylabel('background model', color='r')

            data_background = np.array(data_background)
            data_background_mean = np.array(data_background_mean)
            data_plot_label += ax_sub_background.plot(range(len(data_background)), data_background, '--r', alpha=0.5, label='background model')
            ax_sub_background.fill_between(range(len(data_background)), data_background + data_background_mean, data_background - data_background_mean, facecolor='red', alpha=0.3)

        ax.set_xticks([0, viewpoint_index, len(data)])

        ax.set_xticklabels([upstream_range, referencePointString, downstream_range])
        ax.set_ylabel('Number of interactions')

        # multiple legends in one figure
        data_legend = [label.get_label() for label in data_plot_label]
        ax.legend(data_plot_label, data_legend, loc=0)

        outFileName = '.'.join(args.outFileName.split('.')[:-1])
        fileFormat = args.outFileName.split('.')[-1]

        plt.savefig(outFileName + '_' + header + '.' + fileFormat, dpi=args.dpi)
        plt.close(fig)

        if args.outputViewpointFile:

            interaction_file_data_raw_sorted = sorted(interaction_file_data_raw)
            with open(args.outputViewpointFile, 'w') as output_file:
                output_file.write('#ChrViewpoint\tStart\tEnd\tChrInteraction\tStart\tEnd\tRelative position\tRelative Interactions\tZ-score')

                if args.backgroundModelFile:
                    output_file.write('\tbackground model\tbackground model SEM\n')
                else:
                    output_file.write('\n')
                for key in interaction_file_data_raw_sorted:

                    array_ = interaction_file_data_raw[key]

                    output_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(array_[0], array_[1], array_[2],
                                                                                  array_[3], array_[4], array_[5],
                                                                                  array_[6], array_[7], array_[8]))
                    if args.backgroundModelFile:
                        output_file.write('\t{}\t{}\n'.format(array_[9], array_[10]))
                    else:
                        output_file.write('\n')
