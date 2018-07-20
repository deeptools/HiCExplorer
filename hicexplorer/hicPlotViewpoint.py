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

    parserRequired.add_argument('--matrix', '-m',
                                help='path of the Hi-C matrices to plot',
                                required=True,
                                nargs='+')

    parserRequired.add_argument('--region',
                                help='The format is chr:start-end ',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the image.',
                                required=True)

    parserRequired.add_argument('--referencePoint', '-rp', help='Reference point. Needs to be in the format: \'chr:100\' for a '
                                'single reference point or \'chr:100-200\' for a reference region.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--chromosome', '-C',
                           help='Optional parameter: Only show results for this chromosome.')

    parserOpt.add_argument('--interactionOutFileName', '-i', help='Optional parameter:  If set a bedgraph file with all interaction'
                           ' will be created.',
                           required=False)

    parserOpt.add_argument('--dpi',
                           help='Optional parameter: Resolution for the image in case the'
                           'ouput is a raster graphics image (e.g png, jpg)',
                           type=int,
                           default=300)
    parserOpt.add_argument('--averageContactBin',
                           help='Average the contacts of n bins, written to last column.',
                           type=int,
                           default=0)
    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    return parser


def relabelTicks(pTick):
    if pTick < 1e6:
        xlabels = "{:.2f} Kb".format(int(pTick) / 1e3)
    else:
        xlabels = "{:.2f} Mb".format(int(pTick) / 1e6)
    return xlabels


# def getViewpointValues(pMatrix, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end, pInteractionList=None, pChromosome=None):

#     hic = hm.hiCMatrix(pMatrix)
#     if pChromosome is not None:
#         hic.keepOnlyTheseChr(pChromosome)

#     if len(pReferencePoint) == 2:
#         view_point_start, view_point_end = hic.getRegionBinRange(pReferencePoint[0], int(pReferencePoint[1]), int(pReferencePoint[1]))
#     elif len(pReferencePoint) == 3:
#         view_point_start, view_point_end = hic.getRegionBinRange(pReferencePoint[0], int(pReferencePoint[1]), int(pReferencePoint[2]))
#     else:
#         log.error("No valid reference point given. {}".format(pReferencePoint))
#         exit(1)

#     view_point_range = hic.getRegionBinRange(pChromViewpoint, pRegion_start, pRegion_end)
#     elements_of_viewpoint = view_point_range[1] - view_point_range[0]
#     data_list = np.zeros(elements_of_viewpoint)
#     view_point_start_ = view_point_start
#     interactions_list = None
#     if pInteractionList is not None:
#         interactions_list = []
#     while view_point_start_ <= view_point_end:
#         chrom, start, end, _ = hic.getBinPos(view_point_start_)
#         for j, idx in zip(range(elements_of_viewpoint), range(view_point_range[0], view_point_range[1], 1)):
#             data_list[j] += hic.matrix[view_point_start_, idx]
#             if interactions_list is not None:
#                 chrom_second, start_second, end_second, _ = hic.getBinPos(idx)
#                 interactions_list.append((chrom, start, end, chrom_second, start_second, end_second, hic.matrix[view_point_start_, idx]))
#         view_point_start_ += 1

#     return [view_point_start, view_point_end, view_point_range, data_list, interactions_list]


def main(args=None):
    args = parse_arguments().parse_args(args)
    log.warning('This tool is depricated. Please use chicViewpoint, chicViewpointBackgroundModel and chicPlotViewpoint.')
    if args.region:
        if sys.version_info[0] == 2:
            args.region = args.region.translate(None, ",.;|!{}()").replace("-", ":")
        if sys.version_info[0] == 3:
            args.region = args.region.replace(",", "")
            args.region = args.region.replace(";", "")
            args.region = args.region.replace("!", "")
            args.region = args.region.replace("-", ":")
        region = args.region.split(":")
        if len(region) != 3:
            log.error("Region format is invalid {}".format(args.region))
            exit(0)
        chrom, region_start, region_end = region[0], int(region[1]), int(region[2])
    if sys.version_info[0] == 2:
        args.referencePoint = args.referencePoint.translate(None, ",.;|!{}()").replace("-", ":")
    if sys.version_info[0] == 3:
        args.referencePoint = args.referencePoint.replace(",", "")
        args.referencePoint = args.referencePoint.replace(";", "")
        args.referencePoint = args.referencePoint.replace("!", "")
        args.referencePoint = args.referencePoint.replace("-", ":")
    referencePoint = args.referencePoint.split(":")

    data_list = []
    interactions_list = None
    if args.interactionOutFileName is not None:
        interactions_list = []
    matrix_name_legend = []
    for matrix in args.matrix:
        hic_matrix = hm.hiCMatrix(matrix)
        if args.chromosome is not None:
            hic_matrix.keepOnlyTheseChr(args.chromosome)
        view_point_start, view_point_end, view_point_range, data_list_, interactions_list_ \
            = getViewpointValues(hic_matrix, referencePoint, chrom, region_start, region_end, args.interactionOutFileName, args.chromosome)
        data_list.append(data_list_)
        if args.interactionOutFileName is not None:
            interactions_list.append(interactions_list_)
        matrix_name_legend.append(os.path.basename(matrix))

    fig = plt.figure(figsize=(6.4, 4.8))
    ax = plt.subplot(111)
    matrices_plot_legend = []
    for i, data in enumerate(data_list):
        matrices_plot_legend.append(ax.plot(range(len(data)), data, alpha=0.7, label=matrix_name_legend[i])[0])
    if len(referencePoint) == 2:
        log.debug("Single reference point mode: {}".format(referencePoint))
        log.debug("label 0: {}".format((int(referencePoint[1]) - region_start) * (-1)))
        log.debug("referencePoint[1]: {}".format(referencePoint[1]))
        log.debug("region_start: {}".format(region_start))
        log.debug("label 1: {}".format(referencePoint[0] + ":" + relabelTicks(int(referencePoint[1]))))
        log.debug("label 2: {}".format(region_end - int(referencePoint[1])))

        ax.set_xticks([0, view_point_start - view_point_range[0], view_point_range[1] - view_point_range[0]])
        xticklabels = [None] * 3
        xticklabels[0] = relabelTicks((int(referencePoint[1]) - region_start) * (-1))
        xticklabels[1] = referencePoint[0] + ":" + relabelTicks(int(referencePoint[1]))
        xticklabels[2] = relabelTicks(region_end - int(referencePoint[1]))

    elif len(referencePoint) == 3:
        log.debug("Range mode: {}".format(referencePoint))

        # fit scale: start coordinate is 0 --> view_point_range[0]
        ax.set_xticks([0, view_point_start - view_point_range[0], view_point_end - view_point_range[0], view_point_range[1] - view_point_range[0]])
        xticklabels = [None] * 4
        xticklabels[0] = relabelTicks((int(referencePoint[1]) - region_start) * (-1))
        xticklabels[1] = referencePoint[0] + ":" + relabelTicks(int(referencePoint[1]))
        xticklabels[2] = referencePoint[0] + ":" + relabelTicks(int(referencePoint[2]))
        xticklabels[3] = relabelTicks(region_end - int(referencePoint[1]))

    ax.set_xticklabels(xticklabels)
    ax.set_ylabel('Number of interactions')
    # left, width = .45, .5
    # bottom, height = .25, .7
    # right = left + width
    # top = bottom + height

    plt.legend(handles=matrices_plot_legend)
    plt.savefig(args.outFileName, dpi=args.dpi)
    plt.close(fig)

    if interactions_list is not None:

        # averageContactBin

        if args.averageContactBin > 0:
            average_contacts = []

            window_size = np.int(math.floor(args.averageContactBin / 2))
            window_size_upstream = window_size
            if args.averageContactBin % 2 == 0:
                window_size_upstream -= 1
            for interactions_list_ in interactions_list:
                interactions = np.swapaxes(interactions_list_, 0, 1)[6].astype(float)
                average_contacts_ = np.zeros(len(interactions))

                # add upstream and downstream, handle regular case
                for i in range(window_size_upstream, len(interactions) - window_size):
                    start = i - window_size_upstream
                    end = i + window_size + 1
                    if i == window_size_upstream:
                        log.info('args.averageContactBin {} window_size_upstream {} window_size {}'.format(args.averageContactBin, window_size_upstream, window_size))
                        log.info('start {}, end {} , interactions[start:end] {}'.format(start, end, interactions[start:end]))
                    average_contacts_[i] = np.mean(interactions[start:end])
                # handle border conditions
                for i in range(window_size):
                    start = i - window_size_upstream
                    if start < 0:
                        start = 0
                    end = i + window_size + 1
                    # log.info('start {}, end {}'.format(start, end))
                    # log.info('interactions[start:end]: {}'.format(interactions[start:end]))

                    # log.info('interactions[-end:]: {}'.format(interactions[-end:]))

                    average_contacts_[i] = np.mean(interactions[start:end])
                    average_contacts_[-(i + 1)] = np.mean(interactions[-end:])
                average_contacts.append(average_contacts_)
                log.info(average_contacts_)

            for i, interactions_list_ in enumerate(interactions_list):
                with open(args.interactionOutFileName + '_' + matrix_name_legend[i] + '.bedgraph', 'w') as fh:
                    for j, interaction in enumerate(interactions_list_):
                        fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{:.12f}\t{:.12f}\n".
                                 format(toString(interaction[0]), toString(interaction[1]), toString(interaction[2]),
                                        toString(interaction[3]), toString(interaction[4]),
                                        toString(interaction[5]), float(interaction[6]),
                                        float(average_contacts[i][j])))

                        # print("{}\t{}\t{}\t{}\t{}\t{}\t{:.12f}\t{:.12f}\n".\
                        # format(toString(interaction[0]), toString(interaction[1]), toString(interaction[2]),
                        #         toString(interaction[3]), toString(interaction[4]),
                        #         toString(interaction[5]), float(interaction[6]),
                        #         float(average_contacts[i][j])))
                return

        for i, interactions_list_ in enumerate(interactions_list):
            print(matrix_name_legend[i])
            with open(args.interactionOutFileName + '_' + matrix_name_legend[i] + '.bedgraph', 'w') as fh:
                for interaction in interactions_list_:
                    fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{:.12f}\n".
                             format(toString(interaction[0]), toString(interaction[1]),
                                    toString(interaction[2]), toString(interaction[3]),
                                    toString(interaction[4]), toString(interaction[5]),
                                    float(interaction[6])))
                    # print("{}\t{}\t{}\t{}\t{}\t{}\t{:.12f}\n".\
                    # format(toString(interaction[0]), toString(interaction[1]),
                    #                 toString(interaction[2]), toString(interaction[3]),
                    #                 toString(interaction[4]), toString(interaction[5]),
                    #                 float(interaction[6])))
