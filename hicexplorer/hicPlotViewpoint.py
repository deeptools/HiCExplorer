import argparse

import numpy as np
import hicexplorer.HiCMatrix as hm
import matplotlib.pyplot as plt
import os


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description='Plots the number of interactions around a given reference point in a region.')

    parser.add_argument('--matrix', '-m',
                        help='path of the Hi-C matrices to plot',
                        required=True)

    parser.add_argument('--region',
                        help='The format is chr:start-end ',
                        required=True)

    parser.add_argument('--outFileName', '-out',
                        help='File name to save the image.',
                        required=True)

    parser.add_argument('--chromosome', '-C',
                        help='Only show results for this chromosome.')

    parser.add_argument('--referencePoint', '-rp', help='Reference point.',
                        required=True)

    parser.add_argument('--dpi',
                        help='Resolution for the image in case the'
                             'ouput is a raster graphics image (e.g png, jpg)',
                        type=int,
                        default=300)
    return parser


def relabelTicks(pTick):
    if pTick < 1e6:
        xlabels = "{:.2f} Kb".format(int(pTick) / 1e3)
    else:
        xlabels = "{:.2f} Mb".format(int(pTick) / 1e6)
    return xlabels


def main(args=None):
    args = parse_arguments().parse_args(args)

    hic = hm.hiCMatrix(args.matrix)
    if args.chromosome:
        hic.keepOnlyTheseChr(args.chromosome)
    if args.region:
        args.region = args.region.translate(None, ",.;|!{}()").replace("-", ":")
        region = args.region.split(":")
        if len(region) != 3:
            print("Region format is invalid {}".format(args.region))
            exit(0)
        chrom, start, end = region[0], int(region[1]), int(region[2])

    args.referencePoint = args.referencePoint.translate(None, ",.;|!{}()").replace("-", ":")
    referencePoint = args.referencePoint.split(":")

    view_point = hic.getRegionBinRange(referencePoint[0], int(referencePoint[1]), int(referencePoint[1]))[0]
    view_point_range = hic.getRegionBinRange(chrom, start, end)
    data_list = []
    for idx in range(view_point_range[0], view_point_range[1], 1):
        data = hic.matrix[view_point, idx]
        data_list.append(data)

    ax = plt.subplot(111)
    ax.plot(range(len(data_list)), data_list)
    ax.set_xticks([0, view_point - view_point_range[0], view_point_range[1] - view_point_range[0]])
    xticklabels = [None] * 3
    xticklabels[0] = relabelTicks((int(referencePoint[1]) - start) * (-1))
    xticklabels[1] = referencePoint[0] + ":" + relabelTicks(int(referencePoint[1]))
    xticklabels[2] = relabelTicks(end - int(referencePoint[1]))
    ax.set_xticklabels(xticklabels)
    ax.set_ylabel('Number of interactions')
    left, width = .45, .5
    bottom, height = .25, .7
    right = left + width
    top = bottom + height
    ax.text(right, top, os.path.basename(args.matrix),
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes,
            color='black', fontsize=8)
    plt.savefig(args.outFileName, dpi=args.dpi)
