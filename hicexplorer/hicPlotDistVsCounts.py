#!/usr/bin/env python
#-*- coding: utf-8 -*-
from __future__ import division

import numpy as np
import argparse
from matplotlib import use as mplt_use
mplt_use('Agg')
import matplotlib.pyplot as plt

import hicexplorer.HiCMatrix as HiCMatrix
from hicexplorer.utilities import remove_outliers


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description='This program makes a distance vs. hi-c counts plot per chromosome.')

    # define the arguments
    parser.add_argument('--matrix', '-m',
                        help='HiC normalized (corrected) matrix',
                        required=True)

    parser.add_argument('--skipDiagonal', '-s',
                        help='If set, diagonal counts are not included',
                        action='store_true')

    parser.add_argument('--plotFile', '-o',
                       help='File name to save the file. The given file '
                        'ending will be used '
                        'to determine the image format. '
                        'The available options are: .png, .emf, '
                        '.eps, .pdf and .svg.',
                       type=argparse.FileType('w'),
                       metavar='file name',
                       required=True)

    parser.add_argument('--chromosomeExclude',
                        help='Exclude the given list of chrosomes from ',
                        nargs='+')

    parser.add_argument('--plotTitle', '-t',
                       help='Plot title')

    return parser


def chr_diagonals(matrix_file_name, chromosome_exclude):
    """
    extract the counts by distance per chromosome
    """
    hic_ma = HiCMatrix.hiCMatrix(matrix_file_name)
    print "removing unwanted chrs"
    hic_ma.filterUnwantedChr()
    if chromosome_exclude is None:
        chromosome_exclude = []

    chrtokeep = [x for x in hic_ma.interval_trees.keys()
                 if x not in chromosome_exclude]
    print "Number of contacts {}".format(hic_ma.matrix.sum())
    hic_ma.keepOnlyTheseChr(chrtokeep)
    diagonal_dict = hic_ma.getCountsByDistance(per_chr=True)

    common_dist = None
    max_dist = 0
    chrom_list = hic_ma.getChrNames()
    for chrom in chrom_list:
        try:
            distances = diagonal_dict[chrom].keys()
            distances[0]
        except (KeyError, IndexError):
            continue
        # get list of common distances
        if max(distances) > max_dist:
            max_dist = max(distances)
        if common_dist is None:
            common_dist = set(distances)
        else:
            common_dist = common_dist.intersection(distances)

    return diagonal_dict, chrom_list, list(common_dist), max_dist


def main(args=None):
    """
    for each distance, compare the
    distribution of two samples,
    report number of cases were they differ
    """
    args = parse_arguments().parse_args(args)

    diagonals_dict, chrom_list, common_dist, max_dist = \
        chr_diagonals(args.matrix, args.chromosomeExclude)

    fig = plt.figure(figsize=(5, 4))
    max_mean = 0
    min_mean = 1e6
    dist_list = np.sort(common_dist)
    for chrom in chrom_list:
        values = diagonals_dict[chrom]
        mean_values = [remove_outliers(values[dist]).mean()
                       for dist in dist_list]
        prc_values_90 = [np.percentile(
                remove_outliers(values[dist]), 90)
                         for dist in dist_list]
        prc_values_10 = [np.percentile(
                remove_outliers(values[dist]), 10)
                         for dist in dist_list]
        if max(prc_values_90) > max_mean:
            max_mean = max(prc_values_90)
        if min(mean_values) < min_mean:
            min_mean = min(mean_values)
        label = chrom
        plt.plot(dist_list, mean_values, label=label)
    plt.yscale('log')
    plt.xscale('log')
    plt.legend(prop={'size':'small'})
    plt.suptitle(args.plotTitle)
    plt.xlabel('genomic distance')
    plt.ylabel('corrected Hi-C counts')
    plt.xlim(min(dist_list), 1e6)
    plt.ylim(1, max_mean + max_mean * 0.1)
    plt.savefig(args.plotFile.name, bbox_inches='tight')
