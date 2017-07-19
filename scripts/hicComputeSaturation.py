#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from hicexplorer._version import __version__
# import sys
import argparse
import numpy as np
# import pandas as pd
import copy
import pybedtools as bt
import os
import matplotlib
import matplotlib.pyplot as plt

from hicexplorer import HiCMatrix, hicCorrectMatrix, hicFindTADs
from hicexplorer.iterativeCorrection import iterativeCorrection

matplotlib.use('Agg')


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Computes overlap of TADs called at different subsamples of HiCMatrix')

    # define the arguments
    parser.add_argument('--matrixFile', '-m',
                        help='input uncorrected matrix file ',
                        required=True)

    parser.add_argument('--outPrefix', '-o',
                        help='output file name prefix ',
                        default='hicexplorer')

    parser.add_argument('--numberOfProcessors', '-p',
                        help='Number of processors to use for TAD calling',
                        type=int,
                        default=1)

    parser.add_argument('--step',
                        help='Step size for TAD calling',
                        type=int,
                        metavar='INT bp',
                        default=10000)

    parser.add_argument('--minDepth',
                        help='minimum depth for TAD calling',
                        type=int,
                        metavar='INT bp',
                        default=20000)

    parser.add_argument('--maxDepth',
                        help='maximum depth for TAD calling',
                        type=int,
                        metavar='INT bp',
                        default=60000)

    parser.add_argument('--useLogValues',
                        help='use log values for TAD calling',
                        action='store_true')

    parser.add_argument('--overlapResolution',
                        help='Resolution for overlapping TADs. Default resolution'
                        'is the one you used to make the matrix.',
                        type=int,
                        metavar='INT bp'
                        )

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)

    out = args.outPrefix
    hic_matrix = HiCMatrix.hiCMatrix(matrixFile=args.matrixFile)
    if args.overlapResolution is not None:
        resolution = args.overlapResolution
    else:
        resolution = hic_matrix.getBinSize()

    # filtering params
    filterThreshold = list([-3, 3])
    iterNum = 100
    # compute spectra param
    delta = 0.001
    lookahead = 2

    # normalize, call boundaries, and save as bed

    for i in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:

        print("\n Sampling {}% of data \n".format(i))
        # subsample
        ma = copy.deepcopy(hic_matrix)
        ma.matrix.data = ((hic_matrix.matrix.data) *
                          (float(i) / 100)).astype(int)
        ma.matrix.eliminate_zeros()

        # remove outliers
        outlier_regions = hicCorrectMatrix.filter_by_zscore(
            ma, filterThreshold[0], filterThreshold[1])
        pct_outlier = 100 * float(len(outlier_regions)) / ma.matrix.shape[0]
        ma.printchrtoremove(outlier_regions, label="Bins that are MAD outliers after merge ({:.2f}%) "
                                                   "out of".format(pct_outlier, ma.matrix.shape[0]))

        # mask filtered regions
        ma.maskBins(outlier_regions)
        # total_filtered_out = set(outlier_regions)

        # pre_row_sum = np.asarray(ma.matrix.sum(axis=1)).flatten()

        # correct matrix
        correction_factors = []
        corrected_matrix, correction_factors = iterativeCorrection(
            ma.matrix, iterNum)

        ma.setMatrixValues(corrected_matrix)
        ma.setCorrectionFactors(correction_factors)

        # compute Spectra
        chrom, chr_start, chr_end, matrix = hicFindTADs.compute_spectra_matrix(
            args, matrix=ma)

        # findTADs and save
        min_idx = hicFindTADs.find_consensus_minima(
            matrix, lookahead=lookahead, delta=delta)
        boundaries = np.array([chr_start[idx] for idx in min_idx])
        chrom = chrom[min_idx]

        bound_ext = np.array(
            [boundaries - resolution, boundaries + resolution], dtype=int)
        bound_ext = np.transpose(bound_ext)
        outfile = "{}_{}.bed".format(out, i)
        with open(outfile, 'w') as fh:
            for i in range(bound_ext.shape[0]):
                fh.write("{}\t{}\t{}\n".format(
                    chrom[i], bound_ext[i, 0], bound_ext[i, 1]))

    # Intersect all the beds and plot
    full = bt.bedtool.BedTool("{}_100.bed".format(out))

    isect_all = np.array([])
    for num in ['10', '20', '30', '40', '50', '60', '70', '80', '90', '100']:
        bedfile = "{}_{}.bed".format(out, num)
        bed = bt.bedtool.BedTool(bedfile)
        isect = len(bed.intersect(full))
        isect_all = np.append(isect_all, isect)
        os.remove(bedfile)

    isect_all = isect_all * 100 / isect_all[9]
    x = np.arange(10., 110., 10.0)

    plt.figure(figsize=(8, 8), dpi=300)
    plt.plot(x, isect_all)
    plt.xlabel('Sampling (%)')
    plt.ylabel('TADs called (% w.r.t. total)')
    plt.title('Overlap of TADs called per sample')
    plt.savefig(out, dpi=300)


if __name__ == '__main__':
    main()
