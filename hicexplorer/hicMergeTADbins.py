#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import sys
import argparse
from hicmatrix import HiCMatrix as hm
from past.builtins import zip
from builtins import range
import numpy as np
from hicexplorer.utilities import toString
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Uses a BED file of domains or TAD boundaries to merge
the bin counts of a Hi-C matrix.

In this matrix only contains the total counts in each TAD and
its total contacts of all other TADs.""")

    # define the arguments
    parser.add_argument('--matrix', '-m',
                        help='Path to Hi-C matrix to use',
                        required=True)

    parser.add_argument('--domains',
                        help='Path to a bed file containing the domains.',
                        type=argparse.FileType('r'),
                        required=True)

    parser.add_argument('--outFile', '-o',
                        help='Name for the resulting matrix file.',
                        required=True)
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))
    return parser


def merge_tad_bins(hic, boundary_id_list, filename):
    """
    Reduces the HiCMatrix by merging the counts of tad bins.
    :param hic: HiCMatrix object
    :param boundary_id_list list of tad boundary bin ids
    :param filename Name to save the resulting matrix
    :return: HiCMatrix object
    """

    from hicexplorer.reduceMatrix import reduce_matrix
    hic.restoreMaskedBins()
    ref_name_list, start_list, end_list, coverage_list = zip(*hic.cut_intervals)
    new_bins = []
    bins_to_merge = []
    prev_ref = ref_name_list[0]

    # prepare new intervals
    idx_start = 0
    new_start = start_list[0]
    count = 0
    for idx, ref in enumerate(ref_name_list):
        if (count > 0 and idx in boundary_id_list) or ref != prev_ref:
            coverage = np.mean(coverage_list[idx_start:idx])
            new_bins.append((ref_name_list[idx_start], new_start,
                             end_list[idx - 1], coverage))
            bins_to_merge.append(list(range(idx_start, idx)))
            idx_start = idx
            new_start = start_list[idx]
            count = 0

        prev_ref = ref
        count += 1
    # check that the previous for loop ran, otherwise
    # some variables may not be set
    if len(bins_to_merge) > 0:
        coverage = np.mean(coverage_list[idx_start:])
        new_bins.append((ref, new_start, end_list[idx], coverage))
        bins_to_merge.append(list(range(idx_start, idx + 1)))
        # remove correction factors otherwise they are
        # saved but they no longer correspond to the
        # size of the matrix.
        hic.correction_factors = None

        hic.update_matrix(
            reduce_matrix(hic.matrix, bins_to_merge, diagonal=True), new_bins)

        hic.save(filename)
    else:
        log.info("Nothing to merge.")


def get_boundary_bin_id(hic, bed_fh):
    """
    :param hic: HiCMatrix object
    :param bed_fh: file handle of the bed file
    :return: Sorted list of bin indices.
    """
    line_number = 0
    boundaries = set()
    for line in bed_fh.readlines():
        line_number += 1
        line = toString(line)
        if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
            continue
        try:
            chrom, start, end = line.strip().split('\t')[0:3]
        except Exception as detail:
            msg = 'Could not read line\n{}\n. {}'.format(line, detail)
            log.exception(msg)
            sys.exit()

        try:
            start = int(start)
            end = int(end)
        except ValueError as detail:
            msg = "Error reading line: {}. One of the fields is not " \
                  "an integer.\nError message: {}".format(line_number, detail)
            log.exception(msg)
            sys.exit()

        assert start <= end, "Error in line #{}, end1 larger than start1 in {}".format(line_number, line)

        # check the overlap of the region with the hic matrix bins
        start_bin, end_bin = hic.getRegionBinRange(chrom, start, end)
        boundaries.add(start_bin)
        boundaries.add(end_bin)

    return np.sort(list(boundaries))


def main(args=None):

    args = parse_arguments().parse_args(args)
    hic_ma = hm.hiCMatrix(args.matrix)
    hic_ma.restoreMaskedBins()

    # the bin id of boundary positions
    boundary_id_list = get_boundary_bin_id(hic_ma, args.domains)

    # make a reduce matrix by merging the TAD bins
    log.info("Generating matrix with merged bins")
    merge_tad_bins(hic_ma, boundary_id_list, args.outFile)
