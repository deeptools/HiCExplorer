#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import argparse

from hicexplorer import readBed
import hicexplorer.HiCMatrix as hm


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='This script computes the sum of contacts inside TADs, between TADs and the ratio'
                    'inter TAD / intra TAD. Inter-TAD values are computed for neighboring TADs only.')

    parser.add_argument('--matrix', '-m',
                        help='Corrected Hi-C matrix',
                        required=True)

    parser.add_argument('--tads',
                        help='Path to a TADs bed file. Each entry should contain the chrom, start and '
                             'end position of a TAD',
                        type=argparse.FileType('r'),
                        required=True)

    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)
    hicma = hm.hiCMatrix(args.matrix)

    prev_chrom = None
    prev_start = None

    bed_h = readBed.ReadBed(args.tads)
    sum_tads = 0
    sum_inter = 0
    for bed in bed_h:
        if prev_chrom is None or bed.chromosome != prev_chrom:
            start_bin, end_bin = hicma.getRegionBinRange(bed.chromosome, bed.start, bed.end)
            sum_tads += hicma.matrix[start_bin:end_bin, start_bin:end_bin].sum()
            prev_chrom = bed.chromosome
            prev_start = start_bin
            continue

        start_bin, end_bin = hicma.getRegionBinRange(bed.chromosome, bed.start, bed.end)

        sum_inter += hicma.matrix[prev_start:start_bin, start_bin:end_bin].sum()
        sum_tads += hicma.matrix[start_bin:end_bin, start_bin:end_bin].sum()

    print "sum tads\t{}\nsum inter\t{}\nratio inter/tads\t{:.3f}".format(sum_tads, sum_inter, sum_inter / sum_tads)


if __name__ == "__main__":
    main()
