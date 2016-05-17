#!/usr/bin/env python
#-*- coding: utf-8 -*-

import argparse
import numpy as np
from scipy.sparse import csr_matrix, dia_matrix, coo_matrix
from scipy.sparse import vstack as sparse_vstack
from scipy.sparse import hstack as sparse_hstack
from scipy.sparse import triu, tril

from hicexplorer import HiCMatrix
from hicexplorer._version import __version__

def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Merges multiple hic Matrices of different lengths (eg, from different chromosomes)')

    # define the arguments
    parser.add_argument('--matrixList', '-m',
                        help='List of matrices to merge ',
                        nargs='+',
                        required=True)

    parser.add_argument('--outFileName', '-o',
                        help='output matrix file name ',
                        required=True,
                        default='hicexplorer')

    parser.add_argument('--bplimit', '-b',
                        help='maximum limit (in base pairs) after which the matrix '
                              'will be truncated. i.e. TADs bigger than this size will'
                              'not be shown. For Matrices with very high resolution, '
                              'truncating the matrix after a limit helps in saving memory '
                              'during processing, without much loss of data. You can use '
                              'bplimit of 2 x size of biggest expected TAD. ',
                        type=int,
                        metavar='INT bp',
                        default=100000)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)
    matrixList = args.matrixList
    bplimit = args.bplimit
    filename = args.outFileName

    ## Create empty row, col, value for the matrix
    cut_intervals = []
    row = np.array([])
    col = np.array([])
    values = np.array([])
    ## for each chr, append the row, col, value to the first one. Extend the dim
    size = 0
    nan_bins = []
    for i in range(0, len(matrixList)):
        hicmat = HiCMatrix.hiCMatrix(matrixFile= matrixList[i])
        # trim matrix
        binsize = hicmat.getBinSize()
        limit = int(bplimit / binsize)
        hicmat.matrix = (triu(hicmat.matrix, k=-limit) - triu(hicmat.matrix, k=limit)).tocoo()
        # add data
        row = np.concatenate([row, hicmat.matrix.row + size])
        col = np.concatenate([col, hicmat.matrix.col + size])
        values = np.concatenate([values, hicmat.matrix.data])
        nan_bins = np.concatenate([nan_bins, hicmat.nan_bins + size])

        cut_intervals.append(hicmat.cut_intervals)
        size += hicmat.matrix.shape[0]

    final_mat = coo_matrix((values, (row, col)), shape=(size, size))
    chrNameList, startList, endList, extraList = zip(*cut_intervals[0])

    # save only the upper triangle of the
    # symmetric matrix
    matrix = triu(final_mat, k=0, format='csr')
    try:
        np.savez(
            filename, matrix=matrix, chrNameList=chrNameList,
            startList=startList, endList=endList, extraList=extraList,
            nan_bins=nan_bins)
    except Exception as e:
        print "error saving matrix: {}".format(e)
        try:
            print "Matrix can not be saved because is too big!"
            print "Eliminating entries with only one count."

            # try to remove noise by deleting 1
            matrix.data = matrix.data - 1
            matrix.eliminate_zeros()
            np.savez(
                filename, matrix=matrix, chrNameList=chrNameList,
                startList=startList, endList=endList, extraList=extraList,
                nan_bins=nan_bins)
        except:
            print "Matrix can not be saved because is too big!"
        exit()


if __name__ == '__main__':
    main()
