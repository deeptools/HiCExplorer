from __future__ import division
import sys
import argparse
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
import numpy as np
import time

from utilities import getPearson
from utilities import convertNansToZeros
from collections import OrderedDict
from six import iteritems
from multiprocessing import Process, Queue
from copy import deepcopy

from scipy.linalg import matrix_balance
from scipy.sparse import csr_matrix, lil_matrix
import logging
from scipy import linalg, dot, cov
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
from utilities import getPearson
from utilities import convertNansToZeros
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

import pyBigWig


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Converts the (interaction) matrix to a observed/expected matrix or a pearson_correlated.')

    # define the arguments
    parser.add_argument('--matrix', '-m',
                        help='input file. The computation is done per chromosome.',
                        required=True)

    parser.add_argument('--outFileName', '-o',
                        help='File name to save the exported matrix.',
                        required=True)
    parser.add_argument('--threads', '-t',
                        help='Number of threads for pearson correlation.',
                        required=False,
                        default=4,
                        type=int)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--obs_exp', '-oe', action='store_true')
    group.add_argument('--pearson_corrected', '-pc', action='store_true')
    group.add_argument('--covariance', '-cov', action='store_true')
    

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))
    parser.add_argument('--chromosomes',
                        help='List of chromosomes to be included in the '
                        'correlation.',
                        default=None,
                        nargs='+')
    return parser

def expected_interactions_in_distance(pLength_chromosome_dict, pCopy_submatrix):
    print("expected_interactions_in_distance...")
    expected_interactions = np.zeros(pCopy_submatrix.shape[0])
    for distance in range(pCopy_submatrix.shape[0]):
        row = 0
        col = distance
        sum_distance = 0.0
        while row < pCopy_submatrix.shape[0] and col < pCopy_submatrix.shape[1]:
            sum_distance += pCopy_submatrix[row, col]
            row += 1
            col += 1
        sum_distance_genome = 0.0
        for element in pLength_chromosome_dict:
            sum_distance_genome += pLength_chromosome_dict[element] - distance
        expected_interactions[distance] = sum_distance / sum_distance_genome
    print("expected_interactions_in_distance...Done")
        
    return expected_interactions

def exp_obs_matrix(pSubmatrix, pLength_chromosome_dict):
    copy_submatrix = deepcopy(pSubmatrix)
    pSubmatrix = pSubmatrix.todense().astype(float)
    expected_interactions_in_distance_ = expected_interactions_in_distance(pLength_chromosome_dict, copy_submatrix)
    for row in range(pSubmatrix.shape[0]):
        for col in range(pSubmatrix.shape[1]):
            distance = abs(row - col)
            pSubmatrix[row, col] = pSubmatrix[row, col] / expected_interactions_in_distance_[distance]
    return pSubmatrix



def main(args=None):
    args = parse_arguments().parse_args(args)

    hic_ma = hm.hiCMatrix(matrixFile=args.matrix)
   
    if args.chromosomes:
        hic_ma.keepOnlyTheseChr(args.chromosomes)

    length_chromosome_dict = {}
    for chrname in hic_ma.getChrNames():
        chr_range = hic_ma.getChrBinRange(chrname)
        length_chromosome_dict[chrname] = chr_range[1] - chr_range[0]
    
    if args.obs_exp:
        trasf_matrix = lil_matrix(hic_ma.matrix.shape)

        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
            exp_obs_matrix_ = exp_obs_matrix(submatrix, length_chromosome_dict)
            exp_obs_matrix_ = convertNansToZeros(csr_matrix(exp_obs_matrix_))
            trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = exp_obs_matrix_.tolil()
            
    elif args.pearson_corrected:
        trasf_matrix = lil_matrix(hic_ma.matrix.shape)

        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
           
            pearson_correlation_matrix = np.corrcoef(submatrix.todense())
            pearson_correlation_matrix = convertNansToZeros(csr_matrix(pearson_correlation_matrix))

            trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(pearson_correlation_matrix)
    elif args.covariance:
        trasf_matrix = lil_matrix(hic_ma.matrix.shape)
        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
            corrmatrix = cov(submatrix.todense())
            trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(corrmatrix)


    hic_ma.setMatrix(trasf_matrix.tocsr(), cut_intervals=hic_ma.cut_intervals)
        


    sys.stderr.write('saving...\n')

    hic_ma.save(args.outFileName, pSymmetric=False)
