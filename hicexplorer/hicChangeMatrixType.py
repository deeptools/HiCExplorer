from __future__ import division
import argparse

from scipy.sparse import csr_matrix, lil_matrix
import numpy as np

from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import exp_obs_matrix_lieberman
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros


import logging

logging.basicConfig()
log = logging.getLogger("hicChangeMatrixType")
log.setLevel(logging.WARN)


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
    group.add_argument('--pearson_correlated', '-pc', action='store_true')
    group.add_argument('--covariance', '-cov', action='store_true')

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))
    parser.add_argument('--chromosomes',
                        help='List of chromosomes to be included in the '
                        'correlation.',
                        default=None,
                        nargs='+')
    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)

    hic_ma = hm.hiCMatrix(matrixFile=args.matrix)

    if args.chromosomes:
        hic_ma.keepOnlyTheseChr(args.chromosomes)

    length_chromosome = 0
    chromosome_count = len(hic_ma.getChrNames())
    for chrname in hic_ma.getChrNames():
        chr_range = hic_ma.getChrBinRange(chrname)
        length_chromosome += chr_range[1] - chr_range[0]

    if args.obs_exp:
        trasf_matrix = lil_matrix(hic_ma.matrix.shape)

        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
            exp_obs_matrix_ = exp_obs_matrix_lieberman(submatrix, length_chromosome, chromosome_count)
            exp_obs_matrix_ = convertNansToZeros(ecsr_matrix(exp_obs_matrix_))
            exp_obs_matrix_ = convertInfsToZeros(csr_matrix(exp_obs_matrix_)).todense()

            trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = exp_obs_matrix_.tolil()

    elif args.pearson_correlated:
        trasf_matrix = lil_matrix(hic_ma.matrix.shape)

        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]

            pearson_correlation_matrix = np.corrcoef(submatrix.todense())
            pearson_correlation_matrix = convertNansToZeros(csr_matrix(pearson_correlation_matrix))
            pearson_correlation_matrix = convertInfsToZeros(csr_matrix(pearson_correlation_matrix)).todense() 

            trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(pearson_correlation_matrix)
    elif args.covariance:
        trasf_matrix = lil_matrix(hic_ma.matrix.shape)
        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
            corrmatrix = np.cov(submatrix.todense())
            trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(corrmatrix)

    hic_ma.setMatrix(trasf_matrix.tocsr(), cut_intervals=hic_ma.cut_intervals)

    hic_ma.save(args.outFileName, pSymmetric=False)
