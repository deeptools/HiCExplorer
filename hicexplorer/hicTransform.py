from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse

from scipy.sparse import csr_matrix, lil_matrix
import numpy as np

from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import obs_exp_matrix_lieberman, obs_exp_matrix_norm, obs_exp_matrix_non_zero, obs_exp_matrix
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros


import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='Converts the (interaction) matrix to different types of obs/exp, pearson or covariance matrix.')

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrix', '-m',
                                help='input file. The computation is done per chromosome.',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the exported matrix.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--method', '-me',
                           help='Transformation method to use for input matrix. Transformation is computed per chromosome.'
                           'obs_exp computes the expected matrix as the sum per genomic distance j divided by maximal possible contacts: sum(diagonal(j) / number of elements in diagonal(j) '
                           'obs_exp_lieberman computes the expected matrix as the sum per genomic distance j divided by the : sum(diagonal(j) / (length of chromosome - j))'
                           'obs_exp_non_zero computes the expected matrix as the sum per genomic distance j divided by sum of non-zero contacts: sum(diagonal(j) / number of non-zero elements in diagonal(j)'
                           'obs_exp_norm computes the expected matrix for exp_i,j: sum(diagonal(i-j)) * sum(row(j)) * sum(row(i)) / sum(matrix)'
                           'pearson computes the Pearson correlation matrix on the input matrix: Pearson_i,j = C_i,j / sqrt(C_i,i * C_j,j) and C is the covariance matrix'
                           'covariance computes the Covariance matrix on the input matrix: Cov_i,j = E[M_i, M_j] - my_i * my_j where M is the input matrix and my the mean.',
                           choices=['obs_exp', 'obs_exp_lieberman', 'obs_exp_non_zero', 'obs_exp_norm', 'pearson', 'covariance'],
                           default='obs_exp')

    parserOpt.add_argument('--chromosomes',
                           help='List of chromosomes to be included in the computation.',
                           default=None,
                           nargs='+')
    parserOpt.add_argument('--perChromosome', '-pc',
                           help='Each chromosome is processed individually, inter-chromosomal interactions are ignored. Option not valid for obs_exp_lieberman.',
                           action='store_true')

    parserOpt.add_argument("-help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def _obs_exp_lieberman(pSubmatrix, pLengthChromosome, pChromosomeCount):

    obs_exp_matrix_ = obs_exp_matrix_lieberman(pSubmatrix, pLengthChromosome, pChromosomeCount)
    obs_exp_matrix_ = convertNansToZeros(csr_matrix(obs_exp_matrix_))
    obs_exp_matrix_ = convertInfsToZeros(csr_matrix(obs_exp_matrix_)).todense()
    return obs_exp_matrix_


def _pearson(pSubmatrix):
    pearson_correlation_matrix = np.corrcoef(pSubmatrix)
    pearson_correlation_matrix = convertNansToZeros(csr_matrix(pearson_correlation_matrix))
    pearson_correlation_matrix = convertInfsToZeros(csr_matrix(pearson_correlation_matrix)).todense()
    return pearson_correlation_matrix


def _obs_exp_norm(pSubmatrix):

    obs_exp_matrix_ = obs_exp_matrix_norm(pSubmatrix)
    obs_exp_matrix_ = convertNansToZeros(csr_matrix(obs_exp_matrix_))
    obs_exp_matrix_ = convertInfsToZeros(csr_matrix(obs_exp_matrix_)).todense()
    return obs_exp_matrix_


def _obs_exp(pSubmatrix):

    obs_exp_matrix_ = obs_exp_matrix(pSubmatrix)
    obs_exp_matrix_ = convertNansToZeros(csr_matrix(obs_exp_matrix_))
    obs_exp_matrix_ = convertInfsToZeros(csr_matrix(obs_exp_matrix_)).todense()
    return obs_exp_matrix_


def _obs_exp_non_zero(pSubmatrix):

    obs_exp_matrix_ = obs_exp_matrix_non_zero(pSubmatrix)
    obs_exp_matrix_ = convertNansToZeros(csr_matrix(obs_exp_matrix_))
    obs_exp_matrix_ = convertInfsToZeros(csr_matrix(obs_exp_matrix_)).todense()
    return obs_exp_matrix_


def main(args=None):

    args = parse_arguments().parse_args(args)

    if not args.outFileName.endswith('.h5') and not args.outFileName.endswith('.cool'):
        log.error('Output filetype not known.')
        log.error('It is: {}'.format(args.outFileName))
        log.error('Accepted is .h5 or .cool')
        exit(1)

    if args.matrix.endswith('cool') and args.chromosomes is not None and len(args.chromosomes) == 1:
        hic_ma = hm.hiCMatrix(pMatrixFile=args.matrix, pChrnameList=args.chromosomes)
    else:
        hic_ma = hm.hiCMatrix(pMatrixFile=args.matrix)
        if args.chromosomes:
            hic_ma.keepOnlyTheseChr(args.chromosomes)

    trasf_matrix = lil_matrix(hic_ma.matrix.shape)

    if args.method == 'obs_exp_norm':
        # trasf_matrix = lil_matrix(hic_ma.matrix.shape)
        if args.perChromosome:
            for chrname in hic_ma.getChrNames():
                chr_range = hic_ma.getChrBinRange(chrname)
                submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
                submatrix.astype(float)
                submatrix = _obs_exp_norm(submatrix)
                trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(submatrix)
        else:
            submatrix = _obs_exp_norm(hic_ma.matrix)
            trasf_matrix = csr_matrix(submatrix)

    elif args.method == 'obs_exp':
        if args.perChromosome:

            for chrname in hic_ma.getChrNames():
                chr_range = hic_ma.getChrBinRange(chrname)
                submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
                submatrix.astype(float)
                trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(_obs_exp(submatrix))
        else:
            submatrix = _obs_exp(hic_ma.matrix)
            trasf_matrix = csr_matrix(submatrix)

    elif args.method == 'obs_exp_non_zero':
        if args.perChromosome:

            for chrname in hic_ma.getChrNames():
                chr_range = hic_ma.getChrBinRange(chrname)
                submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
                submatrix.astype(float)
                trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(_obs_exp_non_zero(submatrix))
        else:
            submatrix = _obs_exp(hic_ma.matrix)
            trasf_matrix = csr_matrix(submatrix)
    elif args.method == 'obs_exp_lieberman':
        length_chromosome = 0
        chromosome_count = len(hic_ma.getChrNames())
        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            length_chromosome += chr_range[1] - chr_range[0]
        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
            submatrix.astype(float)
            trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(_obs_exp_lieberman(submatrix, length_chromosome, chromosome_count))
        trasf_matrix = trasf_matrix.tocsr()
        # log.debug('type: {}'.format(type(trasf_matrix)))
    elif args.method == 'pearson':
        if args.perChromosome:
            for chrname in hic_ma.getChrNames():
                chr_range = hic_ma.getChrBinRange(chrname)
                submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]

                submatrix.astype(float)

                trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(_pearson(submatrix.todense()))
        else:
            trasf_matrix = csr_matrix(_pearson(hic_ma.matrix.todense()))

    elif args.method == 'covariance':
        if args.perChromosome:

            for chrname in hic_ma.getChrNames():
                chr_range = hic_ma.getChrBinRange(chrname)
                submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]

                submatrix.astype(float)

                corrmatrix = np.cov(submatrix.todense())
                trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(corrmatrix)
        else:
            corrmatrix = np.cov(hic_ma.matrix.todense())
            trasf_matrix = csr_matrix(corrmatrix)

    if args.perChromosome:
        hic_ma.setMatrix(trasf_matrix.tocsr(), cut_intervals=hic_ma.cut_intervals)
    else:
        hic_ma.setMatrix(trasf_matrix, cut_intervals=hic_ma.cut_intervals)

    hic_ma.save(args.outFileName, pSymmetric=True, pApplyCorrection=False)
