from __future__ import division
import argparse
from os.path import basename, dirname

from scipy.sparse import csr_matrix, lil_matrix
import numpy as np

from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import exp_obs_matrix_lieberman, exp_obs_matrix_norm
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
        description='Converts the (interaction) matrix to a observed/expected matrix or a pearson_correlated.')

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
                           help='Transformation method to use. If the option all is used, all three matrices in '
                           'consecutively way (input -> obs_exp -> pearson -> covariance) are created.',
                           choices=['obs_exp', 'pearson', 'covariance', 'all', 'norm'],
                           default='obs_exp')

    parserOpt.add_argument('--chromosomes',
                           help='List of chromosomes to be included in the '
                           'correlation.',
                           default=None,
                           nargs='+')
   
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads for pearson correlation.',
                           required=False,
                           default=4,
                           type=int)

    parserOpt.add_argument("-help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def __obs_exp(pSubmatrix, pLengthChromosome, pChromosomeCount):

    exp_obs_matrix_ = exp_obs_matrix_lieberman(pSubmatrix, pLengthChromosome, pChromosomeCount)
    exp_obs_matrix_ = convertNansToZeros(csr_matrix(exp_obs_matrix_))
    exp_obs_matrix_ = convertInfsToZeros(csr_matrix(exp_obs_matrix_)).todense()
    return exp_obs_matrix_


def __pearson(pSubmatrix):
    pearson_correlation_matrix = np.corrcoef(pSubmatrix)
    pearson_correlation_matrix = convertNansToZeros(csr_matrix(pearson_correlation_matrix))
    pearson_correlation_matrix = convertInfsToZeros(csr_matrix(pearson_correlation_matrix)).todense()
    return pearson_correlation_matrix

def _obs_exp_norm(pSubmatrix, pLengthChromosome, pChromosomeCount):
    
    exp_obs_matrix_ = exp_obs_matrix_norm(pSubmatrix, pLengthChromosome, pChromosomeCount)
    exp_obs_matrix_ = convertNansToZeros(csr_matrix(exp_obs_matrix_))
    exp_obs_matrix_ = convertInfsToZeros(csr_matrix(exp_obs_matrix_)).todense()
    return exp_obs_matrix_

def main(args=None):

    args = parse_arguments().parse_args(args)

    if not args.outFileName.endswith('.h5') or args.outFileName.endswith('.cool'):
        log.error('Output filetype not known.')
        log.error('It is: {}'.format(args.outFileName))
        log.error('Accepted is .h5 or .cool')
        exit(1)

    hic_ma = hm.hiCMatrix(pMatrixFile=args.matrix)
    log.info("hic_ma.matrix: {}".format(hic_ma.matrix))
    if args.chromosomes:
        hic_ma.keepOnlyTheseChr(args.chromosomes)

    length_chromosome = 0
    chromosome_count = len(hic_ma.getChrNames())
    for chrname in hic_ma.getChrNames():
        chr_range = hic_ma.getChrBinRange(chrname)
        length_chromosome += chr_range[1] - chr_range[0]
    trasf_matrix = lil_matrix(hic_ma.matrix.shape)

    if args.method == 'norm':
        trasf_matrix = lil_matrix(hic_ma.matrix.shape)
        # trasf_matrix_pearson = lil_matrix(hic_ma.matrix.shape)
        # trasf_matrix_corr = lil_matrix(hic_ma.matrix.shape)

        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]

            submatrix.astype(float)
            submatrix = _obs_exp_norm(submatrix, length_chromosome, chromosome_count)

            submatrix = __pearson(submatrix)
            trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(submatrix)

        # hic_ma.setMatrix(trasf_matrix.tocsr(), cut_intervals=hic_ma.cut_intervals)
        # hic_ma.save('obs_norm_pearson.'+ args.outFileName, pSymmetric=False, pApplyCorrection=False)

    elif args.method == 'obs_exp':
        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
            submatrix.astype(float)
            trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(__obs_exp(submatrix, length_chromosome, chromosome_count))

    elif args.method == 'pearson':
        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
            log.debug("shape: {}".format(submatrix.shape))

            submatrix.astype(float)
            log.debug("shape: {}".format(submatrix.shape))

            trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(__pearson(submatrix.todense()))

    elif args.method == 'covariance':
        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
            log.debug("shape: {}".format(submatrix.shape))

            submatrix.astype(float)
            log.debug("shape: {}".format(submatrix.shape))

            corrmatrix = np.cov(submatrix.todense())
            trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(corrmatrix)

    elif args.method == 'all':
        trasf_matrix_obs_exp = lil_matrix(hic_ma.matrix.shape)
        trasf_matrix_pearson = lil_matrix(hic_ma.matrix.shape)
        trasf_matrix_corr = lil_matrix(hic_ma.matrix.shape)

        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]

            submatrix.astype(float)
            submatrix = __obs_exp(submatrix, length_chromosome, chromosome_count)

            trasf_matrix_obs_exp[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(submatrix)
            submatrix = __pearson(submatrix)

            trasf_matrix_pearson[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(submatrix)
            corrmatrix = np.cov(submatrix)
            trasf_matrix_corr[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(corrmatrix)

        hic_ma.setMatrix(trasf_matrix_obs_exp.tocsr(), cut_intervals=hic_ma.cut_intervals)

        basename_outFileName = basename(args.outFileName)
        basename_obs_exp = "obs_exp_" + basename_outFileName
        basename_pearson = "pearson_" + basename_outFileName
        basename_covariance = "covariance_" + basename_outFileName
        path = dirname(args.outFileName)
        if path != '':
            path += '/'

        hic_ma.save(path + basename_obs_exp, pSymmetric=False, pApplyCorrection=False)

        hic_ma.setMatrix(trasf_matrix_pearson.tocsr(), cut_intervals=hic_ma.cut_intervals)
        hic_ma.save(path + basename_pearson, pSymmetric=False, pApplyCorrection=False)

        hic_ma.setMatrix(trasf_matrix_corr.tocsr(), cut_intervals=hic_ma.cut_intervals)
        hic_ma.save(path + basename_covariance, pSymmetric=False, pApplyCorrection=False)

    if not args.method == 'all':
        hic_ma.setMatrix(trasf_matrix.tocsr(), cut_intervals=hic_ma.cut_intervals)
        hic_ma.save(args.outFileName, pSymmetric=False, pApplyCorrection=False)
