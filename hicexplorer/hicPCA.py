from __future__ import division

import argparse

from scipy.sparse import csr_matrix
from scipy import linalg

import numpy as np
import pyBigWig

from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import exp_obs_matrix_lieberman
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros

import logging
log = logging.getLogger(__name__)


def parse_arguments():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve',
        usage="%(prog)s --matrix hic_matrix -o pca1.bedgraph pca2.bedgraph ",
        description="""
Computes PCA eigenvectors for the HiC matrix.

    $ hicPCA --matrix hic_matrix -o pca1.bedgraph pca2.bedgraph

"""
    )

    parser.add_argument('--matrix', '-m',
                        help='HiCExplorer matrix.',
                        required=True)

    parser.add_argument('--outputFileName', '-o',
                        help='File names for the result of the pca. Number of output file '
                             'must match the number of computed eigenvectors.',
                        nargs='+',
                        default=['pca1', 'pca2'],
                        required=True)

    parser.add_argument('--numberOfEigenvectors', '-noe',
                        help='The number of eigenvectors that the PCA should compute.',
                        default=2,
                        type=int,
                        required=False)

    parser.add_argument('--format', '-f',
                        help='output format. Either bedgraph (default) or bigwig.',
                        choices=['bedgraph', 'bigwig'],
                        default='bedgraph',
                        required=False)

    parser.add_argument('--chromosomes',
                        help='List of chromosomes to be included in the '
                        'correlation.',
                        default=None,
                        nargs='+')

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)
    if int(args.numberOfEigenvectors) != len(args.outputFileName):
        log.error("Number of output file names and number of eigenvectors does not match: {} {}".format(len(args.outputFileName), args.numberOfEigenvectors))
        exit(1)

    ma = hm.hiCMatrix(args.matrix)
    ma.maskBins(ma.nan_bins)

    if args.chromosomes:
        ma.keepOnlyTheseChr(args.chromosomes)

    vecs_list = []
    chrom_list = []
    start_list = []
    end_list = []
    # PCA is computed per chromosome
    length_chromosome = 0
    chromosome_count = len(ma.getChrNames())
    for chrname in ma.getChrNames():
        chr_range = ma.getChrBinRange(chrname)
        length_chromosome += chr_range[1] - chr_range[0]
    for chrname in ma.getChrNames():
        chr_range = ma.getChrBinRange(chrname)
        log.debug("Computing pca for chromosome: {}".format(chrname))

        submatrix = ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]

        exp_obs_matrix_ = exp_obs_matrix_lieberman(submatrix, length_chromosome, chromosome_count)
        exp_obs_matrix_ = convertNansToZeros(csr_matrix(exp_obs_matrix_)).todense()
        exp_obs_matrix_ = convertInfsToZeros(csr_matrix(exp_obs_matrix_)).todense()

        pearson_correlation_matrix = np.corrcoef(exp_obs_matrix_)
        pearson_correlation_matrix = convertNansToZeros(csr_matrix(pearson_correlation_matrix)).todense()
        pearson_correlation_matrix = convertInfsToZeros(csr_matrix(pearson_correlation_matrix)).todense()
        corrmatrix = np.cov(pearson_correlation_matrix)

        evals, eigs = linalg.eig(corrmatrix)
        k = args.numberOfEigenvectors

        chrom, start, end, _ = zip(*ma.cut_intervals[chr_range[0]:chr_range[1]])
        vecs_list += eigs[:, :k].tolist()

        chrom_list += chrom
        start_list += start
        end_list += end

    if args.format == 'bedgraph':
        for idx, outfile in enumerate(args.outputFileName):
            assert(len(vecs_list) == len(chrom_list))

            with open(outfile, 'w') as fh:
                for i, value in enumerate(vecs_list):
                    if len(value) == args.numberOfEigenvectors:
                        if isinstance(value[idx], np.complex):
                            value[idx] = value[idx].real
                        fh.write("{}\t{}\t{}\t{}\n".format(chrom_list[i], start_list[i], end_list[i], value[idx]))
    elif args.format == 'bigwig':
        if not pyBigWig.numpy == 1:
            log.error("ERROR: Your version of pyBigWig is not supporting numpy: {}".format(pyBigWig.__file__))
            exit(1)
        old_chrom = chrom_list[0]
        header = []
        for i, chrom_ in enumerate(chrom_list):
            if old_chrom != chrom_:
                header.append((old_chrom, end_list[i - 1]))
            old_chrom = chrom_

        header.append((chrom_list[-1], end_list[-1]))
        for idx, outfile in enumerate(args.outputFileName):
            assert(len(vecs_list) == len(chrom_list))
            values = []

            bw = pyBigWig.open(outfile, 'w')
            # set big wig header
            bw.addHeader(header)
            # create entry lists
            for i, value in enumerate(vecs_list):
                values.append(value[idx].real)
            # write entries

            bw.addEntries(chrom_list, start_list, ends=end_list, values=values)
            bw.close()
    else:
        log.error("Output format not known: {}".format(args.format))
        exit(1)
