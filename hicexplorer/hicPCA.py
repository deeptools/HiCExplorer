from __future__ import division

import sys
import argparse
# from scipy.sparse.linalg import eigs

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
from scipy import cov

import pyBigWig

logging.basicConfig()
log = logging.getLogger("hicPCA")
log.setLevel(logging.WARN)


def parse_arguments():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve',
        usage="%(prog)s --matrix hic_matrix -o pca1.bw pca2.bw ",
        description="""
Computes PCA eigenvectors for the HiC matrix.

    $ hicPCA --matrix hic_matrix -o pca1.bw pca2.bw

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
                        required=True)

    parser.add_argument('--format', '-f',
                        help='output format. Either bigwig (default) or bedgraph.',
                        choices=['bedgraph', 'bigwig'],
                        default='bigwig',
                        required=False)
    parser.add_argument('--chromosomes',
                        help='List of chromosomes to be included in the '
                        'correlation.',
                        default=None,
                        nargs='+')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

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
    if int(args.numberOfEigenvectors) != len(args.outputFileName):
        exit("Number of output file names and number of eigenvectors does not match: {} {}".format(len(args.outputFileName), args.numberOfEigenvectors))

    # normalized contact matrix M* is missing.
    # dividing each entry by the gnome-wide
    # average contact probability for loci at
    # that genomic distance
    ma = hm.hiCMatrix(args.matrix)
    ma.maskBins(ma.nan_bins)

    if args.chromosomes:
        ma.keepOnlyTheseChr(args.chromosomes)

    vecs_list = []
    chrom_list = []
    start_list = []
    end_list = []
    # PCA is computed per chromosome
    length_chromosome_dict = {}
    for chrname in ma.getChrNames():
        chr_range = ma.getChrBinRange(chrname)
        length_chromosome_dict[chrname] = chr_range[1] - chr_range[0]
    for chrname in ma.getChrNames():
        chr_range = ma.getChrBinRange(chrname)

        submatrix = ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]


        exp_obs_matrix_ = exp_obs_matrix(submatrix, length_chromosome_dict)
        exp_obs_matrix_ = convertNansToZeros(csr_matrix(exp_obs_matrix_)).todense()
        pearson_correlation_matrix = np.corrcoef(exp_obs_matrix_)
        pearson_correlation_matrix = convertNansToZeros(csr_matrix(pearson_correlation_matrix)).todense()
        
        # # similar to Lieberman-Aiden 2009
       
        # corrmatrix = convertNansToZeros(csr_matrix(corrmatrix)).todense()
        
        # copymatrix = deepcopy(corrmatrix)
        # for row in range(copymatrix.shape[0]):
        #     # print(corrmatrix[row, :].tolist())
        #     # print("bla", type((corrmatrix[row, :].tolist())) )
        #     # exit()
        #     row_value = float(sum(corrmatrix[row, :].tolist()[0]))
        #     for col in range(copymatrix.shape[1]):
        #         copymatrix[row, col] = float(corrmatrix[row, col]) / (row_value / corrmatrix.shape[0])

        corrmatrix = cov(pearson_correlation_matrix)

        # corrmatrix = np.corrcoef(corrmatrix.todense())

        # corrmatrix = convertNansToZeros(csr_matrix(corrmatrix)).todense()
        evals, eigs = linalg.eig(corrmatrix)
        k = int(args.numberOfEigenvectors)

        # vectors_eigen = [[]] * k
        # for idx in range(k):
        #     for i, value in enumerate(eigs[:, :k]):
        #         vectors_eigen[idx].append(value[idx])
        #     vectors_eigen[idx] = np.gradient(vectors_eigen[idx], 2)
        #     print(len(vectors_eigen[1]))
        #     print(len(eigs[1][idx]))
        #     for i in range(len(vecs_list)):
        #         eigs[i][idx] = vectors_eigen[idx][i]


        chrom, start, end, _ = zip(*ma.cut_intervals[chr_range[0]:chr_range[1]])
        vecs_list += eigs[:, :k].tolist()

        # vecs_list += vectors_eigen[:, :].tolist()
        chrom_list += chrom
        start_list += start
        end_list += end

    # vectors_eigen = [[]] * k
    # for idx in range(k):
    #     for i, value in enumerate(vecs_list):
    #         vectors_eigen[idx].append(value[idx])
    #     vectors_eigen[idx] = np.gradient(vectors_eigen[idx])
    #     for i in range(len(vecs_list)):
    #         vecs_list[i][idx] = vectors_eigen[idx][i]

    if args.format == 'bedgraph':
        for idx, outfile in enumerate(args.outputFileName):
            assert(len(vecs_list) == len(chrom_list))

            with open(outfile, 'w') as fh:
                for i, value in enumerate(vecs_list):
                    if len(value) == int(args.numberOfEigenvectors):
                        fh.write("{}\t{}\t{}\t{}\n".format(chrom_list[i], start_list[i], end_list[i], value[idx]))
    elif args.format == 'bigwig':
        if not pyBigWig.numpy == 1:
            exit("ERROR: Your version of pyBigWig is not supporting numpy: {}".format(pyBigWig.__file__))
        old_chrom = chrom[0]
        header = []
        for i, chrom_ in enumerate(chrom_list):
            if old_chrom != chrom_:
                header.append((old_chrom, end_list[i - 1]))
            old_chrom = chrom_
        header.append((chrom_list[-1], end_list[-1]))

        for idx, outfile in enumerate(args.outputFileName):
            assert(len(vecs_list[:, idx]) == len(chrom_list))
            values = []

            bw = pyBigWig.open(outfile, 'w')
            # set big wig header
            bw.addHeader(header)
            # create entry lists
            for i, value in enumerate(vecs_list[:, idx]):
                values.append(value.real)
            # write entries
            bw.addEntries(list(chrom_list), list(start_list), ends=list(end_list), values=values)
            bw.close()
    else:
        exit("Output format not known: {}".format(args.format))
