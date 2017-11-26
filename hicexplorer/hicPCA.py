import sys
import argparse
# from scipy.sparse.linalg import eigs

from scipy.linalg import matrix_balance
from scipy.sparse import csr_matrix, lil_matrix
import logging

from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
from utilities import getPearson
from utilities import convertNansToZeros
import matplotlib.pyplot as plt
import numpy as np
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

    if len(args.chromosomes) > 1:
        exit("Only one chromosome or all is supported right now.")
    if args.chromosomes:
        ma.keepOnlyTheseChr(args.chromosomes)

    vecs_list = []
    chrom_list = []
    start_list = []
    end_list = []
    # PCA is computed per chromosome
    for chrname in ma.getChrNames():
        chr_range = ma.getChrBinRange(chrname)

        ma.keepOnlyTheseChr(chrname)
        chr_submatrix = ma.convert_to_obs_exp_matrix()
        # chr_submatrix = ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]

        chr_submatrix = getPearson(chr_submatrix)
        # chr_submatrix = matrix_balance(chr_submatrix)
        # eigenvectors and eigenvalues for the from the matrix
        # chr_submatrix = convertNansToZeros(chr_submatrix).todense()
        vals, vecs = np.linalg.eig(chr_submatrix)#, k=int(args.numberOfEigenvectors), which='LM')
        k = int(args.numberOfEigenvectors)
        vals = vals[:k]
        vecs = vecs[:, :k]
        chrom, start, end, _ = zip(*ma.cut_intervals[chr_range[0]:chr_range[1]])
        vecs_list += vecs.tolist()
        chrom_list += chrom
        start_list += start
        end_list += end

    if args.format == 'bedgraph':
        for idx, outfile in enumerate(args.outputFileName):
            assert(len(vecs_list) == len(chrom_list))

            with open(outfile, 'w') as fh:
                for i, value in enumerate(vecs_list):
                    fh.write("{}\t{}\t{}\t{}\n".format(chrom_list[i], start_list[i], end_list[i], value[idx].real))
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
