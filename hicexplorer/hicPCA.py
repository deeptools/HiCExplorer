import sys, argparse
from scipy.sparse.linalg import eigs
import logging

from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__

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
                        help='File names for the result of the pca. Number of output file must match the number of computed eigenvectors.',
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
                        default = 'bigwig',
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
    ma = hm.hiCMatrix(args.matrix)
    if args.chromosomes:
        ma.keepOnlyTheseChr(args.chromosomes)
    ma.maskBins(ma.nan_bins)
    ma.matrix = ma.matrix.asfptype()
    # eigenvectors and eigenvalues for the from the matrix
    vals, vecs = eigs(ma.matrix, k=int(args.numberOfEigenvectors), which='LM')
    
    # save eigenvectors
    chrom, start, end, _ = zip(*ma.cut_intervals)
    
    if args.format == 'bedgraph':
        for idx, outfile in enumerate(args.outputFileName):
            assert(len(vecs[:, idx]) == len(chrom))
            with open(outfile, 'w') as fh:
                for i, value in enumerate(vecs[:, idx]):
                    fh.write("{}\t{}\t{}\t{}\n".format(chrom[i], start[i], end[i], value.real))
    elif args.format == 'bigwig':
        if not pyBigWig.numpy == 1:
            exit("ERROR: Your version of pyBigWig is not supporting numpy: {}".format(pyBigWig.__file__))
        old_chrom = chrom[0]
        header = []
        for i, chrom_ in enumerate(chrom):
            if old_chrom != chrom_:
                header.append((old_chrom, end[i - 1]))
            old_chrom = chrom_
        header.append((chrom[-1], end[-1]))
       
        for idx, outfile in enumerate(args.outputFileName):
            assert(len(vecs[:, idx]) == len(chrom))
            values = []
            
            bw = pyBigWig.open(outfile, 'w')
            # set big wig header
            bw.addHeader(header)
            # create entry lists
            for i, value in enumerate(vecs[:, idx]):
                values.append(value.real)
            # write entires
            bw.addEntries(list(chrom), list(start), ends=list(end), values=values)
            bw.close()
    else:
        exit("Output format not known: {}".format(args.format))
    