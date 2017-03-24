from __future__ import division
import argparse
import numpy as np
import sys
from hicexplorer import HiCMatrix
import hicexplorer.parserCommon
from hicexplorer.utilities import transformMatrix, applyFdr, getPearson
from hicexplorer._version import __version__


def parse_arguments(args=None):
    """
    parse arguments
    """

    parent_parser = hicexplorer.parserCommon.getParentArgParse()
    parser = argparse.ArgumentParser(
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Applies a transformation to a Hi-C matrix.')

    parser.add_argument('--originalMat',
                        help='File name containing the original hic matrix',
                        type=argparse.FileType('r'))

    parser.add_argument('--outFileName', '-o',
                        help='File name to save the resulting matrix ',
                        type=hicexplorer.parserCommon.writableFile,
                        required=True)

    parser.add_argument(
        '--method',
        help='method to transform the matrix values',
        choices=['z-score', 't-score', 'residuals', 'obs/exp', 'p-value',
                 'pearson', 'nbinom-p-value', 'nbinom-expected',
                 'nbinom-est-dist', 'log-norm', 'chi-squared', 'none'],
        required=True)

    parser.add_argument(
        '--applyFdr',
        help='If set the computed values are corrected using fdr. '
        'Should only be used with the methods that produce pvalues.',
        action='store_true')

    parser.add_argument(
        '--perChromosome',
        help='Default is to fit distributions per each distance. Setting this '
        'option will fit distributions per distance per chromosome',
        action='store_true')

    parser.add_argument(
        '--skipDiagonal', '-s',
        help='If set, diagonal counts are not included',
        action='store_true')

    parser.add_argument(
        '--outFormat',
        help='Output format',
        choices=['hdf5', 'dekker'],
        default='hdf5')

    parser.add_argument(
        '--chromosomes',
        help='List of chromosomes to be included in transformation '
        'correction.',
        default=None,
        nargs='+')

    parser.add_argument(
        '--depth',
        help='Depth (in base pairs) up to which the computations will be carried out. A depth of 10.0000 bp '
             'means that any computations involving points that are over 10kbp are not considered.',
        type=int,
        default=None)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def main():
    """
    collects all arguments and executes
    the appropriate functions
    """
    args = parse_arguments().parse_args()

    hic_ma = HiCMatrix.hiCMatrix(args.matrix)
    if args.originalMat:
        orig_ma = HiCMatrix.hiCMatrix(args.originalMat.name)
    else:
        orig_ma = None

    try:
        hic_ma.maskBins(hic_ma.nan_bins)
    except AttributeError:
        pass

    # remove unwanted Chrs or select a given chromosome
    # in case is given
    hic_ma.filterUnwantedChr()
    if args.originalMat:
        orig_ma.filterUnwantedChr()
    if args.chromosomes:
        hic_ma.keepOnlyTheseChr(args.chromosomes)
        if args.originalMat:
            orig_ma.reorderChromosomes(hic_ma.chrBinBoundaries.keys())

    if args.skipDiagonal:
        hic_ma.diagflat()
        if args.originalMat:
            orig_ma.diagflat()

    max_depth_in_bins = None
    if args.depth:
        binsize = hic_ma.getBinSize()
        if args.depth < binsize:
            exit("Please specify a depth larger than bin size ({})".format(binsize))
        max_depth_in_bins = int(args.depth / binsize)
        import scipy.sparse
        # work only with the upper matrix
        # and remove all pixels that are beyond
        # max_depth_in_bis
        # (this is done by subtracting a second sparse matrix
        # that contains only the upper matrix that wants to be removed.
        hic_ma.matrix = scipy.sparse.triu(hic_ma.matrix, k=0, format='csr') - \
            scipy.sparse.triu(hic_ma.matrix, k=max_depth_in_bins, format='csr')
        hic_ma.matrix.eliminate_zeros()

    if args.method == 'obs/exp':
        hic_ma.convert_to_obs_exp_matrix()
        new_ma = hic_ma.matrix
    elif args.method == 'pearson':
        sys.stderr.write("\nComputing observed / expected\n")
        hic_ma.convert_to_obs_exp_matrix()
        sys.stderr.write("\nComputing pearson\n")
        new_ma = getPearson(hic_ma.matrix)
    elif args.method != 'none':
        # check that the normalized and original matrices
        # have the same size
        if orig_ma:
            assert np.all(hic_ma.matrix.shape == orig_ma.matrix.shape), \
                "original and derived matrices do not have same shape"
        new_ma = transformMatrix(hic_ma, args.method,
                                 per_chr=args.perChromosome,
                                 original_matrix=orig_ma,
                                 depth_in_bins=max_depth_in_bins)
    else:
        new_ma = hic_ma.matrix

    if args.applyFdr:
        new_ma = applyFdr(new_ma)

    hic_ma.setMatrixValues(new_ma)
    hic_ma.restoreMaskedBins()
    if args.outFormat == 'dekker':
        hic_ma.save_dekker(args.outFileName)
    else:
        hic_ma.save(args.outFileName)
