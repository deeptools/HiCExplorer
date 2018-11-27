from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
import numpy as np
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     add_help=False,
                                     description=('Takes two matrices as input, normalizes them and applies '
                                                  'the given operation. To normalize the matrices '
                                                  'each element is divided by the sum of the matrix.'))

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='Name of the matrices in .h5 format to use, separated by a space.',
                                metavar='matrix.h5',
                                nargs=2,
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting matrix. The output is '
                                'also a .h5 file.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--operation',
                           help='Operation to apply to the matrices.',
                           choices=['diff', 'ratio', 'log2ratio'],
                           default='log2ratio')

    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)
    if args.operation not in ['diff', 'ratio', 'log2ratio']:
        exit("Operation not found. Please use 'diff', 'ratio' or 'log2ratio'.")

    hic1 = hm.hiCMatrix(args.matrices[0])
    hic2 = hm.hiCMatrix(args.matrices[1])

    if hic1.matrix.shape != hic2.matrix.shape:
        exit("The two matrices have different size. Use matrices having the same resolution and created using"
             "the same parameters. Check the matrix values using the tool `hicInfo`.")

    if hic1.chrBinBoundaries != hic2.chrBinBoundaries:
        exit("The two matrices have different chromosome order. Use the tool `hicExport` to change the order.\n"
             "{}: {}\n"
             "{}: {}".format(args.matrices[0], hic1.chrBinBoundaries.keys(),
                             args.matrices[1], hic2.chrBinBoundaries.keys()))

    # normalize by total matrix sum
    hic1.matrix.data = hic1.matrix.data.astype(float) / hic1.matrix.data.sum()
    hic2.matrix.data = hic2.matrix.data.astype(float) / hic2.matrix.data.sum()

    nan_bins = set(hic1.nan_bins)
    nan_bins = nan_bins.union(hic2.nan_bins)

    if args.operation == 'diff':
        new_matrix = hic1.matrix - hic2.matrix
    elif args.operation == 'ratio' or args.operation == 'log2ratio':
        hic2.matrix.data = float(1) / hic2.matrix.data
        new_matrix = hic1.matrix.multiply(hic2.matrix)
        # just in case
        new_matrix.eliminate_zeros()
        if args.operation == 'log2ratio':
            new_matrix.data = np.log2(new_matrix.data)
            new_matrix.eliminate_zeros()

    hic1.setMatrixValues(new_matrix)
    hic1.maskBins(sorted(nan_bins))
    hic1.save(args.outFileName)
