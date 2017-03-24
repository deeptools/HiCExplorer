import argparse
import numpy as np
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=('Takes two matrices, normalizes them and applies'
                                                  'the given operation. To normalize the matrices '
                                                  'each element is divided by sum of the matrix.'))

    parser.add_argument('--matrices', '-m',
                        help='matrices to use.',
                        metavar='.h5 file format',
                        nargs=2,
                        required=True)

    parser.add_argument('--outFileName', '-o',
                        help='File name to save the resulting matrix. The output is '
                             'also a .h5 file.',
                        required=True)

    parser.add_argument('--operation',
                        help='Operation to apply for the matrices. Options are: diff, ratio, log2ratio',
                        required=True)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def main():

    args = parse_arguments().parse_args()
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
    else:
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
