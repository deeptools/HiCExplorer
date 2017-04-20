import argparse
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=('Prints information about a matrix including size, '
                                                  'number of elements, sum of elements, etc.'))

    parser.add_argument('--matrices', '-m',
                        help='matrices to add. Must have the same shape.',
                        metavar='.h5 of .npz file format',
                        nargs='+',
                        required=True)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def main():

    args = parse_arguments().parse_args()
    for matrix in args.matrices:
        print("File:\t{}".format(matrix))

        hic_ma = hm.hiCMatrix(matrix)
        size = hic_ma.matrix.shape[0]
        num_non_zero = hic_ma.matrix.nnz
        sum_elements = hic_ma.matrix.sum() / 2
        bin_length = hic_ma.getBinSize()
        num_nan_bins = len(hic_ma.nan_bins)
        min_non_zero = hic_ma.matrix.data.min()
        max_non_zero = hic_ma.matrix.data.max()
        chromosomes = hic_ma.chrBinBoundaries.keys()

        print("Size:\t{:,}".format(size))
        print("Sum:\t{:,}".format(sum_elements))
        print("Bin_length:\t{}".format(bin_length))
        print("Chromosomes:\t{}".format(", ".join(chromosomes)))
        print("Non-zero elements:\t{:,}".format(num_non_zero))
        print("Minimum (non zero):\t{}".format(min_non_zero))
        print("Maximum:\t{}".format(max_non_zero))
        print("NaN bins:\t{}".format(num_nan_bins))
        print("")
