from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import toString
from hicmatrix.HiCMatrix import check_cooler
import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
Prints information about a matrix or matrices including matrix size,
number of elements, sum of elements, etc.
An example usage is:
$ hicInfo -m matrix1.h5 matrix2.h5 matrix3.h5
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='The matrix (or multiple matrices) to get information about. '
                                'HiCExplorer supports the following file formats: h5 (native HiCExplorer format) '
                                'and cool.',
                                nargs='+',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--outFileName', '-o',
                           help='File name to save information of the matrix instead of writing it to the bash.'
                           )
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def main():

    args = parse_arguments().parse_args()
    for matrix in args.matrices:
        # if
        hic_ma = hm.hiCMatrix(matrix)
        size = hic_ma.matrix.shape[0]
        num_non_zero = hic_ma.matrix.nnz
        sum_elements = hic_ma.matrix.sum() / 2
        bin_length = hic_ma.getBinSize()
        num_nan_bins = len(hic_ma.nan_bins)
        min_non_zero = hic_ma.matrix.data.min()
        max_non_zero = hic_ma.matrix.data.max()

        chromosomes = list(hic_ma.chrBinBoundaries)

        if args.outFileName:
            with open(args.outFileName, 'w') as file:
                file.write("# Matrix information file. Created with HiCExplorer's hicInfo version {}\n".format(__version__))
                file.write("File:\t{}\n".format(matrix))
                file.write("Size:\t{:,}\n".format(size))
                file.write("Sum:\t{:,}\n".format(sum_elements))
                file.write("Bin_length:\t{}\n".format(bin_length))
                file.write("Chromosomes:\t{}\n".format(", ".join(toString(chromosomes))))
                file.write("Non-zero elements:\t{:,}\n".format(num_non_zero))
                file.write("Minimum (non zero):\t{}\n".format(min_non_zero))
                file.write("Maximum:\t{}\n".format(max_non_zero))
                file.write("NaN bins:\t{}\n".format(num_nan_bins))
                if check_cooler(matrix):
                    file.write('The following columns are available: {}'.format(hic_ma.getInformationCoolerBinNames()))
        else:
            print("File:\t{}".format(matrix))
            print("Size:\t{:,}".format(size))
            print("Sum:\t{:,}".format(sum_elements))
            print("Bin_length:\t{}".format(bin_length))
            print("Chromosomes:\t{}".format(", ".join(toString(chromosomes))))
            print("Non-zero elements:\t{:,}".format(num_non_zero))
            print("Minimum (non zero):\t{}".format(min_non_zero))
            print("Maximum:\t{}".format(max_non_zero))
            print("NaN bins:\t{}".format(num_nan_bins))
            if check_cooler(matrix):
                print('The following columns are available: {}'.format(hic_ma.getInformationCoolerBinNames()))
