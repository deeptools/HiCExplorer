from __future__ import division
import argparse
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Prints information about a matrix or matrices including matrix size,
number of elements, sum of elements, etc.

An example usage is:

$ hicInfo -m matrix1.h5 matrix2.h5 matrix3.h5

""")

    parser.add_argument('--matrices', '-m',
                        help='The matrix (or multiple matrices) to get information about. '
                             'HiCExplorer supports the following file formats: h5 (native HiCExplorer format), '
                             'npz (format used by earlier versions of HiCExplorer), '
                             'dekker (matrix format used in Job Dekker publications), '
                             'and lieberman (format used by Erez Lieberman Aiden). This last formats may change '
                             'in the future.',
                        nargs='+',
                        required=True)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def main():

    args = parse_arguments().parse_args()
    for matrix in args.matrices:

        hic_ma = hm.hiCMatrix(matrix)
        size = hic_ma.matrix.shape[0]
        num_non_zero = hic_ma.matrix.nnz
        sum_elements = hic_ma.matrix.sum() / 2
        bin_length = hic_ma.getBinSize()
        num_nan_bins = len(hic_ma.nan_bins)
        min_non_zero = hic_ma.matrix.data.min()
        max_non_zero = hic_ma.matrix.data.max()
        if not matrix.endswith("lieberman"):
            log.debug("lieberman matrix")
            chromosomes = list(hic_ma.chrBinBoundaries)

        log.info("File:\t{}".format(matrix))
        log.info("Size:\t{:,}".format(size))
        log.info("Sum:\t{:,}".format(sum_elements))
        log.info("Bin_length:\t{}".format(bin_length))
        log.info("Chromosomes:\t{}".format(", ".join(chromosomes)))
        log.info("Non-zero elements:\t{:,}".format(num_non_zero))
        log.info("Minimum (non zero):\t{}".format(min_non_zero))
        log.info("Maximum:\t{}".format(max_non_zero))
        log.info("NaN bins:\t{}".format(num_nan_bins))
        log.info("")
