from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
import logging
log = logging.getLogger(__name__)

import numpy as np


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
Normalizes given matrices either to the smallest given read number of all matrices or to 0 - 1 range. However, it does NOT compute the contact probability.
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='The matrix (or multiple matrices) to get information about. '
                                'HiCExplorer supports the following file formats: h5 (native HiCExplorer format) '
                                'and cool.',
                                nargs='+',
                                required=True)

    parserRequired.add_argument('--normalize', '-n',
                                help='Normalize to a) 0 to 1 range, b) all matrices to the lowest read count of the given matrices.',
                                choices=['norm_range', 'smallest'],
                                default='smallest',
                                required=True)
    parserRequired.add_argument('--outFileName', '-o',
                                help='Output file name for the Hi-C matrix.',
                                metavar='FILENAME',
                                nargs='+',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)
    hic_matrix_list = []
    sum_list = []
    for matrix in args.matrices:
        hic_ma = hm.hiCMatrix(matrix)
        if args.normalize == 'smallest':
            sum_list.append(hic_ma.matrix.sum())
        hic_matrix_list.append(hic_ma)

    if args.normalize == 'norm_range':
        for i, hic_matrix in enumerate(hic_matrix_list):
            hic_matrix.matrix.data = hic_matrix.matrix.data.astype(np.float32)
            mask = np.isnan(hic_matrix.matrix.data)
            hic_matrix.matrix.data[mask] = 0

            mask = np.isinf(hic_matrix.matrix.data)
            hic_matrix.matrix.data[mask] = 0
            min_value = np.min(hic_matrix.matrix.data)
            max_value = np.max(hic_matrix.matrix.data)
            min_max_difference = np.float64(max_value - min_value)

            hic_matrix.matrix.data -= min_value
            hic_matrix.matrix.data /= min_max_difference

            mask = np.isnan(hic_matrix.matrix.data)
            hic_matrix.matrix.data[mask] = 0

            mask = np.isinf(hic_matrix.matrix.data)
            hic_matrix.matrix.data[mask] = 0
            hic_matrix.matrix.eliminate_zeros()

            hic_matrix.save(args.outFileName[i], pApplyCorrection=False)
    elif args.normalize == 'smallest':
        argmin = np.argmin(sum_list)

        for i, hic_matrix in enumerate(hic_matrix_list):
            hic_matrix.matrix.data = hic_matrix.matrix.data.astype(np.float32)
            if i != argmin:
                mask = np.isnan(hic_matrix.matrix.data)
                hic_matrix.matrix.data[mask] = 0

                mask = np.isinf(hic_matrix.matrix.data)
                hic_matrix.matrix.data[mask] = 0
                adjust_factor = sum_list[i] / sum_list[argmin]
                hic_matrix.matrix.data /= adjust_factor
                mask = np.isnan(hic_matrix.matrix.data)

            mask = np.isnan(hic_matrix.matrix.data)
            hic_matrix.matrix.data[mask] = 0

            mask = np.isinf(hic_matrix.matrix.data)
            hic_matrix.matrix.data[mask] = 0
            hic_matrix.matrix.eliminate_zeros()

            hic_matrix.save(args.outFileName[i], pApplyCorrection=False)
