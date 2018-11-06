from __future__ import division
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import toString
from hicmatrix.HiCMatrix import check_cooler
import logging
log = logging.getLogger(__name__)

import numpy as np

def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
Normalizes given matrices either to 0 - 1 range or to the smallest given read number of all matrices.
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
                                default = 'norm_range',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def main():

    args = parse_arguments().parse_args()
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

            min_value = np.min(hic_matrix.matrix.data)
            max_value = np.max(hic_matrix.matrix.data)
            min_max_difference = np.float64(max_value - min_value)
            # log.debug(type(min_max_difference))
            # log.debug(type(hic_matrix.matrix.data[0]))

            hic_matrix.matrix.data -= min_value
            hic_matrix.matrix.data /= min_max_difference
            output_name_array = args.matrices[i].split('.')
            output_name_array[-2] += '_norm'
            output_name = '.'.join(output_name_array)
            hic_matrix.save(output_name, pApplyCorrection=False)
    elif args.normalize == 'smallest':
        argmin = np.argmin(sum_list)
        # log.debug('argmin {}'.format(argmin))
        # log.debug('sum_list {}'.format(sum_list))
        # log.debug('len sum_list {}'.format(len(sum_list)))


        for i, hic_matrix in enumerate(hic_matrix_list):
            hic_matrix.matrix.data = hic_matrix.matrix.data.astype(np.float32)
            if i != argmin:
                adjust_factor = sum_list[i] / sum_list[argmin]
                # log.debug('sum_list[i] {}'.format(sum_list[i]))
                # log.debug('sum_list[argmin]{}'.format(sum_list[argmin]))
                # log.debug('adjust_factor {}'.format(adjust_factor))

                # log.debug('sum of data BEFORE correction {}'.format(np.sum(hic_matrix.matrix.data)))

                hic_matrix.matrix.data /= adjust_factor

                # log.debug('sum of data after correction {}'.format(np.sum(hic_matrix.matrix.data)))
            output_name_array = args.matrices[i].split('.')
            output_name_array[-2] += '_norm'
            output_name = '.'.join(output_name_array)
            hic_matrix.save(output_name, pApplyCorrection=False)


# TODO: remove normalization factors because they no longer hold true