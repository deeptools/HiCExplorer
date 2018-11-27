from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     add_help=False,
                                     description=('Adds Hi-C matrices of the same size. Format '
                                                  'has to be hdf5 (.h5) or npz. In order to minimze the '
                                                  'the loss of information, it is recommended to '
                                                  'to sum uncorrected matrices (before hicCorrectMatrix).'))

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='Space-delimited names of the matrices to add. The matrices must have the same shape/size. '
                                'You can verify their size by using `hicInfo`.',
                                metavar='.h5 or cooler file format',
                                nargs='+',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting matrix. The output is '
                                'also a .h5 file. Please, do not add the .h5 suffix.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument("-h", "--help", action="help", help="show this help message and exit")
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)

    hic = hm.hiCMatrix(args.matrices[0])
    summed_matrix = hic.matrix
    nan_bins = set(hic.nan_bins)
    for matrix in args.matrices[1:]:
        hic_to_append = hm.hiCMatrix(matrix)
        if hic.chrBinBoundaries != hic_to_append.chrBinBoundaries:
            log.error("The two matrices have different chromosome order. Use the tool `hicExport` to change the order.\n"
                      "{}: {}\n"
                      "{}: {}".format(args.matrices[0], list(hic.chrBinBoundaries),
                                      matrix, list(hic_to_append.chrBinBoundaries)))
            exit(1)

        try:
            summed_matrix = summed_matrix + hic_to_append.matrix
            if len(hic_to_append.nan_bins):
                nan_bins = nan_bins.union(hic_to_append.nan_bins)
        except Exception:
            log.exception(
                "\nMatrix {} seems to be corrupted or of different shape".format(matrix))
            exit(1)

    # save only the upper triangle of the
    # symmetric matrix
    hic.setMatrixValues(summed_matrix)
    hic.maskBins(sorted(nan_bins))
    hic.save(args.outFileName)
    return
