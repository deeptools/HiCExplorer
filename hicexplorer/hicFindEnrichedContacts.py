from __future__ import division
import argparse
from hicexplorer import HiCMatrix
import hicexplorer.parserCommon
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """
    parse arguments
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Identifies enriched contacts by computing an observe vs. expected or a z-score matrix.',
        add_help=False)

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='Name of the Hi-C matrix in .h5 format.',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting matrix.',
                                type=hicexplorer.parserCommon.writableFile,
                                required=True)

    parserRequired.add_argument(
        '--method',
        help='Method to transform the matrix values.',
        choices=['z-score', 'obs/exp'],
        required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument(
        '--perchr',
        help='Default is to fit distributions per each distance. Setting this '
        'option will fit distributions per distance per chromosome.',
        action='store_true')

    parserOpt.add_argument(
        '--skipDiagonal', '-s',
        help='If set, diagonal counts are not included.',
        action='store_true')

    parserOpt.add_argument(
        '--depth',
        help='Depth (in base pairs) up to which the computations will be carried out. A depth of 10.0000 bp '
             'means that any computations involving bins that are over 10kbp apart are not considered.',
        type=int,
        default=None)

    parserOpt.add_argument("-h", "--help", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):
    log.debug(args)
    args = parse_arguments().parse_args(args)

    hic_ma = HiCMatrix.hiCMatrix(args.matrix)
    try:
        hic_ma.maskBins(hic_ma.nan_bins)
    except AttributeError:
        pass

    if args.skipDiagonal:
        hic_ma.diagflat()

    if args.method == 'obs/exp':
        hic_ma.convert_to_obs_exp_matrix(maxdepth=args.depth, perchr=args.perchr)
    else:
        hic_ma.convert_to_zscore_matrix(maxdepth=args.depth, perchr=args.perchr)

    hic_ma.save(args.outFileName)
