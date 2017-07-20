from __future__ import division
import argparse
from hicexplorer import HiCMatrix
import hicexplorer.parserCommon
from hicexplorer._version import __version__


def parse_arguments(args=None):
    """
    parse arguments
    """

    parent_parser = hicexplorer.parserCommon.getParentArgParse()
    parser = argparse.ArgumentParser(
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Identifies enriched contacts by computing a observe vs. expected or a z-score matrix')

    parser.add_argument('--outFileName', '-o',
                        help='File name to save the resulting matrix',
                        type=hicexplorer.parserCommon.writableFile,
                        required=True)

    parser.add_argument(
        '--method',
        help='Method to transform the matrix values',
        choices=['z-score', 'obs/exp'],
        required=True)

    parser.add_argument(
        '--perchr',
        help='Default is to fit distributions per each distance. Setting this '
        'option will fit distributions per distance per chromosome',
        action='store_true')

    parser.add_argument(
        '--skipDiagonal', '-s',
        help='If set, diagonal counts are not included',
        action='store_true')

    parser.add_argument(
        '--depth',
        help='Depth (in base pairs) up to which the computations will be carried out. A depth of 10.0000 bp '
             'means that any computations involving bins that are over 10kbp apart are not considered.',
        type=int,
        default=None)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):
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
