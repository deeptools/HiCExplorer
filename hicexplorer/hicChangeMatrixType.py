from __future__ import division
import sys
import argparse
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
import numpy as np

from scipy.sparse import csr_matrix, lil_matrix

from utilities import getPearson
from utilities import convertNansToZeros


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Converts the (interaction) matrix to a observed/expected matrix or a pearson_correlated.')

    # define the arguments
    parser.add_argument('--matrix', '-m',
                        help='input file(s). Could be one or many files. '
                        'Multiple input files are allowed for hicexplorer or lieberman format. '
                        ' In case of multiple input files, they will be combined. ',
                        required=True)

    parser.add_argument('--outFileName', '-o',
                        help='File name to save the exported matrix. In the case of "lieberman" '
                             'output format this should be the path of a folder where the information '
                             'per chromosome is stored.',
                        required=True)

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--obs_exp', '-oe', action='store_true')
    group.add_argument('--pearson_corrected', '-pc', action='store_true')

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))
    parser.add_argument('--chromosomes',
                        help='List of chromosomes to be included in the '
                        'correlation.',
                        default=None,
                        nargs='+')
    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)

    hic_ma = hm.hiCMatrix(matrixFile=args.matrix)
    print(hic_ma.matrix.todense())
    if len(args.chromosomes) > 1:
        exit("Only one chromosome or all is supported right now.")
    if args.chromosomes:
        hic_ma.keepOnlyTheseChr(args.chromosomes)

    if args.obs_exp:
        hic_ma.convert_to_obs_exp_matrix(perchr=True)
        print(hic_ma.matrix.todense())

    elif args.pearson_corrected:
        for chrname in hic_ma.getChrNames():
            pearson_matrix = getPearson(hic_ma.matrix)

            hic_ma.setMatrix(csr_matrix(pearson_matrix), cut_intervals=hic_ma.cut_intervals)

    sys.stderr.write('saving...\n')

    hic_ma.save(args.outFileName, pSymmetric=False)
