from __future__ import division
import argparse
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import toString

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
                    Build a background model for viewpoints.
                    An example usage is:
                    $ hicInfo -m matrix1.h5 matrix2.h5 matrix3.h5
                    """
    )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='The input matrices to build the background model on.',
                                nargs='+',
                                required=True)

    parserRequired.add_argument('--viewpoints', '-b',
                                help='Bed file contains all viewpoints which should be used to build the background model.',
                                nargs='+',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def main():

    args = parse_arguments().parse_args()

    viewpoints = []

    with open(args.viewpoints, 'r') as file:
        for line in fh.readlines():
            line.split('\t')
            
        # parser.read_file(file_h)
    background_model_data = []
    for matrix in args.matrices:

        hic_ma = hm.hiCMatrix(matrix)
