from __future__ import division
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
import numpy as np

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Computes the log2 ratio between two matrices. The larger matrix is scaled down to match the '
                    'total sum of the smaller matrix. ')

    # define the arguments
    parser.add_argument('--treatment', '-t',
                        help='Treatment Hi-C matrix',
                        required=True)

    parser.add_argument('--control', '-c',
                        help='Control Hi-C matrix',
                        required=True)

    parser.add_argument('--outFileName', '-o',
                        help='File name to save the resulting matrix',
                        required=True)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)

    hic_t = hm.hiCMatrix(matrixFile=args.treatment)
    hic_c = hm.hiCMatrix(matrixFile=args.control)

    # scale larger matrix down
    total_t = hic_t.matrix.sum()
    total_c = hic_c.matrix.sum()

    if total_c > total_t:
        scale_factor = [1, float(total_t) / total_c]
    else:
        scale_factor = [float(total_c) / total_t, 1]

    hic_t.matrix.data = hic_t.matrix.data * scale_factor[0]
    hic_c.matrix.data = hic_c.matrix.data * scale_factor[1]

    """
    Uses sparse matrix tricks to convert
    into a vector the matrix values such
    that zero values that appear in only
    one of the matrices is kept. But
    zeros in two matrices are removed

    Requires two sparse matrices as input
    """
    assert hic_t.matrix.shape == hic_c.matrix.shape, log.error("Matrices have different shapes.")

    assert (hic_t.matrix - hic_c.matrix).sum() != 0, log.error("Matrices are identical.")

    # create a new matrix that is the sum of the two
    # matrices to compare. The goal is to have
    # a matrix that contains all the positions
    # that are non-zero in both matrices
    _mat = hic_t.matrix + hic_c.matrix

    # add one to each element in the new matrix
    _mat.data += 1

    # get a vector of the values in hic_t from
    # _mat
    values_t = (_mat - hic_t.matrix).data - 1

    # get a vector of the values in hic_c from
    # _mat
    values_c = (_mat - hic_c.matrix).data - 1

    # compute log2ratio
    _mat.data = np.log2(values_t / values_c)

    hic_t.matrix = _mat
    hic_t.matrix.eliminate_zeros()
    hic_t.save(args.outFileName)
