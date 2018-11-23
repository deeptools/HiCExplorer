from __future__ import division

import sys
import argparse
import numpy as np
import os
from unidecode import unidecode
from hicmatrix import HiCMatrix as hm
from hicexplorer.utilities import convertInfsToZeros_ArrayFloat
import warnings

def parse_arguments():
    parser = argparse.ArgumentParser(description = "Making the expected matrix")
    #Required
    parser.add_argument("-m",
                        dest = "count_matrix",
                        required = True)
    parser.add_argument("-o",
                        dest = "output_matrix",
                        required = True)
    return parser
def expected_interactions(matrix):
    """
        Computes the expected number of interactions per distance
    """

    expected_interactions = np.zeros(matrix.shape[0])
    row, col = matrix.nonzero()
    distance = np.absolute(row - col)
    occurences = np.arange(matrix.shape[0], 0, -1)
    
    for i, distance_ in enumerate(distance):
        expected_interactions[distance_] += matrix.data[i]
    expected_interactions /= occurences

    mask = np.isnan(expected_interactions)
    expected_interactions[mask] = 0
    mask = np.isinf(expected_interactions)
    expected_interactions[mask] = 0
    return expected_interactions


def _expected_matrix(matrix):
    expected_interactions_in_distance_ = expected_interactions(matrix)
    row, col = matrix.nonzero()
    distance = np.absolute(row - col).astype(np.int32)

    if len(matrix.data) > 0:
        data_type = type(matrix.data[0])
        expected = expected_interactions_in_distance_[distance]
        matrix.data = matrix.data.astype(np.float32)
        matrix.data /= expected
        matrix.data = convertInfsToZeros_ArrayFloat(matrix.data).astype(data_type)
    return matrix


def main():
    parser = parse_arguments()
    args = parser.parse_args()
    ma = hm.hiCMatrix(args.count_matrix)
    expected_matrix = _expected_matrix(ma.matrix)
    ma.setMatrixValues(expected_matrix)
    ma.save(args.output_matrix, pApplyCorrection=False)
