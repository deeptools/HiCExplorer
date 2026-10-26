"""hicPlotSVL: a thin layer over the C++ executable hicPlotSVL.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from io import StringIO
from multiprocessing import Process, Queue
import time
import numpy as np
import cooler
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
from scipy.stats import ranksums

from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import toString
from hicmatrix.HiCMatrix import check_cooler
import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """The argparse parser of hicPlotSVL, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicPlotSVL')


def compute_relation_short_long_range(pHiCMatrix, pChromosomes, pDistance, pIsCooler, pQueue):
    # hic_matrix = hm.hiCMatrix(
    #             pMatrixFile=args.matrix, pChrnameList=[args.region])
    svl_relations = []
    sum_smaller = []
    sum_greater = []
    for chromosome in pChromosomes:
        if pIsCooler:
            hic_matrix_obj = hm.hiCMatrix(
                pMatrixFile=pHiCMatrix, pChrnameList=[chromosome])
            max_distance = pDistance / hic_matrix_obj.getBinSize()
            hic_matrix = hic_matrix_obj.matrix

        else:
            # in case it is h5 pHiCMatrix is an HiCMatrix object
            indices_chromosome = pHiCMatrix.getChrBinRange(chromosome)
            hic_matrix = pHiCMatrix.matrix[indices_chromosome[0]:indices_chromosome[1], indices_chromosome[0]:indices_chromosome[1]]
            max_distance = pDistance / pHiCMatrix.getBinSize()

        instances, features = hic_matrix.nonzero()
        distances = np.absolute(instances - features)
        mask = distances <= max_distance

        sum_smaller_max_distance = np.sum(hic_matrix.data[mask])
        sum_greater_max_distance = np.sum(hic_matrix.data[~mask])
        svl_relation = sum_smaller_max_distance / sum_greater_max_distance
        if np.isinf(svl_relation) or np.isnan(svl_relation):
            continue
        svl_relations.append(svl_relation)
        sum_smaller.append(sum_smaller_max_distance)
        sum_greater.append(sum_greater_max_distance)

    pQueue.put([svl_relations, sum_smaller, sum_greater])
    return


def main(args=None):
    """Run the C++ implementation of hicPlotSVL with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicPlotSVL')(args)
