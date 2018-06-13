from __future__ import division
import argparse
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import toString
import math
import logging
log = logging.getLogger(__name__)

from .lib import Viewpoint
import numpy as np


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

    parserRequired.add_argument('--referencePoints', '-b',
                                help='Bed file contains all viewpoints which should be used to build the background model.',
                                type=str,
                                required=True,
                                nargs=2)
    parserRequired.add_argument('--range', '-r',
                                help='Region around the viewpoint to consider.',
                                type=int,
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--averageContactBin',
                           help='Average the contacts of n bins, written to last column.',
                           type=int,
                           default=0)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def smooth_viewpoint():
    pass


def main():

    args = parse_arguments().parse_args()

    viewpointObj = Viewpoint()
    referencePoints = viewpointObj.readReferncePointFile(args.referencePoints)
    interactions = []

    # elements_of_viewpoint = args.range * 2 / hic_ma.bin_size
    relative_counts_conditions = []
    relative_positions = set()
    bin_size = 0

    # - compute for each condition (matrix):
    # - all viewpoints and smooth them: sliding window approach
    # - after smoothing, sum all viewpoints up to one
    # - compute the percentage of each position with respect to the total interaction count
    # for models of all conditions:
    # - compute mean percentage for each bin
    # - compute SEM = Standard deviation / (square root of sample size)

    for matrix in args.matrices:
        hic_ma = hm.hiCMatrix(matrix)
        viewpointObj.setHiCMatrixObj(hic_ma)
        background_model_data = {}
        bin_size = hic_ma.getBinSize()
        for referencePoint in referencePoints:
            region_start, region_end = viewpointObj.calculateViewpointRange(referencePoint, args.range)

            data_list = viewpointObj.computeViewpoint(referencePoint, referencePoint[0], region_start, region_end)

            if args.averageContactBin > 0:
                data_list = viewpointObj.smoothInteractionValues(data_list, args.averageContactBin)

            # set data in relation to viewpoint, upstream are negative values, downstream positive, zero is viewpoint
            view_point_start, _ = viewpointObj.getReferencePointAsMatrixIndices(referencePoint)
            view_point_range_start, view_point_range_end = \
                viewpointObj.getViewpointRangeAsMatrixIndices(referencePoint[0], region_start, region_end)
            for i, data in zip(range(view_point_range_start, view_point_range_end, 1), data_list):
                relative_position = i - view_point_start
                if relative_position in background_model_data:
                    background_model_data[relative_position] += data
                else:
                    background_model_data[relative_position] = data
                    relative_positions.add(relative_position)

        # compute the percentage of each position with respect to the total interaction count
        total_count = 0.0
        for key, value in background_model_data.items():
            total_count += value

        relative_counts = background_model_data

        for key in relative_counts:
            relative_counts[key] /= total_count
        relative_counts_conditions.append(relative_counts)

    # for models of all conditions:
    # - compute mean percentage for each bin
    # - compute SEM = Standard deviation / (square root of sample size)
    mean_percentage = {}
    sem = {}
    relative_positions = sorted(relative_positions)
    for relative_position in relative_positions:
        i = 0
        count = 0
        for condition in relative_counts_conditions:
            if relative_position in condition:
                i += 1
                count += condition[relative_position]
        mean_percentage[relative_position] = count / i
        sem[relative_position] = (count / i) / math.sqrt(i)

    with open('background_model.bed', 'w') as file:
        for relative_position in relative_positions:
            relative_position_in_genomic_scale = relative_position * bin_size
            file.write("{}\t{:.12f}\t{:.12f}\n".format(relative_position_in_genomic_scale, mean_percentage[relative_position],
                                                       sem[relative_position]))
