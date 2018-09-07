from __future__ import division
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
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
                    Build a background model for viewpoints, the viewpoints over all samples (matrices) are used to build the background model.
                    An example usage is:
                    $ chicViewpointBackgroundModel -m matrix1.h5 matrix2.h5 matrix3.h5 -rp referencePointsFile.bed --range 20000 40000 
                    """
    )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='The input matrices to build the background model on.',
                                nargs='+',
                                required=True)

    parserRequired.add_argument('--referencePoints', '-rp',
                                help='Bed file contains all reference points which should be used to build the background model.',
                                type=str,
                                required=True)
    parserRequired.add_argument('--range', '-r',
                                help='Region around the viewpoint to consider in basepairs, first is upstream, second downstream.'
                                ' E.g. \'--range 20000 40000\' for 20kb upstream to 40kb downstream of a viewpoint. Do NOT use a \'-\' prefix for upstream range.',
                                type=int,
                                required=True,
                                nargs=2)
    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--averageContactBin',
                           help='Average the contacts of n bins, written to last column.',
                           type=int,
                           default=0)
    parserOpt.add_argument('--outFileName', '-o',
                           help='The name of the background model file',
                           default='background_model.bed')
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def main():
    args = parse_arguments().parse_args()

    viewpointObj = Viewpoint()
    referencePoints, _ = viewpointObj.readReferencePointFile(args.referencePoints)

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
        viewpointObj.hicMatrix = hic_ma
        background_model_data = {}
        bin_size = hic_ma.getBinSize()
        for referencePoint in referencePoints:

            region_start, region_end, _ = viewpointObj.calculateViewpointRange(referencePoint, args.range)

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

        total_count = sum(background_model_data.values())

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

    # write result to file
    with open(args.outFileName, 'w') as file:
        for relative_position in relative_positions:
            relative_position_in_genomic_scale = relative_position * bin_size
            file.write("{}\t{:.12f}\t{:.12f}\n".format(relative_position_in_genomic_scale, mean_percentage[relative_position],
                                                       sem[relative_position]))
