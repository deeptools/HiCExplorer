from __future__ import division
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
import math
from multiprocessing import Process, Queue
import time
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

    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--averageContactBin',
                           help='Average the contacts of n bins, written to last column.',
                           type=int,
                           default=0)
    parserOpt.add_argument('--outFileName', '-o',
                           help='The name of the background model file',
                           default='background_model.bed')
    parserOpt.add_argument('--threads',
                           help='Number of threads. Using the python multiprocessing module. ',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--fixateRange', '-fs',
                           help='Fixate score of backgroundmodel starting at distance x. E.g. all values greater 500kb are set to the value of the 500kb bin.',
                           required=False,
                           default=500000,
                           type=int
                           )
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def compute_background(pReferencePoints, pViewpointObj, pArgs, pQueue):

    background_model_data = {}
    relative_positions = set()
    for i, referencePoint in enumerate(pReferencePoints):

        region_start, region_end, _ = pViewpointObj.calculateViewpointRange(referencePoint, (pArgs.fixateRange, pArgs.fixateRange))

        data_list = pViewpointObj.computeViewpoint(referencePoint, referencePoint[0], region_start, region_end)

        if pArgs.averageContactBin > 0:
            data_list = pViewpointObj.smoothInteractionValues(data_list, pArgs.averageContactBin)

        # set data in relation to viewpoint, upstream are negative values, downstream positive, zero is viewpoint
        view_point_start, _ = pViewpointObj.getReferencePointAsMatrixIndices(referencePoint)
        view_point_range_start, view_point_range_end = \
            pViewpointObj.getViewpointRangeAsMatrixIndices(referencePoint[0], region_start, region_end)

        for i, data in zip(range(view_point_range_start, view_point_range_end, 1), data_list):
            relative_position = i - view_point_start
            if relative_position in background_model_data:
                background_model_data[relative_position] += data
            else:
                background_model_data[relative_position] = data
                relative_positions.add(relative_position)
    pQueue.put([background_model_data, relative_positions])
    return


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

    referencePointsPerThread = len(referencePoints) // args.threads
    queue = [None] * args.threads
    process = [None] * args.threads
    for matrix in args.matrices:
        hic_ma = hm.hiCMatrix(matrix)
        # hic_ma.
        viewpointObj.hicMatrix = hic_ma
        background_model_data = None
        bin_size = hic_ma.getBinSize()
        all_data_collected = False
        thread_done = [False] * args.threads
        log.debug('matrix read, starting processing')
        for i in range(args.threads):

            if i < args.threads - 1:
                referencePointsThread = referencePoints[i * referencePointsPerThread:(i + 1) * referencePointsPerThread]
            else:
                referencePointsThread = referencePoints[i * referencePointsPerThread:]

            queue[i] = Queue()
            process[i] = Process(target=compute_background, kwargs=dict(
                pReferencePoints=referencePointsThread,
                pViewpointObj=viewpointObj,
                pArgs=args,
                pQueue=queue[i]
            )
            )

            process[i].start()

        while not all_data_collected:
            for i in range(args.threads):
                if queue[i] is not None and not queue[i].empty():
                    foo = queue[i].get()
                    # log.debug('len(queue[i].get() {}'.format(len(foo)))
                    background_model_data_thread, relative_positions_thread = foo
                    if background_model_data is None:
                        background_model_data = background_model_data_thread
                    else:
                        for relativePosition in background_model_data_thread:
                            if relativePosition in background_model_data:
                                background_model_data[relativePosition] += background_model_data_thread[relativePosition]
                            else:
                                background_model_data[relativePosition] = background_model_data_thread[relativePosition]
                        # log.debug('relative_positions_thread {}'.format(relative_positions_thread))

                    relative_positions = relative_positions.union(relative_positions_thread)
                    queue[i] = None
                    process[i].join()
                    process[i].terminate()
                    process[i] = None
                    thread_done[i] = True
            all_data_collected = True
            for thread in thread_done:
                if not thread:
                    all_data_collected = False
            time.sleep(1)

        del hic_ma
        del viewpointObj.hicMatrix

        total_count = sum(background_model_data.values())
        relative_counts = background_model_data
        # for the current condition compute the relative interactions per relative distance
        for key in relative_counts:
            relative_counts[key] /= total_count

        relative_counts_conditions.append(relative_counts)

    # log.debug('background_model_data {}'.format(background_model_data))
    # log.debug('relative_positions {}'.format(relative_positions))

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
            # count the number of relative interactions at 'relative_position' if this position exists in the condition
            if relative_position in condition:
                i += 1
                count += condition[relative_position]
        # mean_percentage is given by number of relative interactions at a relative position divided by the number of conditions
        # log.debug('relative pos: {} count: {}'.format(relative_position, count))
        mean_percentage[relative_position] = count / i
        sem[relative_position] = (count / i) / math.sqrt(i)
    # lower_limit_range = args.range[0] * (-1) / bin_size
    # upper_limit_range = args.range[1]  / bin_size

    # write result to file
    with open(args.outFileName, 'w') as file:
        for relative_position in relative_positions:
            # if lower_limit_range <= relative_position and relative_position <= upper_limit_range:
            relative_position_in_genomic_scale = relative_position * bin_size
            file.write("{}\t{:.12f}\t{:.12f}\n".format(relative_position_in_genomic_scale, mean_percentage[relative_position],
                                                       sem[relative_position]))
