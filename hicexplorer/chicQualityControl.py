import argparse
import math
from multiprocessing import Process, Queue
import time
import logging
log = logging.getLogger(__name__)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from .lib import Viewpoint


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
Computes the sparsity of each viewpoint to determine the quality. A viewpoint is considered to be of bad quality if it is too sparse i.e. if there are too many locations with no interactions recorded.

This script creates three output files: a plot with the sparsity distribution per matrix, a plot with the sparsity distribution as histograms and a filtered reference points file.

An example usage is:

$ chicQualityControl -m matrix1.h5 matrix2.h5 -rp referencePointsFile.bed --range 20000 40000 --sparsity 0.01 -o referencePointFile_QC_passed.bed
"""
    )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='The input matrices to apply the QC on.',
                                nargs='+',
                                required=True)

    parserRequired.add_argument('--referencePoints', '-rp',
                                help='Bed file contains all reference points which are checked for a sufficient number of interactions.',
                                type=str,
                                required=True)
    parserRequired.add_argument('--sparsity', '-s',
                                help='Viewpoints with a sparsity less than the value given are considered of bad quality. If multiple matrices are given, '
                                'the viewpoint is removed as soon as it is of bad quality in at least one matrix.',
                                type=float,
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--outFileName', '-o',
                           help='The output file name of the passed reference points. Used as prefix for the plots as well.',
                           default='new_referencepoints.bed')
    parserOpt.add_argument('--outFileNameHistogram', '-oh',
                           help='The output file for the histogram plot.',
                           default='histogram.png')
    parserOpt.add_argument('--outFileNameSparsity', '-os',
                           help='The output file for the sparsity distribution plot.',
                           default='sparsity.png')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads.',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--fixateRange', '-fs',
                           help='Fixate score of background model starting at distance x. E.g. all values greater than 500kb are set to the value of the 500kb bin.',
                           required=False,
                           default=500000,
                           type=int
                           )
    parserOpt.add_argument('--dpi',
                           help='Optional parameter: Resolution for the image if the'
                           'output is a raster graphics image (e.g png, jpg)',
                           type=int,
                           default=300,
                           required=False)
    parserOpt.add_argument('--help', '-h', action='help',
                           help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def compute_sparsity(pReferencePoints, pViewpointObj, pArgs, pQueue):

    sparsity_list = []

    for i, referencePoint in enumerate(pReferencePoints):
        region_start, region_end, _ = pViewpointObj.calculateViewpointRange(
            referencePoint, (pArgs.fixateRange, pArgs.fixateRange))
        data_list = pViewpointObj.computeViewpoint(
            referencePoint, referencePoint[0], region_start, region_end)
        sparsity = (np.count_nonzero(data_list) / len(data_list))

        sparsity_list.append(sparsity)

    pQueue.put(sparsity_list)
    return


def main(args=None):
    args = parse_arguments().parse_args(args)

    viewpointObj = Viewpoint()
    referencePoints, _ = viewpointObj.readReferencePointFile(
        args.referencePoints)

    # compute for each viewpoint the sparsity and consider these as bad with a sparsity less than given.

    referencePointsPerThread = len(referencePoints) // args.threads
    queue = [None] * args.threads
    process = [None] * args.threads
    sparsity = []

    for j, matrix in enumerate(args.matrices):
        sparsity_local = [None] * args.threads
        hic_ma = hm.hiCMatrix(matrix)
        viewpointObj.hicMatrix = hic_ma

        all_data_collected = False
        thread_done = [False] * args.threads
        for i in range(args.threads):

            if i < args.threads - 1:
                referencePointsThread = referencePoints[i * referencePointsPerThread:(i + 1) * referencePointsPerThread]
            else:
                referencePointsThread = referencePoints[i * referencePointsPerThread:]
            if len(referencePointsThread) == 0:
                process[i] = None
                queue[i] = None
                sparsity_local[i] = []
                continue
            else:
                queue[i] = Queue()
                process[i] = Process(target=compute_sparsity, kwargs=dict(
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
                    sparsity_ = queue[i].get()
                    sparsity_local[i] = sparsity_
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

        # merge sparsity data per matrix from each thread to one list

        sparsity_local = [
            item for sublist in sparsity_local for item in sublist]
        sparsity.append(sparsity_local)

    # change sparsity to sparsity values per viewpoint per matrix: viewpoint = [matrix1, ..., matrix_n]
    sparsity = np.array(sparsity).T

    with open(args.referencePoints, 'r') as reference_file_input:
        with open(args.outFileName + '_raw_filter', 'w') as output_file_raw:
            output_file_raw.write('# Created with chicQualityControl version {}\n'.format(__version__))
            output_file_raw.write('# Chromosome\tStart\tEnd\t')
            for matrix in args.matrices:
                output_file_raw.write('Sparsity {}\t'.format(matrix))
            output_file_raw.write('\n')

            with open(args.outFileName + '_rejected_filter', 'w') as output_file_rejected:
                with open(args.outFileName, 'w') as output_file:
                    for i, line in enumerate(reference_file_input.readlines()):
                        sparsity_str = '\t'.join(str(x) for x in sparsity[i])
                        output_file_raw.write(
                            line.strip() + '\t' + sparsity_str + '\n')
                        count = 0
                        for j in range(len(sparsity[i])):
                            if sparsity[i][j] > args.sparsity:
                                count += 1
                        if count:
                            output_file.write(line)
                        else:
                            output_file_rejected.write(line)
    # output plot of sparsity distribution per sample

    # re-arange values again

    x = [[]] * len(args.matrices)
    y = [[]] * len(args.matrices)

    for i in range(len(args.matrices)):
        y[i] = [i] * len(sparsity)
    sparsity = sparsity.T

    for i in range(len(args.matrices)):
        x[i] = sparsity[i].flatten()

    for i in range(len(args.matrices)):
        plt.plot(x[i], y[i], 'o', mfc='none', markersize=0.3,
                 label=args.matrices[i].split('/')[-1])
    plt.yticks([])
    plt.xlabel("Sparsity level")

    plt.axvline(x=args.sparsity, c='r', label='sparsity threshold', linewidth=0.3)
    plt.xscale('log')
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
    plt.legend(loc='center', bbox_to_anchor=(1.4, 0.5))
    plt.savefig(args.outFileNameSparsity, dpi=args.dpi)

    # plt.xlabel("Length of list (number)")
    # plt.ylabel("Time taken (seconds)")
    plt.close()
    for i in range(len(args.matrices)):
        plt.hist(x[i], bins=100, alpha=0.5, label=args.matrices[i].split('/')[-1])
    plt.xlabel("Sparsity level")
    plt.ylabel("Number of counts")
    # plt.legend(loc='upper right')

    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
    plt.legend(loc='center', bbox_to_anchor=(1.4, 0.5))
    plt.savefig(args.outFileNameHistogram, dpi=args.dpi)
