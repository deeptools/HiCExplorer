from __future__ import division
import argparse
from hicexplorer import HiCMatrix as hm
# from hicexplorer
from hicexplorer._version import __version__
from multiprocessing import Process, Lock, Queue
import time
import os
import pandas as pd
import numpy as np


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=('Adds Hi-C matrices of the same size. Format '
                                                  'has to be hdf5 or npz'))

    parser.add_argument('--matrices', '-m',
                        help='matrices to add. Must have the same shape.',
                        metavar='.h5 file format',
                        nargs='+',
                        required=True)

    parser.add_argument('--outFileName', '-o',
                        help='File name to save the resulting matrix. The output is '
                             'also a .h5 file. But don\'t add the suffix',
                        required=True)
    parser.add_argument('--threads',
                        help='Number of threads. Using the python multiprocessing module. Is only used with \'cool\' matrix format.'
                        ' One master process which is used to read the input file into the buffer and one process which is merging '
                        'the output bam files of the processes into one output bam file. All other threads do the actual computation.',
                        required=False,
                        default=4,
                        type=int
                        )
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def sum_cool_matrix(pBinsList, pMatrixList, pQueue):

    cut_intervals_0 = [tuple(x) for x in pBinsList[0].values]
    cut_intervals_1 = [tuple(x) for x in pBinsList[1].values]
    cut_intervals_tuple = []
    pBinsList[0] = None
    pBinsList[1] = None

    i = 0
    j = 0
    while i < len(cut_intervals_0) and j < len(cut_intervals_1):

        if cut_intervals_0[i][1] == cut_intervals_1[j][1] and \
                cut_intervals_0[i][2] == cut_intervals_1[j][2]:
            # add only if both are not nan.
            # if one is nan, use the other. this is either a number or nan too.
            if not np.isnan(cut_intervals_0[i][3]) and not np.isnan(cut_intervals_1[j][3]):
                cut_intervals_tuple.append((cut_intervals_0[i][0], cut_intervals_0[i][1], cut_intervals_0[i][2],
                                            cut_intervals_0[i][3] + cut_intervals_1[j][3]))
            elif np.isnan(cut_intervals_0[i][3]):
                cut_intervals_tuple.append((cut_intervals_0[i][0], cut_intervals_0[i][1], cut_intervals_0[i][2],
                                            cut_intervals_1[j][3]))
            else:
                cut_intervals_tuple.append((cut_intervals_0[i][0], cut_intervals_0[i][1], cut_intervals_0[i][2],
                                            cut_intervals_0[i][3]))
            cut_intervals_0[i] = None
            cut_intervals_1[j] = None
            i += 1
            j += 1

        elif cut_intervals_0[i][0] == cut_intervals_1[j][0] and \
                cut_intervals_0[i][1] < cut_intervals_1[j][1]:
            cut_intervals_tuple.append(cut_intervals_0[i])
            cut_intervals_0[i] = None
            i += 1
        elif cut_intervals_0[i][0] == cut_intervals_1[j][0] and \
                cut_intervals_0[i][1] > cut_intervals_1[j][1]:
            cut_intervals_tuple.append(cut_intervals_1[j])
            cut_intervals_1[j] = None
            j += 1
        elif cut_intervals_0[i][0] < cut_intervals_1[j][0]:
            cut_intervals_tuple.append(cut_intervals_0[i])
            cut_intervals_0[i] = None
            i += 1
        else:
            cut_intervals_tuple.append(cut_intervals_1[j])
            cut_intervals_1[j] = None
            j += 1

    while i < len(cut_intervals_0):
        cut_intervals_tuple.append(cut_intervals_0[i])
        cut_intervals_0[i] = None
        i += 1
    while j < len(cut_intervals_1):
        cut_intervals_tuple.append(cut_intervals_1[j])
        cut_intervals_1[j] = None
        j += 1

    cool_pandas_bins = pd.DataFrame(cut_intervals_tuple, columns=['chrom', 'start', 'end', 'weight'])
    cut_intervals_tuple = None

    matrix_0 = [tuple(x) for x in pMatrixList[0].values]
    matrix_1 = [tuple(x) for x in pMatrixList[1].values]
    matrix_tuple = []
    pMatrixList[0] = None
    pMatrixList[1] = None
    i = 0
    j = 0

    while i < len(matrix_0) and j < len(matrix_1):
        if matrix_0[i][1] == matrix_1[j][1] and \
                matrix_0[i][2] == matrix_1[j][2]:
            # if value is nan, do not add and take the other one.
            if not np.isnan(matrix_0[i][2]) and not np.isnan(matrix_1[j][2]):
                matrix_tuple.append((matrix_0[i][0], matrix_0[i][1],
                                     matrix_0[i][2] + matrix_1[j][2]))
            elif np.isnan(matrix_0[i][2]):
                matrix_tuple.append((matrix_0[i][0], matrix_0[i][1],
                                     matrix_1[j][2]))
            else:
                matrix_tuple.append((matrix_0[i][0], matrix_0[i][1],
                                     matrix_0[i][2]))
            matrix_0[i] = None
            matrix_1[j] = None
            i += 1
            j += 1
        elif matrix_0[i][0] == matrix_1[j][0] and matrix_0[i][1] < matrix_1[j][1]:
            matrix_tuple.append(matrix_0[i])
            matrix_0[i] = None
            i += 1
        elif matrix_0[i][0] == matrix_1[j][0] and matrix_0[i][1] > matrix_1[j][1]:
            matrix_tuple.append(matrix_1[j])
            cut_intervals_1[j] = None
            j += 1
        elif matrix_0[i][0] < matrix_1[j][0]:
            matrix_tuple.append(matrix_0[i])
            cut_intervals_0[i] = None
            i += 1
        else:
            matrix_tuple.append(matrix_1[j])
            cut_intervals_1[j] = None
            j += 1
    while i < len(matrix_0):
        matrix_tuple.append(matrix_0[i])
        matrix_0[i] = None
        i += 1
    while j < len(matrix_1):
        matrix_tuple.append(matrix_1[j])
        matrix_1[j] = None
        j += 1

    cool_matrix_pixel = pd.DataFrame(matrix_tuple, columns=['bin1_id', 'bin2_id', 'count'])
    matrix_tuple = None
    pQueue.put([cool_pandas_bins, cool_matrix_pixel])


def main(args=None):

    args = parse_arguments().parse_args(args)
    if args.matrices[0].endswith('.cool'):
        hic = hm.hiCMatrix(args.matrices[0], cooler_only_init=True)
        chromosome_list = hic.cooler_file.chromnames
        process = [None] * args.threads
        lock = Lock()
        queue = [None] * args.threads

        for matrix in args.matrices[1:]:
            hic_to_append = hm.hiCMatrix(matrix, cooler_only_init=True)

            chromosome_list_to_append = hic_to_append.cooler_file.chromnames
            if chromosome_list != chromosome_list_to_append:
                exit("The two matrices have different chromosome order. Use the tool `hicExport` to change the order.\n"
                     "{}: {}\n"
                     "{}: {}".format(args.matrices[0], chromosome_list,
                                     matrix, chromosome_list_to_append))

            chr_element = 0
            thread_done = [False] * args.threads
            all_threads_done = False
            first_to_save = True
            dataFrameBins = None
            dataFrameMatrix = None

            while chr_element < len(chromosome_list) or not all_threads_done:
                for i in range(args.threads):
                    if queue[i] is None and chr_element < len(chromosome_list):

                        queue[i] = Queue()
                        chromosome = chromosome_list[chr_element]
                        process[i] = Process(target=sum_cool_matrix, kwargs=dict(
                            pBinsList=[hic.load_cool_bins(chromosome), hic_to_append.load_cool_bins(chromosome)],
                            pMatrixList=[hic.load_cool_matrix(chromosome), hic_to_append.load_cool_matrix(chromosome)],
                            pQueue=queue[i]
                        ))
                        process[i].start()
                        chr_element += 1
                    elif queue[i] is not None and not queue[i].empty():
                        dataFrameBins_, dataFrameMatrix_ = queue[i].get()
                        if dataFrameBins is None:
                            dataFrameBins = dataFrameBins_
                        else:
                            dataFrameBins = pd.concat([dataFrameBins, dataFrameBins_], ignore_index=True)  # .append(dataFrameBins_, ignore_index=True)
                        if dataFrameMatrix is None:
                            dataFrameMatrix = dataFrameMatrix_
                        else:
                            dataFrameMatrix = pd.concat([dataFrameMatrix, dataFrameMatrix_], ignore_index=True)  # .append(dataFrameMatrix_, ignore_index=True)
                        dataFrameBins_ = None
                        dataFrameMatrix_ = None
                        queue[i] = None
                        process[i].join()
                        process[i].terminate()
                        process[i] = None
                        thread_done[i] = True
                    elif chr_element >= len(chromosome_list) and queue[i] is None:
                        thread_done[i] = True
                if chr_element >= len(chromosome_list):
                    all_threads_done = True
                    for thread in thread_done:
                        if not thread:
                            all_threads_done = False

            hic.save_cool_pandas(args.outFileName, dataFrameBins, dataFrameMatrix)
    else:
        hic = hm.hiCMatrix(args.matrices[0])
        summed_matrix = hic.matrix
        nan_bins = set(hic.nan_bins)
        for matrix in args.matrices[1:]:
            hic_to_append = hm.hiCMatrix(matrix)
            if hic.chrBinBoundaries != hic_to_append.chrBinBoundaries:
                exit("The two matrices have different chromosome order. Use the tool `hicExport` to change the order.\n"
                     "{}: {}\n"
                     "{}: {}".format(args.matrices[0], list(hic.chrBinBoundaries),
                                     matrix, list(hic_to_append.chrBinBoundaries)))

            try:
                summed_matrix = summed_matrix + hic_to_append.matrix
                if len(hic_to_append.nan_bins):
                    nan_bins = nan_bins.union(hic_to_append.nan_bins)
            except:
                print("\nMatrix {} seems to be corrupted or of different shape".format(matrix))
                exit(1)

        # save only the upper triangle of the
        # symmetric matrix
        hic.setMatrixValues(summed_matrix)
        hic.maskBins(sorted(nan_bins))
        hic.save(args.outFileName)
