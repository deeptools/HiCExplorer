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

import operator



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


def sum_cool_matrix(pBinsList, pMatrixList, pQueue, pHic):

    cool_pandas_bins = pHic.compute_dataframe_bins(pBinsList, "+")    
    cool_matrix_pixel = pHic.compute_dataframe_matrix(pMatrixList, "+")
    pQueue.put([cool_pandas_bins, cool_matrix_pixel])


def main(args=None):
    args = parse_arguments().parse_args(args)
    if args.matrices[0].endswith('.cool'):
        if args.threads < 2:
            exit("At least two threads are necessary. Given are: {}.".format(args.threads))
        hic = hm.hiCMatrix(args.matrices[0], cooler_only_init=True)
        chromosome_list = hic.cooler_file.chromnames
        args.threads = args.threads - 1
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
                            pQueue=queue[i],
                            pHic=hic
                        ))
                        process[i].start()
                        chr_element += 1
                        thread_done[i] = False
                    elif queue[i] is not None and not queue[i].empty():
                        dataFrameBins_, dataFrameMatrix_ = queue[i].get()
                        if dataFrameBins_ is not None:
                            if dataFrameBins is None:
                                dataFrameBins = dataFrameBins_
                            else:
                                dataFrameBins = pd.concat([dataFrameBins, dataFrameBins_], ignore_index=True)
                        if dataFrameMatrix_ is not None:
                            if dataFrameMatrix is None:
                                dataFrameMatrix = dataFrameMatrix_
                            else:
                                dataFrameMatrix = pd.concat([dataFrameMatrix, dataFrameMatrix_], ignore_index=True)
                        dataFrameBins_ = None
                        dataFrameMatrix_ = None
                        queue[i] = None
                        process[i].join()
                        process[i].terminate()
                        
                        process[i] = None
                        thread_done[i] = True
                    elif chr_element >= len(chromosome_list) and queue[i] is None:
                        thread_done[i] = True
                    else:
                        time.sleep(0.1)
                if chr_element >= len(chromosome_list):
                    all_threads_done = True
                    for thread in thread_done:
                        if not thread:
                            all_threads_done = False
            hic.save_cool_pandas(args.outFileName, dataFrameBins, dataFrameMatrix)
    else:
        if args.threads:
            print("Multiple threads are only used if the matrices are in 'cool'-format.")
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
