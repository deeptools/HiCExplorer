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

from  scipy.sparse import csr_matrix, vstack

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


def sum_cool_matrix(pMatrix, pMatrixToAppend, pNanBins, pNanBinsAppend, pQueue):

    matrix = pMatrix + pMatrixToAppend
    # pNanBins = set(pNanBins)                        
    nan_bins = np.append(pNanBins, pNanBinsAppend)
    # if 
    # print("pMatrixList[0]", pMatrixList[0])
    # print("pMatrixList[0]", pMatrixList[0])
    # print pMatrixList[0].values
    # print csr_matrix(pMatrixList[0], shape=(pHic.cooler_matrix_shape, pHic.cooler_matrix_shape,))

    # cool_pandas_bins = pHic.compute_dataframe_bins(pBinsList, "+")    
    # cool_matrix_pixel = pHic.compute_dataframe_matrix(pMatrixList, "+")
    pQueue.put([matrix, nan_bins])


def main(args=None):
    args = parse_arguments().parse_args(args)
    if args.matrices[0].endswith('.cool'):
        if args.threads < 2:
            exit("At least two threads are necessary. Given are: {}.".format(args.threads))
        hic = hm.hiCMatrix(args.matrices[0], cooler_only_init=True)
        hic.load_cool_cut_intervals()
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
            # dataFrameBins = None
            # dataFrameMatrix = None
            matrix = None
            nan_bins = None
            while chr_element < len(chromosome_list) or not all_threads_done:
                for i in range(args.threads):
                    if queue[i] is None and chr_element < len(chromosome_list):
                        print("Computing: ", chromosome_list[chr_element])
                        
                        queue[i] = Queue()
                        chromosome = chromosome_list[chr_element]
                        matrix_org, nanbins_org = hic.load_cool_per_chr_csr(chromosome)
                        matrix_append, nanbins_append = hic_to_append.load_cool_per_chr_csr(chromosome)
                        process[i] = Process(target=sum_cool_matrix, kwargs=dict(
                            pMatrix=matrix_org,
                            pMatrixToAppend=matrix_append,
                            pNanBins=nanbins_org,
                            pNanBinsAppend=nanbins_append,
                            pQueue=queue[i]
                        ))
                        process[i].start()
                        chr_element += 1
                        thread_done[i] = False
                    elif queue[i] is not None and not queue[i].empty():
                        matrix_, nan_bins_ = queue[i].get()
                        # print("matrix_", matrix_)
                        if matrix is None:
                            matrix = matrix_
                        else:
                            matrix += matrix_
                        if nan_bins is None:
                            nan_bins = nan_bins_
                        else:
                            nan_bins = np.append(nan_bins, nan_bins_)
                      
                        matrix_ = None
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

            # print("matrixFOOF", matrix)
            
            hic.setMatrixValues(matrix)
            nan_bins = set(nan_bins)
            hic.setNanBins(nan_bins)
            hic.save(args.outFileName)
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
