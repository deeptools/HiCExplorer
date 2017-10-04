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
from scipy.sparse import csr_matrix
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
                        'the output bam files of the processes into one output bam file. All other threads do the actual computation.'
                        ' Multithread computation is only used for larger matrices.',
                        required=False,
                        default=2,
                        type=int
                        )
    parser.add_argument('--chunkSize',
                        help='Chunk size defines how many elements per core should be computed. Decrease it to decrease memory usage.'
                        ' If the chunk size is greater as the shape of the matrix, only one compute thread is used (and one master thread). '
                        ' Only used with \'cool\' matrix format.',
                        required=False,
                        default=50000,
                        type=int
                        )
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def sum_cool_matrix(pMatrixList, pStartListElement, pQueue, pHic):

    matrix = pMatrixList[0] + pMatrixList[1]

    pQueue.put([matrix, pStartListElement])


def singleCore(args):
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
    return


def main(args=None):
    args = parse_arguments().parse_args(args)
    if args.matrices[0].endswith('.cool'):

        hic = hm.hiCMatrix(args.matrices[0], cooler_only_init=True)
        hic.load_cool_bins()
        dimension = hic.cooler_file.info['nbins']
        if dimension < 50000:
            print("Single core mode is used, matrix too small to split the computation.")
            hic = None
            singleCore(args)
        else:

            if args.threads < 2:
                exit("At least two threads are necessary. Given are: {}.".format(args.threads))
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

                chunk_size = args.chunkSize
                startX = 0
                startY = 0
                start_list = []
                chunk_size_x = chunk_size
                for i in range(0, dimension, chunk_size):
                    if i + chunk_size > dimension:
                        chunk_size_x = dimension % chunk_size
                    else:
                        chunk_size_x = chunk_size
                    for j in range(0, dimension, chunk_size):
                        if j + chunk_size < dimension:
                            start_list.append((i, j, chunk_size_x, chunk_size))
                        else:
                            start_list.append((i, j, chunk_size_x, dimension % chunk_size))

                col = []
                row = []
                data = []
                while chr_element < len(start_list) or not all_threads_done:
                    for i in range(args.threads):
                        if queue[i] is None and chr_element < len(start_list):
                            print (chr_element, " / ", len(start_list))
                            queue[i] = Queue()
                            process[i] = Process(target=sum_cool_matrix, kwargs=dict(
                                pMatrixList=[hic.load_cool_matrix_csr(start_list[chr_element]), hic_to_append.load_cool_matrix_csr(start_list[chr_element])],
                                pStartListElement=chr_element,
                                pQueue=queue[i],
                                pHic=hic
                            ))
                            process[i].start()
                            chr_element += 1
                            thread_done[i] = False

                        elif queue[i] is not None and not queue[i].empty():
                            matrix_, index = queue[i].get()
                            row.extend(matrix_.nonzero()[0] + start_list[index][0])
                            col.extend(matrix_.nonzero()[1] + start_list[index][1])
                            data.extend(matrix_.data)

                            matrix_ = None
                            queue[i] = None
                            process[i].join()
                            process[i].terminate()

                            process[i] = None
                            thread_done[i] = True
                        elif chr_element >= len(start_list) and queue[i] is None:
                            thread_done[i] = True
                        else:
                            time.sleep(0.01)
                    if chr_element >= len(start_list):
                        all_threads_done = True
                        for thread in thread_done:
                            if not thread:
                                all_threads_done = False

                hic.matrix = csr_matrix((data, (row, col)), shape=(len(row), len(col)))

                data = None
                row = None
                col = None
                hic.save(args.outFileName)
    else:
        singleCore(args)
