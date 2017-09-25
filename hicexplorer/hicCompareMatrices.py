from __future__ import division

import argparse
import numpy as np
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
from multiprocessing import Process
from multiprocessing import Process, Lock, Queue
import pandas as pd
import operator
import time



def parse_arguments(args=None):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=('Takes two matrices, normalizes them and applies'
                                                  'the given operation. To normalize the matrices '
                                                  'each element is divided by sum of the matrix.'))

    parser.add_argument('--matrices', '-m',
                        help='matrices to use.',
                        metavar='.h5 file format',
                        nargs=2,
                        required=True)

    parser.add_argument('--outFileName', '-o',
                        help='File name to save the resulting matrix. The output is '
                             'also a .h5 file.',
                        required=True)

    parser.add_argument('--operation',
                        help='Operation to apply for the matrices. Options are: diff, ratio, log2ratio',
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

def compute_cool_matrix(pBinsList, pMatrixList, pQueue, pHic, pOperation):
   
    # normalize by total matrix sum
    sum_matrix_0 = pMatrixList[0]['count'].sum()
    sum_matrix_1 = pMatrixList[1]['count'].sum()

    pMatrixList[0]['count'] = pMatrixList[0]['count'].apply(lambda x: x / sum_matrix_0)
    pMatrixList[1]['count'] = pMatrixList[1]['count'].apply(lambda x: x / sum_matrix_0)

    cool_pandas_bins = pHic.compute_dataframe_bins(pBinsList, "+")    

    if pMatrixList[0] is None:
        print("First is none")
    if pMatrixList[1] is None:
        print("Second is none")
    new_matrix = None
    if pOperation == 'diff':
        new_matrix = pHic.compute_dataframe_matrix(pMatrixList, "-")
    elif pOperation == 'ratio' or pOperation == 'log2ratio':
        pMatrixList[1]['count'] = pMatrixList[1]['count'].apply(lambda x: 1 / x)
        new_matrix = pHic.compute_dataframe_matrix(pMatrixList, "*")
        if pOperation == 'log2ratio':
            if new_matrix is not None:
                new_matrix['count'] = new_matrix['count'].apply(np.log2)       
    pQueue.put([cool_pandas_bins, new_matrix])
    return



def main(args=None):

    args = parse_arguments().parse_args(args)
    if args.operation not in ['diff', 'ratio', 'log2ratio']:
        exit("Operation not found. Please use 'diff', 'ratio' or 'log2ratio'.")

    if args.matrices[0].endswith('.cool'):
        if args.threads < 2:
            exit("At least two threads are necessary. Given are: {}.".format(args.threads))
        
        args.threads = args.threads - 1
        process = [None] * args.threads
        lock = Lock()
        queue = [None] * args.threads

        hic = hm.hiCMatrix(args.matrices[0], cooler_only_init=True)
        hic_second_matrix = hm.hiCMatrix(args.matrices[1], cooler_only_init=True)        
        chromosome_list = hic.cooler_file.chromnames
        chromosome_list_second_matrix = hic_second_matrix.cooler_file.chromnames
        thread_done = [False] * args.threads
        all_threads_done = False
        
        dataFrameBins = None
        dataFrameMatrix = None
        chr_element = 0
        # TODO: shape check is missing
        if chromosome_list != chromosome_list_second_matrix:
            exit("The two matrices have different chromosome order. Use the tool `hicExport` to change the order.\n"
                    "{}: {}\n"
                    "{}: {}".format(args.matrices[0], chromosome_list,
                                    args.matrices[1], chromosome_list_second_matrix))
        while chr_element < len(chromosome_list) or not all_threads_done:
            for i in range(args.threads):
                if queue[i] is None and chr_element < len(chromosome_list):
                    queue[i] = Queue()
                    chromosome = chromosome_list[chr_element]
                    process[i] = Process(target=compute_cool_matrix, kwargs=dict(
                        pBinsList=[hic.load_cool_bins(chromosome), hic_second_matrix.load_cool_bins(chromosome)],
                        pMatrixList=[hic.load_cool_matrix(chromosome), hic_second_matrix.load_cool_matrix(chromosome)],
                        pQueue=queue[i],
                        pHic=hic,
                        pOperation=args.operation
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
                            dataFrameBins = pd.concat([dataFrameBins, dataFrameBins_], ignore_index=True)  # .append(dataFrameBins_, ignore_index=True)
                    if dataFrameMatrix_ is not None:
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
                else:
                    time.sleep(0.1)
                if chr_element >= len(chromosome_list):
                    all_threads_done = True
                    for thread in thread_done:
                        if not thread:
                            all_threads_done = False
        
        # TODO
        hic.save_cool_pandas(args.outFileName, dataFrameBins, dataFrameMatrix)

                
    else:
        hic1 = hm.hiCMatrix(args.matrices[0])
        hic2 = hm.hiCMatrix(args.matrices[1])

        if hic1.matrix.shape != hic2.matrix.shape:
            exit("The two matrices have different size. Use matrices having the same resolution and created using"
                "the same parameters. Check the matrix values using the tool `hicInfo`.")

        if hic1.chrBinBoundaries != hic2.chrBinBoundaries:
            exit("The two matrices have different chromosome order. Use the tool `hicExport` to change the order.\n"
                "{}: {}\n"
                "{}: {}".format(args.matrices[0], hic1.chrBinBoundaries.keys(),
                                args.matrices[1], hic2.chrBinBoundaries.keys()))

        # normalize by total matrix sum
        hic1.matrix.data = hic1.matrix.data.astype(float) / hic1.matrix.data.sum()
        hic2.matrix.data = hic2.matrix.data.astype(float) / hic2.matrix.data.sum()

        nan_bins = set(hic1.nan_bins)
        nan_bins = nan_bins.union(hic2.nan_bins)

        if args.operation == 'diff':
            new_matrix = hic1.matrix - hic2.matrix
        elif args.operation == 'ratio' or args.operation == 'log2ratio':
            hic2.matrix.data = float(1) / hic2.matrix.data
            new_matrix = hic1.matrix.multiply(hic2.matrix)
            # just in case
            new_matrix.eliminate_zeros()
            if args.operation == 'log2ratio':
                new_matrix.data = np.log2(new_matrix.data)
                new_matrix.eliminate_zeros()

        hic1.setMatrixValues(new_matrix)
        hic1.maskBins(sorted(nan_bins))
        hic1.save(args.outFileName)
