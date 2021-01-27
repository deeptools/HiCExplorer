import argparse
import sys
import os
import errno
import math
from multiprocessing import Process, Queue
import time
import traceback

import logging
log = logging.getLogger(__name__)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

import hicmatrix.HiCMatrix as hm
from hicexplorer import utilities
from hicexplorer._version import __version__
from .lib import Viewpoint

import h5py
import io
import tarfile
from contextlib import closing


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
chicExportData exports the data stored in the intermediate hdf5 files to text files per reference point.

"""
)
    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--file', '-f',
                                help='path to the file which should be used for data export',
                                required=True)


    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--outFileName', '-o',
                           help='Output tar.gz of the files'
                           ' (Default: %(default)s).',
                           required=False,
                           default='data.tar.gz')
  
    parserOpt.add_argument('--fileType',
                        '-ft',
                        help=''
                        ' (Default: %(default)s).',
                        default='interaction',
                        choices=['interaction', 'significant', 'target', 'aggregated']
                        )

    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads (uses the python multiprocessing module)'
                           ' (Default: %(default)s).',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser

def exportData(pFileList, pArgs, pViewpointObj, pQueue):

    file_list = []
    file_content_list = []
    try:
        if pArgs.fileType == 'interaction' or pArgs.fileType == 'significant':
            header_information = '# Chromosome\tStart\tEnd\tGene\tSum of interactions\tRelative position\tRelative Interactions\tp-value\tx-fold\tRaw\n'

            for file in pFileList:
                for sample in file:
                    data = pViewpointObj.readInteractionFile(pArgs.file, sample)
                    
                    file_content_string = header_information
                    key_list = sorted(list(data[1].keys()))
                    for key in key_list:
                        file_content_string += '\t'.join(str(x) for x in data[1][key]) + '\n'
                file_content_list.append(file_content_string)
                file_name = '_'.join(sample) + '.txt'
                file_list.append(file_name)
        elif args.fileType == 'target':
            pass
        elif args.fileType == 'aggregated':
            pass
        
        elif args.fileType == 'differential':
            pass
    
    except Exception as exp:
        log.debug("FAIL: {}".format(str(exp) + traceback.format_exc()))
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
        return
 
    pQueue.put([file_list, file_content_list])
    log.debug('RETRUN')
    return

def main(args=None):
    args = parse_arguments().parse_args(args)
    viewpointObj = Viewpoint()
    
    fileList = []


    ### read hdf file
    fileHDF5Object = h5py.File(args.file, 'r')
    keys_file = list(fileHDF5Object.keys())       

    if args.fileType == 'interaction' or args.fileType == 'significant':
        
        # resolution = interactionFileHDF5Object.attrs['resolution'][()]

        if len(keys_file) > 1:
            for i, sample in enumerate(keys_file):
                    
                    matrix_obj1 = fileHDF5Object[sample]
                    chromosomeList1 = sorted(list(matrix_obj1.keys()))
                    chromosomeList1.remove('genes')
                    for chromosome1 in chromosomeList1:
                        geneList1 = sorted(list(matrix_obj1[chromosome1].keys()))
                        for gene1 in geneList1:
                            fileList.append([[sample,chromosome1, gene1]])

    
     # log.debug(interactionFileList[:10])
    elif args.fileType == 'target':
        
        for outer_matrix in keys_file:
            inner_matrix_object = fileHDF5Object[outer_matrix]
            keys_inner_matrices = list(inner_matrix_object.keys())
            for inner_matrix in keys_inner_matrices:
                inner_object = inner_matrix_object[inner_matrix]
                gene_object = inner_object['genes']
                keys_genes = list(gene_object.keys())
                for gen in keys_genes:
                    fileList.append([outer_matrix, inner_matrix, 'genes', gen])
                    
        
    elif args.fileType == 'aggregated':

        for i, combinationOfMatrix in enumerate(keys_file):
            # log.debug('list(aggregatedFileHDF5Object[combinationOfMatrix].keys()) {}'.format(list(aggregatedFileHDF5Object[combinationOfMatrix].keys())))
            keys_matrix_intern = list(fileHDF5Object[combinationOfMatrix].keys())
            if len(keys_matrix_intern) == 0:
                continue
        # if len(keys_aggregatedFile) > 1:

            log.debug('combinationOfMatrix {} keys_matrix_intern {}'.format(combinationOfMatrix, keys_matrix_intern))
            matrix1 = keys_matrix_intern[0]
            matrix2 = keys_matrix_intern[1]

            matrix_obj1 = aggregatedFileHDF5Object[combinationOfMatrix + '/' + matrix1]
            matrix_obj2 = aggregatedFileHDF5Object[combinationOfMatrix + '/' + matrix2]
            # for 
            # for sample2 in keys_aggregatedFile[i + 1:]:
                
            #     matrix_obj1 = aggregatedFileHDF5Object[sample]
            #     matrix_obj2 = aggregatedFileHDF5Object[sample]

            chromosomeList1 = sorted(list(matrix_obj1.keys()))
            chromosomeList2 = sorted(list(matrix_obj2.keys()))
            chromosomeList1.remove('genes')
            chromosomeList2.remove('genes')
            for chromosome1, chromosome2 in zip(chromosomeList1, chromosomeList2):
                geneList1 = sorted(list(matrix_obj1[chromosome1].keys()))
                geneList2 = sorted(list(matrix_obj2[chromosome2].keys()))

                for gene1, gene2 in zip(geneList1, geneList2):
                    # if gene1 in present_genes[sample][sample2]:
                    fileList.append([[combinationOfMatrix, matrix1, chromosome1, gene1],[combinationOfMatrix, matrix2, chromosome2, gene2]])


    elif args.fileType == 'differential':

        for plotGroup in interactionFileList:
            differential_group = []

            if plotGroup[0][0] in keys_file:
                matrix_object = fileHDF5Object[plotGroup[0][0]]
                if plotGroup[1][0] in matrix_object:

                    matrix1_object = matrix_object[plotGroup[1][0]]
                    if plotGroup[1][1] in matrix1_object:
                        chromosome_object = matrix1_object[plotGroup[1][1]]
                
                        if plotGroup[1][2] in chromosome_object:
                            differential_group = [plotGroup[0][0], plotGroup[1][0], plotGroup[1][1], plotGroup[1][2]]
            log.debug('differential_group {}'.format(differential_group))
            fileList.append(differential_group)


    fileHDF5Object.close()

    filesPerThread = len(fileList) // args.threads
    # highlightSignificantRegionsFileListThread = len(highlightSignificantRegionsFileList) // args.threads

    all_data_collected = False
    stringIO_data = [None] * args.threads
    file_name_list = [None] * args.threads

    queue = [None] * args.threads
    process = [None] * args.threads
    thread_done = [False] * args.threads
    fail_flag = False
    fail_message = ''

    for i in range(args.threads):

        if i < args.threads - 1:
            fileListPerThread = fileList[i * filesPerThread:(i + 1) * filesPerThread]
        else:
            fileListPerThread = fileList[i * filesPerThread:]
        queue[i] = Queue()

        process[i] = Process(target=exportData, kwargs=dict(
            pFileList=fileListPerThread,
            pArgs=args,
            pViewpointObj=viewpointObj,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(args.threads):
            if queue[i] is not None and not queue[i].empty():
                return_content = queue[i].get()
                if 'Fail:' in return_content:
                    fail_flag = True
                    fail_message = return_content[6:]
                    log.debug('fail flag')
                else:
                    file_name_list[i], stringIO_data[i] = return_content
                    log.debug('return content')

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
    if fail_flag:
        log.error(fail_message)
        exit(1)

    stringIO_data = [item for sublist in stringIO_data for item in sublist]
    file_name_list = [item for sublist in file_name_list for item in sublist]
    
    with tarfile.open(args.outFileName, "w:gz") as tar:
        for i, file_content_string in enumerate(stringIO_data):
            # with closing(bufferObject) as fobj:

            tar_info = tarfile.TarInfo(name=file_name_list[i])
            tar_info.mtime=time.time()
            file_content_string = file_content_string.encode('utf-8')
            tar_info.size=len(file_content_string)
            file = io.BytesIO(file_content_string)
            tar.addfile(tarinfo=tar_info, fileobj=file)
               