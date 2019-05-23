import argparse
import math
from multiprocessing import Process, Queue
import time
import logging
log = logging.getLogger(__name__)

import numpy as np
import fit_nbinom
from scipy.stats import nbinom, kstest
from statsmodels.stats.gof import gof_chisquare_discrete
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
    parserRequired.add_argument('--sparsity', '-s',
                                help='Viewpoints with a sparsity less than given are considered of bad quality.',
                                type=float,
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')
    
    parserOpt.add_argument('--averageContactBin',
                           help='Average the contacts of n bins, written to last column.',
                           type=int,
                           default=5)
            
    parserOpt.add_argument('--outFileName', '-o',
                           help='The name of the outputfile',
                           default='new_referencepoints.bed')
    parserOpt.add_argument('--threads', '-t',
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

    test_result = []
    average_count = []
    # undecided = []
    for i, referencePoint in enumerate(pReferencePoints):

        region_start, region_end, _ = pViewpointObj.calculateViewpointRange(referencePoint, (pArgs.fixateRange, pArgs.fixateRange))

        data_list = pViewpointObj.computeViewpoint(referencePoint, referencePoint[0], region_start, region_end)
        
        sparsity = (np.count_nonzero(data_list) / len(data_list))

        # fit_nb_parameters = fit_nbinom.fit(np.array(data_list))
        # fit_nb_parameters_array = [fit_nb_parameters['size'], fit_nb_parameters['prob']]
        # # np.savetxt('foo.data', np.array(data_list))
        # # log.debug('fit_nb_parameters {}'.format(fit_nb_parameters))
        # # nb_distribution = nbinom(fit_nb_parameters['size'], fit_nb_parameters['prob'])
        # test_result.append(gof_chisquare_discrete(nbinom, arg=fit_nb_parameters_array, rvs=data_list, alpha=0.1, msg=referencePoint))
        # # log.debug('test_result {}'.format(test_result))
        # exit()


        test_result.append(sparsity)
        # average_count.append(np.average(data_list))
    # with open('zero_share.txt', 'w') as file:
    #     for i, referencePoint in enumerate(pReferencePoints):
    #         file.write('{} zero share: {} {} {} {} {} \n'.format(referencePoint, test_result[i], test_result[i] > 0.95 , average_count[i], average_count[i] < 0.1, test_result[i] > 0.95 and average_count[i] < 0.1))
    #         if test_result[i] > 0.95:
    #             undecided.append('{} {} {}\n'.format(referencePoint, test_result[i]))
    # with open('undecided.txt', 'w') as file:
    #     for undecided_ in undecided:
    #         file.write(undecided_)
    pQueue.put(test_result)
    return


def main():
    args = parse_arguments().parse_args()

    viewpointObj = Viewpoint()
    referencePoints, _ = viewpointObj.readReferencePointFile(args.referencePoints)

    relative_positions = set()
    bin_size = 0

    # compute for each viewpoint the sparsity and consider these as bad with a sparsity less than given.

    referencePointsPerThread = len(referencePoints) // args.threads
    queue = [None] * args.threads
    process = [None] * args.threads
    background_model_data = None
    sparsity =  []
    
        

    for j, matrix in enumerate(args.matrices):
        sparsity_local = [None] * args.threads
        hic_ma = hm.hiCMatrix(matrix)
        viewpointObj.hicMatrix = hic_ma
        

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
                    sparsity_ = queue[i].get()
                    sparsity_local[i] = sparsity_
                    # relative_positions = relative_positions.union(relative_positions_thread)
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
        # for i in range(len(sparsity)):
        # log.debug('sparsity_local {}'.format(sparsity_local))

        sparsity_local = [item for sublist in sparsity_local for item in sublist]
        # log.debug('sparsity_local {}'.format(sparsity_local))

        sparsity.append(sparsity_local)
    # log.debug('sparsity {}'.format(sparsity))

    # change sparsity to sparsity values per viewpoint per matrix: viewpoint = [matrix1, ..., matrix_n] 
    sparsity = np.array(sparsity).T

    # log.debug('sparsity {}'.format(sparsity))
    with open(args.referencePoints, 'r') as reference_file_input:
        with open(args.outFileName + '_raw_filter', 'w') as output_file_raw:
            with open(args.outFileName, 'w') as output_file:
                for i, line in enumerate(reference_file_input.readlines()):
                    sparsity_str = '\t'.join(str(x) for x in sparsity[i])
                    # f                        sparsity_str += sparsity[j]
                    output_file_raw.write(line.strip() + '\t' + sparsity_str + '\n')
                    count = 0
                    for j in range(len(sparsity[i])):
                        if sparsity[i][j] > args.sparsity:
                            count += 1
                    if count:
                        output_file.write(line)
    # output plot of sparsity distribution per sample
    
    # re-arange values again 

    x = [[]] * len(args.matrices)
    y = [[]] *  len(args.matrices)

    for i in range(len(args.matrices)):
        y[i] = [i] * len(sparsity)
    sparsity = sparsity.T

    for i in range(len(args.matrices)):
        x[i] = sparsity[i].flatten()
    
    # log.debug(x)
    # log.debug(y)

    for i in range(len(args.matrices)):
        plt.plot(x[i], y[i],  'o', mfc='none', markersize=0.3)
    

    # plt.axvline(x=args.sparsity, color='r')
    plt.savefig(args.outFileName + '_sparsity_distribution.png', dpi=300)

    for i in range(len(args.matrices)):
        plt.hist(x[i], bins=100,  alpha=0.5)
        # plt.hist()
    plt.savefig(args.outFileName + '_sparsity_distribution_histogram.png', dpi=300)
    