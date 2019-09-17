# get two pca files
# make correlation to find areas worth of investigation
# apply sliding window approach to find differential expressed regions

import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from io import StringIO
from multiprocessing import Process, Queue
from scipy.stats import anderson_ksamp, ranksums
import time
import cooler
from scipy.stats import pearsonr
import numpy as np
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import toString
from hicmatrix.HiCMatrix import check_cooler
import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
Prints information about a matrix or matrices including matrix size,
number of elements, sum of elements, etc.
An example usage is:
$ hicInfo -m matrix1.h5 matrix2.h5 matrix3.h5
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--pcaFiles', '-m',
                                help='Two PCA files computed by hicPCA.',
                                nargs=2,
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--outFileName', '-o',
                           help='File name to save information of the matrix instead of writing it to the bash.'
                           )
    parserOpt.add_argument('--slidingWindowSize', '-sws',
                           type=int,
                           help='Sliding window size for differential test to detect differential expressed PCa regions',
                           default=5)
    parserOpt.add_argument('--pValue', '-p',
                           type=float,
                           help='P-vlaue for differential expression acceptance.',
                           default=0.01)
    parserOpt.add_argument('--help', '-h', action='help',
                           help='show this help message and exit')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads to use, the parallelization is implemented per chromosome.',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser



def readPCAFile(pFile):
    data = {}
    chromosomeList = []
    chromosome = None
    with open(pFile, 'rb') as file:
        for line in file.readlines():
            line = toString(line)
            fields = line.strip().split('\t')
            try:
                chrom_X, start_X, end_X, value = fields
            except Exception:
                pass
            
            if chromosome == chrom_X:
                data[chrom_X].append(float(value))
            else:
                chromosome = chrom_X
                data[chrom_X] = [float(value)]
                chromosomeList.append(chromosome)
    return data, chromosomeList

def computeDifferentialPCA(pPCADictFile1, pPCADictFile2, pChromosomeList, pPValue, pSlidingWindowSize, pQueue):

    differential_regions_per_chromosome = {}
    for key in pChromosomeList:

        # correlation = np.corrcoef(pPCADictFile1[key],
        #                         pPCADictFile2[key])
        # mask =  -pCorrelationThreshold <= correlation
        # mask2 = correlation <= pCorrelationThreshold

        # mask = np.logical_and(mask, mask2)
        start_i = None
        end_i = None

        differential_expressed_regions = []
        for i in range(0, len(pPCADictFile1[key])-pSlidingWindowSize):
        
            statistic, significance_level_test1 = ranksums(pPCADictFile1[key][i:i+pSlidingWindowSize], pPCADictFile2[key][i:i+pSlidingWindowSize])

            if significance_level_test1 <= pPValue:
                if start_i is None:

                    start_i = i
                else:
                    end_i = i + pSlidingWindowSize
            else:
                if start_i is not None and end_i is not None:
                    differential_expressed_regions.append([start_i, end_i])
                start_i = None
                end_i = None
        
        differential_regions_per_chromosome[key] = differential_expressed_regions

    pQueue.put(differential_regions_per_chromosome)

                


            # no clear correlation given, and is considered as candidate

def main(args=None):

    args = parse_arguments().parse_args(args)
    
    data_file_1, chromosomeList1 = readPCAFile(args.pcaFiles[0])
    data_file_2, chromosomeList2 = readPCAFile(args.pcaFiles[1])

    differential_expressed_threads = [None] * args.threads
    pcaListPerThread = len(chromosomeList1) // args.threads
    all_data_collected = False
    queue = [None] * args.threads
    process = [None] * args.threads
    thread_done = [False] * args.threads
    for i in range(args.threads):

        if i < args.threads - 1:
            pca1ListPerThread = chromosomeList1[i * pcaListPerThread:(i + 1) * pcaListPerThread]
            # pca2ListPerThread = chromosomeList2[i * pcaListPerThread:(i + 1) * pcaListPerThread]

        else:
            pca1ListPerThread = chromosomeList1[i * pcaListPerThread:]
            # pca2ListPerThread = chromosomeList1[i * pcaListPerThread:]

# pPCADictFile1, pPCADictFile2, pPValue, pSlidingWindowSize, pQueue
        queue[i] = Queue()
        process[i] = Process(target=computeDifferentialPCA, kwargs=dict(
            pPCADictFile1=data_file_1,
            pPCADictFile2=data_file_2,
            pChromosomeList=pca1ListPerThread,
            pPValue=args.pValue,
            pSlidingWindowSize=args.slidingWindowSize,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(args.threads):
            if queue[i] is not None and not queue[i].empty():
                differential_expressed_threads[i] = queue[i].get()
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

    
    log.debug('differential_expressed_threads {}'.format(differential_expressed_threads))
 
    with open(args.outFileName, 'w') as file:

        header  = '# Created with HiCExplorer\'s hicDifferentialPCA ' + __version__ + '\n'
        header += '# P-value: ' + args.pValue +'\n'
        file.write(header)

        for thread_data in differential_expressed_threads:
            for chromosome in thread_data:
                if len(thread_data[chromosome]) != 0:
                    for region in thread_data[chromosome]:
                        file.write(data_file_1[chromosome][region[0]], data_file_1[chromosome][region[1]], data_file_1[chromosome][region[2]])
                    # data_file_1[chromosome]
                    # data_file_2[chromosome]