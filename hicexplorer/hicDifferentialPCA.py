# get two pca files
# make correlation to find areas worth of investigation
# apply sliding window approach to find differential expressed regions

import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from io import StringIO

import cooler
from scipy.stats import pearsonr

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

    parserOpt.add_argument('--help', '-h', action='help',
                           help='show this help message and exit')

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
                data[chrom_X] = np.array([float(value)])
                chromosomeList.append(chromosome)
    return data, chromosomeList

def computeDifferentialPCA(pPCADictFile1, pPCADictFile2, pPValue, pCorrelationThreshold, pSlidingWindowSize, pQueue):

    differential_regions_per_chromosome = {}
    for key in pPCADictFile1:

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
                else
                    end_i = i + pSlidingWindowSize
            else:
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

    pca_threads = [None] * args.threads
    pcaListPerThread = len(chromosomes_list) // args.threads
    all_data_collected = False
    queue = [None] * args.threads
    process = [None] * args.threads
    thread_done = [False] * args.threads
    for i in range(args.threads):

        if i < args.threads - 1:
            pcaListPerThread = chromosomes_list[i * chromosomesListPerThread:(i + 1) * chromosomesListPerThread]
        else:
            pcaListPerThread = chromosomes_list[i * chromosomesListPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=computeDifferentialPCA, kwargs=dict(
            pPCAList=pcaListPerThread,
            pPValue=args.pValue
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(args.threads):
            if queue[i] is not None and not queue[i].empty():
                short_v_long_range_matrix_threads[i] = queue[i].get()
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