from __future__ import division
import argparse
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
from multiprocessing import Process, Lock

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
                        help='Number of threads. Using the python multiprocessing module.'
                        ' One master process which is used to read the input file into the buffer and one process which is merging '
                        'the output bam files of the processes into one output bam file. All other threads do the actual computation.'
                        required=False,
                        default=4,
                        type=int
                        )
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser

def sum_matrix(pCoolerMatrix, pCoolerMatrixToAppend, pChr, pLock):

    matrix = pCoolerMatrix.matrix(balance=False, sparse=True).fetch(pChr).tocsr()
    matrix_to_append = pCoolerMatrixToAppend.matrix(balance=False, sparse=True).fetch(pChr).tocsr()
    matrix = matrix + matrix_to_append

    pCoolerMatrix.

def main(args=None):

    args = parse_arguments().parse_args(args)
    if args.format == 'cooler':
        # hic = hm.hiCMatrix(args.matrices[0], cooler_only_init=True)
        cooler_file = cooler.Cooler(pMatrixFile)
        chromosome_list = cooler_file.chromnames
        process = [None] * args.threads
        lock = Lock()
        for matrix in args.matrices[1:]:
            cooler_file_to_append = cooler.Cooler(matrix)
            chromosome_list_to_append = cooler_file_to_append.chromnames
            if chromosome_list != chromosome_list_to_append:
                exit("The two matrices have different chromosome order. Use the tool `hicExport` to change the order.\n"
                    "{}: {}\n"
                    "{}: {}".format(args.matrices[0], chromosome_list,
                                    matrix, chromosome_list_to_append))
            
            for i in range(args.threads):
                for chromosome in chromosome_list:
                
                    process[i] = Process(target=summed_matrix, kwargs=dict(
                                    pCoolerMatrix=cooler_file, 
                                    pCoolerMatrixToAppend=cooler_file_to_append, 
                                    pChr=chromosome, 
                                    pLock=lock
                                    ))
                    process[i].start()
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