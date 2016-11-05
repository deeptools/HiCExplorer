from __future__ import division
import sys
import argparse
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
import numpy as np


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Converts between different matrix file formats')

    # define the arguments
    parser.add_argument('--inFile', '-in',
                        help='input file(s). Could be one or many files. '
                        'Multiple input files are allowed for hicexplorer or lieberman format. '
                        ' In case of multiple input files, they will be combined. ',
                        nargs='+',
                        required=True)

    parser.add_argument('--inputFormat',
                        help='file format for input file. \n'
                             '(options : hicexplorer, lieberman, npz (file format of previous hicexplorer versions),'
                             'dekker.',
                        default='hicexplorer')

    parser.add_argument('--chrNameList',
                        help='list of chromosome names (only if input format is lieberman), eg : 1 2 .',
                        nargs='+',
                        )

    parser.add_argument('--outFileName', '-o',
                        help='File name to save the exported matrix. In the case of "lieberman" '
                             'output format this should be the path of a folder where the information '
                             'per chromosome is stored.',
                        required=True)

    parser.add_argument('--chromosomeOrder',
                        help='Chromosomes and order in which the chromosomes should be saved. If not all chromosomes '
                             'are given, those chromosomes are left out. For example, --chromosomeOrder chrX will '
                             'export a matrix only containing chromosome X',
                        nargs='+')

    parser.add_argument('--bplimit', '-b',
                        help='When merging many matrices : maximum limit (in base pairs) after '
                             'which the matrix will be truncated. i.e. TADs bigger than this '
                             'size will not be shown. For Matrices with very high resolution, '
                             'truncating the matrix after a limit helps in saving memory '
                             'during processing, without much loss of data. You can use '
                             'bplimit of 2 x size of biggest expected TAD. ',
                        type=int,
                        metavar='INT bp',
                        default=None)

    parser.add_argument('--outputFormat',
                        help='Output format. The possibilities are "dekker",  "ren", "hicexplorer, '
                             'npz (former hicexplorer format) and "GInteractoins". '
                             'The dekker format outputs the whole matrix where the '
                             'first column and first row are the bin widths and labels. '
                             'The "ren" format is a list of tuples of the form '
                             'chrom, bin_star, bin_end, values. '
                             'The lieberman format writes separate files for each chromosome,'
                             'with three columns : contact start, contact end, and raw observed score. '
                             'This corresponds to the RawObserved files from lieberman group. The '
                             'hicexplorer format stores the data using a hdf5 format. Optionally, '
                             'the numpy npz format can be used for small datasets (< 4GB).'
                             'The GInteractions format is in the form : Bin1, Bin2 , Interaction,'
                             'where Bin1 and Bin2 are intervals (chr,start,end), seperated by tab.',
                        default='dekker',
                        choices=['dekker', 'ren', 'lieberman', 'hicexplorer', 'npz', 'GInteractions'])

    parser.add_argument('--clearMaskedBins',
                        help='if set, masked bins are removed from the matrix. Masked bins '
                             'are those that do not have any values, mainly because they are'
                             'repetitive regions of the genome',
                        action='store_true')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def combine_matrices(matrix_list, bplimit=None):
    """
    Combines individual matrices, stored per chromosome into
    one matrix

    :param matrix_list: name of the matrices that will be combined into one.
    :param bplimit: To reduce the final file size, counts over the given distance can be removed
    :return: sparse matrix, bin intervals, nan bins, corrrections factors, distance counts
    """

    # Create empty row, col, value for the matrix
    from scipy.sparse import coo_matrix, triu
    new_cut_intervals = []
    row = np.array([]).astype("int")
    col = np.array([]).astype("int")
    values = np.array([])
    new_nan_bins = np.array([]).astype('int')
    new_correction_factors = np.array([])
    new_distance_counts = np.array([])

    # for each chr, append the row, col, value to the first one. Extend the dim
    size = 0
    for i in range(0, len(matrix_list)):
        hic = hm.hiCMatrix(matrix_list[i])

        # trim matrix if bplimit given
        if bplimit is not None:
            limit = int(bplimit / hic.getBinSize())
            matrix = (triu(hic.matrix, k=-limit) - triu(hic.matrix, k=limit)).tocoo()
        else:
            matrix = hic.matrix.tocoo()

        # add data
        row = np.concatenate([row, matrix.row + size])
        col = np.concatenate([col, matrix.col + size])
        values = np.concatenate([values, matrix.data])
        new_nan_bins = np.concatenate([new_nan_bins, hic.nan_bins + size])
        new_cut_intervals.extend(hic.cut_intervals)
        size += matrix.shape[0]

        # also add correction_factors
        if hic.correction_factors is not None:
            new_correction_factors = np.append(new_correction_factors, hic.correction_factors)
        else:
            # add an array with NaNs
            arr = np.empty(matrix.shape[0])
            arr[:] = np.NAN
            new_correction_factors = np.concatenate([new_correction_factors, arr])
        if hic.distance_counts is not None:
            new_distance_counts = np.concatenate([new_distance_counts, hic.distance_counts])

    final_mat = coo_matrix((values, (row, col)), shape=(size, size)).tocsr()

    assert len(new_cut_intervals) == final_mat.shape[0], \
        "Corrupted matrix file. Matrix size and " \
        "matrix bin definitions do not correspond"

    if len(new_distance_counts) == 0:
        new_distance_counts = None
    if len(new_correction_factors) == 0:
        new_correction_factors = None

    return final_mat, new_cut_intervals, new_nan_bins, new_correction_factors, new_distance_counts


def main():
    args = parse_arguments().parse_args()

    # create hiC matrix with given input format
    # additional file needed for lieberman format
    if args.inputFormat == 'lieberman':
        if args.chrNameList is None:
            exit("Error: --chrNameList is required when the input format is lieberman. ")
        else:
            hic_ma = hm.hiCMatrix(matrixFile=args.inFile, file_format='lieberman', chrnameList=args.chrNameList)

    elif args.inputFormat == 'npz' and len(args.inFile) > 1:  # assume hicexplorer_multi format
        if args.bplimit:
            sys.stderr.write("\nCutting maximum matrix depth to {} for saving\n".format(args.bplimit))

        matrix, cut_intervals, nan_bins, corrections_factors, distance_counts = \
            combine_matrices(args.inFile, bplimit=args.bplimit)
        hic_ma = hm.hiCMatrix()
        hic_ma.setMatrix(matrix, cut_intervals=cut_intervals)

        if len(nan_bins):
            hic_ma.nan_bins = nan_bins
        if corrections_factors is not None:
            hic_ma.correction_factors = corrections_factors
        if distance_counts is not None:
            hic_ma.distance_counts = distance_counts

    else:
        hic_ma = hm.hiCMatrix(matrixFile=args.inFile[0], file_format=args.inputFormat)
        if args.bplimit:
            from scipy.sparse import triu
            sys.stderr.write("\nCutting maximum matrix depth to {} for saving\n".format(args.bplimit))

            limit = int(args.bplimit / hic_ma.getBinSize())
            hic_ma.matrix = (triu(hic_ma.matrix, k=-limit) - triu(hic_ma.matrix, k=limit)).tocsr()
            hic_ma.matrix.eliminate_zeros()

    if args.chromosomeOrder:
        hic_ma.keepOnlyTheseChr(args.chromosomeOrder)

    if args.clearMaskedBins:
        hic_ma.maskBins(hic_ma.nan_bins)

    sys.stderr.write('saving...\n')

    if args.outputFormat == 'dekker':
        hic_ma.save_dekker(args.outFileName)
    elif args.outputFormat == 'ren':
        hic_ma.save_bing_ren(args.outFileName)
    elif args.outputFormat == 'lieberman':
        hic_ma.save_lieberman(args.outFileName)
    elif args.outputFormat == 'npz':
        hic_ma.save_npz(args.outFileName)
    elif args.outputFormat == 'GInteractions':
        hic_ma.save_GInteractions(args.outFileName)
    else:
        hic_ma.save(args.outFileName)
