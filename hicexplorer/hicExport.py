from __future__ import division
import argparse
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
import numpy as np

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Conversion of Hi-C matrices between different file formats.',
        add_help=False)

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--inFile', '-in',
                                help='input file(s). Could be one or many files. '
                                'Multiple input files are allowed for hicexplorer or lieberman format. '
                                ' In case of multiple input files, they will be combined. ',
                                nargs='+',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the exported matrix. In the case of "lieberman" '
                                'output format this should be the path of a folder where the information '
                                'per chromosome is stored.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--inputFormat',
                           help='file format for the matrix file. \n'
                           'The following options are available: `hicexplorer` or `h5` (native HiCExplorer '
                           'format based on hdf5 storage format), '
                           '`npz` (format used by earlier versions of HiCExplorer), '
                           '`dekker` (matrix format used in Job Dekker publications), '
                           '`lieberman` (format used by Erez Lieberman Aiden) and '
                           ' `cool`. This last formats may change '
                           'in the future.',
                           choices=['dekker', 'ren', 'lieberman', 'h5',
                                    'npz', 'GInteractions', 'cool', 'hicexplorer'],
                           default='hicexplorer')

    parserOpt.add_argument('--outputFormat',
                           help='Output format. The possibilities are "hicexplorer" or "h5" (native HiCExplorer format), '
                           '"dekker",  "ren",  '
                           'npz (former hicexplorer format), "GInteractoins" and "cool". '
                           'The dekker format outputs the whole matrix where the '
                           'first column and first row are the bin widths and labels. '
                           'The "ren" format is a list of tuples of the form '
                           'chrom, bin_star, bin_end, values. '
                           'The lieberman format writes separate files for each chromosome,'
                           'with three columns: contact start, contact end, and raw observed score. '
                           'This corresponds to the RawObserved files from lieberman group. The '
                           'hicexplorer format stores the data using a hdf5 format. Optionally, '
                           'the numpy npz format can be used for small datasets (< 4GB).'
                           'The GInteractions format is in the form : Bin1, Bin2 , Interaction, '
                           'where Bin1 and Bin2 are intervals (chr,start,end), seperated by tab.',
                           default='dekker',
                           choices=['dekker', 'ren', 'lieberman', 'h5', 'npz', 'GInteractions', 'cool', 'hicexplorer'])

    parserOpt.add_argument('--chrNameList',
                           help='list of chromosome names (only if input format is lieberman), eg : 1 2 .',
                           nargs='+',
                           )

    parserOpt.add_argument('--chromosomeOrder',
                           help='Chromosomes and order in which the chromosomes should be saved. If not all chromosomes '
                           'are given, the missing chromosomes are left out. For example, --chromosomeOrder chrX will '
                           'export a matrix only containing chromosome X.',
                           nargs='+')

    parserOpt.add_argument('--bplimit', '-b',
                           help='When merging many matrices : maximum limit (in base pairs) after '
                           'which the matrix will be truncated. i.e. TADs bigger than this '
                           'size will not be shown. For Matrices with very high resolution, '
                           'truncating the matrix after a limit helps in saving memory '
                           'during processing, without much loss of data. You can use '
                           'bplimit of 2 x size of biggest expected TAD.',
                           type=int,
                           metavar='INT bp',
                           default=None)

    parserOpt.add_argument('--clearMaskedBins',
                           help='if set, masked bins are removed from the matrix. Masked bins '
                           'are those that do not have any values, mainly because they are '
                           'repetitive regions of the genome',
                           action='store_true')

    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
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
            limit = bplimit // hic.getBinSize()
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


def main(args=None):
    log.debug(args)
    args = parse_arguments().parse_args(args)
    are_chrom_reordered = False
    # create hiC matrix with given input format
    # additional file needed for lieberman format
    if args.inputFormat == 'lieberman':
        if args.chrNameList is None:
            log.error("Error: --chrNameList is required when the input format is lieberman.")
            exit()
        else:
            hic_ma = hm.hiCMatrix(matrixFile=args.inFile, file_format='lieberman', chrnameList=args.chrNameList)

    elif args.inputFormat in ['npz', 'hicexplorer', 'h5'] and len(args.inFile) > 1:  # assume hicexplorer_multi format
        if args.bplimit:
            log.info("\nCutting maximum matrix depth to {} for saving\n".format(args.bplimit))

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
        if args.inputFormat == 'cool' and args.chromosomeOrder is not None and len(args.chromosomeOrder) == 1:
            # We have to use == 1 because we can only use the benefits of the cooler format to load the matrix partial
            # if we load one chromosome. More are so far not possible.
            hic_ma = hm.hiCMatrix(matrixFile=args.inFile[0], file_format=args.inputFormat, chrnameList=args.chromosomeOrder)
            are_chrom_reordered = True
        else:
            hic_ma = hm.hiCMatrix(matrixFile=args.inFile[0], file_format=args.inputFormat)

        if args.bplimit:
            from scipy.sparse import triu
            log.info("\nCutting maximum matrix depth to {} for saving\n".format(args.bplimit))

            limit = args.bplimit // hic_ma.getBinSize()
            hic_ma.matrix = (triu(hic_ma.matrix, k=-limit) - triu(hic_ma.matrix, k=limit)).tocsr()
            hic_ma.matrix.eliminate_zeros()

    if args.chromosomeOrder and are_chrom_reordered is False:
        hic_ma.keepOnlyTheseChr(args.chromosomeOrder)

    if args.clearMaskedBins:
        hic_ma.maskBins(hic_ma.nan_bins)

    if not args.outFileName.endswith(args.outputFormat):
        args.outFileName += "."
        args.outFileName += args.outputFormat

    if args.outputFormat == 'dekker':
        log.info('saving as dekker...')
        hic_ma.save_dekker(args.outFileName)
    elif args.outputFormat == 'ren':
        log.info('saving as ren...')
        hic_ma.save_bing_ren(args.outFileName)
    elif args.outputFormat == 'lieberman':
        log.info('saving as lieberman...')
        hic_ma.save_lieberman(args.outFileName)
    elif args.outputFormat == 'npz':
        log.info('saving as npz...')
        hic_ma.save_npz(args.outFileName)
    elif args.outputFormat == 'GInteractions':
        log.info('saving as GInteractions...')
        hic_ma.save_GInteractions(args.outFileName)
    elif args.outputFormat == 'cool':
        log.info('saving as cool...')
        hic_ma.save_cooler(args.outFileName)
    elif args.outputFormat == 'h5':
        log.info('saving as h5...')
        hic_ma.save(args.outFileName)
    else:
        log.error("An error occurred. hicExport aborted!")
        exit()
