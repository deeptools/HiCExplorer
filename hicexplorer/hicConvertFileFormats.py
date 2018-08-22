from __future__ import division
import argparse
# from hicexplorer import HiCMatrix
from hicexplorer._version import __version__
# from hicexplorer.hicMergeMatrixBins import merge_bins
# import numpy as np
import sys
from hic2cool import hic2cool_convert
import logging
log = logging.getLogger(__name__)

from .lib import MatrixFileHandler


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Conversion of Hi-C matrices of different file formats to cool.',
        add_help=False)

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrices', '-m',
                                help='input file(s). Could be one or many files. ',
                                nargs='+',
                                required=True)

    # parserRequired.add_argument('--outFileName', '-o',
    #                             help='File name to save the exported matrix.',
    #                             required=True)

    parserRequired.add_argument('--inputFormat',
                                help='file format of the matrix file. \n'
                                'The following options are available: `h5` (native HiCExplorer '
                                'format based on hdf5 storage format), '
                                ' `cool`, `hic`, `homer`, `hicpro`',
                                choices=['h5', 'cool', 'hic', 'homer', 'hicpro'],
                                default='h5',
                                required=True)

    parserRequired.add_argument('--outputFormat',
                                help='Output format. The following options are available: `h5` (native HiCExplorer '
                                'format based on hdf5 storage format). '
                                ' `cool` and `ginteractions`',
                                default='cool',
                                choices=['cool', 'h5', 'ginteractions'],
                                required=True)

    # parserRequired.add_argument("--modus", "-mo",
    #                             choices=['resolution', 'combineSample', 'resolutionAndCombineSample', 'hic2cool'],
    #                             default='resolution',
    #                             help="Store different sample in one mcool file.")
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument("--removeCorrection", "-rc",
                           action='store_true',
                           help="Do not apply correction factors and store original data. Option only for cool input files.")
    # parserOpt.
    parserOpt.add_argument("--resolutions", '-r',
                           nargs='+',
                           help='List of resolutions that should be added.')
    parserOpt.add_argument("--help", "-h", action="help",
                           help='A conversion is only possible from '
                           'hic to cool, homer to [cool, h5], h5 to cool or cool to h5. '
                           )

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    parserOpt.add_argument('--bedFileHicpro', '-bf',
                           help='Bed file(s) of hicpro file format.',
                           nargs='+',
                           required=False)
    return parser


def main(args=None):
    log.debug(args)
    args = parse_arguments().parse_args(args)

    # parse from hicpro, homer, h5 and hic to cool

    if args.inputFormat == 'hic' and args.outputFormat == 'cool':
        log.info('Converting with hic2cool.')
        for matrix in args.matrices:
            hic2cool_convert(matrix, args.outFileName, 0)
        return
    elif args.inputFormat in ['hicpro', 'homer', 'h5']:  # and args.outputFormat in ['cool':
        if args.inputFormat == 'hicpro':
            if len(args.matrices) != len(args.bedFileHicpro):
                log.error('Number of matrices and associated bed files need to be the same.')
                log.error('Matrices: {}; Bed files: {}'.format(len(args.matrices), len(args.bedFileHicpro)))
                sys.exit(1)

        for i, matrix in enumerate(args.matrices):
            if args.inputFormat == 'hicpro':
                matrixFileHandlerInput = MatrixFileHandler(pFileType=args.inputFormat, pMatrixFile=matrix,
                                                           pBedFileHicPro=args.bedFileHicpro[i])
            else:
                matrixFileHandlerInput = MatrixFileHandler(pFileType=args.inputFormat, pMatrixFile=matrix)

            _matrix, cut_intervals, nan_bins, \
                correction_factors, distance_counts = matrixFileHandlerInput.load()

            matrixFileHandlerOutput = MatrixFileHandler(pFileType=args.outputFormat)

            matrixFileHandlerOutput.set_matrix_variables(_matrix, cut_intervals, nan_bins,
                                                         correction_factors, distance_counts)
            log.debug('Setting done')

            matrixFileHandlerOutput.save(matrix + '.' + args.outputFormat, pSymmetric=True, pApplyCorrection=False)

    # create hiC matrix with given input format
    # additional file needed for lieberman format

    # args.removeCorrection = not args.removeCorrection
    # if args.inputFormat in ['h5', 'cool', 'homer'] and args.outputFormat in ['mcool']:
    #     # create mcool file

    #     # option 1: create out of n files one mcooler, either naming or resolution differs
    #     # option 2: create out of n files one mcooler, with naming and different resolutions
    #     # option 3: create out of one file one mcooler

    #     # option 1
    #     if args.modus == 'resolution':
    #         if len(args.matrix) > 1:
    #             log.error('Please provide only one matrix')
    #             return
    #         if args.resolutions:
    #             hic_matrix = HiCMatrix.hiCMatrix(matrix, pApplyCorrectionCooler=args.removeCorrection)
    #             resolution = hic_matrix.getBinSize()
    #             for resolution_ in args.addResolutions:
    #                 merged_matrix = merge_bins(hic_matrix, float(resolution_) / resolution)
    #                 merged_matrix.save_cooler(args.outFileName + '::/resolutions/' + str(resolution_))
    #         else:
    #             log.error('Please define --resolutions')
    #             return
    #     elif args.modus == 'combineSample':
    #         if len(args.matrix) < 2:
    #             log.error('Please provide more than one matrix')
    #             return

    #         for matrix in args.matrix:
    #             hic_matrix = HiCMatrix.hiCMatrix(matrix, pApplyCorrectionCooler=args.removeCorrection)
    #             hic_matrix.save_cooler(args.outFileName + '::/samples/' + matrix)

    #     elif args.modus == 'resolutionAndCombineSample':
    #         if len(args.matrix) < 2:
    #             log.error('Please provide more than one matrix')
    #             return

    #         for matrix in args.matrix:
    #             hic_matrix = HiCMatrix.hiCMatrix(matrix, pApplyCorrectionCooler=args.removeCorrection)
    #             resolution = hic_matrix.getBinSize()
    #             for resolution_ in args.resolutions:
    #                 merged_matrix = merge_bins(hic_matrix, float(resolution_) / resolution)
    #                 merged_matrix.save_cooler(args.outFileName + '::/samples/' + matrix'/resolutions/' + str(resolution_))

    #     elif args.modus == 'hic2cool':
    #         if (args.inputFormat == 'hic'):
    #             log.info('Converting with hic2cool.')
    #             for matrix in args.matrix:
    #                 hic2cool_convert(matrix, args.outFileName, 0)
    #             return

    # elif args.inputFormat in ['h5', 'cool', 'homer'] and args.outputFormat in ['cool']:
    #     hic_matrix = HiCMatrix.hiCMatrix(args.matrix[0], pApplyCorrectionCooler=args.removeCorrection)
    #     hic_matrix.save_cooler(args.outFileName)
