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

from hicexplorer import hicMergeMatrixBins
from hicmatrix import HiCMatrix


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

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the exported matrix.',
                                required=True,
                                nargs='+')

    parserRequired.add_argument('--inputFormat',
                                help='file format of the matrix file. \n'
                                'The following options are available: `h5` (native HiCExplorer '
                                'format based on hdf5 storage format), '
                                ' `cool`, `hic`, `homer`, `hicpro`',
                                choices=['h5', 'cool', 'hic',
                                         'homer', 'hicpro'],
                                required=True)

    parserRequired.add_argument('--outputFormat',
                                help='Output format. The following options are available: `h5` (native HiCExplorer '
                                'format based on hdf5 storage format). '
                                ' `cool` and `ginteractions`',
                                default='cool',
                                choices=['cool', 'h5',
                                         'ginteractions', 'mcool'],
                                required=True)

    # parserRequired.add_argument("--modus", "-mo",
    #                             choices=['resolution', 'combineSample', 'resolutionAndCombineSample', 'hic2cool'],
    #                             default='resolution',
    #                             help="Store different sample in one mcool file.")
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--correction_name',
                           help='Name of the column which stores the correction factors. The information about the '
                                'column names can be figured out with the tool hicInfo. Option only for cool input files.',
                           default='weight')
    parserOpt.add_argument('--correction_division',
                           help='If set, division is applied for correction. Default is a multiplication. Option only for cool input files.',
                           action='store_true')
    parserOpt.add_argument('--store_applied_correction',
                           help='Store the applied correction and do not set correction factors. Option only for cool input files.',
                           action='store_true')
    parserOpt.add_argument('--chromosome',
                           help='Load only one chromosome. Option only for cool input files.')
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
    args = parse_arguments().parse_args(args)
    log.debug(args)

    # parse from hicpro, homer, h5 and hic to cool
    if args.inputFormat != 'hic' and args.outputFormat != 'mcool':
        if len(args.matrices) != len(args.outFileName):
            log.error(
                'Number of input matrices does not match number output matrices!')
            exit(1)
    if args.inputFormat == 'hic' and args.outputFormat == 'cool':
        log.info('Converting with hic2cool.')
        for matrix in args.matrices:
            hic2cool_convert(matrix, args.outFileName, 0)
        return
    elif args.inputFormat in ['hicpro', 'homer', 'h5', 'cool']:
        if args.inputFormat == 'hicpro':
            if len(args.matrices) != len(args.bedFileHicpro):
                log.error(
                    'Number of matrices and associated bed files need to be the same.')
                log.error('Matrices: {}; Bed files: {}'.format(
                    len(args.matrices), len(args.bedFileHicpro)))
                sys.exit(1)

        for i, matrix in enumerate(args.matrices):
            if args.inputFormat == 'hicpro':
                matrixFileHandlerInput = MatrixFileHandler(pFileType=args.inputFormat, pMatrixFile=matrix,
                                                           pBedFileHicPro=args.bedFileHicpro[i])
            else:
                correction_operator = None

                if args.correction_division:
                    correction_operator = '/'

                chromosomes_to_load = None
                if args.chromosome:
                    chromosomes_to_load = [args.chromosome]
                matrixFileHandlerInput = MatrixFileHandler(pFileType=args.inputFormat, pMatrixFile=matrix,
                                                           pCorrectionFactorTable=args.correction_name,
                                                           pCorrectionOperator=correction_operator,
                                                           pChrnameList=chromosomes_to_load)

            _matrix, cut_intervals, nan_bins, \
                correction_factors, distance_counts = matrixFileHandlerInput.load()

            log.debug('Setting done')

            if args.outputFormat in ['cool', 'h5']:
                matrixFileHandlerOutput = MatrixFileHandler(
                    pFileType=args.outputFormat)
                if args.store_applied_correction:
                    correction_factors = None
                matrixFileHandlerOutput.set_matrix_variables(_matrix, cut_intervals, nan_bins,
                                                             correction_factors, distance_counts)
                matrixFileHandlerOutput.save(
                    args.outFileName[i] + '.' + args.outputFormat, pSymmetric=True, pApplyCorrection=False)
            elif args.outputFormat in ['mcool']:
                log.debug('outformat is mcool')
                if args.resolutions and len(args.matrices) > 1:
                    log.error(
                        'Please define either one matrix and many resolutions which should be created.')
                if args.resolutions:
                    log.info(
                        'Correction factors are removed. They are not valid for any new created resolution')
                    hic_matrix = HiCMatrix.hiCMatrix()
                    hic_matrix.setMatrix(_matrix, cut_intervals)
                    bin_size = hic_matrix.getBinSize()

                    for resolution in args.resolutions:
                        _mergeFactor = int(resolution) // bin_size
                        merged_matrix = hicMergeMatrixBins.merge_bins(
                            hic_matrix, _mergeFactor)
                        matrixFileHandlerOutput = MatrixFileHandler(
                            pFileType='cool')
                        matrixFileHandlerOutput.set_matrix_variables(merged_matrix.matrix,
                                                                     merged_matrix.cut_intervals,
                                                                     merged_matrix.nan_bins,
                                                                     merged_matrix.correction_factors,
                                                                     merged_matrix.distance_counts)
                        matrixFileHandlerOutput.save(args.outFileName[0] + '.mcool' + '::/resolutions/' + str(
                            resolution), pSymmetric=True, pApplyCorrection=False)

                else:
                    hic_matrix = HiCMatrix.hiCMatrix()
                    hic_matrix.setMatrix(_matrix, cut_intervals)
                    bin_size = hic_matrix.getBinSize()
                    matrixFileHandlerOutput = MatrixFileHandler(
                        pFileType='cool')
                    matrixFileHandlerOutput.set_matrix_variables(_matrix, cut_intervals, nan_bins,
                                                                 correction_factors, distance_counts)
                    matrixFileHandlerOutput.save(args.outFileName[0] + '.mcool' + '::/resolutions/' + str(
                        bin_size), pSymmetric=True, pApplyCorrection=False)
