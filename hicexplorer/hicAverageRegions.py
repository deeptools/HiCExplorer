import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
import logging
log = logging.getLogger(__name__)
import numpy as np
from scipy.sparse import csr_matrix, save_npz, lil_matrix


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
       Sums Hi-C contacts around given reference points and computes their average. This tool is useful to detect differences at certain reference points as for example TAD boundaries between samples.

WARNING: This tool can only be used with fixed bin size Hi-C matrices. No guarantees how and if it works on restriction site interaction matrices.
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The matrix to use for the average of TAD regions.',
                                required=True)
    parserRequired.add_argument('--regions', '-r',
                                help='BED file which stores a list of regions that are summed and averaged',
                                required=True)
    parserMutuallyExclusiveGroup = parser.add_mutually_exclusive_group(required=True)
    parserMutuallyExclusiveGroup.add_argument('--range', '-ra',
                                              help='Range of region up- and downstream of each region to include in genomic units.',
                                              nargs=2,
                                              type=int)
    parserMutuallyExclusiveGroup.add_argument('--rangeInBins', '-rib',
                                              help='Range of region up- and downstream of each region to include in bin units.',
                                              nargs=2,
                                              type=int)
    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the average regions TADs matrix.',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--coordinatesToBinMapping', '-cb',
                           help='If the region contains start and end coordinates, define if the start, center (start + (end-start) / 2) or end bin should be used as start for range.'
                           'This parameter is only important to set if the given start and end coordinates are not in the same bin'
                           ' (Default: %(default)s).',
                                choices=['start', 'center', 'end'],
                                default='start')
    parserOpt.add_argument('--considerStrandDirection',
                           help='This parameter specifies if the strand information is taken into account for the aggregation. '
                           'It has the effect that the contacts of a reverse strand region are inverted e.g. [1,2,3] becomes [3,2,1].',
                           action='store_true',
                           required=False)
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def calculateViewpointRange(pHiCMatrix, pViewpoint, pRange, pCoordinatesToBinMapping):
    '''
    This function computes the correct start and end position of a viewpoint given the reference and the range.
    '''
    start_out_of_range = False
    end_out_of_range = False
    max_length = pHiCMatrix.getBinPos(pHiCMatrix.getChrBinRange(pViewpoint[0])[1] - 1)[2]
    # bin_size = pHiCMatrix.getBinSize()
    # _range = [pRange[0], pRange[1]]

    if pCoordinatesToBinMapping == 'start':
        region_start = int(pViewpoint[1]) - pRange[0]
        region_end = int(pViewpoint[1]) + pRange[1]
    elif pCoordinatesToBinMapping == 'end':
        region_start = int(pViewpoint[2]) - pRange[0]
        region_end = int(pViewpoint[2]) + pRange[1]
    elif pCoordinatesToBinMapping == 'center':
        viewpoint_center_value = int(float(pViewpoint[1]) + ((float(pViewpoint[2]) - float(pViewpoint[1])) / 2))
        region_start = viewpoint_center_value - pRange[0]
        region_end = viewpoint_center_value + pRange[1]

    if region_start < 0:
        region_start = 0
        start_out_of_range = True

    if region_end > max_length:
        # -1 is important, otherwise self.hicMatrix.getRegionBinRange will crash
        region_end = max_length - 1
        end_out_of_range = True
    return region_start, region_end, start_out_of_range, end_out_of_range


def getBinIndices(pHiCMatrix, pViewpoint):

    return pHiCMatrix.getRegionBinRange(pViewpoint[0], pViewpoint[1], pViewpoint[2])


def calculateViewpointRangeBins(pHiCMatrix, pViewpoint, pRange, pCoordinatesToBinMapping):
    # if pCoordinatesToBinMapping == 'start_end':
    #     viewpoint_index_start = getBinIndices(pHiCMatrix, pViewpoint)[0]
    #     viewpoint_index_end = getBinIndices(pHiCMatrix, pViewpoint)[1]
    start_out_of_range = False
    end_out_of_range = False
    if pCoordinatesToBinMapping == 'start':
        viewpoint_index = getBinIndices(pHiCMatrix, pViewpoint)[0]
    elif pCoordinatesToBinMapping == 'end':
        viewpoint_index = getBinIndices(pHiCMatrix, pViewpoint)[1]
    else:
        viewpoint_center_value = int(float(pViewpoint[1]) + ((float(pViewpoint[2]) - float(pViewpoint[1])) / 2))
        viewpoint_center = [pViewpoint[0], viewpoint_center_value, viewpoint_center_value]
        viewpoint_index = getBinIndices(pHiCMatrix, viewpoint_center)[1]

    # if pCoordinatesToBinMapping == 'start_end':
    #     start = viewpoint_index_start - pRange[0]
    #     end = viewpoint_index_end + pRange[1]
    # else:
    first_bin, last_bin = pHiCMatrix.getChrBinRange(pViewpoint[0])
    start = viewpoint_index - pRange[0]
    end = viewpoint_index + pRange[1]
    if start < first_bin:
        start = first_bin
        start_out_of_range = True
    if end > last_bin:
        end = last_bin
        end_out_of_range = True
    return start, end, start_out_of_range, end_out_of_range


def main(args=None):

    args = parse_arguments().parse_args(args)

    hic_ma = hm.hiCMatrix(pMatrixFile=args.matrix)
    indices_values = []

    with open(args.regions, 'r') as file:
        for line in file.readlines():
            _line = line.strip().split('\t')
            if len(line) == 0:
                continue
            if len(_line) == 2:
                chrom, start = _line[0], _line[1]

                viewpoint = (chrom, start, start)
            elif len(_line) >= 3:
                chrom, start, end = _line[0], _line[1], _line[2]
                if args.considerStrandDirection and len(_line) < 6:
                    log.error('Strand orientation should be considered but file does not contain the 6th column of the bed file containing this information. Exiting!')
                    exit(1)

                viewpoint = (chrom, start, end)
            if args.range:
                start_range_genomic, end_range_genomic, start_out, end_out = calculateViewpointRange(hic_ma, viewpoint, args.range, args.coordinatesToBinMapping)
                start_bin, end_bin = getBinIndices(hic_ma, (chrom, start_range_genomic, end_range_genomic))
            else:
                start_bin, end_bin, start_out, end_out = calculateViewpointRangeBins(hic_ma, viewpoint, args.rangeInBins, args.coordinatesToBinMapping)
            if args.considerStrandDirection:
                indices_values.append([start_bin, end_bin, start_out, end_out, _line[5]])

            else:
                indices_values.append([start_bin, end_bin, start_out, end_out, None])

    if args.range:
        dimensions_new_matrix = (args.range[0] // hic_ma.getBinSize()) + (args.range[1] // hic_ma.getBinSize())
    elif args.rangeInBins:
        dimensions_new_matrix = args.rangeInBins[0] + args.rangeInBins[1]

    summed_matrix = lil_matrix((dimensions_new_matrix, dimensions_new_matrix), dtype=np.float32)
    count_matrix = np.zeros(shape=(dimensions_new_matrix, dimensions_new_matrix))

    # max_length = hic_ma.matrix.shape[1]
    for start, end, start_out, end_out, orientation in indices_values:
        _start = 0
        _end = summed_matrix.shape[1]
        # if start < 0:
        #     _start = np.absolute(start)
        #     start = 0
        # if end >= max_length:
        #     _end = end
        #     end = max_length
        orig_matrix_length = end - start
        if start_out:
            _start = _end - orig_matrix_length
        if end_out:
            _end = start + orig_matrix_length
        submatrix = hic_ma.matrix[start:end, start:end]
        if summed_matrix.shape != submatrix.shape:
            log.warning('Shape of a submatrix does not match. It is ignored.')
            log.warning('Region: {}'.format(hic_ma.getBinPos(start)))
            continue
        count_matrix[_start:_end, _start:_end] += 1

        if orientation is None or orientation == '+':
            summed_matrix[_start:_end, _start:_end] += hic_ma.matrix[start:end, start:end]
        elif orientation == '-':

            summed_matrix[_start:_end, _start:_end] += hic_ma.matrix[start:end, start:end].T
    summed_matrix /= count_matrix
    summed_matrix = np.array(summed_matrix)
    data = summed_matrix[np.nonzero(summed_matrix)]
    row = np.nonzero(summed_matrix)[0]
    col = np.nonzero(summed_matrix)[1]
    summed_matrix = csr_matrix((data, (row, col)), shape=(dimensions_new_matrix, dimensions_new_matrix))
    save_npz(args.outFileName, summed_matrix)
