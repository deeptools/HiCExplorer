from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
import logging
log = logging.getLogger(__name__)
import numpy as np
from scipy.sparse import csr_matrix, save_npz


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

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def calculateViewpointRange(pHiCMatrix, pViewpoint, pRange):
    '''
    This function computes the correct start and end position of a viewpoint given the reference and the range.
    '''

    max_length = pHiCMatrix.getBinPos(pHiCMatrix.getChrBinRange(pViewpoint[0])[1] - 1)[2]
    bin_size = pHiCMatrix.getBinSize()
    _range = [pRange[0], pRange[1]]
    region_start = int(pViewpoint[1]) - pRange[0]
    if region_start < 0:
        region_start = 0
        _range[0] = int(pViewpoint[1])

    region_end = int(pViewpoint[2]) + pRange[1]
    if region_end > max_length:
        # -1 is important, otherwise self.hicMatrix.getRegionBinRange will crash
        region_end = max_length - 1
        _range[1] = (max_length - int(pViewpoint[2])) + bin_size
    return region_start, region_end, _range


def getBinIndices(pHiCMatrix, pViewpoint):
    return pHiCMatrix.getRegionBinRange(pViewpoint[0], pViewpoint[1], pViewpoint[2])


def calculateViewpointRangeBins(pHiCMatrix, pViewpoint, pRange):

    viewpoint_index = getBinIndices(pHiCMatrix, pViewpoint)[0]
    start = viewpoint_index - pRange[0]
    end = viewpoint_index + pRange[1]

    return start, end


def main(args=None):

    args = parse_arguments().parse_args(args)

    hic_ma = hm.hiCMatrix(pMatrixFile=args.matrix)
    indices_values = []
    with open(args.regions, 'r') as file:
        for line in file.readlines():
            _line = line.strip().split('\t')
            if len(line) == 0:
                continue
            if len(_line) >= 2:
                chrom, start = _line[0], _line[1]

            viewpoint = (chrom, start, start)
            if args.range:
                start_range_genomic, end_range_genomic, _ = calculateViewpointRange(hic_ma, viewpoint, args.range)
                start_bin, end_bin = getBinIndices(hic_ma, (chrom, start_range_genomic, end_range_genomic))
            else:
                start_bin, end_bin = calculateViewpointRangeBins(hic_ma, viewpoint, args.rangeInBins)
            indices_values.append([start_bin, end_bin])

    if args.range:
        dimensions_new_matrix = (args.range[0] // hic_ma.getBinSize()) + (args.range[1] // hic_ma.getBinSize())
    elif args.rangeInBins:
        dimensions_new_matrix = args.rangeInBins[0] + args.rangeInBins[1]
    summed_matrix = csr_matrix((dimensions_new_matrix, dimensions_new_matrix), dtype=np.float32)
    max_length = hic_ma.matrix.shape[1]
    for start, end in indices_values:
        _start = 0
        _end = summed_matrix.shape[1]
        if start < 0:
            _start = np.absolute(start)
            start = 0
        if end >= max_length:
            _end = end
            end = max_length

        summed_matrix[_start:_end, _start:_end] += hic_ma.matrix[start:end, start:end]

    summed_matrix /= len(indices_values)

    save_npz(args.outFileName, summed_matrix)
