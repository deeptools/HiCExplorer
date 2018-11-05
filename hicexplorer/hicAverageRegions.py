from __future__ import division
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import toString
from hicmatrix.HiCMatrix import check_cooler
import logging
log = logging.getLogger(__name__)

from scipy.sparse import csr_matrix, save_npz


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

    parserRequired.add_argument('--matrix', '-m',
                                help='The matrix (or multiple matrices) to get information about. '
                                'HiCExplorer supports the following file formats: h5 (native HiCExplorer format) '
                                'and cool.',
                                required=True)
    parserRequired.add_argument('--regions', '-r',
                                help='BED file which stores a list of regions that are summed and averaged',
                                required=True)
    parserRequired.add_argument('--range', '-ra',
                                help='BED file which stores a list of regions to keep / remove',
                                nargs=2,
                                required=True)
    parserRequired.add_argument('--outFileName', '-out',
                                help='File name to save the adjusted matrix.',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def calculateViewpointRange(self, pHiCMatrix, pViewpoint, pRange):
    '''
    This function computes the correct start and end position of a viewpoint given the viewpoint and the range.
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


def main():

    args = parse_arguments().parse_args()

    hic_ma = hm.hiCMatrix(pMatrixFile=args.matrix)
    indices_values = []

    with open(args.regions, 'r') as file:
        for line in file.readlines():
            _line = line.strip().split('\t')
            if len(line) == 0:
                continue
            if len(_line) == 3:
                chrom, start = _line[0], _line[1]

            viewpoint = (chrom, start, start)
            start_bin, end_bin = calculateViewpointRange(hic_ma, viewpoint, args.range)
            indices_values.append([start_bin, end_bin])

    dimensions_new_matrix = (args.range[0] // hic_ma.getBinSize()) + (args.range[1] // hic_ma.getBinSize())

    summed_matrix = csr_matrix((dimensions_new_matrix, dimensions_new_matrix), dtype=np.float32)

    for start, end in indices_values:
        summed_matrix += hic_ma.matrix[start:end, start:end]

    summed_matrix /= len(indices_values)

    save_npz(args.outFileName, summed_matrix)