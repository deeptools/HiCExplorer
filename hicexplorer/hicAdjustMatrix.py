from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicmatrix.HiCMatrix import check_cooler
import numpy as np
import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""

""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The matrix to adjust. '
                                'HiCExplorer supports the following file formats: h5 (native HiCExplorer format) '
                                'and cool.',
                                required=True)
    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the adjusted matrix.',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--chromosomes', '-c',
                           nargs='+',
                           help='List of chromosomes to keep / remove')
    parserOpt.add_argument('--action',
                           help='Keep, remove or mask the list of specified chromosomes / regions ',
                           default='keep',
                           choices=['keep', 'remove', 'mask']
                           )
    parserOpt.add_argument('--regions', '-r',
                           help='BED file which stores a list of regions to keep / remove')
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)

    if args.chromosomes is not None and args.regions is not None:
        log.error('Please specify either --chromosomes or --regions.')
        exit(1)
    hic_ma = None
    if args.chromosomes:

        if check_cooler(args.matrix) and len(args.chromosomes) == 1 and args.action == 'keep':
            hic_ma = hm.hiCMatrix(args.matrix, pChrnameList=args.chromosomes)
        else:
            hic_ma = hm.hiCMatrix(args.matrix)

        if args.action == 'keep':
            hic_ma.reorderChromosomes(args.chromosomes)
        elif args.action == 'remove':
            chromosomes = list(hic_ma.chrBinBoundaries)
            for chromosome in args.chromosomes:
                if chromosome in chromosomes:
                    chromosomes.remove(chromosome)
            hic_ma.reorderChromosomes(chromosomes)
        elif args.action == 'mask':
            hic_ma.maskChromosomes(args.chromosomes)
    elif args.regions:
        hic_ma = hm.hiCMatrix(args.matrix)
        genomic_regions = []
        with open(args.regions, 'r') as file:
            for line in file.readlines():
                _line = line.strip().split('\t')
                if len(line) == 0:
                    continue
                if len(_line) == 3:
                    chrom, start, end = _line[0], _line[1], int(_line[2]) - 1

                genomic_regions.append((chrom, start, end))

        # log.debug('genomic_regions {}'.format(genomic_regions))
        matrix_indices_regions = []
        for region in genomic_regions:
            _regionBinRange = hic_ma.getRegionBinRange(region[0], region[1], region[2])
            if _regionBinRange is not None:
                start, end = _regionBinRange
                matrix_indices_regions.extend(list(range(start, end)))

        # log.debug('matrix_indices_regions {}'.format(matrix_indices_regions))
        if args.action == 'keep':
            hic_ma.reorderBins(matrix_indices_regions)
        elif args.action == 'mask':
            hic_ma.maskBins(matrix_indices_regions)

        elif args.action == 'remove':

            full_matrix_range = np.array(range(0, max(hic_ma.matrix.shape[0], hic_ma.matrix.shape[1])))
            matrix_indices_regions = np.array(matrix_indices_regions)
            full_matrix_range[matrix_indices_regions] = -1
            mask = full_matrix_range != -1
            full_matrix_range = full_matrix_range[mask]

            hic_ma.reorderBins(full_matrix_range)
    else:
        log.info('No data to adjust given. Please specify either --chromosomes or --region parameter.')

    if hic_ma is not None:
        hic_ma.save(args.outFileName)
