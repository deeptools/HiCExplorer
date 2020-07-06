import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicmatrix.HiCMatrix import check_cooler
import numpy as np
import cooler
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
                                help='The Hi-C matrix to adjust. '
                                'HiCExplorer supports the following file formats: h5 (native HiCExplorer format) '
                                'and cool.',
                                required=True)
    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the adjusted matrix.',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')
    parserMutuallyExclusive = parser.add_mutually_exclusive_group()
    parserMutuallyExclusive.add_argument('--chromosomes', '-c',
                                         nargs='+',
                                         help='List of chromosomes to keep / remove.')
    parserMutuallyExclusive.add_argument('--regions', '-r',
                                         help='BED file which stores a list of regions to keep / remove.')
    parserMutuallyExclusive.add_argument('--maskBadRegions', '-mbr',
                                         help='Bad regions are identified and masked.')
    parserOpt.add_argument('--action',
                           help='Keep, remove or mask the list of specified chromosomes / regions. ',
                           default='keep',
                           choices=['keep', 'remove', 'mask']
                           )

    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def adjustMatrix(pArgs):
    if pArgs.chromosomes is not None and pArgs.regions is not None:
        log.error('Please specify either --chromosomes or --regions.')
        exit(1)
    hic_matrix = None
    if pArgs.chromosomes:

        if check_cooler(pArgs.matrix) and len(pArgs.chromosomes) == 1 and pArgs.action == 'keep':
            chromosomes_list = cooler.Cooler(pArgs.matrix).chromnames
            if pArgs.chromosomes[0] in chromosomes_list:
                hic_matrix = hm.hiCMatrix(pArgs.matrix, pChrnameList=pArgs.chromosomes)
            else:
                log.error('Chromosome not available in matrix: {} {}'.format(pArgs.matrix, pArgs.chromosomes[0]))
                exit(1)
        else:
            hic_matrix = hm.hiCMatrix(pArgs.matrix)

        chromosomes_list = list(hic_matrix.chrBinBoundaries)
        chromosomes_list_to_operate_on = []
        for chromosome in pArgs.chromosomes:
            if chromosome in chromosomes_list:
                chromosomes_list_to_operate_on.append(chromosome)
            else:
                log.warning('Chromosome not available in matrix: {} {}'.format(pArgs.matrix, chromosome))
        if len(chromosomes_list_to_operate_on) == 0:
            log.error('No valid chromosome given: {}. Available: {}'.format(pArgs.chromosomes, chromosomes_list))
            exit(1)
        if pArgs.action == 'keep':
            hic_matrix.reorderChromosomes(chromosomes_list_to_operate_on)
        elif pArgs.action == 'remove':
            # chromosomes = list(hic_matrix.chrBinBoundaries)
            for chromosome in chromosomes_list:
                if chromosome in chromosomes_list_to_operate_on:
                    chromosomes_list.remove(chromosome)
            hic_matrix.reorderChromosomes(chromosomes_list)
        elif pArgs.action == 'mask':
            hic_matrix.maskChromosomes(chromosomes_list_to_operate_on)
    elif pArgs.regions:
        hic_matrix = hm.hiCMatrix(pArgs.matrix)
        chromosomes_list = list(hic_matrix.chrBinBoundaries)
        genomic_regions = []
        with open(pArgs.regions, 'r') as file:
            for line in file.readlines():
                _line = line.strip().split('\t')
                log.debug('_line {}'.format(_line))
                if len(line) < 3:
                    log.warning("An entry shorter than 3 columns has been found!")
                    continue
                if len(_line) >= 3:
                    chrom, start, end = _line[0], int(_line[1]), int(_line[2])
                    log.debug('chrom {}'.format(chrom))
                    if chrom in chromosomes_list:
                        genomic_regions.append((chrom, start, end))
                    else:
                        log.warning('Chromosome not available in matrix, '
                                    'ignoring regions: {} {}'.format(pArgs.matrix, chrom))
        if len(genomic_regions) == 0:
            log.error('No valid chromosome given. Available: {}'.format(chromosomes_list))
            exit(1)
        log.debug('genomic_regions {}'.format(genomic_regions))
        matrix_indices_regions = []
        for region in genomic_regions:
            log.debug('region {}'.format(region))
            _regionBinRange = hic_matrix.getRegionBinRange(region[0], int(region[1]), int(region[2])-1)
            if _regionBinRange is not None:
                start, end = _regionBinRange
                matrix_indices_regions.extend(list(range(start, end)))
        log.debug('foo {}'.format(matrix_indices_regions))

        # log.debug('matrix_indices_regions {}'.format(matrix_indices_regions))
        if pArgs.action == 'keep':
            hic_matrix.reorderBins(matrix_indices_regions)

        elif pArgs.action == 'mask':
            hic_matrix.maskBins(matrix_indices_regions)

        elif pArgs.action == 'remove':

            full_matrix_range = np.array(range(0, max(hic_matrix.matrix.shape[0], hic_matrix.matrix.shape[1])))
            matrix_indices_regions = np.array(matrix_indices_regions)
            full_matrix_range[matrix_indices_regions] = -1
            mask = full_matrix_range != -1
            full_matrix_range = full_matrix_range[mask]

            hic_matrix.reorderBins(full_matrix_range)
    elif pArgs.maskBadRegions:
        if check_cooler(pArgs.matrix) and len(pArgs.chromosomes) == 1 and pArgs.action == 'keep':
            hic_matrix = hm.hiCMatrix(pArgs.matrix, pChrnameList=pArgs.chromosomes)
        else:
            hic_matrix = hm.hiCMatrix(pArgs.matrix)

    else:
        log.info('No data to adjust given. Please specify either --chromosomes or --region parameter.')

    return hic_matrix


def main(args=None):

    args = parse_arguments().parse_args(args)

    hic_matrix = adjustMatrix(args)

    if hic_matrix is not None:
        hic_matrix.save(args.outFileName)
