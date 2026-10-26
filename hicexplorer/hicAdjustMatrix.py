"""hicAdjustMatrix: a thin layer over the C++ executable hicAdjustMatrix.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import check_cooler
import numpy as np
import cooler
import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """The argparse parser of hicAdjustMatrix, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicAdjustMatrix')


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
        log.debug('pArgs.chromosomes {}'.format(pArgs.chromosomes))
        for chromosome in pArgs.chromosomes:
            if chromosome in chromosomes_list:
                chromosomes_list_to_operate_on.append(chromosome)
            else:
                log.warning('Chromosome not available in matrix: {} {}'.format(pArgs.matrix, chromosome))
        log.debug('85: chromosomes_list_to_operate_on {}'.format(chromosomes_list_to_operate_on))
        if len(chromosomes_list_to_operate_on) == 0:
            log.error('No valid chromosome given: {}. Available: {}'.format(pArgs.chromosomes, chromosomes_list))
            exit(1)
        if pArgs.action == 'keep':
            hic_matrix.reorderChromosomes(chromosomes_list_to_operate_on)
        elif pArgs.action == 'remove':
            log.debug('chromosomes_list {}'.format(chromosomes_list))

            for chromosome in chromosomes_list_to_operate_on:
                if chromosome in chromosomes_list:
                    chromosomes_list.remove(chromosome)
            # for chromosome in chromosomes_list:
            #     if chromosome in chromosomes_list_to_operate_on:
            #         chromosomes_list.remove(chromosome)
                    log.debug('chromosome {}'.format(chromosome))
            log.debug('chromosomes_list {}'.format(chromosomes_list))
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

        matrix_indices_regions = []
        for region in genomic_regions:
            log.debug('region {}'.format(region))
            _regionBinRange = hic_matrix.getRegionBinRange(region[0], int(region[1]), int(region[2]))
            if _regionBinRange is not None:
                start, end = _regionBinRange
                matrix_indices_regions.extend(list(range(start, end + 1)))  # end is inclusive, so +1
                if pArgs.action == 'remove':
                    chr_start, chr_end = hic_matrix.getChrBinRange(chrom)
                    if (start > chr_start) and (end < chr_end - 1):
                        log.warning("{}:{}-{} entry may generate discounted regions on a chromosome."
                                    "Please consider using `mask` action to deal with that.".format(chrom, start, end))

        if pArgs.action == 'keep':
            hic_matrix.reorderBins(matrix_indices_regions)

        elif pArgs.action == 'mask':
            hic_matrix.maskBins(matrix_indices_regions)

        elif pArgs.action == 'remove':
            hic_matrix.maskBins(matrix_indices_regions)
            hic_matrix.orig_bin_ids = []
            hic_matrix.orig_cut_intervals = []
            hic_matrix.nan_bins = []

    elif pArgs.maskBadRegions:
        if check_cooler(pArgs.matrix) and len(pArgs.chromosomes) == 1 and pArgs.action == 'keep':
            hic_matrix = hm.hiCMatrix(pArgs.matrix, pChrnameList=pArgs.chromosomes)
        else:
            hic_matrix = hm.hiCMatrix(pArgs.matrix)

    else:
        log.info('No data to adjust given. Please specify either --chromosomes or --region parameter.')

    if pArgs.interIntraHandling is not None:
        if pArgs.interIntraHandling == 'inter':
            # iterate over all given chromosomes
            chromosomes_list = list(hic_matrix.chrBinBoundaries)
            for chromosome in chromosomes_list:
                # get start and end bins of that chromosome
                start, end = hic_matrix.getChrBinRange(chromosome)
                # remove inter region
                hic_matrix.matrix[start:end, end:] = 0

        elif pArgs.interIntraHandling == 'intra':
            # iterate over all given chromosomes
            chromosomes_list = list(hic_matrix.chrBinBoundaries)
            for chromosome in chromosomes_list:
                # get start and end bins of that chromosome
                start, end = hic_matrix.getChrBinRange(chromosome)
                # remove intra region
                hic_matrix.matrix[start:end, start:end] = 0
        else:
            log.warning('Option not valid: {}'.format(pArgs.interIntraHandling))

    return hic_matrix


def main(args=None):
    """Run the C++ implementation of hicAdjustMatrix with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicAdjustMatrix')(args)
