"""hicPCA: a thin layer over the C++ executable hicPCA.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import warnings
import argparse

from scipy.sparse import csr_matrix, lil_matrix
from scipy import linalg
from scipy.stats import pearsonr
import numpy as np
import pyBigWig

from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import obs_exp_matrix_lieberman, obs_exp_matrix_non_zero
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros
from hicexplorer.utilities import enlarge_bins
from hicexplorer.parserCommon import CustomFormatter
from hicexplorer.utilities import toString
from hicexplorer.utilities import opener
from hicmatrix.lib import MatrixFileHandler
from .readBed import ReadBed
import logging
log = logging.getLogger(__name__)

warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)


def parse_arguments(args=None):
    """The argparse parser of hicPCA, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicPCA')


def correlateEigenvectorWithGeneTrack(pMatrix, pEigenvector, pGeneTrack):
    '''
    This function correlates the eigenvectors per chromosome with the gene
    density. If the correlation is negative, the eigenvector values are
    multiplied with -1.
    '''
    log.debug('correlate eigenvector!')
    file_h = opener(pGeneTrack)
    bed = ReadBed(file_h)

    gene_occurrence = np.zeros(len(pMatrix.cut_intervals))
    gene_occurrence_per_chr = {}

    chromosome_list = pMatrix.getChrNames()

    for interval in bed:
        chromosome_name = interval.chromosome
        if chromosome_name not in chromosome_list:
            continue
        if interval.start > pMatrix.get_chromosome_sizes()[chromosome_name]:
            log.warning('Your chromosome sizes do not match the chromosome sizes of the extraTrack data!')
            log.warning('Your chromosome {}; Size {}. ExtraTrack data {} {} {}'.format(chromosome_name, pMatrix.get_chromosome_sizes()[chromosome_name], interval.chromosome,
                                                                                       interval.start, interval.end))
            log.warning('Please create your interaction matrix with a chromosome size file! However, if the sizes are intended and it is accepted that certain regions are not part of the correlation, you can ignore this message.')
            continue
        # in which bin of the Hi-C matrix is the given gene?
        bin_id = pMatrix.getRegionBinRange(interval.chromosome,
                                           interval.start, interval.end)
        if bin_id is None:
            log.warning('Your chromosome sizes do not match the chromosome sizes of the extraTrack data!')
            log.warning('Your chromosome {}; Size {}. ExtraTrack data {} {} {}'.format(chromosome_name, pMatrix.get_chromosome_sizes()[chromosome_name], interval.chromosome,
                                                                                       interval.start, interval.end))
            log.warning('Please create your interaction matrix with a chromosome size file! However, if the sizes are intended and it is accepted that certain regions are not part of the correlation, you can ignore this message.')
            continue
        # add +1 for one gene occurrence in this bin
        gene_occurrence[bin_id[1]] += 1

    for chromosome in chromosome_list:
        # where is the start and the end bin of a chromosome?
        bin_id = pMatrix.getChrBinRange(chromosome)
        gene_occurrence_per_chr[chromosome] = \
            gene_occurrence[bin_id[0]:bin_id[1]]

    # change from [[1,2], [3,4], [5,6]] to [[1,3,5],[2,4,6]]
    pEigenvector = np.array(pEigenvector).real.transpose()

    # correlate gene density and eigenvector values.
    # if positive correlation, do nothing, if negative, flip the values.
    # computed per chromosome
    for chromosome in chromosome_list:
        bin_id = pMatrix.getChrBinRange(chromosome)
        for i, eigenvector in enumerate(pEigenvector):
            _correlation = pearsonr(np.array(eigenvector[bin_id[0]:bin_id[1]]).real,
                                    gene_occurrence_per_chr[chromosome])
            if _correlation[0] < 0:
                eigenvector[bin_id[0]:bin_id[1]] = np.negative(eigenvector[bin_id[0]:bin_id[1]])
            # log.debug('correlated to {}!'.format(_correlation[0]))
    return np.array(pEigenvector).transpose()


def correlateEigenvectorWithHistonMarkTrack(pEigenvector, bwTrack, chromosome,
                                            start, end, pHistonMarkTrack,
                                            pHistonMarkType):
    """
        This function flip the signs only if both compartments exist for the
        given `chromosome`, otherwise it doesn't change the signs and
        the investigation should be done manually.
    """
    for index, vector in enumerate(pEigenvector):
        pos_indices = np.where(vector > 0)[0]
        neg_indices = np.where(vector < 0)[0]
        pos_sum = 0
        neg_sum = 0
        pos_mean = 0
        neg_mean = 0
        if chromosome in bwTrack.chroms().keys():
            for ind in pos_indices:
                if bwTrack.stats(chromosome, start[ind], end[ind])[0]:
                    pos_sum += bwTrack.stats(chromosome, start[ind],
                                             end[ind])[0]
            if pos_sum != 0:
                pos_mean = pos_sum / len(pos_indices)
            for ind in neg_indices:
                if bwTrack.stats(chromosome, start[ind], end[ind])[0]:
                    neg_sum += bwTrack.stats(chromosome, start[ind],
                                             end[ind])[0]
            if neg_sum != 0:
                neg_mean = neg_sum / len(neg_indices)
            if pHistonMarkType == 'active':
                if (pos_mean < neg_mean) and (neg_mean != 0) and (pos_mean != 0):
                    # flip the sign
                    vector[pos_indices] = np.negative(vector[pos_indices])
                    vector[neg_indices] = np.negative(vector[neg_indices])
            else:
                assert (pHistonMarkType == 'inactive')
                if (pos_mean > neg_mean) and (neg_mean != 0) and (pos_mean != 0):
                    # flip the sign
                    vector[pos_indices] = -1 * vector[pos_indices]
                    vector[neg_indices] = -1 * vector[neg_indices]
        pEigenvector[index] = vector


def main(args=None):
    """Run the C++ implementation of hicPCA with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicPCA')(args)
