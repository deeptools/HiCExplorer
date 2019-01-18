from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse

from scipy.sparse import csr_matrix, lil_matrix
from scipy import linalg
from scipy.stats import pearsonr
import numpy as np
import pyBigWig

from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import obs_exp_matrix_lieberman, obs_exp_matrix_norm
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros
from hicexplorer.parserCommon import CustomFormatter
from hicexplorer.utilities import toString
from hicexplorer.utilities import opener
from hicmatrix.lib import MatrixFileHandler
from .readBed import ReadBed
import logging
log = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=CustomFormatter,
        add_help=False,
        conflict_handler='resolve',
        #        usage="%(prog)s --matrix hic_matrix.h5 -o pca1.bedgraph pca2.bedgraph ",
        description="""
Computes PCA eigenvectors for a Hi-C matrix.

    $ hicPCA --matrix hic_matrix.h5 -o pca1.bedgraph pca2.bedgraph

"""
    )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='HiCExplorer matrix in h5 format.',
                                required=True)

    parserRequired.add_argument('--outputFileName', '-o',
                                help='File names for the result of the pca. Number of output file '
                                'must match the number of computed eigenvectors.',
                                nargs='+',
                                default=['pca1', 'pca2'],
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--numberOfEigenvectors', '-noe',
                           help='The number of eigenvectors that the PCA should compute.',
                           default=2,
                           type=int,
                           required=False)

    parserOpt.add_argument('--format', '-f',
                           help='output format. Either bedgraph or bigwig.',
                           choices=['bedgraph', 'bigwig'],
                           default='bigwig',
                           required=False)

    parserOpt.add_argument('--chromosomes',
                           help='List of chromosomes to be included in the '
                           'correlation.',
                           default=None,
                           nargs='+')
    parserOpt.add_argument('--norm',
                           help='Different obs-exp normalization as used by Homer software.',
                           action='store_true')
    parserOpt.add_argument('--geneTrack',
                           help='The gene track is needed to decide if the values of the eigenvector need a sign flip or not.',
                           default=None)
    parserOpt.add_argument('--pearsonMatrix', '-pm',
                           help='Internally the input matrix is converted per chromosome to obs_exp matrix and consecutively to a Pearson matrix.'
                           ' Set this parameter to write the pearson matrix to a file.'
                           )
    parserOpt.add_argument('--obsexpMatrix', '-oem',
                           help='Internally the input matrix is converted per chromosome to obs_exp matrix and consecutively to a Pearson matrix.'
                           ' Set this parameter to write the observe / expected matrix to a file.'
                           )
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def correlateEigenvectorWithGeneTrack(pMatrix, pEigenvector, pGeneTrack):
    '''
    This function correlates the eigenvectors per chromosome with the gene density.
    If the correlation is negative, the eigenvector values are multiplied with -1.
    '''

    file_h = opener(pGeneTrack)
    bed = ReadBed(file_h)

    gene_occurrence = np.zeros(len(pMatrix.cut_intervals))
    gene_occurrence_per_chr = {}

    chromosome_list = pMatrix.getChrNames()

    for interval in bed:
        chromosome_name = interval.chromosome
        if chromosome_name not in chromosome_list:
            continue
        # in which bin of the Hi-C matrix is the given gene?
        bin_id = pMatrix.getRegionBinRange(interval.chromosome, interval.start, interval.end)

        # add +1 for one gene occurrence in this bin
        gene_occurrence[bin_id[1]] += 1

    for chromosome in chromosome_list:
        # where is the start and the end bin of a chromosome?
        bin_id = pMatrix.getChrBinRange(chromosome)
        gene_occurrence_per_chr[chromosome] = gene_occurrence[bin_id[0]: bin_id[1]]

    # change from [[1,2], [3,4], [5,6]] to [[1,3,5],[2,4,6]]
    pEigenvector = np.array(pEigenvector).real.transpose()

    # correlate gene density and eigenvector values.
    # if positive correlation, do nothing, if negative, flip the values.
    # computed per chromosome
    for chromosome in chromosome_list:
        bin_id = pMatrix.getChrBinRange(chromosome)
        for i, eigenvector in enumerate(pEigenvector):
            _correlation = pearsonr(eigenvector[bin_id[0]:bin_id[1]].real, gene_occurrence_per_chr[chromosome])
            if _correlation[0] < 0:
                eigenvector[bin_id[0]:bin_id[1]] = np.negative(eigenvector[bin_id[0]:bin_id[1]])
    return np.array(pEigenvector).transpose()


def main(args=None):
    args = parse_arguments().parse_args(args)
    if int(args.numberOfEigenvectors) != len(args.outputFileName):
        log.error("Number of output file names and number of eigenvectors does not match. Please"
                  "provide the name of each file.\nFiles: {}\nNumber of eigenvectors: {}".format(args.outputFileName,
                                                                                                 args.numberOfEigenvectors))
        exit(1)

    ma = hm.hiCMatrix(args.matrix)
    ma.maskBins(ma.nan_bins)

    if args.chromosomes:
        ma.keepOnlyTheseChr(args.chromosomes)

    vecs_list = []
    chrom_list = []
    start_list = []
    end_list = []
    # PCA is computed per chromosome
    length_chromosome = 0
    chromosome_count = len(ma.getChrNames())
    if args.pearsonMatrix:
        trasf_matrix_pearson = lil_matrix(ma.matrix.shape)

    if args.obsexpMatrix:
        trasf_matrix_obsexp = lil_matrix(ma.matrix.shape)

    for chrname in ma.getChrNames():
        chr_range = ma.getChrBinRange(chrname)
        length_chromosome += chr_range[1] - chr_range[0]
    for chrname in ma.getChrNames():
        chr_range = ma.getChrBinRange(chrname)

        submatrix = ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
        if args.norm:
            obs_exp_matrix_ = obs_exp_matrix_norm(submatrix)

        else:
            obs_exp_matrix_ = obs_exp_matrix_lieberman(submatrix, length_chromosome, chromosome_count)
        obs_exp_matrix_ = convertNansToZeros(csr_matrix(obs_exp_matrix_)).todense()
        obs_exp_matrix_ = convertInfsToZeros(csr_matrix(obs_exp_matrix_)).todense()

        if args.obsexpMatrix:
            trasf_matrix_obsexp[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(obs_exp_matrix_)

        pearson_correlation_matrix = np.corrcoef(obs_exp_matrix_)
        pearson_correlation_matrix = convertNansToZeros(csr_matrix(pearson_correlation_matrix)).todense()
        pearson_correlation_matrix = convertInfsToZeros(csr_matrix(pearson_correlation_matrix)).todense()

        if args.pearsonMatrix:
            trasf_matrix_pearson[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(pearson_correlation_matrix)

        corrmatrix = np.cov(pearson_correlation_matrix)
        corrmatrix = convertNansToZeros(csr_matrix(corrmatrix)).todense()
        corrmatrix = convertInfsToZeros(csr_matrix(corrmatrix)).todense()
        evals, eigs = linalg.eig(corrmatrix)
        k = args.numberOfEigenvectors

        chrom, start, end, _ = zip(*ma.cut_intervals[chr_range[0]:chr_range[1]])
        vecs_list += eigs[:, :k].tolist()

        chrom_list += chrom
        start_list += start
        end_list += end

    if args.pearsonMatrix:
        file_type = 'cool'
        if args.pearsonMatrix.endswith('.h5'):
            file_type = 'h5'
        matrixFileHandlerOutput = MatrixFileHandler(pFileType=file_type)
        matrixFileHandlerOutput.set_matrix_variables(trasf_matrix_pearson.tocsr(),
                                                     ma.cut_intervals,
                                                     ma.nan_bins,
                                                     ma.correction_factors,
                                                     ma.distance_counts)
        matrixFileHandlerOutput.save(args.pearsonMatrix, pSymmetric=True, pApplyCorrection=False)

    if args.obsexpMatrix:
        file_type = 'cool'
        if args.obsexpMatrix.endswith('.h5'):
            file_type = 'h5'
        matrixFileHandlerOutput = MatrixFileHandler(pFileType=file_type)
        matrixFileHandlerOutput.set_matrix_variables(trasf_matrix_obsexp.tocsr(),
                                                     ma.cut_intervals,
                                                     ma.nan_bins,
                                                     ma.correction_factors,
                                                     ma.distance_counts)
        matrixFileHandlerOutput.save(args.obsexpMatrix, pSymmetric=True, pApplyCorrection=False)

    if args.geneTrack:
        vecs_list = correlateEigenvectorWithGeneTrack(ma, vecs_list, args.geneTrack)

    if args.format == 'bedgraph':
        for idx, outfile in enumerate(args.outputFileName):
            assert(len(vecs_list) == len(chrom_list))

            with open(outfile, 'w') as fh:
                for i, value in enumerate(vecs_list):
                    if len(value) == args.numberOfEigenvectors:
                        if isinstance(value[idx], np.complex):
                            value[idx] = value[idx].real
                        fh.write("{}\t{}\t{}\t{:.12f}\n".format(toString(chrom_list[i]), start_list[i], end_list[i], value[idx]))

    elif args.format == 'bigwig':
        if not pyBigWig.numpy == 1:
            log.error("ERROR: Your version of pyBigWig is not supporting numpy: {}".format(pyBigWig.__file__))
            exit(1)
        old_chrom = chrom_list[0]
        header = []
        for i, _chrom in enumerate(chrom_list):
            if old_chrom != _chrom:
                header.append((toString(old_chrom), end_list[i - 1]))
            old_chrom = _chrom

        header.append((toString(chrom_list[-1]), end_list[-1]))
        for idx, outfile in enumerate(args.outputFileName):
            log.debug("bigwig: len(vecs_list) {}".format(len(vecs_list)))
            log.debug("bigwig: len(chrom_list) {}".format(len(chrom_list)))

            assert(len(vecs_list) == len(chrom_list))
            _chrom_list = []
            _start_list = []
            _end_list = []
            values = []

            bw = pyBigWig.open(outfile, 'w')
            # set big wig header
            bw.addHeader(header)
            # create entry lists
            for i, value in enumerate(vecs_list):
                # it can happen that some 'value' is having less dimensions than it should
                if len(value) == args.numberOfEigenvectors:
                    if isinstance(value[idx], np.complex):
                        value[idx] = value[idx].real
                    values.append(value[idx])
                    _chrom_list.append(toString(chrom_list[i]))
                    _start_list.append(start_list[i])
                    _end_list.append(end_list[i])

            # write entries
            bw.addEntries(_chrom_list, _start_list, ends=_end_list, values=values)
            bw.close()
    else:
        log.error("Output format not known: {}".format(args.format))
        exit(1)
