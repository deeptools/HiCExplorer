"""hicexplorer.lib.buildMatrixMethods, implemented in C++.

Every public name keeps the name, signature and behaviour it had in
HiCExplorer 3.7.6 and calls hicexplorer._hicx.lib_buildMatrixMethods
(cpp/python/src/api/lib_buildMatrixMethods.cpp): the read pair
classification, the binning and the whole createMatrix run are the C++ of
the build matrix tools (cpp/core/include/hicx/build_matrix_run.hpp). The
functions that take pysam objects keep taking them; the C++ reads the
attributes 3.7.6 read from the objects it is given.
"""
import logging
log = logging.getLogger(__name__)

import argparse
import numpy as np
from scipy.sparse import coo_matrix, dia_matrix
import time
from os import unlink
import os
from pathlib import Path
from io import StringIO
import traceback
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

from collections import OrderedDict

from copy import deepcopy
from ctypes import Structure, c_uint, c_ushort
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Array, RawArray

from intervaltree import IntervalTree, Interval


# own tools
from hicmatrix import HiCMatrix as hm
from hicexplorer.utilities import getUserRegion, genomicRegion
from hicexplorer._version import __version__
import hicexplorer.hicPrepareQCreport as QC

from hicmatrix.lib import MatrixFileHandler

from hicexplorer import hicMergeMatrixBins
import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import lib_buildMatrixMethods as _impl


class C_Interval(_impl.C_Interval):
    """Struct to map a Interval form intervaltree as a multiprocessing.sharedctype"""


class C_Coverage(_impl.C_Coverage):
    """Struct to model the coverage as a multiprocessing.sharedctype"""


class ReadPositionMatrix(object):
    """A class to check for PCR duplicates.
       A set storing all possible
       start sites. Checks if read is already in the set.
    """

    def __init__(self):
        """
        >>> rp = ReadPositionMatrix()
        >>> rp.is_duplicated('1', 0, '2', 0)
        False
        >>> rp.is_duplicated('1', 0, '2', 0)
        True
        """
        _impl.ReadPositionMatrix___init__(self)

    def is_duplicated(self, chrom1, start1, chrom2, start2):
        return _impl.ReadPositionMatrix_is_duplicated(self, chrom1, start1, chrom2, start2)


def intervalListToIntervalTree(interval_list):
    r"""
    given a dictionary containing tuples of chrom, start, end,
    this is transformed to an interval trees. To each
    interval an id is assigned, this id corresponds to the
    position of the interval in the given array of tuples
    and if needed can be used to identify
    the index of a row/colum in the hic matrix.

    >>> bin_list = [('chrX', 0, 50000), ('chrX', 50000, 100000)]
    >>> res = intervalListToIntervalTree(bin_list)
    >>> sorted(res['chrX'])
    [Interval(0, 50000, 0), Interval(50000, 100000, 1)]
    """
    return _impl.intervalListToIntervalTree(interval_list)


def get_bins(bin_size, chrom_size, region=None):
    r"""
    Split the chromosomes into even sized bins
    of length bin_size.

    Returns a list of tuples containing
    ('chrom', start, end)

    >>> from collections import namedtuple
    >>> Handle = namedtuple('Handle', 'references lengths')
    >>> chrom_size = get_chrom_sizes(Handle(('contig-1', 'contig-2'), (7125, 3345)))
    >>> sorted(get_bins(50000, chrom_size))
    [('contig-1', 0, 7125), ('contig-2', 0, 3345)]
    >>> get_bins(50000, chrom_size, region='contig-1')
    [('contig-1', 0, 7125)]
    """
    return _impl.get_bins(bin_size, chrom_size, region)


def bed2interval_list(bed_file_handler, pChromosomeSize, pRegion):
    r"""
    reads a BED file and returns
    a list of tuples containing
    (chromosome name, start, end)

    >>> import tempfile, os

    Make a temporary BED file
    >>> _file = tempfile.NamedTemporaryFile(delete=False)
    >>> _file.close()
    >>> file_tmp = open(_file.name, 'w')
    >>> foo = file_tmp.write('chr1\t10\t20\tH1\t0\n')
    >>> foo = file_tmp.write("chr1\t60\t70\tH2\t0\n")
    >>> file_tmp.close()
    >>> bed2interval_list(open(_file.name, 'r'), None, None)
    [('chr1', 10, 20), ('chr1', 60, 70)]
    >>> os.remove(_file.name)
    """
    return _impl.bed2interval_list(bed_file_handler, pChromosomeSize, pRegion)


def get_rf_bins(rf_cut_intervals, min_distance=200, max_distance=800):
    r"""
    returns a list of tuples containing a bin having the format:
    ('chrom', start, end)

    The goal is to have bins containing each  a restriction site close to
    the center. Restriction sites that are less than 'min_distance' appart
    are merged. To the left and right of the restriction site
    'max_distance' bp are added unless they clash with other bin.
    Resulting bins are not immediately one after the other.

    The given cut sites list has to be ordered chr, and start site

    :param rf_cut_intervals: list of tuples containing the position of restriction fragment sites
    :param min_distance: min distance between restriction fragment sites.
    :param max_distance: max distance between restriction fragment and bin border.
    :return:

    # two restriction sites of length 10
    >>> rf_cut_interval = [('chr1', 10, 20), ('chr1', 60, 70)]

    The following two bins are close together
    and should be  merged under min_distance = 20
    >>> rf_cut_interval.extend([('chr2', 20, 30), ('chr2', 40, 50),
    ... ('chr2', 70, 80)])
    >>> get_rf_bins(rf_cut_interval, min_distance=10, max_distance=20)
    [('chr1', 0, 40), ('chr1', 40, 90), ('chr2', 0, 60), ('chr2', 60, 100)]
    """
    return _impl.get_rf_bins(rf_cut_intervals, min_distance, max_distance)


def get_chrom_sizes(bam_handle):
    """
    return the list of chromosome names and their
    size from the bam file
    The return value is a list of the form
    [('chr1', 2343434), ('chr2', 43432432)]

    The handle needs the attributes `references` and `lengths`, as a
    pysam.AlignmentFile has them.

    >>> from collections import namedtuple
    >>> Handle = namedtuple('Handle', 'references lengths')
    >>> get_chrom_sizes(Handle(('contig-1', 'contig-2'), (7125, 3345)))
    [('contig-1', 7125), ('contig-2', 3345)]
    """
    return _impl.get_chrom_sizes(bam_handle)


def check_dangling_end(read, dangling_sequences):
    """
    given a pysam read object, this function
    checks if a forward read starts with
    the dangling sequence or if a reverse
    read ends with the dangling sequence.
    """
    return _impl.check_dangling_end(read, dangling_sequences)


def get_supplementary_alignment(read, pysam_obj):
    """Checks if a read has a supplementary alignment
    :param read pysam AlignedSegment
    :param pysam_obj pysam file object

    :return pysam AlignedSegment of the supplementary aligment or None in case of no supplementary alignment
    """
    return _impl.get_supplementary_alignment(read, pysam_obj)


def get_correct_map(primary, supplement_list):
    r"""
    Decides which of the mappings, the primary or supplement, is correct. In the case of
    long reads (eg. 150bp) the restriction enzyme site could split the read into two parts
    but only the mapping corresponding to the start of the read should be considered as
    the correct one.

    For example:

    Forward read:

                          _  <- cut site
    |======================|==========>
                           -
    |----------------------|
        correct mapping


    reverse read:

                          _  <- cut site
    <======================|==========|
                           -
                           |----------|
                           correct mapping


    :param primary: pysam AlignedSegment for primary mapping
    :param supplement_list: list of pysam AlignedSegment for secondary mapping

    :return: pysam AlignedSegment that is mapped correctly

    Examples
    --------
    >>> from collections import namedtuple
    >>> Read = namedtuple('Read', 'qname is_reverse cigartuples pos')
    >>> primary = Read('read', False, [(4, 10), (2, 1), (0, 25)], 32)
    >>> secondary_list = [Read('read', False, [(4, 1), (0, 34)], 87),
    ...                   Read('read', False, [(4, 5), (0, 30)], 119)]
    >>> first_mapped = get_correct_map(primary, secondary_list)

    The first mapped read is the first secondary at position 88 (sam 1-based) = 87 (0-based)
    >>> print(first_mapped.pos)
    87

    """
    return _impl.get_correct_map(primary, supplement_list)


def enlarge_bins(bin_intervals, chrom_sizes):
    r"""
    takes a list of consecutive but not
    one after the other bin intervals
    and joins them such that the
    end and start of consecutive bins
    is the same.

    >>> chrom_sizes = [('chr1', 100), ('chr2', 100)]
    >>> bin_intervals =     [('chr1', 10, 30), ('chr1', 50, 80),
    ... ('chr2', 10, 60), ('chr2', 60, 90)]
    >>> enlarge_bins(bin_intervals, chrom_sizes)
    [('chr1', 0, 40), ('chr1', 40, 100), ('chr2', 0, 60), ('chr2', 60, 100)]
    """
    return _impl.enlarge_bins(bin_intervals, chrom_sizes)


def readBamFiles(pFileOneIterator, pFileTwoIterator, pNumberOfItemsPerBuffer, pSkipDuplicationCheck, pReadPosMatrix, pRefId2name, pMinMappingQuality):
    """Read the two bam input files into n buffers each with pNumberOfItemsPerBuffer
        with n = number of processes. The duplication check is handled here too."""
    return _impl.readBamFiles(pFileOneIterator, pFileTwoIterator, pNumberOfItemsPerBuffer, pSkipDuplicationCheck, pReadPosMatrix, pRefId2name, pMinMappingQuality)


def process_data(pMateBuffer1, pMateBuffer2, pMinMappingQuality,
                 pKeepSelfCircles, pRestrictionSequence, pKeepSelfLigation, pMatrixSize,
                 pRfPositions, pRefId2name,
                 pDanglingSequences, pBinsize, pResultIndex,
                 pQueueOut, pTemplate, pOutputBamSet, pCounter,
                 pSharedBinIntvalTree, pDictBinIntervalTreeIndex, pCoverage, pCoverageIndex,
                 pOutputFileBufferDir, pRow, pCol, pData,
                 pMaxInsertSize, pQuickQCMode):
    """
    This function computes for a given number of elements in pMateBuffer1 and pMaterBuffer2 a partial interaction matrix.
    This function is used by multiple processes to speed up the computation.
    All partial matrices are merged in the end into one interaction matrix.

    Parameters
    ----------
    pMateBuffer1 : List of n reads of type 'pysam.libcalignedsegment.AlignedSegment' of sam input file 1
    pMateBuffer2 : List of n reads of type 'pysam.libcalignedsegment.AlignedSegment' of sam input file 2
    pMinMappingQuality : integer, minimum mapping quality of a read
    pKeepSelfCircles : boolean, if self circles should be kept
    pRestrictionSequence : List of String, the restriction sequence
    pKeepSelfLigation : If self ligations should be removed
    pMatrixSize : integer, the size of the interaction matrix
    pRfPositions : intervalTree, only used if a restriction cut file and not a bin size was defined.
    pRefId2name : Tuple, Maps a reference id to a name
    pDanglingSequences : dict, dict of dangling sequences
    pBinsize : integer, the size of the bins
    pResultIndex : integer, number of processs, range(0, threads). Is returned via the queue to have access to the right row, col and data array after the computation.
    pQueueOut : multiprocessing.Queue, queue to return the computed counting variables:
            one_mate_unmapped, one_mate_low_quality, one_mate_not_unique, dangling_end, self_circle, self_ligation, same_fragment,
            mate_not_close_to_rf, count_inward, count_outward, count_left, count_right, inter_chromosomal, short_range, long_range,
            pair_added, len(pMateBuffer1), pResultIndex, pCounter
    pTemplate : The template for the output bam file
    pOutputBamSet : If a output bam file should be written. Depending on the input parameter '--outBam'
    pOutputName : String, Name of the partial bam file
    pCounter : integer, value which is returned to the main process. The main process can than write a pCounter.bam_done file
                to signal the background process, which is merging the partial bam files into one, that this dataset can be merged.
    pSharedBinIntvalTree : multiprocessing.sharedctype.RawArray of C_Interval, stores the interval tree in a 1D-RawArray.
    pDictBinIntervalTreeIndex : dict, stores the information at which index position a given interval starts and ends in the 1D-array 'pSharedBinIntvalTree'
    pCoverage : multiprocessing.sharedctype.Array of c_uint, Stores the coverage in a 1D-Array
    pCoverageIndex :  multiprocessing.sharedctype.RawArray of C_Coverage, stores the information in the 1D-array 'pCoverage'
    pOutputFileBufferDir : String, the directory where the partial output bam files are buffered. Default is '/dev/shm/'
    pRow : multiprocessing.sharedctype.RawArray of c_uint, Stores the row index information. It is available for all processes and does not need to be copied.
    pCol : multiprocessing.sharedctype.RawArray of c_uint, stores the column index information. It is available for all processes and does not need to be copied.
    pData : multiprocessing.sharedctype.RawArray of c_ushort, stores a 1 for each row - column pair. It is available for all processes and does not need to be copied.
    pMaxInsertSize : maximum illumina insert size
    """
    return _impl.process_data(pMateBuffer1, pMateBuffer2, pMinMappingQuality, pKeepSelfCircles, pRestrictionSequence, pKeepSelfLigation, pMatrixSize, pRfPositions, pRefId2name, pDanglingSequences, pBinsize, pResultIndex, pQueueOut, pTemplate, pOutputBamSet, pCounter, pSharedBinIntvalTree, pDictBinIntervalTreeIndex, pCoverage, pCoverageIndex, pOutputFileBufferDir, pRow, pCol, pData, pMaxInsertSize, pQuickQCMode)


def createMatrix(pOutFileName, pMaxDistance, pMaxLibraryInsertSize, pQCfolder,
                 pThreads, pDanglingSequence, pRestrictionSequence, pSamFiles,
                 pDoTestRun, pOutBam, pChromosomeSizes, pRestrictionCutFile,
                 pRegion, pBinSize, pInputBufferSize,
                 pDoTestRunLines, pSkipDuplicationCheck, pMinMappingQuality,
                 pKeepSelfCircles, pKeepSelfLigation, pMinDistance, pGenomeAssembly):
    return _impl.createMatrix(pOutFileName, pMaxDistance, pMaxLibraryInsertSize, pQCfolder, pThreads, pDanglingSequence, pRestrictionSequence, pSamFiles, pDoTestRun, pOutBam, pChromosomeSizes, pRestrictionCutFile, pRegion, pBinSize, pInputBufferSize, pDoTestRunLines, pSkipDuplicationCheck, pMinMappingQuality, pKeepSelfCircles, pKeepSelfLigation, pMinDistance, pGenomeAssembly)


class Tester(object):
    def __init__(self):
        _impl.Tester___init__(self, __file__)
