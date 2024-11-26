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

import pysam
from collections import OrderedDict

from copy import deepcopy
from ctypes import Structure, c_uint, c_ushort
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Array, RawArray

from intervaltree import IntervalTree, Interval

from Bio.Seq import Seq

# own tools
from hicmatrix import HiCMatrix as hm
from hicexplorer.utilities import getUserRegion, genomicRegion
from hicexplorer._version import __version__
import hicexplorer.hicPrepareQCreport as QC

from hicmatrix.lib import MatrixFileHandler

from hicexplorer import hicMergeMatrixBins
import logging
log = logging.getLogger(__name__)


class C_Interval(Structure):
    """Struct to map a Interval form intervaltree as a multiprocessing.sharedctype"""
    _fields_ = [("begin", c_uint),
                ("end", c_uint),
                ("data", c_uint)]


class C_Coverage(Structure):
    """Struct to model the coverage as a multiprocessing.sharedctype"""

    _fields_ = [("begin", c_uint),
                ("end", c_uint)]


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

        self.pos_matrix = set()

    def is_duplicated(self, chrom1, start1, chrom2, start2):
        if chrom1 < chrom2:
            id_string = "{}-{}".format(chrom1, chrom2)
        else:
            id_string = "{}-{}".format(chrom2, chrom1)

        if start1 < start2:
            id_string += "-{}-{}".format(start1, start2)
        else:
            id_string += "-{}-{}".format(start2, start1)

        if id_string in self.pos_matrix:
            return True
        else:
            self.pos_matrix.add(id_string)
            return False


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
    bin_int_tree = {}

    for intval_id, intval in enumerate(interval_list):
        chrom, start, end = intval[0:3]
        if chrom not in bin_int_tree:
            bin_int_tree[chrom] = IntervalTree()
        bin_int_tree[chrom].add(Interval(start, end, intval_id))

    return bin_int_tree


def get_bins(bin_size, chrom_size, region=None):
    r"""
    Split the chromosomes into even sized bins
    of length bin_size.

    Returns a list of tuples containing
    ('chrom', start, end)

    >>> test = Tester()
    >>> chrom_size = get_chrom_sizes(pysam.Samfile(test.bam_file_1))
    >>> sorted(get_bins(50000, chrom_size))
    [('contig-1', 0, 7125), ('contig-2', 0, 3345)]
    >>> get_bins(50000, chrom_size, region='contig-1')
    [('contig-1', 0, 7125)]
    """
    bin_intvals = []
    start = 0
    if region:
        chrom_size, start, _, _ = \
            getUserRegion(chrom_size, region)

    for chrom, size in chrom_size:
        for interval in range(start, size, bin_size):
            bin_intvals.append((chrom, interval,
                                min(size, interval + bin_size)))
    return bin_intvals


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

    if pRegion is not None:
        chrom_size, start_region, end_region, _ = \
            getUserRegion(pChromosomeSize, pRegion)

    count = 0
    interval_list = []
    for line in bed_file_handler:
        count += 1
        fields = line.strip().split()
        try:
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])

        except IndexError:
            log.error("error reading BED file at line {}".format(count))

        if pRegion is not None:
            if chrom == chrom_size[0] and start_region <= start and end_region <= end:
                interval_list.append((chrom, start, end))
        else:
            interval_list.append((chrom, start, end))

    return interval_list


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
    log.info("Minimum distance considered between "
             "restriction sites is {}\nMax "
             "distance: {}\n".format(min_distance, max_distance))

    chrom, start, end = list(zip(*rf_cut_intervals))

    rest_site_len = end[0] - start[0]

    # find sites that are less than min_distance apart
    to_merge = np.flatnonzero(np.diff(start) - rest_site_len <= min_distance)
    to_merge += 1  # + 1 to account for np.diff index handling
    merge_idx = 0
    # add max_distance to both sides
    start = np.array(start) - max_distance
    end = np.array(end) + max_distance
    new_start = [max(0, start[0])]
    new_end = []
    new_chrom = [chrom[0]]
    for idx in range(1, len(start)):
        # identify end of chromosome
        if chrom[idx] != chrom[idx - 1]:
            new_start.append(max(0, start[idx]))
            new_end.append(end[idx - 1])
            new_chrom.append(chrom[idx])
            merge_idx += 1
            continue

        if merge_idx < len(to_merge) and idx == to_merge[merge_idx]:
            merge_idx += 1
            continue

        # identify overlapping bins
        if chrom[idx] == chrom[idx - 1] and end[idx - 1] > start[idx]:
            middle = start[idx] + int((end[idx - 1] - start[idx]) / 2)
            new_start.append(middle)
            new_end.append(middle)
        else:
            new_start.append(start[idx])
            new_end.append(end[idx - 1])
        new_chrom.append(chrom[idx])

    new_end.append(end[-1])
    assert len(new_chrom) == len(new_start), "error"
    assert len(new_end) == len(new_start), "error"

    intervals = zip(new_chrom, new_start, new_end)
    intervals = [(_chrom, _start, _end) for _chrom, _start,
                 _end in intervals if _end - _start >= min_distance]
    return intervals


def get_chrom_sizes(bam_handle):
    """
    return the list of chromosome names and their
    size from the bam file
    The return value is a list of the form
    [('chr1', 2343434), ('chr2', 43432432)]

    >>> test = Tester()
    >>> get_chrom_sizes(pysam.Samfile(test.bam_file_1, 'rb'))
    [('contig-1', 7125), ('contig-2', 3345)]
    """

    # in some cases there are repeated entries in
    # the bam file. Thus, I first convert to dict,
    # then to list.
    list_chrom_sizes = OrderedDict(zip(bam_handle.references,
                                       bam_handle.lengths))
    return list(list_chrom_sizes.items())


def check_dangling_end(read, dangling_sequences):
    """
    given a pysam read object, this function
    checks if a forward read starts with
    the dangling sequence or if a reverse
    read ends with the dangling sequence.
    """
    ds = dangling_sequences
    # check if keys are existing, return false otherwise
    if 'pat_forw' not in ds or 'pat_rev' not in ds:
        return False
    # skip forward read that stars with the restriction sequence
    if not read.is_reverse and \
            read.seq.upper().startswith(ds['pat_forw']):
        return True

    # skip reverse read that ends with the restriction sequence
    if read.is_reverse and \
            read.seq.upper().endswith(ds['pat_rev']):
        return True

    return False


def get_supplementary_alignment(read, pysam_obj):
    """Checks if a read has a supplementary alignment
    :param read pysam AlignedSegment
    :param pysam_obj pysam file object

    :return pysam AlignedSegment of the supplementary aligment or None in case of no supplementary alignment
    """

    # the SA field contains a list of other alignments as a ';' delimited list in the format
    # rname,pos,strand,CIGAR,mapQ,NM;
    if read.has_tag('SA'):
        # field always ends in ';' thus last element after split is always
        # empty, hence [0:-1]
        other_alignments = read.get_tag('SA').split(";")[0:-1]
        supplementary_alignment = []
        for i in range(len(other_alignments)):
            _sup = next(pysam_obj)
            if _sup.qname == read.qname:
                supplementary_alignment.append(_sup)

        return supplementary_alignment

    else:
        return None


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
    >>> sam_file_name = "/tmp/test.sam"
    >>> sam_file = open(sam_file_name, 'w')
    >>> _ = sam_file.write('''@HD\tVN:1.0
    ... @SQ\tSN:chr1\tLN:1575\tAH:chr1:5000000-5010000
    ... read\t65\tchr1\t33\t20\t10S1D25M\t=\t200\t167\tAGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG\t<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<\t
    ... read\t2113\tchr1\t88\t30\t1S34M\t=\t500\t412\tACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA\t<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<\t
    ... read\t2113\tchr1\t120\t30\t5S30M\t=\t500\t412\tACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA\t<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<''')
    >>> sam_file.close()
    >>> sam = pysam.Samfile(sam_file_name, "r")
    >>> read_list = [read for read in sam]
    >>> primary = read_list[0]
    >>> secondary_list = read_list[1:]
    >>> first_mapped = get_correct_map(primary, secondary_list)

    The first mapped read is the first secondary at position 88 (sam 1-based) = 87 (0-based)
    >>> print(first_mapped.pos)
    87

    """

    for supplement in supplement_list:
        assert primary.qname == supplement.qname, "ERROR, primary " \
            "and supplementary reads do not have the same id. The ids " \
            "are as follows\n{}\n{}".format(primary.qname, supplement.qname)
    read_list = [primary] + supplement_list
    first_mapped = []
    for idx, read in enumerate(read_list):
        if read.is_reverse:
            cigartuples = read.cigartuples[::-1]
        else:
            cigartuples = read.cigartuples[:]

        # For each read in read_list, calculate the position of the first match (operation M in CIGAR string) in the read sequence.
        # The calculation is done by adding up the lengths of all the operations until the first match.
        # CIGAR string is a list of tuples of (operation, length). Match is stored as CMATCH.
        match_sum = 0
        for op, count in cigartuples:
            if op == pysam.CMATCH:
                break
            match_sum += count

        first_mapped.append(match_sum)
    # find which read has a cigar string that maps first than any of the
    # others.
    idx_min = first_mapped.index(min(first_mapped))

    return read_list[idx_min]


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
    # enlarge remaining bins
    chr_start = True
    chrom_sizes_dict = dict(chrom_sizes)
    for idx in range(len(bin_intervals) - 1):
        chrom, start, end = bin_intervals[idx]
        chrom_next, start_next, end_next = bin_intervals[idx + 1]
        if chr_start is True:
            start = 0
            chr_start = False
        if chrom == chrom_next and \
                end != start_next:
            middle = start_next - int((start_next - end) / 2)
            bin_intervals[idx] = (chrom, start, middle)
            bin_intervals[idx + 1] = (chrom, middle, end_next)
        if chrom != chrom_next:
            bin_intervals[idx] = (chrom, start, chrom_sizes_dict[chrom])
            bin_intervals[idx + 1] = (chrom_next, 0, end_next)

    chrom, start, end = bin_intervals[-1]
    bin_intervals[-1] = (chrom, start, chrom_sizes_dict[chrom])

    return bin_intervals


def readBamFiles(pFileOneIterator, pFileTwoIterator, pNumberOfItemsPerBuffer, pSkipDuplicationCheck, pReadPosMatrix, pRefId2name, pMinMappingQuality):
    """Read the two bam input files into n buffers each with pNumberOfItemsPerBuffer
        with n = number of processes. The duplication check is handled here too."""
    buffer_mate1 = []
    buffer_mate2 = []
    duplicated_pairs = 0
    one_mate_unmapped = 0
    one_mate_not_unique = 0
    one_mate_low_quality = 0

    all_data_read = False
    j = 0
    iter_num = 0
    while j < pNumberOfItemsPerBuffer:
        try:
            mate1 = next(pFileOneIterator)
            mate2 = next(pFileTwoIterator)
        except StopIteration:
            all_data_read = True
            break
        iter_num += 1

        # skip 'not primary' alignments
        while mate1.flag & 256 == 256:
            try:
                mate1 = next(pFileOneIterator)
            except StopIteration:
                all_data_read = True
                break
        while mate2.flag & 256 == 256:
            try:
                mate2 = next(pFileTwoIterator)
            except StopIteration:
                all_data_read = True
                break

        assert mate1.qname == mate2.qname, "FATAL ERROR {} {} " \
            "Be sure that the sam files have the same read order " \
            "If using Bowtie2 or Hisat2 add " \
            "the --reorder option".format(mate1.qname, mate2.qname)

        # check for supplementary alignments
        # (needs to be done before skipping any unmapped reads
        # to keep the order of the two bam files in sync)
        mate1_supplementary_list = get_supplementary_alignment(
            mate1, pFileOneIterator)
        mate2_supplementary_list = get_supplementary_alignment(
            mate2, pFileTwoIterator)

        if mate1_supplementary_list:
            mate1 = get_correct_map(mate1, mate1_supplementary_list)

        if mate2_supplementary_list:
            mate2 = get_correct_map(mate2, mate2_supplementary_list)

        # skip if any of the reads is not mapped
        if mate1.flag & 0x4 == 4 or mate2.flag & 0x4 == 4:
            one_mate_unmapped += 1
            continue

        # skip if the read quality is low
        if mate1.mapq < pMinMappingQuality or mate2.mapq < pMinMappingQuality:
            # for bwa other way to test
            # for multi-mapping reads is with a mapq = 0
            # the XS flag is not reliable.
            if mate1.mapq == 0 & mate2.mapq == 0:
                one_mate_not_unique += 1
                continue

            """
            # check if low quality is because of
            # read being repetitive
            # by reading the XS flag.
            # The XS:i field is set by bowtie when a read is
            # multi read and it contains the mapping score of the next
            # best match
            if 'XS' in dict(mate1.tags) or 'XS' in dict(mate2.tags):
                one_mate_not_unique += 1
                continue
            """

            one_mate_low_quality += 1
            continue

        if pSkipDuplicationCheck is False:
            if pReadPosMatrix.is_duplicated(pRefId2name[mate1.rname],
                                            mate1.pos,
                                            pRefId2name[mate2.rname],
                                            mate2.pos):
                duplicated_pairs += 1
                continue
        buffer_mate1.append(mate1)
        buffer_mate2.append(mate2)
        j += 1

    if all_data_read and len(buffer_mate1) != 0 and len(buffer_mate2) != 0:
        return buffer_mate1, buffer_mate2, True, duplicated_pairs, one_mate_unmapped, one_mate_not_unique, one_mate_low_quality, iter_num - len(buffer_mate1)
    if all_data_read and len(buffer_mate1) == 0 or len(buffer_mate2) == 0:
        return None, None, True, duplicated_pairs, one_mate_unmapped, one_mate_not_unique, one_mate_low_quality, iter_num - len(buffer_mate1)
    return buffer_mate1, buffer_mate2, False, duplicated_pairs, one_mate_unmapped, one_mate_not_unique, one_mate_low_quality, iter_num - len(buffer_mate1)


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
    try:

        one_mate_unmapped = 0
        one_mate_low_quality = 0
        one_mate_not_unique = 0
        dangling_end = {}
        if pRestrictionSequence is not None:
            for restrictionSequence in pRestrictionSequence:
                dangling_end[restrictionSequence] = 0
        self_circle = 0
        self_ligation = 0
        same_fragment = 0
        mate_not_close_to_rf = 0

        count_inward = 0
        count_outward = 0
        count_left = 0
        count_right = 0
        inter_chromosomal = 0
        short_range = 0
        long_range = 0

        pair_added = 0

        iter_num = 0
        hic_matrix = None

        out_bam_index_buffer = []

        if pMateBuffer1 is None or pMateBuffer2 is None:

            pQueueOut.put([hic_matrix, [one_mate_unmapped, one_mate_low_quality, one_mate_not_unique, dangling_end, self_circle, self_ligation, same_fragment,
                                        mate_not_close_to_rf, count_inward, count_outward,
                                        count_left, count_right, inter_chromosomal, short_range, long_range, pair_added, iter_num, pResultIndex, out_bam_index_buffer]])
            return

        while iter_num < len(pMateBuffer1) and iter_num < len(pMateBuffer2):
            mate1 = pMateBuffer1[iter_num]
            mate2 = pMateBuffer2[iter_num]
            iter_num += 1

            # check if reads belong to a bin
            #
            # pDictBinInterval stores the start and end position for each chromsome in the array 'pSharedBinIntvalTree'
            # To get to the right interval a binary search is used.
            mate_bins = []
            mate_is_unasigned = False
            for mate in [mate1, mate2]:
                mate_ref = pRefId2name[mate.rname]
                # find the middle genomic position of the read. This is used to
                # find the bin it belongs to.
                read_middle = mate.pos + int(mate.qlen / 2)
                try:
                    start, end = pDictBinIntervalTreeIndex[mate_ref]
                    middle_pos = int((start + end) / 2)
                    mate_bin = None
                    while not start > end:
                        if pSharedBinIntvalTree[middle_pos].begin <= read_middle and read_middle <= pSharedBinIntvalTree[middle_pos].end:
                            mate_bin = pSharedBinIntvalTree[middle_pos]
                            mate_is_unasigned = False
                            break
                        elif pSharedBinIntvalTree[middle_pos].begin > read_middle:
                            end = middle_pos - 1
                            middle_pos = int((start + end) / 2)
                            mate_is_unasigned = True
                        else:
                            start = middle_pos + 1
                            middle_pos = int((start + end) / 2)
                            mate_is_unasigned = True

                except KeyError:
                    # for small contigs it can happen that they are not
                    # in the bin_intval_tree keys if no restriction site is found
                    # on the contig.
                    mate_is_unasigned = True
                    break

                # report no match case
                if mate_bin is None:
                    mate_is_unasigned = True
                    break
                mate_bin_id = mate_bin.data
                mate_bins.append(mate_bin_id)

            # if a mate is unassigned, it means it is not close
            # to a restriction site
            if mate_is_unasigned is True:
                mate_not_close_to_rf += 1
                continue

            # check if mates are in the same chromosome
            if mate1.reference_id != mate2.reference_id:
                orientation = 'diff_chromosome'
            else:
                # to identify 'inward' and 'outward' orientations
                # the order or the mates in the genome has to be
                # known.
                if mate1.pos < mate2.pos:
                    first_mate = mate1
                    second_mate = mate2
                else:
                    first_mate = mate2
                    second_mate = mate1

                """
                outward
                <---------------              ---------------->

                inward
                --------------->              <----------------

                same-strand-right
                --------------->              ---------------->

                same-strand-left
                <---------------              <----------------
                """

                if not first_mate.is_reverse and second_mate.is_reverse:
                    orientation = 'inward'
                elif first_mate.is_reverse and not second_mate.is_reverse:
                    orientation = 'outward'
                elif first_mate.is_reverse and second_mate.is_reverse:
                    orientation = 'same-strand-left'
                else:
                    orientation = 'same-strand-right'

                # check self-circles
                # self circles are defined as outward pairs that do not
                # have a restriction sequence in between. The distance of < 25kb is
                # used to only check close outward pairs as far apart pairs can not be self-circles
                if abs(mate2.pos - mate1.pos) < 25000 and orientation == 'outward':
                    has_rf = []

                    if pRfPositions and pRestrictionSequence:
                        # check if in between the two mate
                        # ends the restriction fragment is found.
                        # check for multiple restriction sequences
                        for restrictionSequence in pRestrictionSequence:
                            # the interval used is:
                            # start of fragment + length of restriction sequence
                            # end of fragment - length of restriction sequence
                            # the restriction sequence length is subtracted
                            # such that only fragments internally containing
                            # the restriction site are identified
                            frag_start = min(mate1.pos, mate2.pos) + \
                                len(restrictionSequence)
                            frag_end = max(mate1.pos + mate1.qlen, mate2.pos + mate2.qlen) - len(restrictionSequence)
                            mate_ref = pRefId2name[mate1.rname]
                            if mate_ref in pRfPositions:
                                has_rf.extend(sorted(
                                    pRfPositions[mate_ref][frag_start: frag_end]))

                            if len(has_rf) == 0:
                                self_circle += 1
                                if not pKeepSelfCircles:
                                    continue
                if abs(mate2.pos - mate1.pos) < pMaxInsertSize and orientation == 'inward':
                    # check for dangling ends if the restriction sequence is known and if they look
                    # like 'same fragment'
                    if pRestrictionSequence:
                        if pDanglingSequences:

                            # check for dangling ends in sequence. Stop check with first match
                            one_match = False
                            for restrictionSequence in pRestrictionSequence:
                                if check_dangling_end(mate1, pDanglingSequences[restrictionSequence]) or \
                                        check_dangling_end(mate2, pDanglingSequences[restrictionSequence]):
                                    dangling_end[restrictionSequence] += 1
                                    one_match = True
                                    break
                            if one_match:
                                continue
                    has_rf = []

                    if pRfPositions and pRestrictionSequence:
                        # check if in between the two mate
                        # ends the restriction fragment is found.

                        # the interval used is:
                        # start of fragment + length of restriction sequence
                        # end of fragment - length of restriction sequence
                        # the restriction sequence length is subtracted
                        # such that only fragments internally containing
                        # the restriction site are identified

                        for restrictionSequence in pRestrictionSequence:
                            frag_start = min(mate1.pos, mate2.pos) + \
                                len(restrictionSequence)
                            frag_end = max(mate1.pos + mate1.qlen, mate2.pos + mate2.qlen) - len(restrictionSequence)
                            mate_ref = pRefId2name[mate1.rname]
                            if mate_ref in pRfPositions:
                                has_rf.extend(sorted(
                                    pRfPositions[mate_ref][frag_start: frag_end]))

                    # case when there is no restriction fragment site between the
                    # mates
                    if len(has_rf) == 0:
                        same_fragment += 1
                        continue

                    self_ligation += 1

                    if not pKeepSelfLigation:
                        # skip self ligations
                        continue

                # set insert size to save bam
                mate1.isize = mate2.pos - mate1.pos
                mate2.isize = mate1.pos - mate2.pos

            # if mate_bins, which is set in the previous section
            # does not have size=2, it means that one
            # of the exceptions happened
            if len(mate_bins) != 2:
                continue

            # count type of pair (distance, orientation)
            if mate1.reference_id != mate2.reference_id:
                inter_chromosomal += 1

            elif abs(mate2.pos - mate1.pos) < 20000:
                short_range += 1
            else:
                long_range += 1

            if orientation == 'inward':
                count_inward += 1
            elif orientation == 'outward':
                count_outward += 1
            elif orientation == 'same-strand-left':
                count_left += 1
            elif orientation == 'same-strand-right':
                count_right += 1

            for mate in [mate1, mate2]:
                # fill in coverage vector
                vec_start = int(max(0, mate.pos - mate_bin.begin) / pBinsize)
                length_coverage = pCoverageIndex[mate_bin_id].end - \
                    pCoverageIndex[mate_bin_id].begin
                vec_end = min(length_coverage, int(
                    vec_start + len(mate.seq) / pBinsize))
                coverage_index = pCoverageIndex[mate_bin_id].begin + vec_start
                coverage_end = pCoverageIndex[mate_bin_id].begin + vec_end
                for i in range(coverage_index, coverage_end, 1):
                    pCoverage[i] += 1

            if not pQuickQCMode:
                pRow[pair_added] = mate_bins[0]
                pCol[pair_added] = mate_bins[1]
                pData[pair_added] = np.uint8(1)

            pair_added += 1
            if pOutputBamSet:

                out_bam_index_buffer.append(iter_num - 1)

        pQueueOut.put([[one_mate_unmapped, one_mate_low_quality, one_mate_not_unique, dangling_end, self_circle, self_ligation, same_fragment,
                        mate_not_close_to_rf, count_inward, count_outward,
                        count_left, count_right, inter_chromosomal, short_range, long_range, pair_added, len(pMateBuffer1), pResultIndex, pCounter, out_bam_index_buffer]])
    except Exception as exp:
        pQueueOut.put('Fail: ' + str(exp) + traceback.format_exc())
        return
    return


def createMatrix(pOutFileName, pMaxDistance, pMaxLibraryInsertSize, pQCfolder,
                 pThreads, pDanglingSequence, pRestrictionSequence, pSamFiles,
                 pDoTestRun, pOutBam, pChromosomeSizes, pRestrictionCutFile,
                 pRegion, pBinSize, pInputBufferSize,
                 pDoTestRunLines, pSkipDuplicationCheck, pMinMappingQuality,
                 pKeepSelfCircles, pKeepSelfLigation, pMinDistance, pGenomeAssembly):
    # pOutFileName.name = pOutFileName.name.strip()
    # log.debug('pOutFileName.name: {}'.format(pOutFileName.name.endswith('.h5')))
    if not pOutFileName.name.endswith('.h5') and not pOutFileName.name.endswith('.cool'):
        if '.mcool' not in pOutFileName.name:
            log.error('Please define the file extension. h5 and cool are supported, or the specializations of cool, mcool. Given input {}'.format(pOutFileName.name))
            exit(1)
    # for backwards compatibility
    if pMaxDistance is not None:
        pMaxLibraryInsertSize = pMaxDistance
    try:
        QC.make_sure_path_exists(pQCfolder)
    except OSError:
        exit("Can't open/create QC folder path: {}. Please check".format(pQCfolder))

    if pThreads < 2:
        pThreads = 2
        warnings.warn(
            "\nAt least two threads need to be defined. Setting --threads = 2!s\n")

    if pDanglingSequence and not pRestrictionSequence:
        exit("\nIf --danglingSequence is set, --restrictionSequence needs to be set too.\n")

    log.info("reading {} and {} to build hic_matrix\n".format(pSamFiles[0].name,
                                                              pSamFiles[1].name))
    str1 = pysam.Samfile(pSamFiles[0].name, 'rb')
    str2 = pysam.Samfile(pSamFiles[1].name, 'rb')

    pSamFiles[0].close()
    pSamFiles[1].close()
    if not pDoTestRun:
        if pOutBam:
            pOutBam.close()
            out_bam_file = pysam.Samfile(pOutBam.name, 'wb', template=str1)

    if pChromosomeSizes is None:
        chrom_sizes = get_chrom_sizes(str1)
    else:
        chrom_sizes = OrderedDict()
        with open(pChromosomeSizes.name, 'r') as file:
            file_ = True
            while file_:
                file_ = file.readline().strip()
                if file_ != '':
                    line_split = file_.split('\t')
                    chrom_sizes[line_split[0]] = int(line_split[1])
        chrom_sizes = list(chrom_sizes.items())

    # log.debug('chrom_sizes {}'.format(chrom_sizes))
    read_pos_matrix = ReadPositionMatrix()

    rf_interval = []
    if pRestrictionCutFile:
        for restrictionCutFile in pRestrictionCutFile:
            rf_interval.extend(bed2interval_list(restrictionCutFile, chrom_sizes, pRegion))

        rf_positions = intervalListToIntervalTree(rf_interval)
        log.debug('rf_positions {}'.format(rf_positions.keys()))
    else:
        rf_positions = None

    if pBinSize:
        bin_intervals = get_bins(pBinSize[0], chrom_sizes, pRegion)
    else:
        bin_intervals = get_rf_bins(rf_interval,
                                    min_distance=pMinDistance,
                                    max_distance=pMaxLibraryInsertSize)

    matrix_size = len(bin_intervals)
    bin_intval_tree = intervalListToIntervalTree(bin_intervals)
    ref_id2name = str1.references

    # build c_type shared memory for the interval tree
    shared_array_list = []
    index_dict = {}
    end = -1
    for i, seq in enumerate(bin_intval_tree):
        start = end + 1
        interval_list = []
        for interval in bin_intval_tree[seq]:
            interval_list.append((interval.begin, interval.end, interval.data))
        end = start + len(bin_intval_tree[seq]) - 1
        index_dict[seq] = (start, end)
        interval_list = sorted(interval_list)
        shared_array_list.extend(interval_list)
    shared_build_intval_tree = RawArray(C_Interval, shared_array_list)
    bin_intval_tree = None
    dangling_sequences = {}
    if pDanglingSequence:
        # build a list of dangling sequences
        if pRestrictionSequence is not None:
            for i in range(0, len(pRestrictionSequence)):
                pRestrictionSequence[i] = pRestrictionSequence[i].upper()
                pDanglingSequence[i] = pDanglingSequence[i].upper()
                dangling_sequences[pRestrictionSequence[i]] = {}
                dangling_sequences[pRestrictionSequence[i]]['pat_forw'] = pDanglingSequence[i]
                dangling_sequences[pRestrictionSequence[i]]['pat_rev'] = str(
                    Seq(pDanglingSequence[i]).reverse_complement())

        log.info("dangling sequences to check "
                 "are {}\n".format(dangling_sequences))

    # initialize coverage vectors that
    # save the number of reads that overlap
    # a bin.
    # To save memory, coverage is not measured by bp
    # but by bins of length 10bp
    binsize = 10
    number_of_elements_coverage = 0
    start_pos_coverage = []
    end_pos_coverage = []

    for chrom, start, end in bin_intervals:
        start_pos_coverage.append(number_of_elements_coverage)

        number_of_elements_coverage += (end - start) // binsize
        end_pos_coverage.append(number_of_elements_coverage - 1)
    pos_coverage = RawArray(C_Coverage, list(zip(
        start_pos_coverage, end_pos_coverage)))
    start_pos_coverage = None
    end_pos_coverage = None
    coverage = Array(c_uint, [0] * number_of_elements_coverage)

    # define global shared ctypes arrays for row, col and data
    pThreads = pThreads - 1
    row = [None] * pThreads
    col = [None] * pThreads
    data = [None] * pThreads
    for i in range(pThreads):
        row[i] = RawArray(c_uint, pInputBufferSize)
        col[i] = RawArray(c_uint, pInputBufferSize)
        data[i] = RawArray(c_ushort, pInputBufferSize)

    start_time = time.time()

    iter_num = 0
    # pair_added = 0
    # hic_matrix = None

    one_mate_unmapped = 0
    one_mate_low_quality = 0
    one_mate_not_unique = 0
    dangling_end = {}
    if pRestrictionSequence is not None:
        for restrictionSequence in pRestrictionSequence:
            dangling_end[restrictionSequence] = 0

    self_circle = 0
    self_ligation = 0
    same_fragment = 0
    mate_not_close_to_rf = 0
    duplicated_pairs = 0

    count_inward = 0
    count_outward = 0
    count_left = 0
    count_right = 0
    inter_chromosomal = 0
    short_range = 0
    long_range = 0

    pair_added = 0

    # input buffer for bam files
    buffer_workers1 = [None] * pThreads
    buffer_workers2 = [None] * pThreads

    # output buffer to write bam with mate1 and mate2 pairs
    process = [None] * pThreads
    all_data_processed = False
    hic_matrix = coo_matrix((matrix_size, matrix_size), dtype='uint32')
    queue = [None] * pThreads

    all_threads_done = False
    thread_done = [False] * pThreads
    count_output = 0
    count_call_of_read_input = 0
    computed_pairs = 0

    if pDoTestRun:
        pInputBufferSize = pDoTestRunLines
    fail_flag = False
    fail_message = ''
    while not all_data_processed or not all_threads_done:

        for i in range(pThreads):
            if queue[i] is None and not all_data_processed:
                count_call_of_read_input += 1

                buffer_workers1[i], buffer_workers2[i], all_data_processed, \
                    duplicated_pairs_, one_mate_unmapped_, one_mate_not_unique_, \
                    one_mate_low_quality_, iter_num_ = readBamFiles(pFileOneIterator=str1,
                                                                    pFileTwoIterator=str2,
                                                                    pNumberOfItemsPerBuffer=pInputBufferSize,
                                                                    pSkipDuplicationCheck=pSkipDuplicationCheck,
                                                                    pReadPosMatrix=read_pos_matrix,
                                                                    pRefId2name=ref_id2name,
                                                                    pMinMappingQuality=pMinMappingQuality
                                                                    )
                duplicated_pairs += duplicated_pairs_
                one_mate_unmapped += one_mate_unmapped_
                one_mate_not_unique += one_mate_not_unique_
                one_mate_low_quality += one_mate_low_quality_
                iter_num += iter_num_
                queue[i] = Queue()
                thread_done[i] = False
                computed_pairs += len(buffer_workers1[i])
                # create process to compute hic matrix for this buffer
                process[i] = Process(target=process_data, kwargs=dict(
                    pMateBuffer1=buffer_workers1[i],
                    pMateBuffer2=buffer_workers2[i],
                    pMinMappingQuality=pMinMappingQuality,
                    pKeepSelfCircles=pKeepSelfCircles,
                    pRestrictionSequence=pRestrictionSequence,
                    pKeepSelfLigation=pKeepSelfLigation,
                    pMatrixSize=matrix_size,
                    pRfPositions=rf_positions,
                    pRefId2name=ref_id2name,
                    pDanglingSequences=dangling_sequences,
                    pBinsize=binsize,
                    pResultIndex=i,
                    pQueueOut=queue[i],
                    pTemplate=str1,
                    pOutputBamSet=pOutBam,
                    pCounter=count_output,
                    pSharedBinIntvalTree=shared_build_intval_tree,
                    pDictBinIntervalTreeIndex=index_dict,
                    pCoverage=coverage,
                    pCoverageIndex=pos_coverage,
                    pOutputFileBufferDir="",
                    pRow=row[i],
                    pCol=col[i],
                    pData=data[i],
                    pMaxInsertSize=pMaxLibraryInsertSize,
                    pQuickQCMode=pDoTestRun
                ))
                process[i].start()
                count_output += 1

            elif queue[i] is not None and not queue[i].empty():
                result = queue[i].get()
                if 'Fail:' in result:
                    fail_flag = True
                    fail_message = result[6:]
                else:
                    if result[0] is not None:
                        elements = result[0][15]
                        if hic_matrix is None:
                            hic_matrix = coo_matrix(
                                (data[i][:elements], (row[i][:elements], col[i][:elements])), shape=(matrix_size, matrix_size))
                        else:
                            hic_matrix += coo_matrix(
                                (data[i][:elements], (row[i][:elements], col[i][:elements])), shape=(matrix_size, matrix_size))

                        for sequence in result[0][3]:
                            dangling_end[sequence] += result[0][3][sequence]
                        self_circle += result[0][4]
                        self_ligation += result[0][5]
                        same_fragment += result[0][6]
                        mate_not_close_to_rf += result[0][7]

                        count_inward += result[0][8]
                        count_outward += result[0][9]
                        count_left += result[0][10]
                        count_right += result[0][11]
                        inter_chromosomal += result[0][12]
                        short_range += result[0][13]
                        long_range += result[0][14]

                        pair_added += result[0][15]
                        iter_num += result[0][16]

                    for bam_index in result[0][19]:
                        mate1 = buffer_workers1[i][bam_index]
                        mate2 = buffer_workers2[i][bam_index]

                        mate1.flag |= 0x1
                        mate2.flag |= 0x1

                        # set one read as the first in pair and the
                        # other as second
                        mate1.flag |= 0x40
                        mate2.flag |= 0x80

                        # set chrom of mate
                        mate1.mrnm = mate2.rname
                        mate2.mrnm = mate1.rname

                        # set position of mate
                        mate1.mpos = mate2.pos
                        mate2.mpos = mate1.pos

                        out_bam_file.write(mate1)
                        out_bam_file.write(mate2)

                buffer_workers1[i] = None
                buffer_workers2[i] = None
                queue[i] = None
                process[i].join()
                process[i].terminate()
                process[i] = None
                thread_done[i] = True

                # caused by the architecture I try to display this output
                # information after +-1e5 of 1e6 reads.
                if iter_num % 1e6 < 100000:
                    elapsed_time = time.time() - start_time
                    log.info("processing {} lines took {:.2f} "
                             "secs ({:.1f} lines per "
                             "second)\n".format(iter_num,
                                                elapsed_time,
                                                iter_num / elapsed_time))
                    log.info("{} ({:.2f}%) valid pairs added to matrix"
                             "\n".format(pair_added, float(100 * pair_added) / iter_num))
                if pDoTestRun and iter_num > pDoTestRunLines:
                    log.debug(
                        "\n## *WARNING*. Early exit because of --doTestRun parameter  ##\n\n")
                    all_data_processed = True
                    thread_done[i] = True
                    break
            elif all_data_processed and queue[i] is None:
                thread_done[i] = True
            else:
                time.sleep(1)

        if all_data_processed:
            all_threads_done = True
            for thread in thread_done:
                if not thread:
                    all_threads_done = False
    if fail_flag:
        log.error(fail_message)
        exit(1)
    else:
        log.debug('Parallel stuff done')
    if not pDoTestRun:
        # the resulting matrix is only filled unevenly with some pairs
        # int the upper triangle and others in the lower triangle. To construct
        # the definite matrix I add the values from the upper and lower triangles
        # and subtract the diagonal to avoid double counting it.
        # The resulting matrix is symmetric.
        if pOutBam:
            out_bam_file.close()

        dia = dia_matrix(([hic_matrix.diagonal()], [0]),
                         shape=hic_matrix.shape)
        hic_matrix = hic_matrix + hic_matrix.T - dia
        # extend bins such that they are next to each other
        bin_intervals = enlarge_bins(bin_intervals[:], chrom_sizes)
        # compute max bin coverage
        bin_max = []

        for cover in pos_coverage:
            max_element = 0
            for i in range(cover.begin, cover.end, 1):
                if coverage[i] > max_element:
                    max_element = coverage[i]
            if max_element == 0:
                bin_max.append(np.nan)
            else:
                bin_max.append(max_element)

        chr_name_list, start_list, end_list = list(zip(*bin_intervals))
        bin_intervals = list(zip(chr_name_list, start_list, end_list, bin_max))
        hic_ma = hm.hiCMatrix()
        hic_ma.setMatrix(hic_matrix, cut_intervals=bin_intervals)

    """
    if pRestrictionCutFile:
        # load the matrix to mask those
        # bins that most likely didn't
        # have a restriction site that was cutted

        # reload the matrix as a Hi-C matrix object
        hic_matrix = hm.hiCMatrix(pOutFileName.name)

        hic_matrix.maskBins(get_poor_bins(bin_max))
        hic_matrix.save(pOutFileName.name)
    """
    if not pKeepSelfLigation:
        msg = " (removed)"
    else:
        msg = " (not removed)"

    mappable_unique_high_quality_pairs = iter_num - \
        (one_mate_unmapped + one_mate_low_quality + one_mate_not_unique)

    intermediate_qc_log = StringIO()

    if pMinDistance:
        intermediate_qc_log.write("""
File\t{}\t\t
Sequenced reads\t{}\t\t
Min rest. site distance\t{}\t\t
Max library insert size\t{}\t\t

""".format(pOutFileName.name, iter_num, pMinDistance, pMaxLibraryInsertSize))
    else:
        intermediate_qc_log.write("""
File\t{}\t\t
Sequenced reads\t{}\t\t
Max library insert size\t{}\t\t

""".format(pOutFileName.name, iter_num, pMaxLibraryInsertSize))

    intermediate_qc_log.write(
        "#\tcount\t(percentage w.r.t. total sequenced reads)\n")

    intermediate_qc_log.write("Pairs mappable, unique and high quality\t{}\t({:.2f})\n".
                              format(mappable_unique_high_quality_pairs,
                                     100 * float(mappable_unique_high_quality_pairs) / iter_num))

    intermediate_qc_log.write("Hi-C contacts\t{}\t({:.2f})\n".
                              format(pair_added, 100 * float(pair_added) / iter_num))

    intermediate_qc_log.write("One mate unmapped\t{}\t({:.2f})\n".
                              format(one_mate_unmapped, 100 * float(one_mate_unmapped) / iter_num))

    intermediate_qc_log.write("One mate not unique\t{}\t({:.2f})\n".
                              format(one_mate_not_unique, 100 * float(one_mate_not_unique) / iter_num))

    intermediate_qc_log.write("Low mapping quality\t{}\t({:.2f})\n".
                              format(one_mate_low_quality, 100 * float(one_mate_low_quality) / iter_num))

    intermediate_qc_log.write(
        "\n#\tcount\t(percentage w.r.t. mappable, unique and high quality pairs)\n")

    if len(dangling_end) > 0:
        log.debug('dangling_sequences {}'.format(dangling_sequences.keys()))
        log.debug('dangling_end {}'.format(dangling_end.keys()))
        # log.debug('dangling_end {}'.format(res.keys()))

        for key in dangling_end:
            # dangling_sequences[pRestrictionSequence[i]]['pat_forw']
            intermediate_qc_log.write("dangling end {} (restriction sequence {})\t{}\t({:.2f})\n".
                                      format(dangling_sequences[key]['pat_forw'], key, dangling_end[key], 100 * float(dangling_end[key]) / mappable_unique_high_quality_pairs))
    # else:
    #     intermediate_qc_log.write("dangling end\t{}\t({:.2f})\n".
    #                               format(0, 100 * float(0) / mappable_unique_high_quality_pairs))
    if pRestrictionCutFile is not None:
        intermediate_qc_log.write("self ligation{}\t{}\t({:.2f})\n".
                                  format(msg, self_ligation, 100 * float(self_ligation) / mappable_unique_high_quality_pairs))

        intermediate_qc_log.write("One mate not close to rest site\t{}\t({:.2f})\n".
                                  format(mate_not_close_to_rf, 100 * float(mate_not_close_to_rf) / mappable_unique_high_quality_pairs))
    intermediate_qc_log.write("same fragment\t{}\t({:.2f})\n".
                              format(same_fragment, 100 * float(same_fragment) / mappable_unique_high_quality_pairs))
    if pRestrictionCutFile is not None:
        intermediate_qc_log.write("self circle\t{}\t({:.2f})\n".
                                  format(self_circle, 100 * float(self_circle) / mappable_unique_high_quality_pairs))

    intermediate_qc_log.write("duplicated pairs\t{}\t({:.2f})\n".
                              format(duplicated_pairs, 100 * float(duplicated_pairs) / mappable_unique_high_quality_pairs))

    if pair_added > 0:
        intermediate_qc_log.write(
            "\n#\tcount\t(percentage w.r.t. total valid pairs used)\n")
        intermediate_qc_log.write("inter chromosomal\t{}\t({:.2f})\n".
                                  format(inter_chromosomal, 100 * float(inter_chromosomal) / pair_added))

        intermediate_qc_log.write("Intra short range (< 20kb)\t{}\t({:.2f})\n".
                                  format(short_range, 100 * float(short_range) / pair_added))

        intermediate_qc_log.write("Intra long range (>= 20kb)\t{}\t({:.2f})\n".
                                  format(long_range, 100 * float(long_range) / pair_added))

        intermediate_qc_log.write("Read pair type: inward pairs\t{}\t({:.2f})\n".
                                  format(count_inward, 100 * float(count_inward) / pair_added))

        intermediate_qc_log.write("Read pair type: outward pairs\t{}\t({:.2f})\n".
                                  format(count_outward, 100 * float(count_outward) / pair_added))

        intermediate_qc_log.write("Read pair type: left pairs\t{}\t({:.2f})\n".
                                  format(count_left, 100 * float(count_left) / pair_added))

        intermediate_qc_log.write("Read pair type: right pairs\t{}\t({:.2f})\n".
                                  format(count_right, 100 * float(count_right) / pair_added))

    log_file_name = os.path.join(pQCfolder, "QC.log")
    log_file = open(log_file_name, 'w')
    log_file.write(intermediate_qc_log.getvalue())
    log_file.close()
    log.debug('log_file_name {}'.format(log_file_name))
    log.debug('pQCfolder {}'.format(pQCfolder))

    QC.main("-l {} -o {}".format(log_file_name, pQCfolder).split())

    pOutFileName.close()
    # removing the empty file. Otherwise the save method
    # will say that the file already exists.
    unlink(pOutFileName.name)

    hic_metadata = {}
    hic_metadata['statistics'] = intermediate_qc_log.getvalue()
    hic_metadata['matrix-generated-by'] = np.string_(
        'HiCExplorer-' + __version__)
    hic_metadata['matrix-generated-by-url'] = np.string_(
        'https://github.com/deeptools/HiCExplorer')
    if pGenomeAssembly:
        hic_metadata['genome-assembly'] = np.string_(pGenomeAssembly)

    intermediate_qc_log.close()
    if pOutFileName.name.endswith('.mcool') and pBinSize is not None and len(pBinSize) > 2:

        matrixFileHandlerOutput = MatrixFileHandler(
            pFileType='cool', pHiCInfo=hic_metadata)
        matrixFileHandlerOutput.set_matrix_variables(hic_ma.matrix,
                                                     hic_ma.cut_intervals,
                                                     hic_ma.nan_bins,
                                                     hic_ma.correction_factors,
                                                     hic_ma.distance_counts)
        matrixFileHandlerOutput.save(pOutFileName.name + '::/resolutions/' + str(
            pBinSize[0]), pSymmetric=True, pApplyCorrection=False)

        for resolution in pBinSize[1:]:
            hic_matrix_to_merge = deepcopy(hic_ma)
            _mergeFactor = int(resolution) // pBinSize[0]
            merged_matrix = hicMergeMatrixBins.merge_bins(
                hic_matrix_to_merge, _mergeFactor)
            matrixFileHandlerOutput = MatrixFileHandler(
                pFileType='cool', pAppend=True, pHiCInfo=hic_metadata)
            matrixFileHandlerOutput.set_matrix_variables(merged_matrix.matrix,
                                                         merged_matrix.cut_intervals,
                                                         merged_matrix.nan_bins,
                                                         merged_matrix.correction_factors,
                                                         merged_matrix.distance_counts)
            matrixFileHandlerOutput.save(pOutFileName.name + '::/resolutions/' + str(
                resolution), pSymmetric=True, pApplyCorrection=False)

    else:
        if not pDoTestRun:
            hic_ma.save(pOutFileName.name, pHiCInfo=hic_metadata)


class Tester(object):
    def __init__(self):
        hic_test_data_dir = os.environ.get('HIC_TEST_DATA_DIR', False)
        if hic_test_data_dir:
            self.root = hic_test_data_dir
        else:
            self.root = os.path.dirname(os.path.dirname(
                os.path.abspath(__file__))) + "/test/test_data/"
        self.bam_file_1 = os.path.join(self.root, "hic.bam")
