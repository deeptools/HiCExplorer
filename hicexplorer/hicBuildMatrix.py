import argparse
import sys
import numpy as np
from scipy.sparse import coo_matrix, dia_matrix
import time
from os import unlink
import os
import shutil
import pysam

from six.moves import xrange

from ctypes import Structure, c_uint, c_ushort
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Array, RawArray

from intervaltree import IntervalTree, Interval

# own tools
from hicexplorer import HiCMatrix as hm
from hicexplorer.utilities import getUserRegion, genomicRegion
from hicexplorer._version import __version__
import hicexplorer.hicPrepareQCreport as QC


class C_Interval(Structure):
    _fields_ = [("begin", c_uint),
                ("end", c_uint),
                ("data", c_uint)]


class C_Coverage(Structure):
    _fields_ = [("begin", c_uint),
                ("end", c_uint)]


class SharedData():
    """A class to wrap the input to process the buffers with the function \'process data\'."""

    def __init__(self, pMinMappingQuality, pRemoveSelfCircles,
                 pRestrictionSequence, pRemoveSelfLigation, pReadPosMatrix, pRfPositions,
                 pRefId2name, pDanglingSequences, pBinsize, pTemplate,
                 pSharedBinIntvalTree, pDictBinIntervalTreeIndex, pCoverage,
                 pCoverageIndex, pOutputFileBufferDir):

        self.min_mapping_quality = pMinMappingQuality
        self.removeSelfCircles = pRemoveSelfCircles
        self.restrictionSequence = pRestrictionSequence
        self.removeSelfLigation = pRemoveSelfLigation
        self.readPosMatrix = pReadPosMatrix
        self.rfPositions = pRfPositions
        self.refId2name = pRefId2name
        self.danglingSequences = pDanglingSequences
        self.binsize = pBinsize
        self.template = pTemplate
        # self.outputName = pOutputName
        # self.counter = pCounter
        self.sharedBinIntvalTree = pSharedBinIntvalTree
        self.dictBinIntervalTreeIndex = pDictBinIntervalTreeIndex
        self.coverage = pCoverage
        self.coverageIndex = pCoverageIndex
        self.outputFileBufferDir = pOutputFileBufferDir


class ProcessData():
    """A class to wrap all variables which are unique for a process. """

    def __init__(self, pBuffer1=None, pBuffer2=None, pRow=None,
                 pCol=None, pData=None, pResultIndex=None,
                 pOutputName=None, pCounter=None, pQueueOut=None):

        self.one_mate_unmapped = 0
        self.one_mate_low_quality = 0
        self.one_mate_not_unique = 0
        self.dangling_end = 0
        self.self_circle = 0
        self.self_ligation = 0
        self.same_fragment = 0
        self.mate_not_close_to_rf = 0
        self.count_inward = 0
        self.count_outward = 0
        self.count_left = 0
        self.count_right = 0
        self.inter_chromosomal = 0
        self.short_range = 0
        self.long_range = 0
        self.pair_added = 0
        self.iter_num = 0
        self.duplicated_pairs = 0
        self.one_mate_unmapped = 0
        self.one_mate_not_unique = 0
        self.one_mate_low_quality = 0

        self.result_index = pResultIndex
        self.output_name = pOutputName
        self.counter = pCounter

        self.row = pRow
        self.col = pCol
        self.data = pData
        self.buffer1 = pBuffer1
        self.buffer2 = pBuffer2
        self.queue_out = pQueueOut

    def synchronize(self, pProcessDataObj):
        self.one_mate_unmapped += pProcessDataObj.one_mate_unmapped
        self.one_mate_low_quality += pProcessDataObj.one_mate_low_quality
        self.one_mate_not_unique += pProcessDataObj.one_mate_not_unique
        self.dangling_end += pProcessDataObj.dangling_end
        self.self_circle += pProcessDataObj.self_circle
        self.self_ligation += pProcessDataObj.self_ligation
        self.same_fragment += pProcessDataObj.same_fragment
        self.mate_not_close_to_rf += pProcessDataObj.mate_not_close_to_rf
        self.count_inward += pProcessDataObj.count_inward
        self.count_outward += pProcessDataObj.count_outward
        self.count_left += pProcessDataObj.count_left
        self.count_right += pProcessDataObj.count_right
        self.inter_chromosomal += pProcessDataObj.inter_chromosomal
        self.short_range += pProcessDataObj.short_range
        self.long_range += pProcessDataObj.long_range
        self.pair_added += pProcessDataObj.pair_added
        self.iter_num += pProcessDataObj.iter_num

    def processing_done(self):
        if self.queue_out is not None and not self.queue_out.empty():
            return True
        return False


class ReadPositionMatrix(object):
    """ class to check for PCR duplicates.
    A sparse matrix having as bins all possible
    start sites (single bp resolution)
    is created. PCR duplicates
    are determined by checking if the matrix
    cell is already filled.

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
        """Checks if some sequence was seen before and is a duplicate.
        Returns true if it is a duplicate, returns false if not and inserts the elements into  a set.
        """
        if chrom1 < chrom2:
            id_string = "%s-%s" % (chrom1, chrom2)
        else:
            id_string = "%s-%s" % (chrom2, chrom1)

        if start1 < start2:
            id_string += "-%s-%s" % (start1, start2)
        else:
            id_string += "-%s-%s" % (start2, start1)

        if id_string in self.pos_matrix:
            return True
        else:
            self.pos_matrix.add(id_string)
            return False


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=('Using an alignment from Bowtie2 where both '
                     'PE reads are mapped using  the --local '
                     'option, this program reads such file and '
                     'creates a matrix of interactions.'))

    # define the arguments
    parser.add_argument('--samFiles', '-s',
                        help='The two sam files to process',
                        metavar='two sam files',
                        nargs=2,
                        type=argparse.FileType('r'),
                        required=True)

    # define the arguments
    parser.add_argument('--outBam', '-b',
                        help='Bam file to process. Optional parameter. Computation will be significant longer if this option is set.',
                        metavar='bam file',
                        type=argparse.FileType('w'),
                        required=True)

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('--binSize', '-bs',
                       help='Size in bp for the bins.',
                       type=int,
                       default=10000)

    group.add_argument('--restrictionCutFile', '-rs',
                       help=('BED file with all restriction cut places. '
                             'Should contain only  mappable '
                             'restriction sites. If given, the bins are '
                             'set to match the restriction fragments (i.e. '
                             'the region between one restriction site and '
                             'the next).'),
                       type=argparse.FileType('r'),
                       metavar='BED file')

    parser.add_argument('--minDistance',
                        help='Minimum distance between restriction sites. '
                        'Restriction sites that are closer that this '
                        'distance are merged into one. This option only '
                        'applies if --restrictionCutFile is given.',
                        type=int,
                        default=300,
                        required=False)

    parser.add_argument('--maxDistance',
                        help='Maximum distance in bp from restriction site '
                        'to read, to consider a read a valid one. This option '
                        'only applies if --restrictionCutFile is given.',
                        type=int,
                        default=800,
                        required=False)

    parser.add_argument('--restrictionSequence', '-seq',
                        help='Sequence of the restriction site. This is used '
                        'to discard reads that end/start with such sequence '
                        'and that are considered un-ligated fragments or '
                        '"dangling-ends". If not given, such statistics will '
                        'not be available.')

    parser.add_argument('--outFileName', '-o',
                        help='Output file name for a matrix',
                        metavar='FILENAME',
                        type=argparse.FileType('w'),
                        required=True)

    parser.add_argument('--QCfolder',
                        help='Path of folder to save the quality control data for the matrix',
                        metavar='FOLDER',
                        required=True)

    parser.add_argument('--region', '-r',
                        help='Region of the genome to limit the operation. '
                        'The format is chr:start-end. Also valid is just to '
                        'specify a chromosome, for example --region chr10',
                        metavar="CHR:START-END",
                        required=False,
                        type=genomicRegion
                        )

    parser.add_argument('--removeSelfLigation',
                        # help='If set, inward facing reads less than 1000 bp apart and having a restriction'
                        #     'site in between are removed. Although this reads do not contribute to '
                        #     'any distant contact, they are useful to account for bias in the data.',
                        help=argparse.SUPPRESS,
                        required=False,
                        default=True
                        # action='store_true'
                        )

    parser.add_argument('--removeSelfCircles',
                        help='If set, outward facing reads, at a distance of less thatn 25kbs are removed.',
                        required=False,
                        action='store_true'
                        )

    parser.add_argument('--minMappingQuality',
                        help='minimun mapping quality for reads to be accepted. Because the restriction '
                             'enzyme site could be located on top of the read, this may reduce the '
                             'reported quality of the read. Thus, this parameter may be adusted if too many '
                             'low quality (but otherwise perfectly valid hic-reads) are found. A good strategy '
                             'is to make a test run (using the --doTestRun), then checking the results to see '
                             'if too many low quality reads are present and then using the bam file generated to '
                             'check if those low quality reads are caused by the read not being mapped entirely.',
                        required=False,
                        default=15,
                        type=int
                        )
    parser.add_argument('--threads',
                        help='Number of threads. Using the python multiprocessing module.'
                        ' One master process which is used to read the input file into the buffer and one process which is merging '
                        'the output bam files of the processes into one output bam file.'
                        ' This means that two processes less as defined do the computation. Minimum value is 3.',
                        required=False,
                        default=4,
                        type=int
                        )
    parser.add_argument('--inputBufferSize',
                        help='Size of the input buffer of each thread. 100,000 elements per input file per thread is the default value.'
                             ' Reduce value to decrease memory usage.',
                        required=False,
                        default=1e5,
                        type=int
                        )
    parser.add_argument('--outputFileBufferDir',
                        help='The location of the output file buffer. Per default /dev/shm/ is used which is in the most Linux systems a RAM disk. '
                        'Please make sure no other instance of hicBuildMatrix is accessing this directory at the same time or that old tmp files, maybe from '
                        'an interupted run of hicBuildMatrix, are stored there. It could cause some non expected behaviour and or results.',
                        required=False,
                        default='/dev/shm/',
                        type=str
                        )
    parser.add_argument('--doTestRun',
                        help='A test run is useful to test the quality of a Hi-C experiment quickly. It works by '
                             'testing only 1,000.000 reads. This option is useful to get an idea of quality control'
                             'values like inter-chromosomal interactins, duplication rates etc.',
                        action='store_true'
                        )

    parser.add_argument('--skipDuplicationCheck',
                        help='Identification of duplicated read pairs is memory consuming. Thus, in case of '
                             'memory errors this check can be skipped. However, consider running a `--doTestRun` '
                             'first to get an estimation of the duplicated reads. ',
                        action='store_true'
                        )

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def intervalListToIntervalTree(interval_list):
    """
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
    """
    Split the chromosomes into even sized bins
    of length bin_size.

    Returns a list of tuples containing
    ('chrom', start, end)

    >>> test = Tester()
    >>> chrom_size = get_chrom_sizes(pysam.Samfile(test.bam_file_1))
    >>> get_bins(50000, chrom_size)
    [('contig-2', 0, 3345), ('contig-1', 0, 7125)]
    >>> get_bins(50000, chrom_size, region='contig-1')
    [('contig-1', 0, 7125)]
    """
    bin_intvals = []
    start = 0
    if region:
        chrom_size, start, _, _ = \
            getUserRegion(chrom_size, region)

    for chrom, size in chrom_size:
        for interval in xrange(start, size, bin_size):
            bin_intvals.append((chrom, interval,
                                min(size, interval + bin_size)))
    return bin_intvals


def bed2interval_list(bed_file_handler):
    r"""
    reads a BED file and returns
    a list of tuples containing
    (chromosome name, start, end)

    >>> import tempfile, os

    Make a temporary BED file
    >>> _file = tempfile.NamedTemporaryFile(delete=False)
    >>> _file.write('chr1\t10\t20\tH1\t0\n')
    >>> _file.write("chr1\t60\t70\tH2\t0\n")
    >>> _file.close()
    >>> bed2interval_list(open(_file.name))
    [('chr1', 10, 20), ('chr1', 60, 70)]
    >>> os.remove(_file.name)
    """
    count = 0
    interval_list = []
    for line in bed_file_handler:
        count += 1
        fields = line.strip().split()
        try:
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        except IndexError:
            sys.stderr.write("error reading BED file at line {}".format(count))

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
    sys.stderr.write("Minimum distance considered between "
                     "restriction sites is {}\nMax "
                     "distance: {}\n".format(min_distance, max_distance))

    chrom, start, end = zip(*rf_cut_intervals)
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
    intervals = [(_chrom, _start, _end) for _chrom, _start, _end in intervals if _end - _start >= min_distance]
    return intervals


def get_chrom_sizes(bam_handle):
    """
    return the list of chromosome names and their
    size from the bam file
    The return value is a list of the form
    [('chr1', 2343434), ('chr2', 43432432)]

    >>> test = Tester()
    >>> get_chrom_sizes(pysam.Samfile(test.bam_file_1, 'rb'))
    [('contig-2', 3345), ('contig-1', 7125)]
    """

    # in some cases there are repeated entries in
    # the bam file. Thus, I first convert to dict,
    # then to list.
    list_chrom_sizes = dict(zip(bam_handle.references,
                                bam_handle.lengths))
    return list_chrom_sizes.items()


def check_dangling_end(read, dangling_sequences):
    """
    given a pysam read object, this function
    checks if a forward read starts with
    the dangling sequence or if a reverse
    read ends with the dangling sequence.
    """
    ds = dangling_sequences
    # skip forward read that stars with the restriction sequence
    if not read.is_reverse and \
            read.seq.upper()[0:len(ds['pat_forw'])] == ds['pat_forw']:
        return True

    # skip reverse read that ends with the restriction sequence
    if read.is_reverse and \
            read.seq.upper()[-len(ds['pat_rev']):] == ds['pat_rev']:
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
        # field always ends in ';' thus last element after split is always empty, hence [0:-1]
        other_alignments = read.get_tag('SA').split(";")[0:-1]
        supplementary_alignment = []
        for i in range(len(other_alignments)):
            _sup = pysam_obj.next()
            if _sup.is_supplementary and _sup.qname == read.qname:
                supplementary_alignment.append(_sup)

        return supplementary_alignment

    else:
        return None


def get_correct_map(primary, supplement_list):
    """
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

        first_mapped.append([x for x, cig in enumerate(cigartuples) if cig[0] == 0][0])
    # find which read has a cigar string that maps first than any of the others.
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
    for idx in xrange(len(bin_intervals) - 1):
        chrom, start, end = bin_intervals[idx]
        chrom_next, start_next, end_next = bin_intervals[idx + 1]
        if chr_start is True:
            start = 0
            chr_start = False
        if chrom == chrom_next and \
                end != start_next:
            middle = start_next - (start_next - end) / 2
            bin_intervals[idx] = (chrom, start, middle)
            bin_intervals[idx + 1] = (chrom, middle, end_next)
        if chrom != chrom_next:
            bin_intervals[idx] = (chrom, start, chrom_sizes_dict[chrom])
            bin_intervals[idx + 1] = (chrom_next, 0, end_next)

    chrom, start, end = bin_intervals[-1]
    bin_intervals[-1] = (chrom, start, chrom_sizes_dict[chrom])

    return bin_intervals


def readBamFiles(pFileOneIterator, pFileTwoIterator, pNumberOfItemsPerBuffer,
                 pSkipDuplicationCheck, pReadPosMatrix, pRefId2name, pMinMappingQuality,
                 pProcessData):
    """Read the two bam input files into n buffers each with pNumberOfItemsPerBuffer
        with n = number of processes. The duplication check is handled here too."""
    buffer_mate1 = []
    buffer_mate2 = []
    # duplicated_pairs = 0
    # one_mate_unmapped = 0
    # one_mate_not_unique = 0
    # one_mate_low_quality = 0

    all_data_read = False
    j = 0
    iter_num = 0
    while j < pNumberOfItemsPerBuffer:
        try:
            mate1 = pFileOneIterator.next()
            mate2 = pFileTwoIterator.next()
        except StopIteration:
            all_data_read = True
            break
        iter_num += 1

        # skip 'not primary' alignments
        while mate1.flag & 256 == 256:
            try:
                mate1 = pFileOneIterator.next()
            except StopIteration:
                all_data_read = True
                break
        while mate2.flag & 256 == 256:
            try:
                mate2 = pFileTwoIterator.next()
            except StopIteration:
                # pTerminateSignal.set()
                all_data_read = True
                break

        assert mate1.qname == mate2.qname, "FATAL ERROR {} {} " \
            "Be sure that the sam files have the same read order " \
            "If using Bowtie2 or Hisat2 add " \
            "the --reorder option".format(mate1.qname, mate2.qname)

        # check for supplementary alignments
        # (needs to be done before skipping any unmapped reads
        # to keep the order of the two bam files in sync)
        mate1_supplementary_list = get_supplementary_alignment(mate1, pFileOneIterator)
        mate2_supplementary_list = get_supplementary_alignment(mate2, pFileTwoIterator)

        if mate1_supplementary_list:
            mate1 = get_correct_map(mate1, mate1_supplementary_list)

        if mate2_supplementary_list:
            mate2 = get_correct_map(mate2, mate2_supplementary_list)

        # skip if any of the reads is not mapped
        if mate1.flag & 0x4 == 4 or mate2.flag & 0x4 == 4:
            pProcessData.one_mate_unmapped += 1
            continue
        # skip if the read quality is low
        if mate1.mapq < pMinMappingQuality or mate2.mapq < pMinMappingQuality:
            # for bwa other way to test
            # for multi-mapping reads is with a mapq = 0
            # the XS flag is not reliable.
            if mate1.mapq == 0 & mate2.mapq == 0:
                pProcessData.one_mate_not_unique += 1
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

            pProcessData.one_mate_low_quality += 1
            continue

        if pSkipDuplicationCheck is False:
            if pReadPosMatrix.is_duplicated(pRefId2name[mate1.rname],
                                            mate1.pos,
                                            pRefId2name[mate2.rname],
                                            mate2.pos):
                pProcessData.duplicated_pairs += 1
                continue
        buffer_mate1.append(mate1)
        buffer_mate2.append(mate2)
        j += 1

    pProcessData.iter_num += iter_num - len(buffer_mate1)
    if all_data_read and len(buffer_mate1) != 0 and len(buffer_mate2) != 0:
        return buffer_mate1, buffer_mate2, True, pProcessData
    if all_data_read and len(buffer_mate1) == 0 or len(buffer_mate2) == 0:
        return None, None, True, pProcessData
    return buffer_mate1, buffer_mate2, False, pProcessData


def process_data(pSharedData, pProcessData, pCol, pRow, pData, pQueue):
    """This function is used to compute the data in parallel. To do this it is called at 
    the process start together with a buffer of input data and some parameters. This all is
    embedded in the parameter pData which is of the datatype \'MultiprocessingData\'."""
    iter_num = 0
    pair_added = 0
    # hic_matrix = None
    out_bam = pysam.Samfile(os.path.join(pSharedData.outputFileBufferDir, pProcessData.output_name), 'wb', template=pSharedData.template)

    if pProcessData.buffer1 is None or pProcessData.buffer2 is None:
        pQueue.put(pair_added)

    while iter_num < len(pProcessData.buffer1) and iter_num < len(pProcessData.buffer2):
        mate1 = pProcessData.buffer1[iter_num]
        mate2 = pProcessData.buffer2[iter_num]
        iter_num += 1

        # check if reads belong to a bin
        #
        # pDictBinInterval stores the start and end position for each chromsome in the array 'pSharedBinIntvalTree'
        # To get to the right interval a binary search is used.
        mate_bins = []
        mate_is_unasigned = False
        for mate in [mate1, mate2]:
            mate_ref = pSharedData.refId2name[mate.rname]
            # find the middle genomic position of the read. This is used to find the bin it belongs to.
            read_middle = mate.pos + int(mate.qlen / 2)
            try:
                start, end = pSharedData.dictBinIntervalTreeIndex[mate_ref]
                middle_pos = int((start + end) / 2)
                mate_bin = None
                while not start > end:
                    if pSharedData.sharedBinIntvalTree[middle_pos].begin <= read_middle and read_middle <= pSharedData.sharedBinIntvalTree[middle_pos].end:
                        mate_bin = pSharedData.sharedBinIntvalTree[middle_pos]
                        mate_is_unasigned = False
                        break
                    elif pSharedData.sharedBinIntvalTree[middle_pos].begin > read_middle:
                        end = middle_pos - 1
                        middle_pos = int((start + end) / 2)
                        mate_is_unasigned = True
                    else:
                        start = middle_pos + 1
                        middle_pos = int((start + end) / 2)
                        mate_is_unasigned = True

            except:
                # for small contigs it can happen that they are not
                # in the bin_intval_tree keys if no restriction site is found on the contig.
                mate_is_unasigned = True
                break

            # report no match case
            if mate_bin is None:
                mate_is_unasigned = True
                break
            mate_bin_id = mate_bin.data
            mate_bins.append(mate_bin_id)

        # if a mate is unassigned, it means it is not close
        # to a restriction sites
        if mate_is_unasigned is True:
            pProcessData.mate_not_close_to_rf += 1
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
            # self circles are defined as pairs within 25kb
            # with 'outward' orientation (Jin et al. 2013. Nature)
            if abs(mate2.pos - mate1.pos) < 25000 and orientation == 'outward':
                pProcessData.self_circle += 1
                if pSharedData.removeSelfCircles:
                    continue

            # check for dangling ends if the restriction sequence
            # is known:
            if pSharedData.restrictionSequence:
                if check_dangling_end(mate1, pSharedData.danglingSequences) or \
                        check_dangling_end(mate2, pSharedData.danglingSequences):
                    pProcessData.dangling_end += 1
                    continue

            if abs(mate2.pos - mate1.pos) < 1000 and orientation == 'inward':
                has_rf = []

                if pSharedData.rfPositions and pSharedData.restrictionSequence:
                    # check if in between the two mate
                    # ends the restriction fragment is found.

                    # the interval used is:
                    # start of fragment + length of restriction sequence
                    # end of fragment - length of restriction sequence
                    # the restriction sequence length is subtracted
                    # such that only fragments internally containing
                    # the restriction site are identified
                    frag_start = min(mate1.pos, mate2.pos) + len(pSharedData.restrictionSequence)
                    frag_end = max(mate1.pos + mate1.qlen, mate2.pos + mate2.qlen) - len(pSharedData.restrictionSequence)
                    mate_ref = pSharedData.refId2name[mate1.rname]
                    has_rf = sorted(pSharedData.rfPositions[mate_ref][frag_start: frag_end])

                # case when there is no restriction fragment site between the mates
                if len(has_rf) == 0:
                    pProcessData.same_fragment += 1
                    continue

                pProcessData.self_ligation += 1

                if pSharedData.removeSelfLigation:
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
            pProcessData.inter_chromosomal += 1

        elif abs(mate2.pos - mate1.pos) < 20000:
            pProcessData.short_range += 1
        else:
            pProcessData.long_range += 1

        if orientation == 'inward':
            pProcessData.count_inward += 1
        elif orientation == 'outward':
            pProcessData.count_outward += 1
        elif orientation == 'same-strand-left':
            pProcessData.count_left += 1
        elif orientation == 'same-strand-right':
            pProcessData.count_right += 1

        for mate in [mate1, mate2]:
            # fill in coverage vector
            vec_start = max(0, mate.pos - mate_bin.begin) / pSharedData.binsize
            length_coverage = pSharedData.coverageIndex[mate_bin_id].end - pSharedData.coverageIndex[mate_bin_id].begin
            vec_end = min(length_coverage, vec_start +
                          len(mate.seq) / pSharedData.binsize)
            coverage_index = pSharedData.coverageIndex[mate_bin_id].begin + vec_start
            coverage_end = pSharedData.coverageIndex[mate_bin_id].begin + vec_end
            for i in xrange(coverage_index, coverage_end, 1):
                pSharedData.coverage[i] += 1

        pRow[pair_added] = mate_bins[0]
        pCol[pair_added] = mate_bins[1]
        pData[pair_added] = np.uint8(1)

        pair_added += 1
        # prepare data for bam output
        # set the flag to point that this data is paired
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
        out_bam.write(mate1)
        out_bam.write(mate2)

    out_bam.close()
    pProcessData.pair_added += pair_added
    pQueue.put([pair_added, pProcessData])
    return


def write_bam(pOutputBam, pTemplate, pOutputFileBufferDir):
    """This function is used by the background process to merge the sub-output bam files
    of each input buffer. If it is done, it is creating a file named \'done_processing\' to signal 
    the main process it can continue."""
    out_bam = pysam.Samfile(os.path.join(pOutputFileBufferDir, pOutputBam), 'wb', template=pTemplate)

    counter = 0
    while not os.path.isfile(os.path.join(pOutputFileBufferDir, 'done_processing')) or os.path.isfile(os.path.join(pOutputFileBufferDir, str(counter), '.bam_done')):
        if os.path.isfile(os.path.join(pOutputFileBufferDir, str(counter), '.bam_done')):
            out_put_threads = pysam.Samfile(os.path.join(pOutputFileBufferDir, str(counter), '.bam'), 'rb')
            while True:
                try:
                    data = out_put_threads.next()
                except StopIteration:
                    break
                out_bam.write(data)
            out_put_threads.close()
            os.remove(os.path.join(pOutputFileBufferDir, str(counter), '.bam'))
            os.remove(os.path.join(pOutputFileBufferDir, str(counter), '.bam_done'))
            counter += 1
        else:
            time.sleep(3)

    out_bam.close()
    os.remove(pOutputFileBufferDir + 'done_processing')


def main(args=None):
    """
    Reads line by line two bam files that are not sorted.
    Each line in the two bam files should correspond
    to the mapped position of the two ends of a Hi-C
    fragment. 

    Each mate pair is assessed to determine if it is
    a valid Hi-C pair, in such case a matrix
    reporting the counts of mates is constructed.

    A bam file containing the valid Hi-C reads
    is also constructed
    """

    args = parse_arguments().parse_args(args)
    if args.threads < 3:
        exit("\nAt least three threads need to be defined.\n")

    sys.stderr.write("reading {} and {} to build hic_matrix\n".format(args.samFiles[0].name,
                                                                      args.samFiles[1].name))
    str1 = pysam.Samfile(args.samFiles[0].name, 'rb')
    str2 = pysam.Samfile(args.samFiles[1].name, 'rb')

    args.samFiles[0].close()
    args.samFiles[1].close()
    outputFileBufferDir = args.outputFileBufferDir

    if args.outBam:
        args.outBam.close()

    chrom_sizes = get_chrom_sizes(str1)

    read_pos_matrix = ReadPositionMatrix()

    # define bins
    rf_positions = None
    if args.restrictionCutFile:
        rf_interval = bed2interval_list(args.restrictionCutFile)
        bin_intervals = get_rf_bins(rf_interval,
                                    min_distance=args.minDistance,
                                    max_distance=args.maxDistance)

        rf_positions = intervalListToIntervalTree(rf_interval)
    else:
        bin_intervals = get_bins(args.binSize, chrom_sizes, args.region)

    matrix_size = len(bin_intervals)
    bin_intval_tree = intervalListToIntervalTree(bin_intervals)
    ref_id2name = str1.references

    # build c_type shared memory for the interval tree
    shared_array_list = []
    index_dict = {}
    j = 0
    interval_list = []
    for i, seq in enumerate(bin_intval_tree):
        start = j
        interval_list = []
        for interval in bin_intval_tree[seq]:
            interval_list.append((interval.begin, interval.end, interval.data))
            j += 1
        end = j - 1
        index_dict[seq] = (start, end)
        interval_list = sorted(interval_list)
        shared_array_list.extend(interval_list)
    shared_build_intval_tree = RawArray(C_Interval, shared_array_list)
    bin_intval_tree = None
    dangling_sequences = dict()
    if args.restrictionSequence:
        # build a list of dangling sequences
        args.restrictionSequence = args.restrictionSequence.upper()
        dangling_sequences['pat_forw'] = args.restrictionSequence[1:]
        dangling_sequences['pat_rev'] = args.restrictionSequence[:-1]
        sys.stderr.write("dangling sequences to check "
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

    for value in bin_intervals:
        chrom, start, end = value
        start_pos_coverage.append(number_of_elements_coverage)

        number_of_elements_coverage += (end - start) / binsize
        end_pos_coverage.append(number_of_elements_coverage - 1)
    pos_coverage_ = zip(start_pos_coverage, end_pos_coverage)
    pos_coverage = RawArray(C_Coverage, pos_coverage_)
    pos_coverage_ = None
    start_pos_coverage = None
    end_pos_coverage = None
    coverage = Array(c_uint, number_of_elements_coverage)
    # foo = RawArray(pysam.libcalignedsegment.AlignedSegment, 50)
    # print foo
    # define global shared ctypes arrays for row, col and data
    args.threads = args.threads - 2
    row = [None] * args.threads
    col = [None] * args.threads
    data = [None] * args.threads
    for i in xrange(args.threads):
        row[i] = RawArray(c_uint, args.inputBufferSize)
        col[i] = RawArray(c_uint, args.inputBufferSize)
        data[i] = RawArray(c_ushort, args.inputBufferSize)

    start_time = time.time()

    # iter_num = 0
    # pair_added = 0
    hic_matrix = None


    # pair_added = 0

    # input buffer for bam files
    buffer_workers1 = [None] * args.threads
    buffer_workers2 = [None] * args.threads
    process_data_obj = [None] * args.threads
    # output buffer to write bam with mate1 and mate2 pairs
    process = [None] * args.threads
    all_data_processed = False
    hic_matrix = coo_matrix((matrix_size, matrix_size), dtype='uint32')
    queue = [None] * args.threads

    all_threads_done = False
    thread_done = [False] * args.threads
    count_output = 0
    count_call_of_read_input = 0
    main_process_data = ProcessData()

    shared_data = SharedData(pMinMappingQuality=args.minMappingQuality,
                             pRemoveSelfCircles=args.removeSelfCircles,
                             pRestrictionSequence=args.restrictionSequence,
                             pRemoveSelfLigation=args.removeSelfLigation,
                             pReadPosMatrix=read_pos_matrix,
                             pRfPositions=rf_positions,
                             pRefId2name=ref_id2name,
                             pDanglingSequences=dangling_sequences,
                             pBinsize=binsize,
                             pTemplate=str1,
                             pSharedBinIntvalTree=shared_build_intval_tree,
                             pDictBinIntervalTreeIndex=index_dict,
                             pCoverage=coverage,
                             pCoverageIndex=pos_coverage,
                             pOutputFileBufferDir=outputFileBufferDir)

    process_write_bam_file = Process(target=write_bam, kwargs=dict(pOutputBam=args.outBam.name, pTemplate=str1, pOutputFileBufferDir=outputFileBufferDir))
    process_write_bam_file.start()
    while not all_data_processed or not all_threads_done:

        for i in xrange(args.threads):
            if queue[i] is None and not all_data_processed:
                count_call_of_read_input += 1

                mate_buffer1, mate_buffer2, all_data_processed, main_process_data = readBamFiles(pFileOneIterator=str1,
                                                                                                 pFileTwoIterator=str2,
                                                                                                 pNumberOfItemsPerBuffer=args.inputBufferSize,
                                                                                                 pSkipDuplicationCheck=args.skipDuplicationCheck,
                                                                                                 pReadPosMatrix=read_pos_matrix,
                                                                                                 pRefId2name=ref_id2name,
                                                                                                 pMinMappingQuality=args.minMappingQuality,
                                                                                                 pProcessData=main_process_data)

                queue[i] = Queue()
                process_data_obj[i] = ProcessData(pBuffer1=mate_buffer1, pBuffer2=mate_buffer2,
                                                  pResultIndex=i, pOutputName=str(count_output) + '.bam',
                                                  pCounter=count_output)
                thread_done[i] = False
                process[i] = Process(target=process_data, kwargs=dict(pSharedData=shared_data,
                                                                      pProcessData=process_data_obj[i],
                                                                      pCol=col[i],
                                                                      pRow=row[i],
                                                                      pData=data[i],
                                                                      pQueue=queue[i]))
                process[i].start()
                count_output += 1
                buffer_workers1[i] = None
                buffer_workers2[i] = None
                process_data_obj[i] = None
                
            elif queue[i] is not None and not queue[i].empty():
                result = queue[i].get()

                if result[0] is not None:
                    elements = result[0]
                    if hic_matrix is None:
                        hic_matrix = coo_matrix((data[i][:elements], (row[i][:elements], col[i][:elements])), shape=(matrix_size, matrix_size), dtype='uint32')
                    else:
                        hic_matrix += coo_matrix((data[i][:elements], (row[i][:elements], col[i][:elements])), shape=(matrix_size, matrix_size), dtype='uint32')

                    main_process_data.synchronize(result[1])

                queue[i] = None
                process[i].join()
                process[i].terminate()
                process[i] = None
                thread_done[i] = True

                open(os.path.join(outputFileBufferDir, str(result[1].counter) + '.bam_done'), 'a').close()
                result = None
                # caused by the architecture I try to display this output information after +-1e5 of 1e6 reads.
                if main_process_data.iter_num % 1e6 < 100000:
                    elapsed_time = time.time() - start_time
                    sys.stderr.write("processing {} lines took {:.2f} "
                                     "secs ({:.1f} lines per "
                                     "second)\n".format(main_process_data.iter_num,
                                                        elapsed_time,
                                                        main_process_data.iter_num / elapsed_time))
                    sys.stderr.write("{} ({:.2f}%) valid pairs added to matrix"
                                     "\n".format(main_process_data.pair_added, float(100 * main_process_data.pair_added) / main_process_data.iter_num))
                if args.doTestRun and main_process_data.iter_num > 1e5:
                    sys.stderr.write("\n## *WARNING*. Early exit because of --doTestRun parameter  ##\n\n")
                    break
            else:
                time.sleep(1)

        if all_data_processed:
            all_threads_done = True
            for thread in thread_done:
                if not thread:
                    all_threads_done = False

    # the resulting matrix is only filled unevenly with some pairs
    # int the upper triangle and others in the lower triangle. To construct
    # the definite matrix I add the values from the upper and lower triangles
    # and subtract the diagonal to avoid double counting it.
    # The resulting matrix is symmetric.
    open(os.path.join(outputFileBufferDir, 'done_processing'), 'a').close()
    print "wait for bam merging process to finish"
    process_write_bam_file.join()
    print "wait for bam merging process to finish...DONE!"

    dia = dia_matrix(([hic_matrix.diagonal()], [0]), shape=hic_matrix.shape)
    hic_matrix = hic_matrix + hic_matrix.T - dia
    # extend bins such that they are next to each other
    bin_intervals = enlarge_bins(bin_intervals[:], chrom_sizes)
    # compute max bin coverage
    bin_max = []

    for cover in pos_coverage:
        max_element = 0
        for i in xrange(cover.begin, cover.end, 1):
            if coverage[i] > max_element:
                max_element = coverage[i]
        if max_element == 0:
            bin_max.append(np.nan)
        else:
            bin_max.append(max_element)

    chr_name_list, start_list, end_list = zip(*bin_intervals)
    bin_intervals = zip(chr_name_list, start_list, end_list, bin_max)
    hic_ma = hm.hiCMatrix()
    hic_ma.setMatrix(hic_matrix, cut_intervals=bin_intervals)

    args.outFileName.close()
    # removing the empty file. Otherwise the save method
    # will say that the file already exists.
    unlink(args.outFileName.name)

    hic_ma.save(args.outFileName.name)

    if args.outBam:
        shutil.move(outputFileBufferDir + args.outBam.name, args.outBam.name)

    """
    if args.restrictionCutFile:
        # load the matrix to mask those
        # bins that most likely didn't
        # have a restriction site that was cutted

        # reload the matrix as a HiCMatrix object
        hic_matrix = hm.hiCMatrix(args.outFileName.name)

        hic_matrix.maskBins(get_poor_bins(bin_max))
        hic_matrix.save(args.outFileName.name)
    """
    if args.removeSelfLigation:
        msg = " (removed)"
    else:
        msg = " (not removed)"

    mappable_pairs = main_process_data.iter_num - main_process_data.one_mate_unmapped

    log_file_name = os.path.join(args.QCfolder, "QC.log")
    log_file = open(log_file_name, 'w')
    log_file.write("""
File\t{}\t\t
Pairs considered\t{}\t\t
Min rest. site distance\t{}\t\t
Max rest. site distance\t{}\t\t

""".format(args.outFileName.name, main_process_data.iter_num, args.minDistance,
           args.maxDistance))

    log_file.write("Pairs used\t{}\t({:.2f})\t({:.2f})\n".format(main_process_data.pair_added,
                                                                 100 * float(main_process_data.pair_added) / main_process_data.iter_num,
                                                                 100 * float(main_process_data.pair_added) / mappable_pairs))
    log_file.write("One mate unmapped\t{}\t({:.2f})\t({:.2f})\n".format(main_process_data.one_mate_unmapped,
                                                                        100 * float(main_process_data.one_mate_unmapped) / main_process_data.iter_num,
                                                                        100 * float(main_process_data.one_mate_unmapped) / mappable_pairs))

    log_file.write("One mate not unique\t{}\t({:.2f})\t({:.2f})\n".format(main_process_data.one_mate_not_unique,
                                                                          100 * float(main_process_data.one_mate_not_unique) / main_process_data.iter_num,
                                                                          100 * float(main_process_data.one_mate_not_unique) / mappable_pairs))

    log_file.write("One mate low quality\t{}\t({:.2f})\t({:.2f})\n".format(main_process_data.one_mate_low_quality,
                                                                           100 * float(main_process_data.one_mate_low_quality) / main_process_data.iter_num,
                                                                           100 * float(main_process_data.one_mate_low_quality) / mappable_pairs))

    log_file.write("dangling end\t{}\t({:.2f})\t({:.2f})\n".format(main_process_data.dangling_end,
                                                                   100 * float(main_process_data.dangling_end) / main_process_data.iter_num,
                                                                   100 * float(main_process_data.dangling_end) / mappable_pairs))

    log_file.write("self ligation{}\t{}\t({:.2f})\t({:.2f})\n".format(msg, main_process_data.self_ligation,
                                                                      100 * float(main_process_data.self_ligation) / main_process_data.iter_num,
                                                                      100 * float(main_process_data.self_ligation) / mappable_pairs))

    log_file.write("One mate not close to rest site\t{}\t({:.2f})\t({:.2f})\n".format(main_process_data.mate_not_close_to_rf,
                                                                                      100 * float(main_process_data.mate_not_close_to_rf) / main_process_data.iter_num,
                                                                                      100 * float(main_process_data.mate_not_close_to_rf) / mappable_pairs))

    log_file.write("same fragment (800 bp)\t{}\t({:.2f})\t({:.2f})\n".format(main_process_data.same_fragment,
                                                                             100 * float(main_process_data.same_fragment) / main_process_data.iter_num,
                                                                             100 * float(main_process_data.same_fragment) / mappable_pairs))
    log_file.write("self circle\t{}\t({:.2f})\t({:.2f})\n".format(main_process_data.self_circle,
                                                                  100 * float(main_process_data.self_circle) / main_process_data.iter_num,
                                                                  100 * float(main_process_data.self_circle) / mappable_pairs))
    log_file.write("duplicated pairs\t{}\t({:.2f})\t({:.2f})\n".format(main_process_data.duplicated_pairs,
                                                                       100 * float(main_process_data.duplicated_pairs) / main_process_data.iter_num,
                                                                       100 * float(main_process_data.duplicated_pairs) / mappable_pairs))
    if main_process_data.pair_added > 0:
        log_file.write("Of pairs used:\n")
        log_file.write("inter chromosomal\t{}\t({:.2f})\n".format(main_process_data.inter_chromosomal, 100 * float(main_process_data.inter_chromosomal) / main_process_data.pair_added))

        log_file.write("short range < 20kb\t{}\t({:.2f})\n".format(main_process_data.short_range, 100 * float(main_process_data.short_range) / main_process_data.pair_added))

        log_file.write("long range\t{}\t({:.2f})\n".format(main_process_data.long_range, 100 * float(main_process_data.long_range) / main_process_data.pair_added))

        log_file.write("inward pairs\t{}\t({:.2f})\n".format(main_process_data.count_inward, 100 * float(main_process_data.count_inward) / main_process_data.pair_added))

        log_file.write("outward pairs\t{}\t({:.2f})\n".format(main_process_data.count_outward, 100 * float(main_process_data.count_outward) / main_process_data.pair_added))

        log_file.write("left pairs\t{}\t({:.2f})\n".format(main_process_data.count_left, 100 * float(main_process_data.count_left) / main_process_data.pair_added))

        log_file.write("right pairs\t{}\t({:.2f})\n".format(main_process_data.count_right, 100 * float(main_process_data.count_right) / main_process_data.pair_added))

    log_file.close()
    QC.main("-l {} -o {}".format(log_file_name, args.QCfolder).split())


class Tester(object):
    def __init__(self):
        hic_test_data_dir = os.environ.get('HIC_TEST_DATA_DIR', False)
        if hic_test_data_dir:
            self.root = hic_test_data_dir
        else:
            self.root = os.path.dirname(os.path.abspath(__file__)) + "/test/test_data/"
        self.bam_file_1 = os.path.join(self.root, "hic.bam")
