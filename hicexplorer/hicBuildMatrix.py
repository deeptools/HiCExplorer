import argparse
import sys
import numpy as np
from scipy.sparse import coo_matrix, dia_matrix, dok_matrix, csr_matrix
import time
from os import unlink
import os
from fs.memoryfs import MemoryFS

# fs.memoryfs.MemoryFS
from itertools import izip
import pysam
# bx python
from intervaltree import IntervalTree, Interval

# own tools
from hicexplorer import HiCMatrix as hm
from hicexplorer.utilities import getUserRegion, genomicRegion
from hicexplorer._version import __version__


import collections  # for the buffer
# import time  # to ease polling
# import threading
import multiprocessing
from multiprocessing import Lock
from multiprocessing.sharedctypes import Array
from ctypes import Structure, c_int, c_bool

from pyhashxx import hashxx
debug = 1
from multiprocessing.managers import SyncManager

from multiprocessing.managers import MakeProxyType

BaseSetProxy = MakeProxyType('BaseSetProxy', [
    '__and__', '__contains__', '__iand__', '__ior__',
    '__isub__', '__ixor__', '__len__', '__or__', '__rand__', '__ror__', '__rsub__',
    '__rxor__', '__sub__', '__xor__', 'add', 'clear', 'copy', 'difference',
    'difference_update', 'discard', 'intersection', 'intersection_update', 'isdisjoint',
    'issubset', 'issuperset', 'pop', 'remove', 'symmetric_difference',
    'symmetric_difference_update', 'union', 'update']
    )


class SetProxy(BaseSetProxy):
    # in-place hooks need to return `self`, specify these manually
    def __iand__(self, value):
        self._callmethod('__iand__', (value,))
        return self
    def __ior__(self, value):
        self._callmethod('__ior__', (value,))
        return self
    def __isub__(self, value):
        self._callmethod('__isub__', (value,))
        return self
    def __ixor__(self, value):
        self._callmethod('__ixor__', (value,))
        return self







# class Duplicate(Structure):
#     _fields_ = [('hash', c_int), ('value', c_bool)]
class ReadPositionMatrix(object):
    """ class to check for PCR duplicates.
    A sparse matrix having as bins all possible
    start sites (single bp resolution)
    is created. PCR duplicates
    are determined by checking if the matrix
    cell is already filled.

    """

    def __init__(self, pManager):
        """
        >>> rp = ReadPositionMatrix([('1', 10), ('2', 10)])
        >>> rp.pos2matrix_bin('1', 0)
        0
        >>> rp.pos2matrix_bin('2', 0)
        10
        >>> rp.is_duplicated('1', 0, '2', 0)
        False
        >>> rp.pos_matrix[0,10]
        True
        >>> rp.pos_matrix[10,0]
        True
        >>> rp.is_duplicated('1', 0, '2', 0)
        True
        """
        # determine number of bins
        # total_size = 0
        # m = multiprocessing.Manager()
        # m.
        # self.chr_start_pos = pManager.dict()
        # for chrom, size in chrom_sizes:
        #     self.chr_start_pos[chrom] = total_size
        #     total_size += size

        self.pos_matrix = pManager.set()
        # self.pos_matrix = dok_matrix((total_size, total_size), dtype=bool)
        # print self.pos_matrix

    def is_duplicated(self, chrom1, start1, chrom2, start2):
       
        # pos1 = hashxx(chrom1, seed=start1)
        # pos2 = hashxx(chrom2, seed=start2)
        
        # chrom = chrom1 if chrom1 < chrom2 else chrom2
        # pos = start1 if start1 < start2 else start2
        # pos3 = hashxx(int(pos1))
        # pos4 = hashxx(int(pos2))
        # pos = pos1 + pos2
        # pos = hashxx(), seed=)
        # hash(start1)
        # hash(start2)
        id_string = "%s%s-%s%s" % (chrom1, start1, chrom2, start2)
        # chrom1 + str(start1) + '-' + str(chrom2) + str(start2)
        if id_string in self.pos_matrix:
            return True
        else:
            self.pos_matrix.add(id_string)
            self.pos_matrix.add("%s%s-%s%s" % (chrom2, start2, chrom1, start1))
            # self.pos_matrix.add(hash(pos2) + pos1)
            
            # self.pos_matrix.add(hashxx(hashxx(chrom2, seed=start2), seed=hashxx(chrom1, seed=start1)))
            
        # pos1 = self.pos2matrix_bin(chrom1, start1)
        # pos2 = self.pos2matrix_bin(chrom2, start2)

        # if str(pos1) + '-' + str(pos2) in self.pos_matrix:
        #     return True
        # else:
        #     self.pos_matrix[str(pos1) + '-' + str(pos2)] = True
        #     self.pos_matrix[str(pos2) + '-' + str(pos1)] = True
        #     return False

    def pos2matrix_bin(self, chrom, start):
        return self.chr_start_pos[chrom] + start

    def add_matrix(self, pReadPosMatrix):
        self.pos_matrix += pReadPosMatrix.pos_matrix
    # def to_coo(self):
    #     self.pos_matrix = self.pos_matrix.tocoo()
    # def to_dok(self):
    #     self.pos_matrix = self.pos_matrix.todok()


# class ReadPositionMatrix(object):
#     """ class to check for PCR duplicates.
#     A sparse matrix having as bins all possible
#     start sites (single bp resolution)
#     is created. PCR duplicates
#     are determined by checking if the matrix
#     cell is already filled.

#     """

#     def __init__(self, chrom_sizes):
#         """
#         >>> rp = ReadPositionMatrix([('1', 10), ('2', 10)])
#         >>> rp.pos2matrix_bin('1', 0)
#         0
#         >>> rp.pos2matrix_bin('2', 0)
#         10
#         >>> rp.is_duplicated('1', 0, '2', 0)
#         False
#         >>> rp.pos_matrix[0,10]
#         True
#         >>> rp.pos_matrix[10,0]
#         True
#         >>> rp.is_duplicated('1', 0, '2', 0)
#         True
#         """
#         # determine number of bins
#         total_size = 0

#         self.chr_start_pos = {}
#         for chrom, size in chrom_sizes:
#             self.chr_start_pos[chrom] = total_size
#             total_size += size

#         # self.pos_matrix = pManager.dict()
#         self.pos_matrix = dok_matrix((total_size, total_size), dtype=bool)
#         # print self.pos_matrix

#     def is_duplicated(self, chrom1, start1, chrom2, start2):
#         pos1 = self.pos2matrix_bin(chrom1, start1)
#         pos2 = self.pos2matrix_bin(chrom2, start2)

#         if self.pos_matrix[pos1, pos2]:
#             return True
#         else:
#             self.pos_matrix[pos1, pos2] = True
#             self.pos_matrix[pos2, pos1] = True
#             return False

#     def pos2matrix_bin(self, chrom, start):
#         return self.chr_start_pos[chrom] + start

#     def add_matrix(self, pReadPosMatrix):
#         self.pos_matrix += pReadPosMatrix.pos_matrix

#     def init_matrix_shared(self, pData, pIndices, pIndptr, pShape):
#         self.pos_matrix = dok_matrix((pData, pIndices, pIndptr), shape=pShape, copy=False)
#     def get_arrays(self):
#         return self.pos_matrix.data, self.pos_matrix.indices, self.pos_matrix.indptr
#     # def to_dok(self):
#     #     self.pos_matrix = self.pos_matrix.todok()


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = ('Using an alignment from Bowtie2 where both '
                     'PE reads are mapped using  the --local '
                     'option, this program reads such file and '
                     'creates a matrix of interactions.'))

    # define the arguments
    parser.add_argument('--samFiles', '-s',
                        help = 'The two sam files to process',
                        metavar = 'two sam files',
                        nargs = 2,
                        type = argparse.FileType('r'),
                        required = True)

    # define the arguments
    parser.add_argument('--outBam', '-b',
                        help = 'Bam file to process',
                        metavar = 'bam file',
                        type = argparse.FileType('w'),
                        required = True)

    group=parser.add_mutually_exclusive_group(required = True)

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
                        help='Number of threads. Using the python multiprocessing module.',
                        required=False,
                        default=4,
                        type=int
                        )
    parser.add_argument('--inputBufferSize',
                        help='Size of the input buffer of each thread. One million elements per input file per thread is the default value.'
                             ' Reduce value to decrease memory usage.',
                        required=False,
                        default=1e6,
                        type=int
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

        for i, _sup in enumerate(pysam_obj):
            if i > len(other_alignments):
                break
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


def readBamFiles(pFileOneIterator, pFileTwoIterator, pNumberOfItemsPerBuffer):
    """Read the two bam input files into n buffers each with pNumberOfItemsPerBuffer with n = number of processes """
    buffer_mate1 = []
    buffer_mate2 = []

    iter_num = 0
    all_data_read = False
    j = 0
    # print pFileOneIterator
    # print pFileTwoIterator
    # for mate1, mate2 in izip(pFileOneIterator, pFileTwoIterator):
    #     if iter_num >= pNumberOfItemsPerBuffer:
    #         break
    #     # if mate1.flag & 256 == 256:
    #     #     for mate1 in pFileOneIterator:
    #     #         if not (mate1.flag & 256 == 256):
    #     #             break
    #     # if mate2.flag & 256 == 256:
    #     #     for mate2 in pFileTwoIterator:
    #     #         if not (mate2.flag & 256 == 256):
    #     #             break

    #     while mate1.flag & 256 == 256:
    #         try:
    #             mate1 = pFileOneIterator.next()
    #         except StopIteration:
    #             # pTerminateSignal.set()
    #             all_data_read = True
    #             print "all_data_processed: ", all_data_read
    #             break
    #     while mate2.flag & 256 == 256:
    #         try:
    #             mate2 = pFileTwoIterator.next()
    #         except StopIteration:
    #             # pTerminateSignal.set()
    #             all_data_read = True
    #             print "all_data_processed: ", all_data_read
    #             break
    #     assert mate1.qname == mate2.qname, "FATAL ERROR {} {} " \
    #         "Be sure that the sam files have the same read order " \
    #         "If using Bowtie2 or Hisat2 add " \
    #         "the --reorder option".format(mate1.qname, mate2.qname)

    #     # check for supplementary alignments
    #     # (needs to be done before skipping any unmapped reads
    #     # to keep the order of the two bam files in sync)
    #     mate1_supplementary_list = get_supplementary_alignment(mate1, pFileOneIterator)
    #     mate2_supplementary_list = get_supplementary_alignment(mate2, pFileTwoIterator)

    #     if mate1_supplementary_list:
    #         mate1 = get_correct_map(mate1, mate1_supplementary_list)

    #     if mate2_supplementary_list:
    #         mate2 = get_correct_map(mate2, mate2_supplementary_list)

    #     buffer_mate1.append(mate1)
    #     buffer_mate2.append(mate2)
    #     iter_num += 1
    # print "iter_num: ", iter_num
    # print "pNumberOfItemsPerBuffer: ", pNumberOfItemsPerBuffer
    # if iter_num < pNumberOfItemsPerBuffer and len(buffer_mate1) != 0:
    #     return buffer_mate1, buffer_mate2, True, iter_num
    # if iter_num == 0:
    #     return None, None, True, 0
    # return buffer_mate1, buffer_mate2, False, iter_num
    while j < pNumberOfItemsPerBuffer:
        iter_num += 1
        try:
            mate1 = pFileOneIterator.next()
            mate2 = pFileTwoIterator.next()
        except StopIteration:
            # pTerminateSignal.set()
            all_data_read = True
            print "all_data_processed: ", all_data_read
            break

        # skip 'not primary' alignments
        while mate1.flag & 256 == 256:
            try:
                mate1 = pFileOneIterator.next()
            except StopIteration:
                # pTerminateSignal.set()
                all_data_read = True
                print "all_data_processed: ", all_data_read
                break
        while mate2.flag & 256 == 256:
            try:
                mate2 = pFileTwoIterator.next()
            except StopIteration:
                # pTerminateSignal.set()
                all_data_read = True
                print "all_data_processed: ", all_data_read
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

        buffer_mate1.append(mate1)
        buffer_mate2.append(mate2)
        j += 1

    if all_data_read and len(buffer_mate1) != 0 and len(buffer_mate2) != 0:
        return buffer_mate1, buffer_mate2, True, iter_num
    if all_data_read and len(buffer_mate1) == 0 or len(buffer_mate2) == 0:
        return None, None, True, 0
    return buffer_mate1, buffer_mate2, False, iter_num


def process_data(pMateBuffer1, pMateBuffer2, pMinMappingQuality, pSkipDuplicationCheck,
                 pRemoveSelfCircles, pRestrictionSequence, pRemoveSelfLigation, pMatrixSize,
                 pReadPosMatrix, pRfPositions, pBinIntvalTree, pRefId2name,
                #  pRfPositions, pBinIntvalTree, pRefId2name,
                 pDanglingSequences, pBinsize, pResultIndex, pBinIntervals,
                 pQueueOut, pTemplate, pOutputName):

                 
                 #, p
                #  ReadPosMatrixData, pReadPosMatrixIndices, pReadPosMatrixIndptr, pShape):
    # pReadPosMatrix = 
    bufferOutputBam = []
    one_mate_unmapped = 0
    one_mate_low_quality = 0
    one_mate_not_unique = 0
    dangling_end = 0
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

    row = []
    col = []
    data = []

    iter_num = 0
    hic_matrix = None
    # mem_fs = MemoryFS()
    # output_bam_memory = mem_fs.createfile(pOutputName)
    # out_bam = pysam.Samfile(output_bam_memory, 'wb', template=pTemplate)

    coverage = []
    binsize = 10
    # print len(pBinIntervals)
    # out_bam_buffer = []
    # print pBinIntvalTree
    for value in pBinIntervals:
        chrom, start, end = value
        # print value
        # coverage[i] = np.zeros((end - start) / binsize, dtype='uint16')
        coverage.append(np.zeros((end - start) / binsize, dtype='uint16'))
    if pMateBuffer1 is None or pMateBuffer2 is None:

        pQueueOut.put([hic_matrix, [one_mate_unmapped, one_mate_low_quality, one_mate_not_unique, dangling_end, self_circle, self_ligation, same_fragment,
                                    mate_not_close_to_rf, duplicated_pairs, count_inward, count_outward,
                                    count_left, count_right, inter_chromosomal, short_range, long_range, pair_added, iter_num, pResultIndex], coverage, bufferOutputBam])
        return

    while iter_num < len(pMateBuffer1) and iter_num < len(pMateBuffer2):
        mate1 = pMateBuffer1[iter_num]
        mate2 = pMateBuffer2[iter_num]
        iter_num += 1

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

        # check if reads belong to a bin
        mate_bins = []
        mate_is_unasigned = False
        for mate in [mate1, mate2]:
            mate_ref = pRefId2name[mate.rname]
            # find the middle genomic position of the read. This is used to find the bin it belongs to.
            read_middle = mate.pos + int(mate.qlen / 2)
            try:
                mate_bin = sorted(pBinIntvalTree[mate_ref][read_middle:read_middle + 1])
            except KeyError:
                # for small contigs it can happen that they are not
                # in the bin_intval_tree keys if no restriction site is found on the contig.
                mate_is_unasigned = True
                break

            # report no match case
            if len(mate_bin) == 0:
                mate_is_unasigned = True
                break
            # take by default only the first match
            # (although always there should be only
            # one match
            mate_bin = mate_bin[0]

            mate_bin_id = mate_bin.data
            mate_bins.append(mate_bin_id)

        # if a mate is unassigned, it means it is not close
        # to a restriction sites
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
            # self circles are defined as pairs within 25kb
            # with 'outward' orientation (Jin et al. 2013. Nature)
            if abs(mate2.pos - mate1.pos) < 25000 and orientation == 'outward':
                self_circle += 1
                if pRemoveSelfCircles:
                    continue

            # check for dangling ends if the restriction sequence
            # is known:
            if pRestrictionSequence:
                if check_dangling_end(mate1, pDanglingSequences) or \
                        check_dangling_end(mate2, pDanglingSequences):
                    dangling_end += 1
                    continue

            if abs(mate2.pos - mate1.pos) < 1000 and orientation == 'inward':
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
                    frag_start = min(mate1.pos, mate2.pos) + len(pRestrictionSequence)
                    frag_end = max(mate1.pos + mate1.qlen, mate2.pos + mate2.qlen) - len(pRestrictionSequence)
                    mate_ref = pRefId2name[mate1.rname]
                    has_rf = sorted(pRfPositions[mate_ref][frag_start: frag_end])
                    # print

                # case when there is no restriction fragment site between the mates
                if len(has_rf) == 0:
                    # print "foo781"

                    same_fragment += 1
                    continue

                self_ligation += 1
                # print "pRemoveSelfLigation: ", pRemoveSelfLigation

                if pRemoveSelfLigation:
                    # print "foo"
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
            vec_start = max(0, mate.pos - mate_bin.begin) / pBinsize
            vec_end = min(len(coverage[mate_bin_id]), vec_start +
                          len(mate.seq) / pBinsize)
            coverage[mate_bin_id][vec_start:vec_end] += 1

        row.append(mate_bins[0])
        col.append(mate_bins[1])
        data.append(1)

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
        # out_bam.write(mate1)
        # out_bam.write(mate2)

        # import copy
        # bufferOutputBam.append(mate1)
        # bufferOutputBam.append(mate2)

        # if iter_num % 5e6 == 0:
        #     # every 5 million iterations append to the matrix
        #     # otherwise the row, col and data vectors continue growing and
        #     # for a large dataset the system could run out of memory
        #     if hic_matrix is None:
        #         hic_matrix = coo_matrix((data, (row, col)), shape=(pMatrixSize, pMatrixSize), dtype='uint16')
        #     else:
        #         hic_matrix += coo_matrix((data, (row, col)), shape=(pMatrixSize, pMatrixSize), dtype='uint16')
        #     row = []
        #     col = []
        #     data = []
        # iter_num += 1

    if hic_matrix is None:
        hic_matrix = coo_matrix((data, (row, col)), shape=(pMatrixSize, pMatrixSize), dtype='uint16')
    else:
        hic_matrix += coo_matrix((data, (row, col)), shape=(pMatrixSize, pMatrixSize), dtype='uint16')
    row = []
    col = []
    data = []
    print "Process x: ", pResultIndex, " DONE!"
    # print(bufferOutputBam.pop().flag)
    # pTerminateOutput.set()
    # pLock.acquire()
    # pReadPosMatrix.to_coo()
    pQueueOut.put([hic_matrix, [one_mate_unmapped, one_mate_low_quality, one_mate_not_unique, dangling_end, self_circle, self_ligation, same_fragment,
                                mate_not_close_to_rf, duplicated_pairs, count_inward, count_outward,
                                count_left, count_right, inter_chromosomal, short_range, long_range, pair_added, iter_num, pResultIndex],
                   coverage])
    # out_bam.close()
    # pLock.release()
    return


def write_output_bam(pOutputBamName, pBufferOutBam, pTemplate):
    print "write output bam"
    # print "pBufferOutBam", pBufferOutBam
    # for i, buffer in enumerate(pBufferOutBam):
    # print buffer
    # print "buffer x: ", i
    # print "type buffer:", type(pBufferOutBam)
    # print "type file: ", type(pOutputBam)
    out_bam = pysam.Samfile(pOutputBamName, 'wb', template=pTemplate)
    while pBufferOutBam:
        try:
            data = pBufferOutBam.pop()  # pop from RIGHT end of buffer
        except IndexError:
            break
        print "typedata: ", type(data)
        print(data.flag)
        # print "write data"
        # out_bam.write(data)
        # print "wrote data"
    pBufferOutBam.clear()
    out_bam.close()


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

    # check pysam version
    from distutils.version import LooseVersion
    if LooseVersion(pysam.__version__) < LooseVersion("0.8.3"):
        exit("\n*ERROR*\n\nVersion of pysam has to be higher than 0.8.3. Current installed version is {}\n".format(pysam.__version__))

    args = parse_arguments().parse_args(args)

    sys.stderr.write("reading {} and {} to build hic_matrix\n".format(args.samFiles[0].name,
                                                                      args.samFiles[1].name))
    str1 = pysam.Samfile(args.samFiles[0].name, 'rb')
    str2 = pysam.Samfile(args.samFiles[1].name, 'rb')

    args.samFiles[0].close()
    args.samFiles[1].close()
    args.outBam.close()

    # lock = threading.Lock()
    out_bam = pysam.Samfile(args.outBam.name, 'wb', template=str1)
    manager = SyncManager()
    manager.start()
    SyncManager.register('set', set, SetProxy)
    
    lock = multiprocessing.Lock()

    chrom_sizes = get_chrom_sizes(str1)
    # print chrom_sizes
    # initialize read start positions matrix
    manager_multiprocessing = multiprocessing.Manager()
    read_pos_matrix = ReadPositionMatrix( manager_multiprocessing)
    # read_pos_matrix_data, read_pos_matrix_indices, read_pos_matrix_indptr = read_pos_matrix.get_arrays()
    # read_pos_matrix.to_coo()
    # print "read_pos_matrix: ", read_pos_matrix.pos_matrix
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

    sys.stderr.write("Matrix size: {}\n".format(len(bin_intervals)))
    matrix_size = len(bin_intervals)
    bin_intval_tree = intervalListToIntervalTree(bin_intervals)
    ref_id2name = str1.references

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
    coverage = []
    binsize = 10
    for value in bin_intervals:
        chrom, start, end = value
        coverage.append(np.zeros((end - start) / binsize, dtype='uint16'))
    start_time = time.time()

    iter_num = 0
    pair_added = 0
    hic_matrix = None

    one_mate_unmapped = 0
    one_mate_low_quality = 0
    one_mate_not_unique = 0
    dangling_end = 0
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
    buffer_workers1 = [collections.deque()] * args.threads
    buffer_workers2 = [collections.deque()] * args.threads

    # output buffer to write bam with mate1 and mate2 pairs
    process = [None] * args.threads
    all_data_processed = False
    hic_matrix = None
    queue = [None] * args.threads
    all_threads_done = False
    thread_done = [False] * args.threads
    count_output = 0
    count_call_of_read_input = 0
    computed_pairs = 0
    out_bam_in_memory = []
    print "args.removeSelfLigation: ", args.removeSelfLigation
    while not all_data_processed or not all_threads_done:
        # out_q = multiprocessing.Queue()

        for i in xrange(args.threads):
            if queue[i] is None and not all_data_processed:
                count_call_of_read_input += 1

                buffer_workers1[i], buffer_workers2[i], all_data_processed, iter_num_ = readBamFiles(pFileOneIterator=str1,
                                                                                                     pFileTwoIterator=str2,
                                                                                                     pNumberOfItemsPerBuffer=args.inputBufferSize)
                # iter_num += iter_num_
                queue[i] = multiprocessing.Queue()
                thread_done[i] = False
                computed_pairs += len(buffer_workers1[i])
                # create n processes to compute hic matrix
                # for i in xrange(args.threads):
                process[i] = multiprocessing.Process(target=process_data, kwargs=dict(
                    pMateBuffer1=buffer_workers1[i],
                    pMateBuffer2=buffer_workers2[i],
                    pMinMappingQuality=args.minMappingQuality,
                    pSkipDuplicationCheck=args.skipDuplicationCheck,
                    pRemoveSelfCircles=args.removeSelfCircles,
                    pRestrictionSequence=args.restrictionSequence,
                    pRemoveSelfLigation=args.removeSelfLigation,
                    pMatrixSize=matrix_size,
                    pReadPosMatrix=read_pos_matrix,
                    pRfPositions=rf_positions,
                    pBinIntvalTree=bin_intval_tree,
                    pRefId2name=ref_id2name,
                    pDanglingSequences=dangling_sequences,
                    pBinsize=binsize,
                    pResultIndex=i,
                    pBinIntervals=bin_intervals,
                    pQueueOut=queue[i],
                    pTemplate=str1,
                    pOutputName=str(count_output) + '.bam'
                ))
                process[i].start()
                print "Start thread: ", i
                count_output += 1

            elif queue[i] is not None and not queue[i].empty():
                result = queue[i].get()

                print "Get result from process: ", i
                if result[0] is not None:
                    if hic_matrix is None:
                        hic_matrix = result[0]  # hicmatrix
                    else:
                        hic_matrix += result[0]

                    # one_mate_unmapped, one_mate_low_quality, one_mate_not_unique, dangling_end, self_circle, self_ligation, same_fragment,
                    # mate_not_close_to_rf, duplicated_pairs, count_inward, count_outward,
                    # count_left, count_right, inter_chromosomal, short_range, long_range, pair_added, iter_num, pResultIndex
                    one_mate_unmapped += result[1][0]
                    one_mate_low_quality += result[1][1]
                    one_mate_not_unique += result[1][2]
                    dangling_end += result[1][3]
                    self_circle += result[1][4]
                    self_ligation += result[1][5]
                    same_fragment += result[1][6]
                    mate_not_close_to_rf += result[1][7]
                    duplicated_pairs += result[1][8]

                    count_inward += result[1][9]
                    count_outward += result[1][10]
                    count_left += result[1][11]
                    count_right += result[1][12]
                    inter_chromosomal += result[1][13]
                    short_range += result[1][14]
                    long_range += result[1][15]

                    pair_added += result[1][16]
                    iter_num += result[1][17]

                    if result[2] is not None:
                        coverage = np.add(coverage, result[2])

                    # out_bam_in_memory.append(result[3])
                    # read_pos_matrix.add_matrix(result[4])
                # if result[3] is not None:
                    # print type(result[3][0])
                    # print result[3][0].flag

                # buffer_workers1[i] = []
                # buffer_workers2[i] = []
                queue[i] = None
                process[i].join()
                process[i].terminate()
                process[i] = None
                thread_done[i] = True
                print "Get result from process: ", i, " DONE!"
        if all_data_processed:
            all_threads_done = True
            for thread in thread_done:
                if not thread:
                    all_threads_done = False

        # if iter_num % 1e6 == 0:
        #     elapsed_time = time.time() - start_time
        #     sys.stderr.write("processing {} lines took {:.2f} "
        #                      "secs ({:.1f} lines per "
        #                      "second)\n".format(iter_num,
        #                                         elapsed_time,
        #                                         iter_num / elapsed_time))
        #     sys.stderr.write("{} ({:.2f}%) valid pairs added to matrix"
        #                      "\n".format(pair_added, float(100 * pair_added) / iter_num))
        # if args.doTestRun and iter_num > 1e5:
        #     sys.stderr.write("\n## *WARNING*. Early exit because of --doTestRun parameter  ##\n\n")
        #     break

        # for result_buffer in buffer_workers_out:
        #     write_output_bam(out_bam, result_buffer)

    # the resulting matrix is only filled unevenly with some pairs
    # int the upper triangle and others in the lower triangle. To construct
    # the definite matrix I add the values from the upper and lower triangles
    # and subtract the diagonal to avoid double counting it.
    # The resulting matrix is symmetric.
    print "iter_num: ", iter_num
    print "computed_pairs: ", computed_pairs
    dia = dia_matrix(([hic_matrix.diagonal()], [0]), shape=hic_matrix.shape)
    hic_matrix = hic_matrix + hic_matrix.T - dia
    # extend bins such that they are next to each other
    bin_intervals = enlarge_bins(bin_intervals[:], chrom_sizes)
    print "count_call_of_read_input: ", count_call_of_read_input
    # compute max bin coverage
    bin_max = []

    for cov in coverage:
        # bin_coverage.append(round(float(len(cov[cov > 0])) / len(cov), 3))
        if len(cov) == 0:
            bin_max.append(np.nan)
        else:
            bin_max.append(max(cov))

    # print "max coverage value: ", max(bin_max)
    chr_name_list, start_list, end_list = zip(*bin_intervals)
    bin_intervals = zip(chr_name_list, start_list, end_list, bin_max)
    hic_ma = hm.hiCMatrix()
    hic_ma.setMatrix(hic_matrix, cut_intervals=bin_intervals)

    args.outFileName.close()
    # removing the empty file. Otherwise the save method
    # will say that the file already exists.
    unlink(args.outFileName.name)

    hic_ma.save(args.outFileName.name)
    # for memory_file in out_bam_in_memory:
    #     out_put_threads = pysam.Samfile(memory_file, 'rb')
    #     while True:
    #         try:
    #             data = out_put_threads.next()
    #         except StopIteration:
    #             break
    #         out_bam.write(data)
    #     out_put_threads.close()
    #     # os.remove(str(i) + '.bam')
    #     memory_file.
    # out_bam.close()

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

    mappable_pairs = iter_num - one_mate_unmapped
    print("""
File\t{}\t\t
Pairs considered\t{}\t\t
Min rest. site distance\t{}\t\t
Max rest. site distance\t{}\t\t

""".format(args.outFileName.name, iter_num, args.minDistance,
           args.maxDistance))

    print("Pairs used\t{}\t({:.2f})\t({:.2f})".format(pair_added, 100 * float(pair_added) / iter_num,
                                                      100 * float(pair_added) / mappable_pairs))
    print("One mate unmapped\t{}\t({:.2f})\t({:.2f})".format(one_mate_unmapped, 100 * float(one_mate_unmapped) / iter_num,
                                                             100 * float(one_mate_unmapped) / mappable_pairs))

    print("One mate not unique\t{}\t({:.2f})\t({:.2f})".format(one_mate_not_unique, 100 * float(one_mate_not_unique) / iter_num,
                                                               100 * float(one_mate_not_unique) / mappable_pairs))
    print("One mate low quality\t{}\t({:.2f})\t({:.2f})".format(one_mate_low_quality, 100 * float(one_mate_low_quality) / iter_num,
                                                                100 * float(one_mate_low_quality) / mappable_pairs))
    print("dangling end\t{}\t({:.2f})\t({:.2f})".format(dangling_end, 100 * float(dangling_end) / iter_num,
                                                        100 * float(dangling_end) / mappable_pairs))
    print("self ligation{}\t{}\t({:.2f})\t({:.2f})".format(msg, self_ligation, 100 * float(self_ligation) / iter_num,
                                                           100 * float(self_ligation) / mappable_pairs))
    print("One mate not close to rest site\t{}\t({:.2f})\t({:.2f})".format(mate_not_close_to_rf, 100 * float(mate_not_close_to_rf) / iter_num,
                                                                           100 * float(mate_not_close_to_rf) / mappable_pairs))
    print("same fragment (800 bp)\t{}\t({:.2f})\t({:.2f})".format(same_fragment, 100 * float(same_fragment) / iter_num,
                                                                  100 * float(same_fragment) / mappable_pairs))
    print("self circle\t{}\t({:.2f})\t({:.2f})".format(self_circle, 100 * float(self_circle) / iter_num,
                                                       100 * float(self_circle) / mappable_pairs))
    print("duplicated pairs\t{}\t({:.2f})\t({:.2f})".format(duplicated_pairs, 100 * float(duplicated_pairs) / iter_num,
                                                            100 * float(duplicated_pairs) / mappable_pairs))
    if pair_added > 0:
        print("Of pairs used:")
        print("inter chromosomal\t{}\t({:.2f})".format(inter_chromosomal, 100 * float(inter_chromosomal) / pair_added))

        print("short range < 20kb\t{}\t({:.2f})".format(short_range, 100 * float(short_range) / pair_added))

        print("long range\t{}\t({:.2f})".format(long_range, 100 * float(long_range) / pair_added))

        print("inward pairs\t{}\t({:.2f})".format(count_inward, 100 * float(count_inward) / pair_added))

        print("outward pairs\t{}\t({:.2f})".format(count_outward, 100 * float(count_outward) / pair_added))

        print("left pairs\t{}\t({:.2f})".format(count_left, 100 * float(count_left) / pair_added))

        print("right pairs\t{}\t({:.2f})".format(count_right, 100 * float(count_right) / pair_added))


class Tester(object):
    def __init__(self):
        import os
        hic_test_data_dir = os.environ.get('HIC_TEST_DATA_DIR', False)
        if hic_test_data_dir:
            self.root = hic_test_data_dir
        else:
            self.root = os.path.dirname(os.path.abspath(__file__)) + "/test/test_data/"
        self.bam_file_1 = self.root + "hic.bam"
