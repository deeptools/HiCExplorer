from __future__ import division

import argparse
import numpy as np
from scipy.sparse import coo_matrix, dia_matrix
import time
from os import unlink
import os
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
# import itertools

import pysam
from collections import OrderedDict

from six.moves import xrange
from future.utils import listitems
from copy import deepcopy
from ctypes import Structure, c_uint, c_ushort
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Array, RawArray

from intervaltree import IntervalTree, Interval

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

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


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description=('Using an alignment from a program that supports '
                     'local alignment (eg. Bowtie2) where both '
                     'PE reads are mapped using  the --local '
                     'option, this program reads such file and '
                     'creates a matrix of interactions.'
                     ))

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--samFiles', '-s',
                                help='The two PE alignment sam files to process',
                                metavar='two sam files',
                                nargs=2,
                                type=argparse.FileType('r'),
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='Output file name for the Hi-C matrix.',
                                metavar='FILENAME',
                                type=argparse.FileType('w'),
                                required=True)

    parserRequired.add_argument('--QCfolder',
                                help='Path of folder to save the quality control data for the matrix. The log files '
                                'produced this way can be loaded into `hicQC` in order to compare the quality of multiple '
                                'Hi-C libraries.',
                                metavar='FOLDER',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--outBam', '-b',
                           help='Output bam file to process. Optional parameter. '
                           'A bam file containing all valid Hi-C reads can be created '
                           'using this option. This bam file could be useful to inspect '
                           'the distribution of valid Hi-C reads pairs or for other '
                           'downstream analyses, but is not used by any HiCExplorer tool. '
                           'Computation will be significantly longer if this option is set.',
                           metavar='bam file',
                           type=argparse.FileType('w'),
                           required=False)

    group = parserOpt.add_mutually_exclusive_group(required=True)

    group.add_argument('--binSize', '-bs',
                       help='Size in bp for the bins. The bin size depends '
                            'on the depth of sequencing. Use a larger bin size for '
                            'libraries sequenced with lower depth. Alternatively, the location of '
                            'the restriction sites can be given (see --restrictionCutFile). '
                            'Optional for mcool file format: Define multiple resolutions which are all a multiple of the first value. '
                            ' Example: --binSize 10000 20000 50000 will create a mcool file formate containing the three defined resolutions.',
                       type=int,
                       nargs='+',
                       default=10000)

    group.add_argument('--restrictionCutFile', '-rs',
                       help=('BED file with all restriction cut places '
                             '(output of "findRestSite" command). '
                             'Should contain only  mappable '
                             'restriction sites. If given, the bins are '
                             'set to match the restriction fragments (i.e. '
                             'the region between one restriction site and '
                             'the next).'),
                       type=argparse.FileType('r'),
                       metavar='BED file')

    parserOpt.add_argument('--minDistance',
                           help='Minimum distance between restriction sites. '
                           'Restriction sites that are closer than this '
                           'distance are merged into one. This option only '
                           'applies if --restrictionCutFile is given.',
                           type=int,
                           default=300,
                           required=False)

    parserOpt.add_argument('--maxDistance',
                           help='This parameter is now obsolete. Use --maxLibraryInsertSize instead',
                           type=int)

    parserOpt.add_argument('--maxLibraryInsertSize',
                           help='The maximum library insert size defines different cut offs based on the maximum expected '
                           'library size. *This is not the average fragment size* but the higher end of the '
                           'the fragment size distribution (obtained using for example a Fragment Analyzer or a Bioanalyzer) '
                           'which usually is between 800 to 1500 bp. If this value if not known use the default of '
                           '1000. The insert value is used to decide if two mates belong to the same fragment (by '
                           'checking if they are within this max insert size) and to decide if a mate is too far '
                           'away from the nearest restriction site.',
                           type=int,
                           default=1000,
                           required=False)

    parserOpt.add_argument('--restrictionSequence', '-seq',
                           help='Sequence of the restriction site.')

    parserOpt.add_argument('--danglingSequence',
                           help='Sequence left by the restriction enzyme after cutting. Each restriction enzyme '
                                'recognizes a different DNA sequence and, after cutting, they leave behind a specific '
                                '"sticky" end or dangling end sequence.  For example, for HindIII the restriction site '
                                'is AAGCTT and the dangling end is AGCT. For DpnII, the restriction site and dangling '
                                'end sequence are the same: GATC. This information is easily found on the description '
                                'of the restriction enzyme. The dangling sequence is used to classify and report reads '
                                'whose 5\' end starts with such sequence as dangling-end reads. A significant portion '
                                'of dangling-end reads in a sample are indicative of a problem with the re-ligation '
                                'step of the protocol.')

    parserOpt.add_argument('--region', '-r',
                           help='Region of the genome to limit the operation to. '
                           'The format is chr:start-end. It is also possible to just '
                           'specify a chromosome, for example --region chr10',
                           metavar="CHR:START-END",
                           required=False,
                           type=genomicRegion
                           )
    # currently not implemented
    parserOpt.add_argument('--removeSelfLigation',
                           # help='If set, inward facing reads less than 1000 bp apart and having a restriction'
                           #     'site in between are removed. Although this reads do not contribute to '
                           #     'any distant contact, they are useful to account for bias in the data.',
                           help=argparse.SUPPRESS,
                           required=False,
                           default=True
                           # action='store_true'
                           )

    parserOpt.add_argument('--keepSelfCircles',
                           help='If set, outward facing reads without any restriction fragment (self circles) are kept. '
                           'They will be counted and shown in the QC plots.',
                           required=False,
                           action='store_true'
                           )

    parserOpt.add_argument('--minMappingQuality',
                           help='minimum mapping quality for reads to be accepted. '
                           'Because the restriction enzyme site could be located '
                           'on top of the read, this may reduce the '
                           'reported quality of the read. Thus, this parameter '
                           'may be adusted if too many low quality '
                           '(but otherwise perfectly valid Hi-C reads) are found. '
                           'A good strategy is to make a test run (using the --doTestRun), '
                           'then checking the results to see if too many low quality '
                           'reads are present and then using the bam file generated to '
                           'check if those low quality reads are caused by the read '
                           'not being mapped entirely.',
                           required=False,
                           default=15,
                           type=int
                           )
    parserOpt.add_argument('--threads',
                           help='Number of threads. Using the python multiprocessing module. '
                           'One master process which is used to read the input file into the buffer and one process which is merging '
                           'the output bam files of the processes into one output bam file. All other threads do the actual computation. '
                           'Minimum value for the \'--thread\' parameter is 2. '
                           'The usage of 8 threads is optimal if you have an HDD. A higher number of threads is only '
                           'useful if you have a fast SSD. Have in mind that the performance of hicBuildMatrix is influenced by '
                           'the number of threads, the speed of your hard drive and the inputBufferSize. To clearify: the peformance '
                           'with a higher thread number is not negative influenced but not positiv too. With a slow HDD and a high number of '
                           'threads many threads will do nothing most of the time. ',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--inputBufferSize',
                           help='Size of the input buffer of each thread. 400,000 read pairs per input file per thread is the default value. '
                           'Reduce this value to decrease memory usage.',
                           required=False,
                           default=400000,
                           type=int
                           )
    parserOpt.add_argument('--doTestRun',
                           help='A test run is useful to test the quality '
                           'of a Hi-C experiment quickly. It works by '
                           'testing only 1,000,000 reads. This option '
                           'is useful to get an idea of quality control '
                           'values like inter-chromosomal interactions, '
                           'duplication rates etc.',
                           action='store_true'
                           )

    parserOpt.add_argument('--skipDuplicationCheck',
                           help='Identification of duplicated read pairs is '
                           'memory consuming. Thus, in case of memory '
                           'errors this check can be skipped. However, '
                           'consider running a `--doTestRun` first to '
                           'get an estimation of the duplicated reads. ',
                           action='store_true'
                           )

    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


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
    >>> _file.close()
    >>> file_tmp = open(_file.name, 'w')
    >>> foo = file_tmp.write('chr1\t10\t20\tH1\t0\n')
    >>> foo = file_tmp.write("chr1\t60\t70\tH2\t0\n")
    >>> file_tmp.close()
    >>> bed2interval_list(open(_file.name, 'r'))
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
            log.error("error reading BED file at line {}".format(count))

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
    return listitems(list_chrom_sizes)


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
    for idx in xrange(len(bin_intervals) - 1):
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
                 pKeepSelfCircles, pRestrictionSequence, pRemoveSelfLigation, pMatrixSize,
                 pRfPositions, pRefId2name,
                 pDanglingSequences, pBinsize, pResultIndex,
                 pQueueOut, pTemplate, pOutputBamSet, pCounter,
                 pSharedBinIntvalTree, pDictBinIntervalTreeIndex, pCoverage, pCoverageIndex,
                 pOutputFileBufferDir, pRow, pCol, pData,
                 pMaxInsertSize):
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
    pRestrictionSequence : String, the restriction sequence
    pRemoveSelfLigation : If self ligations should be removed
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

    one_mate_unmapped = 0
    one_mate_low_quality = 0
    one_mate_not_unique = 0
    dangling_end = 0
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

                    if len(has_rf) == 0:
                        self_circle += 1
                        if not pKeepSelfCircles:
                            continue

            if abs(mate2.pos - mate1.pos) < pMaxInsertSize and orientation == 'inward':
                # check for dangling ends if the restriction sequence is known and if they look
                # like 'same fragment'
                if pRestrictionSequence:
                    if pDanglingSequences:
                        if check_dangling_end(mate1, pDanglingSequences) or \
                                check_dangling_end(mate2, pDanglingSequences):
                            dangling_end += 1
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
                    frag_start = min(mate1.pos, mate2.pos) + \
                        len(pRestrictionSequence)
                    frag_end = max(mate1.pos + mate1.qlen, mate2.pos +
                                   mate2.qlen) - len(pRestrictionSequence)
                    mate_ref = pRefId2name[mate1.rname]
                    has_rf = sorted(
                        pRfPositions[mate_ref][frag_start: frag_end])

                # case when there is no restriction fragment site between the
                # mates
                if len(has_rf) == 0:
                    same_fragment += 1
                    continue

                self_ligation += 1

                if pRemoveSelfLigation:
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
            length_coverage = pCoverageIndex[mate_bin_id].end - pCoverageIndex[mate_bin_id].begin
            vec_end = min(length_coverage, int(vec_start +
                                               len(mate.seq) / pBinsize))
            coverage_index = pCoverageIndex[mate_bin_id].begin + vec_start
            coverage_end = pCoverageIndex[mate_bin_id].begin + vec_end
            for i in xrange(coverage_index, coverage_end, 1):
                pCoverage[i] += 1

        pRow[pair_added] = mate_bins[0]
        pCol[pair_added] = mate_bins[1]
        pData[pair_added] = np.uint8(1)

        pair_added += 1
        if pOutputBamSet:

            out_bam_index_buffer.append(iter_num - 1)

    pQueueOut.put([[one_mate_unmapped, one_mate_low_quality, one_mate_not_unique, dangling_end, self_circle, self_ligation, same_fragment,
                    mate_not_close_to_rf, count_inward, count_outward,
                    count_left, count_right, inter_chromosomal, short_range, long_range, pair_added, len(pMateBuffer1), pResultIndex, pCounter, out_bam_index_buffer]])
    return


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
    is also constructed.
    """

    args = parse_arguments().parse_args(args)

    # for backwards compatibility
    if args.maxDistance is not None:
        args.maxLibraryInsertSize = args.maxDistance
    try:
        QC.make_sure_path_exists(args.QCfolder)
    except OSError:
        exit("Can't open/create QC folder path: {}. Please check".format(args.QCfolder))

    if args.threads < 2:
        args.threads = 2
        warnings.warn("\nAt least two threads need to be defined. Setting --threads = 2!s\n")

    if args.danglingSequence and not args.restrictionSequence:
        exit("\nIf --danglingSequence is set, --restrictonSequence needs to be set too.\n")

    log.info("reading {} and {} to build hic_matrix\n".format(args.samFiles[0].name,
                                                              args.samFiles[1].name))
    str1 = pysam.Samfile(args.samFiles[0].name, 'rb')
    str2 = pysam.Samfile(args.samFiles[1].name, 'rb')

    args.samFiles[0].close()
    args.samFiles[1].close()

    if args.outBam:
        args.outBam.close()
        out_bam_file = pysam.Samfile(args.outBam.name, 'wb', template=str1)

    chrom_sizes = get_chrom_sizes(str1)

    read_pos_matrix = ReadPositionMatrix()

    # define bins
    rf_positions = None
    if args.restrictionCutFile:
        rf_interval = bed2interval_list(args.restrictionCutFile)
        bin_intervals = get_rf_bins(rf_interval,
                                    min_distance=args.minDistance,
                                    max_distance=args.maxLibraryInsertSize)

        rf_positions = intervalListToIntervalTree(rf_interval)
    else:
        bin_intervals = get_bins(args.binSize[0], chrom_sizes, args.region)

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
    dangling_sequences = dict()
    if args.danglingSequence:
        # build a list of dangling sequences
        args.restrictionSequence = args.restrictionSequence.upper()
        args.danglingSequence = args.danglingSequence.upper()
        dangling_sequences['pat_forw'] = args.danglingSequence
        dangling_sequences['pat_rev'] = str(
            Seq(args.danglingSequence, generic_dna).reverse_complement())

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
    args.threads = args.threads - 1
    row = [None] * args.threads
    col = [None] * args.threads
    data = [None] * args.threads
    for i in xrange(args.threads):
        row[i] = RawArray(c_uint, args.inputBufferSize)
        col[i] = RawArray(c_uint, args.inputBufferSize)
        data[i] = RawArray(c_ushort, args.inputBufferSize)

    start_time = time.time()

    iter_num = 0
    # pair_added = 0
    # hic_matrix = None

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
    buffer_workers1 = [None] * args.threads
    buffer_workers2 = [None] * args.threads

    # output buffer to write bam with mate1 and mate2 pairs
    process = [None] * args.threads
    all_data_processed = False
    hic_matrix = coo_matrix((matrix_size, matrix_size), dtype='uint32')
    queue = [None] * args.threads

    all_threads_done = False
    thread_done = [False] * args.threads
    count_output = 0
    count_call_of_read_input = 0
    computed_pairs = 0

    while not all_data_processed or not all_threads_done:

        for i in xrange(args.threads):
            if queue[i] is None and not all_data_processed:
                count_call_of_read_input += 1

                buffer_workers1[i], buffer_workers2[i], all_data_processed, \
                    duplicated_pairs_, one_mate_unmapped_, one_mate_not_unique_, \
                    one_mate_low_quality_, iter_num_ = readBamFiles(pFileOneIterator=str1,
                                                                    pFileTwoIterator=str2,
                                                                    pNumberOfItemsPerBuffer=args.inputBufferSize,
                                                                    pSkipDuplicationCheck=args.skipDuplicationCheck,
                                                                    pReadPosMatrix=read_pos_matrix,
                                                                    pRefId2name=ref_id2name,
                                                                    pMinMappingQuality=args.minMappingQuality
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
                    pMinMappingQuality=args.minMappingQuality,
                    pKeepSelfCircles=args.keepSelfCircles,
                    pRestrictionSequence=args.restrictionSequence,
                    pRemoveSelfLigation=args.removeSelfLigation,
                    pMatrixSize=matrix_size,
                    pRfPositions=rf_positions,
                    pRefId2name=ref_id2name,
                    pDanglingSequences=dangling_sequences,
                    pBinsize=binsize,
                    pResultIndex=i,
                    pQueueOut=queue[i],
                    pTemplate=str1,
                    pOutputBamSet=args.outBam,
                    pCounter=count_output,
                    pSharedBinIntvalTree=shared_build_intval_tree,
                    pDictBinIntervalTreeIndex=index_dict,
                    pCoverage=coverage,
                    pCoverageIndex=pos_coverage,
                    pOutputFileBufferDir="",
                    pRow=row[i],
                    pCol=col[i],
                    pData=data[i],
                    pMaxInsertSize=args.maxLibraryInsertSize
                ))
                process[i].start()
                count_output += 1

            elif queue[i] is not None and not queue[i].empty():
                result = queue[i].get()

                if result[0] is not None:
                    elements = result[0][15]
                    if hic_matrix is None:
                        hic_matrix = coo_matrix(
                            (data[i][:elements], (row[i][:elements], col[i][:elements])), shape=(matrix_size, matrix_size))
                    else:
                        hic_matrix += coo_matrix(
                            (data[i][:elements], (row[i][:elements], col[i][:elements])), shape=(matrix_size, matrix_size))

                    dangling_end += result[0][3]
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
                if args.doTestRun and iter_num > 1e5:
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

    # the resulting matrix is only filled unevenly with some pairs
    # int the upper triangle and others in the lower triangle. To construct
    # the definite matrix I add the values from the upper and lower triangles
    # and subtract the diagonal to avoid double counting it.
    # The resulting matrix is symmetric.
    if args.outBam:
        out_bam_file.close()

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

    chr_name_list, start_list, end_list = list(zip(*bin_intervals))
    bin_intervals = list(zip(chr_name_list, start_list, end_list, bin_max))
    hic_ma = hm.hiCMatrix()
    hic_ma.setMatrix(hic_matrix, cut_intervals=bin_intervals)

    args.outFileName.close()
    # removing the empty file. Otherwise the save method
    # will say that the file already exists.
    unlink(args.outFileName.name)

    if args.outFileName.name.endswith('.cool') and args.binSize is not None and len(args.binSize) > 2:

        matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool')
        matrixFileHandlerOutput.set_matrix_variables(hic_ma.matrix,
                                                     hic_ma.cut_intervals,
                                                     hic_ma.nan_bins,
                                                     hic_ma.correction_factors,
                                                     hic_ma.distance_counts)
        matrixFileHandlerOutput.save(args.outFileName.name + '::/resolutions/' + str(args.binSize[0]), pSymmetric=True, pApplyCorrection=False)

        for resolution in args.binSize[1:]:
            hic_matrix_to_merge = deepcopy(hic_ma)
            _mergeFactor = int(resolution) // args.binSize[0]
            merged_matrix = hicMergeMatrixBins.merge_bins(hic_matrix_to_merge, _mergeFactor)
            matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool', pAppend=True)
            matrixFileHandlerOutput.set_matrix_variables(merged_matrix.matrix,
                                                         merged_matrix.cut_intervals,
                                                         merged_matrix.nan_bins,
                                                         merged_matrix.correction_factors,
                                                         merged_matrix.distance_counts)
            matrixFileHandlerOutput.save(args.outFileName.name + '::/resolutions/' + str(resolution), pSymmetric=True, pApplyCorrection=False)

    else:
        hic_ma.save(args.outFileName.name)

    """
    if args.restrictionCutFile:
        # load the matrix to mask those
        # bins that most likely didn't
        # have a restriction site that was cutted

        # reload the matrix as a Hi-C matrix object
        hic_matrix = hm.hiCMatrix(args.outFileName.name)

        hic_matrix.maskBins(get_poor_bins(bin_max))
        hic_matrix.save(args.outFileName.name)
    """
    if args.removeSelfLigation:
        msg = " (removed)"
    else:
        msg = " (not removed)"

    mappable_unique_high_quality_pairs = iter_num - (one_mate_unmapped + one_mate_low_quality + one_mate_not_unique)

    log_file_name = os.path.join(args.QCfolder, "QC.log")
    log_file = open(log_file_name, 'w')
    log_file.write("""
File\t{}\t\t
Pairs considered\t{}\t\t
Min rest. site distance\t{}\t\t
Max library insert size\t{}\t\t

""".format(args.outFileName.name, iter_num, args.minDistance, args.maxLibraryInsertSize))

    log_file.write("#\tcount\t(percentage w.r.t. total sequenced reads)\n")

    log_file.write("Pairs mappable, unique and high quality\t{}\t({:.2f})\n".
                   format(mappable_unique_high_quality_pairs,
                          100 * float(mappable_unique_high_quality_pairs) / iter_num))

    log_file.write("Pairs used\t{}\t({:.2f})\n".
                   format(pair_added, 100 * float(pair_added) / iter_num))

    log_file.write("One mate unmapped\t{}\t({:.2f})\n".
                   format(one_mate_unmapped, 100 * float(one_mate_unmapped) / iter_num))

    log_file.write("One mate not unique\t{}\t({:.2f})\n".
                   format(one_mate_not_unique, 100 * float(one_mate_not_unique) / iter_num))

    log_file.write("One mate low quality\t{}\t({:.2f})\n".
                   format(one_mate_low_quality, 100 * float(one_mate_low_quality) / iter_num))

    log_file.write("\n#\tcount\t(percentage w.r.t. mappable, unique and high quality pairs)\n")

    log_file.write("dangling end\t{}\t({:.2f})\n".
                   format(dangling_end, 100 * float(dangling_end) / mappable_unique_high_quality_pairs))

    log_file.write("self ligation{}\t{}\t({:.2f})\n".
                   format(msg, self_ligation, 100 * float(self_ligation) / mappable_unique_high_quality_pairs))

    log_file.write("One mate not close to rest site\t{}\t({:.2f})\n".
                   format(mate_not_close_to_rf, 100 * float(mate_not_close_to_rf) / mappable_unique_high_quality_pairs))

    log_file.write("same fragment\t{}\t({:.2f})\n".
                   format(same_fragment, 100 * float(same_fragment) / mappable_unique_high_quality_pairs))

    log_file.write("self circle\t{}\t({:.2f})\n".
                   format(self_circle, 100 * float(self_circle) / mappable_unique_high_quality_pairs))

    log_file.write("duplicated pairs\t{}\t({:.2f})\n".
                   format(duplicated_pairs, 100 * float(duplicated_pairs) / mappable_unique_high_quality_pairs))

    if pair_added > 0:
        log_file.write("\n#\tcount\t(percentage w.r.t. total valid pairs used)\n")
        log_file.write("inter chromosomal\t{}\t({:.2f})\n".
                       format(inter_chromosomal, 100 * float(inter_chromosomal) / pair_added))

        log_file.write("short range < 20kb\t{}\t({:.2f})\n".
                       format(short_range, 100 * float(short_range) / pair_added))

        log_file.write("long range\t{}\t({:.2f})\n".
                       format(long_range, 100 * float(long_range) / pair_added))

        log_file.write("inward pairs\t{}\t({:.2f})\n".
                       format(count_inward, 100 * float(count_inward) / pair_added))

        log_file.write("outward pairs\t{}\t({:.2f})\n".
                       format(count_outward, 100 * float(count_outward) / pair_added))

        log_file.write("left pairs\t{}\t({:.2f})\n".
                       format(count_left, 100 * float(count_left) / pair_added))

        log_file.write("right pairs\t{}\t({:.2f})\n".
                       format(count_right, 100 * float(count_right) / pair_added))

    log_file.close()
    QC.main("-l {} -o {}".format(log_file_name, args.QCfolder).split())


class Tester(object):
    def __init__(self):
        hic_test_data_dir = os.environ.get('HIC_TEST_DATA_DIR', False)
        if hic_test_data_dir:
            self.root = hic_test_data_dir
        else:
            self.root = os.path.dirname(
                os.path.abspath(__file__)) + "/test/test_data/"
        self.bam_file_1 = os.path.join(self.root, "hic.bam")
