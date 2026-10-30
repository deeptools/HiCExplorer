"""hicexplorer.readBed, implemented in C++.

ReadBed keeps the attributes, the methods and the behaviour it had in
HiCExplorer 3.7.6; every method body calls hicexplorer._hicx.readBed with
self (cpp/python/src/api/readBed.cpp, cpp/core/include/hicx/read_bed.hpp).
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import readBed as _impl


class ReadBed(object):
    """
    Reads a bed file. Based on the number of fields
    it tries to guess the type of bed file used. Current options
    are bed3, bed6 and bed12

    Example:
    bed = readBed(open("file.bed", 'r'))
    for interval in bed:
        print bed['start']

    """

    def __init__(self, file_handle):
        """
        :param file_handle: file handle
        :return:
        """
        _impl.ReadBed___init__(self, file_handle)

    def __iter__(self):
        return self

    def get_no_comment_line(self):
        """
        Skips comment lines starting with '#'
        "track" or "browser" in the bed files
        :return:
        """
        return _impl.ReadBed_get_no_comment_line(self)

    def guess_file_type(self, line_values):
        """try to guess type of bed file by counting the fields
        """
        return _impl.ReadBed_guess_file_type(self, line_values)

    def next(self):
        """
        :return: bedInterval object
        """
        return _impl.ReadBed_next(self)

    def __next__(self):
        """
        :return: bedInterval object
        """
        return _impl.ReadBed_next(self)

    def get_bed_interval(self, bed_line):
        r"""
        Processes each bed line from a bed file, casts the values and returns
        a namedtuple object

        >>> bed_line="chr1\t0\t1000\tgene_1\t0.5\t-\t0\t1000\t0\t3\t10,20,100\t20,200,700"
        >>> with open('/tmp/test.bed', 'w') as fh:
        ...     foo = fh.write(bed_line)
        >>> bed_f = ReadBed(open('/tmp/test.bed','r'))
        >>> bed = bed_f.get_bed_interval(bed_line)
        >>> bed.chromosome
        'chr1'
        >>> bed.block_starts
        [20, 200, 700]

        >>> bed_line="chr2\t0\t1000\tgene_1\t0.5\t-\n"
        >>> with open('/tmp/test.bed', 'w') as fh:
        ...     foo = fh.write(bed_line)
        >>> bed_f = ReadBed(open('/tmp/test.bed','r'))
        >>> bed_f.get_bed_interval(bed_line)
        BedInterval(chromosome='chr2', start=0, end=1000, name='gene_1', score=0.5, strand='-')
        """
        return _impl.ReadBed_get_bed_interval(self, bed_line)
