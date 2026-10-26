"""Equivalence cases for hicexplorer.readBed.ReadBed (apieq.py).

ReadBed is a class reading from an open file, so every case is a driver that
opens the case's file, constructs the reader and reads it; what the reader
returns and the state it keeps are compared.
"""
import os

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
T = 'hicexplorer.readBed.ReadBed'


def state(bed):
    return {'file_type': bed.file_type, 'line_number': bed.line_number, 'fields': bed.fields,
            'interval_fields': bed.BedInterval._fields, 'prev_chrom': bed.prev_chrom,
            'prev_start': bed.prev_start, 'prev_line': bed.prev_line}


def iterate(ReadBed, path, mode='r'):
    """All intervals, then the reader's state."""
    with open(path, mode) as handle:
        bed = ReadBed(handle)
        after_init = state(bed)
        intervals = []
        error = None
        try:
            for interval in bed:
                intervals.append(interval)
        except Exception as e:  # compared as part of the result
            error = (type(e).__name__, str(e))
        return {'init': after_init, 'intervals': intervals, 'error': error, 'final': state(bed)}


def next_method(ReadBed, path):
    """The Python 2 style next() method, three times."""
    with open(path) as handle:
        bed = ReadBed(handle)
        return [bed.next() for _ in range(3)], state(bed)


def interval_of(ReadBed, path, line):
    """get_bed_interval on a line the file did not provide."""
    with open(path) as handle:
        bed = ReadBed(handle)
        return bed.get_bed_interval(line), state(bed)


def subclass(ReadBed, path):
    """A subclass that overrides a method sees it called by the others."""
    class Upper(ReadBed):
        def get_bed_interval(self, bed_line):
            return ReadBed.get_bed_interval(self, bed_line.replace('chr', 'CHR'))
    with open(path) as handle:
        return list(Upper(handle))


BED6 = b"chr1\t10\t20\tg1\t0.5\t+\nchr1\t15\t30\tg2\t1\t-\nchr2\t5\t9\tg3\tx\t.\n"
BED12 = (b"chr1\t0\t1000\tgene_1\t0.5\t-\t0\t1000\t0\t3\t10,20,100\t20,200,700\n"
         b"chr1\t5\t900\tgene_2\t1\t+\t10\t800\t255,0,0\t2\t10,20,\t0,100,\n"
         b"chr1\t6\t901\tgene_3\t1\t+\t10\t800\t1,x,2\t2\t10,y\t0,100\n")
BED9 = b"chr1\t0\t10\tn\t1\t+\t0\t10\t0,0,255\nchr1\t2\t12\tn\t1\t1\t0\t10\t1,2\n"


CASES = [
    Case('real_bed6', T, args=(os.path.join(TEST_DATA, 'hicAggregateContacts', 'test_regions_strand.bed'),),
         driver='iterate'),
    Case('real_bed9', T, args=(os.path.join(TEST_DATA, 'domains.bed'),), driver='iterate'),
    Case('real_bed3', T, args=(os.path.join(TEST_DATA, 'regions.bed'),), driver='iterate'),
    Case('real_reference_points', T, args=(os.path.join(TEST_DATA, 'cHi-C', 'referencePoints.bed'),),
         driver='iterate'),
    Case('bed6', T, args=(File('a.bed'),), files={'a.bed': BED6}, driver='iterate'),
    Case('bed12', T, args=(File('a.bed'),), files={'a.bed': BED12}, driver='iterate'),
    Case('bed9', T, args=(File('a.bed'),), files={'a.bed': BED9}, driver='iterate'),
    Case('bed3_comments', T, args=(File('a.bed'),), driver='iterate',
         files={'a.bed': b"#x\ntrack name=a\nbrowser p\n   \nchr1\t1\t2\n\n# y\nchr1\t5\t6\n"}),
    Case('bed3_trailing_comments', T, args=(File('a.bed'),), driver='iterate',
         files={'a.bed': b"chr1\t1\t2\n# end\n#\n"}),
    Case('bedgraph', T, args=(File('a.bed'),), files={'a.bed': b"chr1\t1\t2\t0.5\n"}, driver='iterate'),
    Case('five_fields', T, args=(File('a.bed'),), files={'a.bed': b"chr1\t1\t2\tn\t3\n"}, driver='iterate'),
    Case('seven_fields', T, args=(File('a.bed'),), files={'a.bed': b"chr1\t1\t2\tn\t3\t+\t9\n"}, driver='iterate'),
    Case('unsorted', T, args=(File('a.bed'),), files={'a.bed': b"chr1\t10\t20\nchr1\t5\t6\n"}, driver='iterate'),
    Case('start_not_before_end', T, args=(File('a.bed'),), files={'a.bed': b"chr1\t10\t10\n"}, driver='iterate'),
    Case('not_an_integer', T, args=(File('a.bed'),), files={'a.bed': b"chr1\t1x\t20\n"}, driver='iterate'),
    Case('invalid_strand', T, args=(File('a.bed'),), driver='iterate',
         files={'a.bed': b"chr1\t1\t2\tn\t0\tw\nchr1\t3\t4\tn\t0\t-1\n"}),
    Case('python_int_syntax', T, args=(File('a.bed'),), files={'a.bed': b"chr1\t+1_0\t 20\n"}, driver='iterate'),
    Case('float_scores', T, args=(File('a.bed'),), driver='iterate',
         files={'a.bed': b"chr1\t1\t2\tn\tinf\t+\nchr1\t1\t2\tn\t1e-3\t+\nchr1\t1\t2\tn\tNaN\t+\n"}),
    Case('mixed_field_counts', T, args=(File('a.bed'),), files={'a.bed': b"chr1\t1\t2\tn\t0\t+\nchr1\t3\t4\n"},
         driver='iterate'),
    Case('binary_handle', T, args=(File('a.bed'), 'rb'), files={'a.bed': BED6}, driver='iterate'),
    Case('unicode', T, args=(File('a.bed'),), files={'a.bed': "chré\t1\t2\nchré\t3\t4　\n".encode()},
         driver='iterate'),
    Case('crlf', T, args=(File('a.bed'),), files={'a.bed': b"chr1\t1\t2\r\nchr1\t3\t4\r\n"}, driver='iterate'),
    Case('empty_file', T, args=(File('a.bed'),), files={'a.bed': b""}, driver='iterate'),
    Case('only_comments', T, args=(File('a.bed'),), files={'a.bed': b"# a\n# b\n"}, driver='iterate'),
    Case('next_method', T, args=(File('a.bed'),), files={'a.bed': BED6}, driver='next_method'),
    Case('get_bed_interval_doctest12', T, driver='interval_of', files={'a.bed': BED12},
         args=(File('a.bed'), "chr1\t0\t1000\tgene_1\t0.5\t-\t0\t1000\t0\t3\t10,20,100\t20,200,700")),
    Case('get_bed_interval_doctest6', T, driver='interval_of', files={'a.bed': BED6},
         args=(File('a.bed'), "chr2\t0\t1000\tgene_1\t0.5\t-\n")),
    Case('get_bed_interval_bytes', T, driver='interval_of', files={'a.bed': BED6},
         args=(File('a.bed'), b"chr2\t0\t1000\tgene_1\t0.5\t-\n")),
    Case('get_bed_interval_wrong_count', T, driver='interval_of', files={'a.bed': BED6},
         args=(File('a.bed'), "chr2\t0\t1000\n")),
    Case('subclass_override', T, args=(File('a.bed'),), files={'a.bed': BED6}, driver='subclass'),
]
