"""Equivalence cases for hicexplorer.lib.buildMatrixMethods (apieq.py).

The inputs are the repository's real alignments (hicBuildMatrix/R1_1000.bam
and R2_1000.bam, hic.bam, small_test_R1/R2_unsorted.bam), restriction site
files (DpnII.bed, hicFindRestSite/hindIII.bed) and chromosome sizes, read with
pysam in the drivers exactly as 3.7.6 code reads them, plus small SAM files
for supplementary alignments and stand-in objects for the edge cases.
pysam reads are compared by their SAM text (AlignedSegment.to_string()).

createMatrix has its own cases file, cases_lib_buildMatrixMethods_createMatrix.py.

Deviations from 3.7.6 (cpp/python/src/api/lib_buildMatrixMethods.cpp):
  * types 3.7.6 accepted only by accident raise TypeError naming HiCExplorer
    4: coordinates that are neither int nor numpy.int64 (get_bins,
    get_rf_bins, enlarge_bins), non-integer CIGAR tuples, flags, mapping
    qualities, insert sizes and bin sizes, reference names and sequences
    that are not str;
  * process_data's failure message carries the Python traceback for Python
    errors and a one-frame traceback for errors raised in C++; the message
    before "Traceback" is identical;
  * get_rf_bins returns, for each bin, the chromosome object of the first
    site of that chromosome rather than of the site that opened the bin
    (they differ only for equal chromosome objects of different types);
  * check_dangling_end upper-cases non-ASCII sequences with str.upper().
"""
import os

import numpy as np

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'test_data')
BM = os.path.join(DATA, 'hicBuildMatrix')
R1 = os.path.join(BM, 'R1_1000.bam')
R2 = os.path.join(BM, 'R2_1000.bam')
HIC_BAM = os.path.join(DATA, 'hic.bam')
DPNII = os.path.join(BM, 'DpnII.bed')
HINDIII = os.path.join(DATA, 'hicFindRestSite', 'hindIII.bed')
HINDIII_SMALL = os.path.join(DATA, 'hicFindRestSite', 'hindIII_DpnII.bed')
SIZES = os.path.join(BM, 'dm3.chrom.sizes')
SMALL_R1 = os.path.join(DATA, 'small_test_R1_unsorted.bam')
SMALL_R2 = os.path.join(DATA, 'small_test_R2_unsorted.bam')
M = 'hicexplorer.lib.buildMatrixMethods.'


def module_of(target):
    import sys
    name = target.__module__ if not isinstance(target, type(os)) else target.__name__
    return sys.modules[name]


def sam(read):
    return None if read is None else read.to_string()


def sams(reads):
    return None if reads is None else [sam(r) for r in reads]


# ---------------------------------------------------------------------------
# small inputs

SA_SAM = (
    "@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:1575\n"
    "read\t65\tchr1\t33\t20\t10S1D25M\t=\t200\t167\tAGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG\t"
    "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<\tSA:Z:chr1,88,+,1S34M,30,0;chr1,120,+,5S30M,30,0;\n"
    "read\t2113\tchr1\t88\t30\t1S34M\t=\t500\t412\tACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA\t"
    "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<\n"
    "read\t2113\tchr1\t120\t30\t5S30M\t=\t500\t412\tACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA\t"
    "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<\n"
    "other\t65\tchr1\t300\t30\t35M\t=\t200\t167\tAGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG\t"
    "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<\tSA:Z:chr1,400,+,35M,30,0;\n"
    "stranger\t2113\tchr1\t400\t30\t35M\t=\t200\t167\tAGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG\t"
    "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<\n"
    "last\t0\tchr1\t500\t30\t35M\t*\t0\t0\tAGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG\t"
    "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<\n").encode()

REVERSE_SAM = (
    "@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:1575\n"
    "r\t81\tchr1\t33\t20\t20M15S\t=\t200\t167\tAGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG\t*\n"
    "r\t2129\tchr1\t88\t30\t15M2I18S\t=\t500\t412\tACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA\t*\n"
    "r\t2129\tchr1\t120\t30\t30M5S\t=\t500\t412\tACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA\t*\n"
    "q\t2129\tchr1\t120\t30\t30M5S\t=\t500\t412\tACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA\t*\n").encode()


class FakeRead(object):
    def __init__(self, **fields):
        self.__dict__.update(fields)

    def has_tag(self, tag):
        return tag in getattr(self, 'tags', {})

    def get_tag(self, tag):
        return self.tags[tag]


class FakeBam(object):
    def __init__(self, references, lengths):
        self.references = references
        self.lengths = lengths


def open_reads(path):
    import pysam
    return pysam.AlignmentFile(path, 'r' if path.endswith('.sam') else 'rb')


# ---------------------------------------------------------------------------
# drivers

def struct_layout(cls):
    import ctypes
    from multiprocessing.sharedctypes import RawArray
    value = cls(*range(1, len(cls._fields_) + 1))
    array = RawArray(cls, [tuple(range(10 * i, 10 * i + len(cls._fields_))) for i in range(3)])
    return ([(name, kind.__name__) for name, kind in cls._fields_], ctypes.sizeof(cls),
            issubclass(cls, ctypes.Structure), [getattr(value, name) for name, _ in cls._fields_],
            [[getattr(item, name) for name, _ in cls._fields_] for item in array], cls.__name__, cls.__doc__)


def position_matrix(ReadPositionMatrix, calls):
    rp = ReadPositionMatrix()
    results = []
    for call in calls:
        try:
            results.append(rp.is_duplicated(*call))
        except Exception as error:
            results.append((type(error).__name__, str(error)))
    return results, sorted(rp.pos_matrix)


def trees(target, *args):
    result = target(*args)
    return {key: sorted((iv.begin, iv.end, iv.data) for iv in tree) for key, tree in result.items()}, \
        [type(tree).__name__ for tree in result.values()]


def bins_of_bam(target, bin_size, path, *args):
    module = module_of(target)
    return target(bin_size, module.get_chrom_sizes(open_reads(path)), *args)


def bed_list(target, path, sizes, region, mode='r'):
    with open(path, mode) as handle:
        result = target(handle, sizes, region)
        return result, handle.tell()


def bed_from_list(target, lines, sizes, region):
    return target(iter(lines), sizes, region)


def rf_bins_of_bed(target, path, *args):
    module = module_of(target)
    with open(path) as handle:
        sites = module.bed2interval_list(handle, None, None)
    return target(sites, *args)


def chrom_sizes_of(target, path):
    return target(open_reads(path))


def chrom_sizes_fake(target, references, lengths):
    return target(FakeBam(references, lengths))


def dangling_reads(target, path, ds, count=200):
    reads = open_reads(path)
    out = []
    for _ in range(count):
        read = next(reads)
        if read.seq is None:
            continue
        out.append(target(read, ds))
    return out


def dangling_fake(target, ds, **fields):
    return target(FakeRead(**fields), ds)


def supplementary(target, path, skip=0):
    reads = open_reads(path)
    for _ in range(skip):
        next(reads)
    read = next(reads)
    result = target(read, reads)
    try:
        following = sam(next(reads))
    except StopIteration:
        following = None
    return sams(result), following


def supplementary_fake(target, tags, names):
    read = FakeRead(qname='a', tags=tags)
    others = iter([FakeRead(qname=n) for n in names])
    result = target(read, others)
    return None if result is None else [r.qname for r in result], list(r.qname for r in others)


def correct_map(target, path, as_tuple=False):
    reads = list(open_reads(path))
    primary = reads[0]
    supplements = [r for r in reads[1:] if r.qname == primary.qname]
    if as_tuple:
        supplements = tuple(supplements)
    result = target(primary, supplements)
    return sam(result), [i for i, r in enumerate([primary] + list(supplements)) if r is result]


def correct_map_mixed(target, path):
    reads = list(open_reads(path))
    return sam(target(reads[0], reads[1:]))


def correct_map_fake(target, primary, others):
    reads = [FakeRead(qname='x', is_reverse=rev, cigartuples=cig) for rev, cig in [primary] + others]
    result = target(reads[0], reads[1:])
    return reads.index(result)


def call_and_list_state(target, intervals, sizes):
    result = target(intervals, sizes)
    return result, result is intervals


def enlarge_partial(target, intervals, sizes):
    try:
        target(intervals, sizes)
    except KeyError as error:
        return 'KeyError', str(error), intervals


def enlarge_rf(target, path, sizes_path):
    module = module_of(target)
    with open(path) as handle:
        sites = module.bed2interval_list(handle, None, None)
    bins = module.get_rf_bins(sites, 150, 1500)
    sizes = module.get_chrom_sizes(open_reads(sizes_path))
    return target(bins, sizes)


def read_buffers(target, path1, path2, size, skip, minimum, rounds=50):
    """readBamFiles until it reports that all data was read."""
    module = module_of(target)
    file1 = open_reads(path1)
    file2 = open_reads(path2)
    positions = module.ReadPositionMatrix()
    results = []
    for _ in range(rounds):
        result = target(file1, file2, size, skip, positions, file1.references, minimum)
        results.append((sams(result[0]), sams(result[1])) + tuple(result[2:]))
        if result[2]:
            break
    return results, len(positions.pos_matrix)


def read_buffers_error(target, path1, path2):
    module = module_of(target)
    file1 = open_reads(path1)
    file2 = open_reads(path2)
    return target(file1, file2, 100, False, module.ReadPositionMatrix(), file1.references, 15)


# process_data: the inputs createMatrix builds for it

def process_inputs(module, path1, path2, bin_size, rf_paths, restriction, dangling, buffer_size):
    from multiprocessing.sharedctypes import Array, RawArray
    from ctypes import c_uint, c_ushort
    file1 = open_reads(path1)
    file2 = open_reads(path2)
    chrom_sizes = module.get_chrom_sizes(file1)
    rf_interval = []
    for path in rf_paths:
        with open(path) as handle:
            rf_interval.extend(module.bed2interval_list(handle, chrom_sizes, None))
    rf_positions = module.intervalListToIntervalTree(rf_interval) if rf_paths else None
    if bin_size:
        bin_intervals = module.get_bins(bin_size, chrom_sizes)
    else:
        bin_intervals = module.get_rf_bins(rf_interval, min_distance=150, max_distance=1500)
    tree = module.intervalListToIntervalTree(bin_intervals)
    shared, index, end = [], {}, -1
    for seq in tree:
        start = end + 1
        items = sorted((iv.begin, iv.end, iv.data) for iv in tree[seq])
        end = start + len(tree[seq]) - 1
        index[seq] = (start, end)
        shared.extend(items)
    shared = RawArray(module.C_Interval, shared)
    starts, ends, n = [], [], 0
    for _, s, e in bin_intervals:
        starts.append(n)
        n += (e - s) // 10
        ends.append(n - 1)
    coverage_index = RawArray(module.C_Coverage, list(zip(starts, ends)))
    coverage = Array(c_uint, [0] * n)
    dangling_sequences = {}
    if dangling:
        for r, d in zip(restriction, dangling):
            dangling_sequences[r] = {'pat_forw': d, 'pat_rev': d[::-1].translate(str.maketrans('ACGT', 'TGCA'))}
    buffers = module.readBamFiles(file1, file2, buffer_size, False, module.ReadPositionMatrix(), file1.references, 15)
    return dict(file1=file1, rf_positions=rf_positions, shared=shared, index=index, coverage=coverage,
                coverage_index=coverage_index, dangling_sequences=dangling_sequences, buffers=buffers,
                row=RawArray(c_uint, buffer_size), col=RawArray(c_uint, buffer_size),
                data=RawArray(c_ushort, buffer_size))


def process(target, path1, path2, bin_size=5000, rf_paths=(DPNII,), restriction=('GATC',), dangling=('GATC',),
            buffer_size=2000, keep_circles=False, keep_ligation=False, quick=False, out_bam=False,
            coverage_bin=10, none_buffer=False, max_insert=1000):
    import queue
    module = module_of(target)
    inputs = process_inputs(module, path1, path2, bin_size, list(rf_paths), list(restriction) if restriction else None,
                            list(dangling) if dangling else None, buffer_size)
    b1, b2 = inputs['buffers'][0], inputs['buffers'][1]
    if none_buffer:
        b1 = None
    q = queue.Queue()
    target(pMateBuffer1=b1, pMateBuffer2=b2, pMinMappingQuality=15, pKeepSelfCircles=keep_circles,
           pRestrictionSequence=list(restriction) if restriction else None, pKeepSelfLigation=keep_ligation,
           pMatrixSize=len(inputs['shared']), pRfPositions=inputs['rf_positions'], pRefId2name=inputs['file1'].references,
           pDanglingSequences=inputs['dangling_sequences'], pBinsize=coverage_bin, pResultIndex=1, pQueueOut=q,
           pTemplate=inputs['file1'], pOutputBamSet=out_bam, pCounter=7, pSharedBinIntvalTree=inputs['shared'],
           pDictBinIntervalTreeIndex=inputs['index'], pCoverage=inputs['coverage'],
           pCoverageIndex=inputs['coverage_index'], pOutputFileBufferDir='', pRow=inputs['row'], pCol=inputs['col'],
           pData=inputs['data'], pMaxInsertSize=max_insert, pQuickQCMode=quick)
    result = q.get()
    if isinstance(result, str):
        return result.split('Traceback')[0]
    counts = result[0] if result[0] is not None else result[1]
    added = counts[15]
    isizes = None if b1 is None else [(r.isize, s.isize) for r, s in zip(b1, b2)]
    return (result[0] is None, counts, list(inputs['row'][:added]), list(inputs['col'][:added]),
            list(inputs['data'][:added]), [(i, v) for i, v in enumerate(inputs['coverage'][:]) if v], isizes)


# ---------------------------------------------------------------------------

CASES = [
    # C_Interval, C_Coverage ----------------------------------------------------
    Case('C_Interval/layout', M + 'C_Interval', driver='struct_layout'),
    Case('C_Coverage/layout', M + 'C_Coverage', driver='struct_layout'),

    # ReadPositionMatrix ---------------------------------------------------------
    Case('ReadPositionMatrix/doctest', M + 'ReadPositionMatrix',
         args=([('1', 0, '2', 0), ('1', 0, '2', 0)],), driver='position_matrix'),
    Case('ReadPositionMatrix/orders', M + 'ReadPositionMatrix',
         args=([('chr2', 500, 'chr1', 100), ('chr1', 100, 'chr2', 500), ('chr1', 500, 'chr2', 100),
                ('chr2', 100, 'chr1', 500), ('chrX', 7, 'chrX', 3), ('chrX', 3, 'chrX', 7), ('chrX', 3, 'chrX', 3),
                ('a-b', 1, 'c', 2), ('a', 1, 'b-c', 2), (1, 2, 3, 4), (np.int64(5), 6, np.int64(5), 6)],),
         driver='position_matrix'),
    Case('ReadPositionMatrix/type_error', M + 'ReadPositionMatrix',
         args=([('chr1', 1, 2, 3), ('chr1', 1, 'chr2', None)],), driver='position_matrix'),

    # intervalListToIntervalTree -------------------------------------------------
    Case('intervalListToIntervalTree/doctest', M + 'intervalListToIntervalTree',
         args=([('chrX', 0, 50000), ('chrX', 50000, 100000)],), driver='trees'),
    Case('intervalListToIntervalTree/chromosomes', M + 'intervalListToIntervalTree',
         args=([('chr2', 0, 10), ('chr1', 5, 9, 'extra'), ['chr2', 10, 20], ('chr1', 0, 5)],), driver='trees'),
    Case('intervalListToIntervalTree/empty', M + 'intervalListToIntervalTree', args=([],), driver='trees'),
    Case('intervalListToIntervalTree/null_interval', M + 'intervalListToIntervalTree',
         args=([('chr1', 0, 10), ('chr1', 5, 5)],), driver='trees'),
    Case('intervalListToIntervalTree/short_tuple', M + 'intervalListToIntervalTree',
         args=([('chr1', 0)],), driver='trees'),
    Case('intervalListToIntervalTree/numpy', M + 'intervalListToIntervalTree',
         args=lambda: ([('chr1', np.int64(0), np.int64(10)), ('chr1', 10, np.int64(20))],), driver='trees'),

    # get_bins --------------------------------------------------------------------
    Case('get_bins/doctest', M + 'get_bins', args=(50000, HIC_BAM), driver='bins_of_bam'),
    Case('get_bins/doctest_region', M + 'get_bins', args=(50000, HIC_BAM, 'contig-1'), driver='bins_of_bam'),
    Case('get_bins/real', M + 'get_bins', args=(1000000, R1), driver='bins_of_bam'),
    Case('get_bins/region_interval', M + 'get_bins', args=(1000, HIC_BAM, 'contig-1:2000:5500'), driver='bins_of_bam'),
    Case('get_bins/region_unknown', M + 'get_bins', args=(1000, HIC_BAM, 'chr9'), driver='bins_of_bam'),
    Case('get_bins/sizes', M + 'get_bins', args=(30, [('a', 100), ('b', 30), ('c', 0), ('d', 31)])),
    Case('get_bins/numpy_sizes', M + 'get_bins', args=lambda: (30, [('a', np.int64(100)), ('b', np.int64(29))])),
    Case('get_bins/numpy_bin_size', M + 'get_bins', args=lambda: (np.int64(40), [('a', 100), ('b', 80)])),
    Case('get_bins/region_empty_string', M + 'get_bins', args=(30, [('a', 100)], '')),
    Case('get_bins/zero', M + 'get_bins', args=(0, [('a', 100)])),
    Case('get_bins/negative', M + 'get_bins', args=(-10, [('a', 100)])),
    Case('get_bins/float', M + 'get_bins', args=(10.0, [('a', 100)])),
    Case('get_bins/bad_entry', M + 'get_bins', args=(10, [('a', 100, 3)])),

    # bed2interval_list -------------------------------------------------------------
    Case('bed2interval_list/doctest', M + 'bed2interval_list', args=(File('a.bed'), None, None),
         files={'a.bed': b'chr1\t10\t20\tH1\t0\nchr1\t60\t70\tH2\t0\n'}, driver='bed_list'),
    Case('bed2interval_list/real', M + 'bed2interval_list', args=(HINDIII, None, None), driver='bed_list'),
    Case('bed2interval_list/bytes', M + 'bed2interval_list', args=(HINDIII_SMALL, None, None, 'rb'),
         driver='bed_list'),
    Case('bed2interval_list/short_line', M + 'bed2interval_list', args=(File('a.bed'), None, None),
         files={'a.bed': b'chr1 10 20\n\nchr2\t5\nchr3\t1\t2\n'}, driver='bed_list'),
    Case('bed2interval_list/short_first', M + 'bed2interval_list', args=(File('a.bed'), None, None),
         files={'a.bed': b'chr1\t10\nchr1\t1\t2\n'}, driver='bed_list'),
    Case('bed2interval_list/not_int', M + 'bed2interval_list', args=(File('a.bed'), None, None),
         files={'a.bed': b'chr1\t10\t20\nchr1\tx\t30\n'}, driver='bed_list'),
    Case('bed2interval_list/int_forms', M + 'bed2interval_list', args=(File('a.bed'), None, None),
         files={'a.bed': b'chr1\t+10\t1_000\nchr1\t007\t30\n'}, driver='bed_list'),
    Case('bed2interval_list/region', M + 'bed2interval_list',
         args=(File('a.bed'), [('chr1', 100), ('chr2', 50)], 'chr1'),
         files={'a.bed': b'chr1\t10\t20\nchr1\t60\t100\nchr2\t1\t2\n'}, driver='bed_list'),
    Case('bed2interval_list/region_short_first', M + 'bed2interval_list',
         args=(File('a.bed'), [('chr1', 100)], 'chr1'), files={'a.bed': b'\nchr1\t1\t2\n'}, driver='bed_list'),
    Case('bed2interval_list/region_unknown', M + 'bed2interval_list',
         args=(File('a.bed'), [('chr1', 100)], 'chr5'), files={'a.bed': b'chr1\t1\t2\n'}, driver='bed_list'),
    Case('bed2interval_list/list_of_lines', M + 'bed2interval_list',
         args=(['chr1 1 2\n', '  chr2 3\t4  \n'], None, None), driver='bed_from_list'),

    # get_rf_bins ---------------------------------------------------------------------
    Case('get_rf_bins/doctest', M + 'get_rf_bins',
         args=([('chr1', 10, 20), ('chr1', 60, 70), ('chr2', 20, 30), ('chr2', 40, 50), ('chr2', 70, 80)],),
         kwargs={'min_distance': 10, 'max_distance': 20}),
    Case('get_rf_bins/defaults', M + 'get_rf_bins',
         args=([('chr1', 1000, 1004), ('chr1', 1100, 1104), ('chr1', 5000, 5004), ('chr2', 900, 904)],)),
    Case('get_rf_bins/real', M + 'get_rf_bins', args=(HINDIII, 150, 1500), driver='rf_bins_of_bed'),
    Case('get_rf_bins/real_defaults', M + 'get_rf_bins', args=(HINDIII_SMALL,), driver='rf_bins_of_bed'),
    Case('get_rf_bins/numpy', M + 'get_rf_bins',
         args=lambda: ([('c', np.int64(500), np.int64(504)), ('c', 2000, 2004)], np.int64(10), 300)),
    Case('get_rf_bins/empty', M + 'get_rf_bins', args=([],)),
    Case('get_rf_bins/two_columns', M + 'get_rf_bins', args=([('c', 1), ('c', 2)],)),
    Case('get_rf_bins/four_columns', M + 'get_rf_bins', args=([('c', 1, 2, 3)],)),
    Case('get_rf_bins/min_none', M + 'get_rf_bins', args=([('c', 1, 2)], None, 3)),

    # get_chrom_sizes ------------------------------------------------------------------
    Case('get_chrom_sizes/doctest', M + 'get_chrom_sizes', args=(HIC_BAM,), driver='chrom_sizes_of'),
    Case('get_chrom_sizes/real', M + 'get_chrom_sizes', args=(R1,), driver='chrom_sizes_of'),
    Case('get_chrom_sizes/repeated', M + 'get_chrom_sizes',
         args=(('chr1', 'chr2', 'chr1', 'chr3'), (10, 20, 30, 40)), driver='chrom_sizes_fake'),
    Case('get_chrom_sizes/zip_shorter', M + 'get_chrom_sizes', args=(('a', 'b', 'c'), [1, 2]),
         driver='chrom_sizes_fake'),

    # check_dangling_end ---------------------------------------------------------------
    Case('check_dangling_end/real', M + 'check_dangling_end', args=(R1, {'pat_forw': 'GATC', 'pat_rev': 'GATC'}),
         driver='dangling_reads'),
    Case('check_dangling_end/real_agct', M + 'check_dangling_end',
         args=(R2, {'pat_forw': 'AGCT', 'pat_rev': 'AGCT'}), driver='dangling_reads'),
    Case('check_dangling_end/missing_key', M + 'check_dangling_end', args=({'pat_forw': 'A'},),
         kwargs={'is_reverse': False, 'seq': 'AC'}, driver='dangling_fake'),
    Case('check_dangling_end/forward', M + 'check_dangling_end', args=({'pat_forw': 'GATC', 'pat_rev': 'TT'},),
         kwargs={'is_reverse': False, 'seq': 'gatcAA'}, driver='dangling_fake'),
    Case('check_dangling_end/reverse', M + 'check_dangling_end', args=({'pat_forw': 'GATC', 'pat_rev': 'TTA'},),
         kwargs={'is_reverse': 1, 'seq': 'GATCtta'}, driver='dangling_fake'),
    Case('check_dangling_end/reverse_no', M + 'check_dangling_end', args=({'pat_forw': 'GATC', 'pat_rev': 'TTA'},),
         kwargs={'is_reverse': True, 'seq': 'GATC'}, driver='dangling_fake'),
    Case('check_dangling_end/tuple', M + 'check_dangling_end', args=({'pat_forw': ('X', 'GA'), 'pat_rev': 'T'},),
         kwargs={'is_reverse': False, 'seq': 'GATC'}, driver='dangling_fake'),
    Case('check_dangling_end/empty_pattern', M + 'check_dangling_end', args=({'pat_forw': '', 'pat_rev': ''},),
         kwargs={'is_reverse': False, 'seq': ''}, driver='dangling_fake'),
    Case('check_dangling_end/unicode', M + 'check_dangling_end', args=({'pat_forw': 'SS', 'pat_rev': ''},),
         kwargs={'is_reverse': False, 'seq': 'ß!'}, driver='dangling_fake'),
    Case('check_dangling_end/none_seq', M + 'check_dangling_end', args=({'pat_forw': 'A', 'pat_rev': 'A'},),
         kwargs={'is_reverse': False, 'seq': None}, driver='dangling_fake'),
    Case('check_dangling_end/int_pattern', M + 'check_dangling_end', args=({'pat_forw': 5, 'pat_rev': 'A'},),
         kwargs={'is_reverse': False, 'seq': 'A'}, driver='dangling_fake'),

    # get_supplementary_alignment -----------------------------------------------------
    Case('get_supplementary_alignment/real_none', M + 'get_supplementary_alignment', args=(R1,),
         driver='supplementary'),
    Case('get_supplementary_alignment/sam', M + 'get_supplementary_alignment', args=(File('s.sam'),),
         files={'s.sam': SA_SAM}, driver='supplementary'),
    Case('get_supplementary_alignment/other_name', M + 'get_supplementary_alignment', args=(File('s.sam'), 3),
         files={'s.sam': SA_SAM}, driver='supplementary'),
    Case('get_supplementary_alignment/fake', M + 'get_supplementary_alignment',
         args=({'SA': 'a;b;c;'}, ['a', 'b', 'a', 'z']), driver='supplementary_fake'),
    Case('get_supplementary_alignment/empty_tag', M + 'get_supplementary_alignment',
         args=({'SA': ''}, ['a']), driver='supplementary_fake'),
    Case('get_supplementary_alignment/exhausted', M + 'get_supplementary_alignment',
         args=({'SA': 'x;y;'}, ['a']), driver='supplementary_fake'),
    Case('get_supplementary_alignment/int_tag', M + 'get_supplementary_alignment',
         args=({'SA': 3}, ['a']), driver='supplementary_fake'),

    # get_correct_map ------------------------------------------------------------------
    Case('get_correct_map/doctest', M + 'get_correct_map', args=(File('s.sam'),), files={'s.sam': SA_SAM},
         driver='correct_map'),
    Case('get_correct_map/reverse', M + 'get_correct_map', args=(File('s.sam'),), files={'s.sam': REVERSE_SAM},
         driver='correct_map'),
    Case('get_correct_map/tuple', M + 'get_correct_map', args=(File('s.sam'), True), files={'s.sam': SA_SAM},
         driver='correct_map'),
    Case('get_correct_map/other_name', M + 'get_correct_map', args=(File('s.sam'),), files={'s.sam': REVERSE_SAM},
         driver='correct_map_mixed'),
    Case('get_correct_map/fake_ties', M + 'get_correct_map',
         args=((False, [(4, 5), (0, 10)]), [(False, [(4, 5), (0, 3)]), (True, [(0, 1), (4, 5)]), (False, [(4, 2)])]),
         driver='correct_map_fake'),
    Case('get_correct_map/fake_no_match_op', M + 'get_correct_map',
         args=((False, [(4, 50)]), [(False, [(5, 3), (1, 2), (0, 9)])]), driver='correct_map_fake'),
    Case('get_correct_map/none_cigar', M + 'get_correct_map', args=((False, None), []), driver='correct_map_fake'),

    # enlarge_bins ----------------------------------------------------------------------
    Case('enlarge_bins/doctest', M + 'enlarge_bins',
         args=([('chr1', 10, 30), ('chr1', 50, 80), ('chr2', 10, 60), ('chr2', 60, 90)],
               [('chr1', 100), ('chr2', 100)]), driver='call_and_list_state'),
    Case('enlarge_bins/touching_first', M + 'enlarge_bins',
         args=([('a', 10, 30), ('a', 30, 80), ('a', 85, 99), ('b', 5, 9), ('b', 20, 30), ('c', 3, 4)],
               [('a', 100), ('b', 50), ('c', 7)]), driver='call_and_list_state'),
    Case('enlarge_bins/single', M + 'enlarge_bins', args=([('a', 10, 30)], [('a', 100)]),
         driver='call_and_list_state'),
    Case('enlarge_bins/lists_and_numpy', M + 'enlarge_bins',
         args=lambda: ([['a', np.int64(10), np.int64(30)], ['a', np.int64(41), 80], ('a', 90, np.int64(95)),
                        ['b', 1, 2]], [('a', np.int64(100)), ('b', 9)]), driver='call_and_list_state'),
    Case('enlarge_bins/overlap', M + 'enlarge_bins', args=([('a', 10, 50), ('a', 41, 80)], [('a', 100)]),
         driver='call_and_list_state'),
    Case('enlarge_bins/real', M + 'enlarge_bins', args=(HINDIII_SMALL, HIC_BAM), driver='enlarge_rf'),
    Case('enlarge_bins/missing_size', M + 'enlarge_bins',
         args=([('a', 10, 30), ('a', 40, 80), ('b', 1, 2), ('b', 3, 4)], [('b', 10)]), driver='enlarge_partial'),
    Case('enlarge_bins/missing_last', M + 'enlarge_bins', args=([('a', 10, 30), ('a', 40, 80)], [('b', 10)]),
         driver='enlarge_partial'),
    Case('enlarge_bins/empty', M + 'enlarge_bins', args=([], [('a', 1)])),

    # readBamFiles ----------------------------------------------------------------------
    Case('readBamFiles/real_one_buffer', M + 'readBamFiles', args=(R1, R2, 400000, False, 15),
         driver='read_buffers'),
    Case('readBamFiles/real_buffers', M + 'readBamFiles', args=(R1, R2, 70, False, 15), driver='read_buffers'),
    Case('readBamFiles/exact_multiple', M + 'readBamFiles', args=(R1, R2, 50, False, 15), driver='read_buffers'),
    Case('readBamFiles/skip_check_none', M + 'readBamFiles', args=(R1, R2, 400000, None, 0), driver='read_buffers'),
    Case('readBamFiles/skip_check_true', M + 'readBamFiles', args=(SMALL_R1, SMALL_R2, 3000, True, 30, 3),
         driver='read_buffers'),
    Case('readBamFiles/duplicates', M + 'readBamFiles', args=(SMALL_R1, SMALL_R2, 5000, False, 15, 2),
         driver='read_buffers'),
    Case('readBamFiles/min_quality_high', M + 'readBamFiles', args=(R1, R2, 400000, False, 255),
         driver='read_buffers'),
    Case('readBamFiles/order', M + 'readBamFiles', args=(R1, HIC_BAM), driver='read_buffers_error'),

    # process_data ---------------------------------------------------------------------------
    Case('process_data/bins', M + 'process_data', args=(R1, R2), driver='process'),
    Case('process_data/rf_bins', M + 'process_data', args=(R1, R2), kwargs={'bin_size': None}, driver='process'),
    Case('process_data/no_restriction', M + 'process_data', args=(R1, R2),
         kwargs={'rf_paths': (), 'restriction': None, 'dangling': None}, driver='process'),
    Case('process_data/two_enzymes', M + 'process_data', args=(R1, R2),
         kwargs={'rf_paths': (DPNII, HINDIII), 'restriction': ('GATC', 'AAGCTT'), 'dangling': ('GATC', 'AGCT')},
         driver='process'),
    Case('process_data/keep', M + 'process_data', args=(R1, R2),
         kwargs={'keep_circles': True, 'keep_ligation': True, 'out_bam': True}, driver='process'),
    Case('process_data/quick_qc', M + 'process_data', args=(R1, R2), kwargs={'quick': True, 'max_insert': 5000},
         driver='process'),
    Case('process_data/coverage_bin', M + 'process_data', args=(R1, R2), kwargs={'coverage_bin': 3},
         driver='process'),
    Case('process_data/coverage_bin_zero', M + 'process_data', args=(R1, R2), kwargs={'coverage_bin': 0},
         driver='process'),
    Case('process_data/none_buffer', M + 'process_data', args=(R1, R2), kwargs={'none_buffer': True},
         driver='process'),
    Case('process_data/small_test', M + 'process_data', args=(SMALL_R1, SMALL_R2),
         kwargs={'buffer_size': 20000}, driver='process'),

]
