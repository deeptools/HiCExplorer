"""Equivalence cases for hicexplorer.lib.buildMatrixMethods.createMatrix
(apieq.py), apart from cases_lib_buildMatrixMethods.py because every case
runs a whole hicBuildMatrix (3.7.6 needs about ten seconds for one).

The inputs are the repository's real alignments (hicBuildMatrix/R1_1000.bam
and R2_1000.bam, small_test_R1/R2_unsorted.bam), restriction site files
(DpnII.bed, hicFindRestSite/hindIII.bed) and dm3.chrom.sizes, passed as open
file objects as hicBuildMatrix's argparse passes them.

createMatrix writes a matrix and a QC folder: the drivers return QC.log, the
datasets and attributes of the matrix file, the QC folder listing, the
--outBam records, the state of the file objects and lists passed in, and the
warnings raised. Two kinds of log records are not compared:
  * the timing of "processing N lines took X secs (Y lines per second)",
    which a driver filter replaces by placeholders; with several workers the
    numbers in both progress records depend on the order in which 3.7.6's
    worker processes finished, so there the whole message is masked;
  * the debug records of hicexplorer.hicPrepareQCreport: createMatrix calls
    hicPrepareQCreport.main, which in 3.7.6 ran in-process and in HiCExplorer
    4 runs the hicPrepareQCreport executable.
The files hicPrepareQCreport draws (PNG, HTML) are compared by name only.

Deviations from 3.7.6 (cpp/python/src/api/lib_buildMatrixMethods_run.hpp):
  * the schedule of 3.7.6's worker processes is emulated as every worker
    holding a buffer and results taken in reading order, which decides the
    --outBam record order with several buffers and workers (cpp/tools/
    hicBuildMatrix.cpp Q8) and how many buffers a --doTestRun counts there;
  * BAM files htslib cannot open raise RuntimeError with htslib's message
    instead of pysam's OSError;
  * a restriction sequence given twice is classified as 3.7.6 classifies it
    and its QC rows are merged as 3.7.6's dictionary merges them;
  * the timing of the progress records: 3.7.6 waits a second whenever no
    worker has finished.
"""
import os
import re

import numpy as np

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'test_data')
BM = os.path.join(DATA, 'hicBuildMatrix')
R1 = os.path.join(BM, 'R1_1000.bam')
R2 = os.path.join(BM, 'R2_1000.bam')
DPNII = os.path.join(BM, 'DpnII.bed')
HINDIII = os.path.join(DATA, 'hicFindRestSite', 'hindIII.bed')
SIZES = os.path.join(BM, 'dm3.chrom.sizes')
SMALL_R1 = os.path.join(DATA, 'small_test_R1_unsorted.bam')
SMALL_R2 = os.path.join(DATA, 'small_test_R2_unsorted.bam')
M = 'hicexplorer.lib.buildMatrixMethods.'


# createMatrix

class _Timing(object):
    """Masks the timing (or, with several buffers, the numbers) of the
    progress records; drops hicPrepareQCreport's records."""

    def __init__(self, mask_all):
        self.mask_all = mask_all

    def filter(self, record):
        message = record.getMessage()
        if message.startswith('processing ') and 'lines took' in message:
            message = re.sub(r'took [0-9.]+ secs \([0-9.e+]+ lines per second\)',
                             'took <secs> secs (<rate> lines per second)', message)
            if self.mask_all:
                message = 'processing <lines> lines took <secs>'
            record.msg, record.args = message, None
        elif self.mask_all and message.endswith('valid pairs added to matrix\n'):
            record.msg, record.args = '<pairs> valid pairs added to matrix', None
        return True


def _drop(record):
    return False


def matrix_file(path):
    """The datasets and attributes of an h5 or cool file."""
    import h5py
    if not os.path.exists(path):
        return None
    out = []

    def attributes(obj):
        items = []
        for key in sorted(obj.attrs):
            if key in ('creation-date', 'CLASS', 'VERSION', 'TITLE', 'PYTABLES_FORMAT_VERSION', 'FILTERS'):
                continue
            value = obj.attrs[key]
            if isinstance(value, bytes):
                value = value.decode()
            if isinstance(value, np.ndarray):
                value = value.tolist()
            if isinstance(value, np.generic):
                value = value.item()
            items.append((key, value))
        return items

    with h5py.File(path, 'r') as f:
        out.append(('/', attributes(f)))

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                data = obj[()]
                if isinstance(data, np.ndarray):
                    data = [v.decode() if isinstance(v, bytes) else v for v in data.tolist()]
                elif isinstance(data, bytes):
                    data = data.decode()
                out.append((name, str(obj.dtype) if obj.dtype.kind != 'S' else 'S', data, attributes(obj)))
            else:
                out.append((name, attributes(obj)))
        f.visititems(visit)
    return out


def bam_records(path):
    if not path or not os.path.exists(path):
        return None
    if os.path.getsize(path) == 0:
        return []
    import pysam
    with pysam.AlignmentFile(path, 'rb') as f:
        return [r.to_string() for r in f]


def create(target, directory, out_name='matrix.h5', sam_files=(R1, R2), bin_size=(5000,), rf=(DPNII,),
           restriction=('GATC',), dangling=('GATC',), threads=4, out_bam=None, chromosome_sizes=None,
           region=None, buffer_size=400000, test_run=False, test_lines=1000000, skip_dup=False, min_mapq=15,
           keep_circles=False, keep_ligation=False, min_distance=300, max_distance=None, max_insert=1000,
           assembly=None, qc='QC', mask_all=False):
    import logging
    import warnings
    out_path = os.path.join(directory, out_name)
    qc_folder = os.path.join(directory, qc)
    handles = {'out': open(out_path, 'w'), 'sam': [open(p, 'rb') for p in sam_files]}
    handles['bam'] = open(os.path.join(directory, out_bam), 'w') if out_bam else None
    handles['sizes'] = open(chromosome_sizes) if chromosome_sizes else None
    handles['rf'] = [open(p) for p in rf] if rf is not None else None
    restriction_list = list(restriction) if restriction is not None else None
    dangling_list = list(dangling) if dangling is not None else None
    timing = _Timing(mask_all)
    own = logging.getLogger('hicexplorer.lib.buildMatrixMethods')
    qc_logger = logging.getLogger('hicexplorer.hicPrepareQCreport')
    own.addFilter(timing)
    qc_logger.addFilter(_drop)
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            try:
                result = target(pOutFileName=handles['out'], pMaxDistance=max_distance, pMaxLibraryInsertSize=max_insert,
                                pQCfolder=qc_folder, pThreads=threads, pDanglingSequence=dangling_list,
                                pRestrictionSequence=restriction_list, pSamFiles=handles['sam'],
                                pDoTestRun=test_run, pOutBam=handles['bam'], pChromosomeSizes=handles['sizes'],
                                pRestrictionCutFile=handles['rf'], pRegion=region,
                                pBinSize=list(bin_size) if bin_size is not None else None,
                                pInputBufferSize=buffer_size, pDoTestRunLines=test_lines,
                                pSkipDuplicationCheck=skip_dup, pMinMappingQuality=min_mapq,
                                pKeepSelfCircles=keep_circles, pKeepSelfLigation=keep_ligation,
                                pMinDistance=min_distance, pGenomeAssembly=assembly)
                outcome = ('returned', result)
            except SystemExit as error:
                outcome = ('SystemExit', error.code)
            except Exception as error:
                outcome = (type(error).__name__, str(error))
            messages = [(w.category.__name__, str(w.message)) for w in caught
                        if issubclass(w.category, UserWarning) and 'threads' in str(w.message)]
    finally:
        own.removeFilter(timing)
        qc_logger.removeFilter(_drop)
    qc_log = None
    if os.path.exists(os.path.join(qc_folder, 'QC.log')):
        with open(os.path.join(qc_folder, 'QC.log')) as f:
            qc_log = f.read()
    listing = sorted(os.listdir(qc_folder)) if os.path.isdir(qc_folder) else None
    closed = {'out': handles['out'].closed, 'sam': [h.closed for h in handles['sam']],
              'bam': None if handles['bam'] is None else handles['bam'].closed,
              'sizes': None if handles['sizes'] is None else handles['sizes'].closed}
    files = {'matrix': matrix_file(out_path)}
    files['bam'] = bam_records(os.path.join(directory, out_bam)) if out_bam else None
    for handle in [handles['out'], handles['bam'], handles['sizes']] + handles['sam'] + (handles['rf'] or []):
        if handle is not None:
            handle.close()
    return (outcome, messages, qc_log, listing, closed, os.path.exists(out_path), files,
            restriction_list, dangling_list)


CASES = [
    # createMatrix ------------------------------------------------------------------------------
    Case('createMatrix/h5', M + 'createMatrix', args=(File(),), driver='create'),
    Case('createMatrix/cool_assembly', M + 'createMatrix', args=(File(),),
         kwargs={'out_name': 'matrix.cool', 'assembly': 'dm3'}, driver='create'),
    Case('createMatrix/mcool', M + 'createMatrix', args=(File(),),
         kwargs={'out_name': 'matrix.mcool', 'bin_size': (5000, 10000, 20000)}, driver='create'),
    Case('createMatrix/mcool_two_resolutions', M + 'createMatrix', args=(File(),),
         kwargs={'out_name': 'matrix.mcool', 'bin_size': (5000, 10000)}, driver='create'),
    Case('createMatrix/rf_bins', M + 'createMatrix', args=(File(),),
         kwargs={'bin_size': None, 'min_distance': 150, 'max_insert': 1500}, driver='create'),
    Case('createMatrix/two_enzymes_lowercase', M + 'createMatrix', args=(File(),),
         kwargs={'rf': (DPNII, HINDIII), 'restriction': ('gatc', 'aagctt'), 'dangling': ('gatc', 'agct')},
         driver='create'),
    Case('createMatrix/region', M + 'createMatrix', args=(File(),),
         kwargs={'rf': (DPNII, HINDIII), 'restriction': ('GATC', 'AAGCTT'), 'dangling': ('GATC', 'AGCT'),
                 'region': 'chr3R'}, driver='create'),
    Case('createMatrix/out_bam_sizes', M + 'createMatrix', args=(File(),),
         kwargs={'out_bam': 'out.bam', 'chromosome_sizes': SIZES, 'keep_circles': True, 'keep_ligation': True,
                 'skip_dup': True, 'max_distance': 1500}, driver='create'),
    Case('createMatrix/buffers', M + 'createMatrix', args=(File(),),
         kwargs={'buffer_size': 40, 'threads': 1}, driver='create'),
    Case('createMatrix/buffers_parallel', M + 'createMatrix', args=(File(),),
         kwargs={'buffer_size': 70, 'threads': 4, 'mask_all': True}, driver='create'),
    Case('createMatrix/test_run', M + 'createMatrix', args=(File(),),
         kwargs={'test_run': True, 'test_lines': 100, 'threads': 2, 'buffer_size': 70}, driver='create'),
    Case('createMatrix/test_run_out_bam', M + 'createMatrix', args=(File(),),
         kwargs={'test_run': True, 'out_bam': 'out.bam'}, driver='create'),
    Case('createMatrix/micro_c', 'hicexplorer.hicBuildMatrixMicroC.createMatrix', args=(File(),),
         kwargs={'rf': None, 'restriction': None, 'dangling': None, 'min_distance': None, 'keep_ligation': None,
                 'bin_size': (10000,), 'qc': 'QC/nested/'}, driver='create'),
    Case('createMatrix/small_test', M + 'createMatrix', args=(File(),),
         kwargs={'sam_files': (SMALL_R1, SMALL_R2), 'out_bam': 'out.bam'}, driver='create'),
    Case('createMatrix/bad_extension', M + 'createMatrix', args=(File(),), kwargs={'out_name': 'matrix.txt'},
         driver='create'),
    Case('createMatrix/dangling_without_restriction', M + 'createMatrix', args=(File(),),
         kwargs={'restriction': None}, driver='create'),
    Case('createMatrix/restriction_without_dangling', M + 'createMatrix', args=(File(),),
         kwargs={'dangling': None}, driver='create'),
    Case('createMatrix/short_dangling', M + 'createMatrix', args=(File(),),
         kwargs={'restriction': ('gatc', 'aagctt'), 'rf': (DPNII, HINDIII)}, driver='create'),
    Case('createMatrix/no_pair_passes', M + 'createMatrix', args=(File(),), kwargs={'min_mapq': 255},
         driver='create'),
    Case('createMatrix/exact_multiple_buffers', M + 'createMatrix', args=(File(),),
         kwargs={'buffer_size': 50, 'threads': 2}, driver='create'),
]
