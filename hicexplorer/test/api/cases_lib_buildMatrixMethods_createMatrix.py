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
    4 runs the hicPrepareQCreport executable;
  * the debug records of hicexplorer.reduceMatrix (timings) and
    hicexplorer.hicMergeMatrixBins ("few bins" notes) that 3.7.6 wrote while
    reducing the matrix to the further --binSize resolutions of an mcool; the
    reduction runs in C++ and writes the same matrices.
The files hicPrepareQCreport draws (PNG, HTML) are compared by name only.

Fix 3 of HiCExplorer 4 (the project owner's decision of 2026-11-06) changes
what is written, so the cases it touches carry a `deviation` (see
test_api_equivalence.py) in place of the comparison with 3.7.6:
  * the provenance attributes (matrix-generated-by, matrix-generated-by-url,
    genome-assembly) are plain strings, where 3.7.6 stored the repr
    "b'HiCExplorer-3.7.6'", and every mcool resolution carries them, where
    3.7.6 gave them to the first resolution only. They are removed from both
    records and checked on the new one, as is the "hicexplorer_qc" key the
    v4 tools add to the metadata of a cool file (3.7.6 writes none); everything
    else still equals 3.7.6;
  * an .mcool name writes one resolution per bin size from one bin size on
    (3.7.6: only for three or more; two bin sizes gave a single cooler at the
    file root);
  * --doTestRun with an .mcool name and three or more bin sizes no longer ends
    in the UnboundLocalError the 3.7.6 code path raised.

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
    qc_loggers = [logging.getLogger(n) for n in ('hicexplorer.hicPrepareQCreport', 'hicexplorer.reduceMatrix',
                                                 'hicexplorer.hicMergeMatrixBins')]
    own.addFilter(timing)
    for qc_logger in qc_loggers:
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
        for qc_logger in qc_loggers:
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


# Fix 3 of HiCExplorer 4, see the module docstring.
_PROVENANCE = ('matrix-generated-by', 'matrix-generated-by-url', 'genome-assembly')
# The QC statistics the v4 tools write into the metadata of a cool file (3.7.6 writes none).
_V4_ONLY = ('hicexplorer_qc',)


def _is_provenance(item):
    return (isinstance(item, tuple) and len(item) == 2 and item[0] == 'tuple' and len(item[1]) == 2
            and item[1][0][0] == 'str' and item[1][0][1] in _PROVENANCE)


def _strip(value):
    """The canonical value without the provenance attributes and without the
    provenance keys of the JSON `metadata` attribute."""
    import json
    if isinstance(value, tuple) and len(value) == 2 and value[0] in ('list', 'tuple'):
        items = []
        for item in value[1]:
            if _is_provenance(item):
                continue
            if (isinstance(item, tuple) and item[0] == 'tuple' and len(item[1]) == 2
                    and item[1][0] == ('str', 'metadata') and item[1][1][0] == 'str'):
                metadata = json.loads(item[1][1][1])
                metadata = {k: v for k, v in metadata.items() if k not in _PROVENANCE + _V4_ONLY}
                item = ('tuple', [item[1][0], ('str', json.dumps(metadata, sort_keys=True))])
            items.append(_strip(item))
        return (value[0], items)
    if isinstance(value, tuple) and len(value) == 2 and value[0] == 'dict':
        return ('dict', [(k, _strip(v)) for k, v in value[1]])
    return value


def _provenance(value, found):
    """Every (key, value) of the provenance in attributes and metadata JSON."""
    import json
    if isinstance(value, tuple) and len(value) == 2 and value[0] in ('list', 'tuple'):
        if _is_provenance(value):
            found.append(('attribute', value[1][0][1], value[1][1][1]))
        elif len(value[1]) == 2 and value[1][0] == ('str', 'metadata') and value[1][1][0] == 'str':
            for key, content in json.loads(value[1][1][1]).items():
                if key in _PROVENANCE:
                    found.append(('metadata', key, content))
        for item in value[1]:
            _provenance(item, found)
    elif isinstance(value, tuple) and len(value) == 2 and value[0] == 'dict':
        for _, item in value[1]:
            _provenance(item, found)
    return found


def _plain(record, resolutions, assembly=None):
    found = _provenance(record['return'][1][6], [])
    versions = [v for kind, key, v in found if key == 'matrix-generated-by']
    if len(versions) < resolutions or any(v != 'HiCExplorer-4.0.0' for v in versions):
        return 'matrix-generated-by is not a plain string on every resolution: {}'.format(versions)
    urls = [v for kind, key, v in found if key == 'matrix-generated-by-url']
    if any(v != 'https://github.com/deeptools/HiCExplorer' for v in urls):
        return 'matrix-generated-by-url is not a plain string: {}'.format(urls)
    if assembly is not None:
        values = [v for kind, key, v in found if key == 'genome-assembly']
        if not values or any(v != assembly for v in values):
            return 'genome-assembly is not the plain string {!r}: {}'.format(assembly, values)
    return None


def provenance(resolutions, assembly=None):
    def check(expected, actual):
        import apieq
        message = _plain(actual, resolutions, assembly)
        if message:
            return message
        exp, act = dict(expected), dict(actual)
        exp['return'] = ('tuple', [_strip(x) for x in expected['return'][1]])
        act['return'] = ('tuple', [_strip(x) for x in actual['return'][1]])
        return apieq.compare_records(exp, act, _NoTolerance)
    return check


class _NoTolerance(object):
    rtol = 0.0
    atol = 0.0


def two_resolutions(expected, actual):
    """An .mcool of two bin sizes: both resolutions exist (3.7.6 wrote a
    single cooler at the file root); the rest of the record equals 3.7.6."""
    import apieq
    message = _plain(actual, 2)
    if message:
        return message
    names = [item[1][0][1] for item in actual['return'][1][6][1][0][1][1]
             if isinstance(item, tuple) and item[0] == 'tuple' and item[1][0][0] == 'str']
    for resolution in ('resolutions/5000', 'resolutions/10000'):
        if resolution not in names:
            return 'no group {} in {}'.format(resolution, names)
    exp, act = dict(expected), dict(actual)
    exp['return'] = ('tuple', [x if i != 6 else ('const', 'None') for i, x in enumerate(expected['return'][1])])
    act['return'] = ('tuple', [x if i != 6 else ('const', 'None') for i, x in enumerate(actual['return'][1])])
    return apieq.compare_records(exp, act, _NoTolerance)


def do_test_run_returns(expected, actual):
    """--doTestRun with an .mcool name: 3.7.6 ended in an UnboundLocalError
    ('hic_ma'); now the call returns like the one with any other name."""
    outcome = actual['return'][1][0]
    if outcome != ('tuple', [('str', 'returned'), ('const', 'None')]):
        return 'the call did not return: {}'.format(outcome)
    if actual['return'][1][2][0] != 'str' or 'Sequenced reads' not in actual['return'][1][2][1]:
        return 'no QC.log was written'
    return None


for _case in CASES:
    if _case.id == 'createMatrix/mcool':
        _case.deviation = provenance(3)
    elif _case.id == 'createMatrix/cool_assembly':
        _case.deviation = provenance(1, 'dm3')
    elif _case.id == 'createMatrix/mcool_two_resolutions':
        _case.deviation = two_resolutions
    else:
        continue
    _case.note = 'Fix 3 of HiCExplorer 4 (project owner 2026-11-06), see the module docstring'

CASES.append(Case('createMatrix/mcool_do_test_run', M + 'createMatrix', args=(File(),),
                  kwargs={'out_name': 'matrix.mcool', 'bin_size': (5000, 10000, 20000), 'test_run': True},
                  driver='create'))
CASES[-1].deviation = do_test_run_returns
CASES[-1].note = 'Fix 3 of HiCExplorer 4 (project owner 2026-11-06), see the module docstring'
