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

Fixes Q1 to Q6 of hicBuildMatrix (cpp/tools/hicBuildMatrix.cpp) change the
numbers of the QC log and the matrix, so every case that classifies read pairs
carries the `deviation` `fixed_counts` (see test_api_equivalence.py) in place
of the comparison with 3.7.6:
  * Q1: a pair is "One mate not unique" when both mates have quality 0 (3.7.6:
    when the first has), so "Low mapping quality" gains what "One mate not
    unique" loses;
  * Q2: self circles are removed unless --keepSelfCircles and counted once, so
    "Hi-C contacts", the outward pairs and the short range pairs lose the
    removed pairs, the matrix loses them too, and --outBam holds their records
    no more;
  * Q3 and Q4: the coverage column of the bin table (intervals/extra_list) is
    not compared;
  * Q5: the --outBam records carry the insert size of the pair;
  * Q6: with --region the restriction cut sites of the chromosome are kept, so
    pairs move between "same fragment", "self ligation" and the dangling end
    counts (the `region` case compares the sums, not the single counts).
Everything else still equals 3.7.6: the other QC counts, the shape and names of
every dataset, the bin table, the pixels' sum (3.7.6 minus the removed pairs),
the QC folder listing, the log records (apart from the pair counts) and the
state of the file objects.

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


# Fixes Q1 to Q6 of hicBuildMatrix, see the module docstring.
_QC_COUNT = re.compile(r'^([^\t\n]+)\t(-?[0-9]+)\t', re.M)
_CHANGING = ('Hi-C contacts', 'One mate not unique', 'Low mapping quality', 'self circle',
             'Intra short range (< 20kb)', 'Intra long range (>= 20kb)',
             'Read pair type: outward pairs')
# the labels the --region fix (Q6) moves pairs between
_REGION_CHANGING = ('self ligation (removed)', 'self ligation (not removed)', 'same fragment')
_MASKED_DATA = ('matrix/data', 'matrix/indices', 'pixels/bin1_id', 'pixels/bin2_id', 'pixels/count')
_LENGTH_ONLY = ('matrix/indptr', 'indexes/bin1_offset', 'intervals/extra_list')
_MASKED_ATTRIBUTES = ('nnz', 'sum')


def _qc_counts(text):
    return [(label, int(count)) for label, count in _QC_COUNT.findall(text)]


def _check_qc(old, new, region):
    """(None, pairs lost) when the QC.log texts differ only as the fixes
    demand, else (reason, 0)."""
    old_counts, new_counts = _qc_counts(old), _qc_counts(new)
    if [label for label, _ in old_counts] != [label for label, _ in new_counts]:
        return 'the QC.log rows changed: {} against {}'.format(
            [label for label, _ in new_counts], [label for label, _ in old_counts]), 0
    o, n = dict(old_counts), dict(new_counts)
    for label in o:
        if label in _CHANGING or (region and (label in _REGION_CHANGING or label.startswith('dangling end'))):
            continue
        if o[label] != n[label]:
            return 'the count of {!r} changed: {} against {}'.format(label, n[label], o[label]), 0
    if o['One mate not unique'] + o['Low mapping quality'] != n['One mate not unique'] + n['Low mapping quality']:
        return 'the not unique and low quality pairs do not add up to the 3.7.6 sum', 0
    if n['One mate not unique'] > o['One mate not unique']:
        return 'more pairs are not unique than in 3.7.6 (Q1 looks at both mates now)', 0
    removed = o['Hi-C contacts'] - n['Hi-C contacts']
    if region:
        return None, removed
    if removed < 0:
        return 'the run kept more pairs than 3.7.6', 0
    if 'self circle' in n:
        if n['self circle'] > o['self circle']:
            return 'self circles are counted more often than in 3.7.6 (Q2)', 0
        if removed > n['self circle']:
            return 'pairs lost {} but the self circle count is {} (Q2)'.format(removed, n['self circle']), 0
    elif removed:
        return 'pairs were lost without a restriction site (Q2)', 0
    if n['Read pair type: outward pairs'] != o['Read pair type: outward pairs'] - removed:
        return 'the outward pairs do not lose exactly the removed pairs', 0
    lost = (o['Intra short range (< 20kb)'] - n['Intra short range (< 20kb)'] +
            o['Intra long range (>= 20kb)'] - n['Intra long range (>= 20kb)'])
    if lost != removed:
        return 'short and long range pairs do not lose exactly the removed pairs', 0
    return None, removed


def _files(record):
    """The files dict of a createMatrix record, {name: canonical value}."""
    value = record['return'][1][6]
    if isinstance(value, tuple) and value[0] == 'dict':
        return {key[1]: item for key, item in value[1]}
    return None


def _sums(matrix):
    """Sum of the pixel counts per resolution (the prefix of the dataset name)."""
    sums = {}
    if matrix is None or matrix[0] != 'list':
        return sums
    for entry in matrix[1]:
        if entry[0] == 'tuple' and len(entry[1]) == 4:
            name = entry[1][0][1]
            for suffix in ('matrix/data', 'pixels/count'):
                if name.endswith(suffix):
                    sums[name[:-len(suffix)]] = sum(v[1] for v in entry[1][2][1])
    return sums


def _masked_attributes(attributes):
    items = []
    for item in attributes[1]:
        if item[0] == 'tuple' and item[1][0][1] in _MASKED_ATTRIBUTES:
            item = ('tuple', [item[1][0], ('const', 'masked')])
        items.append(item)
    return (attributes[0], items)


def _masked_matrix(value):
    """The matrix datasets without the numbers the fixes change (the pixels
    and the coverage column), their lengths and names kept."""
    items = []
    for entry in value[1]:
        if entry[0] == 'tuple' and len(entry[1]) == 4:
            name, dtype, data, attributes = entry[1]
            if any(name[1].endswith(suffix) for suffix in _MASKED_DATA):
                data = ('const', 'masked')
            elif any(name[1].endswith(suffix) for suffix in _LENGTH_ONLY):
                data = ('int', len(data[1]))
            entry = ('tuple', [name, dtype, data, _masked_attributes(attributes)])
        elif entry[0] == 'tuple' and len(entry[1]) == 2:
            entry = ('tuple', [entry[1][0], _masked_attributes(entry[1][1])])
        items.append(entry)
    return (value[0], items)


def _sam_columns(text):
    columns = text.split('\t')
    del columns[8]  # the insert size, Q5
    return columns


def _check_bam(old, new):
    """The records of --outBam: the new ones are the 3.7.6 ones without the
    removed pairs and carry the insert size of the pair (Q5)."""
    if old is None or new is None or old[0] != 'list' or new[0] != 'list':
        return None if old == new else 'the --outBam records appeared or vanished'
    position = 0
    for text in new[1]:
        want = _sam_columns(text[1])
        while position < len(old[1]) and _sam_columns(old[1][position][1]) != want:
            position += 1
        if position == len(old[1]):
            return 'a --outBam record is not one of 3.7.6 apart from its insert size'
        position += 1
    for first, second in zip(new[1][0::2], new[1][1::2]):
        a, b = first[1].split('\t'), second[1].split('\t')
        if a[2] == b[2] and (int(a[8]) != int(b[3]) - int(a[3]) or int(b[8]) != int(a[3]) - int(b[3])):
            return 'the insert size of {} is not the distance of the mates (Q5)'.format(a[0])
    return None


def _masked(record, matrix_too):
    ret = list(record['return'][1])
    ret[2] = ('const', 'masked')
    if matrix_too and isinstance(ret[6], tuple) and ret[6][0] == 'dict':
        files = []
        for key, value in ret[6][1]:
            if key[1] == 'matrix' and value[0] == 'list':
                value = _masked_matrix(value)
            elif key[1] == 'bam':
                value = ('const', 'masked')
            files.append((key, value))
        ret[6] = ('dict', files)
    elif not matrix_too:
        ret[6] = ('const', 'masked')
    out = dict(record)
    out['return'] = ('tuple', ret)
    return _masked_logs(out)


def _masked_logs(record):
    out = dict(record)
    if 'logs' in record:
        out['logs'] = [(entry[0], entry[1], '<pairs> valid pairs added to matrix')
                       if entry[2].endswith('valid pairs added to matrix\n') else entry
                       for entry in record['logs']]
    return out


def _pair_counts_of_logs(old, new):
    for was, now in zip(old.get('logs', []), new.get('logs', [])):
        if was[2].endswith('valid pairs added to matrix\n'):
            if int(now[2].split()[0]) > int(was[2].split()[0]):
                return 'a progress record counts more pairs than 3.7.6'
    return None


def fixed_counts(region=False, matrix=True, prepare=None):
    """The deviation of the cases that classify read pairs (Q1 to Q6);
    `prepare(expected, actual)` returns the records to compare and a message
    (the provenance cases of Fix 3 strip their attributes first)."""
    def check(expected, actual):
        import apieq
        if prepare is not None:
            expected, actual, message = prepare(expected, actual)
            if message:
                return message
        if 'raise' in expected or 'raise' in actual:
            message = _pair_counts_of_logs(expected, actual)
            if message:
                return message
            return apieq.compare_records(_masked_logs(expected), _masked_logs(actual), _NoTolerance)
        old_ret, new_ret = expected['return'][1], actual['return'][1]
        removed = 0
        if old_ret[2][0] == 'str' and new_ret[2][0] == 'str':
            message, removed = _check_qc(old_ret[2][1], new_ret[2][1], region)
            if message:
                return message
        if matrix:
            old_files, new_files = _files(expected), _files(actual)
            if old_files is not None and new_files is not None:
                old_sums, new_sums = _sums(old_files.get('matrix')), _sums(new_files.get('matrix'))
                for prefix in old_sums:
                    lost = old_sums[prefix] - new_sums.get(prefix, 0)
                    if lost != removed:
                        return 'the pixels of {!r} lose {} counts, not the {} removed pairs'.format(prefix, lost, removed)
                message = _check_bam(old_files.get('bam'), new_files.get('bam'))
                if message:
                    return message
        message = _pair_counts_of_logs(expected, actual)
        if message:
            return message
        return apieq.compare_records(_masked(expected, matrix), _masked(actual, matrix), _NoTolerance)
    return check


def _provenance_stripped(resolutions, assembly=None):
    def prepare(expected, actual):
        message = _plain(actual, resolutions, assembly)
        if message:
            return None, None, message
        exp, act = dict(expected), dict(actual)
        exp['return'] = ('tuple', [_strip(x) for x in expected['return'][1]])
        act['return'] = ('tuple', [_strip(x) for x in actual['return'][1]])
        return exp, act, None
    return prepare


def _two_resolutions_prepared(expected, actual):
    """The checks of `two_resolutions` (both groups exist, plain provenance),
    without its comparison of the rest, which fixed_counts does."""
    message = _plain(actual, 2)
    if message:
        return None, None, message
    names = [item[1][0][1] for item in actual['return'][1][6][1][0][1][1]
             if isinstance(item, tuple) and item[0] == 'tuple' and item[1][0][0] == 'str']
    for resolution in ('resolutions/5000', 'resolutions/10000'):
        if resolution not in names:
            return None, None, 'no group {} in {}'.format(resolution, names)
    exp, act = dict(expected), dict(actual)
    exp['return'] = ('tuple', [x if i != 6 else ('const', 'None') for i, x in enumerate(expected['return'][1])])
    act['return'] = ('tuple', [x if i != 6 else ('const', 'None') for i, x in enumerate(actual['return'][1])])
    return exp, act, None


# Cases whose numbers do not change under Fix 3 alone keep that deviation
# and add the one of Q1 to Q6 on top.
_FIXES = 'Fixes Q1 to Q6 of hicBuildMatrix (see the module docstring)'
_AFFECTED = ('createMatrix/h5', 'createMatrix/rf_bins', 'createMatrix/two_enzymes_lowercase',
             'createMatrix/out_bam_sizes', 'createMatrix/buffers', 'createMatrix/buffers_parallel',
             'createMatrix/test_run', 'createMatrix/micro_c', 'createMatrix/small_test',
             'createMatrix/restriction_without_dangling', 'createMatrix/exact_multiple_buffers')
for _case in CASES:
    if _case.id == 'createMatrix/mcool':
        _case.deviation = fixed_counts(prepare=_provenance_stripped(3))
        _case.note = 'Fix 3 and ' + _FIXES
    elif _case.id == 'createMatrix/cool_assembly':
        _case.deviation = fixed_counts(prepare=_provenance_stripped(1, 'dm3'))
        _case.note = 'Fix 3 and ' + _FIXES
    elif _case.id == 'createMatrix/mcool_two_resolutions':
        _case.deviation = fixed_counts(matrix=False, prepare=_two_resolutions_prepared)
        _case.note = 'Fix 3 and ' + _FIXES
    elif _case.id == 'createMatrix/region':
        _case.deviation = fixed_counts(region=True)
        _case.note = _FIXES
    elif _case.id in _AFFECTED:
        _case.deviation = fixed_counts()
        _case.note = _FIXES

CASES.append(Case('createMatrix/mcool_do_test_run', M + 'createMatrix', args=(File(),),
                  kwargs={'out_name': 'matrix.mcool', 'bin_size': (5000, 10000, 20000), 'test_run': True},
                  driver='create'))
CASES[-1].deviation = do_test_run_returns
CASES[-1].note = 'Fix 3 of HiCExplorer 4 (project owner 2026-11-06), see the module docstring'
