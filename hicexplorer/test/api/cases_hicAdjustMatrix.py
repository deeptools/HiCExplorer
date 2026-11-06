"""Equivalence cases for hicexplorer.hicAdjustMatrix (apieq.py).

adjustMatrix gets the argparse namespace the tool builds, on the matrices
and region files of the hicAdjustMatrix tests (small_test_matrix_50kb_res
as h5 and cool, gm12878_1_2_3.cool with weights). It returns a hiCMatrix,
which has no value representation; the drivers return its state (matrix,
bin table, nan_bins, correction factors, the masked bin bookkeeping,
chrBinBoundaries, interval trees and the file handler's state that save()
reads) and, for 'saved', the state of the matrix written by save() and
loaded again. exit(1) is compared as the SystemExit code.

Cases of defects that HiCExplorer 4 fixes on purpose (the same fixes as the
command line tool) carry a `deviation`, a callable(expected, actual) that
replaces the comparison with the 3.7.6 record and checks what the fix leaves
unchanged (see test_api_equivalence.py): --maskBadRegions masks the bins of
its BED file, an unknown chromosome name and a call without a selection exit
1, --action remove of --regions removes only the bins of the regions, and the
correction factors of an h5 matrix are kept and follow their bins.

Deviation (hic_object.hpp): the cool metadata of the file handler holds
str() values, so only the three entries save() reads are compared.
"""
import os

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
A = 'hicexplorer.hicAdjustMatrix.adjustMatrix'
H5 = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5')
COOL = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.cool')
GM = os.path.join(TEST_DATA, 'hicAdjustMatrix', 'gm12878_1_2_3.cool')
REGIONS = os.path.join(TEST_DATA, 'regions.bed')
REGIONS_XFAIL = os.path.join(TEST_DATA, 'regions_xfail.bed')
MASK = os.path.join(TEST_DATA, 'hicAdjustMatrix', 'mask.bed')
LI = os.path.join(TEST_DATA, 'Li_et_al_2015.h5')


def namespace(matrix, chromosomes=None, regions=None, action='keep', maskBadRegions=None, interIntraHandling=None):
    def make():
        import argparse
        return (argparse.Namespace(matrix=matrix, chromosomes=chromosomes, regions=regions, action=action,
                                   maskBadRegions=maskBadRegions, interIntraHandling=interIntraHandling,
                                   outFileName='unused.cool'),)
    return make


def state(hic):
    if hic is None:
        return None
    trees = getattr(hic, 'interval_trees', None)
    boundaries = getattr(hic, 'chrBinBoundaries', None)
    handler = hic.matrixFileHandler
    handler_state = None
    if handler is not None:
        f = handler.matrixFile
        handler_state = (type(f).__name__, getattr(f, 'matrixFileName', None), getattr(f, 'chrnameList', None),
                         getattr(f, 'correctionOperator', None), getattr(f, 'hic2cool_version', None),
                         getattr(f, 'hicmatrix_version', None),
                         [str((getattr(f, 'hic_metadata', None) or {}).get(k)) for k in
                          ('genome-assembly', 'matrix-generated-by', 'matrix-generated-by-url')])
    return {
        'matrix': hic.matrix,
        'cut_intervals': hic.cut_intervals,
        'nan_bins': hic.nan_bins,
        'correction_factors': hic.correction_factors,
        'distance_counts': hic.distance_counts,
        'orig_bin_ids': hic.orig_bin_ids,
        'orig_cut_intervals': hic.orig_cut_intervals,
        'chrBinBoundaries': None if boundaries is None else list(boundaries.items()),
        'interval_trees': None if trees is None else
        [(name, sorted((iv.begin, iv.end, iv.data) for iv in tree)) for name, tree in trees.items()],
        'handler': handler_state,
        'bin_size': hic.bin_size,
    }


def run(target, args):
    try:
        return state(target(args))
    except SystemExit as error:
        return 'SystemExit', error.code


def run_with_regions(target, regions, matrix, action):
    import argparse
    return run(target, argparse.Namespace(matrix=matrix, chromosomes=None, regions=regions, action=action,
                                          maskBadRegions=None, interIntraHandling=None))


def run_with_mask_bad_regions(target, bed, matrix, chromosomes=None):
    import argparse
    return run(target, argparse.Namespace(matrix=matrix, chromosomes=chromosomes, regions=None, action='keep',
                                          maskBadRegions=bed, interIntraHandling=None))


def saved(target, args, out):
    from hicmatrix import HiCMatrix as hm
    hic = target(args)
    hic.save(out)
    written = out if os.path.exists(out) else out + '.h5'
    reloaded = hm.hiCMatrix(written)
    result = state(reloaded)
    result.pop('handler')
    return os.path.basename(written), result


R = dict(driver='run', compare_args=False)
CASES = []
for name, path in (('h5', H5), ('cool', COOL)):
    for action in ('keep', 'remove', 'mask'):
        CASES += [
            Case('{}/chromosomes_{}'.format(name, action), A, args=namespace(path, ['chrX', 'chr3R'], action=action), **R),
            Case('{}/one_chromosome_{}'.format(name, action), A, args=namespace(path, ['chr3R'], action=action), **R),
            Case('{}/regions_{}'.format(name, action), A, args=namespace(path, regions=REGIONS, action=action), **R),
            Case('{}/regions_xfail_{}'.format(name, action), A, args=namespace(path, regions=REGIONS_XFAIL,
                                                                               action=action), **R),
            Case('{}/unknown_chromosome_{}'.format(name, action), A, args=namespace(path, ['chr10'], action=action),
                 **R),
            Case('{}/unknown_and_known_{}'.format(name, action), A,
                 args=namespace(path, ['chr10', 'chrX'], action=action), **R),
        ]
CASES += [
    Case('gm/keep_1', A, args=namespace(GM, ['1']), **R),
    Case('gm/keep_2_3', A, args=namespace(GM, ['3', '2']), **R),
    Case('gm/mask_2', A, args=namespace(GM, ['2'], action='mask'), **R),
    Case('gm/inter', A, args=namespace(GM, ['1', '2', '3'], interIntraHandling='inter'), **R),
    Case('gm/intra', A, args=namespace(GM, ['1', '2', '3'], interIntraHandling='intra'), **R),
    Case('gm/invalid_handling', A, args=namespace(GM, ['1', '2', '3'], interIntraHandling='both'), **R),
    Case('gm/mask_bad_regions', A, args=(File('bad.bed'), GM, []), driver='run_with_mask_bad_regions',
         compare_args=False, files={'bad.bed': b'1\t0\t5000000\n2\t10000000\t12000000\n'}),
    Case('gm/mask_bad_regions_no_chromosomes', A, args=(File('bad.bed'), GM), driver='run_with_mask_bad_regions',
         compare_args=False, files={'bad.bed': b'1\t0\t5000000\n2\t10000000\t12000000\n'}),
    Case('h5/mask_bad_regions', A, args=(MASK, H5), driver='run_with_mask_bad_regions', compare_args=False),
    Case('cool/mask_bad_regions', A, args=(MASK, COOL), driver='run_with_mask_bad_regions', compare_args=False),
    Case('li/keep_regions_with_factors', A, args=(File('part.bed'), LI, 'keep'), driver='run_with_regions',
         compare_args=False, files={'part.bed': b'X\t1000000\t3000000\n'}),
    Case('li/remove_regions_with_factors', A, args=(File('part.bed'), LI, 'remove'), driver='run_with_regions',
         compare_args=False, files={'part.bed': b'X\t1000000\t3000000\n'}),
    Case('h5/inter', A, args=namespace(H5, ['chrX', 'chr2L'], interIntraHandling='inter'), **R),
    Case('cool/keep_one_unknown', A, args=namespace(COOL, ['chrZ']), **R),
    Case('both_given', A, args=namespace(H5, ['chrX'], regions=REGIONS), **R),
    Case('nothing_given', A, args=namespace(H5), **R),
    Case('nothing_given_handling', A, args=namespace(H5, interIntraHandling='inter'), **R),
    Case('short_region_lines', A, args=(File('r.bed'), H5, 'remove'), driver='run_with_regions', compare_args=False,
         files={'r.bed': b'chrX\t10000000\t12000000\n\nab\nchr3R\t0\t1000000\tx\n'}),
    Case('region_lines_mask', A, args=(File('r.bed'), COOL, 'mask'), driver='run_with_regions', compare_args=False,
         files={'r.bed': b'chrX\t10000000\t12000000\nchr2L\t500000\t900000\n'}),
    Case('bad_region_int', A, args=(File('r.bed'), H5, 'keep'), driver='run_with_regions', compare_args=False,
         files={'r.bed': b'chrX\tten\t12000000\n'}),
    Case('missing_region_file', A, args=(File('missing.bed'), H5, 'keep'), driver='run_with_regions',
         compare_args=False),
    Case('saved/h5_to_cool_name', A, args=lambda: namespace(H5, ['chrX'])() + (File('out.cool'),),
         driver='saved', compare_args=False),
    Case('saved/cool_to_cool', A, args=lambda: namespace(COOL, ['chrX', 'chr2L'], action='remove')() +
         (File('out.cool'),), driver='saved', compare_args=False),
    Case('saved/gm_weights', A, args=lambda: namespace(GM, ['2'])() + (File('out.cool'),), driver='saved',
         compare_args=False),
    Case('saved/gm_to_h5', A, args=lambda: namespace(GM, ['1', '3'])() + (File('out.h5'),), driver='saved',
         compare_args=False),
]


# ---------------------------------------------------------------------------
# deviations (fixes of HiCExplorer 4)

def _fields(record):
    """The state of a returned matrix as {name: canonical value}, or None."""
    value = record.get('return')
    if not (isinstance(value, tuple) and len(value) == 2 and value[0] == 'dict'):
        return None
    return {key[1]: item for key, item in value[1]}


def _items(canonical):
    """The elements of a canonical list/tuple/array, or None."""
    if isinstance(canonical, tuple) and len(canonical) == 2 and isinstance(canonical[1], list):
        return canonical[1]
    return None


def _exit_one(record):
    from apieq import canon
    return 'return' in record and record['return'] == canon(('SystemExit', 1), ())


def _logged(record, level, text):
    return any(lvl == level and text in message for _, lvl, message in record['logs'])


def exits_one(text):
    """The call exits 1 and logs an error that contains `text`."""
    def check(expected, actual):
        if _exit_one(expected) and _exit_one(actual) and not text:
            return None
        if not _exit_one(actual):
            return 'expected exit(1), got {!r}'.format(actual.get('return', actual.get('raise')))
        if not _logged(actual, 'ERROR', text):
            return 'no error log containing {!r}: {!r}'.format(text, actual['logs'])
        return None
    return check


def _same_apart(expected, actual, ignore):
    """Everything but the fields in `ignore` and the discounted-regions
    warnings is the record of 3.7.6."""
    from apieq import compare
    if 'raise' in actual or 'raise' in expected or _fields(actual) is None or _fields(expected) is None:
        return 'unexpected record: {!r} against {!r}'.format(actual.get('raise', actual.get('return')),
                                                           expected.get('raise', expected.get('return')))
    new, old = _fields(actual), _fields(expected)
    if sorted(new) != sorted(old):
        return 'the state has other fields: {} against {}'.format(sorted(new), sorted(old))
    for key in sorted(old):
        if key in ignore:
            continue
        message = compare(old[key], new[key], 0.0, 0.0, key)
        if message:
            return message
    keep = lambda logs: [r for r in logs if 'discounted regions' not in r[2]]  # noqa: E731
    return compare(keep(expected['logs']), keep(actual['logs']), 0.0, 0.0, 'logs')


def removes_only_the_regions(expected, actual):
    """--regions --action remove: the bins of 3.7.6 are all in the result, the
    other bins it dropped are exactly the NaN bins, the logs are those of
    3.7.6 apart from the chromosome in the discounted-regions warning."""
    if _exit_one(expected):
        return exits_one('')(expected, actual)
    if 'raise' in expected:
        # 3.7.6 looked up the chromosome of the last BED line (a leaked loop
        # variable) and raised when the matrix does not have it
        return None if _fields(actual) is not None else 'the call failed: {!r}'.format(actual.get('raise'))
    message = _same_apart(expected, actual, ('matrix', 'cut_intervals', 'nan_bins', 'chrBinBoundaries',
                                             'interval_trees', 'orig_bin_ids', 'orig_cut_intervals',
                                             'correction_factors', 'distance_counts'))
    if message:
        return message
    old = [repr(i) for i in _items(_fields(expected)['cut_intervals'])]
    new = [repr(i) for i in _items(_fields(actual)['cut_intervals'])]
    position = 0
    for interval in old:
        while position < len(new) and new[position] != interval:
            position += 1
        if position == len(new):
            return 'a bin that 3.7.6 kept is missing from the result'
        position += 1
    if len(new) <= len(old):
        return 'the result is not larger than the 3.7.6 one although the matrix has NaN bins'
    return None


def masks_the_regions(expected, actual):
    """--maskBadRegions: the bins of the file are masked; 3.7.6 returned the
    matrix unchanged (or raised for a cool matrix)."""
    new = _fields(actual)
    if new is None:
        return 'the call failed: {!r}'.format(actual.get('raise', actual.get('return')))
    shown = _items(new['cut_intervals'])
    original = _items(new['orig_cut_intervals'])
    if not shown or not original or len(shown) >= len(original):
        return 'no bin was masked'
    if 'raise' not in expected and _fields(expected) is not None:
        before = _items(_fields(expected)['cut_intervals'])
        if len(before) != len(original):
            return 'the unmasked bin table is not the one of 3.7.6'
        if _fields(expected)['handler'] != new['handler']:
            return 'the file handler differs from 3.7.6'
    return None


def keeps_the_factors(expected, actual):
    """An h5 matrix keeps /correction_factors in correction_factors (3.7.6 left
    them in distance_counts), and the factors follow their bins."""
    new = _fields(actual)
    if new is None:
        return 'the call failed: {!r}'.format(actual.get('raise', actual.get('return')))
    factors = new['correction_factors']
    rows = len(_items(new['cut_intervals']))
    if not (isinstance(factors, tuple) and factors[0] == 'ndarray'):
        return 'the correction factors are gone: {!r}'.format(factors)
    if factors[3] != (rows,):
        return 'the factors do not match the {} bins: {!r}'.format(rows, factors[3])
    if new['distance_counts'] != ('const', 'None'):
        return 'distance_counts is not None: it still holds the factors'
    return None


def _deviation(case_id, check, note):
    for case in CASES:
        if case.id == case_id:
            case.deviation = check
            case.note = note
            return
    raise KeyError(case_id)


for _name in ('h5', 'cool'):
    for _action in ('keep', 'remove', 'mask'):
        _deviation('{}/unknown_chromosome_{}'.format(_name, _action), exits_one('chr10'),
                   'The unknown chromosome is reported as an error (3.7.6: a warning and then the same error).')
        _deviation('{}/unknown_and_known_{}'.format(_name, _action), exits_one('chr10'),
                   'An unknown chromosome name exits 1 (HiCExplorer 4); 3.7.6 only warned.')
    _deviation('{}/regions_remove'.format(_name), removes_only_the_regions,
               'Fixed in HiCExplorer 4: --action remove removes only the bins of the regions.')
    _deviation('{}/regions_xfail_remove'.format(_name), removes_only_the_regions,
               'Fixed in HiCExplorer 4: --action remove removes only the bins of the regions.')
_deviation('short_region_lines', removes_only_the_regions,
           'Fixed in HiCExplorer 4: --action remove removes only the bins of the regions; the warning names '
           'the chromosome of its region.')
for _id in ('gm/mask_bad_regions', 'gm/mask_bad_regions_no_chromosomes', 'h5/mask_bad_regions',
            'cool/mask_bad_regions'):
    _deviation(_id, masks_the_regions, 'Fixed in HiCExplorer 4: --maskBadRegions masks the bins of its BED file.')
_deviation('nothing_given', exits_one('No data to adjust'),
           'Fixed in HiCExplorer 4: no selection exits 1 (3.7.6 logged and returned None).')
_deviation('nothing_given_handling', exits_one('No data to adjust'),
           'Fixed in HiCExplorer 4: no selection exits 1 (3.7.6 raised an AttributeError).')
for _id in ('li/keep_regions_with_factors', 'li/remove_regions_with_factors'):
    _deviation(_id, keeps_the_factors,
               'Fixed in HiCExplorer 4: the factors of an h5 matrix stay in correction_factors and follow '
               'their bins (3.7.6 left them in distance_counts and did not reorder them).')
