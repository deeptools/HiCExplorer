"""Equivalence cases for hicexplorer.hicDifferentialTAD (apieq.py).

computeDifferentialTADs is the worker main() of 3.7.6 starts per share of a
chromosome's TADs. The worker cases build its arguments the way main() did
(readDomainBoundaries, .values.tolist(), the TADs of one chromosome, the
share of a thread) on the real GSM2644945 / GSM2644947 matrices of the test
data, through both of main()'s routes: cool file names (one region loaded
per TAD) and hiCMatrix objects loaded from the h5 files. What the worker put
into its queue is compared; for a failure only the 'Fail: <message>' part,
because the traceback after it is Python's in 3.7.6 and names the C++
worker here.

HiCExplorer 4 cuts the left and right inter-TAD blocks from one TAD geometry for
cool and h5 (the project owner's decision of 2026-11-06, the same as in
hicInterIntraTAD): every TAD but the first has a left block (previous TAD x this
TAD) and every TAD but the last a right block (this TAD x next TAD). 3.7.6 dropped
the last bin of the next TAD from the cool right block, gave the last TAD a stale
left block (NaN on h5) and the second TAD of a two-TAD chromosome none. The worker
cases of these shapes carry a `deviation` (see test_api_equivalence.py) in place
of the comparison with 3.7.6: the intra-TAD columns and the TAD rows must still
equal 3.7.6, and the new blocks must tile (the right test of a TAD is the left test
of its successor on the same cells, so statistic and p-value agree).

Deviations from 3.7.6, kept on purpose (cpp/python/src/api/hicDifferentialTAD.cpp):

* The traceback text after 'Fail: <message>' of a failed worker.
* readDomainBoundaries reads the file with hicx::TextTable, which raises
  pandas' ParserError also for a row with fewer fields than the first, where
  pandas fills the missing fields with NaN.
"""
import os

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'test_data', 'hicDifferentialTAD')
M = 'hicexplorer.hicDifferentialTAD.'
TARGET = os.path.join(DATA, 'GSM2644945_Untreated-R1.100000_chr1.cool')
CONTROL = os.path.join(DATA, 'GSM2644947_Auxin2days-R1.100000_chr1.cool')
TARGET_H5 = os.path.join(DATA, 'GSM2644945_Untreated-R1.100000_chr1.h5')
CONTROL_H5 = os.path.join(DATA, 'GSM2644947_Auxin2days-R1.100000_chr1.h5')
TARGET_2 = os.path.join(DATA, 'GSM2644945_Untreated-R1.100000_chr1_chr2.cool')
CONTROL_2 = os.path.join(DATA, 'GSM2644947_Auxin2days-R1.100000_chr1_chr2.cool')
DOMAINS = os.path.join(DATA, 'untreated_R1_domains.bed')
DOMAINS_2 = os.path.join(DATA, 'untreated_R1_domains_chr1_chr2.bed')
FOUR = os.path.join(DATA, 'domains_chr1_4.bed')


def frame(df):
    """A DataFrame as plain data: its dtypes, labels and rows."""
    return (type(df).__name__, [str(d) for d in df.dtypes], list(df.columns), type(df.columns).__name__,
            type(df.index).__name__, list(df.index), df.values.tolist())


def read_frame(target, *args):
    return frame(target(*args))


def chromosomes(domains):
    """The per chromosome lists of main()."""
    out = []
    for row in domains:
        if out and out[-1][0][0] == row[0]:
            out[-1].append(row)
        else:
            out.append([row])
    return out


def worker(target, target_matrix, control_matrix, domains_file, cool, p_value, thread_id, chromosome=0,
           share=None):
    """One worker of main(): the TADs of one chromosome, or a share of them."""
    import queue
    import hicexplorer.hicDifferentialTAD as module
    domains = module.readDomainBoundaries(domains_file).values.tolist()
    rows = chromosomes(domains)[chromosome]
    if share is not None:
        rows = rows[share[0]:share[1]]
    if not cool:
        from hicmatrix import HiCMatrix as hm
        target_matrix = hm.hiCMatrix(target_matrix)
        control_matrix = hm.hiCMatrix(control_matrix)
    q = queue.Queue()
    result = target(target_matrix, control_matrix, rows, cool, p_value, thread_id, q)
    data = q.get()
    if isinstance(data, str):
        data = data.split('Traceback')[0]
    return result, data


TABLE = (b"chr1\t0\t100\tID_1\t-0.5\t.\nchr1\t100\t250\tID_2\t0.25\t.\n"
         b"chr2\t0\t80\tID_3\tnan\t.\nchr2\t80\t90\t4\t1e-3\t.\n")

CASES = [
    Case('read_domains', M + 'readDomainBoundaries', args=(DOMAINS,), driver='read_frame'),
    Case('read_domains_types', M + 'readDomainBoundaries', args=(File('d.bed'),), files={'d.bed': TABLE},
         driver='read_frame'),
    Case('read_domains_nine_columns', M + 'readDomainBoundaries', args=(FOUR,), driver='read_frame'),
    Case('read_domains_five_columns', M + 'readDomainBoundaries', args=(File('d.bed'),),
         files={'d.bed': b"chr1\t0\t100\tID_1\t-0.5\n"}, driver='read_frame'),
    Case('read_domains_three_columns', M + 'readDomainBoundaries', args=(File('d.bed'),),
         files={'d.bed': b"chr1\t0\t100\n"}, driver='read_frame'),
    Case('read_domains_more_fields', M + 'readDomainBoundaries', args=(File('d.bed'),),
         files={'d.bed': TABLE + b"chr3\t0\t5\tx\t1\t.\textra\n"}, driver='read_frame'),
    Case('read_domains_empty', M + 'readDomainBoundaries', args=(File('d.bed'),), files={'d.bed': b""},
         driver='read_frame'),
    Case('read_domains_missing', M + 'readDomainBoundaries', args=(File('missing.bed'),), driver='read_frame'),
    Case('worker_cool_all', M + 'computeDifferentialTADs', args=(TARGET, CONTROL, DOMAINS, True, 0.05, ''),
         driver='worker'),
    Case('worker_cool_first_share', M + 'computeDifferentialTADs',
         args=(TARGET, CONTROL, DOMAINS, True, 0.05, None), kwargs={'share': (0, 36)}, driver='worker'),
    Case('worker_cool_middle_share', M + 'computeDifferentialTADs',
         args=(TARGET, CONTROL, DOMAINS, True, 0.01, True), kwargs={'share': (34, 71)}, driver='worker'),
    Case('worker_cool_last_share', M + 'computeDifferentialTADs',
         args=(TARGET, CONTROL, DOMAINS, True, 0.05, False), kwargs={'share': (69, 140)}, driver='worker'),
    Case('worker_cool_four_tads', M + 'computeDifferentialTADs', args=(TARGET, CONTROL, FOUR, True, 0.05, ''),
         driver='worker'),
    Case('worker_cool_second_chromosome', M + 'computeDifferentialTADs',
         args=(TARGET_2, CONTROL_2, DOMAINS_2, True, 0.05, ''), kwargs={'chromosome': 1}, driver='worker'),
    Case('worker_h5_all', M + 'computeDifferentialTADs', args=(TARGET_H5, CONTROL_H5, DOMAINS, False, 0.05, ''),
         driver='worker'),
    Case('worker_h5_first_share', M + 'computeDifferentialTADs',
         args=(TARGET_H5, CONTROL_H5, DOMAINS, False, 0.05, None), kwargs={'share': (0, 50)}, driver='worker'),
    Case('worker_cool_unknown_chromosome', M + 'computeDifferentialTADs',
         args=(TARGET, CONTROL, DOMAINS_2, True, 0.05, ''), kwargs={'chromosome': 1}, driver='worker'),
]


# Cases whose floats differ from 3.7.6 in the last bits since the special
# functions (normal tail, incomplete beta and gamma functions, log-gamma) are
# Boost.Math's and no longer the routines scipy uses: the p-values agree
# to a relative 2e-14 or better. The cases compare every float at rtol 1e-12;
# integers, strings, coordinates and the structure of the result are still
# compared exactly, so a changed call, a moved TAD or a different row count
# fails.
_BOOST_LAST_BITS = {
    'worker_cool_all',
    'worker_cool_first_share',
    'worker_cool_four_tads',
    'worker_cool_last_share',
    'worker_cool_middle_share',
    'worker_cool_second_chromosome',
    'worker_h5_all',
    'worker_h5_first_share',
}
for _case in CASES:
    if _case.id in _BOOST_LAST_BITS:
        _case.rtol = 1e-12
        _case.note = (_case.note + '; ' if _case.note else '') + (
            'Boost.Math special functions round the last bits of the p-values differently from scipy 1.14')


# Fix 2 of HiCExplorer 4 (the TAD geometry). The six lists of what the worker
# puts into its queue: 0 the statistics and 1 the p-values, one row per TAD
# with [left, right, intra]; 2, 3, 4 the rejection flags of the left, right
# and intra test; 5 the TAD rows.
def _lists(record):
    data = record['return'][1][1]
    return data[1] if data[0] == 'list' else None


def tads_tile(expected, actual):
    from apieq import compare
    new = _lists(actual)
    if new is None:
        return 'the worker failed: {!r}'.format(actual['return'])
    old = _lists(expected)
    if old is not None:
        if len(old) != len(new):
            return 'the worker put {} lists into its queue, 3.7.6 {}'.format(len(new), len(old))
        for index in (4, 5):
            message = compare(old[index], new[index], 1e-12, 0.0, 'list {}'.format(index))
            if message:
                return message
        for index in (0, 1):
            if len(old[index][1]) != len(new[index][1]):
                return 'list {} has another number of TADs'.format(index)
            for row, (was, now) in enumerate(zip(old[index][1], new[index][1])):
                message = compare(was[1][2], now[1][2], 1e-12, 0.0, 'intra of TAD {} in list {}'.format(row, index))
                if message:
                    return message
    for index in (0, 1):
        rows = new[index][1]
        for i in range(len(rows) - 1):
            message = compare(rows[i][1][1], rows[i + 1][1][0], 1e-12, 0.0,
                              'right test of TAD {} against left test of TAD {} (list {})'.format(i, i + 1, index))
            if message:
                return message
    return None


for _case in CASES:
    if _case.id in _BOOST_LAST_BITS:
        _case.deviation = tads_tile
        _case.note = (_case.note + '; ' if _case.note else '') + 'Fix 2 of HiCExplorer 4 (TAD geometry), see the module docstring'
