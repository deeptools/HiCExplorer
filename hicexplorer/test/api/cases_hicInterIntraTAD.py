"""Equivalence cases for hicexplorer.hicInterIntraTAD (apieq.py).

computeInterIntraTADs is the worker main() of 3.7.6 starts per share of a
chromosome's TADs. The worker cases build its arguments the way main() did
on the real GSM2644945 / GSM2644947 matrices of the test data, through cool
file names and through hiCMatrix objects loaded from the h5 files, and
compare what the worker put into its queue, with the type of every value:
the literal 0 of a missing neighbour, the numpy scalars of the sums and
ratios, the Python ints and floats of the counts and densities. On an h5
matrix the last TAD of a chromosome with three or more TADs had an empty
left block and 3.7.6 failed with ZeroDivisionError; HiCExplorer 4 gives it
the block (fix 2). For a failure only the 'Fail: <message>' part is compared.

HiCExplorer 4 cuts the left and right inter-TAD blocks from one TAD geometry for
cool and h5 (the project owner's decision of 2026-11-06): every TAD but the first
has a left block (previous TAD x this TAD) and every TAD but the last a right
block (this TAD x next TAD). 3.7.6 dropped the last bin of the next TAD from the
cool right block, gave the last TAD a stale left block (empty on h5, so that the
density divided by zero) and the second TAD of a two-TAD chromosome none. The
worker cases of these shapes therefore carry a `deviation` (see
test_api_equivalence.py) in place of the comparison with 3.7.6: the intra-TAD
lists and the TAD rows must still equal 3.7.6, the inter lists keep its types,
and the new blocks must tile (the right block of a TAD is the left block of its
successor, the same cells, so sums, counts and densities agree).

Deviations from 3.7.6, kept on purpose (cpp/python/src/api/hicInterIntraTAD.cpp):

* The traceback text after 'Fail: <message>' of a failed worker.
* readDomainBoundaries reads the file with hicx::TextTable, which raises
  pandas' ParserError also for a row with fewer fields than the first, where
  pandas fills the missing fields with NaN.
"""
import os

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'test_data', 'hicDifferentialTAD')
M = 'hicexplorer.hicInterIntraTAD.'
COOL = os.path.join(DATA, 'GSM2644947_Auxin2days-R1.100000_chr1.cool')
H5 = os.path.join(DATA, 'GSM2644945_Untreated-R1.100000_chr1.h5')
COOL_2 = os.path.join(HERE, '..', 'test_data', 'hicInterIntraTAD', 'GSM2644947_Auxin2days-R1.100000_chr1_chr2.cool')
DOMAINS = os.path.join(DATA, 'untreated_R1_domains.bed')
DOMAINS_2 = os.path.join(HERE, '..', 'test_data', 'hicInterIntraTAD', 'untreated_R1_domains_chr1_chr2.bed')
FOUR = os.path.join(DATA, 'domains_chr1_4.bed')


def frame(df):
    return (type(df).__name__, [str(d) for d in df.dtypes], list(df.columns), type(df.columns).__name__,
            type(df.index).__name__, list(df.index), df.values.tolist())


def read_frame(target, *args):
    return frame(target(*args))


def chromosomes(domains):
    out = []
    for row in domains:
        if out and out[-1][0][0] == row[0]:
            out[-1].append(row)
        else:
            out.append([row])
    return out


def worker(target, matrix, domains_file, cool, thread_id, chromosome=0, share=None):
    import queue
    import hicexplorer.hicInterIntraTAD as module
    domains = module.readDomainBoundaries(domains_file).values.tolist()
    rows = chromosomes(domains)[chromosome]
    if share is not None:
        rows = rows[share[0]:share[1]]
    if not cool:
        from hicmatrix import HiCMatrix as hm
        matrix = hm.hiCMatrix(matrix)
    q = queue.Queue()
    result = target(matrix, rows, cool, thread_id, q)
    data = q.get()
    if isinstance(data, str):
        data = data.split('Traceback')[0]
    return result, data


CASES = [
    Case('read_domains', M + 'readDomainBoundaries', args=(DOMAINS_2,), driver='read_frame'),
    Case('read_domains_five_columns', M + 'readDomainBoundaries', args=(File('d.bed'),),
         files={'d.bed': b"chr1\t0\t100\tID_1\t-0.5\n"}, driver='read_frame'),
    Case('worker_cool_all', M + 'computeInterIntraTADs', args=(COOL, DOMAINS, True, ''), driver='worker'),
    Case('worker_cool_first_share', M + 'computeInterIntraTADs', args=(COOL, DOMAINS, True, None),
         kwargs={'share': (0, 36)}, driver='worker'),
    Case('worker_cool_middle_share', M + 'computeInterIntraTADs', args=(COOL, DOMAINS, True, True),
         kwargs={'share': (34, 71)}, driver='worker'),
    Case('worker_cool_last_share', M + 'computeInterIntraTADs', args=(COOL, DOMAINS, True, False),
         kwargs={'share': (69, 140)}, driver='worker'),
    Case('worker_cool_four_tads', M + 'computeInterIntraTADs', args=(COOL, FOUR, True, ''), driver='worker'),
    Case('worker_cool_two_chromosomes', M + 'computeInterIntraTADs', args=(COOL_2, DOMAINS_2, True, ''),
         kwargs={'chromosome': 1}, driver='worker'),
    Case('worker_h5_zero_division', M + 'computeInterIntraTADs', args=(H5, DOMAINS, False, ''), driver='worker'),
    Case('worker_h5_first_share', M + 'computeInterIntraTADs', args=(H5, DOMAINS, False, None),
         kwargs={'share': (0, 60)}, driver='worker'),
    Case('worker_h5_two_tads', M + 'computeInterIntraTADs', args=(H5, DOMAINS, False, ''),
         kwargs={'share': (10, 12)}, driver='worker'),
]


# Fix 2 of HiCExplorer 4 (the TAD geometry). The 15 lists of what the worker
# puts into its queue: 0 left sum, 1 right sum, 2 left density, 3 right
# density, 4 left contacts, 5 right contacts, 6 left nnz, 7 right nnz,
# 8 intra sum, 9 intra contacts, 10 intra nnz, 11 intra density,
# 12 left/intra ratio, 13 right/intra ratio, 14 left/right/intra ratio,
# then the TAD rows.
_INTER = (0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14)
_TILES = ((0, 1), (2, 3), (4, 5), (6, 7))


def _lists(record):
    data = record['return'][1][1]
    return data[1] if data[0] == 'list' else None


def _shape(value):
    """A canonical value without its number: the kind stays."""
    if isinstance(value, tuple) and value and value[0] in ('float', 'npscalar'):
        return value[:2] if value[0] == 'npscalar' else ('float',)
    if isinstance(value, tuple) and value and value[0] == 'int':
        return ('int',)
    if isinstance(value, tuple) and len(value) == 2 and value[0] in ('list', 'tuple'):
        return (value[0], [_shape(v) for v in value[1]])
    return value


def _number(value):
    return value[2][1] if value[0] == 'npscalar' else value[1][1] if value[0] == 'float' else value[1]


def tads_tile(expected, actual):
    """The deviation of the worker cases: see the module docstring."""
    from apieq import compare
    new = _lists(actual)
    if new is None:
        return 'the worker failed: {!r}'.format(actual['return'])
    old = _lists(expected)
    if old is not None:
        if len(old) != len(new):
            return 'the worker put {} lists into its queue, 3.7.6 {}'.format(len(new), len(old))
        for index in range(len(new)):
            if index in _INTER:
                # same kinds of value, except where a block that 3.7.6 left out (the integer 0) now exists
                for position, (was, now) in enumerate(zip(old[index][1], new[index][1])):
                    if _shape(was) != _shape(now) and was != ('int', 0):
                        return 'list {}[{}] changed its type: {} against {}'.format(index, position, now, was)
            else:
                message = compare(old[index], new[index], 0.0, 0.0, 'list {}'.format(index))
                if message:
                    return message
    # the blocks tile: right block of TAD i == left block of TAD i + 1
    for left, right in _TILES:
        lefts, rights = new[left][1], new[right][1]
        for i in range(len(lefts) - 1):
            message = compare(rights[i], lefts[i + 1], 1e-12, 0.0, 'right of TAD {} against left of TAD {}'.format(i, i + 1))
            if message:
                return '{} (list {}/{})'.format(message, right, left)
    # no empty block: every left block that exists has contacts
    if any(c != ('int', 0) and _number(c) <= 0 for c in new[4][1]):
        return 'a left block has no contacts'
    return None


def _deviation(case_id, check):
    for case in CASES:
        if case.id == case_id:
            case.deviation = check
            case.note = 'Fix 2 of HiCExplorer 4 (TAD geometry), see the module docstring'
            return
    raise KeyError(case_id)


for _id in ('worker_cool_all', 'worker_cool_first_share', 'worker_cool_middle_share', 'worker_cool_last_share',
            'worker_cool_four_tads', 'worker_cool_two_chromosomes', 'worker_h5_zero_division', 'worker_h5_two_tads'):
    _deviation(_id, tads_tile)
