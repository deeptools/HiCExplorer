"""Equivalence cases for hicexplorer.hicInterIntraTAD (apieq.py).

computeInterIntraTADs is the worker main() of 3.7.6 starts per share of a
chromosome's TADs. The worker cases build its arguments the way main() did
on the real GSM2644945 / GSM2644947 matrices of the test data, through cool
file names and through hiCMatrix objects loaded from the h5 files, and
compare what the worker put into its queue, with the type of every value:
the literal 0 of a missing neighbour, the numpy scalars of the sums and
ratios, the Python ints and floats of the counts and densities. On an h5
matrix the last TAD of a chromosome with three or more TADs has an empty
left block and 3.7.6 fails with ZeroDivisionError (pinned); for a failure
only the 'Fail: <message>' part is compared.

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
