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
