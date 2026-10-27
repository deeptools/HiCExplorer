"""Equivalence cases for hicexplorer.hicFindRestSite.find_pattern (apieq.py).

find_pattern writes into the file handle it is given and closes it; the
driver returns the return value, whether the handle was closed and the text
written.

Deviations from 3.7.6 (cpp/python/src/api/hicFindRestSite.cpp):
  * 3.7.6 left a NamedTemporaryFile behind in the temporary directory on every
    call; the C++ sorts in memory and leaves none.
  * patterns other than letters and '.' are matched with std::regex
    (ECMAScript), so a construct only Python's re has behaves differently and
    an invalid pattern raises ValueError instead of re.error. A pattern with a
    character class is mirrored by the reverse complement into an invalid one
    in 3.7.6 already (A[CG]T becomes A]CG[T), so only palindromic regular
    expressions are covered here.
  * a pattern that is not a str raises TypeError naming HiCExplorer 4.
"""
import os
import pathlib

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'test_data', 'hicFindRestSite')
CHRM = os.path.join(DATA, 'dm3_chrM.fasta')
CHRM_GZ = os.path.join(DATA, 'dm3_chrM.fasta.gz')
T = 'hicexplorer.hicFindRestSite.find_pattern'

DOCTEST_FASTA = b">chr1\nCTACGGTACGAACGTACGGTACGcgtaCGNAGTCATG\n"
TWO_RECORDS = (b">chrB description text\nGATCnnGATCgatc\r\nGA TC\n"
               b">chrA\nAAGCTTGATC\n>chr10\nGATC\n>chr2\nTTGATCAA\n")


def write_to(target, patterns, fasta, directory):
    out_path = os.path.join(directory, 'out.bed')
    handle = open(out_path, 'w')
    result = target(patterns, fasta, handle)
    closed = handle.closed
    with open(out_path) as written:
        return result, closed, written.read()


def fasta_handle(target, patterns, fasta, directory):
    """fasta_file given as an open file, whose name is used."""
    with open(fasta) as handle:
        return write_to(target, patterns, handle, directory)


CASES = [
    Case('real/hindIII', T, args=(['AAGCTT'], CHRM, File()), driver='write_to'),
    Case('real/dpnII', T, args=(['GATC'], CHRM, File()), driver='write_to'),
    Case('real/two_patterns_collide', T, args=(['AAGCTT', 'AAGC'], CHRM, File()), driver='write_to'),
    Case('real/two_patterns_other_order', T, args=(['AAGC', 'AAGCTT'], CHRM, File()), driver='write_to'),
    Case('real/gzip', T, args=(['GATC', 'AAGCTT'], CHRM_GZ, File()), driver='write_to'),
    Case('real/lowercase_pattern', T, args=(['aagctt'], CHRM, File()), driver='write_to'),
    Case('real/fasta_handle', T, args=(['GATC'], CHRM, File()), driver='fasta_handle'),
    Case('real/pathlib', T, args=lambda: (['GATC'], pathlib.Path(CHRM)), kwargs=None,
         driver='write_to_pathlib'),
    Case('doctest/palindrome', T, args=(['GTAC'], File('t.fa'), File()), files={'t.fa': DOCTEST_FASTA},
         driver='write_to'),
    Case('doctest/regexp', T, args=(['CG.AG'], File('t.fa'), File()), files={'t.fa': DOCTEST_FASTA},
         driver='write_to'),
    Case('regexp/quantifier', T, args=(['GA+TC'], File('t.fa'), File()), files={'t.fa': TWO_RECORDS},
         driver='write_to'),
    Case('records/sorted_by_name', T, args=(['GATC', 'TTGA'], File('t.fa'), File()), files={'t.fa': TWO_RECORDS},
         driver='write_to'),
    Case('patterns/str_is_iterated', T, args=('GA', File('t.fa'), File()), files={'t.fa': TWO_RECORDS},
         driver='write_to'),
    Case('patterns/tuple', T, args=(('AAGCTT', 'GATC'), File('t.fa'), File()), files={'t.fa': TWO_RECORDS},
         driver='write_to'),
    Case('patterns/empty', T, args=([], File('missing.fa'), File()), driver='write_to'),
    Case('patterns/no_match', T, args=(['CCCCCCCC'], File('t.fa'), File()), files={'t.fa': TWO_RECORDS},
         driver='write_to'),
    Case('fasta/empty', T, args=(['GATC'], File('t.fa'), File()), files={'t.fa': b''}, driver='write_to'),
    Case('fasta/missing', T, args=(['GATC'], File('missing.fa'), File()), driver='write_to'),
    Case('fasta/directory', T, args=(['GATC'], File(), File()), driver='write_to'),
]


def write_to_pathlib(target, patterns, fasta):
    import tempfile
    directory = tempfile.mkdtemp()
    try:
        return write_to(target, patterns, fasta, directory)
    finally:
        import shutil
        shutil.rmtree(directory)
