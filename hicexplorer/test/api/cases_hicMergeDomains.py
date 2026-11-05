"""Equivalence cases for hicexplorer.hicMergeDomains (apieq.py).

The chain cases run the functions the way main() of 3.7.6 did on the real
domain files of the test data (10, 50 and 100 kb hicFindTADs domains, CTCF
peaks named like the domain files' chromosomes): read and filter each file,
merge, number, relate, write both lists and draw the relation trees. The
trees are drawn in the 'plain' and 'svg' formats, whose output dot writes
without a time stamp, so the files can be compared byte for byte. The other
cases call each function on its own, including the in place changes, the
aliasing of row objects that merge_list creates and the exceptions.

Deviations from 3.7.6, kept on purpose (cpp/python/src/api/hicMergeDomains.cpp):

* create_tree runs graphviz's dot with the command graphviz 0.20.3 runs, on
  the same DOT source, without the graphviz package; patching
  graphviz.Digraph does not reach it.
* merge_protein raises ValueError for a bin size that is not positive, where
  3.7.6 loops forever (no case can pin an endless loop).
* Row fields that are not str raise TypeError naming HiCExplorer 4.
"""
import contextlib
import io
import os

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'test_data', 'hicMergeDomains')
TEN = os.path.join(DATA, '10kbtad_domains.bed')
FIFTY = os.path.join(DATA, '50kbtad_domains.bed')
HUNDRED = os.path.join(DATA, '100kbtad_domains.bed')
CTCF = os.path.join(DATA, 'ctcf_sorted_nochr.bed')
M = 'hicexplorer.hicMergeDomains.'


def files_of(directory):
    out = {}
    for name in sorted(os.listdir(directory)):
        path = os.path.join(directory, name)
        if os.path.isfile(path):
            with open(path, 'rb') as handle:
                out[name] = handle.read()
    return out


def chain(target, directory, domain_files, protein, min_peak=1, value=5000, percent=0.5, fmt='plain'):
    """main() of 3.7.6 with every step a public function."""
    import hicexplorer.hicMergeDomains as m
    proteins = m.read_protein(protein) if protein else None
    lists = [m.create_list_with_protein(f, min_peak, proteins) for f in domain_files]
    merged = lists[0]
    for other in lists[1:]:
        merged = m.merge_list(merged, other, value)
    m.add_id(merged)
    m.write_in_file(merged, os.path.join(directory, 'merged.bed'))
    printed = io.StringIO()
    if len(domain_files) > 1:
        relations = m.create_relationsship_list(merged, percent)
        m.write_in_file(relations, os.path.join(directory, 'relations.txt'), True)
        with contextlib.redirect_stdout(printed):
            m.create_tree(relations, merged, os.path.join(directory, 'tree'), fmt)
    return {'files': files_of(directory), 'printed': printed.getvalue(), 'rows': len(merged)}


def aliasing(target, *paths, value=5000):
    """merge_list keeps the row objects: which result slots are the same object."""
    import hicexplorer.hicMergeDomains as m
    lists = [m.create_list_of_file(p, 0)[0][:300] for p in paths]
    merged = target(lists[0], lists[1], value)
    ids = {}
    pattern = [ids.setdefault(id(row), len(ids)) for row in merged]
    sources = [[(k, i) for k, lst in enumerate(lists) for i, r in enumerate(lst) if r is row][:1] for row in merged]
    return merged, pattern, sources


def reversed_order(target, path_a, path_b):
    """The 3.7.6 defect pinned by the tool's tests: merge_list appends one row
    object twice and add_id then numbers it twice."""
    import hicexplorer.hicMergeDomains as m
    a = m.create_list_of_file(path_a, 0)[0]
    b = m.create_list_of_file(path_b, 0)[0]
    merged = target(b, a, 5000)
    m.add_id(merged)
    ids = {}
    return [r[3] for r in merged][:50], [r[3] for r in merged][-50:], len(merged), \
        [ids.setdefault(id(row), len(ids)) for row in merged][-50:]


def protein_filter(target, domain_path, protein_path, min_peak, para=None):
    import hicexplorer.hicMergeDomains as m
    rows, bin_size = m.create_list_of_file(domain_path, 0)
    binned = m.merge_protein(m.read_protein(protein_path), bin_size, min_peak)
    before = len(rows)
    result = target(rows, binned) if para is None else target(rows, binned, para)
    return before, len(result), result is rows, result[:40], binned[0][:20]


def tree(target, directory, relations, domains, fmt):
    printed = io.StringIO()
    with contextlib.redirect_stdout(printed):
        result = target(relations, domains, os.path.join(directory, 'tree'), fmt)
    return result, printed.getvalue(), files_of(directory)


def write(target, directory, rows, name, *args):
    result = target(rows, os.path.join(directory, name), *args)
    return result, files_of(directory)


def in_place(target, *args):
    result = target(*args)
    return result, args[0]


def rows(text):
    return [line.split('\t') for line in text.strip('\n').split('\n')]


SMALL = ("1\t0\t100000\tID_1\t-1.5\t.\n1\t100000\t250000\tID_2\t-0.5\t.\n1\t250000\t300000\tID_3\t-0.45\t.\n"
         "1\t300000\t600000\tID_4\t-0.40\t.\n2\t0\t50000\tID_5\t0.1\t.\n2\t50000\t200000\tID_6\t0.2\t.\n")
OVERLAP = ("1\t0\t500000\tx\t1\n1\t0\t100000\tx\t1\n1\t80000\t300000\tx\t1\n1\t260000\t500000\tx\t1\n"
           "1\t400000\t520000\tx\t1\n2\t0\t90000\tx\t1\n2\t0\t40000\tx\t1\n2\t30000\t90000\tx\t1\n")
# start coordinates without a 0: create_list_of_file takes its bin size from
# them, and a start of 0 gives the bin size 0 on which merge_protein of 3.7.6
# never returns
SHIFTED = ("1\t20000\t100000\tID_1\t-1.5\t.\n1\t100000\t250000\tID_2\t-0.5\t.\n1\t250000\t300000\tID_3\t-0.45\t.\n"
           "1\t300000\t600000\tID_4\t-0.40\t.\n2\t10000\t50000\tID_5\t0.1\t.\n2\t50000\t200000\tID_6\t0.2\t.\n")
PEAKS = ("1\t5000\t5100\n1\t7000\t7100\n1\t260000\t260200\n1\t270000\t270100\n1\t650000\t650100\n"
         "2\t1000\t1100\n2\t120000\t120100\n")


def numbered(text):
    out = rows(text)
    for i, r in enumerate(out):
        r[3] = 'ID_{}'.format(i + 1)
    return out


CASES = [
    Case('chain_two_files', M + 'main', args=(File(), [TEN, FIFTY], None), driver='chain'),
    Case('chain_two_files_protein', M + 'main', args=(File(), [TEN, FIFTY], CTCF), driver='chain'),
    Case('chain_three_files_svg', M + 'main', args=(File(), [TEN, FIFTY, HUNDRED], CTCF, 3, 10000, 0.3, 'svg'),
         driver='chain'),
    Case('chain_one_file_protein', M + 'main', args=(File(), [FIFTY], CTCF, 2), driver='chain'),
    Case('create_list_of_file', M + 'create_list_of_file', args=(HUNDRED, 0)),
    Case('create_list_of_file_small', M + 'create_list_of_file', args=(File('a.bed'), []),
         files={'a.bed': SMALL.encode()}),
    Case('create_list_of_file_crlf_trailing_space', M + 'create_list_of_file', args=(File('a.bed'), 0),
         files={'a.bed': '1\t1200\t2000\tx \t1.0\t. \r\n1\t2000\t3500\tx\t2\t.\u3000\n'.encode()}),
    Case('create_list_of_file_no_tab', M + 'create_list_of_file', args=(File('a.bed'), 0),
         files={'a.bed': b'1 0 100\n'}),
    Case('create_list_of_file_not_a_number', M + 'create_list_of_file', args=(File('a.bed'), 0),
         files={'a.bed': b'1\tx10\t100\n'}),
    Case('create_list_of_file_missing', M + 'create_list_of_file', args=(File('missing.bed'), 0)),
    Case('merge_list_aliasing', M + 'merge_list', args=(TEN, FIFTY), driver='aliasing'),
    Case('merge_list_aliasing_value', M + 'merge_list', args=(FIFTY, HUNDRED), kwargs={'value': 60000},
         driver='aliasing'),
    Case('merge_list_reversed_order', M + 'merge_list', args=(TEN, FIFTY), driver='reversed_order'),
    Case('merge_list_empty_second', M + 'merge_list', args=(lambda: (rows(SMALL), [], 5000))),
    Case('add_id', M + 'add_id', args=lambda: (rows(SMALL),)),
    Case('add_id_short_row', M + 'add_id', args=lambda: ([['1', '2', '3']],)),
    Case('relations_small', M + 'create_relationsship_list', args=lambda: (numbered(OVERLAP),)),
    Case('relations_percent', M + 'create_relationsship_list', args=lambda: (numbered(OVERLAP), 0.05)),
    Case('add_relation_new', M + 'add_relation_to_list', args=lambda: ([], '1', 'ID_1', 'ID_2'), driver='in_place'),
    Case('add_relation_existing', M + 'add_relation_to_list',
         args=lambda: ([['1', 'ID_1', ['ID_2']], ['1', 'ID_4', ['ID_5']]], '1', 'ID_1', 'ID_3'), driver='in_place'),
    Case('add_relation_later', M + 'add_relation_to_list',
         args=lambda: ([['1', 'ID_1', ['ID_2']], ['1', 'ID_4', ['ID_5']]], '1', 'ID_7', 'ID_8'), driver='in_place'),
    Case('add_relation_out_of_range', M + 'add_relation_to_list',
         args=lambda: ([['1', 'ID_5', ['ID_6']], ['1', 'ID_4', ['ID_7']]], '1', 'ID_2', 'ID_3'), driver='in_place'),
    Case('write_domains', M + 'write_in_file', args=(File(), rows(SMALL), 'out.bed'), driver='write'),
    Case('write_relations', M + 'write_in_file',
         args=(File(), [['1', 'ID_1', ['ID_2', 'ID_3']], ['2', 'ID_5', ['ID_6']]], 'rel.txt', True), driver='write'),
    Case('write_empty_row', M + 'write_in_file', args=(File(), [['a', 'b'], []], 'out.bed'), driver='write'),
    Case('write_unwritable', M + 'write_in_file', args=(File(), rows(SMALL), 'missing/out.bed'), driver='write'),
    Case('tree_plain', M + 'create_tree',
         args=(File(), [['1', 'ID_1', ['ID_2', 'ID_3']], ['1', 'ID_3', ['ID_4']], ['2', 'ID_5', ['ID_6']]],
               numbered(SMALL), 'plain'), driver='tree'),
    Case('tree_dot_source', M + 'create_tree',
         args=(File(), [['1', 'ID_2', ['ID_3']], ['2', 'ID_5', ['ID_6']]], numbered(SMALL), 'DOT'), driver='tree'),
    Case('tree_unknown_format', M + 'create_tree',
         args=(File(), [['1', 'ID_1', ['ID_2']]], numbered(SMALL), 'PDF2'), driver='tree'),
    Case('tree_no_relations', M + 'create_tree', args=(File(), [], numbered(SMALL), 'plain'), driver='tree'),
    Case('small_list', M + 'create_small_list', args=lambda: (numbered(SMALL),)),
    Case('small_list_empty', M + 'create_small_list', args=([],)),
    Case('read_protein', M + 'read_protein', args=(File('p.bed'),), files={'p.bed': PEAKS.encode()}),
    Case('read_protein_real', M + 'read_protein', args=(CTCF,), driver=None, compare_args=True),
    Case('read_protein_empty', M + 'read_protein', args=(File('p.bed'),), files={'p.bed': b''}),
    Case('merge_protein', M + 'merge_protein', args=lambda: ([[r[0:3] for r in rows(PEAKS) if r[0] == '1'],
                                                             [r[0:3] for r in rows(PEAKS) if r[0] == '2']], 100000, 1)),
    Case('merge_protein_min_peak', M + 'merge_protein', args=lambda: ([[r[0:3] for r in rows(PEAKS)]], 100000, 2)),
    Case('compare_boundaries_real', M + 'compare_boundaries_protein', args=(FIFTY, CTCF, 1), driver='protein_filter'),
    Case('compare_boundaries_real_para', M + 'compare_boundaries_protein', args=(TEN, CTCF, 2, 0.5),
         driver='protein_filter'),
    Case('compare_boundaries_small', M + 'compare_boundaries_protein',
         args=lambda: (rows(SMALL), [[['1', 0, 100000, 2], ['1', 200000, 300000, 2]], [['2', 100000, 200000, 1]]]),
         driver='in_place'),
    Case('compare_boundaries_small_para', M + 'compare_boundaries_protein',
         args=lambda: (rows(SMALL), [[['1', 350000, 400000, 2]], [['2', 100000, 200000, 1]]], 0.9),
         driver='in_place'),
    Case('compare_boundaries_empty_chromosome', M + 'compare_boundaries_protein',
         args=lambda: (rows(SMALL), [[], [['2', 100000, 200000, 1]]]), driver='in_place'),
    Case('list_with_protein', M + 'create_list_with_protein', args=(File('a.bed'), 1, None),
         files={'a.bed': SMALL.encode()}),
    Case('list_with_protein_peaks', M + 'create_list_with_protein',
         args=lambda: (File('a.bed'), 1, [[r[0:3] for r in rows(PEAKS) if r[0] == '1'],
                                          [r[0:3] for r in rows(PEAKS) if r[0] == '2']]),
         files={'a.bed': SHIFTED.encode()}),
]
