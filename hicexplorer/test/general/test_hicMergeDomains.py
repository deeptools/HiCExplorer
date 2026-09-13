import os
import sys
import hashlib
import shutil
from tempfile import NamedTemporaryFile
from tempfile import mkdtemp
from psutil import virtual_memory
import subprocess
import pytest
import logging
log = logging.getLogger(__name__)

import graphviz
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from matplotlib.testing.exceptions import ImageComparisonFailure
from hicexplorer import hicMergeDomains
from hicexplorer.test.test_compute_function import compute

mem = virtual_memory()
memory = mem.total / 2**30

# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 4
HIGH_MEMORY = 120

REMOVE_OUTPUT = True

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")

# The relation trees are rendered by graphviz's `dot` executable. Without it on
# PATH, graphviz raises ExecutableNotFound from inside compute(), whose generic
# Exception does not match xfail(raises=ImageComparisonFailure), so the two tests
# below were reported as failures of the tool. They are skipped instead.
DOT_MISSING = shutil.which('dot') is None


def are_files_equal(file1, file2, delta=None):
    equal = True
    if delta:
        mismatches = 0
    with open(file1) as textfile1, open(file2) as textfile2:
        for x, y in zip(textfile1, textfile2):
            if x.startswith('File'):
                continue
            if x != y:
                if delta:
                    mismatches += 1
                    if mismatches > delta:
                        equal = False
                        break
                else:
                    equal = False
                    break
    return equal


@pytest.mark.xfail(raises=SystemExit)
def test_main_one_file():
    outfile_domains = NamedTemporaryFile(suffix='.txt', delete=True)
    args = "-d {} -om {} ".format(
        ROOT + 'hicMergeDomains/10kbtad_domains.bed',
        outfile_domains.name).split()
    # hicMergeDomains.main(args)
    compute(hicMergeDomains.main, args, 5)


def test_main_one_file_protein():
    outfile_domains = NamedTemporaryFile(suffix='.txt', delete=True)

    args = "-d {} -om {} -p {} ".format(
        ROOT + 'hicMergeDomains/10kbtad_domains.bed',
        outfile_domains.name, ROOT + 'hicMergeDomains/ctcf_sorted.bed').split()
    # hicMergeDomains.main(args)
    compute(hicMergeDomains.main, args, 5)

    assert are_files_equal(outfile_domains.name, ROOT + 'hicMergeDomains/one_file', delta=2)


@pytest.mark.skipif(DOT_MISSING, reason='graphviz dot executable not on PATH')
@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
def test_main_two_file_protein():
    outfile_domains = NamedTemporaryFile(suffix='.txt', delete=True)
    outfile_ctcf_relation = NamedTemporaryFile(suffix='.txt', delete=True)
    plot_folder = mkdtemp(prefix="plot_relations")

    args = "-d {} {} -om {} -or {} -ot {} -of {} -p {} ".format(
        ROOT + 'hicMergeDomains/10kbtad_domains.bed', ROOT + 'hicMergeDomains/50kbtad_domains.bed',
        outfile_domains.name, outfile_ctcf_relation.name, plot_folder + '/two_files_plot_ctcf', 'png', ROOT + 'hicMergeDomains/ctcf_sorted.bed').split()
    # hicMergeDomains.main(args)
    compute(hicMergeDomains.main, args, 5)
    assert are_files_equal(outfile_domains.name, ROOT + 'hicMergeDomains/two_files_ctcf.bed', delta=2)
    assert are_files_equal(outfile_ctcf_relation.name, ROOT + 'hicMergeDomains/two_files_relation_ctcf.txt', delta=2)

    list_data = list(range(1, 23, 1))
    list_data = ["{}".format(x) for x in list_data]
    list_data.append('X')
    for i in list_data:
        res = compare_images(ROOT + '/hicMergeDomains/ctcf_plot/two_files_plot_ctcf_' + i + '.png', plot_folder + '/two_files_plot_ctcf_' + i + '.png', tol=40)
        assert res is None, res


@pytest.mark.skipif(DOT_MISSING, reason='graphviz dot executable not on PATH')
@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
def test_main_two_file_no_protein():
    outfile_domains = NamedTemporaryFile(suffix='.txt', delete=True)
    outfile_ctcf_relation = NamedTemporaryFile(suffix='.txt', delete=True)
    plot_folder = mkdtemp(prefix="plot_relations")

    args = "-d {} {} -om {} -or {} -ot {} -of {}".format(
        ROOT + 'hicMergeDomains/10kbtad_domains.bed', ROOT + 'hicMergeDomains/50kbtad_domains.bed',
        outfile_domains.name, outfile_ctcf_relation.name, plot_folder + '/two_files_plot', 'png').split()
    # hicMergeDomains.main(args)
    compute(hicMergeDomains.main, args, 5)
    assert are_files_equal(outfile_domains.name, ROOT + 'hicMergeDomains/two_files.bed', delta=2)
    assert are_files_equal(outfile_ctcf_relation.name, ROOT + 'hicMergeDomains/two_files_relation.txt', delta=2)

    list_data = list(range(1, 23, 1))
    list_data = ["{}".format(x) for x in list_data]
    list_data.append('X')
    for i in list_data:
        res = compare_images(ROOT + '/hicMergeDomains/no_ctcf_plot/two_files_plot_' + i + '.png', plot_folder + '/two_files_plot_' + i + '.png', tol=40)
        assert res is None, res


# ---------------------------------------------------------------------------
# Characterization tests for the C++ port (cpp/tools/hicMergeDomains.cpp).
#
# The tests above compare with zip() and allow two differing lines, and three of
# them are xfail, so they cannot constrain the output. The tests below pin the
# Python's exact output on the real domain files. They need no graphviz: the
# rendering is replaced by a recorder that keeps every Digraph's filename and
# source, which is also what pins the DOT sources the port has to reproduce.
#
# ctcf_sorted.bed names chromosomes chr1..chrX while the domain files use 1..X,
# so the protein filter never acts on it. ctcf_sorted_nochr.bed is the same file
# with the prefix removed (sed 's/^chr//'), which makes --proteinFile and
# --minimumNumberOfPeaks observable.

DATA = ROOT + 'hicMergeDomains/'
TEN = DATA + '10kbtad_domains.bed'
FIFTY = DATA + '50kbtad_domains.bed'
HUNDRED = DATA + '100kbtad_domains.bed'
CTCF = DATA + 'ctcf_sorted.bed'
CTCF_NOCHR = DATA + 'ctcf_sorted_nochr.bed'


def _sha256(path):
    with open(path, 'rb') as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def _line_count(path):
    with open(path, 'rb') as handle:
        return handle.read().count(b'\n')


def _run(tmp_path, monkeypatch, args):
    """Run main in tmp_path with rendering recorded instead of executed."""
    graphs = []

    def recorder(self, filename=None, directory=None, view=False, cleanup=False, **kwargs):
        name = filename if filename is not None else self.filename
        graphs.append((name, self.source))
        return name

    monkeypatch.setattr(graphviz.Digraph, 'render', recorder)
    monkeypatch.chdir(tmp_path)
    hicMergeDomains.main(args + ['-om', 'merged.bed', '-or', 'relations.txt', '-ot', 'tree'])
    return graphs


def _read_recording(path):
    with open(path) as handle:
        content = handle.read()
    graphs = []
    for record in content.split('### ')[1:]:
        name, _, source = record.partition('\n')
        graphs.append((name, source))
    return graphs


def _file_bytes(path):
    with open(path, 'rb') as handle:
        return handle.read()


@pytest.mark.parametrize('args, merged_sha, merged_lines, relations_sha, relations_lines', [
    (['-d', TEN, '-p', CTCF_NOCHR],
     'ec33af6c3ddb5901b6faac9fc5b3f45fe4b16f3e3a6d6f08bf5e40d68e1fd0ba', 8007, None, None),
    (['-d', TEN, FIFTY],
     '75b51f0a9b7fade83bca39bd37e8415a17fb94f0440f64b8cfa6de74e3a4b6c4', 13329,
     'b08b1b5c0b56ce365bf278dc7fa38b107af4d9e4af12678d2df3754d53143339', 7152),
    (['-d', TEN, FIFTY, '-p', CTCF_NOCHR],
     'f5efa6570107a5caf54410c3404d730d3bfa7b5b47d23d9f8d50234df81d1429', 11966,
     '66042c1ff68e061c330760c014a1441f5c3c2c2c5e4c9d150e010dcac4fa678b', 6002),
    (['-d', TEN, FIFTY, '-p', CTCF_NOCHR, '-m', '3'],
     'aece7f8d347e615c45035384a4762588d8ce13240739d2cc00e607c2bea3263d', 11410,
     '82b4cc0c92d24680ef6c780673a7d77976b131450ac4d34051dc6ab605949c34', 5503),
    (['-d', TEN, FIFTY, '-v', '20000'],
     '932fcaa3fbc48a74e0962152d066ca20ce270fceb48fca29991ff4cbe3966444', 13038,
     '21d3e8b2f42066aaa2558eb516c3d3aba51f1a6db634458971733f1898d82cf0', 6847),
    (['-d', TEN, FIFTY, '-pe', '0.1'],
     '75b51f0a9b7fade83bca39bd37e8415a17fb94f0440f64b8cfa6de74e3a4b6c4', 13329,
     'c78277d431c609487a6a88021590254e2da24784ce1180c12bf9053f05770a48', 8315),
    (['-d', TEN, FIFTY, HUNDRED],
     '14c0b781be5db8ffc88bc311cd675d82201324450b456fdcbc406934a3f1ffcb', 15402,
     '92b002b8464cd310c4dabcd76bdbeb8581e2829be843ed1a3fa4295416d14204', 18108),
    (['-d', TEN, FIFTY, HUNDRED, '-p', CTCF_NOCHR],
     'aeca666a03587c3be6b96e272359d529914d4d7a05948573884ce13baa6f95d3', 13978,
     '7bb2c56e41846c42e09dfb5f47e58df27fb4f554ee9d0554326bc970f17bf001', 15635),
    (['-d', FIFTY, TEN],
     '63d3a8ebf81d2e2412031beb40d8e6d22bcf2e43171d54726d749afcc87ddb3f', 13330,
     'bee8136c6bf7d1a914da26d66ba8ee4d5422bc945d5e6821d63015479fe2147f', 7173),
], ids=['one_file_protein', 'two_files', 'two_files_protein', 'minimumNumberOfPeaks_3',
        'value_20000', 'percent_0.1', 'three_files', 'three_files_protein', 'reversed_order'])
def test_merged_and_relation_lists_are_pinned(tmp_path, monkeypatch, args, merged_sha,
                                              merged_lines, relations_sha, relations_lines):
    _run(tmp_path, monkeypatch, args)
    assert _line_count(tmp_path / 'merged.bed') == merged_lines
    assert _sha256(tmp_path / 'merged.bed') == merged_sha
    if relations_sha is None:
        assert not (tmp_path / 'relations.txt').exists()
    else:
        assert _line_count(tmp_path / 'relations.txt') == relations_lines
        assert _sha256(tmp_path / 'relations.txt') == relations_sha


def test_two_files_match_the_committed_references_exactly(tmp_path, monkeypatch):
    _run(tmp_path, monkeypatch, ['-d', TEN, FIFTY])
    assert _file_bytes(tmp_path / 'merged.bed') == _file_bytes(DATA + 'two_files.bed')
    assert _file_bytes(tmp_path / 'relations.txt') == _file_bytes(DATA + 'two_files_relation.txt')


def test_protein_file_with_other_chromosome_names_changes_nothing(tmp_path, monkeypatch):
    (tmp_path / 'one').mkdir()
    _run(tmp_path / 'one', monkeypatch, ['-d', TEN, '-p', CTCF])
    assert _file_bytes(tmp_path / 'one' / 'merged.bed') == _file_bytes(DATA + 'one_file')
    (tmp_path / 'two').mkdir()
    _run(tmp_path / 'two', monkeypatch, ['-d', TEN, FIFTY, '-p', CTCF])
    assert _file_bytes(tmp_path / 'two' / 'merged.bed') == _file_bytes(DATA + 'two_files_ctcf.bed')
    assert _file_bytes(DATA + 'two_files_ctcf.bed') == _file_bytes(DATA + 'two_files.bed')
    assert _file_bytes(tmp_path / 'two' / 'relations.txt') == \
        _file_bytes(DATA + 'two_files_relation_ctcf.txt')


@pytest.mark.parametrize('args, recording', [
    (['-d', TEN, FIFTY], 'two_files.txt'),
    (['-d', TEN, FIFTY, HUNDRED, '-p', CTCF_NOCHR], 'three_files_protein.txt'),
    (['-d', FIFTY, TEN], 'reversed_order.txt'),
], ids=['two_files', 'three_files_protein', 'reversed_order'])
def test_relation_tree_sources_are_pinned(tmp_path, monkeypatch, args, recording):
    # The same recordings pin the C++ DOT writer (cpp/tests/test_merge_domains.cpp).
    graphs = _run(tmp_path, monkeypatch, args)
    expected = _read_recording(DATA + 'tree_sources/' + recording)
    assert [name for name, _ in graphs] == [name for name, _ in expected]
    assert graphs == expected
    # hicMergeDomains.py:271 builds the first Digraph with strict=True, :297 the
    # others without it.
    assert graphs[0][1].startswith('strict digraph {\n')
    assert all(source.startswith('digraph {\n') for _, source in graphs[1:])


def test_reversed_order_writes_one_tad_twice_under_one_id(tmp_path, monkeypatch):
    # merge_list appends the last TAD of X twice as the same list object and
    # add_id numbers that object twice, so both lines carry the later ID.
    _run(tmp_path, monkeypatch, ['-d', FIFTY, TEN])
    with open(tmp_path / 'merged.bed') as handle:
        lines = handle.readlines()
    duplicated = sorted({line for line in lines if lines.count(line) > 1})
    assert duplicated == ['X\t154100000\t154500000\tID_13237\t-0.397066580167\t.\t154100000\t154500000\t51,160,44\n']


def test_one_domain_file_without_protein_file_exits_1(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit) as raised:
        hicMergeDomains.main(['-d', TEN, '-om', 'merged.bed'])
    assert raised.value.code == 1
    assert not (tmp_path / 'merged.bed').exists()


@pytest.mark.skipif(DOT_MISSING, reason='graphviz dot executable not on PATH')
def test_trees_are_rendered_and_their_sources_removed(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    hicMergeDomains.main(['-d', TEN, FIFTY, '-om', 'merged.bed', '-or', 'relations.txt',
                          '-ot', 'trees/tree', '-of', 'png'])
    rendered = sorted(os.listdir(tmp_path / 'trees'))
    expected = sorted('tree_{}.png'.format(chrom) for chrom in [str(i) for i in range(1, 23)] + ['X'])
    assert rendered == expected
    for name in rendered:
        assert _file_bytes(tmp_path / 'trees' / name)[:8] == b'\x89PNG\r\n\x1a\n'
