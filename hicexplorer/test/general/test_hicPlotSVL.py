import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicPlotSVL
from hicmatrix import HiCMatrix as hm

from tempfile import NamedTemporaryFile
import os
import numpy.testing as nt
import numpy as np
import pytest
from hicexplorer.test.test_compute_function import compute

# import pyBigWig
from matplotlib.testing.compare import compare_images

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")

import logging
log = logging.getLogger(__name__)


def are_files_equal(file1, file2, delta=1, skip=0):
    equal = True
    if delta:
        mismatches = 0
    with open(file1) as textfile1, open(file2) as textfile2:
        for i, (x, y) in enumerate(zip(textfile1, textfile2)):
            # if x.startswith('File'):
            #     continue
            if i < skip:
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


def test_plotSVL():
    plot = NamedTemporaryFile(suffix='.png', delete=False)
    outputFileName = NamedTemporaryFile(suffix='.txt', delete=False)
    outputFileNameData = NamedTemporaryFile(suffix='.txt', delete=False)

    plot.close()
    outputFileName.close()
    outputFileNameData.close()
    matrix = ROOT + "small_test_matrix_50kb_res.h5"
    matrix2 = ROOT + "small_test_matrix_50kb_res.h5"

    args = "--matrices {} {} --plotFileName {} --outFileName {} --outFileNameData {} --dpi 300"\
        .format(matrix, matrix2, plot.name, outputFileName.name, outputFileNameData.name).split()
    # hicPlotSVL.main(args)
    compute(hicPlotSVL.main, args, 5)

    log.debug('--plotFileName {} --outFileName {} --outFileNameData {} '.format(plot.name, outputFileName.name, outputFileNameData.name))
    log.debug('matrix {} {}'.format(matrix, matrix2))
    log.debug('matrix {} {}'.format(matrix.split('/')[-1], matrix2.split('/')[-1]))

    assert are_files_equal(ROOT + 'hicPlotSVL/data.txt', outputFileNameData.name, delta=2)
    assert are_files_equal(ROOT + 'hicPlotSVL/p_values.txt', outputFileName.name, delta=2)
    # res = compare_images(ROOT + 'hicPlotSVL/plot.png', plot.name, tol=50)
    # assert res is None, res
    os.unlink(plot.name)
    os.unlink(outputFileName.name)
    os.unlink(outputFileNameData.name)


# ---------------------------------------------------------------------------
# Characterization tests for the C++ port (cpp/tools/hicPlotSVL.cpp).
#
# test_plotSVL allows two differing lines out of each file and compares with
# zip(), which stops at the shorter file. The tests below compare every data row
# exactly and cover the options the Python suite never ran: --distance,
# --chromosomes, --threads and --colorList. The two header lines that carry the
# version and the matrix paths are left out, everything else is compared.

H5 = ROOT + 'small_test_matrix_50kb_res.h5'
COOL = ROOT + 'small_test_matrix_50kb_res.cool'


def _run(tmp_path, monkeypatch, args):
    # hicPlotSVL draws on pyplot's current figure and never closes it, so a
    # second main() in the same process would draw over the first plot.
    import matplotlib.pyplot as plt
    plt.close('all')
    monkeypatch.chdir(tmp_path)
    hicPlotSVL.main(args + ['--outFileName', 'p_values.txt', '--outFileNameData', 'data.txt'])


def _rows(path):
    with open(path) as handle:
        return [line.rstrip('\n') for line in handle if not line.startswith('#')]


def _p_values(path):
    with open(path) as handle:
        return [line.rstrip('\n').split('\t')[2] for line in handle if not line.startswith('#')]


def _chromosome_sums(chromosome, distance=2000000):
    matrix = hm.hiCMatrix(H5)
    first, last = matrix.getChrBinRange(chromosome)
    block = matrix.matrix[first:last, first:last]
    rows, columns = block.nonzero()
    short = np.absolute(rows - columns) <= distance / matrix.getBinSize()
    return np.sum(block.data[short]), np.sum(block.data[~short])


def test_reference_rows_are_misaligned_after_a_skipped_chromosome(tmp_path, monkeypatch):
    _run(tmp_path, monkeypatch, ['--matrices', H5, H5])
    rows = _rows(tmp_path / 'data.txt')
    assert rows == _rows(ROOT + 'hicPlotSVL/data.txt')
    assert len(rows) == 12
    assert sum(1 for row in rows if row.split('\t')[1] != '') == 7
    # chr3RHet has no contact beyond 2 Mb, so its ratio is inf and it is dropped;
    # the row named chr3RHet then shows the sums of the next kept chromosome.
    assert _chromosome_sums('chr3RHet')[1] == 0
    small, large = _chromosome_sums('chr3L')
    assert rows[1].split('\t')[:4] == ['chr3RHet', str(small / large), str(small), str(large)]
    assert _p_values(tmp_path / 'p_values.txt') == ['1.0']


def test_rows_are_named_after_the_last_matrix(tmp_path, monkeypatch):
    _run(tmp_path, monkeypatch, ['--matrices', H5, COOL])
    assert _rows(tmp_path / 'data.txt') == [
        'chr2RHet\t43.166666666666664\t259\t6\t43.166666666666664\t259\t6\t',
        'chr3RHet\t4.33640350877193\t9887\t2280\t4.33640350877193\t9887\t2280\t',
        'chr2LHet\t4.2679445350734095\t10465\t2452\t4.2679445350734095\t10465\t2452\t',
        'chr4\t6.5625\t105\t16\t6.5625\t105\t16\t',
        'chrYHet\t4.818731117824774\t3190\t662\t4.818731117824774\t3190\t662\t',
        'chr3L\t4.604556074766355\t7883\t1712\t4.604556074766355\t7883\t1712\t',
        'chr2L\t4.192916352675208\t11128\t2654\t4.192916352675208\t11128\t2654\t',
        'chrU\t\t6.0\t12\t2\t',
        'chrX\t\t\t', 'chrXHet\t\t\t', 'chr2R\t\t\t', 'chr3R\t\t\t', 'chrUextra\t\t\t',
        'chrM\t\t\t', 'chr3LHet\t\t\t']
    assert _p_values(tmp_path / 'p_values.txt') == ['0.8621866944007625']


def test_distance(tmp_path, monkeypatch):
    _run(tmp_path, monkeypatch, ['--matrices', H5, '--distance', '100000', '--threads', '20'])
    with open(tmp_path / 'data.txt') as handle:
        assert '# Short range contacts: <= 100000\n' in handle.readlines()
    assert _rows(tmp_path / 'data.txt') == [
        'chr2RHet\t0.9485294117647058\t129\t136\t',
        'chr3RHet\t1.325\t106\t80\t',
        'chr2LHet\t4.0\t8\t2\t',
        'chr4\t2.3580246913580245\t382\t162\t',
        'chr3L\t0.9454748960665175\t5913\t6254\t',
        'chr2L\t0.8177596397410639\t5811\t7106\t',
        'chrU\t1.1607142857142858\t65\t56\t',
        'chrX\t0.9183266932270916\t1844\t2008\t',
        'chrXHet\t2.5\t5\t2\t',
        'chr2R\t0.9808009909165979\t4751\t4844\t',
        'chr3R\t0.8695062398263701\t6410\t7372\t',
        'chr3LHet\t0.6544117647058824\t89\t136\t']
    assert not (tmp_path / 'p_values.txt').exists()


def test_chromosomes_keep_the_given_order(tmp_path, monkeypatch):
    _run(tmp_path, monkeypatch, ['--matrices', H5, '--chromosomes', 'chr2L', 'chr3R', 'chrX',
                                 '--threads', '1'])
    assert _rows(tmp_path / 'data.txt') == [
        'chr2L\t4.2679445350734095\t10465\t2452\t',
        'chr3R\t4.192916352675208\t11128\t2654\t',
        'chrX\t4.818731117824774\t3190\t662\t']


@pytest.mark.parametrize('threads', ['1', '3', '20'])
def test_threads_do_not_change_the_output(tmp_path, monkeypatch, threads):
    (tmp_path / 'default').mkdir()
    _run(tmp_path / 'default', monkeypatch, ['--matrices', COOL, H5])
    (tmp_path / 'threads').mkdir()
    _run(tmp_path / 'threads', monkeypatch, ['--matrices', COOL, H5, '--threads', threads])
    for name in ('data.txt', 'p_values.txt'):
        with open(tmp_path / 'default' / name) as a, open(tmp_path / 'threads' / name) as b:
            assert a.read() == b.read()


def test_negative_threads_evaluate_no_chromosome(tmp_path, monkeypatch):
    _run(tmp_path, monkeypatch, ['--matrices', COOL, '--threads', '-1'])
    rows = _rows(tmp_path / 'data.txt')
    assert len(rows) == 15
    assert all(row.endswith('\t\t') and row.count('\t') == 2 for row in rows)


@pytest.mark.parametrize('option', [['--colorList', 'r', 'k'], ['--dpi', '50']],
                         ids=['colorList', 'dpi'])
def test_plot_options_only_change_the_plot(tmp_path, monkeypatch, option):
    (tmp_path / 'default').mkdir()
    _run(tmp_path / 'default', monkeypatch, ['--matrices', H5, COOL])
    (tmp_path / 'changed').mkdir()
    _run(tmp_path / 'changed', monkeypatch, ['--matrices', H5, COOL] + option)
    for name in ('data.txt', 'p_values.txt'):
        with open(tmp_path / 'default' / name) as a, open(tmp_path / 'changed' / name) as b:
            assert a.read() == b.read()
    with open(tmp_path / 'default' / 'plot.png', 'rb') as a, \
            open(tmp_path / 'changed' / 'plot.png', 'rb') as b:
        assert a.read() != b.read()


def test_float64_sums(tmp_path, monkeypatch):
    _run(tmp_path, monkeypatch, ['--matrices', ROOT + 'Li_et_al_2015.h5', '--distance', '50000'])
    assert _rows(tmp_path / 'data.txt') == ['X\t2.090142873782224\t20618387.033307083\t9864582.60434466\t']


def test_float32_sums_print_as_float64_repr(tmp_path, monkeypatch):
    # np.sum over a float32 matrix is a float32, but '{}'.format of a numpy
    # float32 goes through float.__format__ and prints the widened value.
    _run(tmp_path, monkeypatch, ['--matrices', ROOT + 'hicNormalize/smallest_one.h5'])
    assert _rows(tmp_path / 'data.txt')[0] == 'chr2RHet\t48.16666793823242\t289.0\t6.0\t'
