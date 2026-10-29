"""Equivalence cases for hicexplorer.chicDifferentialTest (apieq.py).

The tests run on the aggregate file of the cHi-C test data
(chicAggregateStatistic/aggregate.hdf5), sample against sample per gene, as
chicDifferentialTest's main() runs them, and on constructed tables for the
cases scipy rejects (zeros, negative and missing values, other shapes). The
result files are compared by content (text) or by an h5py dump of every
attribute and dataset.
"""
import argparse
import os

import numpy as np

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
CHIC = os.path.join(HERE, '..', 'test_data', 'cHi-C')
M = 'hicexplorer.chicDifferentialTest.'

AGGREGATE = os.path.join(CHIC, 'chicAggregateStatistic', 'aggregate.hdf5')
AGGREGATE_TARGET = os.path.join(CHIC, 'chicAggregateStatistic', 'aggregate_target.hdf')
COMBINATION = 'FL-E13-5_chr1_MB-E10-5_chr1'


def pairs(genes, combination=COMBINATION, first='FL-E13-5_chr1', second='MB-E10-5_chr1', chromosome='chr1'):
    return [[[combination, first, chromosome, gene], [combination, second, chromosome, gene]] for gene in genes]


GENES = pairs(['Eya1', 'Sox17', 'Tfap2d'])


def dump_h5(path):
    """Every group, attribute and dataset of an HDF5 file with h5py."""
    import h5py
    out = []

    def attributes(obj):
        return sorted((k, repr(v), type(v).__name__) for k, v in obj.attrs.items())

    with h5py.File(path, 'r') as f:
        out.append(('/', attributes(f)))

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                value = obj[()]
                if isinstance(value, np.ndarray):
                    value = value.tolist()
                out.append((name, 'dataset', obj.dtype.str, obj.shape, obj.compression, obj.compression_opts,
                            obj.chunks, value, attributes(obj)))
            else:
                out.append((name, 'group', attributes(obj)))
        f.visititems(visit)
    return out


def aggregated(path, sample):
    from hicexplorer.lib.viewpoint import Viewpoint
    return Viewpoint().readAggregatedFileHDF(path, sample)


def test_samples(target, path, first, second, alpha):
    """chisquare_test / fisher_exact_test on two samples of the aggregate file."""
    _, data1 = aggregated(path, first)
    _, data2 = aggregated(path, second)
    return target(data1, data2, alpha)


def run(target, path, files, test, alpha, use_queue=True):
    import queue
    from hicexplorer.lib.viewpoint import Viewpoint
    q = queue.Queue() if use_queue else None
    result = target(files, argparse.Namespace(aggregatedFile=path, statisticTest=test, alpha=alpha), Viewpoint(), q)
    if q is None:
        return result
    return result, q.get()


def _run_results(path, files, test, alpha):
    import queue
    from hicexplorer import chicDifferentialTest
    from hicexplorer.lib.viewpoint import Viewpoint
    q = queue.Queue()
    chicDifferentialTest.run_statistical_tests(files, argparse.Namespace(aggregatedFile=path, statisticTest=test,
                                                                         alpha=alpha), Viewpoint(), q)
    return q.get()


def write_text(target, out, files, test, alpha, which):
    accepted, rejected, all_results = _run_results(AGGREGATE, files, test, alpha)
    data = {'accepted': accepted, 'rejected': rejected, 'all': all_results}[which]
    for item in data:
        target(out, item, '# header old\nsecond line', '# header new\nx', alpha, test)
    if not data:
        target(out, [], '# header old', '# header new', alpha, test)
    with open(out, 'rb') as handle:
        return handle.read()


def write_text_files(target, out, first, second, test, alpha):
    """writeResult on two text interaction files, as 3.x wrote them."""
    from hicexplorer import chicDifferentialTest as d
    header1, lines1, data1 = d.readInteractionFile(first)
    header2, lines2, data2 = d.readInteractionFile(second)
    test_function = d.chisquare_test if test == 'chi2' else d.fisher_exact_test
    test_result, accepted, rejected = test_function(data1, data2, alpha)
    rows = [[lines1[i], lines2[i], result, data1[i], data2[i]] for i, result in enumerate(test_result)]
    target(out, rows, header1, header2, alpha, test)
    with open(out, 'rb') as handle:
        return handle.read()


def write_hdf(target, out, files, test, alpha):
    accepted, rejected, all_results = _run_results(AGGREGATE, files, test, alpha)
    target(out, accepted, rejected, all_results, files, alpha, test)
    return dump_h5(out)


def write_hdf_twice(target, out, files, test, alpha):
    accepted, rejected, all_results = _run_results(AGGREGATE, files, test, alpha)
    target(out, accepted + accepted, rejected + rejected, all_results + all_results, files + files, alpha, test)
    return dump_h5(out)


INTERACTION_TEXT = (b"# Interaction file, created with HiCExplorer's chicViewpoint version 3.2\n"
                    b"# matrix.cool\tchr1_100_100\t1.0Mb\t2.0Mb\tgeneA\tSum of interactions in fixate range: 1234.5\n"
                    b"# Chromosome\tStart\tEnd\n"
                    b"chr1\t100\t200\tgeneA\t1234.5\t-1000\t0.1\t0.5\t1.2\t3.0\n"
                    b"# comment\n"
                    b"\n"
                    b"chr1\t200\t300\tgeneA\t1234.5\t0\t0.2\t0.01\t2.5\t7.5\n"
                    b"single\n")

CASES = [
    # readInteractionFile --------------------------------------------------------
    Case('readInteractionFile/text', M + 'readInteractionFile', args=(File('i.txt'),),
         files={'i.txt': INTERACTION_TEXT}),
    Case('readInteractionFile/header_only', M + 'readInteractionFile', args=(File('i.txt'),),
         files={'i.txt': b"# first\n# a\tSum of interactions in fixate range: 7\n"}),
    Case('readInteractionFile/no_sum', M + 'readInteractionFile', args=(File('i.txt'),),
         files={'i.txt': b"# first\n# a\tb\n"}),
    Case('readInteractionFile/empty', M + 'readInteractionFile', args=(File('i.txt'),), files={'i.txt': b""}),
    Case('readInteractionFile/missing', M + 'readInteractionFile', args=(File('missing.txt'),)),

    # chisquare_test / fisher_exact_test ------------------------------------------
    Case('chisquare_test/eya1', M + 'chisquare_test', driver='test_samples',
         args=(AGGREGATE, GENES[0][0], GENES[0][1], 0.05)),
    Case('chisquare_test/tfap2d_alpha', M + 'chisquare_test', driver='test_samples',
         args=(AGGREGATE, GENES[2][0], GENES[2][1], 0.5)),
    Case('chisquare_test/aggregate_target', M + 'chisquare_test', driver='test_samples',
         args=(AGGREGATE_TARGET, ['c_adj_norm_t_adj_norm', 'c_adj_norm', 'chr19', 'A1BG'],
               ['c_adj_norm_t_adj_norm', 't_adj_norm', 'chr19', 'A1BG'], 0.05)),
    Case('chisquare_test/constructed', M + 'chisquare_test',
         args=([[10.0, 3.2], [0.0, 0.0], [5.0, -1.0], [float('nan'), 1.0], [1e6, 3.0], [4, 4]],
               [[5.0, 7.0], [3.0, 0.0], [5.0, 2.0], [2.0, 3.0], [1e6, 30.0], [4, 4]], 0.05)),
    Case('chisquare_test/ragged', M + 'chisquare_test', args=([[1.0, 2.0]], [[1.0]], 0.05)),
    Case('chisquare_test/empty', M + 'chisquare_test', args=([], [], 0.05)),
    Case('fisher_exact_test/eya1', M + 'fisher_exact_test', driver='test_samples',
         args=(AGGREGATE, GENES[0][0], GENES[0][1], 0.05)),
    Case('fisher_exact_test/sox17_alpha', M + 'fisher_exact_test', driver='test_samples',
         args=(AGGREGATE, GENES[1][0], GENES[1][1], 0.5)),
    Case('fisher_exact_test/aggregate_target', M + 'fisher_exact_test', driver='test_samples',
         args=(AGGREGATE_TARGET, ['c_adj_norm_t_adj_norm', 'c_adj_norm', 'chr10', 'A1CF'],
               ['c_adj_norm_t_adj_norm', 't_adj_norm', 'chr10', 'A1CF'], 0.05)),
    Case('fisher_exact_test/constructed', M + 'fisher_exact_test',
         args=([[10.0, 3.2], [0.0, 0.0], [5.0, -1.0], [float('nan'), 1.0], [3.0, 3.0], [100.0, 1.0], [1, 2, 3]],
               [[5.0, 7.0], [3.0, 0.0], [5.0, 2.0], [2.0, 3.0], [3.0, 3.0], [1.0, 100.0], [4, 5, 6]], 0.05)),
    Case('fisher_exact_test/mode_table', M + 'fisher_exact_test',
         args=([[5.0, 5.0], [12.0, 3.0]], [[5.0, 5.0], [11.0, 4.0]], 1.0)),

    # writeResult / writeResultHDF -----------------------------------------------
    Case('writeResult/fisher_all', M + 'writeResult', driver='write_text',
         args=(File('r.txt'), GENES[:1], 'fisher', 0.05, 'all')),
    Case('writeResult/chi2_rejected', M + 'writeResult', driver='write_text',
         args=(File('r.txt'), GENES[1:2], 'chi2', 0.5, 'rejected')),
    Case('writeResult/text_files_fisher', M + 'writeResult', driver='write_text_files',
         args=(File('r.txt'), File('a.txt'), File('b.txt'), 'fisher', 0.05),
         files={'a.txt': INTERACTION_TEXT, 'b.txt': INTERACTION_TEXT.replace(b'\t3.0\n', b'\t30.0\n')}),
    Case('writeResult/text_files_chi2', M + 'writeResult', driver='write_text_files',
         args=(File('r.txt'), File('a.txt'), File('b.txt'), 'chi2', 0.5),
         files={'a.txt': INTERACTION_TEXT, 'b.txt': INTERACTION_TEXT.replace(b'\t7.5\n', b'\t0.0\n')}),
    Case('writeResult/empty', M + 'writeResult', driver='write_text',
         args=(File('r.txt'), [], 'fisher', 0.05, 'all')),
    Case('writeResultHDF/fisher', M + 'writeResultHDF', driver='write_hdf',
         args=(File('d.hdf5'), GENES, 'fisher', 0.05)),
    Case('writeResultHDF/chi2', M + 'writeResultHDF', driver='write_hdf',
         args=(File('d.hdf5'), GENES[1:], 'chi2', 0.5)),
    Case('writeResultHDF/twice', M + 'writeResultHDF', driver='write_hdf_twice',
         args=(File('d.hdf5'), GENES[:1], 'fisher', 0.05)),
    Case('writeResultHDF/missing_directory', M + 'writeResultHDF', driver='write_hdf',
         args=(File('no/such/d.hdf5'), GENES[:1], 'fisher', 0.05)),

    # run_statistical_tests ------------------------------------------------------
    Case('run_statistical_tests/fisher', M + 'run_statistical_tests', driver='run',
         args=(AGGREGATE, GENES, 'fisher', 0.05)),
    Case('run_statistical_tests/chi2', M + 'run_statistical_tests', driver='run',
         args=(AGGREGATE, GENES, 'chi2', 0.05)),
    Case('run_statistical_tests/no_queue', M + 'run_statistical_tests', driver='run',
         args=(AGGREGATE, GENES[:1], 'chi2', 0.05, False)),
    Case('run_statistical_tests/unknown_test', M + 'run_statistical_tests', driver='run',
         args=(AGGREGATE, GENES[:1], 'ttest', 0.05)),
    Case('run_statistical_tests/missing_gene', M + 'run_statistical_tests', driver='run',
         args=(AGGREGATE, pairs(['Nope']), 'fisher', 0.05)),
    Case('run_statistical_tests/missing_gene_no_queue', M + 'run_statistical_tests', driver='run',
         args=(AGGREGATE, pairs(['Nope']), 'fisher', 0.05, False)),
    Case('run_statistical_tests/empty', M + 'run_statistical_tests', driver='run',
         args=(AGGREGATE, [], 'fisher', 0.05)),
]


# Cases whose floats differ from 3.7.6 in the last bits since the special
# functions (normal tail, incomplete beta and gamma functions, log-gamma) are
# Boost.Math's and no longer the routines scipy uses: the p-values agree
# to a relative 2e-14 or better. The cases compare every float at rtol 1e-12;
# integers, strings, coordinates and the structure of the result are still
# compared exactly, so a changed call, a moved TAD or a different row count
# fails.
_BOOST_LAST_BITS = {
    'chisquare_test/aggregate_target',
    'chisquare_test/constructed',
    'chisquare_test/eya1',
    'chisquare_test/tfap2d_alpha',
    'run_statistical_tests/chi2',
    'writeResultHDF/chi2',
}
for _case in CASES:
    if _case.id in _BOOST_LAST_BITS:
        _case.rtol = 1e-12
        _case.note = (_case.note + '; ' if _case.note else '') + (
            'Boost.Math special functions round the last bits of the p-values differently from scipy 1.14')
