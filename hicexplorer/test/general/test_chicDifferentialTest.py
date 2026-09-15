from hicexplorer import chicDifferentialTest
from tempfile import NamedTemporaryFile, mkdtemp
import os
import pytest
import h5py
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

ROOT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "test_data/cHi-C/")


def are_files_equal(file1, file2, delta=2, skip=0):
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


def test_regular_mode_fisher():

    outfile_differential = NamedTemporaryFile(suffix='.hdf5', delete=False)
    outfile_differential.close()

    args = "--aggregatedFile {} --alpha {} --statisticTest {} --outFileName {} -t {}\
        ".format(ROOT + 'chicAggregateStatistic/aggregate.hdf5',
                 0.5, 'fisher',
                 outfile_differential.name, 1).split()
    chicDifferentialTest.main(args)

    differentialFileH5Object = h5py.File(outfile_differential.name, 'r')
    assert 'FL-E13-5_chr1' in differentialFileH5Object
    assert len(differentialFileH5Object) == 1

    assert 'MB-E10-5_chr1' in differentialFileH5Object['FL-E13-5_chr1']

    assert 'genes' in differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1']

    assert differentialFileH5Object.attrs['type'] == 'differential'
    assert differentialFileH5Object.attrs['alpha'] == 0.5
    assert differentialFileH5Object.attrs['test'] == 'fisher'

    for chromosome in differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1']:

        assert len(differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1'][chromosome]) == 3

        for gene in differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1'][chromosome]:
            assert len(differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1'][chromosome][gene]) == 3
            for data in differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1'][chromosome][gene]:
                assert len(differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1'][chromosome][gene][data]) == 10
                for status in differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1'][chromosome][gene][data]:

                    assert status in ['chromosome', 'end_list', 'gene', 'pvalue_list', 'raw_target_list_1', 'raw_target_list_2',
                                      'relative_distance_list', 'start_list', 'sum_of_interactions_1', 'sum_of_interactions_2']

    differentialFileH5Object.close()


def _category_counts(path):
    counts = {'accepted': 0, 'rejected': 0, 'all': 0}
    with h5py.File(path, 'r') as handle:
        def visit(name, obj):
            if name.endswith('pvalue_list') and '/genes/' not in name:
                counts[name.split('/')[-2]] += obj.shape[0]
        handle.visititems(visit)
    return counts


def test_characterization_calls_and_failures(tmp_path):
    # Characterization (cpp/AGENTS_CONTRACT.md rule 1), pinned before the C++
    # port: the calls of both tests on the committed aggregate file, and the
    # failures the tests above do not reach.
    aggregate = ROOT + 'chicAggregateStatistic/aggregate.hdf5'
    for test, counts, first in (('fisher', {'accepted': 7, 'rejected': 18, 'all': 25}, 0.00813445170440013),
                                ('chi2', {'accepted': 6, 'rejected': 19, 'all': 25}, 0.00866090872217125)):
        out = str(tmp_path / (test + '.hdf5'))
        chicDifferentialTest.main(["--aggregatedFile", aggregate, "--alpha", "0.5", "--statisticTest", test,
                                   "--outFileName", out, "-t", "1"])
        assert _category_counts(out) == counts
        with h5py.File(out, 'r') as handle:
            pvalues = handle['FL-E13-5_chr1']['MB-E10-5_chr1']['chr1']['Eya1']['all']['pvalue_list'][:]
            assert pvalues[0] == pytest.approx(first, rel=1e-12)
    out = str(tmp_path / 'failures.hdf5')
    with pytest.raises(ZeroDivisionError):
        chicDifferentialTest.main(["--aggregatedFile", aggregate, "--alpha", "0.5", "--outFileName", out, "-t", "0"])
    # A negative count starts no worker; writing the first reference point
    # then raises IndexError after its groups have been created.
    with pytest.raises(IndexError):
        chicDifferentialTest.main(["--aggregatedFile", aggregate, "--alpha", "0.5", "--outFileName", out, "-t", "-1"])
    with h5py.File(out, 'r') as handle:
        assert list(handle['FL-E13-5_chr1']['MB-E10-5_chr1']['chr1']['Eya1'].keys()) == ['accepted', 'all', 'rejected']
    with pytest.raises(SystemExit):
        chicDifferentialTest.main(["--aggregatedFile", ROOT + 'chicViewpoint/two_matrices.hdf5', "--alpha", "0.5",
                                   "--outFileName", str(tmp_path / 'wrong.hdf5'), "-t", "1"])


def test_regular_mode_chi2():

    outfile_differential = NamedTemporaryFile(suffix='.hdf5', delete=False)
    outfile_differential.close()

    args = "--aggregatedFile {} --alpha {} --statisticTest {} --outFileName {} -t {}\
        ".format(ROOT + 'chicAggregateStatistic/aggregate.hdf5',
                 0.5, 'chi2',
                 outfile_differential.name, 1).split()
    chicDifferentialTest.main(args)

    differentialFileH5Object = h5py.File(outfile_differential.name, 'r')
    assert 'FL-E13-5_chr1' in differentialFileH5Object
    assert len(differentialFileH5Object) == 1

    assert 'MB-E10-5_chr1' in differentialFileH5Object['FL-E13-5_chr1']

    assert 'genes' in differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1']

    assert differentialFileH5Object.attrs['type'] == 'differential'
    assert differentialFileH5Object.attrs['alpha'] == 0.5
    assert differentialFileH5Object.attrs['test'] == 'chi2'

    for chromosome in differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1']:

        assert len(differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1'][chromosome]) == 3

        for gene in differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1'][chromosome]:
            assert len(differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1'][chromosome][gene]) == 3
            for data in differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1'][chromosome][gene]:
                assert len(differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1'][chromosome][gene][data]) == 10
                for status in differentialFileH5Object['FL-E13-5_chr1']['MB-E10-5_chr1'][chromosome][gene][data]:

                    assert status in ['chromosome', 'end_list', 'gene', 'pvalue_list', 'raw_target_list_1', 'raw_target_list_2',
                                      'relative_distance_list', 'start_list', 'sum_of_interactions_1', 'sum_of_interactions_2']

    differentialFileH5Object.close()
