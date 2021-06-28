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
