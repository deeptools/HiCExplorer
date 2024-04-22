from hicexplorer import chicAggregateStatistic
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


def test_regular_mode():
    outfile_aggregate = NamedTemporaryFile(suffix='.hdf5', delete=False)
    outfile_aggregate.close()
    args = "--interactionFile {} --targetFile {} --outFileName {} \
            -t {}".format(ROOT + 'chicViewpoint/two_matrices.hdf5',
                          ROOT + 'chicSignificantInteractions/targetFile_dual.hdf5',
                          outfile_aggregate.name, 1).split()
    chicAggregateStatistic.main(args)

    aggregateFileH5Object = h5py.File(outfile_aggregate.name, 'r')
    assert 'FL-E13-5_chr1_MB-E10-5_chr1' in aggregateFileH5Object
    assert 'FL-E13-5_chr1' in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']
    assert 'MB-E10-5_chr1' in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']

    assert 'genes' in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['FL-E13-5_chr1']
    assert 'genes' in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['MB-E10-5_chr1']
    assert len(aggregateFileH5Object) == 1
    assert aggregateFileH5Object.attrs['type'] == 'aggregate'

    for chromosome in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['FL-E13-5_chr1']:

        assert len(aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['FL-E13-5_chr1'][chromosome]) == 3

        for gene in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['FL-E13-5_chr1'][chromosome]:
            assert len(aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['FL-E13-5_chr1'][chromosome][gene]) == 7
            for data in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['FL-E13-5_chr1'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene_name', 'raw_target_list', 'relative_distance_list', 'start_list', 'sum_of_interactions']

    for chromosome in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['MB-E10-5_chr1']:

        assert len(aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['MB-E10-5_chr1'][chromosome]) == 3

        for gene in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['MB-E10-5_chr1'][chromosome]:
            assert len(aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['MB-E10-5_chr1'][chromosome][gene]) == 7
            for data in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['MB-E10-5_chr1'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene_name', 'raw_target_list', 'relative_distance_list', 'start_list', 'sum_of_interactions']

    aggregateFileH5Object.close()


def test_regular_mode_threads():
    outfile_aggregate = NamedTemporaryFile(suffix='.hdf5', delete=False)
    outfile_aggregate.close()
    args = "--interactionFile {} --targetFile {} --outFileName {} \
            -t {}".format(ROOT + 'chicViewpoint/two_matrices.hdf5',
                          ROOT + 'chicSignificantInteractions/targetFile_dual.hdf5',
                          outfile_aggregate.name, 10).split()
    chicAggregateStatistic.main(args)

    aggregateFileH5Object = h5py.File(outfile_aggregate.name, 'r')
    assert 'FL-E13-5_chr1_MB-E10-5_chr1' in aggregateFileH5Object
    assert 'FL-E13-5_chr1' in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']
    assert 'MB-E10-5_chr1' in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']

    assert 'genes' in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['FL-E13-5_chr1']
    assert 'genes' in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['MB-E10-5_chr1']
    assert len(aggregateFileH5Object) == 1
    assert aggregateFileH5Object.attrs['type'] == 'aggregate'

    for chromosome in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['FL-E13-5_chr1']:

        assert len(aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['FL-E13-5_chr1'][chromosome]) == 3

        for gene in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['FL-E13-5_chr1'][chromosome]:
            assert len(aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['FL-E13-5_chr1'][chromosome][gene]) == 7
            for data in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['FL-E13-5_chr1'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene_name', 'raw_target_list', 'relative_distance_list', 'start_list', 'sum_of_interactions']

    for chromosome in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['MB-E10-5_chr1']:

        assert len(aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['MB-E10-5_chr1'][chromosome]) == 3

        for gene in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['MB-E10-5_chr1'][chromosome]:
            assert len(aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['MB-E10-5_chr1'][chromosome][gene]) == 7
            for data in aggregateFileH5Object['FL-E13-5_chr1_MB-E10-5_chr1']['MB-E10-5_chr1'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene_name', 'raw_target_list', 'relative_distance_list', 'start_list', 'sum_of_interactions']

    aggregateFileH5Object.close()


def test_target_list_bed3():
    outfile_aggregate = NamedTemporaryFile(suffix='.hdf5', delete=False)
    outfile_aggregate.close()
    args = "--interactionFile {} --targetFile {} --outFileName {} \
            -t {}".format(ROOT + 'chicViewpoint/two_matrices_custom_keys.hdf5',
                          ROOT + 'chicAggregateStatistic/target_list_3col.bed',
                          outfile_aggregate.name, 10).split()
    chicAggregateStatistic.main(args)

    aggregateFileH5Object = h5py.File(outfile_aggregate.name, 'r')
    assert 'c_adj_norm_t_adj_norm' in aggregateFileH5Object
    assert 'c_adj_norm' in aggregateFileH5Object['c_adj_norm_t_adj_norm']
    assert 't_adj_norm' in aggregateFileH5Object['c_adj_norm_t_adj_norm']

    assert 'genes' in aggregateFileH5Object['c_adj_norm_t_adj_norm']['c_adj_norm']
    assert 'genes' in aggregateFileH5Object['c_adj_norm_t_adj_norm']['t_adj_norm']
    assert len(aggregateFileH5Object) == 1
    assert aggregateFileH5Object.attrs['type'] == 'aggregate'

    for chromosome in aggregateFileH5Object['c_adj_norm_t_adj_norm']['c_adj_norm']:

        assert len(aggregateFileH5Object['c_adj_norm_t_adj_norm']['c_adj_norm'][chromosome]) == 3

        for gene in aggregateFileH5Object['c_adj_norm_t_adj_norm']['c_adj_norm'][chromosome]:
            assert len(aggregateFileH5Object['c_adj_norm_t_adj_norm']['c_adj_norm'][chromosome][gene]) == 7
            for data in aggregateFileH5Object['c_adj_norm_t_adj_norm']['c_adj_norm'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene_name', 'raw_target_list', 'relative_distance_list', 'start_list', 'sum_of_interactions']

    for chromosome in aggregateFileH5Object['c_adj_norm_t_adj_norm']['t_adj_norm']:

        assert len(aggregateFileH5Object['c_adj_norm_t_adj_norm']['t_adj_norm'][chromosome]) == 3

        for gene in aggregateFileH5Object['c_adj_norm_t_adj_norm']['t_adj_norm'][chromosome]:
            assert len(aggregateFileH5Object['c_adj_norm_t_adj_norm']['t_adj_norm'][chromosome][gene]) == 7
            for data in aggregateFileH5Object['c_adj_norm_t_adj_norm']['t_adj_norm'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene_name', 'raw_target_list', 'relative_distance_list', 'start_list', 'sum_of_interactions']

    aggregateFileH5Object.close()


def test_target_list_bed4():
    outfile_aggregate = NamedTemporaryFile(suffix='.hdf5', delete=False)
    outfile_aggregate.close()
    args = "--interactionFile {} --targetFile {} --outFileName {} \
            -t {}".format(ROOT + 'chicViewpoint/two_matrices_custom_keys.hdf5',
                          ROOT + 'chicAggregateStatistic/target_list_4col.bed',
                          outfile_aggregate.name, 10).split()
    chicAggregateStatistic.main(args)

    aggregateFileH5Object = h5py.File(outfile_aggregate.name, 'r')
    assert 'c_adj_norm_t_adj_norm' in aggregateFileH5Object
    assert 'c_adj_norm' in aggregateFileH5Object['c_adj_norm_t_adj_norm']
    assert 't_adj_norm' in aggregateFileH5Object['c_adj_norm_t_adj_norm']

    assert 'genes' in aggregateFileH5Object['c_adj_norm_t_adj_norm']['c_adj_norm']
    assert 'genes' in aggregateFileH5Object['c_adj_norm_t_adj_norm']['t_adj_norm']
    assert len(aggregateFileH5Object) == 1
    assert aggregateFileH5Object.attrs['type'] == 'aggregate'

    for chromosome in aggregateFileH5Object['c_adj_norm_t_adj_norm']['c_adj_norm']:

        assert len(aggregateFileH5Object['c_adj_norm_t_adj_norm']['c_adj_norm'][chromosome]) == 3

        for gene in aggregateFileH5Object['c_adj_norm_t_adj_norm']['c_adj_norm'][chromosome]:
            assert len(aggregateFileH5Object['c_adj_norm_t_adj_norm']['c_adj_norm'][chromosome][gene]) == 7
            for data in aggregateFileH5Object['c_adj_norm_t_adj_norm']['c_adj_norm'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene_name', 'raw_target_list', 'relative_distance_list', 'start_list', 'sum_of_interactions']

    for chromosome in aggregateFileH5Object['c_adj_norm_t_adj_norm']['t_adj_norm']:

        assert len(aggregateFileH5Object['c_adj_norm_t_adj_norm']['t_adj_norm'][chromosome]) == 3

        for gene in aggregateFileH5Object['c_adj_norm_t_adj_norm']['t_adj_norm'][chromosome]:
            assert len(aggregateFileH5Object['c_adj_norm_t_adj_norm']['t_adj_norm'][chromosome][gene]) == 7
            for data in aggregateFileH5Object['c_adj_norm_t_adj_norm']['t_adj_norm'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene_name', 'raw_target_list', 'relative_distance_list', 'start_list', 'sum_of_interactions']

    aggregateFileH5Object.close()