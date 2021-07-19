from hicexplorer import chicSignificantInteractions
from tempfile import NamedTemporaryFile, mkdtemp
import os
import pytest
import h5py
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

ROOT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "test_data/cHi-C/")


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


def test_xFold():
    outfile_significant = NamedTemporaryFile(suffix='.hdf5', delete=False)
    outfile_target = NamedTemporaryFile(suffix='.hdf5', delete=False)

    # output_folder = mkdtemp(prefix="output_")
    # output_folder_target = mkdtemp(prefix="output_target_")

    outfile_significant.close()
    outfile_target.close()

    args = "--interactionFile {} --backgroundModelFile {} --range {} {} --outFileNameSignificant {} --outFileNameTarget {} --xFoldBackground {} --pValue {} -t {} --combinationMode dual".format(ROOT + 'chicViewpoint/two_matrices.hdf5',
                                                                                                                                                                                                 ROOT + 'background.txt',
                                                                                                                                                                                                 200000, 200000, outfile_significant.name,
                                                                                                                                                                                                 outfile_target.name, 1.5, 0.2, 1).split()
    chicSignificantInteractions.main(args)

    significantFileH5Object = h5py.File(outfile_significant.name, 'r')
    assert 'FL-E13-5_chr1' in significantFileH5Object
    assert 'MB-E10-5_chr1' in significantFileH5Object
    assert 'genes' in significantFileH5Object['FL-E13-5_chr1']
    assert 'genes' in significantFileH5Object['MB-E10-5_chr1']
    assert len(significantFileH5Object) == 2
    assert significantFileH5Object.attrs['type'] == 'significant'
    assert len(significantFileH5Object.attrs['range']) == 2
    assert significantFileH5Object.attrs['range'][0] == 200000
    assert significantFileH5Object.attrs['range'][1] == 200000

    # assert significantFileH5Object.attrs['averageContactBin'] == 5
    assert significantFileH5Object.attrs['fixateRange'] == 500000
    assert significantFileH5Object.attrs['mode_preselection'] == 'xfold'
    assert significantFileH5Object.attrs['mode_preselection_value'] == 1.5

    assert significantFileH5Object.attrs['pvalue'] == 0.2
    assert significantFileH5Object.attrs['combinationMode'] == 'dual'
    assert significantFileH5Object.attrs['truncateZeroPvalues'] == False
    assert significantFileH5Object.attrs['peakInteractionsThreshold'] == 5

    for chromosome in significantFileH5Object['FL-E13-5_chr1']:
        assert len(significantFileH5Object['FL-E13-5_chr1'][chromosome]) == 3
        for gene in significantFileH5Object['FL-E13-5_chr1'][chromosome]:
            assert len(significantFileH5Object['FL-E13-5_chr1'][chromosome][gene]) == 12
            for data in significantFileH5Object['FL-E13-5_chr1'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene', 'interaction_data_list', 'pvalue', 'raw', 'reference_point_end', 'reference_point_start', 'relative_position_list', 'start_list', 'sum_of_interactions', 'xfold']

    for chromosome in significantFileH5Object['MB-E10-5_chr1']:
        assert len(significantFileH5Object['MB-E10-5_chr1'][chromosome]) == 3
        for gene in significantFileH5Object['MB-E10-5_chr1'][chromosome]:
            assert len(significantFileH5Object['MB-E10-5_chr1'][chromosome][gene]) == 12
            for data in significantFileH5Object['MB-E10-5_chr1'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene', 'interaction_data_list', 'pvalue', 'raw', 'reference_point_end', 'reference_point_start', 'relative_position_list', 'start_list', 'sum_of_interactions', 'xfold']

    significantFileH5Object.close()

    targetFileH5Object = h5py.File(outfile_target.name, 'r')
    assert 'FL-E13-5_chr1' in targetFileH5Object
    assert 'MB-E10-5_chr1' in targetFileH5Object['FL-E13-5_chr1']

    # assert 'genes' in targetFileH5Object['FL-E13-5_chr1']
    assert len(targetFileH5Object) == 1
    assert len(targetFileH5Object['FL-E13-5_chr1']) == 1

    assert targetFileH5Object.attrs['type'] == 'target'
    assert len(targetFileH5Object.attrs['range']) == 2
    assert targetFileH5Object.attrs['range'][0] == 200000
    assert targetFileH5Object.attrs['range'][1] == 200000

    # assert targetFileH5Object.attrs['averageContactBin'] == 5
    assert targetFileH5Object.attrs['fixateRange'] == 500000
    assert targetFileH5Object.attrs['mode_preselection'] == 'xfold'
    assert targetFileH5Object.attrs['mode_preselection_value'] == 1.5

    assert targetFileH5Object.attrs['pvalue'] == 0.2
    assert targetFileH5Object.attrs['combinationMode'] == 'dual'
    assert targetFileH5Object.attrs['truncateZeroPvalues'] == False
    assert targetFileH5Object.attrs['peakInteractionsThreshold'] == 5

    for matrix in targetFileH5Object['FL-E13-5_chr1']:
        assert len(targetFileH5Object['FL-E13-5_chr1'][matrix]) == 2

        for chromosome in targetFileH5Object['FL-E13-5_chr1'][matrix]:

            assert len(targetFileH5Object['FL-E13-5_chr1'][matrix][chromosome]) == 3

        for gene in targetFileH5Object['FL-E13-5_chr1'][matrix][chromosome]:
            assert len(targetFileH5Object['FL-E13-5_chr1'][matrix][chromosome][gene]) == 5
            for data in targetFileH5Object['FL-E13-5_chr1'][matrix][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene', 'interaction_data_list', 'pvalue', 'raw', 'reference_point_end', 'reference_point_start', 'relative_position_list', 'start_list', 'sum_of_interactions', 'xfold']

    targetFileH5Object.close()


def test_loose_pvalue():
    outfile_significant = NamedTemporaryFile(suffix='.hdf5', delete=False)
    outfile_target = NamedTemporaryFile(suffix='.hdf5', delete=False)

    # output_folder = mkdtemp(prefix="output_")
    # output_folder_target = mkdtemp(prefix="output_target_")

    outfile_significant.close()
    outfile_target.close()

    args = "--interactionFile {} --backgroundModelFile {} --range {} {} --outFileNameSignificant {} --outFileNameTarget {} --loosePValue {} --pValue {} -t {} --combinationMode dual".format(ROOT + 'chicViewpoint/two_matrices.hdf5',
                                                                                                                                                                                             ROOT + 'background.txt',
                                                                                                                                                                                             200000, 200000, outfile_significant.name,
                                                                                                                                                                                             outfile_target.name, 0.5, 0.2, 1).split()
    chicSignificantInteractions.main(args)

    significantFileH5Object = h5py.File(outfile_significant.name, 'r')
    assert 'FL-E13-5_chr1' in significantFileH5Object
    assert 'MB-E10-5_chr1' in significantFileH5Object
    assert 'genes' in significantFileH5Object['FL-E13-5_chr1']
    assert 'genes' in significantFileH5Object['MB-E10-5_chr1']
    assert len(significantFileH5Object) == 2
    assert significantFileH5Object.attrs['type'] == 'significant'
    assert len(significantFileH5Object.attrs['range']) == 2
    assert significantFileH5Object.attrs['range'][0] == 200000
    assert significantFileH5Object.attrs['range'][1] == 200000

    # assert significantFileH5Object.attrs['averageContactBin'] == 5
    assert significantFileH5Object.attrs['fixateRange'] == 500000
    assert significantFileH5Object.attrs['mode_preselection'] == 'loosePValue'
    assert significantFileH5Object.attrs['mode_preselection_value'] == 0.5

    assert significantFileH5Object.attrs['pvalue'] == 0.2
    assert significantFileH5Object.attrs['combinationMode'] == 'dual'
    assert significantFileH5Object.attrs['truncateZeroPvalues'] == False
    assert significantFileH5Object.attrs['peakInteractionsThreshold'] == 5

    for chromosome in significantFileH5Object['FL-E13-5_chr1']:
        assert len(significantFileH5Object['FL-E13-5_chr1'][chromosome]) == 3
        for gene in significantFileH5Object['FL-E13-5_chr1'][chromosome]:
            assert len(significantFileH5Object['FL-E13-5_chr1'][chromosome][gene]) == 12
            for data in significantFileH5Object['FL-E13-5_chr1'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene', 'interaction_data_list', 'pvalue', 'raw', 'reference_point_end', 'reference_point_start', 'relative_position_list', 'start_list', 'sum_of_interactions', 'xfold']

    for chromosome in significantFileH5Object['MB-E10-5_chr1']:
        assert len(significantFileH5Object['MB-E10-5_chr1'][chromosome]) == 3
        for gene in significantFileH5Object['MB-E10-5_chr1'][chromosome]:
            assert len(significantFileH5Object['MB-E10-5_chr1'][chromosome][gene]) == 12
            for data in significantFileH5Object['MB-E10-5_chr1'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene', 'interaction_data_list', 'pvalue', 'raw', 'reference_point_end', 'reference_point_start', 'relative_position_list', 'start_list', 'sum_of_interactions', 'xfold']

    significantFileH5Object.close()

    targetFileH5Object = h5py.File(outfile_target.name, 'r')
    assert 'FL-E13-5_chr1' in targetFileH5Object
    assert 'MB-E10-5_chr1' in targetFileH5Object['FL-E13-5_chr1']

    # assert 'genes' in targetFileH5Object['FL-E13-5_chr1']
    assert len(targetFileH5Object) == 1
    assert len(targetFileH5Object['FL-E13-5_chr1']) == 1

    assert targetFileH5Object.attrs['type'] == 'target'
    assert len(targetFileH5Object.attrs['range']) == 2
    assert targetFileH5Object.attrs['range'][0] == 200000
    assert targetFileH5Object.attrs['range'][1] == 200000

    # assert targetFileH5Object.attrs['averageContactBin'] == 5
    assert targetFileH5Object.attrs['fixateRange'] == 500000
    assert targetFileH5Object.attrs['mode_preselection'] == 'loosePValue'
    assert targetFileH5Object.attrs['mode_preselection_value'] == 0.5

    assert targetFileH5Object.attrs['pvalue'] == 0.2
    assert targetFileH5Object.attrs['combinationMode'] == 'dual'
    assert targetFileH5Object.attrs['truncateZeroPvalues'] == False
    assert targetFileH5Object.attrs['peakInteractionsThreshold'] == 5

    for matrix in targetFileH5Object['FL-E13-5_chr1']:
        assert len(targetFileH5Object['FL-E13-5_chr1'][matrix]) == 2

        for chromosome in targetFileH5Object['FL-E13-5_chr1'][matrix]:

            assert len(targetFileH5Object['FL-E13-5_chr1'][matrix][chromosome]) == 3

        for gene in targetFileH5Object['FL-E13-5_chr1'][matrix][chromosome]:
            assert len(targetFileH5Object['FL-E13-5_chr1'][matrix][chromosome][gene]) == 5
            for data in targetFileH5Object['FL-E13-5_chr1'][matrix][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene', 'interaction_data_list', 'pvalue', 'raw', 'reference_point_end', 'reference_point_start', 'relative_position_list', 'start_list', 'sum_of_interactions', 'xfold']

    targetFileH5Object.close()


def test_loose_pvalue_single():
    outfile_significant = NamedTemporaryFile(suffix='.hdf5', delete=False)
    outfile_target = NamedTemporaryFile(suffix='.hdf5', delete=False)

    # output_folder = mkdtemp(prefix="output_")
    # output_folder_target = mkdtemp(prefix="output_target_")

    outfile_significant.close()
    outfile_target.close()

    args = "--interactionFile {} --backgroundModelFile {} --range {} {} --outFileNameSignificant {} --outFileNameTarget {} --loosePValue {} --pValue {} -t {} --combinationMode single".format(ROOT + 'chicViewpoint/two_matrices.hdf5',
                                                                                                                                                                                               ROOT + 'background.txt',
                                                                                                                                                                                               200000, 200000, outfile_significant.name,
                                                                                                                                                                                               outfile_target.name, 0.5, 0.2, 12).split()
    chicSignificantInteractions.main(args)

    significantFileH5Object = h5py.File(outfile_significant.name, 'r')
    assert 'FL-E13-5_chr1' in significantFileH5Object
    assert 'MB-E10-5_chr1' in significantFileH5Object
    assert 'genes' in significantFileH5Object['FL-E13-5_chr1']
    assert 'genes' in significantFileH5Object['MB-E10-5_chr1']
    assert len(significantFileH5Object) == 2
    assert significantFileH5Object.attrs['type'] == 'significant'
    assert len(significantFileH5Object.attrs['range']) == 2
    assert significantFileH5Object.attrs['range'][0] == 200000
    assert significantFileH5Object.attrs['range'][1] == 200000

    # assert significantFileH5Object.attrs['averageContactBin'] == 5
    assert significantFileH5Object.attrs['fixateRange'] == 500000
    assert significantFileH5Object.attrs['mode_preselection'] == 'loosePValue'
    assert significantFileH5Object.attrs['mode_preselection_value'] == 0.5

    assert significantFileH5Object.attrs['pvalue'] == 0.2
    assert significantFileH5Object.attrs['combinationMode'] == 'single'
    assert significantFileH5Object.attrs['truncateZeroPvalues'] == False
    assert significantFileH5Object.attrs['peakInteractionsThreshold'] == 5

    for chromosome in significantFileH5Object['FL-E13-5_chr1']:
        assert len(significantFileH5Object['FL-E13-5_chr1'][chromosome]) == 3
        for gene in significantFileH5Object['FL-E13-5_chr1'][chromosome]:
            assert len(significantFileH5Object['FL-E13-5_chr1'][chromosome][gene]) == 12
            for data in significantFileH5Object['FL-E13-5_chr1'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene', 'interaction_data_list', 'pvalue', 'raw', 'reference_point_end', 'reference_point_start', 'relative_position_list', 'start_list', 'sum_of_interactions', 'xfold']

    for chromosome in significantFileH5Object['MB-E10-5_chr1']:
        assert len(significantFileH5Object['MB-E10-5_chr1'][chromosome]) == 3
        for gene in significantFileH5Object['MB-E10-5_chr1'][chromosome]:
            assert len(significantFileH5Object['MB-E10-5_chr1'][chromosome][gene]) == 12
            for data in significantFileH5Object['MB-E10-5_chr1'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene', 'interaction_data_list', 'pvalue', 'raw', 'reference_point_end', 'reference_point_start', 'relative_position_list', 'start_list', 'sum_of_interactions', 'xfold']

    significantFileH5Object.close()

    targetFileH5Object = h5py.File(outfile_target.name, 'r')
    assert 'FL-E13-5_chr1' in targetFileH5Object
    assert 'MB-E10-5_chr1' in targetFileH5Object

    # assert 'MB-E10-5_chr1' in targetFileH5Object['FL-E13-5_chr1']

    # assert 'genes' in targetFileH5Object['FL-E13-5_chr1']
    assert len(targetFileH5Object) == 2
    assert len(targetFileH5Object['FL-E13-5_chr1']) == 2
    assert len(targetFileH5Object['MB-E10-5_chr1']) == 2

    assert targetFileH5Object.attrs['type'] == 'target'
    assert len(targetFileH5Object.attrs['range']) == 2
    assert targetFileH5Object.attrs['range'][0] == 200000
    assert targetFileH5Object.attrs['range'][1] == 200000

    # assert targetFileH5Object.attrs['averageContactBin'] == 5
    assert targetFileH5Object.attrs['fixateRange'] == 500000
    assert targetFileH5Object.attrs['mode_preselection'] == 'loosePValue'
    assert targetFileH5Object.attrs['mode_preselection_value'] == 0.5

    assert targetFileH5Object.attrs['pvalue'] == 0.2
    assert targetFileH5Object.attrs['combinationMode'] == 'single'
    assert targetFileH5Object.attrs['truncateZeroPvalues'] == False
    assert targetFileH5Object.attrs['peakInteractionsThreshold'] == 5

    # for matrix in targetFileH5Object['FL-E13-5_chr1']:
    #     assert len(targetFileH5Object['FL-E13-5_chr1'][matrix]) == 2

    for chromosome in targetFileH5Object['FL-E13-5_chr1']:

        assert len(targetFileH5Object['FL-E13-5_chr1'][chromosome]) == 3

        for gene in targetFileH5Object['FL-E13-5_chr1'][chromosome]:
            assert len(targetFileH5Object['FL-E13-5_chr1'][chromosome][gene]) == 5
            for data in targetFileH5Object['FL-E13-5_chr1'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene', 'interaction_data_list', 'pvalue', 'raw', 'reference_point_end', 'reference_point_start', 'relative_position_list', 'start_list', 'sum_of_interactions', 'xfold']

    for chromosome in targetFileH5Object['MB-E10-5_chr1']:

        assert len(targetFileH5Object['MB-E10-5_chr1'][chromosome]) == 3

        for gene in targetFileH5Object['MB-E10-5_chr1'][chromosome]:
            assert len(targetFileH5Object['MB-E10-5_chr1'][chromosome][gene]) == 5
            for data in targetFileH5Object['MB-E10-5_chr1'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene', 'interaction_data_list', 'pvalue', 'raw', 'reference_point_end', 'reference_point_start', 'relative_position_list', 'start_list', 'sum_of_interactions', 'xfold']

    targetFileH5Object.close()
