from hicexplorer import chicViewpoint
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


def test_two_matrices():
    outfile = NamedTemporaryFile(suffix='.hdf5', delete=False)
    outfile.close()
    args = "--matrices {} {} --referencePoints {} --backgroundModel {} --range {} {} -o {} -t {}".format(ROOT + 'FL-E13-5_chr1.cool',
                                                                                                         ROOT + 'MB-E10-5_chr1.cool',
                                                                                                         ROOT + 'referencePoints_chicViewpoint.bed',
                                                                                                         ROOT + 'background.txt',
                                                                                                         200000, 200000,
                                                                                                         outfile.name, 1).split()
    chicViewpoint.main(args)

    interactionFileH5Object = h5py.File(outfile.name, 'r')
    assert 'FL-E13-5_chr1' in interactionFileH5Object
    assert 'MB-E10-5_chr1' in interactionFileH5Object
    assert 'genes' in interactionFileH5Object['FL-E13-5_chr1']
    assert 'genes' in interactionFileH5Object['MB-E10-5_chr1']
    assert len(interactionFileH5Object) == 2
    assert interactionFileH5Object.attrs['type'] == 'interactions'
    assert len(interactionFileH5Object.attrs['range']) == 2
    assert interactionFileH5Object.attrs['range'][0] == 200000
    assert interactionFileH5Object.attrs['range'][1] == 200000
    assert interactionFileH5Object.attrs['averageContactBin'] == 5
    assert interactionFileH5Object.attrs['fixateRange'] == 500000
    for chromosome in interactionFileH5Object['FL-E13-5_chr1']:
        assert len(interactionFileH5Object['FL-E13-5_chr1'][chromosome]) == 3
        for gene in interactionFileH5Object['FL-E13-5_chr1'][chromosome]:
            assert len(interactionFileH5Object['FL-E13-5_chr1'][chromosome][gene]) == 12
            for data in interactionFileH5Object['FL-E13-5_chr1'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene', 'interaction_data_list', 'pvalue', 'raw', 'reference_point_end', 'reference_point_start', 'relative_position_list', 'start_list', 'sum_of_interactions', 'xfold']

    for chromosome in interactionFileH5Object['MB-E10-5_chr1']:
        assert len(interactionFileH5Object['MB-E10-5_chr1'][chromosome]) == 3
        for gene in interactionFileH5Object['MB-E10-5_chr1'][chromosome]:
            assert len(interactionFileH5Object['MB-E10-5_chr1'][chromosome][gene]) == 12
            for data in interactionFileH5Object['MB-E10-5_chr1'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene', 'interaction_data_list', 'pvalue', 'raw', 'reference_point_end', 'reference_point_start', 'relative_position_list', 'start_list', 'sum_of_interactions', 'xfold']


def test_one_matrix():
    outfile = NamedTemporaryFile(suffix='.hdf5', delete=False)
    outfile.close()
    args = "--matrices {} --referencePoints {} --backgroundModel {} --range {} {} -o {} -t {}".format(ROOT + 'FL-E13-5_chr1.cool',
                                                                                                      ROOT + 'referencePoints_chicViewpoint.bed',
                                                                                                      ROOT + 'background.txt',
                                                                                                      200000, 200000,
                                                                                                      outfile.name, 1).split()
    chicViewpoint.main(args)

    interactionFileH5Object = h5py.File(outfile.name, 'r')
    assert 'FL-E13-5_chr1' in interactionFileH5Object
    assert len(interactionFileH5Object) == 1

    assert 'genes' in interactionFileH5Object['FL-E13-5_chr1']

    assert interactionFileH5Object.attrs['type'] == 'interactions'
    assert len(interactionFileH5Object.attrs['range']) == 2
    assert interactionFileH5Object.attrs['range'][0] == 200000
    assert interactionFileH5Object.attrs['range'][1] == 200000
    assert interactionFileH5Object.attrs['averageContactBin'] == 5
    assert interactionFileH5Object.attrs['fixateRange'] == 500000
    for chromosome in interactionFileH5Object['FL-E13-5_chr1']:
        assert len(interactionFileH5Object['FL-E13-5_chr1'][chromosome]) == 3
        for gene in interactionFileH5Object['FL-E13-5_chr1'][chromosome]:
            assert len(interactionFileH5Object['FL-E13-5_chr1'][chromosome][gene]) == 12
            for data in interactionFileH5Object['FL-E13-5_chr1'][chromosome][gene]:
                assert data in ['chromosome', 'end_list', 'gene', 'interaction_data_list', 'pvalue', 'raw', 'reference_point_end', 'reference_point_start', 'relative_position_list', 'start_list', 'sum_of_interactions', 'xfold']
