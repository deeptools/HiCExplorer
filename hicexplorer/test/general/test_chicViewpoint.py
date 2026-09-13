from hicexplorer import chicViewpoint
from tempfile import NamedTemporaryFile, mkdtemp
import os
import pytest
import h5py
import warnings

import numpy as np
from hicmatrix import HiCMatrix as hm
from scipy.special import betainc

from hicexplorer._version import __version__
from hicexplorer.lib import Viewpoint
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


# ---------------------------------------------------------------------------
# Characterization tests added for the C++ port (cpp/AGENTS_CONTRACT.md rule 1).
#
# What the two tests above actually constrain: the names of the groups and of
# the twelve datasets per viewpoint, and the root attributes. Not one value of
# any dataset is read, so a file with the right layout and arbitrary numbers
# passes both. The p-values, the x-fold, the relative interactions, the
# coordinates and the reference point near a chromosome end
# (adjustViewpointData) are untested.

MATRICES = [ROOT + 'FL-E13-5_chr1.cool', ROOT + 'MB-E10-5_chr1.cool']
BIN_SIZE = 1000
CHR1_LENGTH = 197195432


def run_viewpoint(pMatrices, pReferencePoints, pRange=(200000, 200000), pThreads=1):
    outfile = os.path.join(mkdtemp(prefix='chicViewpoint_'), 'viewpoints.hdf5')
    args = "--matrices {} --referencePoints {} --backgroundModelFile {} --range {} {} " \
        "-o {} -t {}".format(' '.join(pMatrices), pReferencePoints, ROOT + 'background.txt',
                             pRange[0], pRange[1], outfile, pThreads).split()
    chicViewpoint.main(args)
    return outfile


def write_reference_points(pLines):
    path = os.path.join(mkdtemp(prefix='chicViewpoint_rp_'), 'referencePoints.bed')
    with open(path, 'w') as handle:
        handle.write(''.join(line + '\n' for line in pLines))
    return path


def read_background_model(pPath):
    """{genomic relative position: (size, prob, mean value)}"""
    model = {}
    with open(pPath) as handle:
        handle.readline()
        for line in handle:
            fields = line.split('\t')
            model[int(fields[0])] = (float(fields[1]), float(fields[2]), float(fields[4]))
    return model


def expected_pvalues(pRaw, pReferenceIndex, pModel):
    """1 - cnb.cdf with the distribution Viewpoint.pvalues looks up.

    The lookup mixes units, and that is pinned here as the Python behaves:
    pvalues() computes the relative distance in bins (viewpoint.py:887), but the
    background model is keyed in base pairs (chicViewpointBackgroundModel.py:255
    writes relative_position * bin_size). With 1 kb bins the only bin distance
    that is also a key is 0, so the reference point uses its own distribution
    and every other position falls back to the model at the smallest or the
    largest key, that is at -fixateRange or +fixateRange.
    """
    result = np.empty(len(pRaw))
    for i, raw in enumerate(pRaw):
        if i == pReferenceIndex:
            key = 0
        elif i < pReferenceIndex:
            key = min(pModel)
        else:
            key = max(pModel)
        size, prob, _ = pModel[key]
        result[i] = 1.0 if raw == 0.0 else 1 - betainc(size, raw + 1, prob)
    result[~np.isfinite(result)] = 1.0
    return result


def dense_rows(pMatrix, pFirstBin, pLastBin, pLow, pHigh):
    """Rows pFirstBin..pLastBin over columns pLow..pHigh, the viewpoint bins summed."""
    rows = np.asarray(pMatrix.matrix[pFirstBin:pLastBin + 1, pLow:pHigh + 1].todense())
    summed = rows.sum(axis=0)
    return np.concatenate([summed[:pFirstBin - pLow],
                           [summed[pFirstBin - pLow:pLastBin - pLow + 1].sum()],
                           summed[pLastBin - pLow + 1:]])


def all_datasets(pFile):
    names = []
    pFile.visit(names.append)
    return sorted(name for name in names if isinstance(pFile[name], h5py.Dataset))


@pytest.fixture(scope='module')
def two_matrices_output():
    return run_viewpoint(MATRICES, ROOT + 'referencePoints_chicViewpoint.bed')


def test_datasets_match_the_reference_file(two_matrices_output):
    """Every dataset but the p-values, exactly; the p-values are tested below.

    The committed two_matrices.hdf5 predates the background model now in
    background.txt, and its p-values differ, so they are recomputed instead.
    """
    with h5py.File(ROOT + 'chicViewpoint/two_matrices.hdf5', 'r') as expected, \
            h5py.File(two_matrices_output, 'r') as observed:
        assert observed.attrs['type'] == 'interactions'
        assert observed.attrs['version'] == __version__
        assert list(observed.attrs['range']) == [200000, 200000]
        assert observed.attrs['range'].dtype == np.int64
        assert observed.attrs['averageContactBin'] == 5
        assert observed.attrs['fixateRange'] == 500000
        assert observed.attrs['resolution'] == 1000
        assert sorted(observed.attrs.keys()) == \
            ['averageContactBin', 'fixateRange', 'range', 'resolution', 'type', 'version']

        assert all_datasets(observed) == all_datasets(expected)
        compared = 0
        for name in all_datasets(expected):
            left = expected[name]
            right = observed[name]
            assert left.dtype == right.dtype, name
            assert left.shape == right.shape, name
            assert (right.compression, right.compression_opts) == \
                (left.compression, left.compression_opts), name
            if left.shape != ():
                assert (right.compression, right.compression_opts) == ('gzip', 9), name
            if name.endswith('/pvalue'):
                continue
            if left.dtype.kind == 'f':
                np.testing.assert_array_equal(right[()], left[()], err_msg=name)
            else:
                assert np.array_equal(left[()], right[()]), name
            compared += 1
        assert compared == 2 * 3 * 11

        for matrix in ('FL-E13-5_chr1', 'MB-E10-5_chr1'):
            assert sorted(observed[matrix].keys()) == ['chr1', 'genes']
            assert sorted(observed[matrix + '/genes'].keys()) == ['Eya1', 'Sox17', 'Tfap2d']
            for gene in ('Eya1', 'Sox17', 'Tfap2d'):
                # a hard link to the same group, not a copy
                assert observed[matrix + '/genes/' + gene] == observed[matrix + '/chr1/' + gene]


def test_pvalues_follow_the_background_lookup_rule(two_matrices_output):
    model = read_background_model(ROOT + 'background.txt')
    between = 0
    with h5py.File(two_matrices_output, 'r') as observed:
        for matrix in ('FL-E13-5_chr1', 'MB-E10-5_chr1'):
            for gene in ('Eya1', 'Sox17', 'Tfap2d'):
                group = observed[matrix + '/chr1/' + gene]
                raw = group['raw'][()]
                reference_index = int(np.flatnonzero(group['relative_position_list'][()] == 0)[0])
                assert reference_index == 200
                pvalues = group['pvalue'][()]
                np.testing.assert_allclose(pvalues, expected_pvalues(raw, reference_index, model),
                                           rtol=1e-12, atol=0)
                between += np.count_nonzero((pvalues > 0.0) & (pvalues < 1.0))
    assert between > 0


def test_derived_datasets_are_consistent(two_matrices_output):
    """Everything in a viewpoint group recomputed from the matrix itself."""
    model = read_background_model(ROOT + 'background.txt')
    reference_points = {'Sox17': 4487435, 'Eya1': 14300280, 'Tfap2d': 19093103}
    with h5py.File(two_matrices_output, 'r') as observed:
        for matrix_file, matrix in zip(MATRICES, ('FL-E13-5_chr1', 'MB-E10-5_chr1')):
            hic = hm.hiCMatrix(matrix_file)
            for gene, position in reference_points.items():
                group = observed[matrix + '/chr1/' + gene]
                bin_id = position // BIN_SIZE
                assert group['chromosome'][()].decode() == 'chr1'
                assert group['gene'][()].decode() == gene
                assert group['reference_point_start'][()] == position
                assert group['reference_point_end'][()] == position
                assert group['reference_point_start'].dtype == np.int64

                # the denominator runs over fixateRange, the data over --range
                fixate = dense_rows(hic, bin_id, bin_id, bin_id - 500, bin_id + 500)
                assert group['sum_of_interactions'][()] == np.sum(fixate)
                window = dense_rows(hic, bin_id, bin_id, bin_id - 200, bin_id + 200)
                raw = Viewpoint().smoothInteractionValues(window, 5)
                np.testing.assert_array_equal(group['raw'][()], raw)
                interaction = raw / np.sum(fixate)
                np.testing.assert_array_equal(group['interaction_data_list'][()], interaction)

                bins = np.arange(bin_id - 200, bin_id + 201)
                np.testing.assert_array_equal(group['start_list'][()], bins * BIN_SIZE)
                np.testing.assert_array_equal(group['end_list'][()], (bins + 1) * BIN_SIZE)
                relative = (bins - bin_id) * BIN_SIZE
                np.testing.assert_array_equal(group['relative_position_list'][()], relative)
                means = np.array([model[key][2] for key in relative])
                np.testing.assert_array_equal(group['xfold'][()], interaction / means)


def test_wide_and_chromosome_end_reference_points():
    """A reference point over four bins, and the adjustViewpointData path.

    NearEnd sits 50 kb before the end of chr1, at an offset of 200 bp inside
    its bin. calculateViewpointRange then gives a range whose background keys
    reach +51 kb, one bin further than the matrix does, so the lengths differ
    and adjustViewpointData runs. It walks range(start, end) with the end
    exclusive (chicViewpoint.py:91-96), which drops the last bin, and it
    returns the data as float32 (chicViewpoint.py:104).
    """
    reference_points = write_reference_points([
        'chr1\t24202042\t24205042\tWide',
        'chr1\t197145200\t197145200\tNearEnd',
        'chr1\t150500\t150500\tNearStart',
        'chr1\t197195000\tThreeColumns'])
    outfile = run_viewpoint(MATRICES[:1], reference_points)
    hic = hm.hiCMatrix(MATRICES[0])
    model = read_background_model(ROOT + 'background.txt')
    with h5py.File(outfile, 'r') as observed:
        assert sorted(observed['FL-E13-5_chr1/chr1'].keys()) == \
            ['NearEnd', 'NearStart', 'ThreeColumns', 'Wide']

        wide = observed['FL-E13-5_chr1/chr1/Wide']
        assert wide['reference_point_start'][()] == 24202042
        assert wide['reference_point_end'][()] == 24205042
        fixate = dense_rows(hic, 24202, 24205, 23702, 24705)
        assert wide['sum_of_interactions'][()] == np.sum(fixate)
        window = dense_rows(hic, 24202, 24205, 24002, 24405)
        assert len(window) == 401
        np.testing.assert_array_equal(wide['raw'][()], Viewpoint().smoothInteractionValues(window, 5))
        starts = wide['start_list'][()]
        relative = wide['relative_position_list'][()]
        # the summed element carries the coordinates of the first bin only
        assert (starts[199], starts[200], starts[201]) == (24201000, 24202000, 24206000)
        assert wide['end_list'][()][200] == 24203000
        assert (relative[199], relative[200], relative[201]) == (-1000, 0, 1000)
        np.testing.assert_array_equal(wide['pvalue'][()],
                                      expected_pvalues(wide['raw'][()], 200, model))

        near_end = observed['FL-E13-5_chr1/chr1/NearEnd']
        raw = near_end['raw'][()]
        assert len(raw) == 250
        fixate = dense_rows(hic, 197145, 197145, 196645, 197195)
        total = np.sum(fixate)
        assert total > 0
        assert near_end['sum_of_interactions'][()] == total
        window = dense_rows(hic, 197145, 197145, 196945, 197195)
        assert len(window) == 251
        expected_raw = Viewpoint().smoothInteractionValues(window[:250].astype(np.float32), 5)
        np.testing.assert_array_equal(raw, expected_raw)
        assert np.count_nonzero(raw) > 0
        assert np.all(raw == raw.astype(np.float32))
        np.testing.assert_array_equal(near_end['start_list'][()], np.arange(196945, 197195) * BIN_SIZE)
        relative = near_end['relative_position_list'][()]
        assert (relative[0], relative[200], relative[-1]) == (-200000, 0, 49000)
        interaction = raw / total
        np.testing.assert_array_equal(near_end['interaction_data_list'][()], interaction)
        means = np.array([model[key][2] for key in relative])
        np.testing.assert_array_equal(near_end['xfold'][()], interaction / means)
        np.testing.assert_array_equal(near_end['pvalue'][()], expected_pvalues(raw, 200, model))

        # No contacts at all near the start of chr1: the denominator is 0, which
        # computeRelativeValues treats as absent and divides by the sum of the
        # data instead, 0 as well, so every relative value is NaN.
        near_start = observed['FL-E13-5_chr1/chr1/NearStart']
        assert len(near_start['raw'][()]) == 351
        assert near_start['relative_position_list'][()][0] == -150000
        assert near_start['start_list'][()][0] == 0
        assert near_start['sum_of_interactions'][()] == 0.0
        assert np.all(np.isnan(near_start['interaction_data_list'][()]))
        assert np.all(near_start['pvalue'][()] == 1.0)

        three = observed['FL-E13-5_chr1/chr1/ThreeColumns']
        assert three['reference_point_start'][()] == three['reference_point_end'][()] == 197195000
        relative = three['relative_position_list'][()]
        assert len(relative) == 200
        assert relative[-1] == -1000
        assert not np.any(relative == 0)


def test_thread_count_does_not_change_the_output():
    reference_points = write_reference_points([
        'chr1\t4487435\t4487435\tSox17',
        'chr1\t24202042\t24205042\tWide',
        'chr1\t197145200\t197145200\tNearEnd',
        'chr1\t14300280\t14300280\tEya1',
        'chr1\t19093103\t19093103\tTfap2d'])
    single = run_viewpoint(MATRICES, reference_points, pThreads=1)
    several = run_viewpoint(MATRICES, reference_points, pThreads=3)
    with h5py.File(single, 'r') as left, h5py.File(several, 'r') as right:
        assert all_datasets(left) == all_datasets(right)
        for name in all_datasets(left):
            assert left[name].dtype == right[name].dtype
            if left[name].dtype.kind == 'f':
                np.testing.assert_array_equal(left[name][()], right[name][()], err_msg=name)
            else:
                assert np.array_equal(left[name][()], right[name][()]), name


def test_a_bad_reference_point_exits_non_zero():
    """A start after the end on chr1, and a chromosome the matrix does not have."""
    outfile = os.path.join(mkdtemp(prefix='chicViewpoint_'), 'viewpoints.hdf5')
    args = "--matrices {} --referencePoints {} --backgroundModelFile {} --range 200000 200000 " \
        "-o {} -t 1".format(MATRICES[0], ROOT + 'referencePoints_chicViewpoint_crash.bed',
                            ROOT + 'background.txt', outfile).split()
    with pytest.raises(SystemExit) as exit_info:
        chicViewpoint.main(args)
    assert exit_info.value.code == 1
