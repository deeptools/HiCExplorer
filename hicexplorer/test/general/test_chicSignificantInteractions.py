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


def _common_args():
    return ["--interactionFile", ROOT + 'chicViewpoint/two_matrices.hdf5',
            "--backgroundModelFile", ROOT + 'background.txt', "--range", "200000", "200000",
            "--outFileNameSignificant", "significant.hdf5", "--outFileNameTarget", "target.hdf5"]


def _relative_positions(path):
    with h5py.File(path, 'r') as handle:
        return {sample: {gene: list(handle[sample]['chr1'][gene]['relative_position_list'][:])
                         for gene in handle[sample]['chr1']}
                for sample in handle}


def test_characterization_truncate_zero_pvalues(tmp_path, monkeypatch):
    # Characterization (cpp/AGENTS_CONTRACT.md rule 1): --truncateZeroPvalues
    # and a peak threshold on the merged raw sums, pinned before the C++ port.
    monkeypatch.chdir(tmp_path)
    chicSignificantInteractions.main(_common_args() + [
        "--loosePValue", "0.5", "--pValue", "0.2", "-t", "2", "--combinationMode", "single",
        "--truncateZeroPvalues", "--peakInteractionsThreshold", "10"])
    assert _relative_positions("significant.hdf5") == {
        'FL-E13-5_chr1': {'Eya1': [-1000], 'Sox17': [0], 'Tfap2d': [-1000]},
        'MB-E10-5_chr1': {'Eya1': [-1000], 'Sox17': [0], 'Tfap2d': [-15000, 0]}}
    with h5py.File("significant.hdf5", 'r') as handle:
        raw = handle['MB-E10-5_chr1']['chr1']['Tfap2d']['raw'][:]
        assert raw[0] == pytest.approx(15.0) and raw[1] == pytest.approx(276.4)
        assert handle.attrs['truncateZeroPvalues']
    with h5py.File("target.hdf5", 'r') as handle:
        assert [x.decode() for x in handle['MB-E10-5_chr1']['chr1']['Tfap2d']['start_list'][:]] == ['19077000']
    assert not os.path.exists("errorLog.txt")


def test_characterization_without_preselection(tmp_path, monkeypatch):
    # Without --xFoldBackground and --loosePValue every position is accepted
    # on its stored p-value, and --peakInteractionsThreshold is compared with
    # the x-fold (filter_by_pvalue), so positions with fewer than 5 raw
    # interactions are accepted. Sox17 has no accepted position, which
    # appends to errorLog.txt in the working directory.
    monkeypatch.chdir(tmp_path)
    chicSignificantInteractions.main(_common_args() + [
        "--pValue", "0.2", "-t", "1", "--combinationMode", "single"])
    positions = _relative_positions("significant.hdf5")
    assert positions['FL-E13-5_chr1']['Eya1'] == [-156000, -103000, -102000, -101000, 192000, 193000,
                                                   194000, 195000, 196000, 199000, 200000]
    assert positions['MB-E10-5_chr1'] == {'Eya1': [199000, 200000],
                                          'Tfap2d': [105000, 124000, 125000, 126000, 130000]}
    with h5py.File("significant.hdf5", 'r') as handle:
        assert max(handle['FL-E13-5_chr1']['chr1']['Eya1']['raw'][:]) < 5
    with open("errorLog.txt") as log:
        assert log.read() == ("Failed for: [['FL-E13-5_chr1', 'chr1', 'Sox17']].\n"
                              "Failed for: [['MB-E10-5_chr1', 'chr1', 'Sox17']].\n")


def test_characterization_reference_points_follow_computation_order(tmp_path, monkeypatch):
    # The significant file is written in np.unique order of the sample
    # triplets, but the reference points are indexed in computation order.
    monkeypatch.chdir(tmp_path)
    chicSignificantInteractions.main(_common_args() + [
        "--xFoldBackground", "1.5", "--pValue", "0.2", "-t", "1", "--combinationMode", "dual"])
    with h5py.File("significant.hdf5", 'r') as handle:
        sox17 = handle['FL-E13-5_chr1']['chr1']['Sox17']
        assert sox17['reference_point_start'][()] == 14300280
    with h5py.File("target.hdf5", 'r') as handle:
        sox17 = handle['FL-E13-5_chr1']['MB-E10-5_chr1']['chr1']['Sox17']
        assert sox17['reference_point_start'][()] == 4487435


def test_characterization_threshold_file_fails_when_written(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(TypeError):
        chicSignificantInteractions.main(_common_args() + [
            "--pValue", ROOT + 'thresholdFile_pValue.txt',
            "--loosePValue", ROOT + 'thresholdFile_loose_pValue.txt', "-t", "1"])
    with h5py.File("significant.hdf5", 'r') as handle:
        assert sorted(handle.attrs.keys()) == ['type', 'version']
        assert len(handle) == 0
    assert not os.path.exists("target.hdf5")


def test_characterization_thread_counts(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ZeroDivisionError):
        chicSignificantInteractions.main(_common_args() + [
            "--xFoldBackground", "1.5", "--pValue", "0.2", "-t", "0"])
    chicSignificantInteractions.main(_common_args() + [
        "--xFoldBackground", "1.5", "--pValue", "0.2", "-t", "-1"])
    for name in ("significant.hdf5", "target.hdf5"):
        with h5py.File(name, 'r') as handle:
            assert len(handle) == 0 and handle.attrs['pvalue'] == 0.2


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
