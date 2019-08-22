import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp

from hicexplorer import chicSignificantInteractions
ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/cHi-C/")


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


def test_xFold_one_interaction_file():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    output_folder = mkdtemp(prefix="output_")
    output_folder_target = mkdtemp(prefix="output_target_")

    outfile.close()
    args = "--interactionFile {} --backgroundModelFile {} --range {} {} -suffix {} --outputFolder {} --xFoldBackground {} --pValue {} --targetFolder {}".format(ROOT + 'chicViewpoint/output_1/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.bed',
                                                                                                                                                                ROOT + 'background.bed',
                                                                                                                                                                200000, 200000,
                                                                                                                                                                'output_significant.bed',
                                                                                                                                                                output_folder,
                                                                                                                                                                1.5, 0.2, output_folder_target).split()
    chicSignificantInteractions.main(args)

    assert are_files_equal(ROOT + "chicSignificantInteractions/output_1/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_output_significant.bed", output_folder + '/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_output_significant.bed')

    assert set(os.listdir(ROOT + "chicSignificantInteractions/output_1/")) == set(os.listdir(output_folder))
    assert set(os.listdir(ROOT + "chicSignificantInteractions/output_1_target/")) == set(os.listdir(output_folder_target))


def test_loose_pValue_one_interaction_file():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    output_folder = mkdtemp(prefix="output_")
    output_folder_target = mkdtemp(prefix="output_target_")

    outfile.close()
    args = "--interactionFile {} --backgroundModelFile {} --range {} {} -suffix {} --outputFolder {} --loosePValue {} --pValue {} --targetFolder {}".format(ROOT + 'chicViewpoint/output_1/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.bed',
                                                                                                                                                            ROOT + 'background.bed',
                                                                                                                                                            200000, 200000,
                                                                                                                                                            'output_significant.bed',
                                                                                                                                                            output_folder,
                                                                                                                                                            0.5, 0.2, output_folder_target).split()
    chicSignificantInteractions.main(args)

    assert are_files_equal(ROOT + "chicSignificantInteractions/output_2/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_output_significant.bed", output_folder + '/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_output_significant.bed')

    assert set(os.listdir(ROOT + "chicSignificantInteractions/output_2/")) == set(os.listdir(output_folder))
    assert set(os.listdir(ROOT + "chicSignificantInteractions/output_2_target/")) == set(os.listdir(output_folder_target))


def test_loose_pValue_multiple_interaction_files():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    output_folder = mkdtemp(prefix="output_")
    output_folder_target = mkdtemp(prefix="output_target_")

    outfile.close()
    args = "--interactionFile {} {} --backgroundModelFile {} --range {} {} -suffix {} --outputFolder {} --xFoldBackground {} --pValue {} --targetFolder {} -csn 1".format(ROOT + 'chicViewpoint/output_1/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.bed',
                                                                                                                                                                          ROOT + 'chicViewpoint/output_1/MB-E10-5_chr1_chr1_14300280_14300280_Eya1.bed',
                                                                                                                                                                          ROOT + 'background.bed',
                                                                                                                                                                          200000, 200000,
                                                                                                                                                                          'output_significant.bed',
                                                                                                                                                                          output_folder,
                                                                                                                                                                          1.5, 0.2, output_folder_target).split()
    chicSignificantInteractions.main(args)

    assert are_files_equal(ROOT + "chicSignificantInteractions/output_3/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_output_significant.bed", output_folder + '/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_output_significant.bed')
    assert are_files_equal(ROOT + "chicSignificantInteractions/output_3/MB-E10-5_chr1_chr1_14300280_14300280_Eya1_output_significant.bed", output_folder + '/MB-E10-5_chr1_chr1_14300280_14300280_Eya1_output_significant.bed')

    assert set(os.listdir(ROOT + "chicSignificantInteractions/output_3/")) == set(os.listdir(output_folder))
    assert set(os.listdir(ROOT + "chicSignificantInteractions/output_3_target/")) == set(os.listdir(output_folder_target))


def test_batchMode_xFold_csn_1():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile_target_list = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile_significant_files_list = NamedTemporaryFile(suffix='.txt', delete=False)

    output_folder = mkdtemp(prefix="output_")
    output_folder_target = mkdtemp(prefix="output_target_")

    outfile.close()
    args = "--interactionFile {} --interactionFileFolder {} --backgroundModelFile {} --range {} {} \
        -suffix {} --outputFolder {} --xFoldBackground {} --pValue {} --targetFolder {} --writeFileNamesToFile {} --targetFileList {} -bm -csn 1 ".format(ROOT + 'chicViewpoint/fileNames_one_matrix.txt',
                                                                                                                                                          ROOT + 'chicViewpoint/output_2',
                                                                                                                                                          ROOT + 'background.bed',
                                                                                                                                                          200000, 200000,
                                                                                                                                                          'output_significant.bed',
                                                                                                                                                          output_folder,
                                                                                                                                                          1.5, 0.2, output_folder_target,
                                                                                                                                                          outfile_significant_files_list.name,
                                                                                                                                                          outfile_target_list.name).split()
    chicSignificantInteractions.main(args)

    assert are_files_equal(ROOT + "chicSignificantInteractions/output_4/FL-E13-5_chr1_chr1_4487435_4487435_Sox17_output_significant.bed", output_folder + '/FL-E13-5_chr1_chr1_4487435_4487435_Sox17_output_significant.bed')
    assert are_files_equal(ROOT + "chicSignificantInteractions/output_4/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_output_significant.bed", output_folder + '/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_output_significant.bed')
    assert are_files_equal(ROOT + "chicSignificantInteractions/output_4/FL-E13-5_chr1_chr1_19093103_19093103_Tfap2d_output_significant.bed", output_folder + '/FL-E13-5_chr1_chr1_19093103_19093103_Tfap2d_output_significant.bed')

    assert are_files_equal(ROOT + "chicSignificantInteractions/output_4_target_list.txt", outfile_target_list.name)
    assert are_files_equal(ROOT + "chicSignificantInteractions/output_4_significant_files.txt", outfile_significant_files_list.name)

    assert set(os.listdir(ROOT + "chicSignificantInteractions/output_4/")) == set(os.listdir(output_folder))
    assert set(os.listdir(ROOT + "chicSignificantInteractions/output_4_target/")) == set(os.listdir(output_folder_target))


def test_batchMode_loose_pValue_csn_2():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile_target_list = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile_significant_files_list = NamedTemporaryFile(suffix='.txt', delete=False)

    output_folder = mkdtemp(prefix="output_")
    output_folder_target = mkdtemp(prefix="output_target_")

    outfile.close()
    args = "--interactionFile {} --interactionFileFolder {} --backgroundModelFile {} --range {} {} \
        -suffix {} --outputFolder {} --loosePValue {} --pValue {} --targetFolder {} --writeFileNamesToFile {} --targetFileList {} -bm -csn 2 ".format(ROOT + 'chicViewpoint/fileNames_two_matrices.txt',
                                                                                                                                                      ROOT + 'chicViewpoint/output_1',
                                                                                                                                                      ROOT + 'background.bed',
                                                                                                                                                      200000, 200000,
                                                                                                                                                      'output_significant.bed',
                                                                                                                                                      output_folder,
                                                                                                                                                      0.5, 0.2, output_folder_target,
                                                                                                                                                      outfile_significant_files_list.name,
                                                                                                                                                      outfile_target_list.name).split()
    chicSignificantInteractions.main(args)

    assert are_files_equal(ROOT + "chicSignificantInteractions/output_5/FL-E13-5_chr1_chr1_4487435_4487435_Sox17_output_significant.bed", output_folder + '/FL-E13-5_chr1_chr1_4487435_4487435_Sox17_output_significant.bed')
    assert are_files_equal(ROOT + "chicSignificantInteractions/output_5/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_output_significant.bed", output_folder + '/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_output_significant.bed')
    assert are_files_equal(ROOT + "chicSignificantInteractions/output_5/FL-E13-5_chr1_chr1_19093103_19093103_Tfap2d_output_significant.bed", output_folder + '/FL-E13-5_chr1_chr1_19093103_19093103_Tfap2d_output_significant.bed')
    assert are_files_equal(ROOT + "chicSignificantInteractions/output_5/MB-E10-5_chr1_chr1_4487435_4487435_Sox17_output_significant.bed", output_folder + '/MB-E10-5_chr1_chr1_4487435_4487435_Sox17_output_significant.bed')
    assert are_files_equal(ROOT + "chicSignificantInteractions/output_5/MB-E10-5_chr1_chr1_14300280_14300280_Eya1_output_significant.bed", output_folder + '/MB-E10-5_chr1_chr1_14300280_14300280_Eya1_output_significant.bed')
    assert are_files_equal(ROOT + "chicSignificantInteractions/output_5/MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_output_significant.bed", output_folder + '/MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_output_significant.bed')

    assert are_files_equal(ROOT + "chicSignificantInteractions/output_5_target/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_target.bed", output_folder_target + '/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_target.bed', skip=4)
    assert are_files_equal(ROOT + "chicSignificantInteractions/output_5_target/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_target.bed", output_folder_target + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_target.bed', skip=4)
    assert are_files_equal(ROOT + "chicSignificantInteractions/output_5_target/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_target.bed", output_folder_target + '/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_target.bed', skip=4)

    assert are_files_equal(ROOT + "chicSignificantInteractions/output_5_target_list.txt", outfile_target_list.name)
    assert are_files_equal(ROOT + "chicSignificantInteractions/output_5_significant_files.txt", outfile_significant_files_list.name)

    assert set(os.listdir(ROOT + "chicSignificantInteractions/output_5/")) == set(os.listdir(output_folder))
    assert set(os.listdir(ROOT + "chicSignificantInteractions/output_5_target/")) == set(os.listdir(output_folder_target))
