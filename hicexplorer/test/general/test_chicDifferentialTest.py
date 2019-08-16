import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp

from hicexplorer import chicDifferentialTest
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


def test_regular_mode_fisher():

    output_folder = mkdtemp(prefix="output_")

    args = "--interactionFile {} {} --alpha {} --statisticTest {} --outputFolder {}\
        ".format(ROOT + 'chicAggregateStatistic/batch_mode/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_aggregated.bed',
                 ROOT + 'chicAggregateStatistic/batch_mode/MB-E10-5_chr1_chr1_14300280_14300280_Eya1_aggregated.bed ',
                 0.5, 'fisher',
                 output_folder).split()
    chicDifferentialTest.main(args)

    assert are_files_equal(ROOT + "chicDifferentialTest/regular_mode_fisher/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_accepted.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_accepted.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/regular_mode_fisher/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_rejected.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_rejected.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/regular_mode_fisher/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_results.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_results.bed')
    assert set(os.listdir(ROOT + "chicDifferentialTest/regular_mode_fisher/")) == set(os.listdir(output_folder))


def test_regular_mode_chi2():

    output_folder = mkdtemp(prefix="output_")

    args = "--interactionFile {} {} --alpha {} --statisticTest {} --outputFolder {}\
        ".format(ROOT + 'chicAggregateStatistic/batch_mode/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_aggregated.bed',
                 ROOT + 'chicAggregateStatistic/batch_mode/MB-E10-5_chr1_chr1_14300280_14300280_Eya1_aggregated.bed ',
                 0.5, 'chi2',
                 output_folder).split()
    chicDifferentialTest.main(args)

    assert are_files_equal(ROOT + "chicDifferentialTest/regular_mode_chi2/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_accepted.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_accepted.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/regular_mode_chi2/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_rejected.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_rejected.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/regular_mode_chi2/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_results.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_results.bed')
    assert set(os.listdir(ROOT + "chicDifferentialTest/regular_mode_chi2/")) == set(os.listdir(output_folder))


def test_batch_mode_fisher():

    output_folder = mkdtemp(prefix="output_")

    args = "--interactionFile {} -iff {} --alpha {} --statisticTest {} --outputFolder {} -bm\
        ".format(ROOT + 'chicAggregateStatistic/batch_mode_file_names.txt',
                 ROOT + 'chicAggregateStatistic/batch_mode',
                 0.5, 'fisher',
                 output_folder).split()
    chicDifferentialTest.main(args)

    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_accepted.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_accepted.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_rejected.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_rejected.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_results.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_results.bed')

    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_H0_accepted.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_H0_accepted.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_H0_rejected.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_H0_rejected.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_results.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_results.bed')

    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_H0_accepted.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_H0_accepted.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_H0_rejected.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_H0_rejected.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_results.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_results.bed')

    assert set(os.listdir(ROOT + "chicDifferentialTest/batch_mode_fisher/")) == set(os.listdir(output_folder))


def test_batch_mode_chi2():

    output_folder = mkdtemp(prefix="output_")

    args = "--interactionFile {} -iff {} --alpha {} --statisticTest {} --outputFolder {} -bm\
        ".format(ROOT + 'chicAggregateStatistic/batch_mode_file_names.txt',
                 ROOT + 'chicAggregateStatistic/batch_mode',
                 0.5, 'chi2',
                 output_folder).split()
    chicDifferentialTest.main(args)

    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_chi2/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_accepted.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_accepted.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_chi2/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_rejected.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_rejected.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_chi2/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_results.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_results.bed')

    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_chi2/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_H0_accepted.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_H0_accepted.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_chi2/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_H0_rejected.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_H0_rejected.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_chi2/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_results.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_results.bed')

    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_chi2/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_H0_accepted.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_H0_accepted.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_chi2/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_H0_rejected.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_H0_rejected.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_chi2/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_results.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_results.bed')

    assert set(os.listdir(ROOT + "chicDifferentialTest/batch_mode_chi2/")) == set(os.listdir(output_folder))


def test_batch_mode_fisher():

    output_folder = mkdtemp(prefix="output_")

    args = "--interactionFile {} -iff {} --alpha {} --statisticTest {} --outputFolder {} -bm --rejectedFileNamesToFile {}\
        ".format(ROOT + 'chicAggregateStatistic/batch_mode_file_names.txt',
                 ROOT + 'chicAggregateStatistic/batch_mode',
                 0.5, 'fisher',
                 output_folder,
                 ).split()
    chicDifferentialTest.main(args)

    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher_outfile/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_accepted.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_accepted.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher_outfile/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_rejected.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_rejected.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher_outfile/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_results.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_results.bed')

    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher_outfile/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_H0_accepted.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_H0_accepted.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher_outfile/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_H0_rejected.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_H0_rejected.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher_outfile/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_results.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_results.bed')

    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher_outfile/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_H0_accepted.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_H0_accepted.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher_outfile/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_H0_rejected.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_H0_rejected.bed')
    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher_outfile/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_results.bed", output_folder + '/FL-E13-5_MB-E10-5_chr1_chr1_4487435_4487435_Sox17_results.bed')

    assert are_files_equal(ROOT + "chicDifferentialTest/batch_mode_fisher_outfile/rejectedFilesList.txt", output_folder + '/rejectedFilesList.txt')

    assert set(os.listdir(ROOT + "chicDifferentialTest/batch_mode_fisher_outfile/")) == set(os.listdir(output_folder))
