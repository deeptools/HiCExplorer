import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp

from hicexplorer import chicAggregateStatistic
ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/cHi-C/")


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
    output_folder = mkdtemp(prefix="output_")

    args = "--interactionFile {} {} --targetFile {} --outFileNameSuffix {} \
            --outputFolder {}".format(ROOT + 'chicViewpoint/output_1/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.bed',
                                      ROOT + 'chicViewpoint/output_1/MB-E10-5_chr1_chr1_14300280_14300280_Eya1.bed ',
                                      ROOT + 'chicSignificantInteractions/output_5_target/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_target.bed',
                                      'aggregated.bed',
                                      output_folder).split()
    chicAggregateStatistic.main(args)

    assert are_files_equal(ROOT + "chicAggregateStatistic/regular_mode/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_aggregated.bed", output_folder + '/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_aggregated.bed')
    assert are_files_equal(ROOT + "chicAggregateStatistic/regular_mode/MB-E10-5_chr1_chr1_14300280_14300280_Eya1_aggregated.bed", output_folder + '/MB-E10-5_chr1_chr1_14300280_14300280_Eya1_aggregated.bed')

    assert set(os.listdir(ROOT + "chicAggregateStatistic/regular_mode/")) == set(os.listdir(output_folder))


def test_batch_mode():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    output_folder = mkdtemp(prefix="output_")

    outfile.close()
    args = "--interactionFile {} --targetFile {} --outFileNameSuffix {} \
        --outputFolder {} -iff {} -tff {} -w {} -bm".format(ROOT + 'chicViewpoint/fileNames_two_matrices.txt',
                                                            ROOT + 'chicSignificantInteractions/output_5_target_list.txt',
                                                            'aggregated.bed',
                                                            output_folder,
                                                            ROOT + 'chicViewpoint/output_1',
                                                            ROOT + 'chicSignificantInteractions/output_5_target',
                                                            outfile.name).split()
    chicAggregateStatistic.main(args)
    assert are_files_equal(ROOT + "chicAggregateStatistic/batch_mode_file_names.txt", outfile.name)

    assert are_files_equal(ROOT + "chicAggregateStatistic/batch_mode/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_aggregated.bed", output_folder + '/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_aggregated.bed')
    assert are_files_equal(ROOT + "chicAggregateStatistic/batch_mode/MB-E10-5_chr1_chr1_14300280_14300280_Eya1_aggregated.bed", output_folder + '/MB-E10-5_chr1_chr1_14300280_14300280_Eya1_aggregated.bed')

    assert are_files_equal(ROOT + "chicAggregateStatistic/batch_mode/FL-E13-5_chr1_chr1_4487435_4487435_Sox17_aggregated.bed", output_folder + '/FL-E13-5_chr1_chr1_4487435_4487435_Sox17_aggregated.bed')
    assert are_files_equal(ROOT + "chicAggregateStatistic/batch_mode/MB-E10-5_chr1_chr1_4487435_4487435_Sox17_aggregated.bed", output_folder + '/MB-E10-5_chr1_chr1_4487435_4487435_Sox17_aggregated.bed')

    assert are_files_equal(ROOT + "chicAggregateStatistic/batch_mode/FL-E13-5_chr1_chr1_19093103_19093103_Tfap2d_aggregated.bed", output_folder + '/FL-E13-5_chr1_chr1_19093103_19093103_Tfap2d_aggregated.bed')
    assert are_files_equal(ROOT + "chicAggregateStatistic/batch_mode/MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_aggregated.bed", output_folder + '/MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d_aggregated.bed')

    assert set(os.listdir(ROOT + "chicAggregateStatistic/batch_mode/")) == set(os.listdir(output_folder))
