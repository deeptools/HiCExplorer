from hicexplorer import hicValidateLocations
import numpy.testing as nt

from tempfile import NamedTemporaryFile
import os
from hicexplorer.test.test_compute_function import compute


ROOT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "test_data/hicValidateLocations/")


def are_files_equal(file1, file2, delta=1, skip=0):

    lines_file1_dict = {}
    mismatches = 0
    matches = 0
    line_count_file1 = 0
    with open(file1, 'r') as textfile1:
        file_content = textfile1.readlines()

        for i, line in enumerate(file_content):
            if i < skip:
                continue
            lines_file1_dict[line] = True
            line_count_file1 += 1
    with open(file2, 'r') as textfile2:

        file_content = textfile2.readlines()
        for i, line in enumerate(file_content):
            if i < skip:
                continue
            if line in lines_file1_dict:
                matches += 1
            else:
                mismatches += 1
    if mismatches < delta and line_count_file1 - delta <= matches:
        return True
    else:
        return False


def test_loop_narrow_peak():
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()

    args = "--data {} --validationData {} --validationType {} --method {} --outFileName {} -r {} --chrPrefixLoops {} ".format(ROOT + 'loops_1.bedgraph',
                                                                                                                              ROOT + 'GSM935376_hg19_Gm12878_Smc3.narrowPeak', 'bed',
                                                                                                                              'loops', outfile.name, 10000, 'add').split()
    compute(hicValidateLocations.main, args, 5)
    assert are_files_equal(
        ROOT + 'overlap_smc3_matched_locations', outfile.name + '_matched_locations')
    assert are_files_equal(ROOT + 'overlap_smc3_statistics',
                           outfile.name + '_statistics', skip=3)


def test_loop_broad_peak():
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()

    args = "--data {} --validationData {} --validationType {} --method {} --outFileName {} -r {} --chrPrefixProtein {} ".format(ROOT + 'loops_1.bedgraph',
                                                                                                                                ROOT + 'GSM733752_hg19_ctcf_GM12878.broadPeak', 'bed',
                                                                                                                                'loops', outfile.name, 10000, 'remove').split()
    compute(hicValidateLocations.main, args, 5)

    print(open(outfile.name + '_matched_locations', "r").read(1000))
    print(open(outfile.name + '_statistics', "r").read(1000))
    assert are_files_equal(
        ROOT + 'overlap_ctcf_matched_locations', outfile.name + '_matched_locations')
    assert are_files_equal(ROOT + 'overlap_ctcf_statistics',
                           outfile.name + '_statistics', skip=3)


def test_loop_cool():
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()

    args = "--data {} --validationData {} --validationType {} --method {} --outFileName {} -r {} --chrPrefixLoop {} ".format(ROOT + 'loops_1.bedgraph',
                                                                                                                             ROOT + 'GSM1436265_RAD21_ENCFF002EMQ_10kb.cool', 'cool',
                                                                                                                             'loops', outfile.name, 10000, 'add').split()
    compute(hicValidateLocations.main, args, 5)

    print(open(outfile.name + '_matched_locations', "r").read(1000))
    print(open(outfile.name + '_statistics', "r").read(1000))
    assert are_files_equal(
        ROOT + 'overlap_rad21_cool_matched_locations', outfile.name + '_matched_locations')
    assert are_files_equal(ROOT + 'overlap_rad21_cool_statistics',
                           outfile.name + '_statistics', skip=3)
