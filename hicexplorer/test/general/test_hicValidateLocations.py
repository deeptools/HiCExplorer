from hicexplorer import hicValidateLocations
import numpy.testing as nt

from tempfile import NamedTemporaryFile
import os


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

    args = "--data {} --protein {} --method {} --outFileName {} -r {} --addChrPrefixLoops".format(ROOT + 'loops_1.bedgraph',
                                                                                                  ROOT + 'GSM935376_hg19_Gm12878_Smc3.narrowPeak',
                                                                                                  'loops', outfile.name, 10000).split()
    hicValidateLocations.main(args)

    assert are_files_equal(
        ROOT + 'overlap_smc3_matched_locations', outfile.name + '_matched_locations')
    assert are_files_equal(ROOT + 'overlap_smc3_statistics',
                           outfile.name + '_statistics', skip=3)


def test_loop_broad_peak():
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()

    args = "--data {} --protein {} --method {} --outFileName {} -r {} --addChrPrefixLoops".format(ROOT + 'loops_1.bedgraph',
                                                                                                  ROOT + 'GSM733752_hg19_ctcf_GM12878.broadPeak',
                                                                                                  'loops', outfile.name, 10000).split()
    hicValidateLocations.main(args)

    # print(outfile.name + '_matched_locations')
    print(open(outfile.name + '_matched_locations', "r").read(1000))
    print(open(outfile.name + '_statistics', "r").read(1000))
    assert are_files_equal(
        ROOT + 'overlap_ctcf_matched_locations', outfile.name + '_matched_locations')
    assert are_files_equal(ROOT + 'overlap_ctcf_statistics',
                           outfile.name + '_statistics', skip=3)
