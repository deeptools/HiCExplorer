from hicexplorer import hicValidateLocations
import numpy.testing as nt

from tempfile import NamedTemporaryFile
import os


ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/hicValidateLocations/")


def are_files_equal(file1, file2, delta=1):
    equal = True
    if delta:
        mismatches = 0
    with open(file1) as textfile1, open(file2) as textfile2:
        for x, y in zip(textfile1, textfile2):
            # if x.startswith('File'):
            #     continue
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


def test_loop_narrow_peak():
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()

    args = "--data {} --protein {} --method {} --outFileName {} -r {}".format(ROOT + 'loops_1.bedgraph',
                                                                              ROOT + 'GSM935376_hg19_Gm12878_Smc3.narrowPeak',
                                                                              'loops', outfile.name, 10000).split()
    hicValidateLocations.main(args)

    # print(outfile.name + '_matched_locations')
    assert are_files_equal(ROOT + 'overlap_smc3_matched_locations', outfile.name + '_matched_locations')
    assert are_files_equal(ROOT + 'overlap_smc3_statistics', outfile.name + '_statistics')


def test_loop_broad_peak():
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()

    args = "--data {} --protein {} --method {} --outFileName {} -r {}".format(ROOT + 'loops_1.bedgraph',
                                                                              ROOT + 'GSM733752_hg19_ctcf_GM12878.broadPeak',
                                                                              'loops', outfile.name, 10000).split()
    hicValidateLocations.main(args)

    # print(outfile.name + '_matched_locations')
    assert are_files_equal(ROOT + 'overlap_ctcf_matched_locations', outfile.name + '_matched_locations')
    assert are_files_equal(ROOT + 'overlap_ctcf_statistics', outfile.name + '_statistics')
