from hicexplorer import hicMergeLoops
import numpy.testing as nt

from tempfile import NamedTemporaryFile
import os


ROOT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "test_data/hicMergeLoops/")


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

    args = "-i {} {} {} -o {} -r {}".format(ROOT + 'gm12878_10kb.bedgraph', ROOT + 'gm12878_25kb.bedgraph',
                                            ROOT + 'gm12878_5kb.bedgraph', outfile.name, 5000).split()
    hicMergeLoops.main(args)

    assert are_files_equal(ROOT + 'gm12878_all.bedgraph', outfile.name, delta=2)
