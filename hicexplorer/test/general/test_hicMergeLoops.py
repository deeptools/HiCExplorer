from hicexplorer import hicMergeLoops
import numpy.testing as nt

from tempfile import NamedTemporaryFile
import os


ROOT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "test_data/hicMergeLoops/")


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


def test_loop_narrow_peak():
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()

    args = "-i {} {} {} -o {} -r {}".format(ROOT + 'gm12878_10kb.bedgraph', ROOT + 'gm12878_25kb.bedgraph',
                                            ROOT + 'gm12878_5kb.bedgraph', outfile.name, 5000).split()
    hicMergeLoops.main(args)

    assert are_files_equal(ROOT + 'gm12878_all.bedgraph', outfile.name)
