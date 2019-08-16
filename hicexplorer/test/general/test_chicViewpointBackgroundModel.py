import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile

from hicexplorer import chicViewpointBackgroundModel

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


def test_compute_background():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile.close()
    args = "--matrices {} {} --referencePoints {} -o {}".format(ROOT + 'FL-E13-5_chr1.cool', ROOT + 'MB-E10-5_chr1.cool',
                                                                ROOT + 'referencePoints.bed', outfile.name).split()
    chicViewpointBackgroundModel.main(args)

    assert are_files_equal(ROOT + 'background.bed', outfile.name, skip=0)
