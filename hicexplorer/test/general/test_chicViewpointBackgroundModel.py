import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile
from sys import platform

import numpy as np
from hicexplorer import chicViewpointBackgroundModel

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/cHi-C/")


def are_files_equal(file1, file2, delta=1, skip=0, eps=0.1):

    mismatches = 0
    with open(file1, 'r') as textfile1:
        with open(file2, 'r') as textfile2:

            file1_content = textfile1.readlines()
            file2_content = textfile2.readlines()

            for i, (line1, line2) in enumerate(zip(file1_content, file2_content)):
                if i < skip:
                    continue
                line1_list = np.array(line1.split('\t'))
                line2_list = np.array(line2.split('\t'))

                line1_list = line1_list.astype(np.float)
                line2_list = line2_list.astype(np.float)

                for value1, value2 in zip(line1_list, line2_list):
                    if np.abs(value1 - value2) < eps:
                        continue
                    else:
                        mismatches += 1

    if mismatches < delta:
        return True
    else:
        return False


def test_compute_background_functional():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile.close()
    args = "--matrices {} {} --referencePoints {} -o {}".format(ROOT + 'FL-E13-5_chr1.cool', ROOT + 'MB-E10-5_chr1.cool',
                                                                ROOT + 'referencePoints.bed', outfile.name).split()
    chicViewpointBackgroundModel.main(args)

    assert are_files_equal(ROOT + 'background.bed', outfile.name, delta=20, skip=1)


def test_compute_background_number_of_lines():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile.close()
    args = "--matrices {} {} --referencePoints {} -o {}".format(ROOT + 'FL-E13-5_chr1.cool', ROOT + 'MB-E10-5_chr1.cool',
                                                                ROOT + 'referencePoints.bed', outfile.name).split()
    chicViewpointBackgroundModel.main(args)

    length_background = 0
    length_background_outfile = 0

    with open(ROOT + 'background.bed') as textfile:
        file_content = textfile.readlines()
        length_background = len(file_content)
    with open(outfile.name) as textfile:
        file_content = textfile.readlines()
        length_background_outfile = len(file_content)

    assert np.abs(length_background - length_background_outfile) < 1
