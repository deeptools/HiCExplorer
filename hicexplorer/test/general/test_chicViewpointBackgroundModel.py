import logging
from hicexplorer import chicViewpointBackgroundModel
import numpy as np
from sys import platform
from tempfile import NamedTemporaryFile
import os
import pytest
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)


ROOT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "test_data/cHi-C/")
log = logging.getLogger(__name__)


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
                        log.debug('{}'.format(line1_list))
                        log.debug('{}'.format(line2_list))

                        mismatches += 1

    if mismatches < delta:
        return True
    else:
        log.debug('mismatches: {}'.format(mismatches))
        return False


def test_compute_background_functional():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile.close()
    args = "--matrices {} {} --referencePoints {} -o {} -t {}".format(ROOT + 'FL-E13-5_chr1.cool', ROOT + 'MB-E10-5_chr1.cool',
                                                                      ROOT + 'referencePoints.bed', outfile.name, 1).split()
    chicViewpointBackgroundModel.main(args)

    assert are_files_equal(ROOT + 'background.txt',
                           outfile.name, delta=700, skip=1)


def test_compute_background_functional_truncate_zeros():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile.close()
    args = "--matrices {} {} --referencePoints {} -o {} -t {} --truncateZeros".format(ROOT + 'FL-E13-5_chr1.cool', ROOT + 'MB-E10-5_chr1.cool',
                                                                                      ROOT + 'referencePoints.bed', outfile.name, 1).split()
    chicViewpointBackgroundModel.main(args)

    assert are_files_equal(ROOT + 'background_truncateZeros.txt',
                           outfile.name, delta=1000, skip=1)


def test_compute_background_number_of_lines():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile.close()
    args = "--matrices {} {} --referencePoints {} -o {} -t {}".format(ROOT + 'FL-E13-5_chr1.cool', ROOT + 'MB-E10-5_chr1.cool',
                                                                      ROOT + 'referencePoints.bed', outfile.name, 1).split()
    chicViewpointBackgroundModel.main(args)

    length_background = 0
    length_background_outfile = 0

    with open(ROOT + 'background.txt') as textfile:
        file_content = textfile.readlines()
        length_background = len(file_content)
    with open(outfile.name) as textfile:
        file_content = textfile.readlines()
        length_background_outfile = len(file_content)

    assert np.abs(length_background - length_background_outfile) < 1


@pytest.mark.xfail
def test_compute_background_functional_fail():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile.close()
    args = "--matrices {} {} --referencePoints {} -o {} -t {} --truncateZeros".format(ROOT + 'FL-E13-5_chr1.cool', ROOT + 'MB-E10-5_chr1.cool',
                                                                                      ROOT + 'referencePoints_qc.bed', outfile.name, 1).split()
    chicViewpointBackgroundModel.main(args)

    assert are_files_equal(ROOT + 'background.txt',
                           outfile.name, delta=700, skip=1)
