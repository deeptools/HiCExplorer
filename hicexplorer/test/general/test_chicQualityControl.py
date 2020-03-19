from hicexplorer import chicQualityControl
from matplotlib.testing.compare import compare_images
import matplotlib as mpl
from tempfile import NamedTemporaryFile, mkdtemp
import os
import pytest
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "test_data/cHi-C/")


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


def test_two_matrices():

    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile_histogram = NamedTemporaryFile(suffix='.png', delete=False)
    outfile_sparsity = NamedTemporaryFile(suffix='.png', delete=False)

    outfile.close()
    args = "--matrices {} {} --referencePoints {} --sparsity {} --outFileName {} --outFileNameHistogram {} --outFileNameSparsity  {} -t {}".format(ROOT + 'FL-E13-5_chr1.cool',
                                                                                                                                                   ROOT + 'MB-E10-5_chr1.cool',
                                                                                                                                                   ROOT + 'referencePoints.bed',
                                                                                                                                                   0.05,
                                                                                                                                                   outfile.name, outfile_histogram.name, outfile_sparsity.name, 1).split()
    chicQualityControl.main(args)

    assert are_files_equal(
        ROOT + "chicQualityControl/new_referencepoints.bed", outfile.name)
    assert are_files_equal(
        ROOT + "chicQualityControl/new_referencepoints.bed_raw_filter", outfile.name + '_raw_filter', skip=4)
    assert are_files_equal(ROOT + "chicQualityControl/new_referencepoints.bed_rejected_filter",
                           outfile.name + '_rejected_filter', skip=2)
    assert are_files_equal(ROOT + "chicQualityControl/new_referencepoints.bed_report",
                           outfile.name + '_rejected_filter', skip=2)

    res = compare_images(
        ROOT + "chicQualityControl/histogram.png", outfile_histogram.name, 50)
    assert res is None, res

    res = compare_images(
        ROOT + "chicQualityControl/sparsity.png", outfile_sparsity.name, 50)
    assert res is None, res
