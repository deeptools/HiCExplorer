import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import os
from tempfile import NamedTemporaryFile, mkdtemp

from hicexplorer import chicViewpoint
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

def test_two_matrices():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    output_folder = mkdtemp(prefix="output_")
    outfile.close()
    args = "--matrices {} {} --referencePoints {} --backgroundModel {} --range {} {} -o {}".format(ROOT+ 'FL-E13-5_chr1.cool', 
                                                                                                    ROOT + 'MB-E10-5_chr1.cool', 
                                                                                                    ROOT + 'referencePoints_chicViewpoint.bed',
                                                                                                    ROOT + 'background.bed',
                                                                                                    200000, 200000,
                                                                                                    output_folder).split()
    chicViewpoint.main(args)

    assert are_files_equal(ROOT + "chicViewpoint/output_1/FL-E13-5_chr1_chr1_4487435_4487435_Sox17.bed", output_folder + '/FL-E13-5_chr1_chr1_4487435_4487435_Sox17.bed')
    assert are_files_equal(ROOT + "chicViewpoint/output_1/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.bed", output_folder + '/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.bed')
    assert are_files_equal(ROOT + "chicViewpoint/output_1/FL-E13-5_chr1_chr1_19093103_19093103_Tfap2d.bed", output_folder + '/FL-E13-5_chr1_chr1_19093103_19093103_Tfap2d.bed')
    assert are_files_equal(ROOT + "chicViewpoint/output_1/MB-E10-5_chr1_chr1_4487435_4487435_Sox17.bed", output_folder + '/MB-E10-5_chr1_chr1_4487435_4487435_Sox17.bed')
    assert are_files_equal(ROOT + "chicViewpoint/output_1/MB-E10-5_chr1_chr1_14300280_14300280_Eya1.bed", output_folder + '/MB-E10-5_chr1_chr1_14300280_14300280_Eya1.bed')
    assert are_files_equal(ROOT + "chicViewpoint/output_1/MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d.bed", output_folder + '/MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d.bed')

    assert set(os.listdir(ROOT + "chicViewpoint/output_1/")) == set(os.listdir(output_folder))



def test_one_matrix():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    output_folder = mkdtemp(prefix="output_")
    outfile.close()
    args = "--matrices {} --referencePoints {} --backgroundModel {} --range {} {} -o {}".format(ROOT+ 'FL-E13-5_chr1.cool',
                                                                                                    ROOT + 'referencePoints_chicViewpoint.bed',
                                                                                                    ROOT + 'background.bed',
                                                                                                    200000, 200000,
                                                                                                    output_folder).split()
    chicViewpoint.main(args)

    assert are_files_equal(ROOT + "chicViewpoint/output_2/FL-E13-5_chr1_chr1_4487435_4487435_Sox17.bed", output_folder + '/FL-E13-5_chr1_chr1_4487435_4487435_Sox17.bed')
    assert are_files_equal(ROOT + "chicViewpoint/output_2/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.bed", output_folder + '/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.bed')
    assert are_files_equal(ROOT + "chicViewpoint/output_2/FL-E13-5_chr1_chr1_19093103_19093103_Tfap2d.bed", output_folder + '/FL-E13-5_chr1_chr1_19093103_19093103_Tfap2d.bed')
   
    assert set(os.listdir(ROOT + "chicViewpoint/output_2/")) == set(os.listdir(output_folder))
    
def test_two_matrices_writeFileName():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile_name_list = NamedTemporaryFile(suffix='.txt', delete=False)

    output_folder = mkdtemp(prefix="output_")
    outfile.close()
    args = "--matrices {} {} --referencePoints {} --backgroundModel {} --range {} {} -o {} -w {}".format(ROOT+ 'FL-E13-5_chr1.cool', 
                                                                                                    ROOT + 'MB-E10-5_chr1.cool', 
                                                                                                    ROOT + 'referencePoints_chicViewpoint.bed',
                                                                                                    ROOT + 'background.bed',
                                                                                                    200000, 200000,
                                                                                                    output_folder, outfile_name_list.name).split()
    chicViewpoint.main(args)
    assert are_files_equal(ROOT + "chicViewpoint/fileNames_two_matrices.txt", outfile_name_list.name)

    assert are_files_equal(ROOT + "chicViewpoint/output_3/FL-E13-5_chr1_chr1_4487435_4487435_Sox17.bed", output_folder + '/FL-E13-5_chr1_chr1_4487435_4487435_Sox17.bed')
    assert are_files_equal(ROOT + "chicViewpoint/output_3/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.bed", output_folder + '/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.bed')
    assert are_files_equal(ROOT + "chicViewpoint/output_3/FL-E13-5_chr1_chr1_19093103_19093103_Tfap2d.bed", output_folder + '/FL-E13-5_chr1_chr1_19093103_19093103_Tfap2d.bed')
    assert are_files_equal(ROOT + "chicViewpoint/output_3/MB-E10-5_chr1_chr1_4487435_4487435_Sox17.bed", output_folder + '/MB-E10-5_chr1_chr1_4487435_4487435_Sox17.bed')
    assert are_files_equal(ROOT + "chicViewpoint/output_3/MB-E10-5_chr1_chr1_14300280_14300280_Eya1.bed", output_folder + '/MB-E10-5_chr1_chr1_14300280_14300280_Eya1.bed')
    assert are_files_equal(ROOT + "chicViewpoint/output_3/MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d.bed", output_folder + '/MB-E10-5_chr1_chr1_19093103_19093103_Tfap2d.bed')

    assert set(os.listdir(ROOT + "chicViewpoint/output_3/")) == set(os.listdir(output_folder))

def test_one_matrix_writeFileName():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile_name_list = NamedTemporaryFile(suffix='.txt', delete=False)

    output_folder = mkdtemp(prefix="output_")
    outfile.close()
    args = "--matrices {} --referencePoints {} --backgroundModel {} --range {} {} -o {} -w {}".format(ROOT+ 'FL-E13-5_chr1.cool', 
                                                                                                    ROOT + 'referencePoints_chicViewpoint.bed',
                                                                                                    ROOT + 'background.bed',
                                                                                                    200000, 200000,
                                                                                                    output_folder, outfile_name_list.name).split()
    chicViewpoint.main(args)
    assert are_files_equal(ROOT + "chicViewpoint/fileNames_one_matrix.txt", outfile_name_list.name)

    assert are_files_equal(ROOT + "chicViewpoint/output_4/FL-E13-5_chr1_chr1_4487435_4487435_Sox17.bed", output_folder + '/FL-E13-5_chr1_chr1_4487435_4487435_Sox17.bed')
    assert are_files_equal(ROOT + "chicViewpoint/output_4/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.bed", output_folder + '/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.bed')
    assert are_files_equal(ROOT + "chicViewpoint/output_4/FL-E13-5_chr1_chr1_19093103_19093103_Tfap2d.bed", output_folder + '/FL-E13-5_chr1_chr1_19093103_19093103_Tfap2d.bed')
   
    assert set(os.listdir(ROOT + "chicViewpoint/output_4/")) == set(os.listdir(output_folder))