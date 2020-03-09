from hicexplorer import chicPlotViewpoint
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

tolerance = 50


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

# one viewpoint


def test_one_viewpoint():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)
    outfile.close()
    args = "-if {} --range {} {} -o {} -t {}".format(
        ROOT + 'chicViewpoint/output_1/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.txt', 200000, 200000, outfile.name, 1).split()
    chicPlotViewpoint.main(args)

    res = compare_images(
        ROOT + 'chicPlotViewpoint/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.png', outfile.name, tolerance)
    assert res is None, res

# n - viewpoints, jpg


def test_two_viewpoint_fileformat_dpi():
    outfile = NamedTemporaryFile(suffix='.pdf', delete=False)
    outfile.close()
    args = "-if {} {} --range {} {} -o {} --dpi {} --outputFormat {} -t {}".format(
        ROOT + 'chicViewpoint/output_1/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.txt',
        ROOT + 'chicViewpoint/output_1/MB-E10-5_chr1_chr1_14300280_14300280_Eya1.txt',
        200000, 200000, outfile.name, 100, 'pdf', 1).split()
    chicPlotViewpoint.main(args)

    res = compare_images(
        ROOT + 'chicPlotViewpoint/two_viewpoint.pdf', outfile.name, tolerance)
    assert res is None, res
# additional background model


def test_two_viewpoint_background():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)
    outfile.close()
    args = "-if {} {} --range {} {} -o {} -bmf {} -t {}".format(
        ROOT + 'chicViewpoint/output_1/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.txt',
        ROOT + 'chicViewpoint/output_1/MB-E10-5_chr1_chr1_14300280_14300280_Eya1.txt',
        200000, 200000, outfile.name, ROOT + 'background.txt', 1).split()
    chicPlotViewpoint.main(args)

    res = compare_images(
        ROOT + 'chicPlotViewpoint/two_viewpoint_background.png', outfile.name, tolerance)
    assert res is None, res


def test_two_viewpoint_background_significant():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)
    outfile.close()
    args = "-if {} {} --range {} {} -o {} -bmf {} --significantInteractions {} {} -t {}".format(
        ROOT + 'chicViewpoint/output_1/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.txt',
        ROOT + 'chicViewpoint/output_1/MB-E10-5_chr1_chr1_14300280_14300280_Eya1.txt',
        200000, 200000, outfile.name, ROOT + 'background.txt',
        ROOT + 'chicSignificantInteractions/output_3/FL-E13-5_chr1_chr1_14300280_14300280_Eya1_output_significant.txt',
        ROOT + 'chicSignificantInteractions/output_3/MB-E10-5_chr1_chr1_14300280_14300280_Eya1_output_significant.txt', 1).split()
    chicPlotViewpoint.main(args)

    res = compare_images(
        ROOT + 'chicPlotViewpoint/two_viewpoint_background_significant.png', outfile.name, tolerance)
    assert res is None, res


def test_two_viewpoint_background_differential():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)
    outfile.close()
    args = "-if {} {} --range {} {} -o {} -bmf {} --differentialTestResult {} -t {}".format(
        ROOT + 'chicViewpoint/output_1/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.txt',
        ROOT + 'chicViewpoint/output_1/MB-E10-5_chr1_chr1_14300280_14300280_Eya1.txt',
        200000, 200000, outfile.name, ROOT + 'background.txt',
        ROOT + 'chicDifferentialTest/batch_mode_fisher/FL-E13-5_MB-E10-5_chr1_chr1_14300280_14300280_Eya1_H0_rejected.txt', 1).split()
    chicPlotViewpoint.main(args)

    res = compare_images(
        ROOT + 'chicPlotViewpoint/two_viewpoint_background_differential.png', outfile.name, tolerance)
    assert res is None, res


def test_one_viewpoint_colormap_pvalue():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)
    outfile.close()
    args = "-if {} --range {} {} --pValue --colorMapPvalue {} -o {} --pValue --colorMapPvalue plasma --maxPValue 0.5 --minPValue 0.1 -t {}".format(
        ROOT + 'chicViewpoint/output_1/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.txt', 200000, 200000, 'plasma', outfile.name, 1).split()
    chicPlotViewpoint.main(args)

    res = compare_images(
        ROOT + 'chicPlotViewpoint/two_viewpoint_pvalue.png', outfile.name, tolerance)
    assert res is None, res


def test_one_viewpoint_per_file_batch_mode():
    output_folder = mkdtemp(prefix="output_")

    args = "-if {} -iff {} --range {} {} --outputFolder {} -psn {} -bm -t {}".format(ROOT + 'chicViewpoint/fileNames_one_matrix.txt',
                                                                                     ROOT + 'chicViewpoint/output_4',
                                                                                     200000, 200000,
                                                                                     output_folder, 1, 1).split()
    chicPlotViewpoint.main(args)

    res = compare_images(ROOT + 'chicPlotViewpoint/batchMode/one/Eya1_FL-E13-5_chr1_chr1_14300280.png',
                         output_folder + '/Eya1_FL-E13-5_chr1_chr1_14300280.png', tolerance)
    assert res is None, res
    res = compare_images(ROOT + 'chicPlotViewpoint/batchMode/one/Sox17_FL-E13-5_chr1_chr1_4487435.png',
                         output_folder + '/Sox17_FL-E13-5_chr1_chr1_4487435.png', tolerance)
    assert res is None, res
    res = compare_images(ROOT + 'chicPlotViewpoint/batchMode/one/Tfap2d_FL-E13-5_chr1_chr1_19093103.png',
                         output_folder + '/Tfap2d_FL-E13-5_chr1_chr1_19093103.png', tolerance)
    assert res is None, res
    assert set(os.listdir(ROOT + "chicPlotViewpoint/batchMode/one")
               ) == set(os.listdir(output_folder))


def test_two_viewpoints_per_file_batch_mode_significances_differential_pvalue():
    output_folder = mkdtemp(prefix="output_")

    args = "-if {} -iff {} --range {} {} --outputFolder {} -psn {} -bm  --differentialTestResult {} -diff {} -si {} -siff {} -p -xf {} -bmf {} -t {}".format(
        ROOT + 'chicViewpoint/fileNames_two_matrices.txt',
        ROOT + 'chicViewpoint/output_1',
        200000, 200000,
        output_folder, 2,
        ROOT + 'chicDifferentialTest/rejectedFilesList.txt',
        ROOT + 'chicDifferentialTest/batch_mode_fisher_outfile/',
        ROOT + 'chicSignificantInteractions/output_5_significant_files.txt',
        ROOT + 'chicSignificantInteractions/output_5/',
        1.5,
        ROOT + 'background.txt', 1).split()
    chicPlotViewpoint.main(args)

    res = compare_images(ROOT + 'chicPlotViewpoint/batchMode/two/Eya1_FL-E13-5_MB-E10-5_chr1_chr1_14300280.png',
                         output_folder + '/Eya1_FL-E13-5_MB-E10-5_chr1_chr1_14300280.png', tolerance)
    assert res is None, res
    res = compare_images(ROOT + 'chicPlotViewpoint/batchMode/two/Sox17_FL-E13-5_MB-E10-5_chr1_chr1_4487435.png',
                         output_folder + '/Sox17_FL-E13-5_MB-E10-5_chr1_chr1_4487435.png', tolerance)
    assert res is None, res
    res = compare_images(ROOT + 'chicPlotViewpoint/batchMode/two/Tfap2d_FL-E13-5_MB-E10-5_chr1_chr1_19093103.png',
                         output_folder + '/Tfap2d_FL-E13-5_MB-E10-5_chr1_chr1_19093103.png', tolerance)
    assert res is None, res
    assert set(os.listdir(ROOT + "chicPlotViewpoint/batchMode/two")
               ) == set(os.listdir(output_folder))
