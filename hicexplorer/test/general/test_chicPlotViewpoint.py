from hicexplorer import chicPlotViewpoint
from matplotlib.testing.compare import compare_images
from matplotlib.testing.exceptions import ImageComparisonFailure

import matplotlib as mpl
import matplotlib.pyplot as plt


from tempfile import NamedTemporaryFile, mkdtemp
import os
import pytest
import tarfile
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


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
def test_one_gene():
    outfile = NamedTemporaryFile(suffix='.tar.gz', delete=False)
    outfile.close()
    args = "-if {} --range {} {} -o {} -t {} --combinationMode oneGene --combinationName Eya1".format(
        ROOT + 'chicViewpoint/two_matrices.hdf5', 200000, 200000, outfile.name, 1).split()
    chicPlotViewpoint.main(args)

    file_obj = tarfile.open(outfile.name, "r")

    namelist = file_obj.getnames()
    assert len(namelist) == 3
    assert 'FL-E13-5_chr1_Eya1.png' in namelist

    output_folder = mkdtemp(prefix="output_")

    file_obj.extractall(output_folder)
    res = compare_images(
        ROOT + 'chicPlotViewpoint/FL-E13-5_chr1_chr1_14300280_14300280_Eya1.png', output_folder + '/FL-E13-5_chr1_Eya1.png', tolerance)
    assert res is None, res


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
def test_one_gene_background_differential_significant():
    outfile = NamedTemporaryFile(suffix='.tar.gz', delete=False)
    outfile.close()
    args = "-if {} --range {} {} -o {} -t {} --combinationMode dual --backgroundModelFile {} --differentialTestResult {} --significantInteractions {} --plotSignificantInteractions".format(
        ROOT + 'chicViewpoint/two_matrices.hdf5', 200000, 200000, outfile.name, 1, ROOT + 'background.txt', ROOT + 'chicDifferentialTest/differential.hdf5', ROOT + 'chicSignificantInteractions/significantInteractions_dual.hdf5').split()
    chicPlotViewpoint.main(args)

    file_obj = tarfile.open(outfile.name, "r")

    namelist = file_obj.getnames()
    assert len(namelist) == 3
    assert 'FL-E13-5_chr1_MB-E10-5_chr1_Eya1.png' in namelist
    assert 'FL-E13-5_chr1_MB-E10-5_chr1_Sox17.png' in namelist
    assert 'FL-E13-5_chr1_MB-E10-5_chr1_Tfap2d.png' in namelist

    output_folder = mkdtemp(prefix="output_")

    file_obj.extractall(output_folder)
    res = compare_images(
        ROOT + 'chicPlotViewpoint/background_diff_significant/FL-E13-5_chr1_MB-E10-5_chr1_Eya1.png', output_folder + '/FL-E13-5_chr1_MB-E10-5_chr1_Eya1.png', tolerance)
    assert res is None, res
    res = compare_images(
        ROOT + 'chicPlotViewpoint/background_diff_significant/FL-E13-5_chr1_MB-E10-5_chr1_Sox17.png', output_folder + '/FL-E13-5_chr1_MB-E10-5_chr1_Sox17.png', tolerance)
    assert res is None, res
    res = compare_images(
        ROOT + 'chicPlotViewpoint/background_diff_significant/FL-E13-5_chr1_MB-E10-5_chr1_Tfap2d.png', output_folder + '/FL-E13-5_chr1_MB-E10-5_chr1_Tfap2d.png', tolerance)
    assert res is None, res


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
def test_single_significant_pvalue():
    outfile = NamedTemporaryFile(suffix='.tar.gz', delete=False)
    outfile.close()
    args = "-if {} --range {} {} -o {} -t {} --combinationMode single --backgroundModelFile {} --significantInteractions {} --pValue --plotSignificantInteractions".format(
        ROOT + 'chicViewpoint/two_matrices.hdf5', 200000, 200000, outfile.name, 1, ROOT + 'background.txt', ROOT + 'chicSignificantInteractions/significantInteractions_dual.hdf5').split()
    chicPlotViewpoint.main(args)

    file_obj = tarfile.open(outfile.name, "r")

    namelist = file_obj.getnames()
    assert len(namelist) == 6
    assert 'FL-E13-5_chr1_Eya1.png' in namelist
    assert 'MB-E10-5_chr1_Eya1.png' in namelist
    assert 'FL-E13-5_chr1_Sox17.png' in namelist
    assert 'MB-E10-5_chr1_Sox17.png' in namelist
    assert 'FL-E13-5_chr1_Tfap2d.png' in namelist
    assert 'MB-E10-5_chr1_Tfap2d.png' in namelist

    output_folder = mkdtemp(prefix="output_")

    file_obj.extractall(output_folder)
    res = compare_images(
        ROOT + 'chicPlotViewpoint/background_single_signficant_pvalues/FL-E13-5_chr1_Eya1.png', output_folder + '/FL-E13-5_chr1_Eya1.png', tolerance)
    assert res is None, res
    res = compare_images(
        ROOT + 'chicPlotViewpoint/background_single_signficant_pvalues/MB-E10-5_chr1_Eya1.png', output_folder + '/MB-E10-5_chr1_Eya1.png', tolerance)
    assert res is None, res
    res = compare_images(
        ROOT + 'chicPlotViewpoint/background_single_signficant_pvalues/FL-E13-5_chr1_Sox17.png', output_folder + '/FL-E13-5_chr1_Sox17.png', tolerance)
    assert res is None, res
    res = compare_images(
        ROOT + 'chicPlotViewpoint/background_single_signficant_pvalues/MB-E10-5_chr1_Sox17.png', output_folder + '/MB-E10-5_chr1_Sox17.png', tolerance)
    assert res is None, res
    res = compare_images(
        ROOT + 'chicPlotViewpoint/background_single_signficant_pvalues/FL-E13-5_chr1_Tfap2d.png', output_folder + '/FL-E13-5_chr1_Tfap2d.png', tolerance)
    assert res is None, res
    res = compare_images(
        ROOT + 'chicPlotViewpoint/background_single_signficant_pvalues/MB-E10-5_chr1_Tfap2d.png', output_folder + '/MB-E10-5_chr1_Tfap2d.png', tolerance)

    assert res is None, res


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
def test_allGenes():
    outfile = NamedTemporaryFile(suffix='.tar.gz', delete=False)
    outfile.close()
    args = "-if {} --range {} {} -o {} -t {} --combinationMode allGenes --backgroundModelFile {} --significantInteractions {} --pValue --plotSignificantInteractions".format(
        ROOT + 'chicViewpoint/two_matrices.hdf5', 200000, 200000, outfile.name, 1, ROOT + 'background.txt', ROOT + 'chicSignificantInteractions/significantInteractions_dual.hdf5').split()
    chicPlotViewpoint.main(args)

    file_obj = tarfile.open(outfile.name, "r")

    namelist = file_obj.getnames()
    # assert len(namelist) == 3
    assert 'FL-E13-5_chr1_MB-E10-5_chr1_Eya1.png' in namelist
    assert 'FL-E13-5_chr1_MB-E10-5_chr1_Sox17.png' in namelist
    assert 'FL-E13-5_chr1_MB-E10-5_chr1_Tfap2d.png' in namelist

    output_folder = mkdtemp(prefix="output_")

    file_obj.extractall(output_folder)
    res = compare_images(
        ROOT + 'chicPlotViewpoint/allGenes/FL-E13-5_chr1_MB-E10-5_chr1_Eya1.png', output_folder + '/FL-E13-5_chr1_MB-E10-5_chr1_Eya1.png', tolerance)
    assert res is None, res
    res = compare_images(
        ROOT + 'chicPlotViewpoint/allGenes/FL-E13-5_chr1_MB-E10-5_chr1_Sox17.png', output_folder + '/FL-E13-5_chr1_MB-E10-5_chr1_Sox17.png', tolerance)
    assert res is None, res
    res = compare_images(
        ROOT + 'chicPlotViewpoint/allGenes/FL-E13-5_chr1_MB-E10-5_chr1_Tfap2d.png', output_folder + '/FL-E13-5_chr1_MB-E10-5_chr1_Tfap2d.png', tolerance)

    assert res is None, res
