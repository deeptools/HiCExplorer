import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from tempfile import NamedTemporaryFile

import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
import os.path
import pytest
from psutil import virtual_memory
mem = virtual_memory()
memory = mem.total / 2 ** 30


import hicexplorer.hicPlotMatrix
tolerance = 30  # default matplotlib pixed difference tolerance
ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")

# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 7
HIGH_MEMORY = 200

REMOVE_OUTPUT = True
# DIFF = 60


@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_region_region2_log1p_clearMaskedBins_and_bigwig():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_h5_', delete=False)

    args = "--matrix {0}/Li_et_al_2015.h5 --region chrX:3000000-3500000 --region2 chrX:3100000-3600000 " \
        "--outFileName  {1} --log1p --clearMaskedBins --bigwig {2} ".format(ROOT, outfile.name,
                                                                            ROOT + "bigwig_chrx_2e6_5e6.bw").split()
    test_image_path = ROOT + "hicPlotMatrix" + '/Li_chrX30-35-chrX31-36_log1p_clearmaskedbins.png'

    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_region_region2_log1p_clearMaskedBins_and_bigwig_vmin_vmax():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_h5_', delete=False)

    args = "--matrix {0}/Li_et_al_2015.h5 --region chrX:3000000-3500000 --region2 chrX:3100000-3600000 " \
        "--outFileName  {1} --log1p --clearMaskedBins --bigwig {2} --vMinBigwig {3} --vMaxBigwig {4}".format(ROOT, outfile.name,
                                                                                                             ROOT + "bigwig_chrx_2e6_5e6.bw", 0, 1).split()
    test_image_path = ROOT + "hicPlotMatrix" + '/Li_chrX30-35-chrX31-36_log1p_clearmaskedbins_vbigwigmin_vbigwigmax.png'

    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

#     if REMOVE_OUTPUT:
#         os.remove(outfile.name)


@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_region_region2_log_no_clearMaskedBins_and_bigwig():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_h5_', delete=False)

    args = "--matrix {0}/Li_et_al_2015.h5 --region chrX:3000000-3500000 --region2 chrX:3100000-3600000 " \
        "--outFileName  {1} --log --bigwig {2} ".format(ROOT, outfile.name,
                                                        ROOT + "bigwig_chrx_2e6_5e6.bw").split()
    test_image_path = ROOT + "hicPlotMatrix" + '/Li_chrX30-35-chrX31-36_log_no_clearmasked.png'

    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_region_region2_no_clearMaskedBins():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_h5_', delete=False)

    args = "--matrix {0}/Li_et_al_2015.h5 --region chrX:3000000-3500000 --region2 chrX:3100000-3600000 " \
        "--outFileName  {1} --clearMaskedBins".format(ROOT, outfile.name).split()
    test_image_path = ROOT + "hicPlotMatrix" + '/Li_chrX30-35-chrX31-36_no_clearmasked.png'

    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_region_region2_no_clearMaskedBins_title():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_h5_', delete=False)
    title = 'Li_chrX:3000000-3500000_chrX:3100000-3600000'
    args = "--matrix {0}/Li_et_al_2015.h5 --region chrX:3000000-3500000 --region2 chrX:3100000-3600000 " \
        "--outFileName  {1} --clearMaskedBins --title {2}".format(ROOT, outfile.name, title).split()
    test_image_path = ROOT + "hicPlotMatrix" + '/Li_chrX30-35-chrX31-36_no_clearmasked_title.png'

    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_cool_region1_region2():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool', delete=False)

    args = "--matrix {0}/Li_et_al_2015.cool --region chrX:3000000-3500000 --region2 chrX:3100000-3600000 " \
           "--outFileName  {1} ".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/Li_chrX30-35-chrX31-36_cool.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(LOW_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_cool_region1():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool', delete=False)

    args = "--matrix {0}/Li_et_al_2015.cool --region X:3000000-3500000 " \
           "--outFileName  {1} ".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/Li_chrX30-35_cool.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(LOW_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_h5_region1():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool', delete=False)

    args = "--matrix {0}/Li_et_al_2015.h5 --region X:3000000-3500000 " \
           "--outFileName  {1} ".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/Li_chrX30-35_cool.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(HIGH_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_cool_log1p():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool', delete=False)

    args = "--matrix {0}/Li_et_al_2015.cool --log1p " \
           "--outFileName  {1} ".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/Li_cool_log1p.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(HIGH_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_cool_log():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool', delete=False)

    args = "--matrix {0}/Li_et_al_2015.cool --log " \
           "--outFileName  {1} ".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/Li_cool_log.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(HIGH_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_cool_full():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool', delete=False)

    args = "--matrix {0}/Li_et_al_2015.cool " \
           "--outFileName  {1} ".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/Li_cool.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(HIGH_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_h5_log1p():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_h5', delete=False)

    args = "--matrix {0}/Li_et_al_2015.h5 --log1p " \
           "--outFileName  {1} ".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/Li_h5_log1p.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(HIGH_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_h5_log():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_h5', delete=False)

    args = "--matrix {0}/Li_et_al_2015.h5 --log " \
           "--outFileName  {1} ".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/Li_h5_log.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(HIGH_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_h5_full():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_h5', delete=False)

    args = "--matrix {0}/Li_et_al_2015.h5 " \
           "--outFileName  {1} ".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/Li_h5.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(LOW_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_cool_log_region1_region2():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool', delete=False)

    args = "--matrix {0}/Li_et_al_2015.cool --region chrX:3000000-3500000 --region2 chrX:3100000-3600000 " \
           "--outFileName  {1} --log ".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/Li_chrX30-35-chrX31-36_cool_log.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


# @pytest.mark.xfail
@pytest.mark.skipif(LOW_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_cool_log_region1_region2_without_cool_suffix():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool', delete=False)

    args = "--matrix {0}/Li_et_al_2015 --region chrX:3000000-3500000 --region2 chrX:3100000-3600000 " \
           "--outFileName  {1} --log ".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/Li_chrX30-35-chrX31-36_cool_log.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(LOW_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_cool_log1p_region1_region2():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool', delete=False)

    args = "--matrix {0}/Li_et_al_2015.cool --region chrX:3000000-3500000 --region2 chrX:3100000-3600000 " \
           "--outFileName  {1} --log1p ".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/Li_chrX30-35-chrX31-36_cool_log1p.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(LOW_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_perChr():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test', delete=False)

    args = "--matrix {0}/small_test_matrix_50kb_res.h5 --perChr --disable_tight_layout " \
           "--outFileName  {1} ".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/small_test_matrix_50kb_res_perChr.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(LOW_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_cool_perChr_log1p():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test', delete=False)

    args = "--matrix {0}/small_test_matrix_50kb_res.h5 --perChr  --disable_tight_layout " \
           "--outFileName  {1} --log1 --vMax 10 ".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/small_test_matrix_50kb_res_perChr_log.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


# @pytest.mark.xfail
@pytest.mark.skipif(LOW_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_h5_perChr_log1p_chromosomeOrder():
    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test', delete=False)

    args = "--matrix {0}/small_test_matrix_50kb_res.h5 --perChr  --disable_tight_layout " \
           "--outFileName  {1} --log --chromosomeOrder chr2L chr3L chr3R chr2R".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/small_test_matrix_perChr_log1p_chromosomeOrder_disable_tight_layout.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(LOW_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_cool_perChr_log1p_chromosomeOrder():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test', delete=False)

    args = "--matrix {0}/small_test_matrix_50kb_res.cool --perChr " \
           "--outFileName  {1} --log1p --chromosomeOrder chr2L chr3L chr3R chr2R".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/small_test_matrix_perChr_log1p_chromosomeOrder.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(LOW_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_perChr_pca1_bigwig():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test', delete=False)

    args = "--matrix {0}/hicTransform/pearson_perChromosome.h5 --perChr  --disable_tight_layout " \
           "--outFileName  {1} --bigwig {2}".format(ROOT, outfile.name, ROOT + "hicPCA/pca1.bw").split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/small_matrix_50kb_pearson_pca1_plot.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(LOW_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_region_pca1_colormap_bigwig():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test', delete=False)

    args = "--matrix {0}/hicTransform/pearson_perChromosome.h5 --region chr2L " \
           "--outFileName  {1} --bigwig {2} --colorMap hot".format(ROOT, outfile.name, ROOT + "hicPCA/pca1.bw").split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/small_test_50kb_pearson_pca1_plot_region__colormap_hot_chr2L_bw.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)


@pytest.mark.skipif(LOW_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_region_start_end_pca1_colormap_bigwig():

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test', delete=False)

    args = "--matrix {0}/hicTransform/pearson_perChromosome.h5 --region chr2L:15000000-20000000 " \
           "--outFileName  {1} --bigwig {2} --colorMap hot".format(ROOT, outfile.name, ROOT + "hicPCA/pca1.bw").split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + "hicPlotMatrix" + '/small_test_50kb_pearson_pca1_plot_region__colormap_hot_chr2L_15mb-20mb_bw.png', outfile.name, tol=tolerance)
    assert res is None, res
    if REMOVE_OUTPUT:
        os.remove(outfile.name)
