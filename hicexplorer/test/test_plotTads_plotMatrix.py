import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
import os.path
import pytest
import sys

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"

h5_browser_track = """

[x-axis]

[hic Kc]
file = Li_et_al_2015.h5
title = Kc DpnII (Li et al. 2015)
colormap = RdYlBu_r
depth = 200000
transform = log1p
x labels = yes
boundaries_file = domains.bed

"""

cool_browser_track = """
[x-axis]

[hic Kc]
file = Li_et_al_2015.cool
title = Kc DpnII (Li et al. 2015)
colormap = RdYlBu_r
depth = 200000
transform = log1p
x labels = yes
boundaries_file = domains.bed

"""

browser_tracks = """

[spacer]
width = 0.05

[tad state]
file = tad_classification.bed
width = 0.5
title = TAD state
display = collapsed
labels = off

[tad-separation score]
file = tad_score.gz
width = 10
type = lines
title= TAD separation score (Ramirez et al.)
file_type = bedgraph_matrix


[spacer]
width = 1

[test bedgraph]
file = bedgraph_chrx_2e6_5e6.bg
color = blue
width = 4
title = bedgraph

[test bigwig]
file = bigwig_chrx_2e6_5e6.bw
color = blue
width = 4
title = rep 1 test fill

[test bigwig lines]
file = bigwig_chrx_2e6_5e6.bw
color = red
width = 4
type = line
title = rep 1 test line

[test bigwig lines]
file = bigwig_chrx_2e6_5e6.bw
color = red
width = 4
type = line:0.2
title = rep 1 test lw=0.1

[test bigwig points]
file = bigwig_chrx_2e6_5e6.bw
color = black
width = 4
type = points:0.5
title = rep 1 test point:0.5

[spacer]
width = 0.5

[genes 2]
file = dm3_genes.bed.gz
width = 5
title = genes
fontsize = 10

[spacer]
width = 1

[test gene rows]
file = dm3_genes.bed.gz
width = 3
title = max num rows 3
fontsize = 8
gene rows = 3

[spacer]
width = 1

[test bed6]
file = dm3_genes.bed6.gz
width = 20
title = bed6 global max row
fontsize = 10
file_type = bed
global max row = yes

[vlines]
file = domains.bed
type = vlines

"""
from tempfile import NamedTemporaryFile
from sys import platform
from sys import version_info
with open(ROOT + "browser_tracks.ini", 'w') as fh:
    fh.write(h5_browser_track)
    fh.write(browser_tracks)

with open(ROOT + "browser_tracks_cool.ini", 'w') as fh:
    fh.write(cool_browser_track)
    fh.write(browser_tracks)

tolerance = 13  # default matplotlib pixed difference tolerance

@pytest.mark.skipif(sys.platform == 'darwin' and sys.version_info[0] == 3,
                    reason="Travis has too less memory to run it.")
def test_hicPlotTads():
    import hicexplorer.hicPlotTADs

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_h5_', delete=False)

    if platform == "darwin" and version_info[0] == 3:
        args = "--tracks {0}/browser_tracks.ini --region chrX:300000-350000  " \
            "--outFileName  {1}".format(ROOT, outfile.name).split()
        test_image_path = ROOT + '/plot_matrix_python3_osx/tads_plot.png'
    else:
        args = "--tracks {0}/browser_tracks.ini --region chrX:3000000-3500000  " \
            "--outFileName  {1}".format(ROOT, outfile.name).split()
        test_image_path = ROOT + '/master_TADs_plot.png'

    hicexplorer.hicPlotTADs.main(args)

    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)

@pytest.mark.skipif(sys.platform == 'darwin' and sys.version_info[0] == 3,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix():
    import hicexplorer.hicPlotMatrix

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_h5_', delete=False)

    if platform == "darwin" and version_info[0] == 3:
        args = "--matrix {0}/Li_et_al_2015.h5 --region chrX:300000-350000 --region2 chrX:310000-360000 " \
            "--outFileName  {1} --log1p --clearMaskedBins".format(ROOT, outfile.name).split()
        test_image_path = ROOT + '/plot_matrix_python3_osx/matrix_plot.png'

    else:
        args = "--matrix {0}/Li_et_al_2015.h5 --region chrX:3000000-3500000 --region2 chrX:3100000-3600000 " \
            "--outFileName  {1} --log1p --clearMaskedBins".format(ROOT, outfile.name).split()
        test_image_path = ROOT + '/master_matrix_plot.png'

    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


# add test cases for all parameters

# # chromosome order
    # outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool', delete=False)
    # args = "--matrix {0}/Li_et_al_2015.cool" \
    #        "--outFileName  {1}".format(ROOT, outfile.name).split()
    # hicexplorer.hicPlotMatrix.main(args)
    # res = compare_images(ROOT + '/master_matrix_plot_cool_chromosomeOrder.png', outfile.name, tolerance)
    # assert res is None, res
    # os.remove(outfile.name)

    # # log1p
    # outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool', delete=False)
    # args = "--matrix {0}/Li_et_al_2015.cool" \
    #        "--outFileName  {1}".format(ROOT, outfile.name).split()
    # hicexplorer.hicPlotMatrix.main(args)
    # res = compare_images(ROOT + '/master_matrix_plot_cool_log1p.png', outfile.name, tolerance)
    # assert res is None, res
    # os.remove(outfile.name)

    # # log
    # outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool', delete=False)
    # args = "--matrix {0}/Li_et_al_2015.cool" \
    #        "--outFileName  {1}".format(ROOT, outfile.name).split()
    # hicexplorer.hicPlotMatrix.main(args)
    # res = compare_images(ROOT + '/master_matrix_plot_cool_log.png', outfile.name, tolerance)
    # assert res is None, res
    # os.remove(outfile.name)

    # # color map
    # outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool', delete=False)
    # args = "--matrix {0}/Li_et_al_2015.cool" \
    #        "--outFileName  {1}".format(ROOT, outfile.name).split()
    # hicexplorer.hicPlotMatrix.main(args)
    # res = compare_images(ROOT + '/master_matrix_plot_cool_colormap.png', outfile.name, tolerance)
    # assert res is None, res
    # os.remove(outfile.name)

@pytest.mark.skipif(sys.platform == 'darwin' and sys.version_info[0] == 3,
                    reason="Travis has too less memory to run it.")
def test_hicPlotTads_cool():
    import hicexplorer.hicPlotTADs

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool_', delete=False)
    args = "--tracks {0}/browser_tracks_cool.ini --region X:3000000-3500000  " \
           "--outFileName  {1}".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotTADs.main(args)

    res = compare_images(ROOT + '/master_TADs_plot_cool_partial_load.png', outfile.name, tol=40)
    assert res is None, res

    os.remove(outfile.name)

@pytest.mark.skipif(sys.platform == 'darwin' and sys.version_info[0] == 3,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_cool():
    import hicexplorer.hicPlotMatrix

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool', delete=False)

    args = "--matrix {0}/Li_et_al_2015.cool --region chrX:3000000-3500000 --region2 chrX:3100000-3600000 " \
           "--outFileName  {1} --log1p --clearMaskedBins".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotMatrix.main(args)
    res = compare_images(ROOT + '/master_matrix_plot.png', outfile.name, tol=40)
    assert res is None, res
    os.remove(outfile.name)
