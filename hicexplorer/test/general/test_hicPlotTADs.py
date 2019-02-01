import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
import os.path
import pytest
ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")

h5_browser_track = """

[x-axis]

[hic Kc]
file = Li_et_al_2015.h5
title = Kc DpnII (Li et al. 2015)
colormap = RdYlBu_r
depth = 200000
transform = log1p
x labels = yes
file_type = hic_matrix
height = 5
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
file_type = hic_matrix
height = 5
"""

browser_tracks = """
[tads]
file = domains.bed
file_type = domains
border color = black
color = none
height = 5
line width = 1.5
overlay previous = share-y
show data range = no

[spacer]
height = 0.05

[bed]
file = tad_classification.bed
file_type = bed
height = 0.5
title = TAD state
display = collapsed
labels = off

[tad-separation score]
file = tad_score.gz
height = 5
type = lines
plot horizontal lines=False
title= TAD separation score (Ramirez et al.)
file_type = bedgraph_matrix


[spacer]
height = 1

[test bedgraph]
file = bedgraph_chrx_2e6_5e6.bg
color = blue
height = 2
title = bedgraph

[test bigwig]
file = bigwig_chrx_2e6_5e6.bw
color = blue
height = 2
title = rep 1 test fill

[test bigwig lines]
file = bigwig_chrx_2e6_5e6.bw
color = red
height = 2
type = line
title = rep 1 test line

[test bigwig lines]
file = bigwig_chrx_2e6_5e6.bw
color = red
height = 2
type = line:0.2
title = rep 1 test lw=0.1

[test bigwig points]
file = bigwig_chrx_2e6_5e6.bw
color = black
height = 2
type = points:0.5
title = rep 1 test point:0.5

[spacer]
height = 0.5

[genes 2]
file = dm3_genes.bed.gz
height = 3
title = genes
fontsize = 10

[spacer]
height = 1

[test gene rows]
file = dm3_genes.bed.gz
height = 1.5
title = max num rows 3
fontsize = 8
gene rows = 3

[spacer]
height = 1

[test bed6]
file = dm3_genes.bed6.gz
height = 10
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


@pytest.mark.skipif(platform == 'darwin' and version_info[0] == 3,
                    reason="Travis has too less memory to run it.")
def test_hicPlotTads():
    import hicexplorer.hicPlotTADs

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_h5_', delete=False)

    args = "--tracks {0}/browser_tracks.ini --region chrX:3000000-3500000  " \
        "--outFileName  {1}".format(ROOT, outfile.name).split()
    test_image_path = ROOT + '/hicPlotTADs/pygenometracks.png'

    hicexplorer.hicPlotTADs.main(args)

    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


@pytest.mark.skipif(platform == 'darwin' and version_info[0] == 3,
                    reason="Travis has too less memory to run it.")
def test_hicPlotTads_cool():
    import hicexplorer.hicPlotTADs

    outfile = NamedTemporaryFile(suffix='.png', prefix='hicexplorer_test_cool_', delete=False)
    args = "--tracks {0}/browser_tracks_cool.ini --region X:3000000-3500000  " \
           "--outFileName  {1}".format(ROOT, outfile.name).split()
    hicexplorer.hicPlotTADs.main(args)

    res = compare_images(ROOT + '/hicPlotTADs/pygenometracks.png', outfile.name, tol=40)
    assert res is None, res

    os.remove(outfile.name)
