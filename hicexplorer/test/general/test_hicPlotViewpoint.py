import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicPlotViewpoint
from tempfile import NamedTemporaryFile
import os

import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
import os.path


ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")


def are_files_equal(file1, file2):
    equal = True
    with open(file1) as textfile1, open(file2) as textfile2:
        for x, y in zip(textfile1, textfile2):
            if x.startswith('File'):
                continue
            if x != y:
                equal = False
                break
    return equal


def test_plot_single_point():

    outfile = NamedTemporaryFile(suffix='.png', prefix='viewpoint1', delete=False)
    matrix = ROOT + 'Li_et_al_2015.h5'
    args = "--matrix {} --region X:3000000-3500000 -rp X:3200000 --outFileName {} --dpi 300".format(matrix, outfile.name).split()
    hicPlotViewpoint.main(args)

    res = compare_images(ROOT + '/hicPlotViewpoint/li_viewpoint_32Mb.png', outfile.name, tol=40)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_single_point_two_matrices():

    outfile = NamedTemporaryFile(suffix='.png', prefix='viewpoint1', delete=False)
    matrix = ROOT + 'Li_et_al_2015.h5' + ' ' + ROOT + 'Li_et_al_2015_twice.h5'
    args = "--matrix {} --region X:3000000-3500000 -rp X:3200000 --outFileName {} --dpi 300".format(matrix, outfile.name).split()
    hicPlotViewpoint.main(args)

    res = compare_images(ROOT + '/hicPlotViewpoint/li_viewpoint_32Mb_twice.png', outfile.name, tol=40)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_single_point_interaction_file():

    outfile = NamedTemporaryFile(suffix='.png', prefix='viewpoint2', delete=False)
    outfile_interactions = NamedTemporaryFile(suffix='.bedgraph', prefix='viewpoint_interactons', delete=False)
    matrix = ROOT + 'Li_et_al_2015.h5'

    args = "--matrix {} --region X:3000000-3500000 -rp X:3200000 --outFileName {} -i {} --dpi 300".format(matrix, outfile.name, outfile_interactions.name).split()
    hicPlotViewpoint.main(args)

    res = compare_images(ROOT + '/hicPlotViewpoint/li_viewpoint_32Mb.png', outfile.name, tol=40)
    assert res is None, res
    assert are_files_equal(ROOT + '/hicPlotViewpoint/li_32mb_interactions.bedgraph', outfile_interactions.name)
    os.remove(outfile.name)
    os.remove(outfile_interactions.name)


def test_plot_single_point_interaction_file_two_matrices():

    outfile = NamedTemporaryFile(suffix='.png', prefix='viewpoint2', delete=False)
    outfile_interactions_one = NamedTemporaryFile(suffix='.bedgraph', prefix='viewpoint_interactons_Li_et_al_2015.h5', delete=False)
    outfile_interactions_two = NamedTemporaryFile(suffix='.bedgraph', prefix='viewpoint_interactons_Li_et_al_2015_twice.h5', delete=False)
    matrix = ROOT + 'Li_et_al_2015.h5' + ' ' + ROOT + 'Li_et_al_2015_twice.h5'

    args = "--matrix {} --region X:3000000-3500000 -rp X:3200000 --outFileName {} -i {} --dpi 300".format(matrix, outfile.name, 'viewpoint_interactons').split()
    hicPlotViewpoint.main(args)

    res = compare_images(ROOT + '/hicPlotViewpoint/li_viewpoint_32Mb_twice.png', outfile.name, tol=40)
    assert res is None, res
    assert are_files_equal(ROOT + '/hicPlotViewpoint/li_32mb_interactions_one.bedgraph', outfile_interactions_one.name)
    assert are_files_equal(ROOT + '/hicPlotViewpoint/li_32mb_interactions_two.bedgraph', outfile_interactions_two.name)

    os.remove(outfile.name)
    os.remove(outfile_interactions_one.name)
    os.remove(outfile_interactions_two.name)


def test_plot_region():

    outfile = NamedTemporaryFile(suffix='.png', prefix='viewpoint3', delete=False)
    matrix = ROOT + 'Li_et_al_2015.h5'

    args = "--matrix {} --region X:3000000-3500000 -rp X:3200000-3300000 --outFileName {} --dpi 300".format(matrix, outfile.name).split()
    hicPlotViewpoint.main(args)

    res = compare_images(ROOT + '/hicPlotViewpoint/li_viewpoint_32-33Mb.png', outfile.name, tol=40)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_region_interaction_file():

    outfile = NamedTemporaryFile(suffix='.png', prefix='viewpoint4', delete=False)

    outfile_interactions = NamedTemporaryFile(suffix='.bedgraph', prefix='viewpoint_interactons', delete=False)
    matrix = ROOT + 'Li_et_al_2015.h5'

    args = "--matrix {} --region X:3000000-3500000 -rp X:3200000-3300000  --outFileName {} -i {} --dpi 300".format(matrix, outfile.name, outfile_interactions.name).split()
    hicPlotViewpoint.main(args)

    res = compare_images(ROOT + 'hicPlotViewpoint/li_viewpoint_32-33Mb.png', outfile.name, tol=40)
    assert res is None, res
    assert are_files_equal(ROOT + 'hicPlotViewpoint/li_32-33mb_interactions.bedgraph', outfile_interactions.name)
    os.remove(outfile.name)
    os.remove(outfile_interactions.name)
