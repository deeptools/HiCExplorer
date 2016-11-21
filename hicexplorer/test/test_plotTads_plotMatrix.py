import filecmp
import sys
import matplotlib as mpl
from tempfile import NamedTemporaryFile
import hicexplorer
import os.path

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"

def compare_svg(file1, file2):
    """
    svg files usually differ on randomly assigned ids and xlink:href tags
    This code compares the files ignoring the lines that contain ids

    :return: bool True if files are similar
    """
    f1 = NamedTemporaryFile(suffix='.svg', delete=False)
    f2 = NamedTemporaryFile(suffix='.svg', delete=False)
    # remove xlink:href, id and url attributes
    os.system('cat {} | perl -lane \'s/xlink:href=".+?"//g; s/id=".+?"//g; s/"url\(.+?\)"//g; print $_\' > {}'.format(file1, f1.name))
    os.system('cat {} | perl -lane \'s/xlink:href=".+?"//g; s/id=".+?"//g; s/"url\(.+?\)"//g; print $_\' > {}'.format(file2, f2.name))
    res = filecmp.cmp(f1.name, f2.name)
    os.remove(f1.name)
    os.remove(f2.name)
    return res


class TestPlottingPrograms(object):
    def __init__(self):
        self.run_image_tests = True

    def setUp(self):
        # the tests based on images were done with
        # matplotlib 1.5.3 and will fail if other
        # version is used
        if mpl.__version__ != '1.5.3':
            sys.stderr.write("\nTests based on images are skipped because of "
                             "different matplotlib version ({}) != 1.5.0\n".format(mpl.__version__))
            exit(1)
            self.run_image_tests = False
        if sys.version_info[0] != 2:
            self.run_image_tests = False

    def test_hicPlotTads(self):
        import hicexplorer.hicPlotTADs

        outfile = NamedTemporaryFile(suffix='.svg', delete=False)
        args = "--tracks {0}/browser_tracks.ini --region chrX:3000000-3500000   " \
               "--outFileName  {1}".format(ROOT, outfile.name).split()
        hicexplorer.hicPlotTADs.main(args)
        assert compare_svg(ROOT + '/master_TADs_plot.svg', outfile.name) is True
        os.remove(outfile.name) 

