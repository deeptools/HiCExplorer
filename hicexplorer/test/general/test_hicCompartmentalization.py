import warnings
import os
from matplotlib.testing.compare import compare_images
from hicexplorer import hicCompartmentalization
from tempfile import NamedTemporaryFile

warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
tolerance = 60  # default matplotlib
ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")


def test_compartmentalization():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)
    outfile.close()

    args = " -m {} --pca {} -o {} --outliers 0.0 --quantile 30".format(ROOT + "hicPCA/obsexp_norm.h5",
                                                                       ROOT + "hicCompartmentalization/pca1.bedgraph",
                                                                       outfile.name).split()
    hicCompartmentalization.main(args)
    test = ROOT + "hicCompartmentalization/compartmentalizationRatio.png"
    res = compare_images(test, outfile.name, tolerance)
    assert res is None, res

    os.unlink(outfile.name)
