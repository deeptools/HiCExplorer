import warnings
import os
from matplotlib.testing.compare import compare_images
from hicexplorer import hicCompartmentsPolarization
from tempfile import NamedTemporaryFile

warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
tolerance = 60  # default matplotlib
ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")


def test_compartments_polarization():
    outfile = NamedTemporaryFile(suffix='.png', delete=False)
    outfile.close()

    args = " -m {} --pca {} -o {} ".format(ROOT + "hicPCA/obsexp_norm.h5",
                                           ROOT + "hicCompartmentsPolarization/pca1.bedgraph",
                                           outfile.name).split()
    hicCompartmentsPolarization.main(args)
    test = ROOT + "hicCompartmentsPolarization/compartmentsPolarizationRatio.png"
    res = compare_images(test, outfile.name, tolerance)
    assert res is None, res

    os.unlink(outfile.name)


def main(args=None):
    test_compartments_polarization()


if __name__ == '__main__':
    main()
