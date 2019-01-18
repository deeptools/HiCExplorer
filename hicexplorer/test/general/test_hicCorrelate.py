import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicCorrelate
from tempfile import NamedTemporaryFile
import os
from matplotlib.testing.compare import compare_images


ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")


def test_correlate():
    outfile_heatmap = NamedTemporaryFile(suffix='heatmap.png', prefix='hicexplorer_test', delete=False)
    outfile_scatter = NamedTemporaryFile(suffix='scatter.png', prefix='hicexplorer_test', delete=False)

    args = "--matrices {} {} --labels 'first' 'second' " \
        " --method spearman --log1p --colorMap jet "\
        "--outFileNameHeatmap {} --outFileNameScatter {}".format(ROOT + "hicCorrectMatrix/small_test_matrix_corrected_chrUextra_chr3LHet.h5",
                                                                 ROOT + "hicCorrectMatrix/small_test_matrix_corrected_chrUextra_chr3LHet.h5",
                                                                 outfile_heatmap.name, outfile_scatter.name).split()
    hicCorrelate.main(args)

    res = compare_images(ROOT + "hicCorrelate" + '/heatmap.png', outfile_heatmap.name, tol=40)
    assert res is None, res

    res = compare_images(ROOT + "hicCorrelate" + '/scatter.png', outfile_scatter.name, tol=40)
    assert res is None, res
    os.remove(outfile_heatmap.name)
    os.remove(outfile_scatter.name)


def test_correlate_chromosomes():
    outfile_heatmap = NamedTemporaryFile(suffix='heatmap.png', prefix='hicexplorer_test', delete=False)
    outfile_scatter = NamedTemporaryFile(suffix='scatter.png', prefix='hicexplorer_test', delete=False)

    args = "--matrices {} {} --labels 'first' 'second' " \
        " --method spearman --log1p --colorMap jet "\
        "--outFileNameHeatmap {} --outFileNameScatter {} " \
        "--chromosomes chrUextra chr3LHet".format(ROOT + "hicCorrectMatrix/small_test_matrix_corrected_chrUextra_chr3LHet.h5",
                                                  ROOT + "hicCorrectMatrix/small_test_matrix_corrected_chrUextra_chr3LHet.h5",
                                                  outfile_heatmap.name, outfile_scatter.name).split()
    hicCorrelate.main(args)

    res = compare_images(ROOT + "hicCorrelate" + '/heatmap_chrom.png', outfile_heatmap.name, tol=40)
    assert res is None, res

    res = compare_images(ROOT + "hicCorrelate" + '/scatter_chrom.png', outfile_scatter.name, tol=40)
    assert res is None, res
    os.remove(outfile_heatmap.name)
    os.remove(outfile_scatter.name)
