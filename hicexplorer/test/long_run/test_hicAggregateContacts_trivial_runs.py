import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
from tempfile import NamedTemporaryFile
import os
from psutil import virtual_memory

import hicexplorer.hicAggregateContacts

mem = virtual_memory()
memory = mem.total / 2**30

# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 7
HIGH_MEMORY = 200

REMOVE_OUTPUT = True


# Some definitions needed for tests
ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")
# test_AggregateContacts
matrix = ROOT + 'Li_et_al_2015.h5'
# matrix = ROOT + 'R1_R2_1000.h5'
BED = ROOT + 'hicAggregateContacts/test_regions.bed'
BED2 = ROOT + 'hicAggregateContacts/test_regions.bed'
outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)
diagnosticHeatmapFile = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_heatmap_', delete=False)


@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
@pytest.mark.parametrize("matrix", [matrix])  # required
@pytest.mark.parametrize("outFileName", [outfile_aggregate_plots])  # required
@pytest.mark.parametrize("BED", [BED])  # required
@pytest.mark.parametrize("ran", ['50000:900000'])  # required
@pytest.mark.parametrize("BED2", [BED2])
@pytest.mark.parametrize("numberOfBins", [30])
@pytest.mark.parametrize("transform", ['total-counts', 'z-score', 'obs/exp', 'none'])
@pytest.mark.parametrize("avgType", ['mean', 'median'])
@pytest.mark.parametrize("outFilePrefixMatrix", ['outFilePrefix'])
@pytest.mark.parametrize("outFileContactPairs", ['outFileContactPairs'])
@pytest.mark.parametrize("diagnosticHeatmapFile", [diagnosticHeatmapFile])
@pytest.mark.parametrize("kmeans", [4])
@pytest.mark.parametrize("hclust", [4])
@pytest.mark.parametrize("howToCluster", ['full', 'center', 'diagonal'])
@pytest.mark.parametrize("chromosomes", ['X'])
@pytest.mark.parametrize("colorMap", ['RdYlBu_r'])
@pytest.mark.parametrize("plotType", ['2d', '3d'])
@pytest.mark.parametrize("vMin", [0.01])
@pytest.mark.parametrize("vMax", [1.0])
def test_aggregate_contacts(capsys, matrix, outFileName, BED, ran, BED2, numberOfBins, transform,
                            avgType, outFilePrefixMatrix, outFileContactPairs,
                            diagnosticHeatmapFile, kmeans, hclust, howToCluster,
                            chromosomes, colorMap, plotType, vMin, vMax):
    """
        Test will run all configurations defined by the parametrized option.
    """
    # test outFilePrefixMatrix
    args = "--matrix {} --outFileName {} --BED {} --range {} --BED2 {} " \
           "--numberOfBins {} --transform {} --avgType {} --outFilePrefixMatrix {} " \
           "--kmeans {} --hclust {} " \
           "--howToCluster {} --chromosomes {} --colorMap {} --plotType {} --vMin {} " \
           "--vMax {} --disable_bbox_tight".format(matrix, outFileName.name, BED, ran,
                                                   BED2, numberOfBins, transform, avgType,
                                                   outFilePrefixMatrix,
                                                   kmeans, hclust,
                                                   howToCluster, chromosomes, colorMap,
                                                   plotType, vMin, vMax).split()
    hicexplorer.hicAggregateContacts.main(args)
    os.remove(outFileName.name)
