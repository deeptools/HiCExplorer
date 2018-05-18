"""
    Testsuite for testing if all components of hiCexplorer are working without raising
    Exceptions.

    This means all tests are designed to fail if there occur any Exceptions in trivial
    calls of the code using the argparser`s arguments.
"""

import pytest
from tempfile import NamedTemporaryFile, mkdtemp
import os
from psutil import virtual_memory

from hicexplorer.utilities import genomicRegion
from hicexplorer import hicBuildMatrix as hicBuildMatrix
import hicexplorer.hicAggregateContacts

mem = virtual_memory()
memory = mem.total / 2**30

# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 7
HIGH_MEMORY = 200

REMOVE_OUTPUT = True


# Some definitions needed for tests
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
# test_build_matrix
sam_R1 = ROOT + "small_test_R1_unsorted.bam"
sam_R2 = ROOT + "small_test_R2_unsorted.bam"
dpnii_file = ROOT + "DpnII.bed"
outFile = NamedTemporaryFile(suffix='.h5', delete=False)
qc_folder = mkdtemp(prefix="testQC_")
# test_AggregateContacts
matrix = ROOT + 'Li_et_al_2015.h5'
BED = ROOT + 'hicAggregateContacts/test_regions.bed'
BED2 = ROOT + 'hicAggregateContacts/test_regions.bed'
outfile_aggregate_plots = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_test_', delete=False)
diagnosticHeatmapFile = NamedTemporaryFile(suffix='.png', prefix='hicaggregate_heatmap_', delete=False)


@pytest.mark.parametrize("sam1", [sam_R1])  # required
@pytest.mark.parametrize("sam2", [sam_R2])  # required
@pytest.mark.parametrize("outFile", [outFile])  # required
@pytest.mark.parametrize("qcFolder", [qc_folder])  # required
@pytest.mark.parametrize("outBam", ['/tmp/test.bam'])
@pytest.mark.parametrize("binSize", [5000])  # required | restrictionCutFile
@pytest.mark.parametrize("restrictionCutFile", [dpnii_file])  # required | binSize
@pytest.mark.parametrize("minDistance", [150])
@pytest.mark.parametrize("maxDistance", [1500])
@pytest.mark.parametrize("maxLibraryInsertSize", [1500])
@pytest.mark.parametrize("restrictionSequence", ['GATC'])
@pytest.mark.parametrize("danglingSequence", ['GATC'])
@pytest.mark.parametrize("region", ["ChrX"])  # region does not work!!
@pytest.mark.parametrize("removeSelfLigation", [True])
@pytest.mark.parametrize("minMappingQuality", [15])
@pytest.mark.parametrize("threads", [4])
@pytest.mark.parametrize("inputBufferSize", [400000])
def test_build_matrix(sam1, sam2, outFile, qcFolder, outBam, binSize, restrictionCutFile,
                      minDistance, maxDistance, maxLibraryInsertSize, restrictionSequence,
                      danglingSequence, region, removeSelfLigation,
                      minMappingQuality, threads, inputBufferSize):
    """
    This test runs buildMatrix with all command line args for one time to ensure all args
    being ok.

    Note: Test will take some time.

    Note: parameters can be expanded (the values in the list) so that many combinations of
          command line args can be tested.
    """
    # test required arguments
    # IMPORTANT: one of binSize or restrictionCutFile is required!
    # test restrictionCutFile
    args = "-s {} {} --restrictionCutFile {} --outFileName {} " \
           "--QCfolder {} ".format(sam1, sam2,
                                   restrictionCutFile,
                                   outFile.name,
                                   qcFolder).split()

    hicBuildMatrix.main(args)

    # test binSize instead of restrictionCutFile
    args = "-s {} {} --binSize {} --outFileName {}  --QCfolder {} ".format(sam1, sam2,
                                                                           binSize,
                                                                           outFile.name,
                                                                           qcFolder).split()

    hicBuildMatrix.main(args)

    # test more args for restrictionCutFile option
    # Note: Something wrong with region argument, test fails due to region param.
    with pytest.raises(TypeError):
        region = genomicRegion(region)

    args = "-s {} {} --restrictionCutFile {} --outFileName {} --QCfolder {} " \
           "--restrictionSequence {} " \
           "--danglingSequence {} " \
           "--minDistance {} " \
           "--maxLibraryInsertSize {} --threads {} " \
           "--region {} --removeSelfLigation {} ".format(sam_R1, sam_R2,
                                                         restrictionCutFile, outFile.name,
                                                         qcFolder,
                                                         restrictionSequence,
                                                         danglingSequence,
                                                         minDistance,
                                                         maxLibraryInsertSize,
                                                         threads, region,
                                                         removeSelfLigation).split()
    with pytest.raises(SystemExit):
        hicBuildMatrix.main(args)

    # test more params with restrictionCutFile (now without region param)
    args = "-s {} {} --restrictionCutFile {} --outFileName {} --QCfolder {} " \
           "--restrictionSequence {} " \
           "--danglingSequence {} " \
           "--minDistance {} " \
           "--maxLibraryInsertSize {} --threads {} " \
           "--removeSelfLigation {} ".format(sam_R1, sam_R2,
                                             restrictionCutFile, outFile.name, qcFolder,
                                             restrictionSequence, danglingSequence,
                                             minDistance, maxLibraryInsertSize, threads,
                                             removeSelfLigation).split()

    hicBuildMatrix.main(args)

    # test more params with restrictionCutFile (now without region param)
    args = "-s {} {} --restrictionCutFile {} --outFileName {} --QCfolder {} " \
           "--restrictionSequence {} " \
           "--danglingSequence {} " \
           "--minDistance {} " \
           "--maxLibraryInsertSize {} --threads {} " \
           "--removeSelfLigation {} --keepSelfCircles ".format(sam_R1, sam_R2,
                                                               restrictionCutFile,
                                                               outFile.name, qcFolder,
                                                               restrictionSequence,
                                                               danglingSequence,
                                                               minDistance,
                                                               maxLibraryInsertSize,
                                                               threads,
                                                               removeSelfLigation).split()

    hicBuildMatrix.main(args)

    # added minMappingQuality
    args = "-s {} {} --restrictionCutFile {} --outFileName {} --QCfolder {} " \
           "--restrictionSequence {} " \
           "--danglingSequence {} " \
           "--minDistance {} " \
           "--maxLibraryInsertSize {} --threads {} " \
           "--removeSelfLigation {} --keepSelfCircles " \
           "--minMappingQuality {} ".format(sam_R1, sam_R2,
                                            restrictionCutFile, outFile.name, qcFolder,
                                            restrictionSequence, danglingSequence,
                                            minDistance, maxLibraryInsertSize, threads,
                                            removeSelfLigation, minMappingQuality).split()

    hicBuildMatrix.main(args)

    # added inputBufferSize
    args = "-s {} {} --restrictionCutFile {} --outFileName {} --QCfolder {} " \
           "--restrictionSequence {} " \
           "--danglingSequence {} " \
           "--minDistance {} " \
           "--maxLibraryInsertSize {} --threads {} " \
           "--removeSelfLigation {} --keepSelfCircles " \
           "--minMappingQuality {} --inputBufferSize {} ".format(sam_R1, sam_R2,
                                                                 restrictionCutFile,
                                                                 outFile.name, qcFolder,
                                                                 restrictionSequence,
                                                                 danglingSequence,
                                                                 minDistance,
                                                                 maxLibraryInsertSize,
                                                                 threads,
                                                                 removeSelfLigation,
                                                                 minMappingQuality,
                                                                 inputBufferSize).split()

    hicBuildMatrix.main(args)

    # added doTestRun
    args = "-s {} {} --restrictionCutFile {} --outFileName {} --QCfolder {} " \
           "--restrictionSequence {} " \
           "--danglingSequence {} " \
           "--minDistance {} " \
           "--maxLibraryInsertSize {} --threads {} " \
           "--removeSelfLigation {} --keepSelfCircles " \
           "--minMappingQuality {} --inputBufferSize {} " \
           "--doTestRun ".format(sam_R1, sam_R2,
                                 restrictionCutFile, outFile.name, qcFolder,
                                 restrictionSequence, danglingSequence,
                                 minDistance, maxLibraryInsertSize, threads,
                                 removeSelfLigation, minMappingQuality,
                                 inputBufferSize).split()

    hicBuildMatrix.main(args)

    # added skipDuplicationCheck
    args = "-s {} {} --restrictionCutFile {} --outFileName {} --QCfolder {} " \
           "--restrictionSequence {} " \
           "--danglingSequence {} " \
           "--minDistance {} " \
           "--maxLibraryInsertSize {} --threads {} " \
           "--removeSelfLigation {} --keepSelfCircles " \
           "--minMappingQuality {} --inputBufferSize {} " \
           "--doTestRun --skipDuplicationCheck ".format(sam_R1, sam_R2,
                                                        restrictionCutFile, outFile.name,
                                                        qcFolder, restrictionSequence,
                                                        danglingSequence, minDistance,
                                                        maxLibraryInsertSize, threads,
                                                        removeSelfLigation,
                                                        minMappingQuality,
                                                        inputBufferSize).split()

    hicBuildMatrix.main(args)


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

    # test outFileContactPairs^
    args = "--matrix {} --outFileName {} --BED {} --range {} --BED2 {} " \
           "--numberOfBins {} --transform {} --avgType {} --outFileContactPairs {} " \
           "--kmeans {} --hclust {} " \
           "--howToCluster {} --chromosomes {} --colorMap {} --plotType {} --vMin {} " \
           "--vMax {} --disable_bbox_tight".format(matrix, outFileName.name, BED, ran,
                                                   BED2, numberOfBins, transform, avgType,
                                                   outFileContactPairs,
                                                   kmeans, hclust,
                                                   howToCluster, chromosomes, colorMap,
                                                   plotType, vMin, vMax).split()
    hicexplorer.hicAggregateContacts.main(args)
    os.remove(outFileName.name)

    # test diagnosticHeatmapFile
    # first test with all parameters failed due to unknown error.
    args = "--matrix {} --BED {} " \
           "--outFileName {out_agg} --numberOfBins 30 --range 50000:900000 --hclust 4 " \
           "--diagnosticHeatmapFile {out_heat} --howToCluster diagonal  --disable_bbox_tight " \
           "--BED2 {}".format(matrix, BED, BED2, out_agg=outFileName.name,
                              out_heat=diagnosticHeatmapFile.name)

    hicexplorer.hicAggregateContacts.main(args.split())
