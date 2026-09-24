# Tools

Tools grouped by task, each with its command-line reference.

## Pre-processing

| Tool | Description |
|---|---|
| [hicFindRestSite](hicFindRestSite.md) | Identifies the genomic locations of restriction sites. |
| [hicBuildMatrix](hicBuildMatrix.md) | Creates a Hi-C matrix from the aligned BAM/SAM files (or a .pairs file) of the Hi-C sequencing reads. |
| [hicBuildMatrixMicroC](hicBuildMatrixMicroC.md) | Creates a Hi-C matrix from Micro-C data (hicBuildMatrix without the restriction-site machinery). |
| [hicSumMatrices](hicSumMatrices.md) | Adds Hi-C matrices of the same size. |
| [hicMergeMatrixBins](hicMergeMatrixBins.md) | Merges consecutive bins on a Hi-C matrix to reduce resolution. |
| [hicCorrectMatrix](hicCorrectMatrix.md) | Uses iterative correction (ICE) or Knight-Ruiz to remove biases from a Hi-C matrix. |
| [hicNormalize](hicNormalize.md) | Normalizes matrices to a 0-1 range or to the smallest total read count. |

## Quality control

| Tool | Description |
|---|---|
| [hicQuickQC](hicQuickQC.md) | Estimates the quality of a Hi-C dataset from a sample of read pairs. |
| [hicQC](hicQC.md) | Plots QC measures from the log files produced by hicBuildMatrix. |
| [hicPrepareQCreport](hicPrepareQCreport.md) | Tabulates and plots QC measures from multiple hicBuildMatrix log files into one HTML report. |
| [hicCorrelate](hicCorrelate.md) | Computes and visualizes the correlation of Hi-C matrices. |
| [hicPlotDistVsCounts](hicPlotDistVsCounts.md) | Plots the decay in interaction frequency with genomic distance. |
| [hicInfo](hicInfo.md) | Prints information about one or more Hi-C matrices (bins, bin size, sum, min, max, and so on). |

## Analysis

| Tool | Description |
|---|---|
| [hicCompareMatrices](hicCompareMatrices.md) | Computes the difference, ratio or log2 ratio between two matrices. |
| [hicPCA](hicPCA.md) | Computes the eigenvectors for A/B compartments. |
| [hicTransform](hicTransform.md) | Computes an observed/expected matrix (Lieberman-Aiden 2009), a Pearson correlation matrix and/or a covariance matrix. |
| [hicAverageRegions](hicAverageRegions.md) | Averages a set of given locations across one or more matrices, usually TAD regions. |
| [hicDetectLoops](hicDetectLoops.md) | Detects enriched interaction regions (peaks/loops). |
| [hicDetectStripes](hicDetectStripes.md) | Detects architectural stripes, a C++ port of Stripenn's detection method. New in v4. |
| [hicValidateLocations](hicValidateLocations.md) | Compares called loops with known protein peak positions. |
| [hicMergeLoops](hicMergeLoops.md) | Merges loop calls made at different resolutions. |
| [hicHyperoptDetectLoops](hicHyperoptDetectLoops.md) | Searches for the best hicDetectLoops parameter setting for a given dataset. |
| [hicHyperoptDetectLoopsHiCCUPS](hicHyperoptDetectLoopsHiCCUPS.md) | Searches for the best Juicer HiCCUPS parameter setting for a given dataset. |
| [hicCompartmentalization](hicCompartmentalization.md) | Computes the global compartmentalization (polarization) signal from PCA eigenvectors. |
| [hicPlotSVL](hicPlotSVL.md) | Computes and plots short-range versus long-range contacts. |
| [hicDifferentialAnalysis](hicDifferentialAnalysis.md) | Replicate-aware, count-based differential analysis of TADs/boundaries, loops or compartments. New in v4, not in the Python HiCExplorer. |

## TADs

| Tool | Description |
|---|---|
| [hicFindTADs](hicFindTADs.md) | Identifies Topologically Associating Domains (TADs). |
| [hicMergeDomains](hicMergeDomains.md) | Merges TAD domains called at different resolutions and their hierarchical relation. |
| [hicDifferentialTAD](hicDifferentialTAD.md) | Identifies differential TADs between two Hi-C matrices. |
| [hicMergeTADbins](hicMergeTADbins.md) | Uses a BED file of domains or TAD boundaries to merge the bin counts of a Hi-C matrix. |
| [hicInterIntraTAD](hicInterIntraTAD.md) | Computes and plots the inter-TAD versus intra-TAD contact ratio. |
| [hicTADClassifier](hicTADClassifier.md) | Predicts TAD boundaries with a supervised classifier. |
| [hicTrainTADClassifier](hicTrainTADClassifier.md) | Trains and tests the classifiers used by hicTADClassifier. |

## Visualization

| Tool | Description |
|---|---|
| [hicPlotMatrix](hicPlotMatrix.md) | Plots a Hi-C matrix as a heatmap. |
| [hicPlotViewpoint](hicPlotViewpoint.md) | Plots the interactions around a reference point or region. |
| [hicAggregateContacts](hicAggregateContacts.md) | Plots aggregated Hi-C sub-matrices for a list of positions. |
| [hicPlotAverageRegions](hicPlotAverageRegions.md) | Visualizes the output of hicAverageRegions. |

## Matrix handling

| Tool | Description |
|---|---|
| [hicConvertFormat](hicConvertFormat.md) | Converts between Hi-C interaction matrix formats (hic, cool, mcool, h5, homer, HiC-Pro, ginteractions). |
| [hicAdjustMatrix](hicAdjustMatrix.md) | Keeps, removes or masks specified regions of a matrix. |
| [hicCreateThresholdFile](hicCreateThresholdFile.md) | Writes a distance-threshold file (one row per bin over a range) for hicDetectLoops. |

## Capture Hi-C

| Tool | Description |
|---|---|
| [chicQualityControl](chicQualityControl.md) | Quality control for capture Hi-C viewpoints: checks sparsity and removes viewpoints that are too sparse. |
| [chicViewpointBackgroundModel](chicViewpointBackgroundModel.md) | Computes the background model used by the other capture Hi-C tools. |
| [chicViewpoint](chicViewpoint.md) | Computes one viewpoint file per sample per reference point, based on the background model. |
| [chicSignificantInteractions](chicSignificantInteractions.md) | Detects significant interactions per viewpoint from the background model. |
| [chicAggregateStatistic](chicAggregateStatistic.md) | Aggregates viewpoint data of two samples into targets for the differential test. |
| [chicDifferentialTest](chicDifferentialTest.md) | Tests for differential interactions between two samples (chi2 or Fisher). |
| [chicPlotViewpoint](chicPlotViewpoint.md) | Plots a viewpoint together with its background model, significant and differential regions. |
| [chicExportData](chicExportData.md) | Exports the data stored in the capture Hi-C intermediate HDF5 files to text or bigWig. |
| [chicChicagoBackgroundModel](chicChicagoBackgroundModel.md) | Estimates the CHiCAGO background model from `.chinput` files or matrices. |
| [chicChicagoScores](chicChicagoScores.md) | Computes CHiCAGO p-values and weighted scores per interaction. |
| [chicChicagoSignificantInteractions](chicChicagoSignificantInteractions.md) | Selects the interactions above the CHiCAGO score threshold. |
| [chicChicagoPlotViewpoint](chicChicagoPlotViewpoint.md) | Plots CHiCAGO viewpoints as scatter or arcs, or writes a pyGenomeTracks links file. |
