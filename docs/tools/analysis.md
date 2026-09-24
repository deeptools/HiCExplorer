# Analysis

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
| [hicHyperoptDetectLoops](hicHyperoptDetectLoops.md) | Searches for the best hicDetectLoops parameter setting for a given dataset. Not yet ported to C++. |
| [hicHyperoptDetectLoopsHiCCUPS](hicHyperoptDetectLoopsHiCCUPS.md) | Searches for the best Juicer HiCCUPS parameter setting for a given dataset. Not yet ported to C++. |
| [hicCompartmentalization](hicCompartmentalization.md) | Computes the global compartmentalization (polarization) signal from PCA eigenvectors. |
| [hicPlotSVL](hicPlotSVL.md) | Computes and plots short-range versus long-range contacts. |
| [hicDifferentialAnalysis](hicDifferentialAnalysis.md) | Replicate-aware, count-based differential analysis of TADs/boundaries, loops or compartments. New in v4, not in the Python HiCExplorer. |
