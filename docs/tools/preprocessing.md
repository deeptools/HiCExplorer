# Pre-processing

| Tool | Description |
|---|---|
| [hicFindRestSite](hicFindRestSite.md) | Identifies the genomic locations of restriction sites. |
| [hicBuildMatrix](hicBuildMatrix.md) | Creates a Hi-C matrix from the aligned BAM/SAM files (or a .pairs file) of the Hi-C sequencing reads. |
| [hicBuildMatrixMicroC](hicBuildMatrixMicroC.md) | Creates a Hi-C matrix from Micro-C data (hicBuildMatrix without the restriction-site machinery). |
| [hicSumMatrices](hicSumMatrices.md) | Adds Hi-C matrices of the same size. |
| [hicMergeMatrixBins](hicMergeMatrixBins.md) | Merges consecutive bins on a Hi-C matrix to reduce resolution. |
| [hicCorrectMatrix](hicCorrectMatrix.md) | Uses iterative correction (ICE) or Knight-Ruiz to remove biases from a Hi-C matrix. |
| [hicNormalize](hicNormalize.md) | Normalizes matrices to a 0-1 range or to the smallest total read count. |
