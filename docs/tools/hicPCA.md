# hicPCA

Computes the eigenvectors for A/B compartments.

```text
usage: hicPCA --matrix MATRIX --outputFileName OUTPUTFILENAME
              [OUTPUTFILENAME ...]
              [--whichEigenvectors WHICHEIGENVECTORS [WHICHEIGENVECTORS ...]]
              [--format {bedgraph,bigwig}]
              [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]
              [--method {dist_norm,lieberman}] [--ligation_factor]
              [--extraTrack EXTRATRACK] [--histonMarkType HISTONMARKTYPE]
              [--pearsonMatrix PEARSONMATRIX] [--obsexpMatrix OBSEXPMATRIX]
              [--ignoreMaskedBins] [--compatMode {v3,v4}]
              [--eigenSolver {dense,lanczos}] [--threads THREADS]
              [--help] [--version]
```

Computes PCA eigenvectors for a Hi-C matrix.

```bash
$ hicPCA --matrix hic_matrix.h5 -o pca1.bedgraph pca2.bedgraph
```

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrix MATRIX, -m MATRIX` | HiCExplorer matrix in h5 format. |
| `--outputFileName OUTPUTFILENAME [OUTPUTFILENAME ...], -o ...` | File names for the result of the pca. Number of output files must match the number of computed eigenvectors. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--whichEigenvectors WHICHEIGENVECTORS [...], -we ...` | The list of eigenvectors that the PCA should compute e.g. 1 2 5 will return the first, second and fifth eigenvector. (Default: 1 2). |
| `--format {bedgraph,bigwig}, -f {bedgraph,bigwig}` | Output format. (Default: bigwig). |
| `--chromosomes CHROMOSOMES [CHROMOSOMES ...]` | List of chromosomes to be included in the correlation. |
| `--method {dist_norm,lieberman}` | Method used to build the obs-exp matrix. (Default: dist_norm). |
| `--ligation_factor` | Multiply a scaling factor to each entry of the expected matrix, as the Homer software does. Only effective with the dist_norm method. |
| `--extraTrack EXTRATRACK` | A gene track (bed) or a histone mark coverage file (bigwig) used to decide the sign of the eigenvector. |
| `--histonMarkType HISTONMARKTYPE` | active or inactive. (Default: active). |
| `--pearsonMatrix PEARSONMATRIX, -pm PEARSONMATRIX` | Write the intermediate Pearson matrix to this file. |
| `--obsexpMatrix OBSEXPMATRIX, -oem OBSEXPMATRIX` | Write the intermediate observed/expected matrix here. |
| `--ignoreMaskedBins` | Remove the masked bins before the PCA is computed. |
| `--compatMode {v3,v4}` | v3 reproduces scipy.linalg.eig, the general solver, in LAPACK's own column order and sign, which is what the Python does. v4 uses the symmetric solver dsyevr for the requested eigenvectors only, sorted by descending eigenvalue with a deterministic sign. Not a Python option; see cpp/PLAN.md 5.8. (Default: v3). |
| `--eigenSolver {dense,lanczos}` | dense forms each chromosome's covariance and solves it as --compatMode says. lanczos computes the requested eigenvectors of the same covariance from the sparse obs/exp matrix without forming it, ordered by eigenvalue, which is not always the Python's order. Not a Python option; see cpp/PLAN.md tier 11. (Default: dense). |
| `--threads THREADS` | Workers for chromosomes, covariance rows and sparse products. The output is byte-identical for any value. (Default: 4). |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |
