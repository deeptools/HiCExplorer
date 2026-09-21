# hicCompartmentalization

Computes the global compartmentalization (polarization) signal from PCA eigenvectors.

```text
usage: hicCompartmentalization --obsexp_matrices OBSEXP_MATRICES
                               [OBSEXP_MATRICES ...] --pca PCA
                               --outputFileName OUTPUTFILENAME
                               [--quantile QUANTILE] [--outliers OUTLIERS]
                               [--outputMatrix OUTPUTMATRIX]
                               [--offset OFFSET [OFFSET ...]] [--noPlot] [-h]
                               [--version]
```

!!! note "Required arguments must be present for `--help`/`-h` to print"
    Like several other tools in this rewrite, `hicCompartmentalization` validates its required
    arguments before honoring `-h`; a bare `hicCompartmentalization -h` reports the missing required
    arguments instead of the usage text. Pass `-h` together with dummy values for the required options
    to see the full help text, or refer to the tables below.

Rearranges the average interaction frequencies using the first PC values to represent the global
compartmentalization signal. To the original authors' knowledge, this was first introduced and
implemented by Wibke Schwarzer et al. 2017 (Nature. 2017 Nov 2; 551(7678): 51-56).

```bash
hicCompartmentalization --obsexp_matrices obsExpMatrix.h5 --pca pc1.bedgraph -o \
  global_signal.png --noPlot
```

## Required arguments

| Flag | Meaning |
|---|---|
| `--obsexp_matrices OBSEXP_MATRICES [...], -m OBSEXP_MATRICES ...` | HiCExplorer matrices in h5/cool format. |
| `--pca PCA` | A PCA vector as a bedgraph file with no header. |
| `--outputFileName OUTPUTFILENAME, -o OUTPUTFILENAME` | Plot to represent the polarization of A/B compartments. The ratios are written to `OUTPUTFILENAME_dat`. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--quantile QUANTILE, -q QUANTILE` | Number of quantiles (default: 30). |
| `--outliers OUTLIERS` | Percentage of outliers to remove (default: 0). |
| `--outputMatrix OUTPUTMATRIX` | Output `.npz` file including all the generated matrices. |
| `--offset OFFSET [OFFSET ...]` | Set NaN for the distances mentioned as offset from the main diagonal; only positive values are accepted. |
| `--noPlot` | C++ port only: write the numeric outputs (`OUTPUTFILENAME_dat` and `--outputMatrix`) without the plot. |
| `--plotData FILE` | C++ port only: write the data of the plot as JSON to FILE instead of drawing it. |
| `-h` | show the help message and exit. |
| `--version` | show program's version number and exit |

C++ port: the ratios and matrices are computed in C++, and the plot is drawn by the `hicexplorer_plot`
drawing layer with the matplotlib calls of the Python tool (`HICX_PLOT_PYTHON` names the interpreter).

## Notes

### PCA to compute the global compartmentalization signal

A global (genome-wide) strength for compartmentalization is computed as `(AA + BB) / (AB + BA)` after
rearranging the bins of the observed/expected matrix based on their corresponding PC1 values: PC1 values
are first reordered incrementally, then the same order of bins is used to rearrange the bins of the
observed/expected matrix.

```bash
hicCompartmentalization --obsexp_matrices obsExpMatrix.h5 --pca pc1.bedgraph -o global_signal.png
```
