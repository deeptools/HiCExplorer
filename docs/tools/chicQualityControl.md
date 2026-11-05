# chicQualityControl

Quality control for capture Hi-C viewpoints: checks sparsity and removes viewpoints that are too sparse.

```text
usage: chicQualityControl --matrices MATRICES [MATRICES ...] --referencePoints
                          REFERENCEPOINTS --sparsity SPARSITY
                          [--outFileName OUTFILENAME]
                          [--outFileNameHistogram OUTFILENAMEHISTOGRAM]
                          [--outFileNameSparsity OUTFILENAMESPARSITY]
                          [--threads THREADS] [--fixateRange FIXATERANGE]
                          [--dpi DPI] [--help] [--version]
```

Computes the sparsity of each viewpoint to determine the quality. A viewpoint is considered to be of bad quality if it is too sparse i.e. if there are too many locations with no interactions recorded.

This script creates three output files: a plot with the sparsity distribution per matrix, a plot with the sparsity distribution as histograms and a filtered reference points file.

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]` | The input matrices to apply the QC on. |
| `--referencePoints REFERENCEPOINTS, -rp REFERENCEPOINTS` | Bed file contains all reference points which are checked for a sufficient number of interactions. |
| `--sparsity SPARSITY, -s SPARSITY` | Viewpoints with a sparsity less than the value given are considered of bad quality. If multiple matrices are given, the viewpoint is removed as soon as it is of bad quality in at least one matrix. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--outFileName OUTFILENAME, -o OUTFILENAME` | The output file name of the passed reference points. Used as prefix for the plots as well (Default: new_referencepoints.bed). |
| `--outFileNameHistogram OUTFILENAMEHISTOGRAM, -oh OUTFILENAMEHISTOGRAM` | The output file for the histogram plot (Default: histogram.png). |
| `--outFileNameSparsity OUTFILENAMESPARSITY, -os OUTFILENAMESPARSITY` | The output file for the sparsity distribution plot (Default: sparsity.png). |
| `--threads THREADS, -t THREADS` | Number of threads (Default: 4). |
| `--fixateRange FIXATERANGE, -fs FIXATERANGE` | Fixate score of background model starting at distance x. E.g. all values greater than 500kb are set to the value of the 500kb bin (Default: 500000). |
| `--dpi DPI` | Optional parameter: Resolution for the image if theoutput is a raster graphics image (e.g png, jpg) (Default: 300). |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## An example usage is

$ chicQualityControl -m matrix1.cool matrix2.cool -rp referencePointsFile.bed --range 20000 40000 --sparsity 0.01 -o referencePointFile_QC_passed.bed

## Additional notes

C++ port: the sparsity and the reference point files are computed in C++, and the two figures are drawn by the hicexplorer_plot drawing layer with the matplotlib calls of the Python tool (HICX_PLOT_PYTHON names the interpreter). The C++-only option --plotData FILE writes the data of the figures as JSON to FILE instead of drawing them.
