# hicPlotViewpoint

Plots the interactions around a reference point or region.

```text
usage: hicPlotViewpoint --matrix MATRIX [MATRIX ...] --region REGION
                        --outFileName OUTFILENAME --referencePoint
                        REFERENCEPOINT [--chromosome CHROMOSOME]
                        [--interactionOutFileName INTERACTIONOUTFILENAME]
                        [--dpi DPI] [--version] [--help]
```

Plots the number of interactions around a given reference point in a region.

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrix MATRIX [MATRIX ...], -m MATRIX [MATRIX ...]` | Hi-C matrix to plot. |
| `--region REGION` | The format is chr:start-end. |
| `--outFileName OUTFILENAME, -o OUTFILENAME` | File name of the image to save. |
| `--referencePoint REFERENCEPOINT, -rp REFERENCEPOINT` | Reference point. Needs to be in the format: 'chr:100' for a single reference point or 'chr:100-200' for a reference region. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--chromosome CHROMOSOME, -C CHROMOSOME` | Optional parameter: Only show results for this chromosome. |
| `--interactionOutFileName INTERACTIONOUTFILENAME, -i INTERACTIONOUTFILENAME` | Optional parameter: If set, a bedgraph file with all interaction will be created. |
| `--dpi DPI` | Optional parameter: Resolution for the image in case theouput is a raster graphics image (e.g png, jpg) (Default: 300). |
| `--version` | show program's version number and exit |
| `--help, -h` | show this help message and exit |

## Additional notes

C++ port: the interactions are summed in C++ and the figure is drawn by the hicexplorer_plot drawing layer with the matplotlib calls of the Python tool (HICX_PLOT_PYTHON names the interpreter). The C++-only option --plotData FILE writes the data of the figure as JSON to FILE instead of drawing it.
