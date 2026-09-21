# hicPlotAverageRegions

Visualizes the output of hicAverageRegions.

```text
usage: hicPlotAverageRegions --matrix MATRIX --outputFile OUTPUTFILE [--log1p]
                             [--log] [--colorMap COLORMAP] [--vMin VMIN]
                             [--vMax VMAX] [--dpi DPI] [--help] [--version]
```

hicPlotAverage regions plots the data computed by hicAverageRegions. It shows the summed up and averaged regions around all given reference points. This tool is useful to plot differences at certain reference points as for example TAD boundaries between samples.

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrix MATRIX, -m MATRIX` | The averaged regions file computed by hicAverageRegions (npz file). |
| `--outputFile OUTPUTFILE, -o OUTPUTFILE` | The averaged regions plot. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--log1p` | Plot log1p of the matrix values. |
| `--log` | Plot log of the matrix values. |
| `--colorMap COLORMAP` | Color map to use for the heatmap. Available values can be seen here: http://matplotlib.org/examples/color/col ormaps_reference.html (Default: hot_r). |
| `--vMin VMIN` | Minimum score value. |
| `--vMax VMAX` | Maximum score value. |
| `--dpi DPI` | Resolution of image ifouput is a raster graphics image (e.g png, jpg) (Default: 300). |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## Additional notes

C++ port: the figure is drawn by the hicexplorer_plot drawing layer with the calls of the Python tool (HICX_PLOT_PYTHON names the interpreter). The C++-only option --plotData FILE writes the data of the figure as JSON to FILE instead of drawing it.

## Notes

See [hicAverageRegions](hicAverageRegions.md) for an example usage.
