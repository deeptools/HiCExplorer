# hicPlotMatrix

Plots a Hi-C matrix as a heatmap.

```text
usage: hicPlotMatrix --matrix MATRIX --outFileName OUTFILENAME
                     [--title TITLE] [--scoreName SCORENAME]
                     [--perChromosome] [--clearMaskedBins]
                     [--chromosomeOrder CHROMOSOMEORDER [CHROMOSOMEORDER ...]]
                     [--region REGION] [--region2 REGION2] [--log1p] [--log]
                     [--colorMap COLORMAP] [--vMin VMIN] [--vMax VMAX]
                     [--dpi DPI] [--bigwig BIGWIG [BIGWIG ...]]
                     [--bigwigAdditionalVerticalAxis]
                     [--vMinBigwig VMINBIGWIG] [--vMaxBigwig VMAXBIGWIG]
                     [--flipBigwigSign]
                     [--scaleFactorBigwig SCALEFACTORBIGWIG]
                     [--fontsize FONTSIZE] [--rotationX ROTATIONX]
                     [--rotationY ROTATIONY]
                     [--increaseFigureWidth INCREASEFIGUREWIDTH]
                     [--increaseFigureHeight INCREASEFIGUREHEIGHT]
                     [--loops LOOPS]
                     [--loopLargeRegionsOperation {first,last,center}]
                     [--tads TADS] [--help] [--version]
```

Creates a heatmap of a Hi-C matrix.

The C++ port computes the matrices, positions and extents of the figure, then draws it through the
`hicexplorer_plot` drawing layer with the same matplotlib calls as the Python tool
(`HICX_PLOT_PYTHON` names the interpreter to use for drawing); see the Python tool's own documentation
for every plotting/styling option (`--title`, `--colorMap`, `--vMin`/`--vMax`, `--log1p`/`--log`,
`--region`/`--region2`, `--perChromosome`, the bigwig overlay options, and so on). The C++-only option
`--plotData FILE` writes the data of the figure as JSON to FILE (and its matrices as `FILE.<n>.npy`)
instead of drawing it.

## Notes

### Details

`hicPlotMatrix` takes a Hi-C matrix and plots the interactions of all or some chromosomes.

### Examples

Hi-C data from wild-type *D. melanogaster* embryos:

![Hi-C matrix heatmap](../images/hicPlotMatrix.png)

This plot shows all contacts of a Hi-C matrix. Its bins were merged into 25 kb bins using
[hicMergeMatrixBins](hicMergeMatrixBins.md). Alternatively, chromosomes can be plotted separately:

![Hi-C matrix heatmap, per chromosome](../images/hicPlotMatrix_perChromosome.png)

```bash
hicPlotMatrix -m Dmel.h5 -o hicPlotMatrix.png \
  -t 'D.melanogaster (--perChromosome)' --log1p \
  --clearMaskedBins --chromosomeOrder 2L 2R 3L 3R X --perChromosome
```
