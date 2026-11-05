# chicPlotViewpoint

Plots a viewpoint together with its background model, significant and differential regions.

```text
usage: chicPlotViewpoint --interactionFile INTERACTIONFILE --range RANGE RANGE
                         [--backgroundModelFile BACKGROUNDMODELFILE]
                         [--differentialTestResult DIFFERENTIALTESTRESULT]
                         [--significantInteractions SIGNIFICANTINTERACTIONS]
                         [--plotSignificantInteractions]
                         [--outFileName OUTFILENAME]
                         [--outputFormat OUTPUTFORMAT] [--dpi DPI]
                         [--combinationMode {dual,single,allGenes,oneGene}]
                         [--combinationName COMBINATIONNAME]
                         [--colorMapPvalue COLORMAPPVALUE]
                         [--maxPValue MAXPVALUE] [--minPValue MINPVALUE]
                         [--pValue]
                         [--pValueSignificanceLevels PVALUESIGNIFICANCELEVELS [PVALUESIGNIFICANCELEVELS ...]]
                         [--xFold XFOLD] [--truncateZeroPvalues]
                         [--colorList COLORLIST [COLORLIST ...]]
                         [--threads THREADS] [--help] [--version]
```

chicPlotViewpoint plots one or many viewpoints with the average background model and the computed p-value per sample. In addition, it can highlight differential interactions of two samples and/or significant regions.

## Required arguments

| Flag | Meaning |
|---|---|
| `--interactionFile INTERACTIONFILE, -if INTERACTIONFILE` | path to the interaction files which should be used for plotting |
| `--range RANGE RANGE` | Defines the region upstream and downstream of a reference point which should be included. Format is --region upstream downstream, e.g.: --region 500000 500000 plots 500kb up- and 500kb downstream. This value should not exceed the range used in the other chic-tools. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--backgroundModelFile BACKGROUNDMODELFILE, -bmf BACKGROUNDMODELFILE` | path to the background file which should be used for plotting |
| `--differentialTestResult DIFFERENTIALTESTRESULT, -dif DIFFERENTIALTESTRESULT` | Path to the H0 rejected files to highlight the regions in the plot. |
| `--significantInteractions SIGNIFICANTINTERACTIONS, -si SIGNIFICANTINTERACTIONS` | Path to the files with detected significant interactions to highlight the regions in the plot. |
| `--plotSignificantInteractions, -psi` | Highlights the significant interactions in the plot itself. If not set, only the p-values are updated |
| `--outFileName OUTFILENAME, -o OUTFILENAME` | Output tar.gz of the files (Default: plots.tar.gz). |
| `--outputFormat OUTPUTFORMAT, -format OUTPUTFORMAT` | Output format of the plot (Default: png). |
| `--dpi DPI` | Optional parameter: Resolution for the image, ifoutput is a raster graphics image (e.g png, jpg) (Default: 300). |
| `--combinationMode {dual,single,allGenes,oneGene}, -cm {dual,single,allGenes,oneGene}` | This option defines how the interaction data should be computed and combined: dual: Combines as follows: [[matrix1_gene1, matrix2_gene1], [matrix2_gene1, matrix3_gene1],[matrix1_gene2, matrix2_gene2], ...]single: Combines as follows: [matrix1_gene1, matrix1_gene2, matrix2_gene1, ...], allGenes: Combines as follows: [[matrix1_gene1, matrix2_gene1, matrix2_gene1], [matrix1_gene2, matrix2_gene2, matrix3_gene2], ...]oneGene: Computes all data of one gene, please specify '--'. If a gene is not unique, each viewpoint is treated independently. (Default: dual). |
| `--combinationName COMBINATIONNAME, -cn COMBINATIONNAME` | Gene name or file name for modes 'oneGene' or 'file' of parameter '--combinationMode' (Default: None). |
| `--colorMapPvalue COLORMAPPVALUE` | Color map to use for the p-value. Available values can be seen here: http://matplotlib.org/examples/color/col ormaps_reference.html (Default: RdYlBu). |
| `--maxPValue MAXPVALUE, -map MAXPVALUE` | Maximal value for p-value. Values above this threshold are set to this value (Default: 0.1). |
| `--minPValue MINPVALUE, -mp MINPVALUE` | Minimal value for p-value. Values below this threshold are set to this value (Default: 0.0). |
| `--pValue, -p` | Plot p-values as a colorbar |
| `--pValueSignificanceLevels PVALUESIGNIFICANCELEVELS [PVALUESIGNIFICANCELEVELS ...], -psl PVALUESIGNIFICANCELEVELS [PVALUESIGNIFICANCELEVELS ...]` | Highlight the p-values by the defined significance levels. |
| `--xFold XFOLD, -xf XFOLD` | Plot x-fold region for the mean background. |
| `--truncateZeroPvalues, -tzpv` | Sets all p-values which are equal to zero to one. |
| `--colorList COLORLIST [COLORLIST ...], -cl COLORLIST [COLORLIST ...]` | Colorlist for the viewpoint lines (Default g b c m y k). |
| `--threads THREADS, -t THREADS` | Number of threads (uses the python multiprocessing module) (Default: 4). |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## Additional notes

C++ port: the plot data is read and computed in C++, and the figures and the archive are written by the hicexplorer_plot drawing layer with the calls of the Python tool (HICX_PLOT_PYTHON names the interpreter). The C++-only option --plotData FILE writes the data of the figures as JSON to FILE instead of drawing them.
