# chicChicagoPlotViewpoint

Plots CHiCAGO viewpoints from the output of [chicChicagoScores](chicChicagoScores.md): one bait, or every bait in a genomic region, as a scatter plot or as arcs. It can also write the interactions as a pyGenomeTracks links file.

This tool exists only in the C++ rewrite. It has no Python HiCExplorer counterpart. [chicPlotViewpoint](chicPlotViewpoint.md) reads the HDF5 files of the chicViewpoint workflow and does not read the output of chicChicagoScores.

```text
usage: chicChicagoPlotViewpoint --scores SCORES (--baitID BAITID |
                                --region CHROM START END)
                                [--style {scatter,arcs}] [--baitmap BAITMAP] [--rmap RMAP]
                                [--backgroundModel BACKGROUNDMODEL]
                                [--range RANGE RANGE] [--plevel1 PLEVEL1]
                                [--plevel2 PLEVEL2] [--keepBait2bait] [--onlySignificant]
                                [--outFileName OUTFILENAME] [--outputFormat OUTPUTFORMAT]
                                [--dpi DPI] [--linksFile FILE] [--plotData FILE] [--help]
                                [--version]
```

## Scope and style

`--baitID` plots one bait. `--region CHROM START END` plots every bait that `--baitmap` places in the region, in absolute genomic coordinates; the other end of an interaction lies at the bait midpoint plus `distSign`. A bait-to-bait pair with both baits inside the region is drawn once.

`--style scatter` reproduces R Chicago's `plotBaits()`: signed distance from the bait on the x axis and the observed count N on the y axis. The score only colors a point, in three tiers: black below `--plevel2` (default 3), blue from `--plevel2`, red from `--plevel1` (default 5). A gray vertical line marks the bait. With `--backgroundModel` the Brownian mean and its dashed upper 95% band, `Bmean + 1.96 * sqrt(Bmean + Bmean^2 / dispersion)`, are drawn as well. The scatter style needs `--baitID`.

`--style arcs` draws one arc per interaction between the bait and the other end, with the height `log1p(N)` and the same three color tiers. It works with `--baitID` and with `--region`.

![Scatter viewpoint of one bait](../images/chicago-viewpoint-scatter.png)

*Scatter viewpoint of bait RBM38 on GM12878 chr20 (PCHiCdata example data), with the Brownian mean and upper band.*

![Arc viewpoint of one bait, significant interactions only](../images/chicago-viewpoint-arcs.png)

*The same bait as arcs with `--onlySignificant`.*

![Arcs of all baits in a region](../images/chicago-region-arcs.png)

*All baits in chr20:0-3,000,000 with `--region 20 0 3000000 --style arcs --onlySignificant`. Grey ticks mark the baits.*

## Links file for pyGenomeTracks

`--linksFile FILE` writes the kept interactions as a pyGenomeTracks `links` file (`chr1 start1 end1 chr2 start2 end2 score`, tab separated). It requires `--baitmap`; `--rmap` supplies the real fragment span of the other end, otherwise a 1 bp placeholder at its position is written. Combined with bigWig, gene and BED tracks in a `tracks.ini`, the links are drawn next to histone marks, genes, promoters and known enhancers by [hicPlotTADs](hicPlotTADs.md). See the [CHiCAGO tutorial](../example-usage/chicago-tutorial.md#combine-the-interactions-with-genome-tracks).

## Required arguments

| Flag | Meaning |
|---|---|
| `--scores SCORES` | Output file of chicChicagoScores. |
| `--baitID BAITID` | The baitID to plot. Mutually exclusive with `--region`; exactly one of the two is required. |
| `--region CHROM START END` | Plot every bait in this region. Requires `--baitmap`. Mutually exclusive with `--baitID`. Not valid with `--style scatter`. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--style {scatter,arcs}` | Plot style (Default: scatter). |
| `--baitmap BAITMAP` | CHiCAGO `.baitmap` file. Required for `--region` and `--linksFile`. With `--baitID` it supplies the plot title and enables the removal of bait-to-bait rows. |
| `--rmap RMAP` | CHiCAGO `.rmap` file, for the other-end fragment spans of `--linksFile`. |
| `--backgroundModel BACKGROUNDMODEL` | Output of chicChicagoBackgroundModel. Enables the Brownian overlay of the scatter style. |
| `--range RANGE RANGE` | Upstream and downstream distance from each bait to plot (Default: 1000000 1000000). |
| `--plevel1 PLEVEL1` | Score threshold of the more significant color tier (Default: 5). |
| `--plevel2 PLEVEL2` | Score threshold of the less significant color tier (Default: 3). |
| `--keepBait2bait` | Keep bait-to-bait rows (Default: removed, as in R). Only has an effect with `--baitmap`. |
| `--onlySignificant` | Drop rows with a score below `--plevel2`, so that only the two significant tiers are drawn. |
| `--outFileName OUTFILENAME, -o OUTFILENAME` | Name of the plot file (Default: chicago_viewpoint.png). |
| `--outputFormat OUTPUTFORMAT, -format OUTPUTFORMAT` | Output format of the plot (Default: png). |
| `--dpi DPI` | Resolution for raster output (Default: 300). |
| `--linksFile FILE` | Write the kept interactions as a pyGenomeTracks links file. |
| `--plotData FILE` | C++-only. Write the data the figure is drawn from as JSON and do not draw it. |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

The figure is drawn with matplotlib. Like the other plotting tools, the tool checks the matplotlib version before it reads any input; `HICX_PLOT_PYTHON` selects the interpreter.

## Example

```bash
# one bait, scatter with the background model
chicChicagoPlotViewpoint --scores scores.txt --baitID 417632 \
    --baitmap design.baitmap --backgroundModel background_model.txt \
    -o viewpoint_scatter.png

# all baits of a region, significant interactions as arcs, plus a links file
chicChicagoPlotViewpoint --scores scores.txt --region 20 0 3000000 \
    --baitmap design.baitmap --rmap design.rmap \
    --style arcs --onlySignificant --linksFile chicago.links \
    -o region_arcs.png
```
