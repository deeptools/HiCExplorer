# chicChicagoScores

Applies a fitted CHiCAGO background model to the interaction counts of one sample and writes CHiCAGO's log p-value and weighted score for every (bait, other end) pair. Second step of the CHiCAGO pipeline, after [chicChicagoBackgroundModel](chicChicagoBackgroundModel.md).

This tool is new in version 4.

```text
usage: chicChicagoScores --backgroundModel BACKGROUNDMODEL
                         (--chinput CHINPUT | --matrices MATRICES [MATRICES ...])
                         --baitmap BAITMAP --rmap RMAP [--outFileName OUTFILENAME]
                         [--minFragLen MINFRAGLEN] [--maxFragLen MAXFRAGLEN]
                         [--minNPerBait MINNPERBAIT] [--maxLBrownEst MAXLBROWNEST]
                         [--noRemoveAdjacent] [--weightAlpha WEIGHTALPHA]
                         [--weightBeta WEIGHTBETA] [--weightGamma WEIGHTGAMMA]
                         [--weightDelta WEIGHTDELTA] [--threads THREADS]
                         [--help] [--version]
```

For each pair the expected count is the Brownian mean (bait factor, other-end factor and the fitted distance function) plus the technical-noise mean of the pair's pool. The p-value is the survival function of a Delaporte distribution (negative binomial plus Poisson) at the observed count. The score is the distance-weighted `-log(p)` of CHiCAGO's `weightedRelative` method; R's default significance threshold is a score of 5.

Pairs whose other end was dropped by the model's trans-count pooling, and baits without a fitted factor, are skipped, as in R. The input has to be read with the same `--minFragLen`, `--maxFragLen`, `--minNPerBait`, `--maxLBrownEst` and `--noRemoveAdjacent` values as the model.

The output is a tab-separated file with the columns `baitID`, `otherEndID`, `N`, `distSign` (`NA` for trans pairs), `Bmean`, `Tmean`, `log_p` and `score`. Its row order follows the input and does not depend on `--threads`.

## Required arguments

| Flag | Meaning |
|---|---|
| `--backgroundModel BACKGROUNDMODEL` | Output file of chicChicagoBackgroundModel. |
| `--chinput CHINPUT` | One `.chinput` interaction count file. Mutually exclusive with `--matrices`; exactly one of the two is required. |
| `--matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]` | One or more Hi-C matrices (cool, h5 or `.hic`) to derive the counts from; cis pairs within `--maxLBrownEst` only. Several matrices are summed per fragment pair. Mutually exclusive with `--chinput`. |
| `--baitmap BAITMAP` | CHiCAGO `.baitmap` file. |
| `--rmap RMAP` | CHiCAGO `.rmap` file. Together with the `.baitmap` it gives the distance-weighting normalization constant. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--outFileName OUTFILENAME, -o OUTFILENAME` | Name of the scores file (Default: chicago_scores.txt). |
| `--minFragLen MINFRAGLEN` | Minimum other-end fragment length (Default: 150). |
| `--maxFragLen MAXFRAGLEN` | Maximum other-end fragment length (Default: 40000). |
| `--minNPerBait MINNPERBAIT` | Minimum total read count per bait (Default: 250). |
| `--maxLBrownEst MAXLBROWNEST` | Maximum distance the Brownian component was estimated over (Default: 1500000). Has to match the model. |
| `--noRemoveAdjacent` | Keep interactions with fragments immediately adjacent to their bait. Has to match the model. |
| `--weightAlpha WEIGHTALPHA`, `--weightBeta WEIGHTBETA`, `--weightGamma WEIGHTGAMMA`, `--weightDelta WEIGHTDELTA` | Parameters of CHiCAGO's distance-weighting curve (Default: R Chicago's own defaults). |
| `--threads THREADS` | Number of threads (Default: 1). Parsing, filtering and scoring run in parallel. |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## Example

```bash
chicChicagoScores --backgroundModel background_model.txt --chinput sample.chinput \
    --baitmap design.baitmap --rmap design.rmap --threads 16 -o scores.txt
```

```text
baitID  otherEndID  N  distSign  Bmean             Tmean               log_p               score
403482  403461      1  -104841   1.49456573187732  0.00313682023102442 -0.367374072616871  0.352650309270243
403482  403462      2  -101616   1.53635641587557  0.00313682023102442 -0.885730815015161  0.872142535750181
```

## Performance

On the real full-genome design of [chicChicagoBackgroundModel](chicChicagoBackgroundModel.md#performance) with `--threads 16`, scoring takes 90.4 s at 20.7 GB peak memory from the `.chinput` file and 19.1 s at 3.9 GB from the equivalent matrix. R Chicago needs 229.3 s for the scoring stage. See [Benchmarks](../benchmarks.md#chicago).
