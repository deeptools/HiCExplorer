# chicChicagoBackgroundModel

Estimates the CHiCAGO background model (Cairns et al. 2016) from capture Hi-C interaction counts. First step of the CHiCAGO pipeline; the output feeds [chicChicagoScores](chicChicagoScores.md).

This tool exists only in the C++ rewrite. It has no Python HiCExplorer counterpart.

```text
usage: chicChicagoBackgroundModel --rmap RMAP --baitmap BAITMAP
                                  (--chinput CHINPUT | --matrices MATRICES [MATRICES ...])
                                  --nperbin NPERBIN --nbaitsperbin NBAITSPERBIN
                                  --proxOE PROXOE [--outFileName OUTFILENAME]
                                  [--minFragLen MINFRAGLEN] [--maxFragLen MAXFRAGLEN]
                                  [--minNPerBait MINNPERBAIT] [--maxLBrownEst MAXLBROWNEST]
                                  [--binsize BINSIZE] [--noRemoveAdjacent]
                                  [--tlbFilterTopPercent TLBFILTERTOPPERCENT]
                                  [--tlbMinProxOEPerBin TLBMINPROXOEPERBIN]
                                  [--tlbMinProxB2BPerBin TLBMINPROXB2BPERBIN]
                                  [--techNoiseMinBaitsPerBin TECHNOISEMINBAITSPERBIN]
                                  [--brownianNoiseSubset BROWNIANNOISESUBSET] [--threads THREADS]
                                  [--help] [--version]
```

The model has three parts: a cubic distance function of the expected count against distance from the bait, bait and other-end scaling factors together with a negative-binomial dispersion (the Brownian component), and a Poisson technical-noise term estimated from trans counts per (bait, other-end) pool. The statistics follow R Chicago's default pipeline (`readSample`, `addTLB`, `normaliseBaits`, `normaliseOtherEnds`, `estimateTechnicalNoise`, `estimateDistFun`, `estimateBrownianComponent`) and were checked against the R package on the PCHiCdata GM12878 and mouse ES data.

The counts come from one of two sources:

- `--chinput`: a CHiCAGO `.chinput` file with one row per (bait, other end) pair, cis and trans. Several replicates have to be summed per pair beforehand; the tool reads one file.
- `--matrices`: one or more Hi-C matrices in cool, h5 or `.hic` format, one bin per restriction fragment. The (bait, other end) counts are derived from the matrix, one chromosome at a time. Only cis pairs within `--maxLBrownEst` are derived, so trans-based technical-noise estimation has no input on this path and a `.chinput` file is needed when trans counts matter. A `.chinput` file can be written as a matrix with [hicConvertFormat](hicConvertFormat.md) `--inputFormat chinput`.

The design tables `--nperbin`, `--nbaitsperbin` and `--proxOE` are the `.npb`, `.nbpb` and `.poe` files that `makeDesignFiles.py` from chicagoTools produces from the `.rmap` and `.baitmap`. See [File formats](../file-formats.md#chicago-formats).

`--brownianNoiseSubset` is reported only. R subsamples baits above that count and averages several stochastic fits; this tool always fits the dispersion on the full design and writes `subsetWouldTriggerInR` to the model file when R would have subsampled.

## Required arguments

| Flag | Meaning |
|---|---|
| `--rmap RMAP` | CHiCAGO `.rmap` restriction fragment file. |
| `--baitmap BAITMAP` | CHiCAGO `.baitmap` baited fragment file. |
| `--chinput CHINPUT` | One `.chinput` interaction count file. Mutually exclusive with `--matrices`; exactly one of the two is required. |
| `--matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]` | One or more Hi-C matrices (cool, h5 or `.hic`) to derive the per (bait, other end) counts from. Several matrices are summed per fragment pair. Mutually exclusive with `--chinput`. |
| `--nperbin NPERBIN` | `.npb` NPerBin design table. |
| `--nbaitsperbin NBAITSPERBIN` | `.nbpb` NBaitsPerBin design table. |
| `--proxOE PROXOE` | `.poe` ProxOE design table. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--outFileName OUTFILENAME, -o OUTFILENAME` | Name of the background model file (Default: chicago_background_model.txt). |
| `--minFragLen MINFRAGLEN` | Minimum other-end fragment length (Default: 150). |
| `--maxFragLen MAXFRAGLEN` | Maximum other-end fragment length (Default: 40000). |
| `--minNPerBait MINNPERBAIT` | Minimum total read count per bait (Default: 250). |
| `--maxLBrownEst MAXLBROWNEST` | Maximum distance the Brownian component is estimated over (Default: 1500000). |
| `--binsize BINSIZE` | Distance bin width (Default: 20000). |
| `--noRemoveAdjacent` | Keep interactions with fragments immediately adjacent to their bait (Default: they are removed). |
| `--tlbFilterTopPercent TLBFILTERTOPPERCENT` | Percent of other ends with the highest trans-counts dropped before pooling (Default: 0.01). |
| `--tlbMinProxOEPerBin TLBMINPROXOEPERBIN` | Minimum non-bait2bait other ends per trans-count pool (Default: 50000). |
| `--tlbMinProxB2BPerBin TLBMINPROXB2BPERBIN` | Minimum bait2bait other ends per trans-count pool (Default: 2500). |
| `--techNoiseMinBaitsPerBin TECHNOISEMINBAITSPERBIN` | Minimum baits per technical-noise pool (Default: 1000). |
| `--brownianNoiseSubset BROWNIANNOISESUBSET` | Reported only: the bait count above which R subsamples (Default: 1000). |
| `--threads THREADS` | Number of threads (Default: 1). Parsing, filtering, the aggregations and the dispersion sums run in parallel. The model does not depend on the thread count beyond floating-point rounding of the dispersion sums. |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## Example

```bash
chicChicagoBackgroundModel --rmap design.rmap --baitmap design.baitmap \
    --chinput sample.chinput \
    --nperbin design.npb --nbaitsperbin design.nbpb --proxOE design.poe \
    --threads 16 -o background_model.txt
```

The same model from a matrix instead of a `.chinput` file:

```bash
chicChicagoBackgroundModel --rmap design.rmap --baitmap design.baitmap \
    --matrices sample.cool \
    --nperbin design.npb --nbaitsperbin design.nbpb --proxOE design.poe \
    --threads 16 -o background_model.txt
```

## Performance

On a real full-genome design (837,161 HindIII fragments, 22,076 baits, 119.8M `.chinput` rows) with `--threads 16`, the model takes 120.7 s at 20.7 GB peak memory from the `.chinput` file and 29.4 s at 3.9 GB from the equivalent matrix. R Chicago needs 160.1 s for the model stage (40.4 GB peak for the whole R run). See [Benchmarks](../benchmarks.md#chicago).
