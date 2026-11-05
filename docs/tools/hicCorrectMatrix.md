# hicCorrectMatrix

Uses iterative correction (ICE) or Knight-Ruiz to remove biases from a Hi-C matrix.

```text
usage: hicCorrectMatrix [-h] [--version]  ...
```

This tool provides two balancing methods that can be applied to a raw matrix:

1. **KR**: balances a matrix using the fast balancing algorithm introduced by Knight and Ruiz (2012).
2. **ICE**: iterative correction of a Hi-C matrix (see Imakaev et al. 2012, Nature Methods, for details).

`hicCorrectMatrix` has two subcommands, `correct` and `diagnostic_plot`.

## correct

Runs Knight-Ruiz matrix balancing (KR) or iterative matrix correction (ICE).

```text
usage: hicCorrectMatrix correct --matrix MATRIX --outFileName OUTFILENAME
                                [--correctionMethod STR]
                                [--filterThreshold FILTERTHRESHOLD FILTERTHRESHOLD]
                                [--iterNum INT] [--inflationCutoff INFLATIONCUTOFF]
                                [--transCutoff TRANSCUTOFF]
                                [--sequencedCountCutoff SEQUENCEDCOUNTCUTOFF]
                                [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]
                                [--skipDiagonal] [--perchr] [--filteredBed FILTEREDBED]
                                [--verbose] [--compatMode {v3,v4}] [--threads INT]
                                [--help]
```

### Required arguments

| Flag | Meaning |
|---|---|
| `--matrix MATRIX, -m MATRIX` | Name of the Hi-C matrix to correct in .h5 format. |
| `--outFileName OUTFILENAME, -o OUTFILENAME` | File name to save the resulting matrix. The output is a `.h5` file. |

### Optional arguments

| Flag | Meaning |
|---|---|
| `--correctionMethod STR` | Method used for matrix correction: `KR` or `ICE` (default: KR). |
| `--filterThreshold FILTERTHRESHOLD FILTERTHRESHOLD, -t ...` | Removes bins of low or large coverage. Applied only for ICE. |
| `--iterNum INT, -n INT` | Number of iterations to compute. Only for ICE (default: 500). |
| `--inflationCutoff INFLATIONCUTOFF` | Maximum number of times a bin can be scaled up during the iterative correction. Only for ICE. |
| `--transCutoff TRANSCUTOFF, -transcut TRANSCUTOFF` | Clip high counts in the top `-transcut` trans regions. Only for ICE. |
| `--sequencedCountCutoff SEQUENCEDCOUNTCUTOFF` | Discard bins covered by fewer reads than this fraction. Only for ICE. |
| `--chromosomes CHROMOSOMES [CHROMOSOMES ...]` | List of chromosomes to include in the correction. |
| `--skipDiagonal, -s` | If set, diagonal counts are not included. Only for ICE. |
| `--perchr` | Normalize each chromosome separately. |
| `--filteredBed FILTEREDBED` | Print bins filtered out by `--filterThreshold` to this file. |
| `--verbose` | Print processing status. |
| `--compatMode {v3,v4}` | `v4` (default) balances in float64. `v3` reproduces krbalancing's float32 input rounding and float32 rescaling accumulators. Only for KR; not a Python option. |
| `--threads INT` | Worker threads. The result does not depend on this (default: 4). Not a Python option. |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## diagnostic_plot

Plots a histogram of the coverage per bin together with the modified z-score, to help choose
`--filterThreshold` for ICE.

```text
usage: hicCorrectMatrix diagnostic_plot --matrix hic_matrix.h5 -o file.png
```

### Required arguments

| Flag | Meaning |
|---|---|
| `--matrix MATRIX, -m MATRIX` | Name of the Hi-C matrix to correct in .h5 format. |
| `--plotName PLOTNAME, -o PLOTNAME` | File name to save the diagnostic plot. |

### Optional arguments

| Flag | Meaning |
|---|---|
| `--chromosomes CHROMOSOMES [CHROMOSOMES ...]` | List of chromosomes to include in the iterative correction. The order given is kept for the resulting corrected matrix. |
| `--xMax XMAX` | Max value for the x-axis in counts per bin. |
| `--perchr` | Compute the histogram per chromosome. For samples from cells with an uneven number of chromosomes and/or translocations, it is advisable to check the histograms per chromosome to find the most conservative `filterThreshold`. |
| `--verbose` | Print processing status. |
| `--help, -h` | show this help message and exit |

C++ port: the coverage per bin and its MAD are computed in C++, and the histogram is drawn by the
`hicexplorer_plot` drawing layer with the matplotlib calls of the Python tool (`HICX_PLOT_PYTHON` names
the interpreter). The C++-only option `--plotData FILE` writes the data of the figure as JSON to FILE
instead of drawing it.

!!! note "diagnostic_plot is implemented, not stubbed out"
    The top-level `hicCorrectMatrix --help` synopsis line for this subcommand still reads "Not
    implemented in the C++ port, see cpp/PLAN.md tier 7"; that description text is stale relative to the
    built binary. `diagnostic_plot`'s own `--help` (reproduced above) shows it is implemented end to end,
    including the C++-only `--plotData` output, delegating only the actual figure drawing to the Python
    plotting layer like the other plotting tools. This page describes the real, current binary.

## Notes

### Knight-Ruiz correction

Alongside the classic iterative correction, HiCExplorer offers the Knight-Ruiz balancing algorithm:

```bash
hicCorrectMatrix correct --matrix matrix.cool --correctionMethod KR \
  --chromosomes chrUextra chr3LHet --outFileName corrected_KR.cool
```

Iterative correction is used via:

```bash
hicCorrectMatrix correct --matrix matrix.cool --correctionMethod ICE \
  --chromosomes chrUextra chr3LHet --iterNum 500 \
  --outFileName corrected_ICE.cool --filterThreshold -1.5 5.0
```

!!! note "Runtime and memory figures below are from HiCExplorer 3"
    The following per-resolution runtime and peak memory figures were measured with HiCExplorer 3 on
    Rao 2014 GM12878 primary + replicate data and do not describe version 4. See
    [Benchmarks](../benchmarks.md) for measured numbers of version 4.

    - KR on 25kb: 165 GB, 1:08 h
    - ICE on 25kb: 224 GB, 3:10 h
    - KR on 10kb: 228 GB, 1:42 h
    - ICE on 10kb: 323 GB, 4:51 h
    - KR on 1kb: 454 GB, 16:50 h
    - ICE on 1kb: over 600 GB, over 2.5 days (the Python authors interrupted the computation and
      strongly recommended using KR at this resolution instead)
