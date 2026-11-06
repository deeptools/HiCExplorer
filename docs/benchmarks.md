# Benchmarks

HiCExplorer 4.0 computes the results of HiCExplorer 3.7.6 two to ten times faster and with two to
four times less memory. The new methods of version 4 also run faster than the reference software they
reproduce. All measurements below are on real data. Wall-clock time and peak resident memory were
recorded with `/usr/bin/time -v`.

## HiCExplorer 4.0 against 3.7.6

| Task | 3.7.6 | 4.0 | Speedup | Peak memory 3.7.6 | Peak memory 4.0 |
|---|---|---|---|---|---|
| `hicBuildMatrix`, 1 thread | 42 min 18 s | 16 min 22 s | 2.6x | 14.9 GB | 6.1 GB |
| `hicBuildMatrix`, 8 threads | 21 min 03 s | 9 min 22 s | 2.2x | 16.4 GB | 6.1 GB |
| `hicDetectLoops` | 20.6 s | 2.1 s | 10.0x | 624 MB | 230 MB |
| `hicFindTADs` | 33.0 s | 6.3 s | 5.2x | 6.0 GB | 1.5 GB |
| `hicPCA` | 24.2 s | 3.8 s | 6.4x | 493 MB | 147 MB |

`hicBuildMatrix` built a 10 kb matrix from 92.7 million sequenced reads (DpnII). Both versions report
the same 73,506,771 usable pairs and 57,566,986 Hi-C contacts. `hicDetectLoops` and `hicFindTADs` ran
on chromosome 1 of GM12878 at 10 kb (24,926 bins, 61.8 million non-zero pixels), single-threaded;
`hicFindTADs` called the same 644 domains in both versions. `hicPCA` ran on the same chromosome at
100 kb, the resolution at which compartments are analysed, with the default dense eigensolver. For
large matrices, `hicPCA --eigenSolver lanczos` is faster still.

### hicBuildMatrix and threads

| Threads | 3.7.6 | 4.0 | Speedup | Peak memory 3.7.6 | Peak memory 4.0 |
|---|---|---|---|---|---|
| 1 | 42 min 18 s | 16 min 22 s | 2.6x | 14.9 GB | 6.1 GB |
| 4 | 23 min 14 s | 10 min 38 s | 2.2x | 15.7 GB | 6.1 GB |
| 8 | 21 min 03 s | 9 min 22 s | 2.2x | 16.4 GB | 6.1 GB |

The memory of version 4 does not grow with the number of threads.

## The new methods against their reference software

| Task | Software | Time | Peak memory |
|---|---|---|---|
| CHiCAGO, 22,076 baits | R CHiCAGO 1.38.0 | 437.5 s | 40.4 GB |
| | HiCExplorer 4, from `.chinput` | 236.6 s | 20.7 GB |
| | HiCExplorer 4, from a cool matrix | 48.5 s | 3.9 GB |
| Stripes, 22 GM12878 autosomes at 10 kb | Stripenn 1.1.65.22 | 3,015 s | 7.2 GB |
| | HiCExplorer 4 (`hicDetectStripes`) | 210 s | 2.2 GB |
| `.pairs` to matrix, 136 million pairs | cooler | 148.0 s | 3.7 GB |
| | HiCExplorer 4 (`hicBuildMatrix`) | 29.5 s | 2.6 GB |
| `.hic` to cool, GM12878 at 25 kb | hic2cool 1.0.1 | 1,951 s | 5.3 GB |
| | hictk 2.2.0, 32 threads | 166 s | 0.8 GB |
| | HiCExplorer 4, 1 thread | 162 s | 0.5 GB |
| | HiCExplorer 4, 32 threads | 34 s | 1.4 GB |
| cool to `.hic`, GM12878 at 25 kb | hictk 2.2.0, 32 threads | 997 s | 9.8 GB |
| | HiCExplorer 4, 32 threads | 56 s | 3.2 GB |

The CHiCAGO tools ran with 16 threads on a genome-wide capture Hi-C design (837,161 HindIII
fragments, 22,076 baits, 119.8 million `.chinput` rows); R CHiCAGO runs as a single process. The
matrix row starts from the cool file that `hicConvertFormat --inputFormat chinput` writes from the
same `.chinput` file in 51 s. Stripes were called with eight threads by both tools. The `.pairs` row
gives CPU time. The GM12878 matrix for the format conversion has 899 million non-zero pixels.

## The same results

Faster tools are only useful if they give the same results. Every tool of version 4 is checked
against HiCExplorer 3.7.6 on real and synthetic data by the equivalence suite of the repository
(`cpp/scripts/equiv.py`). The new methods are checked against their reference software:

- The CHiCAGO tools call the same 1,169 significant interactions as R CHiCAGO on the human data of the
  PCHiCdata package.
- `hicDetectStripes` generates the same 38,535 stripe candidates as Stripenn, with identical p-values,
  and calls the same 887 stripes on the 22 GM12878 autosomes.
- The SCC, GenomeDISCO, HiC-Spector, QuASAR-Rep and ENT3C scores of `hicCorrelate` agree with their
  reference implementations within 10<sup>-14</sup>. For ENT3C this holds on every chromosome
  whose submatrices have at least three bins; with two bins every entropy is zero up to rounding
  and neither implementation gives a meaningful score.
- `.hic` files of versions 6 to 8 convert to the same cool files as with hic2cool, and a matrix
  converted to `.hic` and back returns every pixel unchanged.

Tools with a `--threads` option give the same output for every number of threads.
