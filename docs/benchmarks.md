# Benchmarks

On real data, HiCExplorer 4.0 runs the analysis tools of HiCExplorer 3.7.6 three to eight times
faster with two to four times less peak memory, and `hicBuildMatrix` needs about a twelfth of the
CPU time. Wall-clock time, CPU time and peak resident memory were recorded with `/usr/bin/time -v`;
memory is given in decimal units (1 GB = 10<sup>9</sup> bytes). The script that produced the
numbers is `cpp/scripts/bench/remeasure.sh`, and every run kept its log. Both versions ran one
after the other on the same workstation (32 cores, 125 GB memory), each run starting at a load
average below 4, with the same options.

## HiCExplorer 4.0 against 3.7.6

| Task | 3.7.6 | 4.0 | Speedup | Peak memory 3.7.6 | Peak memory 4.0 |
|---|---|---|---|---|---|
| `hicDetectLoops` | 20.0 s | 6.3 s | 3.2x | 636 MB | 263 MB |
| `hicFindTADs` | 32.3 s | 4.0 s | 8.2x | 6.25 GB | 1.52 GB |
| `hicPCA` | 11.1 s | 4.0 s | 2.8x | 554 MB | 166 MB |

`hicDetectLoops` and `hicFindTADs` ran on chromosome 1 of GM12878 at 10 kb (24,926 bins, 61.8
million non-zero pixels, corrected with `hicCorrectMatrix`), `hicPCA` on the same chromosome merged
to 100 kb, the resolution at which compartments are analysed. Each tool ran three times; the table
gives the median. `hicDetectLoops` ran with `--threads 1` and its default of four threads per
chromosome, `hicFindTADs` with `--numberOfProcessors 1`. `hicPCA` gives the same eigenvector in
both versions. `hicFindTADs` calls 744 domains (3.7.6: 743) and `hicDetectLoops` 1,221 loops
(3.7.6: 859); both differences come from corrected computations of 3.7.6, the expected values of
the observed/expected matrix for the loops and the boundary search for the domains, which the
[news](news.md) describe.

### hicBuildMatrix

| Threads | Wall clock 3.7.6 | Wall clock 4.0 | CPU time 3.7.6 | CPU time 4.0 | Peak memory 3.7.6 | Peak memory 4.0 |
|---|---|---|---|---|---|---|
| 1 | 65 min 42 s | 46 min 57 s | 2,245 s | 197 s | 15.3 GB | 6.3 GB |
| 4 | 25 min 19 s | 14 min 53 s | 2,264 s | 177 s | 16.0 GB | 6.3 GB |
| 8 | 25 min 23 s | 11 min 20 s | 2,536 s | 177 s | 17.1 GB | 6.3 GB |

`hicBuildMatrix` built a 10 kb matrix from 92.7 million sequenced read pairs (SRR927086, DpnII),
given as two SAM files of 58 GB in total on network storage. Version 4 needs 11 to 14 times less CPU
time. Its wall clock is set by reading the input: the process used 7 to 25 percent of one core
while it waited for data, so the wall clock reflects the throughput of the storage more than the
tool. The memory of version 4 does not grow with the number of threads.

Both versions report 73,506,771 usable pairs. Version 4 counts 55,243,683 Hi-C contacts against
57,566,986 in 3.7.6. The difference is the 2,323,303 self circles, which 3.7.6 counted but kept in
the matrix; version 4 removes them unless `--keepSelfCircles` is set. Version 4 classifies a pair as
"one mate not unique" only when both mates have mapping quality 0 (866,537 pairs; 3.7.6 tested the
first mate only and counted 3,941,063).

## The new methods

| Task | Software | Time | Peak memory |
|---|---|---|---|
| CHiCAGO, 22,076 baits, from `.chinput` | HiCExplorer 4 | 278 s | 21.2 GB |
| CHiCAGO, 22,076 baits, from a cool matrix | HiCExplorer 4 | 47.5 s | 3.5 GB |
| Stripes, 22 GM12878 autosomes at 10 kb | Stripenn 1.1.65.22 | 3,015 s | 7.5 GB |
| | HiCExplorer 4 (`hicDetectStripes`) | 210 s | 2.3 GB |

The CHiCAGO runs scored a genome-wide capture Hi-C design (837,161 HindIII fragments, 22,076
baits, 119.8 million `.chinput` rows) with 16 threads; the time is the sum of
`chicChicagoBackgroundModel` and `chicChicagoScores`. The matrix row starts from the cool file
that `hicConvertFormat --inputFormat chinput` writes from the same `.chinput` file in 48 s.
Stripes were called with eight threads by both tools; the table gives the first of two runs each
(second runs: 3,053 s and 209 s). Stripenn runs worker processes, and `/usr/bin/time` counts only
the main process, so its memory is a lower bound.

## The same results

Faster tools are only useful if they give the same results. Every tool of version 4 is checked
against HiCExplorer 3.7.6 on real and synthetic data by the equivalence suite of the repository
(`cpp/scripts/equiv.py`). Where 3.7.6 was wrong, version 4 corrects it and the suite checks the
corrected result instead; these differences are listed in the [news](news.md). The new methods are
checked against their reference software:

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
