# hicDifferentialAnalysis

Replicate-aware, count-based differential analysis of TADs/boundaries, loops or compartments.

```text
usage: hicDifferentialAnalysis [-h] [--version] {tads,loops,compartments} ...
```

!!! note "New in this C++ rewrite"
    `hicDifferentialAnalysis` has no Python HiCExplorer counterpart; it is a new tool added directly in
    C++ (`cpp/PLAN.md` section 9.7).

Replicate-aware, count-based differential analysis of Hi-C contact matrices: TADs and TAD boundaries
(`tads`), loops (`loops`) and A/B compartments (`compartments`). Counts are modelled with a negative
binomial GLM, the dispersion is estimated from the replicates with quasi-likelihood empirical Bayes
moderation, offsets take library size and distance decay into account, bins filtered in any sample are
masked in all samples, a minimum fold change is tested (TREAT), and p-values are adjusted with
Benjamini-Hochberg.

All three subcommands share this model and most of its options (`--blocks`, `--exploratory`, `--fdr`,
`--minFoldChange`, `--filterThreshold`, `--chromosomes`, `--threads`, and the calibration options
`--splitReplicates`, `--plantRegions`, `--plantFold`, `--plantSeed`, used to validate the statistical
test itself rather than for routine use).

## tads

Differential TADs and TAD boundaries. Each TAD is tested on its intra-TAD contacts per distance stratum
and in total (Simes' combination); each boundary between two adjacent TADs is tested on the contacts
crossing it, relative to its two flanks. Writes `<prefix>_tads.tsv` and `<prefix>_boundaries.tsv`.

```text
usage: hicDifferentialAnalysis tads --conditionA MATRIX [MATRIX ...]
                                    --conditionB MATRIX [MATRIX ...]
                                    --domains DOMAINS --outFilePrefix PREFIX
                                    [--blocks LABEL [LABEL ...]] [--exploratory]
                                    [--fdr FDR] [--minFoldChange FOLD]
                                    [--filterThreshold LOWER UPPER]
                                    [--chromosomes CHROM [CHROM ...]]
                                    [--boundaryWindow BP] [--threads THREADS]
                                    [--splitReplicates SEED]
                                    [--plantRegions BED] [--plantFold FOLD]
                                    [--plantSeed SEED] [-h]
```

### Required arguments

| Flag | Meaning |
|---|---|
| `--conditionA MATRIX [...], -a MATRIX [...]` | Cool files of condition A, one per replicate. |
| `--conditionB MATRIX [...], -b MATRIX [...]` | Cool files of condition B, one per replicate. All files must share one bin table; raw counts are used, a weight column is ignored. |
| `--domains DOMAINS, -d DOMAINS` | TAD domains, a BED file such as hicFindTADs' `domains.bed`. |
| `--outFilePrefix PREFIX, -o PREFIX` | Prefix of the output files. |

### Optional arguments

| Flag | Meaning |
|---|---|
| `--blocks LABEL [...]` | One block label per sample, condition A first, for an additive block factor (for example a batch or a paired replicate). The block must not be confounded with the condition. |
| `--exploratory` | Allow a condition with a single sample. The dispersion is then estimated from all samples ignoring the conditions, and every output file is labelled exploratory. Without this option such a comparison is refused. |
| `--fdr FDR` | Benjamini-Hochberg FDR at which a unit is called differential (default: 0.05). |
| `--minFoldChange FOLD` | Minimum fold change the test is against (TREAT); 1 tests for any change (default: 1.1). |
| `--filterThreshold LOWER UPPER` | MAD z-score limits on the cis coverage of a bin, per sample and chromosome; a bin outside them in any sample is masked in all (default: -1.5 5.0). |
| `--chromosomes CHROM [...]` | Chromosomes to analyse (default: all). |
| `--threads THREADS, -t THREADS` | Worker threads; the result does not depend on the number (default: 4). |
| `--boundaryWindow BP` | Width of the windows on either side of a boundary, capped by the two TADs (default: 500000). |

### Calibration arguments

| Flag | Meaning |
|---|---|
| `--splitReplicates SEED` | Split every count of each file into two binomial halves: the first halves form condition A and the second halves condition B. `--conditionA` and `--conditionB` must list the same files in the same order. |
| `--plantRegions BED` | Regions whose contacts are thinned to `1 / --plantFold` in one condition: chrom, start, end and an optional fourth column A or B (default: A). |
| `--plantFold FOLD` | Fold difference of the planted regions (default: 2.0). |
| `--plantSeed SEED` | Seed of the planted thinning (default: 0). |
| `-h, --help` | show this help message and exit |

## loops

Differential loops. The loop calls of all samples are united: positions are mapped to bins of the
matrices, and calls within one bin of each other (in both anchors) are merged into one loop at their
rounded mean position. Each loop is tested on the contacts of its peak square (the loop pixel +/-
`--peakWidth` bins) against the contacts of its local background as the offset: the square of +/-
`--windowSize` bins without the peak square, and the ring of one bin around it. A change is a change of
the loop's enrichment over its background. Writes `<prefix>_loops.tsv`.

```text
usage: hicDifferentialAnalysis loops --conditionA MATRIX [MATRIX ...]
                                     --conditionB MATRIX [MATRIX ...]
                                     --loops LOOPS [LOOPS ...]
                                     --outFilePrefix PREFIX
                                     [--blocks LABEL [LABEL ...]] [--exploratory]
                                     [--fdr FDR] [--minFoldChange FOLD]
                                     [--filterThreshold LOWER UPPER]
                                     [--chromosomes CHROM [CHROM ...]]
                                     [--peakWidth BINS] [--windowSize BINS]
                                     [--threads THREADS]
                                     [--splitReplicates SEED]
                                     [--plantRegions BED] [--plantFold FOLD]
                                     [--plantSeed SEED] [-h]
```

### Required arguments

| Flag | Meaning |
|---|---|
| `--conditionA MATRIX [...], -a MATRIX [...]` | Cool files of condition A, one per replicate. |
| `--conditionB MATRIX [...], -b MATRIX [...]` | Cool files of condition B, one per replicate. All files must share one bin table; raw counts are used, a weight column is ignored. |
| `--loops LOOPS [...]` | Loop calls, one file or more (for example one per sample): BEDPE-like, the first six columns `chrom1 start1 end1 chrom2 start2 end2`, as `hicDetectLoops` writes them. Inter-chromosomal calls are ignored. |
| `--outFilePrefix PREFIX, -o PREFIX` | Prefix of the output file. |

### Optional arguments

| Flag | Meaning |
|---|---|
| `--blocks LABEL [...]` | One block label per sample, condition A first (see `tads` above). |
| `--exploratory` | Allow a single-sample condition (see `tads` above). |
| `--fdr FDR` | Benjamini-Hochberg FDR at which a unit is called differential (default: 0.05). |
| `--minFoldChange FOLD` | Minimum fold change the test is against, TREAT (default: 1.1). |
| `--filterThreshold LOWER UPPER` | MAD z-score limits on cis coverage (default: -1.5 5.0). |
| `--chromosomes CHROM [...]` | Chromosomes to analyse (default: all). |
| `--threads THREADS, -t THREADS` | Worker threads (default: 4). |
| `--peakWidth BINS` | Half width of the peak square in bins (default: 1). |
| `--windowSize BINS` | Half width of the background square in bins; at least `--peakWidth + 2` (default: 5). |

### Calibration arguments

Same as `tads` above: `--splitReplicates`, `--plantRegions`, `--plantFold`, `--plantSeed`.

## compartments

Differential compartment preference. A consensus compartment track is computed per chromosome from the
pooled samples: the Pearson correlation of the pooled observed-over-expected matrix, of its three
leading eigenvectors the one that correlates most with GC content, oriented so that A (positive) is GC
rich. Every bin is then tested on its contacts with A bins, with its contacts with B bins as the offset
(both at a distance of at least `--minDistance`): a change is a change of the bin's preference for the A
over the B compartment. The per-sample score is log2 of observed-over-expected contacts with A over the
same with B. Writes `<prefix>_compartments.tsv`. A planted BED region thins the region's contacts with
the A bins, except the contacts between region bins of opposite compartments. This is dense per
chromosome: use a resolution of 50 kb or coarser.

```text
usage: hicDifferentialAnalysis compartments --conditionA MATRIX [MATRIX ...]
                                            --conditionB MATRIX [MATRIX ...]
                                            --gcTrack BEDGRAPH
                                            --outFilePrefix PREFIX
                                            [--blocks LABEL [LABEL ...]]
                                            [--exploratory] [--fdr FDR]
                                            [--minFoldChange FOLD]
                                            [--filterThreshold LOWER UPPER]
                                            [--chromosomes CHROM [CHROM ...]]
                                            [--minDistance BP] [--threads THREADS]
                                            [--splitReplicates SEED]
                                            [--plantRegions BED] [--plantFold FOLD]
                                            [--plantSeed SEED] [-h]
```

### Required arguments

| Flag | Meaning |
|---|---|
| `--conditionA MATRIX [...], -a MATRIX [...]` | Cool files of condition A, one per replicate. |
| `--conditionB MATRIX [...], -b MATRIX [...]` | Cool files of condition B, one per replicate. |
| `--gcTrack BEDGRAPH` | GC content: chrom, start, end, value; averaged over each bin by overlap. |
| `--outFilePrefix PREFIX, -o PREFIX` | Prefix of the output file. |

### Optional arguments

| Flag | Meaning |
|---|---|
| `--blocks LABEL [...]` | One block label per sample (see `tads` above). |
| `--exploratory` | Allow a single-sample condition (see `tads` above). |
| `--fdr FDR` | Benjamini-Hochberg FDR at which a unit is called differential (default: 0.05). |
| `--minFoldChange FOLD` | Minimum fold change the test is against, TREAT (default: 1.1). |
| `--filterThreshold LOWER UPPER` | MAD z-score limits on cis coverage (default: -1.5 5.0). |
| `--chromosomes CHROM [...]` | Chromosomes to analyse (default: all). |
| `--threads THREADS, -t THREADS` | Worker threads (default: 4). |
| `--minDistance BP` | Least distance of the contacts that are counted (default: 200000). |

### Calibration arguments

Same as `tads` above: `--splitReplicates`, `--plantRegions`, `--plantFold`, `--plantSeed`.

## Calibration

`--splitReplicates`, `--plantRegions`, `--plantFold` and `--plantSeed` are calibration tools used to
validate the statistical test itself (null and planted-effect runs), documented per subcommand above;
they are not needed for a routine differential comparison.
