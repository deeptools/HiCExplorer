# hicDetectStripes

Detects architectural stripes with the method of Stripenn. New in v4.

```text
usage: hicDetectStripes --matrix MATRIX --outFileName OUTFILENAME
                        [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]
                        [--norm NORM] [--canny CANNY] [--minL MINL]
                        [--maxWidth MAXWIDTH]
                        [--maxPixelPercentiles P [P ...]]
                        [--pValue PVALUE] [--mask MASK]
                        [--blurFilter BLURFILTER] [--seed SEED]
                        [--fdr FDR] [--unfilteredFileName UNFILTERED]
                        [--threads THREADS] [--help] [--version]
```

Detects architectural stripes on the given contact matrix with the method of Stripenn (Yoon et al., Stripenn detects architectural stripes from chromatin conformation data using computer vision, Nature Communications 2022; Stripenn 1.1.65.22). HiCExplorer 3.7 has no stripe caller. The tool is Stripenn's algorithm in C++, with its defaults, and gives the same stripes, p-values and stripiness scores as Stripenn on the same matrix.

## Method

The contact image of each chromosome is saturated at several percentiles of its contact values (`--maxPixelPercentiles`) and scanned in overlapping frames of 200 bins. In every frame, Canny edge detection (`--canny`, the sigma of its Gaussian filter), a vertical line filter and a run-length scan of the edge columns (`--minL`) find pairs of edges at most `--maxWidth` bins apart, the sides of a candidate stripe. Overlapping candidates are reduced to the more elongated one. Each candidate gets an empirical p-value against a background distribution that is sampled with the seeded random generator (`--seed`) from all chromosomes of the matrix, and a stripiness score from its observed over expected contacts. A stripe is called when its p-value is below `--pValue`. No multiple-testing correction is applied, as in Stripenn; `--fdr` is an option.

The matrix is read with its balancing weights (`--norm weight`, a cool file with a `weight` column) as Stripenn reads it. As in Stripenn, the background distribution is estimated from all chromosomes of the matrix except unplaced contigs, chrM and chrY, whichever chromosomes are searched.

## Output

One line per called stripe, in the order of Stripenn's `result_filtered.tsv` (descending stripiness):

| Column | Meaning |
|---|---|
| `chrom` | Chromosome. |
| `anchor_start`, `anchor_end` | The narrow axis of the stripe (Stripenn's pos1 - 1 and pos2). |
| `orientation` | `horizontal` when the stripe box starts on the diagonal (Stripenn's downward stripe), `vertical` otherwise. |
| `extent_start`, `extent_end` | The span of the stripe box along the diagonal axis. |
| `enrichment` | Mean observed over expected contacts of the box. |
| `pvalue`, `qvalue` | Stripenn's empirical p-value and its Benjamini-Hochberg adjustment over all candidates. |

`--unfilteredFileName` writes all candidates with the columns of Stripenn's `result_unfiltered.tsv` (`chr`, `pos1` to `pos4`, `length`, `width`, `Mean`, `maxpixel`, `pvalue` and `Stripiness`).

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrix MATRIX, -m MATRIX` | The matrix to compute the stripe detection on (cool, mcool or h5). |
| `--outFileName OUTFILENAME, -o OUTFILENAME` | Outfile name to store the called stripes (tab separated, Stripenn's result_filtered order). |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--chromosomes CHROMOSOMES [CHROMOSOMES ...]` | Chromosomes to search for stripes (Stripenn's chrom). If not set, all chromosomes are used. The background distribution always uses all chromosomes of the matrix. |
| `--norm NORM` | Normalisation: 'weight' (the balancing weights of the cool file), the name of another column of its bin table, or 'None' for the raw counts (Stripenn's norm). (Default: weight). |
| `--canny CANNY` | Canny edge detection parameter, the sigma of its Gaussian filter (Stripenn's canny). (Default: 2.0). |
| `--minL MINL` | Minimum length of a stripe in bins (Stripenn's minL). (Default: 10). |
| `--minStripeLength MINSTRIPELENGTH` | The same in base pairs; when given it replaces --minL by MINSTRIPELENGTH divided by the bin size. |
| `--maxWidth MAXWIDTH` | Maximum width of a stripe in bins (Stripenn's maxW). (Default: 8). |
| `--maxPixelPercentiles P [P ...]` | Percentiles of the chromosome's contact values at which the image is saturated (Stripenn's maxpixel). (Default: 0.95 0.96 0.97 0.98 0.99). |
| `--pValue PVALUE` | P-value cutoff for a called stripe (Stripenn's pvalue). No multiple-testing correction is applied, as in Stripenn. (Default: 0.1). |
| `--mask MASK` | Column coordinates to be masked in the stripiness score, e.g. chr9:12345678-12345789 (Stripenn's mask). (Default: 0, none). |
| `--blurFilter BLURFILTER` | Mean filter size, an odd number (Stripenn's bfilter). (Default: 3). |
| `--seed SEED` | Seed of the pseudo random number generator that samples the background (Stripenn's seed, Python's random module). (Default: 123456789). |
| `--sequentialRandom` | Draw the background samples of all chromosomes from one generator in sequence. Stripenn does so when it runs on a single core; with several cores (and here by default) every chromosome starts from the seeded initial state. |
| `--fdr FDR` | Benjamini-Hochberg q-value threshold, applied instead of --pValue when given. Off by default (Stripenn applies none). |
| `--unfilteredFileName UNFILTEREDFILENAME` | Also write all candidate stripes with Stripenn's result_unfiltered.tsv columns. |
| `--threads THREADS, -t THREADS` | Number of threads to use. (Default: 4). |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## Agreement with Stripenn

The candidate stripes, their p-values and their stripiness scores are those of Stripenn 1.1.65.22 run on the same matrix with the same parameters and seed. Each stage (maxpixel quantiles, expected values, background distribution, edge images, boxes, redundancy removal, p-values and scores) is compared with Stripenn's own intermediate results by the scripts in `cpp/scripts/stripenn_oracle`. The order of the floating point operations of numpy, scipy, scikit-image and OpenCV is reproduced, so that the values are bit for bit equal; the remaining cases are listed in `cpp/STATUS.md`.

## Reference

Yoon S, Vo C, Chandra V, et al. Stripenn detects architectural stripes from chromatin conformation data using computer vision. Nature Communications 13, 2022.
