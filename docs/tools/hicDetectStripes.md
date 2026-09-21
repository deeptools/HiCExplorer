# hicDetectStripes

Detects architectural stripes, a C++ port of Stripenn's detection method. New in v4.

```text
usage: hicDetectStripes --matrix MATRIX --outFileName OUTFILENAME
                        [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]
                        [--minStripeLength MINSTRIPELENGTH]
                        [--maxWidth MAXWIDTH] [--canny CANNY]
                        [--blurFilter BLURFILTER]
                        [--maxPixelPercentiles P [P ...]]
                        [--backgroundSamples BACKGROUNDSAMPLES]
                        [--pValue PVALUE] [--fdr FDR] [--seed SEED]
                        [--threads THREADS] [--help] [--version]
```

Detects architectural stripes on the given contact matrix. New in HiCExplorer v4 (PLAN.md tier 9 section 9.3): a faithful C++ reimplementation of Stripenn 1.1.65.22's detection method (Yoon et al. 2022), not an invented stand-in.

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrix MATRIX, -m MATRIX` | The matrix to compute the stripe detection on. |
| `--outFileName OUTFILENAME, -o OUTFILENAME` | Outfile name to store the detected stripes (tab separated). |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--chromosomes CHROMOSOMES [CHROMOSOMES ...]` | Chromosomes to include in the analysis. If not set, all chromosomes are included. |
| `--minStripeLength MINSTRIPELENGTH` | Shortest candidate stripe length, in base pairs (Stripenn's minL, converted to bins). (Default: 100000). |
| `--maxWidth MAXWIDTH` | Maximum stripe width, in bins (Stripenn's maxW). (Default: 8). |
| `--canny CANNY` | Canny edge detection sigma (Stripenn's canny). (Default: 2.0). |
| `--blurFilter BLURFILTER` | Mean filter size, an odd number (Stripenn's bfilter). (Default: 3). |
| `--maxPixelPercentiles P [P ...]` | Percentiles of the contact frequency data to saturate the image (Stripenn's maxpixel). (Default: 0.95 0.96 0.97 0.98 0.99). |
| `--backgroundSamples BACKGROUNDSAMPLES` | Random samples per chromosome for the background model the p-value is ranked against. (Default: 200000). |
| `--pValue PVALUE` | Raw p-value cutoff for the final call set (Stripenn's own pvalue, and its own practice: no multiple-testing correction by default -- see PLAN.md 9.3 for the real-data evidence this is based on). (Default: 0.1). |
| `--fdr FDR` | Benjamini-Hochberg q-value threshold, applied instead of the raw --pValue cutoff when given. Off by default: on real GM12878 data, even Stripenn's own p-values do not survive q <= 0.05 after correction (PLAN.md 9.3), so this option exists for a user who wants that stricter, opt-in guarantee, not as the default. (Default: off). |
| `--seed SEED` | Seed for the background sample. (Default: 20260915). |
| `--threads THREADS, -t THREADS` | Number of threads to use. (Default: 4). |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |
