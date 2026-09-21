# hicDetectLoops

Detects enriched interaction regions (peaks/loops).

```text
usage: hicDetectLoops --matrix MATRIX --outFileName OUTFILENAME
                      [--peakWidth PEAKWIDTH] [--windowSize WINDOWSIZE]
                      [--pValuePreselection PVALUEPRESELECTION]
                      [--peakInteractionsThreshold PEAKINTERACTIONSTHRESHOLD]
                      [--obsExpThreshold OBSEXPTHRESHOLD] [--pValue PVALUE]
                      [--maxLoopDistance MAXLOOPDISTANCE]
                      [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]
                      [--threads THREADS]
                      [--threadsPerChromosome THREADSPERCHROMOSOME]
                      [--expected {mean,mean_nonzero,mean_nonzero_ligation}]
                      [--help] [--version]
```

Computes enriched regions (peaks) or long range contacts on the given contact matrix.

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrix MATRIX, -m MATRIX` | The matrix to compute the loop detection on. |
| `--outFileName OUTFILENAME, -o OUTFILENAME` | Outfile name to store the detected loops. The file will in bedgraph format. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--peakWidth PEAKWIDTH, -pw PEAKWIDTH` | The width of the peak region in bins. (Default: 2). |
| `--windowSize WINDOWSIZE, -w WINDOWSIZE` | The window size for the neighborhood region the peak is located in. (Default: 5). |
| `--pValuePreselection PVALUEPRESELECTION, -pp PVALUEPRESELECTION` | Only candidates with p-values less the given threshold will be considered as candidates. Can a single value or a threshold file created by hicCreateThresholdFile (Default: 0.1). |
| `--peakInteractionsThreshold PEAKINTERACTIONSTHRESHOLD, -pit PEAKINTERACTIONSTHRESHOLD` | The minimum number of interactions a detected peaks needs to have to be considered (Default: 10). |
| `--obsExpThreshold OBSEXPTHRESHOLD, -oet OBSEXPTHRESHOLD` | The minimum number of obs/exp interactions a detected peaks needs to have to be considered (Default: 1.5). |
| `--pValue PVALUE, -p PVALUE` | Rejection level for Anderson-Darling or Wilcoxon-rank sum test for H0. (Default: 0.025). |
| `--maxLoopDistance MAXLOOPDISTANCE` | Maximum genomic distance of a loop (Default: 2000000). |
| `--chromosomes CHROMOSOMES [CHROMOSOMES ...]` | Chromosomes to include in the analysis. If not set, all chromosomes are included. |
| `--threads THREADS, -t THREADS` | Number of threads to use (Default: 4). |
| `--threadsPerChromosome THREADSPERCHROMOSOME, -tpc THREADSPERCHROMOSOME` | Number of threads to use per parallel thread processing a chromosome (Default: 4). |
| `--expected {mean,mean_nonzero,mean_nonzero_ligation}, -exp {mean,mean_nonzero,mean_nonzero_ligation}` | Method to compute the expected value per distance (Default: mean). |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |
