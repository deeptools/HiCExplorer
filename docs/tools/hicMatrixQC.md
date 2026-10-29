# hicMatrixQC

Quality-control report for cool and mcool matrices. New in version 4; there is no counterpart in version 3.

```text
usage: hicMatrixQC --matrices MATRICES [MATRICES ...] --outFileName OUTFILENAME
                   [--tableFileName TABLEFILENAME] [--labels LABELS [LABELS ...]]
                   [--resolutions RESOLUTIONS [RESOLUTIONS ...]]
                   [--shortRangeCutoff BP] [--decayRange LOW HIGH]
                   [--mapResolutionContacts N] [--mapResolutionFraction F]
                   [--madThreshold LOWER UPPER] [--transThreshold X]
                   [--maxChromosomes N] [--weightColumn NAME] [--threads INT]
                   [--help] [--version]
```

`hicMatrixQC` does for a finished matrix what the QC of `hicBuildMatrix` does for the BAM files. It writes an HTML
report and a tab separated table of all statistics with one column per sample. An mcool file gives one sample per
resolution. Every matrix is read once, in chunks, so memory use depends on the number of bins and not on the number
of contacts; each thread holds two vectors with one value per bin. The report is plain HTML that needs no script and
reads at narrow widths; Galaxy's HTML sanitizer leaves it intact.

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]` | The cool or mcool matrices. |
| `--outFileName OUTFILENAME, -o OUTFILENAME` | File name of the HTML report. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--tableFileName TABLEFILENAME, -t TABLEFILENAME` | File name of the table of all statistics. Default: the report's name with the extension `.tsv`. |
| `--labels LABELS [LABELS ...]` | One label per matrix. Default: the file names. |
| `--resolutions RESOLUTIONS [RESOLUTIONS ...], -r RESOLUTIONS [RESOLUTIONS ...]` | The resolutions of an mcool file to analyse. Default: all resolutions present. Ignored for a cool file. |
| `--shortRangeCutoff BP` | Cis contacts closer than this distance (bp) are short range (default: 20000). |
| `--decayRange LOW HIGH` | Distances (bp) of the decay exponent fit (default: 100000 1000000). |
| `--mapResolutionContacts N` | Minimum number of contacts of a bin for the map resolution (default: 1000). |
| `--mapResolutionFraction F` | Fraction of the bins that must reach that number (default: 0.8). |
| `--madThreshold LOWER UPPER` | Modified z-scores of the coverage outside this range are removed by the MAD filter of hicCorrectMatrix (default: -1.5 5.0). |
| `--transThreshold X` | Chromosome pairs with an observed over expected ratio of trans contacts above X are flagged (default: 2.0). |
| `--maxChromosomes N` | The N longest chromosomes get rows of their own; the genome-wide statistics use all (default: 30). |
| `--weightColumn NAME` | The bin column of the balancing weights (default: weight). |
| `--threads INT, -p INT` | Number of threads (default: 4). |
| `--help, -h` | Show this help message and exit. |
| `--version` | Show program's version number and exit. |

## An example usage is

```bash
$ hicMatrixQC --matrices matrix.cool --outFileName qc.html
$ hicMatrixQC --matrices sample1.mcool sample2.mcool --labels control treated --resolutions 10000 \
    --outFileName qc.html --tableFileName qc.tsv
```

## Statistics

Counts are the counts of the stored pixels (raw, not balanced). The distance of a pixel is the distance between the
starts of its two bins. A pixel at offset d bins holds read pairs between d - 1 and d + 1 bin sizes apart, so every
distance threshold is exact only up to one bin size, and the match to the read-level numbers of `QC.log` is closer
for finer bins. The pairs inside one bin (the main diagonal) are closer than one bin size and are counted as short
range.

1. **Library statistics**, when `hicBuildMatrix` stored them: every count of `QC.log`, the parameters of the run, the
   distribution of the intra-chromosomal read pairs over distance, and the same split by read pair type (inward,
   outward, left, right). The distance is the distance between the read start positions, measured before binning.
2. **Contacts**: total, cis and trans fractions, cis short range and long range, the fractions of the main diagonal
   and of the first off-diagonal.
3. **Contact distance curve P(s)**: the contacts at each distance divided by the number of bin pairs at that distance,
   genome-wide, on the log-spaced distance bins of the QC tables of `hicBuildMatrix` (ten per decade). The decay
   exponent is the slope of log10 P(s) against log10 s between the two distances of `--decayRange`, fitted to ten
   groups of distances equally spaced on the log scale. It is given genome-wide and per chromosome.
4. **Map resolution** after Rao et al. (2014): the smallest bin size at which `--mapResolutionFraction` of the bins
   hold at least `--mapResolutionContacts` contacts (the row sum of the symmetric matrix). For an mcool file the
   resolutions present are used, for a single cool file the bins are aggregated to multiples of the bin size inside
   each chromosome. The fraction is taken over all bins and, in a second column, over the bins with a contact. When the
   smallest bin size of the file already reaches the fraction, the map resolution is that bin size or finer.
5. **Coverage per bin**: quantiles over the non-empty bins, the fraction of empty bins, the fraction of bins that the
   filter of `hicCorrectMatrix` (modified z-score of the coverage without the diagonal) would remove with
   `--madThreshold`, and the coverage along each chromosome as a row of coloured cells.
6. **Balancing**: whether the weight column exists, the attributes cooler stores on it (`converged`,
   `ignore_diags`, and so on), the fraction of bins with a NaN weight, and the coefficient of variation of the
   balanced marginals over the bins with a finite weight.
7. **Chromosome pairs**: the observed trans contacts of each pair of chromosomes divided by the expected ones. The
   expectation distributes the trans contacts between the listed chromosomes in proportion to the product of the two
   chromosomes' trans marginals and is normalised to the observed total. Pairs above `--transThreshold` are flagged.
8. **Further checks**: the ratio of the mean coverage per non-empty bin of chromosome X to that of the autosomes,
   whether chromosome Y and chromosome M exist and their share of the coverage. Names are matched with or without the
   prefix `chr`.

A matrix with bins of different sizes (for example restriction fragments) gets the library statistics only.

## Notes

The table has one row per statistic, named in groups (`contacts.*`, `decay.*`, `coverage.*`, `map_resolution.*`,
`balance.*`, `trans.*`, `checks.*`, `qc.*`), and an empty cell for a statistic that does not exist for a sample.
