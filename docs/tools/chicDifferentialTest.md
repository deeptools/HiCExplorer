# chicDifferentialTest

Tests for differential interactions between two samples (chi2 or Fisher).

```text
usage: chicDifferentialTest --aggregatedFile AGGREGATEDFILE --alpha ALPHA
                            [--outFileName OUTFILENAME]
                            [--statisticTest {fisher,chi2}]
                            [--threads THREADS]
                            [--correctForMultipleTesting {none,fdr,bonferroni}]
                            [--help] [--version]
```

chicDifferentialTest tests if two locations under consideration of the reference point have a different interaction count. For this either Fisher's test or the chi2 contingency test can be used. The file that is accepted for this test can be created with `chicAggregateStatistic`. H0 assumes the interactions are not different. Therefore the differential interaction counts are all where H0 was rejected.

## Required arguments

| Flag | Meaning |
|---|---|
| `--aggregatedFile AGGREGATEDFILE, -af AGGREGATEDFILE` | path to the aggregated files which should be used for the differential test. |
| `--alpha ALPHA, -a ALPHA` | define a significance level (alpha) for accepting samples |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--outFileName OUTFILENAME, -o OUTFILENAME` | Output file for the differential results (Default: differentialResults.hdf5). |
| `--statisticTest {fisher,chi2}` | Type of test used: fisher's exact test or chi2 contingency (Default: fisher). |
| `--threads THREADS, -t THREADS` | Number of threads (uses the python multiprocessing module) (Default: 4). |
| `--correctForMultipleTesting {none,fdr,bonferroni}` | Adjust the p-values of all tested locations of all reference points, Benjamini-Hochberg (fdr) or Bonferroni, and reject where the adjusted value is at most --alpha; the result groups gain pvalue_adjusted_list. Not in the Python tool; none gives its output (Default: none). |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |
