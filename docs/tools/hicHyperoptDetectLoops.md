# hicHyperoptDetectLoops

Searches for the best [hicDetectLoops](hicDetectLoops.md) parameter setting for a given dataset, using
Bayesian hyperparameter optimization (the Python `hyperopt` library) against a set of known protein peak
locations. Because `hicDetectLoops` has many parameters and finding a good setting by hand is difficult,
HiCExplorer added this tool (and its HiCCUPS counterpart,
[hicHyperoptDetectLoopsHiCCUPS](hicHyperoptDetectLoopsHiCCUPS.md)) in version 3.5.

```text
usage: hicHyperoptDetectLoops --matrix MATRIX --proteinFile PROTEINFILE
                              --maximumNumberOfLoops MAXIMUMNUMBEROFLOOPS
                              [--outputFileName OUTPUTFILENAME] [--resolution RESOLUTION]
                              [--chrPrefixLoops {None,add,remove}] [--threads THREADS]
                              [--runs RUNS] [--help] [--version]
```

The tool runs in Python. Each evaluation of a parameter setting runs [hicDetectLoops](hicDetectLoops.md)
and [hicValidateLocations](hicValidateLocations.md).

| Flag | Meaning |
|---|---|
| `--matrix MATRIX, -m MATRIX` | The matrix to compute the loops on. |
| `--proteinFile PROTEINFILE, -p PROTEINFILE` | The protein file to validate the detected loops. |
| `--maximumNumberOfLoops MAXIMUMNUMBEROFLOOPS, -ml MAXIMUMNUMBEROFLOOPS` | The maximum number of loops that should be used for the optimization. |
| `--outputFileName OUTPUTFILENAME, -o OUTPUTFILENAME` | File name for the result of the optimization (Default: hyperopt_result.txt). |
| `--resolution RESOLUTION, -re RESOLUTION` | Resolution of the matrix (Default: 10000). |
| `--chrPrefixLoops {None,add,remove}, -cl {None,add,remove}` | Add, remove or keep a 'chr' prefix on the chromosome names of the loops. |
| `--threads THREADS, -t THREADS` | Number of threads (Default: 4). |
| `--runs RUNS, -r RUNS` | Number of hyperopt runs (Default: 100). |
