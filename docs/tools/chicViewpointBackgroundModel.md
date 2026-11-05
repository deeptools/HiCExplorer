# chicViewpointBackgroundModel

Computes the background model used by the other capture Hi-C tools.

```text
usage: chicViewpointBackgroundModel --matrices MATRICES [MATRICES ...]
                                    --referencePoints REFERENCEPOINTS
                                    [--averageContactBin AVERAGECONTACTBIN]
                                    [--truncateZeros]
                                    [--outFileName OUTFILENAME]
                                    [--threads THREADS]
                                    [--fixateRange FIXATERANGE] [--help]
                                    [--version]
```

chicViewpointBackgroundModel computes a background model for all given samples with all reference points. For all relative distances to a reference point a negative binomial distribution is fitted. In addition, for each relative distance to a reference point the average value for this location is computed. Both background models are used, the first one for p-value and significance computation, the second one to filter out interactions with a smaller x-fold over the mean.

The background distributions are fixed at `--fixateRange`, i.e. all distances lower or higher than this value use the fixed background distribution.

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]` | The input matrices (samples) to build the background model on. |
| `--referencePoints REFERENCEPOINTS, -rp REFERENCEPOINTS` | Bed file contains all reference points which should be used to build the background model. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--averageContactBin AVERAGECONTACTBIN` | Average the contacts of n bins via a sliding window approach (Default: 5). |
| `--truncateZeros, -tz` | Truncates the zeros before the distributions are fitted. Use it in case you observe an over dispersion. |
| `--outFileName OUTFILENAME, -o OUTFILENAME` | The name of the background model file (Default: background_model.txt). |
| `--threads THREADS, -t THREADS` | Number of threads (uses the python multiprocessing module) (Default: 4). |
| `--fixateRange FIXATERANGE, -fs FIXATERANGE` | Fixate score of backgroundmodel starting at distance x. E.g. all values greater 500kb are set to the value of the 500kb bin (Default: 500000). |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## An example usage is

$ chicViewpointBackgroundModel --matrices matrix1.cool matrix2.cool matrix3.cool --referencePoints referencePointsFile.bed --range 20000 40000 --outFileName background_model.bed
