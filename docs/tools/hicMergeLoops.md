# hicMergeLoops

Merges loop calls made at different resolutions.

```text
usage: hicMergeLoops --inputFiles INPUTFILES [INPUTFILES ...] --outFileName
                     OUTFILENAME --lowestResolution LOWESTRESOLUTION [--help]
                     [--version]
```

This script merges the locations of loops detected at several resolutions.

## Required arguments

| Flag | Meaning |
|---|---|
| `--inputFiles INPUTFILES [INPUTFILES ...], -i INPUTFILES [INPUTFILES ...]` | The loop files from hicDetectLoops. To use files from other sources, please follow 'chr start end chr start end' format and remove any header. |
| `--outFileName OUTFILENAME, -o OUTFILENAME` | The name of the merged loop file. |
| `--lowestResolution LOWESTRESOLUTION, -r LOWESTRESOLUTION` | The lowest resolution of all loop files, i.e. 5kb, 10kb and 25kb, please use 25000. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## Loops in the inputFiles need to have the following format

chr start end chr start end

## Additional notes

Loops are merged if the x and y position of a loop overlap with the x and y position of another loop; all loops are considered as an overlap within +/- the bin size of the lowest resolution. I.e. for a loop with coordinates x and y, the overlap with all other loops is checked for (x - lowest resolution) and (y + lowest resolution). If two or more locations are to be merged, the loop at the lowest resolution is taken as the merged loop.
