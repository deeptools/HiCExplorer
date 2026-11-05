# hicValidateLocations

Compares called loops with known protein peak positions.

```text
usage: hicValidateLocations --data DATA --validationData VALIDATIONDATA
                            [--validationType {bed,cool}]
                            [--method {loops,tad}] --resolution RESOLUTION
                            [--outFileName OUTFILENAME]
                            [--chrPrefixLoops {None,add,remove}]
                            [--chrPrefixProtein {None,add,remove}] [--help]
                            [--version]
```

This script overlaps the loop locations with protein locations to determine the accuracy of the loop detection.

## Required arguments

| Flag | Meaning |
|---|---|
| `--data DATA, -d DATA` | The loop file from hicDetectLoops. To use files from other sources, please follow 'chr start end chr start end' format. For TAD data use the boundaries.bed file and not the domains file! |
| `--validationData VALIDATIONDATA, -vd VALIDATIONDATA` | The data file to validate the given locations. Can be narrowPeak, broadPeak (both in bed), or cool |
| `--validationType {bed,cool}, -vt {bed,cool}` | The type of the validation data. Can be bed, or cool format |
| `--method {loops,tad}, -m {loops,tad}` | The method used (for the moment only loop is possible) (Default: loops). |
| `--resolution RESOLUTION, -r RESOLUTION` | The used resolution of the Hi-C interaction matrix. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--outFileName OUTFILENAME, -o OUTFILENAME` | The prefix name of the output files. Two file are written: output_matched_locations and output_statistics.First file contains all loop locations with protein location matches, second file contains statistics about this matching. |
| `--chrPrefixLoops {None,add,remove}, -cl {None,add,remove}` | Adding / removing / do nothing a 'chr'-prefix to chromosome name of the loops. |
| `--chrPrefixProtein {None,add,remove}, -cp {None,add,remove}` | Adding / removing / do nothing a 'chr'-prefix to chromosome name of the protein. |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## Loops need to have format as follows

`chr start end chr start end`

## Additional notes

The protein peaks need to be in narrowPeaks or broadPeak format.

A protein match is successfull if at the bin of the x and y location a protein peak is overlapped. A bin is assumed to have a protein if one or more protein peaks falling within the bin region. The value of the protein is not considered, only match or non-match.
