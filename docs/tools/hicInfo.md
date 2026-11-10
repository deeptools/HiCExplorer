# hicInfo

Prints information about one or more Hi-C matrices (bins, bin size, sum, min, max, and so on).

```text
usage: hicInfo --matrices MATRICES [MATRICES ...] [--outFileName OUTFILENAME]
               [--no_metadata] [--help] [--version]
```

Prints information about a matrix or matrices including matrix size, number of elements, sum of elements, etc.

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]` | The matrix (or multiple matrices) to get information about. HiCExplorer supports the following file formats: h5 (native HiCExplorer format) and cool. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--outFileName OUTFILENAME, -o OUTFILENAME` | File name to save information of the matrix instead of writing it to the bash. |
| `--no_metadata, -nm` | Do not use meta data from cooler file to display information. This method is slower and was the default until version 2.2 of HiCExplorer. H5 files always use this parameter. |
| `--help, -h` | Show this help message and exit. |
| `--version, -v` | show program's version number and exit |

## An example usage is

$ hicInfo -m matrix1.h5 matrix2.h5 matrix3.h5

## Notes

An example output looks like this:

```INI
File:	hic_matrix.h5
WARNING: bin size is not homogeneous. Median 434
Size:	339,684
Sum:	243,779,407
Bin length:	434
Chromosomes:	2L, 2R, 3L, 3R, 4, X
Non-zero elements:	229,122,848
Minimum:	0
Maximum:	2116
NaN bins:	0
```
With HiCExplorer version 3 we support metadata in cooler files. The default behavior is to use the metadata:

```bash
$ hicInfo -m matrix.cool
```
To use the old method (and default for h5) please add the parameter:

```bash
$ hicInfo -m matrix.cool --no_metadata
```
