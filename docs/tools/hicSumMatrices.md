# hicSumMatrices

Adds Hi-C matrices of the same size.

```text
usage: hicSumMatrices --matrices .h5 or cooler file format
                      [.h5 or cooler file format ...] --outFileName
                      OUTFILENAME [-h] [--version]
```

Adds Hi-C matrices of the same size. Format has to be hdf5 (.h5) or npz. In order to minimize the the loss of information, it is recommended to to sum uncorrected matrices (before hicCorrectMatrix).

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrices .h5 or cooler file format [.h5 or cooler file format ...], -m .h5 or cooler file format [.h5 or cooler file format ...]` | Space-delimited names of the matrices to add. The matrices must have the same shape/size. You can verify their size by using `hicInfo`. (default: None) |
| `--outFileName OUTFILENAME, -o OUTFILENAME` | File name to save the resulting matrix. The output is from the same file type as the input. Please add the file ending suffix (either .h5 or .cool), if it is not given, there will be no output. (default: None) |

## Optional arguments

| Flag | Meaning |
|---|---|
| `-h, --help` | show this help message and exit |
| `--version` | show program's version number and exit |
