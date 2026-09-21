# hicMergeTADbins

Uses a BED file of domains or TAD boundaries to merge the bin counts of a Hi-C matrix.

```text
usage: hicMergeTADbins [-h] --matrix MATRIX --domains DOMAINS --outFile
                       OUTFILE [--version]
```

Uses a BED file of domains or TAD boundaries to merge the bin counts of a Hi-C matrix per TAD.

The output matrix contains the total counts per TAD and the total contacts with all other TADs.

## options

| Flag | Meaning |
|---|---|
| `-h, --help` | show this help message and exit |
| `--matrix MATRIX, -m MATRIX` | Path to Hi-C matrix to use. |
| `--domains DOMAINS` | Path to a bed file containing the domains. |
| `--outFile OUTFILE, -o OUTFILE` | Name for the resulting matrix file. |
| `--version` | show program's version number and exit |
