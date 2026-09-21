# hicQuickQC

Estimates the quality of a Hi-C dataset from a sample of read pairs.

```text
usage: hicQuickQC --samFiles two sam files two sam files --QCfolder FOLDER
                  --restrictionCutFile BED file [BED file ...]
                  --restrictionSequence RESTRICTIONSEQUENCE
                  [RESTRICTIONSEQUENCE ...] --danglingSequence DANGLINGSEQUENCE
                  [DANGLINGSEQUENCE ...] [--lines LINES] [--help] [--version]
```

The tool hicQuickQC considers the first n lines of two bam/sam files to get a first estimate of the quality of the data. It is highly recommended to set the restriction enzyme and dangling end parameter to get a good quality report.

## Required arguments

| Flag | Meaning |
|---|---|
| `--samFiles two sam files two sam files, -s two sam files two sam files` | The two PE alignment sam files to process. |
| `--QCfolder FOLDER` | Path of folder to save the quality control data of the matrix. |
| `--restrictionCutFile BED file [BED file ...], -rs BED file [BED file ...]` | BED file(s) with all restriction cut places. |
| `--restrictionSequence RESTRICTIONSEQUENCE [RESTRICTIONSEQUENCE ...], -seq ...` | Sequence of the restriction site. |
| `--danglingSequence DANGLINGSEQUENCE [DANGLINGSEQUENCE ...]` | Sequence left by the restriction enzyme after cutting. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--lines LINES` | Number of lines to consider for the QC test run (Default: 1000000). |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## Notes

For more information see [encode ](https://www.encodeproject.org/documents/75926e4b-77aa-4959-8ca7-87efcba39d79/@@download/attachment/comp_doc_7july2018_final.pdf) .
