# chicChicagoSignificantInteractions

Filters the output of [chicChicagoScores](chicChicagoScores.md) by CHiCAGO's score threshold and writes the accepted interactions with the genomic coordinates of both fragments. Third step of the CHiCAGO pipeline.

This tool is new in version 4. It is unrelated to [chicSignificantInteractions](chicSignificantInteractions.md), which works on the HDF5 files of the chicViewpoint workflow.

```text
usage: chicChicagoSignificantInteractions --scores SCORES --rmap RMAP
                                          --baitmap BAITMAP
                                          [--outFileName OUTFILENAME]
                                          [--scoreThreshold SCORETHRESHOLD]
                                          [--threads THREADS] [--help] [--version]
```

Cairns et al. 2016 call an interaction significant at a score of at least 5, which is the default here. The output has one row per accepted interaction with the columns `bait_chr`, `bait_start`, `bait_end`, `bait_name`, `otherEnd_chr`, `otherEnd_start`, `otherEnd_end`, `N`, `distSign` and `score`. Coordinates are those of the `.rmap` and `.baitmap`.

## Required arguments

| Flag | Meaning |
|---|---|
| `--scores SCORES` | Output file of chicChicagoScores. |
| `--rmap RMAP` | CHiCAGO `.rmap` file (other-end coordinates). |
| `--baitmap BAITMAP` | CHiCAGO `.baitmap` file (bait coordinates and names). |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--outFileName OUTFILENAME, -o OUTFILENAME` | Name of the significant-interactions file (Default: chicago_significant_interactions.txt). |
| `--scoreThreshold SCORETHRESHOLD` | Minimum CHiCAGO score to call an interaction significant (Default: 5). |
| `--threads THREADS` | Number of threads (Default: 1). The output does not depend on it. |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## Example

```bash
chicChicagoSignificantInteractions --scores scores.txt --rmap design.rmap \
    --baitmap design.baitmap -o significant.txt
```

```text
bait_chr  bait_start  bait_end  bait_name           otherEnd_chr  otherEnd_start  otherEnd_end  N   distSign  score
20        325598      341635    NRSN2;RP5-1103G7.4  20            459227          463619        17  127806    13.7655046998409
20        325598      341635    NRSN2;RP5-1103G7.4  20            463620          473617        8   135002    5.81135847913031
```
