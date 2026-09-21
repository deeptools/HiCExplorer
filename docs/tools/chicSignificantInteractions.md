# chicSignificantInteractions

Detects significant interactions per viewpoint from the background model.

```text
usage: chicSignificantInteractions --interactionFile INTERACTIONFILE --pValue
                                   PVALUE
                                   [--xFoldBackground XFOLDBACKGROUND | --loosePValue LOOSEPVALUE]
                                   --backgroundModelFile BACKGROUNDMODELFILE
                                   --range RANGE RANGE
                                   [--outFileNameSignificant OUTFILENAMESIGNIFICANT]
                                   [--outFileNameTarget OUTFILENAMETARGET]
                                   [--combinationMode {dual,single}]
                                   [--threads THREADS] [--truncateZeroPvalues]
                                   [--fixateRange FIXATERANGE]
                                   [--peakInteractionsThreshold PEAKINTERACTIONSTHRESHOLD]
                                   [--correctForMultipleTesting {none,fdr,bonferroni}]
                                   [--help] [--version]
```

Per viewpoint the significant interactions are detected based on the background model. For each viewpoint file, an output file is created with all recorded significant interactions and a target file. The target file is especially useful in the batch mode context; for two consecutive listed control and treatment viewpoints it merges the significant interactions which can then be used to test for a differential interaction scheme.

chicSignificantInteractions supports two modes to detect significant interactions, either by an x-fold over the average background or a loose p-value. In both cases neighboring significant peaks are merged together and an additional p-value is computed based on the sum of interactions for this neighborhood. Only interactions with a higher p-value (as specified by the threshold `--pValue`) are accepted as a significant interaction.

## options

| Flag | Meaning |
|---|---|
| `--xFoldBackground XFOLDBACKGROUND, -xf XFOLDBACKGROUND` | Filter x-fold over background. Used to merge neighboring bins with a broader peak but less significant interactions to a single peak with high significance. Used only for pValue option. |
| `--loosePValue LOOSEPVALUE, -lp LOOSEPVALUE` | loose p-value threshold to filter target regions in a first round. Used to merge neighboring bins with a broader peak but less significant interactions to a single peak with high significance. Used only for pValue option. |

## Required arguments

| Flag | Meaning |
|---|---|
| `--interactionFile INTERACTIONFILE, -if INTERACTIONFILE` | path to the interaction file (HDF5) which should be used for aggregation of the statistics. |
| `--pValue PVALUE, -p PVALUE` | p-value threshold to filter target regions for inclusion in differential analysis. |
| `--backgroundModelFile BACKGROUNDMODELFILE, -bmf BACKGROUNDMODELFILE` | path to the background file. |
| `--range RANGE RANGE` | Defines the region upstream and downstream of a reference point which should be included. Format is --region upstream downstream, e.g. --region 500000 500000 plots 500kb up- and 500kb downstream. This value should not exceed the range used in the other chic-tools. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--outFileNameSignificant OUTFILENAMESIGNIFICANT, -os OUTFILENAMESIGNIFICANT` | File name suffix to save the results; prefix is the input file name (Default: significantInteractions.hdf5). |
| `--outFileNameTarget OUTFILENAMETARGET, -ot OUTFILENAMETARGET` | The file to store the target data (Default: targetFile.hdf5). |
| `--combinationMode {dual,single}, -cm {dual,single}` | This option defines how the interaction data should be computed and combined: dual: Combines as follows: [[matrix1_gene1, matrix2_gene1], [matrix2_gene1, matrix3_gene1],[matrix1_gene2, matrix2_gene2], ...]single: Combines as follows: [matrix1_gene1, matrix1_gene2, matrix2_gene1, ...], (Default: dual). |
| `--threads THREADS, -t THREADS` | Number of threads (uses the python multiprocessing module) (Default: 4). |
| `--truncateZeroPvalues, -tzpv` | Sets all p-values which are equal to zero to one. This has the effect that the associated positions are not part of the significance decision. |
| `--fixateRange FIXATERANGE, -fs FIXATERANGE` | Fixate range of backgroundmodel starting at distance x. E.g. all values greater than 500kb are set to the value of the 500kb bin (Default: 500000). |
| `--peakInteractionsThreshold PEAKINTERACTIONSTHRESHOLD, -pit PEAKINTERACTIONSTHRESHOLD` | The minimum number of interactions a detected peak needs to have to be considered (Default: 5). |
| `--correctForMultipleTesting {none,fdr,bonferroni}` | Adjust the tested p-values across all viewpoints of all samples, Benjamini-Hochberg (fdr) or Bonferroni, and apply --pValue to the adjusted values; the significant file gains pvalue_adjusted. Not in the Python tool; none gives its output (Default: none). |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |
