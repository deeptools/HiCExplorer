# hicPrepareQCreport

Tabulates and plots QC measures from multiple hicBuildMatrix log files into one HTML report.

```text
usage: hicPrepareQCreport --logfiles matrix1_QCfolder/QC.log matrix2_QCfolder/QC.log --labels "sample 1" "sample 2" --outputFolder QC_all_samples
```

Tabulates and plots QC measures from hicBuildMatrix log files within an HTML output

## Required arguments

| Flag | Meaning |
|---|---|
| `--logfiles LOGFILES [LOGFILES ...], -l LOGFILES [LOGFILES ...]` | Path to the log files to be processed |
| `--labels LABELS [LABELS ...]` | Label to assign to each log file. Each label should be separated by a space. Quote labels that contain spaces: E.g. --labels label1 "labels 2" |
| `--outputFolder OUTPUTFOLDER, -o OUTPUTFOLDER` | Several files with be saved under this folder: A table containing the results and a html file with several images. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--dpi DPI` | Image resolution. By default high resolution png images with a 200 dpi are created. |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## Additional notes

C++ port: the tables are computed in C++, and the figures and hicQC.html are drawn by the hicexplorer_plot drawing layer with the calls of the Python tool (HICX_PLOT_PYTHON names the interpreter). The C++-only option --plotData FILE writes the tables and the data of the figures to FILE instead of drawing them.
