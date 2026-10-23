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
| `--plotData FILE` | C++ only: write the tables and the data of the figures as JSON to FILE instead of drawing them. |
| `--prettyReport` | C++ only: also write OUTPUTFOLDER/hicQC_report.html, a card-based report that judges five metrics per sample as PASS, WARN or FAIL. See "PASS, WARN and FAIL" below. |
| `--minValidPairsPercent PCT` | C++ only, read with --prettyReport. Valid pairs as a percentage of sequenced reads: PASS at or above PCT, WARN below PCT down to 0.8 * PCT, FAIL below 0.8 * PCT. Default 25.0. |
| `--maxInterChromosomalPercent PCT` | C++ only, read with --prettyReport. Inter-chromosomal contacts as a percentage of valid pairs: PASS at or below PCT, WARN above PCT up to 1.2 * PCT, FAIL above 1.2 * PCT. Default 10.0. |
| `--maxDuplicatedPairsPercent PCT` | C++ only, read with --prettyReport. Duplicated pairs as a percentage of mappable pairs, judged as for --maxInterChromosomalPercent. Default 20.0. |
| `--maxDanglingEndPercent PCT` | C++ only, read with --prettyReport. Dangling ends as a percentage of mappable pairs, judged as for --maxInterChromosomalPercent. Default 10.0. |
| `--maxSelfLigationPercent PCT` | C++ only, read with --prettyReport. Self-ligated pairs as a percentage of mappable pairs, judged as for --maxInterChromosomalPercent. Default 10.0. |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## Additional notes

C++ port: the tables are computed in C++, and the figures and hicQC.html are drawn by the hicexplorer_plot drawing layer with the calls of the Python tool (HICX_PLOT_PYTHON names the interpreter). The C++-only option --plotData FILE writes the tables and the data of the figures to FILE instead of drawing them.

## PASS, WARN and FAIL (--prettyReport)

Each metric is compared with its cutoff. A value on the right side of the
cutoff, or exactly on it, is PASS. A value past the cutoff by at most 20 % of
the cutoff is WARN, and a value further past it is FAIL. For the default
inter-chromosomal cutoff of 10 % this gives PASS up to 10 %, WARN above 10 % up
to 12 % and FAIL above 12 %; for the default valid pairs cutoff of 25 % it gives
PASS from 25 %, WARN from 20 % up to 25 % and FAIL below 20 %. The 20 % margin
is fixed and is not a command line option.

The valid pairs and inter-chromosomal cutoffs follow HiCExplorer's usage
guidance ("usually 25%-40% of the reads are valid"; "good Hi-C libraries have
lower than 10% inter chromosomal contacts"). The duplicated pairs, dangling end
and self-ligation cutoffs have no published HiCExplorer standard behind them
and are conservative starting points. hicQC.html and the QC tables are written
as without the option.

The report has a light and a dark theme (the button next to the title; the
query `?theme=light` or `?theme=dark` selects one when the page is opened) and
three colour schemes for the charts. The charts are drawn with Plotly 2.35.2,
loaded from cdn.plot.ly when the page is opened. Without network access or
JavaScript (for example where a server strips scripts from HTML outputs) the
charts are missing, while the verdict cards and every table still show.
