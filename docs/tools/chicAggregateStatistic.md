# chicAggregateStatistic

Aggregates viewpoint data of two samples into targets for the differential test.

```text
usage: chicAggregateStatistic --interactionFile INTERACTIONFILE
                              [--targetFile TARGETFILE]
                              [--outFileName OUTFILENAME] [--threads THREADS]
                              [--help] [--version]
```

chicAggregateStatistic is a preprocessing tool for chicDifferentialTest. It takes two consecutive viewpoint files and one target file and creates one file containing all locations which should be tested for differential interactions. Either one target file for two consecutive viewpoint files or one target file for all viewpoints is accepted.

## Required arguments

| Flag | Meaning |
|---|---|
| `--interactionFile INTERACTIONFILE, -if INTERACTIONFILE` | path to the interaction files which should be used for aggregation of the statistics. |
| `--targetFile TARGETFILE, -tf TARGETFILE` | path to the target files which contains the target regions to prepare data for differential analysis. This is either the target file in the hdf format created by chicSignificantInteractions or a regular, three column bed file. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--outFileName OUTFILENAME, -o OUTFILENAME` | File name to save the result (Default: aggregate_target.hdf). |
| `--threads THREADS, -t THREADS` | Number of threads (uses the python multiprocessing module)ist (Default: 4). |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |
