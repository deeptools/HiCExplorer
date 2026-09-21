# chicViewpoint

Computes one viewpoint file per sample per reference point, based on the background model.

```text
usage: chicViewpoint --matrices MATRICES [MATRICES ...] --range RANGE RANGE
                     --referencePoints REFERENCEPOINTS --backgroundModelFile
                     BACKGROUNDMODELFILE [--outFileName OUTFILENAME]
                     [--threads THREADS]
                     [--averageContactBin AVERAGECONTACTBIN]
                     [--fixateRange FIXATERANGE] [--help] [--version]
```

Computes per input matrix all viewpoints which are defined in the reference points file.

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]` | Path to the Hi-C matrices which store the captured Hi-C data per sample. |
| `--range RANGE RANGE` | Defines the region upstream and downstream of a reference point which should be considered in the analysis. Please remember to use the same fixate range setting as for the background model computation and that distances of the range larger than the fixate range use the background model of those.Format is --region upstream downstream |
| `--referencePoints REFERENCEPOINTS, -rp REFERENCEPOINTS` | Reference point file. Needs to be in the format: 'chr 100' for a single reference point or 'chr 100 200' for a reference region and with a single reference point per line |
| `--backgroundModelFile BACKGROUNDMODELFILE, -bmf BACKGROUNDMODELFILE` | path to the background file computed by chicViewpointBackgroundModel |
| `--outFileName OUTFILENAME, -o OUTFILENAME` | This hdf5 file contains all created viewpoint files. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--threads THREADS, -t THREADS` | Number of threads (uses the python multiprocessing module) (Default: 4). |
| `--averageContactBin AVERAGECONTACTBIN` | Average the contacts of n bins via a sliding window approach to smooth the values and be less sensitive for outliers (Default: 5). |
| `--fixateRange FIXATERANGE, -fs FIXATERANGE` | Fixate range of background model starting at distance x. E.g. all values greater 500kb are set to the value of the 500kb bin (Default: 500000). |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |
