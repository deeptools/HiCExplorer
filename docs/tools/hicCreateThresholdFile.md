# hicCreateThresholdFile

Writes a distance-threshold file (one row per bin over a range) for hicDetectLoops.

```text
usage: hicCreateThresholdFile --thresholdValue THRESHOLDVALUE --range RANGE
                              RANGE [--resolution RESOLUTION] --outFileName
                              OUTFILENAME
```

!!! note "No --help / -h at all"
    `hicCreateThresholdFile` reproduces the Python tool's argument parser, which is built with
    `add_help=False` and never adds a help option back: neither `-h` nor `--help` is recognized, and
    either one exits with an "unrecognised argument" style error (after first reporting any missing
    required arguments, since the reference parser validates what it parsed before it reports
    leftovers). There is consequently no `--help` output to generate this reference from; the table
    below is written directly from the C++ source (`cpp/tools/hicCreateThresholdFile.cpp`) and the
    Python tool it ports.

Writes a two-line header and one row per bin over the requested range. The whole tool is about 15 lines
of Python and there is nothing numeric to it beyond a threshold value repeated per bin, but two details
matter for a faithful port:

- the loop runs to `range[1] + resolution`, exclusive, so the last row is at `range[1]` and is included;
  stopping at `range[1]` itself would drop it.
- the threshold is written with Python's `'{}'.format(float)`, which is `repr(float)`: the shortest
  string that round-trips. `-tv 1` prints as `1.0` and `-tv 0.00001` as `1e-05`, neither of which a
  plain `%g` printf format gives; the C++ port reproduces this with its own float-repr routine.

## Required arguments

| Flag | Meaning |
|---|---|
| `--thresholdValue THRESHOLDVALUE, -tv THRESHOLDVALUE` | Standard threshold value for all relative distances. |
| `--range RANGE RANGE` | Defines the region upstream and downstream of a reference point that should be included. |
| `--outFileName OUTFILENAME, -o OUTFILENAME` | The name and path of the created threshold file. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--resolution RESOLUTION, -r RESOLUTION` | Resolution of the bin in genomic units (default: 1000). |
