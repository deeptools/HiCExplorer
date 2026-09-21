# Usage

## General conventions

A typical HiCExplorer command looks like this:

```bash
hicPlotMatrix -m myHiCmatrix.cool -o myHiCmatrix.pdf \
  --clearMaskedBins --region chrX:10,000,000-15,000,000 --vMin -4 --vMax 4
```

Every tool prints its full option list with `--help` (a few tools with only required arguments, such as
`hicCompartmentalization` or `hicCreateThresholdFile`, need those required arguments present before they
will honor `--help`; this is inherited from the Python argparse behaviour they reproduce). Each tool's
own page under [Tools](tools/index.md) lists the real, current flag set of the built C++ binary.

- The output format of a plot is set by the file name's extension, for example `myPlot.pdf` for a PDF
  file or `myPlot.png` for a PNG file.
- Most tools that produce a plot can also write the underlying numeric data instead of, or in addition
  to, the figure (`--plotData FILE` on the plotting tools that delegate drawing to the Python layer, see
  below), so the data can be re-plotted with another tool if the built-in visualization is not what you
  need.
- Native matrix files are `.cool`, `.mcool` or the legacy `.h5` format; see
  [File formats](file-formats.md).

## The Hi-C data pipeline

Building and analyzing a Hi-C contact matrix generally follows four steps: map the reads, build and
filter the matrix, correct it, then analyze or visualize the corrected matrix. See
[Example usage](example-usage.md) for a full worked walkthrough with real commands at each step, and
[Tools](tools/index.md) for the complete tool list grouped by task (pre-processing, quality control,
analysis, TAD calling, visualization, matrix format handling, and capture Hi-C).

## The C++ plotting tools and the Python drawing layer

The C++ rewrite computes every plot's underlying data (matrices, positions, statistics) natively. A
subset of the plotting tools (`hicPlotMatrix`, `hicPlotTADs`, `hicAggregateContacts`,
`hicPlotAverageRegions`, `hicCompartmentalization`, and the `diagnostic_plot` subcommand of
`hicCorrectMatrix`) then hand that data to the original Python/matplotlib (and, for `hicPlotTADs`,
pyGenomeTracks) drawing code through a small layer named `hicexplorer_plot`, rather than reimplementing
matplotlib's rendering in C++. The Python interpreter to use for this is set with the `HICX_PLOT_PYTHON`
environment variable. Passing `--plotData FILE` on these tools skips the drawing step and writes the
figure's data as JSON (plus `.npy` matrices where relevant) instead, which is useful in a pipeline that
does not want a Python dependency at all, or that wants to re-plot with different styling.

If `HICX_PLOT_PYTHON` is not set or the interpreter it names cannot draw (missing matplotlib, for
example), these tools refuse to draw and exit with a non-zero status rather than silently skip the
figure; their numeric outputs, if requested with `--plotData` or the tool's own data-output options, are
still written.

## Threads and determinism

Tools that accept a `--threads` / `--numberOfProcessors` option are deterministic regardless of the
thread count: the number of threads changes only wall-clock time and peak memory, never the numeric
result. This was checked as part of each tool's equivalence validation (see `cpp/STATUS.md` in the
repository for the per-tool ledger).
