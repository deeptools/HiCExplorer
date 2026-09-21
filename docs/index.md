# HiCExplorer

**HiCExplorer** is a suite of command-line tools to process, normalize, analyze and visualize Hi-C,
Micro-C and capture Hi-C (cHi-C) data. This is a from-scratch C++ rewrite of the original Python
HiCExplorer, validated tool by tool against the Python reference implementation (see
[Benchmarks](benchmarks.md) for what that validation covers).

## Why the C++ rewrite

- **Same results.** Every ported tool is checked against the Python HiCExplorer's own output on real
  and synthetic data before it is considered done; equivalence classes and known, documented CLI
  deviations are tracked per tool.
- **Faster, much less memory.** On a real single-chromosome Hi-C matrix (24,926 bins, 61.8M nonzero
  pixels), the C++ tools measured 5 to 70 times faster than Python HiCExplorer 3.7.6, single-threaded,
  at 3 to 9 times less peak memory. See [Benchmarks](benchmarks.md) for the numbers.
- **No Python runtime dependency for the core pipeline.** The tools are native binaries built against
  htslib, HDF5 and a small set of C++ libraries; no conda environment with SciPy/pandas/matplotlib is
  required to build a matrix, correct it, call loops or TADs, or convert formats. A handful of the
  plotting tools still delegate figure drawing to the original Python/matplotlib code through a small
  drawing layer (`HICX_PLOT_PYTHON`), since reproducing matplotlib's rendering pixel-for-pixel in C++
  was judged not worth it; the data these plots draw from is computed in C++.

## At a glance

```bash
# build a Hi-C matrix from two independently mapped SAM/BAM files
hicBuildMatrix --samFiles mate_R1.bam mate_R2.bam \
  --binSize 10000 --restrictionSequence GATC --danglingSequence GATC \
  --restrictionCutFile cut_sites.bed --outBam hic.bam \
  -o hic_matrix.cool --QCfolder ./hicQC

# correct the matrix (Knight-Ruiz, the default)
hicCorrectMatrix correct -m hic_matrix.cool -o hic_corrected.cool

# call TADs
hicFindTADs -m hic_corrected.cool --outPrefix hic_corrected
```

## Where to go next

- [Installation](installation.md): building from source and what changed from the Python package.
- [Usage](usage.md): the general command-line conventions shared by all tools.
- [Tools](tools/index.md): the full list of tools, grouped by task, each with its own CLI reference.
- [File formats](file-formats.md): the native h5/cool matrix format and the capture Hi-C HDF containers.
- [Example usage](example-usage.md): a worked pre-processing-to-visualization walkthrough.
- [Benchmarks](benchmarks.md): runtime and memory versus the Python implementation.
- [News](news.md): the Python HiCExplorer release history this rewrite builds on.
- [Citation](citation.md): how to cite HiCExplorer.
- [Support](support.md): where to get help.
