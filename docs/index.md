# HiCExplorer

**HiCExplorer** is a suite of command-line tools to process, normalize, analyze and visualize Hi-C,
Micro-C and capture Hi-C (cHi-C) data. This is version 4.

## What is new in version 4

- **Faster, much less memory.** On a real single-chromosome Hi-C matrix (24,926 bins, 61.8M nonzero
  pixels), the tools of version 4 measured 5 to 70 times faster than HiCExplorer 3.7.6,
  single-threaded, at 3 to 9 times less peak memory. Most tools also run multithreaded. See
  [Benchmarks](benchmarks.md) for the numbers.
- **The same results.** Every tool is checked against the output of HiCExplorer 3.7.6 on real and
  synthetic data. Command lines, file formats and outputs stay compatible, and known differences are
  listed on the page of the tool.
- **CHiCAGO scoring.** The `chicChicago*` tools score capture Hi-C interactions with the CHiCAGO
  method from cool, h5 or `.hic` matrices, and plot them as viewpoints, arcs and genome tracks. See the
  [CHiCAGO tutorial](example-usage/chicago-tutorial.md).
- **New analysis tools.** `hicDifferentialAnalysis` (replicate-aware differential TADs, loops and
  compartments) and `hicDetectStripes`.
- **One package.** `conda install hicexplorer` installs the tools together with the Python API. Existing
  scripts that call the commands or import `hicexplorer` keep working.

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

- [Installation](installation.md): installing with conda, pip or from source.
- [Usage](usage.md): the general command-line conventions shared by all tools.
- [Tools](tools/index.md): the full list of tools, grouped by task, each with its own CLI reference.
- [File formats](file-formats.md): the native h5/cool matrix format and the capture Hi-C HDF containers.
- [Example usage](example-usage.md): a worked pre-processing-to-visualization walkthrough.
- [Benchmarks](benchmarks.md): runtime and memory compared with HiCExplorer 3.7.6.
- [News](news.md): the release history.
- [Citation](citation.md): how to cite HiCExplorer.
- [Support](support.md): where to get help.
