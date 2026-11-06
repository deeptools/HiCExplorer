# HiCExplorer

**HiCExplorer** is a suite of command-line tools to process, normalize, analyze and visualize Hi-C,
Micro-C and capture Hi-C (cHi-C) data. This is version 4.

## What is new in version 4

- **Faster, with less memory.** On real data, version 4 computes the results of HiCExplorer 3.7.6
  two to ten times faster with two to four times less peak memory, and most tools run
  multithreaded:

    | Task | 3.7.6 | 4.0 | Speedup |
    |---|---|---|---|
    | `hicBuildMatrix`, 92.7M reads, 10 kb | 42 min 18 s | 16 min 22 s | 2.6x |
    | `hicDetectLoops`, GM12878 chr1, 10 kb | 20.6 s | 2.1 s | 10.0x |
    | `hicFindTADs`, GM12878 chr1, 10 kb | 33.0 s | 6.3 s | 5.2x |
    | `hicPCA`, GM12878 chr1, 100 kb | 24.2 s | 3.8 s | 6.4x |

    CHiCAGO scoring of 22,076 baits takes 48.5 s instead of 437.5 s in R CHiCAGO, and stripe
    calling 210 s instead of 3,015 s in Stripenn. See [Benchmarks](benchmarks.md) for memory,
    threads and the data.

- **The same results.** Every tool is checked against the output of HiCExplorer 3.7.6 on real and
  synthetic data. Command lines, file formats and outputs stay compatible, and known differences are
  listed on the page of the tool.
- **CHiCAGO scoring.** The `chicChicago*` tools score capture Hi-C interactions with the CHiCAGO
  method from cool, h5 or `.hic` matrices, and plot them as viewpoints, arcs and genome tracks. See the
  [CHiCAGO tutorial](example-usage/chicago-tutorial.md).
- **New analysis tools.** `hicDifferentialAnalysis` (differential TADs, loops and compartments
  between conditions with replicates), `hicDetectStripes` (architectural stripes with the method of
  Stripenn), `hicArrowhead` (contact domains with the Arrowhead algorithm), nested domains and a CTCF
  motif filter in `hicFindTADs`, and the five reproducibility measures SCC, GenomeDISCO, HiC-Spector,
  QuASAR-Rep and ENT3C in `hicCorrelate`.
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
