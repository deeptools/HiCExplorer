# hicBuildMatrix

Creates a Hi-C matrix from the aligned BAM/SAM files (or a .pairs file) of the Hi-C sequencing reads.

```text
usage: hicBuildMatrix (--samFiles two sam files two sam files | --pairsFile pairs file)
                      --outFileName FILENAME --QCfolder FOLDER
                      [--restrictionCutFile BED file [BED file ...]]
                      [--restrictionSequence RESTRICTIONSEQUENCE [RESTRICTIONSEQUENCE ...]]
                      [--danglingSequence DANGLINGSEQUENCE [DANGLINGSEQUENCE ...]]
                      [--outBam bam file]
                      [--binSize BINSIZE [BINSIZE ...]] [--minDistance MINDISTANCE]
                      [--maxDistance MAXDISTANCE]
                      [--maxLibraryInsertSize MAXLIBRARYINSERTSIZE]
                      [--genomeAssembly GENOMEASSEMBLY] [--region CHR:START-END]
                      [--keepSelfLigation] [--keepSelfCircles]
                      [--minMappingQuality MINMAPPINGQUALITY] [--threads THREADS]
                      [--inputBufferSize INPUTBUFFERSIZE] [--doTestRun]
                      [--doTestRunLines DOTESTRUNLINES] [--skipDuplicationCheck]
                      [--chromosomeSizes txt file] [--noPlot] [--help] [--version]
```

Using an alignment from a program that supports local alignment (eg. Bowtie2) where both PE reads are mapped using the --local option, this program reads such file and creates a matrix of interactions.

## Required arguments

| Flag | Meaning |
|---|---|
| `--samFiles two sam files two sam files, -s two sam files two sam files` | The two PE alignment sam files to process |
| `--pairsFile pairs file` | C++ only, instead of --samFiles: a 4DN or pairtools .pairs file, plain, gzip or bgzip compressed. See "Pairs input" below. |
| `--outFileName FILENAME, -o FILENAME` | Output file name for the Hi-C matrix. |
| `--QCfolder FOLDER` | Path of folder to save the quality control data for the matrix. |
| `--restrictionCutFile BED file [BED file ...], -rs BED file [BED file ...]` | BED file(s) with all restriction cut sites. Required with --samFiles. |
| `--restrictionSequence RESTRICTIONSEQUENCE [RESTRICTIONSEQUENCE ...], -seq ...` | Sequence of the restriction site. Required with --samFiles. |
| `--danglingSequence DANGLINGSEQUENCE [DANGLINGSEQUENCE ...]` | Sequence left by the restriction enzyme after cutting. Required with --samFiles. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--outBam bam file, -b bam file` | Output bam file with all valid Hi-C reads. |
| `--binSize BINSIZE [BINSIZE ...], -bs BINSIZE [BINSIZE ...]` | Size in bp for the bins. (default: None) |
| `--minDistance MINDISTANCE` | Minimum distance between restriction sites. (default: 300) |
| `--maxDistance MAXDISTANCE` | Obsolete. Use --maxLibraryInsertSize instead. |
| `--maxLibraryInsertSize MAXLIBRARYINSERTSIZE` | The maximum library insert size. (default: 1000) |
| `--genomeAssembly GENOMEASSEMBLY, -ga GENOMEASSEMBLY` | The genome the reads were mapped to. |
| `--region CHR:START-END, -r CHR:START-END` | Region of the genome to limit the operation to. |
| `--keepSelfLigation` | Keep self ligations. (default: False) |
| `--keepSelfCircles` | Keep self circles. (default: False) |
| `--minMappingQuality MINMAPPINGQUALITY` | Minimum mapping quality. (default: 15) |
| `--threads THREADS` | Number of threads. (default: 4) |
| `--inputBufferSize INPUTBUFFERSIZE` | Size of the input buffer of each thread. (default: 400000) |
| `--doTestRun` | Test only --doTestRunLines reads. (default: False) |
| `--doTestRunLines DOTESTRUNLINES` | Number of lines for the qc test run. (default: 1000000) |
| `--skipDuplicationCheck` | Skip the identification of duplicated read pairs. |
| `--chromosomeSizes txt file, -cs txt file` | File with the chromosome sizes for your genome. |
| `--noPlot` | C++ only: write QC.log and the QC tables without the QC figures and hicQC.html. |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## Pairs input (--pairsFile, C++ only)

Positions are the one-based 5' ends of the reads; a pair is binned at position - 1 into the half-open bin [start, end), as cooler cload pairs bins it. Chromosome sizes come from --chromosomeSizes or the #chromsize header lines. Bins come from --binSize, or without it from

| Flag | Meaning |
|---|---|
| `--restrictionCutFile, --minDistance and --maxLibraryInsertSize.` |  |
| `--minMappingQuality needs the columns mapq1 and mapq2; without them it is` | refused when given and not applied by default. A pair_type letter N, X or W and the chromosome '!' count as unmapped, M as not unique, DD as a duplicate. The duplicate check runs in file order on a file declared '#sorted: chr1-chr2-pos1-pos2' and '#shape: upper triangle', and with a hash set on any other file. Not allowed, since a .pairs file carries no read sequence and no alignment span: --outBam, --restrictionSequence, |
| `--danglingSequence, --keepSelfLigation, --keepSelfCircles. The coverage` | column of an h5 bin table is NaN, and the QC report has the layout of hicBuildMatrixMicroC. |

## Notes

Please note that the file type extension for the output matrix (`--outFileName`) must be given! This can be **.h5**, **.cool** or the specializations of cool 
**.mcool**, if the path is given! For **.scool** files please create one **.cool** file per cell and merge it together with scHiCExplorer's scHicMergeToSCool.

### Building multicooler matrices

`hicBuildMatrix` supports building multicooler matrices which are for example needed for visualization with [HiGlass ](https://higlass.io/).
To do so, use as out file format either .cool or .mcool and define the desired resolutions as `--binSize`.
`hicBuildMatrix` builds the interaction matrix for the highest resolution and merges the bins for the lower resolutions.
The lower resolutions need to be an integer multiplicative of the highest resolution.

```bash
$ hicBuildMatrix -s forward.bam reverse.bam -o multi_resolution.cool 
  --binSize 10000 20000 50000 100000 --QCfolder QC
```
Introducing with version 3.5 we support multiple restriction and dangling end sequences, and multiple restriction cut site files. 
Hi-C protocols that use multiple restriction cut enzymes benefit from this and get now an improved QC report.
Version 3.5 adds also the support for a chromosome size file which can help to get interaction matrices with a predefined size. Capture Hi-C or 
single-cell Hi-C data, where it is not guaranteed that reads from all areas of the chromosome are present benefit from this latest improvement.

### Missing scaffolds or contigs in the Hi-C matrix

Restriction enzymes cut the DNA at their specific restriction cut site sequences. It can occur for scaffolds or contigs (less likely it happens for chromosomes) that it does not contain any of these cut sites.
In this case, the reads from this scaffold or contig are considered as invalid and are part of the 'same restriction fragment'-statistics in the QC report.
