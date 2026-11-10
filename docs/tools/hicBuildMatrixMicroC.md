# hicBuildMatrixMicroC

Creates a Hi-C matrix from Micro-C data (hicBuildMatrix without the restriction-site machinery).

```text
usage: hicBuildMatrixMicroC --samFiles two sam files two sam files
                            --outFileName FILENAME --QCfolder FOLDER
                            [--outBam bam file] --binSize BINSIZE
                            [BINSIZE ...]
                            [--maxLibraryInsertSize MAXLIBRARYINSERTSIZE]
                            [--genomeAssembly GENOMEASSEMBLY]
                            [--region CHR:START-END] [--keepSelfCircles]
                            [--minMappingQuality MINMAPPINGQUALITY]
                            [--threads THREADS]
                            [--inputBufferSize INPUTBUFFERSIZE] [--doTestRun]
                            [--doTestRunLines DOTESTRUNLINES]
                            [--skipDuplicationCheck]
                            [--chromosomeSizes txt file] [--help] [--version]
```

Using an alignment from a program that supports local alignment (eg. Bowtie2) where both PE reads are mapped using the --local option, this program reads such file and creates a matrix of interactions.

## Required arguments

| Flag | Meaning |
|---|---|
| `--samFiles two sam files two sam files, -s two sam files two sam files` | The two PE alignment sam files to process |
| `--outFileName FILENAME, -o FILENAME` | Output file name for the Hi-C matrix. |
| `--QCfolder FOLDER` | Path of folder to save the quality control data for the matrix. |
| `--binSize BINSIZE [BINSIZE ...], -bs BINSIZE [BINSIZE ...]` | Size in bp for the bins. (default: None) |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--outBam bam file, -b bam file` | Output bam file with all valid Hi-C reads. |
| `--maxLibraryInsertSize MAXLIBRARYINSERTSIZE` | The maximum library insert size. (default: 1000) |
| `--genomeAssembly GENOMEASSEMBLY, -ga GENOMEASSEMBLY` | The genome the reads were mapped to. |
| `--region CHR:START-END, -r CHR:START-END` | Region of the genome to limit the operation to. |
| `--keepSelfCircles` | Keep self circles. (default: False) |
| `--minMappingQuality MINMAPPINGQUALITY` | Minimum mapping quality. (default: 15) |
| `--threads THREADS` | Number of threads. (default: 4) |
| `--inputBufferSize INPUTBUFFERSIZE` | Size of the input buffer of each thread. (default: 400000) |
| `--doTestRun` | Test only --doTestRunLines reads. (default: False) |
| `--doTestRunLines DOTESTRUNLINES` | Number of lines for the qc test run. (default: 1000000) |
| `--skipDuplicationCheck` | Skip the identification of duplicated read pairs. |
| `--chromosomeSizes txt file, -cs txt file` | File with the chromosome sizes for your genome. |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |
