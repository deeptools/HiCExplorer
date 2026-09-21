# chicExportData

Exports the data stored in the capture Hi-C intermediate HDF5 files to text or bigWig.

```text
usage: chicExportData --file FILE [--outFileName OUTFILENAME]
                      [--outputFileType {txt,bigwig}]
                      [--outputMode {all,geneName}]
                      [--outputModeName OUTPUTMODENAME]
                      [--decimalPlaces DECIMALPLACES]
                      [--chromosomeSizes txt file]
                      [--backgroundModelFile BACKGROUNDMODELFILE]
                      [--oneTargetFile] [--range RANGE RANGE]
                      [--outputValueBigwig {relative-interactions,p-value,x-fold,raw}]
                      [--threads THREADS] [--help] [--version]
```

chicExportData exports the data stored in the intermediate hdf5 files to text files per reference point.

## Required arguments

| Flag | Meaning |
|---|---|
| `--file FILE, -f FILE` | path to the file which should be used for data export |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--outFileName OUTFILENAME, -o OUTFILENAME` | Output tar.gz of the files. In case of --outputMode == geneName it is ignored. (Default: data.tar.gz). |
| `--outputFileType {txt,bigwig}, -oft {txt,bigwig}` | Output file type can be set for all file types to txt; except 'interaction' supports also bigwig (Default: txt). |
| `--outputMode {all,geneName}, -om {all,geneName}` | Output mode: Either all date is written or a gene name must be specified. (Default: all). |
| `--outputModeName OUTPUTMODENAME, -omn OUTPUTMODENAME` | ONLY valid if --outputMode geneName! Define the name of the gene |
| `--decimalPlaces DECIMALPLACES` | Decimal places for all output floating numbers in the viewpoint files (Default: 12). |
| `--chromosomeSizes txt file, -cs txt file` | File with the chromosome sizes for your genome. A tab- delimited two column layout "chr_name size" is expectedUsually the sizes can be determined from the SAM/BAM input files, however, for cHi-C or scHi-C it can be that at the start or end no data is present. Please consider that this option causes that only reads are considered which are on the listed chromosomes.Use this option to guarantee fixed sizes. An example file is available via UCSC: http://hgdownlo ad.soe.ucsc.edu/goldenPath/dm3/bigZips/dm3.chrom.sizes |
| `--backgroundModelFile BACKGROUNDMODELFILE, -bmf BACKGROUNDMODELFILE` | Path to the background model file. Required only for fileType=interactions and outputFileTypeBigwig. |
| `--oneTargetFile, -otf` | Compile all target files to one. Applies only if --fileType is target |
| `--range RANGE RANGE` | Defines the region upstream and downstream of a reference point which should be included. Format is --range upstream downstream, e.g.: --range 500000 500000 plots 500kb up- and 500kb downstream. This value should not exceed the range used in the other chic-tools. Applies only for interaction files in the combination with bigwig and a background model file! |
| `--outputValueBigwig {relative-interactions,p-value,x-fold,raw}, -ovb {relative-interactions,p-value,x-fold,raw}` | Select which value the bigwig file should contain: 'relative-interactions', 'p-value', 'x-fold', 'raw' (Default: relative-interactions). |
| `--threads THREADS, -t THREADS` | Number of threads (uses the python multiprocessing module) (Default: 4). |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |
