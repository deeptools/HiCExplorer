# hicConvertFormat

Converts between Hi-C interaction matrix formats (hic, cool, mcool, h5, homer, HiC-Pro, ginteractions).

```text
usage: hicConvertFormat --matrices MATRICES [MATRICES ...] --outFileName
                        OUTFILENAME [OUTFILENAME ...] --inputFormat
                        {h5,cool,hic,homer,hicpro,2D-text,chinput} --outputFormat
                        {cool,h5,homer,ginteractions,mcool,hicpro,hic}
                        [--correction_name CORRECTION_NAME]
                        [--correction_division] [--store_applied_correction]
                        [--chromosome CHROMOSOME] [--enforce_integer]
                        [--load_raw_values] [--resolutions RESOLUTIONS [RESOLUTIONS ...]]
                        [--help] [--chromosomeSizes txt file] [--version]
                        [--bedFileHicpro BEDFILEHICPRO [BEDFILEHICPRO ...]]
                        [--hicVersion {8,9}]
                        [--hicNormalizations {VC,VC_SQRT,KR,SCALE,none} [{VC,VC_SQRT,KR,SCALE,none} ...]]
                        [--threads THREADS]
```

Conversion of Hi-C matrices of different file formats. We support the conversion of hic to cool format via hic2cool, and homer, HicPro, h5 and cool format to h5, cool, homer or ginteractions format. Moreover, hicConvertFormat accepts multiple input files from one format with different resolutions and creates a mcool file. Each original file is stored under the path e.g. ::/resolutions/10000. A batch computation is possible, the number of input files and output files needs to match, all input files need to be of the same format type and all output files too. For input and output of cooler files special options are available, for all other formats they will be ignored. HiCPro file format needs an additional bed file as input.

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]` | input file(s). Could be one or many files. |
| `--outFileName OUTFILENAME [OUTFILENAME ...], -o OUTFILENAME [OUTFILENAME ...]` | File name to save the exported matrix. |
| `--inputFormat {h5,cool,hic,homer,hicpro,2D-text,chinput}` | File format of the input matrix file. `chinput` is C++ only and needs `--rmap`. |
| `--outputFormat {cool,h5,homer,ginteractions,mcool,hicpro,hic}` | Output format. (Default: cool). |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--correction_name CORRECTION_NAME` | Name of the column which stores the correction factors. Option only for cool input files. (Default: weight). |
| `--correction_division` | If set, division is applied for correction. |
| `--store_applied_correction` | Store the applied correction and do not set correction factors. |
| `--chromosome CHROMOSOME` | Load only one chromosome. |
| `--enforce_integer` | Enforce datatype of counts to integer. |
| `--load_raw_values` | Load only 'count' data and do not apply a correction. |
| `--resolutions RESOLUTIONS [RESOLUTIONS ...], -r RESOLUTIONS [RESOLUTIONS ...]` | List of resolutions that should be added. |
| `--help, -h` | show this help message and exit. |
| `--chromosomeSizes txt file, -cs txt file` | For the input format `2D-text` only. |
| `--version` | show program's version number and exit |
| `--bedFileHicpro BEDFILEHICPRO [BEDFILEHICPRO ...], -bf BEDFILEHICPRO [BEDFILEHICPRO ...]` | Bed file(s) of hicpro file format. |
| `--hicVersion {8,9}` | Version of a .hic output file. (Default: 8). |
| `--hicNormalizations {VC,VC_SQRT,KR,SCALE,none} [{VC,VC_SQRT,KR,SCALE,none} ...]` | Normalizations a .hic output file stores, computed as Juicer tools addNorm does. (Default: VC VC_SQRT KR SCALE). |
| `--rmap RMAP [RMAP ...]` | C++-only. CHiCAGO `.rmap` file(s) that define the bins of a `chinput` input, one per input. |
| `--threads THREADS` | Threads for compressing a .hic output file; the file does not depend on the number. (Default: 1). |

## Notes

#### Background

To reproduce analyses and to compare and use different Hi-C analysis software, conversion between interaction matrix file formats is crucial.
However, most Hi-C softwares are only supporting their own data format which makes the exchange difficult. HiCExplorer supports a range of 
interaction matrices, both for import and for export. 

Import:
    - chinput (C++ only, needs `--rmap`; bins are the restriction fragments and both directions of a bait-to-bait pair are summed)
    - hic
    - cool
    - h5
    - homer
    - HicPro

Export:
    - cool
    - mcool
    - h5
    - homer
    - ginteractions

With HiCExplorer version 2.2 hicConvertFormat and hicAdjustMatrix replace hicExport from HiCExplorer 2.1 and older versions.

#### Usage example

##### hic2cool

HiCExplorer uses the library [hic2cool](https://github.com/4dn-dcic/hic2cool)  to convert **.hic** interaction matrix files to the cool format. Usually .hic files 
have the three correction factors **KR**, **VC** or **VC_SQRT**; however these cannot be applied natively by HiCExplorer tools because 
HiCExplorer expects the correction values to be stored in the column **weight**.
To work with corrected data the correction factors need to applied separately, see section cool to cool.

The following example will convert a hic file which contains the resolution of 1000 to a cool file with 10kb resolution. The desired 
resolution needs to be existing in the hic file. If no resolution parameter is defined, a mcool file with all available resolutions is created.

```bash
$ hicConvertFormat -m matrix.hic --inputFormat hic --outputFormat cool -o matrix.cool --resolutions 10000
```
It is only possible to convert from hic to cool format, no other formats are supported.

##### cool to cool

The cool file format is developed and maintained by the [Mirny](https://github.com/mirnylab/cooler) lab and allows to access interaction matrices in an easy to use data format.

The cool data format allows to use the following options:

- correction_name: In case correction factors are not stored in 'weight' the correct column name can be defined using this parameter and the resulting matrix will store the values in 'weight'.
- correction_division: Correction factors can be applied by a multiplication or a division. The default behaviour is to use the multiplication, in case the correction factors are inverted, set this parameter.
- store_applied_correction: Set this parameter if correction factors should be applied on the data and should be written back to colum 'counts' in the corrected form and not as raw. Default: not set.
- chromosomes: Define a list of chromosomes which should be included in the output matrix. All chromosomes which are not defined are not part of the new matrix. This parameter can speed up the processing especially if only one chromosome is used.
- enforce_integer: Raw interaction data is stored as integers, after the correction is applied the data is a float. Set a this parameter to enforce integer values in the new matrix.
- load_raw_values: Set this parameter if the interaction data should not be loaded with the correction factors.

### Example usage

```bash
$ hicConvertFormat -m matrix.cool --inputFormat cool --outputFormat cool -o matrix.cool --correction_name KR
```
##### Homer

[Homer](http://homer.ucsd.edu/homer/index.html) is a software for 'motif discovery and next generation sequencing analysis' and supports Hi-C. HiCExplorer is able to read and write the interaction matrix from Homer. Homer stores the interaction matrix in a simple text file as a dense matrix. To write 
large matrices in Homer format needs a lot of space and can take a few ours to days. 

### Example usage

```bash
$ hicConvertFormat -m matrix.homer --inputFormat homer --outputFormat cool -o matrix.cool
```
##### Hic-Pro

[HiC-Pro](https://github.com/nservant/HiC-Pro) file format needs an additional bed file as input:

### Example usage

```bash
$ hicConvertFormat -m matrix.hicpro --bedFileHicpro hicpro.bed --inputFormat hicpro --outputFormat cool -o matrix.cool
```
##### Create a mcool file

With HiCExplorer it is possible to create a multiple cool (mcool) file. These mcool files can be used e.g. with [HiGlass](http://higlass.io/).

To create an mcool file, use as input either one matrix in one of the supported read formats and define the desired resolutions or define
multiple input matrices. In the second case, the matrices should all have different resolutions.

### Example usage

The resolutions need to be a multiple of the input matrix i.e. matrix with 10kb, 20kb and 30kb are possible but not 35kb.

```bash
$ hicConvertFormat -m matrix.cool --inputFormat cool --outputFormat mcool
   -o multi_matrix.mcool --resolutions 20000 40000 70000 120000 500000
```
```bash
$ hicConvertFormat -m matrix10kb.cool matrix20kb.cool matrix30kb.cool 
    --inputFormat cool --outputFormat mcool -o multi_matrix.mcool
```
The mcool matrix contains the individual matrices as follows:

```
multi_matrix.mcool::/resolutions/10000
multi_matrix.mcool::/resolutions/40000
multi_matrix.mcool::/resolutions/70000
multi_matrix.mcool::/resolutions/120000
multi_matrix.mcool::/resolutions/500000
```
