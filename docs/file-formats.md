# File formats

HiCExplorer has a native interaction matrix format, h5 (an HDF5 container; `.cool`/`.mcool` are also
supported natively and are the recommended format, see below), and several native HDF5 containers for
capture Hi-C analysis.

!!! note "Python code below is the format's reference implementation"
    The `h5` layout was originally defined and is still read/written by the Python `HiCMatrix` library;
    the code samples in this section are that reference implementation, kept as the most precise
    description of the on-disk layout. The tools read and write the same h5 and cool layouts natively
    (see `cpp/core` in the repository), with no dependency on this Python code at run time; the samples are useful if you need to read or write these files from your own
    Python scripts, or as a byte-for-byte specification of the format.

HiCExplorer has a native interaction file format, h5, and several native formats for capture Hi-C analysis.

### Interaction matrix format h5

The h5 format is implemented as a hdf5 container and can be accessed either via HiCMatrix, the hdf5 API or a graphical user interface like hdfview.
The data in the hdf5 format is read and written via the PyTables interface.

The format is structured as follows:

```bash
- intervals
   - chr_list (CArray)
   - end_list (CArray)
   - extra_list (CArray)
   - start_list (CArray)
- matrix
   - data (CArray)
   - indices (CArray)
   - indptr (CArray)
   - shape (tupel)
 - distance counts (CArray) [optional]
 - nan_bins (CArray) [optional]
 - correction_factors (CArray) [optional]
```
The group 'matrix' contains the data which is necessary to create a scipy sparse csr matrix: 

```python
import tables
with tables.open_file(self.matrixFileName, 'r') as f:
        parts = {}
        try:
            for matrix_part in ('data', 'indices', 'indptr', 'shape'):
                parts[matrix_part] = getattr(f.root.matrix, matrix_part).read()
        except Exception as e:
            log.info('No h5 file. Please check parameters concerning the file type!')
            e
        matrix = csr_matrix(tuple([parts['data'], parts['indices'], parts['indptr']]),
                            shape=parts['shape'])
```
The group 'intervals' contains the information of the genomic regions associated to a bin in the matrix via its index. To retrieve the data run:

```python
with tables.open_file(self.matrixFileName, 'r') as f:
    intvals = {}
    for interval_part in ('chr_list', 'start_list', 'end_list', 'extra_list'):
        if toString(interval_part) == toString('chr_list'):
            chrom_list = getattr(f.root.intervals, interval_part).read()
            intvals[interval_part] = toString(chrom_list)
        else:
            intvals[interval_part] = getattr(f.root.intervals, interval_part).read()

    cut_intervals = list(zip(intvals['chr_list'], intvals['start_list'], intvals['end_list'], intvals['extra_list']))
```
To write the format, the above structure needs to be created via PyTables. For example:

```python
with tables.open_file(filename, mode="w", title="HiCExplorer matrix") as h5file:
    matrix_group = h5file.create_group("/", "matrix", )
    for matrix_part in ('data', 'indices', 'indptr', 'shape'):
            arr = np.array(getattr(matrix, matrix_part))
            atom = tables.Atom.from_dtype(arr.dtype)
            ds = h5file.create_carray(matrix_group, matrix_part, atom,
                                      shape=arr.shape,
                                      filters=filters)
```
The matrix object is a scipy `csr_matrix`. 

Please refer to HiCMatrix for a reference implementation: https://github.com/deeptools/HiCMatrix/blob/master/hicmatrix/lib/h5.py

To open a cool / h5 file with the HiCMatrix library, use the following code:

```python
import hicmatrix.HiCMatrix as hm

hic_ma = hm.hiCMatrix("/path/to/matrix")

# csr_matrix
hic_ma.matrix

# the corresponding genomic regions as a list
hic_ma.cut_intervals
```
To write a h5 use the following code. Please consider, once the datatype of a matrix is specified, for example by reading a h5 matrix or by first time writing a h5 matrix, it cannot be changed to a cool file anymore (except by manipulating the internal matrixFilehandler object). 

```python
import hicmatrix.HiCMatrix as hm
from scipy.sparse import csr_matrix

# create a HiCMatrix object
hic_ma = hm.hiCMatrix()

# set important data structures

hic_ma.matrix = csr_matrix()
hic_ma.cut_intervals = list[(chr1, start1, end1, 1), (chr1, start2, end2, 1), ..., (chr1, startN, endN, 1)]

# to store a h5 matrix
hic_ma.save('/path/to/storage/matrix.h5')
# to store a cool matrix
hic_ma.save('/path/to/storage/matrix.cool')
```
Alternatively, the matrixFileHandler object can be accessed directly:

```python
from hicmatrix.lib import MatrixFileHandler

# Load a matrix via the MatrixFileHandler class

matrixFileHandlerInput = MatrixFileHandler(pFileType=args.inputFormat, pMatrixFile=matrix,
                                                        pCorrectionFactorTable=args.correction_name,
                                                        pCorrectionOperator=correction_operator,
                                                        pChrnameList=chromosomes_to_load,
                                                        pEnforceInteger=args.enforce_integer,
                                                        pApplyCorrectionCoolerLoad=applyCorrectionCoolerLoad)

# Load data
_matrix, cut_intervals, nan_bins, \
   distance_counts, correction_factors = matrixFileHandlerInput.load()

# create matrixFileHandler output object
matrixFileHandlerOutput = MatrixFileHandler(pFileType=args.outputFormat, pEnforceInteger=args.enforce_integer, pFileWasH5=format_was_h5, pHic2CoolVersion=hic2CoolVersion, pHiCInfo=cool_metadata)

# set the variables
matrixFileHandlerOutput.set_matrix_variables(_matrix, cut_intervals, nan_bins,
                                                          correction_factors, distance_counts)

matrixFileHandlerOutput.save(args.outFileName, pSymmetric=True, pApplyCorrection=applyCorrection)
```
### CHiCAGO formats

The `chicChicago*` tools read the input files of R Chicago (Cairns et al. 2016). All files are tab
separated text. The examples are the first lines of the fixture in
`hicexplorer/test/test_data/chicago/` (GM12878, chromosomes 20 and 21, HindIII, hg19).

#### Restriction fragment map (`.rmap`)

One line per restriction fragment of the genome, sorted by position: chromosome, start, end and the
fragment ID. Coordinates are 1-based and inclusive, and the fragments of a chromosome tile it without gaps.
The IDs are unique integers, in order of position.

```text
20	1	60024	403460
20	60025	62654	403461
20	62655	66475	403462
20	66476	71941	403463
```

#### Bait map (`.baitmap`)

One line per baited fragment: chromosome, start, end, fragment ID and a bait name (for example the gene
symbol of the promoter). The coordinates and IDs are those of the `.rmap`, and every bait is a fragment of
the `.rmap`.

```text
20	66476	71941	403463	DEFB125
20	119103	138049	403475	DEFB126
20	161620	170741	403482	DEFB128
```

#### Design tables (`.npb`, `.nbpb`, `.poe`)

The three tables are computed from the `.rmap` and the `.baitmap` with `makeDesignFiles.py` of
[chicagoTools](https://github.com/RegulatoryGenomicsGroup/chicago). They depend only on the design and on the
settings of the run (minimum and maximum fragment length, `maxLBrownEst`, the distance bin size, and
whether bait-to-bait pairs and adjacent fragments are removed), so one set serves all samples of a design.
The first line of each table, starting with `#`, records the settings and the two input files. The
tools use the tables of the same settings as the run.

The distance from a bait to the other end is split into bins of `binsize` (20,000 bp by default) up to
`maxLBrownEst` (1,500,000 bp by default), which gives 75 bins.

`.npb` (NPerBin) has one line per bait. The first column is the bait ID and the following 75 columns give
the number of other-end fragments that lie in each distance bin around the bait and pass the fragment length
filters.

```text
#	minFragLen=150	maxFragLen=40000	maxLBrownEst=1500000	binsize=20000	removeb2b=True	removeAdjacent=True	rmapfile=./h19_chr20and21.rmap	baitmapfile=./h19_chr20and21.baitmap
403463	4	4	3	4	1	4	4	6	6	1	1	3	...
```

`.nbpb` (NBaitsPerBin) has one line per restriction fragment. The first column is the fragment ID and the
following 75 columns give the number of baits for which the fragment lies in that distance bin.

```text
#	maxLBrownEst=1500000	binsize=20000	rmapfile=./h19_chr20and21.rmap	baitmapfile=./h19_chr20and21.baitmap
403460	0	1	0	0	1	0	1	0	1	0	1	0	1	1	...
403461	1	0	0	1	0	1	0	1	1	0	1	0	1	1	...
```

`.poe` (ProxOE) lists the proximal pairs: one line per bait and other end within `maxLBrownEst` after the
filters, with bait ID, other-end ID and the distance in base pairs.

```text
#	minFragLen=150	maxFragLen=40000	maxLBrownEst=1500000	binsize=20000	removeb2b=True	removeAdjacent=True	rmapfile=./h19_chr20and21.rmap	baitmapfile=./h19_chr20and21.baitmap
403463	403461	7869
403463	403465	5827
403463	403466	11771
```

#### Interaction counts (`.chinput`)

One sample per file. The first line starts with `##` and records the sample name and the design files.
The second line names the columns. Each following line is one pair of a bait and an other end: bait ID,
other-end ID, the read count N, the length of the other-end fragment and the signed distance
(negative when the other end lies upstream of the bait). Pairs on different chromosomes (trans) have the
literal `NA` as distance. Several replicates are summed per pair beforehand.

```text
##	samplename=GM_rep1	bamname=...	baitmapfile=...	digestfile=...
baitID	otherEndID	N	otherEndLen	distSign
403463	403461	2	2629	-7869
403463	403462	1	3820	-4644
```

`bam2chicago.sh` of chicagoTools writes the file from the aligned reads. `hicConvertFormat --inputFormat
chinput --rmap FILE.rmap` writes it as a matrix; both directions of a reciprocal bait-to-bait pair are
summed in the matrix.

#### Output of the tools

`chicChicagoBackgroundModel` writes the background model and `chicChicagoScores` the scores, both as plain
text tables. `chicChicagoScores` has one line per interaction with the columns `baitID`, `otherEndID`, `N`,
`distSign`, `Bmean` (expected count of the Brownian component), `Tmean` (expected count of the technical
noise), `log_p` (natural logarithm of the p-value) and `score`. `chicChicagoSignificantInteractions` writes the
interactions above the score threshold with the coordinates: `bait_chr`, `bait_start`, `bait_end`,
`bait_name`, `otherEnd_chr`, `otherEnd_start`, `otherEnd_end`, `N`, `distSign` and `score`. See the
[CHiCAGO tutorial](example-usage/chicago-tutorial.md).

### Capture Hi-C HDF containers

The capture Hi-C data analysis creates for the scripts `chicViewpoint`, `chicSignificantInteractions`, `chicAggregateStatistic` and `chicDifferentialTest` individual HDF containers to store the processed data.

##### chicViewpoint

```bash
- averageContactBin (int)
- fixateRange (int)
- range (int, int)
- resolution (int)
- type='interactions' (string)
- matrix 1 
 - chromosome 1
    - gene name 1
        - chromosome (String)
        - end_list (array)
        - gene (String)
        - interaction_data_list (array)
        - pvalue (array)
        - raw (array)
        - reference_point_end (int)
        - reference_point_start (int)
        - relative_position_list (array)
        - start_list (array)
        - sum_of_interactions (float)
        - xfold (array)
    - gene name 2
        - ...
        - ...
     - ...
 - chromosome 2
    - gene 1
       - ...
       - ...
    - ...
 - ...
 - ...
 - genes
   - gene name 1 (link to matrix 1 / chromosome 1 / gene name 1
   - ...
- matrix 1
  - chromosome 1
     - gene name 1
         - chromosome (String)
         - end_list (array)
         - gene (String)
         - interaction_data_list (array)
         - pvalue (array)
         - raw (array)
         - reference_point_end (int)
         - reference_point_start (int)
         - relative_position_list (array)
         - start_list (array)
         - sum_of_interactions (float)
         - xfold (array)
     - gene name 2
         - ...
         - ...
     - ...
  - chromosome 2
     - gene 1
        - ...
        - ...
  - ...
  - ...

 - genes
   - gene name 1 (link to matrix 2 / chromosome 1 / gene name 1
   - ...
```
##### chicSignificantInteractions

`chicSignificantInteractions` creates two files: a target file and a file containing the significant interactions:

Depending on the combination mode (single / dual) the structure is slightly different for the target file:

```bash
- combinationMode = 'single' (String)
- fixateRange (int)
- mode_preselection (String)
- mode_preselection_value (float)
- peakInteractionsThreshold (float)
- pvalue (float)
- range (int, int)
- truncateZeroPvalues (Boolean)
- type='target' (String)
- matrix 1 
 - chromosome 1
    - gene name 1
        - chromosome (String)
        - end_list (array)
        - reference_point_end (int)
        - reference_point_start (int)
        - start_list (array)
    - gene name 2
        - ...
        - ...
     - ...
 - chromosome 2
    - gene 1
       - ...
       - ...
 - ...
 - ...
 - genes
   - gene name 1 (link to matrix 1 / chromosome 1 / gene name 1
   - ...
- matrix 1
  - chromosome 1
     - gene name 1
        - chromosome (String)
        - end_list (array)
        - reference_point_end (int)
        - reference_point_start (int)
        - start_list (array)
     - gene name 2
         - ...
         - ...
      - ...
  - chromosome 2
     - gene 1
        - ...
        - ...
  - ...
  - ...

 - genes
   - gene name 1 (link to matrix 2 / chromosome 1 / gene name 1
   - ...
```
Combination mode 'dual' combines the target regions of one viewpoint region of two matrices.

```bash
- combinationMode = 'dual' (String)
- fixateRange (int)
- mode_preselection (String)
- mode_preselection_value (float)
- peakInteractionsThreshold (float)
- pvalue (float)
- range (int, int)
- truncateZeroPvalues (Boolean)
- type='target' (String)
- matrix 1
   - matrix 2
      - chromosome 1
         - gene name 1
            - chromosome (String)
            - end_list (array)
            - reference_point_end (int)
            - reference_point_start (int)
            - start_list (array)
         - gene name 2
            - ...
            - ...
         - ...
      - chromosome 2
         - gene 1
            - ...
            - ...
      - genes
         - gene name 1 (link to matrix 1 / matrix 2 / chromosome 1 / gene name 1
         - ...
      - matrix 3
         - chromosome 1
            - gene name 1 
               - ...
         - ...
         - genes
            - gene name 1 (link to matrix 1 / matrix 3 / chromosome 1 / gene name 1
            - ...
- matrix 2
   - matrix 3
      - chromosome 1
         - gene name 1
            - ...
         - ...
      - ...
   - ...
      - genes
         - gene name 1 (link to matrix 2 / matrix 3 / chromosome 1 / gene name 1
         - ...
- ...
```
Significant interactions file in the 'single' and 'dual' mode don't have a difference in their structure:

```bash
- combinationMode  = 'single' (String) / 'dual' (String)
- fixateRange (int)
- mode_preselection (String)
- mode_preselection_value (float)
- peakInteractionsThreshold (float)
- pvalue (float)
- range (int, int)
- truncateZeroPvalues (Boolean)
- type='significant' (String)
- matrix 1 
 - chromosome 1
    - gene name 1
         - chromosome (String)
         - end_list (array)
         - gene (String)
         - interaction_data_list (array)
         - pvalue (array)
         - raw (array)
         - reference_point_end (int)
         - reference_point_start (int)
         - relative_position_list (array)
         - start_list (array)
         - sum_of_interactions (float)
         - xfold (array)
    - gene name 2
        - ...
        - ...
     . ---
 - chromosome 2
    - gene 1
       - ...
       - ...
 - ...
 - ...
 - genes
   - gene name 1 (link to matrix 1 / chromosome 1 / gene name 1
   - ...
- matrix 1
  - chromosome 1
     - gene name 1
         - chromosome (String)
         - end_list (array)
         - gene (String)
         - interaction_data_list (array)
         - pvalue (array)
         - raw (array)
         - reference_point_end (int)
         - reference_point_start (int)
         - relative_position_list (array)
         - start_list (array)
         - sum_of_interactions (float)
         - xfold (array)
     - gene name 2
         - ...
         - ...
      . ---
  - chromosome 2
     - gene 1
        - ...
        - ...
  - ...
  - ...

 - genes
   - gene name 1 (link to matrix 2 / chromosome 1 / gene name 1
   - ...
```
##### chicAggregateStatistic

```bash
- type='aggregate' (String)
- matrix_1_matrix_2 
   - matrix 1
      - chromosome 1
         - gene name 1
               - chromosome (String)
               - end_list (array)
               - gene_name (String)
               - interaction_data_list (array)
               - raw_target_list (array)
               - relative_distance_list (array)
               - start_list (array)
               - sum_of_interactions (float)
         - gene name 2
        - ...
        - ...
      - genes
         - gene name 1 (link to matrix_1_matrix_2 / matrix 1 / chromosome 1 / gene name 1
         - ...
   - matrix 2
      - chromosome 1
            - gene name 1
                  - chromosome (String)
                  - end_list (array)
                  - gene_name (String)
                  - interaction_data_list (array)
                  - raw_target_list (array)
                  - relative_distance_list (array)
                  - start_list (array)
                  - sum_of_interactions (float)
            - gene name 2
         - ...
         - ...
      - genes
         - gene name 1 (link to matrix_1_matrix_2/ matrix 2 / chromosome 1 / gene name 1
         - ...
- matrix_2_matrix_2
   - matrix 2
      - ...
   - matrix 3
      - ...
```
##### chicDifferentialTest

```bash
- type='differential' (String)
- alpha (float)
- test: fisher / chi2 (String)
- matrix_1
   - matrix 2
      - chromosome 1
         - gene name 1
               - accepted
                  - chromosome (String)
                  - end_list (array)
                  - gene (String)
                  - interaction_data_list (array)
                  - pvalue (array)
                  - raw_target_list_1 (array)
                  - raw_target_list_2 (array)
                  - relative_distance_list (array)
                  - start_list (array)
                  - sum_of_interactions_1 (float)
                  - sum_of_interactions_2 (float)
               - all
                  - chromosome (String)
                  - end_list (array)
                  - gene (String)
                  - interaction_data_list (array)
                  - pvalue (array)
                  - raw_target_list_1 (array)
                  - raw_target_list_2 (array)
                  - relative_distance_list (array)
                  - start_list (array)
                  - sum_of_interactions_1 (float)
                  - sum_of_interactions_2 (float)
               - rejected 
                  - chromosome (String)
                  - end_list (array)
                  - gene (String)
                  - interaction_data_list (array)
                  - pvalue (array)
                  - raw_target_list_1 (array)
                  - raw_target_list_2 (array)
                  - relative_distance_list (array)
                  - start_list (array)
                  - sum_of_interactions_1 (float)
                  - sum_of_interactions_2 (float)
         - gene name 2
               - accepted
               - all
               - rejected
         - ...
      - genes
         - gene name 1 (link to matrix_1 / matrix 2 / chromosome 1 / gene name 1
         - ...
   - matrix 3
      - chromosome 1
            - gene name 1
                  - accepted
                  - all
                  - rejected
            - gene name 2
         - ...
         - ...
      - genes
         - gene name 1 (link to matrix_1 / matrix 3 / chromosome 1 / gene name 1
         - ...
   - ...
- matrix 2
   - matrix 3
      - ...
```
