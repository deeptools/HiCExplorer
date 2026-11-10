# hicNormalize

Normalizes matrices to a 0-1 range or to the smallest total read count.

```text
usage: hicNormalize --matrices MATRICES [MATRICES ...] --normalize
                    {norm_range,smallest,multiplicative} --outFileName
                    FILENAME [FILENAME ...]
                    [--multiplicativeValue MULTIPLICATIVEVALUE]
                    [--setToZeroThreshold SETTOZEROTHRESHOLD] [--help]
                    [--version]
```

Normalizes given matrices either to the smallest given read number of all matrices or to 0 - 1 range. However, it does NOT compute the contact probability.

We recommend to compute first the normalization (with hicNormalize) and correct the data (with hicCorrectMatrix) in a second step.

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]` | The matrix (or multiple matrices) to get information about. HiCExplorer supports the following file formats: h5 (native HiCExplorer format) and cool. |
| `--normalize {norm_range,smallest,multiplicative}, -n {norm_range,smallest,multiplicative}` | Normalize to a) 0 to 1 range, b) all matrices to the lowest read count of the given matrices (Default: smallest). |
| `--outFileName FILENAME [FILENAME ...], -o FILENAME [FILENAME ...]` | Output file name for the Hi-C matrix. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--multiplicativeValue MULTIPLICATIVEVALUE, -mv MULTIPLICATIVEVALUE` | Value to multiply if --normalize is set to multiplicative. (Default: 1). |
| `--setToZeroThreshold SETTOZEROTHRESHOLD, -sz SETTOZEROTHRESHOLD` | A threshold: after normalization every value smaller than the threshold is set to 0. The default is 0, which has no effect. In the mode "smallest" the counts of the larger matrices are scaled down, so many small values fall below a threshold such as 1.0 and are removed. This increases the sparsity of the matrix. A matrix with sparse inter-chromosomal contacts, which mostly consist of single counts, then loses most of these contacts, except in the smallest matrix, which is not scaled. Set a threshold only when this is intended, and compare the number of inter-chromosomal contacts of the matrices afterwards. |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## Notes

#### Background

To be able to compare different Hi-C interaction matrices the matrices need to be normalized to a equal level of read coverage or
value ranges. hicNormalize accomplish this by offering two modes: 0-1 range normalization or a read count normalization.

#### Usage example

##### Normalize to 0-1 range

norm_range mode of hicNormalize normalizes all reads to the 0 to 1 range i.e. the maximum value of the interaction matrix
becomes 1 and the minimum value 0.

```bash
$ hicNormalize -m matrix.cool --normalize norm_range -o matrix_0_1_range.cool
```
##### Normalize to smallest read count

All matrices are normalized in the way the total read count of each matrix is equal to the read count
of the matrix with the smallest read count of all input matrices.

### Example

- matrix.cool with a read count of 10000
- matrix2.cool with a read count of 12010
- matrix3.cool with a read count of 11000

In this example each entry in matrix2.cool and matrix3.cool are normalized with a factor of 12010 / 10000 respective with 11000 / 10000.

```bash
$ hicNormalize -m matrix.cool matrix2.cool matrix3.cool --normalize smallest 
  -o matrix_normalized.cool matrix2_normalized.cool matrix3_normalized.cool
```
