# hicAdjustMatrix

Keeps, removes or masks specified regions of a matrix.

```text
usage: hicAdjustMatrix --matrix MATRIX --outFileName OUTFILENAME
                       [--chromosomes CHROMOSOMES [CHROMOSOMES ...] |
                       --regions REGIONS | --maskBadRegions MASKBADREGIONS]
                       [--action {keep,remove,mask}]
                       [--interIntraHandling {inter,intra}] [--help]
                       [--version]
```

This tool adjusts hic matrices by keeping, removing or masking a given list of regions or chromosmes.

## Required arguments

| Flag | Meaning |
|---|---|
| `--matrix MATRIX, -m MATRIX` | The Hi-C matrix to adjust. HiCExplorer supports the following file formats: h5 (native HiCExplorer format) and cool. |
| `--outFileName OUTFILENAME, -o OUTFILENAME` | File name to save the adjusted matrix. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--action {keep,remove,mask}, -a {keep,remove,mask}` | Keep, remove or mask the list of specified chromosomes/regions. keep/remove: These options keep/remove bins of matrix by deleting them. This may cause issue plotting the matrix if several parts of a single chromosome are going to be deleted. In that case, one may consider using the mask option (Default: keep). |
| `--interIntraHandling {inter,intra}, -iih {inter,intra}` | Remove the inter- or intra-chromosomal contacts of the given chromosomes. (Default: None). |
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |
| `--chromosomes CHROMOSOMES [CHROMOSOMES ...], -c CHROMOSOMES [CHROMOSOMES ...]` | List of chromosomes to keep/remove. |
| `--regions REGIONS, -r REGIONS` | BED file which stores a list of regions to keep/remove. |
| `--maskBadRegions MASKBADREGIONS, -mbr MASKBADREGIONS` | Bad regions are identified and masked. |

## Notes

hicAdjustMatrix can mask, remove or keep defined regions from a BED file or given chromosomes. This can be useful to create a smaller Hi-C matrix of a single chromosome for e.g. testing purposes or to remove repetitive chromosome ends before calculating A/B compartments using [hicPCA](hicPCA.md).

### Example usages

```bash
$ hicAdjustMatrix -m matrix.cool --action keep --chromosomes chr1 -o matrix_chr1.cool
```
```bash
$ hicAdjustMatrix -m matrix.cool --action mask --regions mask_regions.bed -o matrix_masked.cool
```
mask_regions.bed

```
chr1    10  30
chr1    50  300
```
