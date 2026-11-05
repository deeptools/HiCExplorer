# hicFindRestSite

Identifies the genomic locations of restriction sites.

```text
usage: hicFindRestSite --fasta mm10.fa --searchPattern AAGCTT -o rest_site_positions.bed
```

Identifies the genomic locations of restriction sites.

## Required arguments

| Flag | Meaning |
|---|---|
| `--fasta FASTA, -f FASTA` | Path to fasta file for the organism genome. |
| `--searchPattern SEARCHPATTERN [SEARCHPATTERN ...], -p SEARCHPATTERN [SEARCHPATTERN ...]` | Search pattern. For example, for HindIII this pattern is "AAGCTT". Both, forward and reverse strand are searched for a match. The pattern is a regexp and can contain regexp specif syntax (see https://docs.python.org/2/library/re.html). For example the patternCG..GC will find all occurrence of CG followed by any two bases and then GC. |
| `--outFile OUTFILE, -o OUTFILE` | Name for the resulting bed file. |

## Optional arguments

| Flag | Meaning |
|---|---|
| `--help, -h` | show this help message and exit |
| `--version` | show program's version number and exit |

## Notes

#### Further usage

In case multiple restriction enzymes are used in one experiment, `hicFindRestSite` can be used to find restriction sites individually per enzyme. Afterwards, all output bed files should be combined. However, it should be noted that the QC report will not be correct for this specific usage.
