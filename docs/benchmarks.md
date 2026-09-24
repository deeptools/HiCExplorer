# Benchmarks

HiCExplorer 4 was benchmarked against HiCExplorer 3.7.6 on real data: a full
single-chromosome Hi-C matrix (`hicTADClassifier/gm12878_chr1.cool`, 24,926 bins, 61.8M nonzero pixels)
for `hicDetectLoops`, `hicFindTADs` and `hicPCA`, a real paired-end Hi-C sequencing run for
`hicBuildMatrix` (92.7M sequenced reads, DpnII digestion, 10 kb bins) and a real full-genome capture
Hi-C dataset for CHiCAGO. Wall clock and peak resident set size (RSS) were captured with
`/usr/bin/time -v`. The rows of the table below ran single-threaded (`--threads 1` /
`--numberOfProcessors 1` where the tool has a thread flag).

## Runtime and peak memory

| Tool | C++ v4 | Python 3.7.6 | Speedup | C++ v4 peak RSS | Python peak RSS | Memory |
|---|---|---|---|---|---|---|
| hicDetectLoops | 2.05 s | 20.59 s | 10.0x | 230 MB | 624 MB | 2.7x less |
| hicFindTADs | 6.33 s | 33.00 s | 5.2x | 1,475 MB | 6,021 MB | 4.1x less |
| hicPCA (100 kb) | 3.76 s | 24.20 s | 6.4x | 147 MB | 493 MB | 3.4x less |

`hicPCA` is measured on the same chromosome merged to 100 kb bins (2,253 bins), not the 10 kb matrix
the other rows use. A/B compartment calling, what `hicPCA` computes, is a large-scale feature normally
analyzed at 100 kb to 1 Mb resolution, never at 10 kb; an initial run at 10 kb was abandoned after both
implementations were still running past 1.5 hours (dense, full-spectrum eigendecomposition of a roughly
25,000 by 25,000 matrix is expensive at any implementation and thread count, and is not how this tool
is actually used). See the C++ source's own extensive comments in `cpp/tools/hicPCA.cpp` for why the
dense solver runs single-threaded by design rather than in parallel: multithreaded OpenBLAS was
measured to change the bit pattern of its output, including which eigenvector is returned first when
eigenvalues are nearly equal, which would break byte-identical reproducibility against the Python
reference. A faster, lower-memory `--eigenSolver lanczos` mode exists for large matrices, but it is a
different algorithm from Python's dense path and therefore not directly comparable here.

## hicBuildMatrix

Input: two SAM files with 92,684,960 sequenced reads, DpnII restriction sites, 10 kb bins. The times are
the logged wall clock of one run per setting.

| Threads | C++ v4 | Python 3.7.6 | Speedup | C++ v4 peak RSS | Python peak RSS |
|---|---|---|---|---|---|
| 1 | 16 min 22 s | 42 min 18 s | 2.6x | 6.1 GB | 14.9 GB |
| 4 | 10 min 38 s | 23 min 14 s | 2.2x | 6.1 GB | 15.7 GB |
| 8 | 9 min 22 s | 21 min 03 s | 2.2x | 6.1 GB | 16.4 GB |

The QC report is identical between the two implementations at every thread count: 92,684,960 sequenced
reads, 73,506,771 usable pairs and 57,566,986 Hi-C contacts. An earlier version of this page listed
122.0 s against 700.1 s for a 130.6M-pair run. That run has no retained log and is not reproducible from
the recorded inputs, so the row was removed.

## CHiCAGO

Input: a full-genome capture Hi-C design (837,161 HindIII fragments, 22,076 baits) and one sample with
119.8M `.chinput` rows, run with 16 threads. R Chicago runs as a single process.

| Stage | R Chicago | C++ from `.chinput` | C++ from matrix |
|---|---|---|---|
| Read input | 48.2 s | included below | included below |
| Background model | 160.1 s | 120.7 s (20.7 GB) | 29.4 s (3.9 GB) |
| Scores | 229.3 s | 90.4 s (20.7 GB) | 19.1 s (3.9 GB) |
| Significant interactions | included above | 25.5 s (18.1 GB) | included above |
| Total | 437.5 s (40.4 GB) | 236.6 s | 48.5 s |

The matrix column starts from the cool file that `hicConvertFormat --inputFormat chinput` writes from the
same `.chinput` file, with both directions of a bait-to-bait pair summed. It derives cis pairs within the
Brownian estimation distance only. On a 949-bait subset (17.4M rows), R needs 256.8 s at 5.4 GB and the
C++ tools about 72 s with 8 threads. On the two-chromosome test fixture, R needs 4.5 s and the C++ tools
about 2.2 s.

## Accuracy

hicFindTADs was additionally checked for agreement between the two runs on this benchmark's own real data. It matched in boundary and domain count (644 domains,
both runs). The hicBuildMatrix QC agreement is listed above. This is a runtime sanity check on these benchmark runs, not a full correctness validation;
the project's actual correctness validation is the per-tool equivalence-class work recorded in
`cpp/STATUS.md` in the repository, which checks every ported tool's output against the Python reference
implementation on real and synthetic data before it is considered done.
