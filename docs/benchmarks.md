# Benchmarks

HiCExplorer v4 (C++) was benchmarked against Python HiCExplorer 3.7.6 on real data: a full
single-chromosome Hi-C matrix (`hicTADClassifier/gm12878_chr1.cool`, 24,926 bins, 61.8M nonzero pixels)
for `hicDetectLoops` and `hicFindTADs`, and a real small paired BAM dataset for `hicBuildMatrix`. Wall
clock and peak resident set size (RSS) were captured with `/usr/bin/time -v`, both tools run
single-threaded (`--threads 1` / `--numberOfProcessors 1` where the tool has a thread flag).

## Runtime and peak memory

| Tool | C++ v4 | Python 3.7.6 | Speedup | C++ v4 peak RSS | Python peak RSS | Memory |
|---|---|---|---|---|---|---|
| hicBuildMatrix | 0.23 s | 16.06 s | 69.8x | 98 MB | 900 MB | 9.2x less |
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

## Accuracy

hicBuildMatrix and hicFindTADs were additionally checked for output size and domain count between the
two runs (matrix dimensions and file size for hicBuildMatrix, boundary and domain count for
hicFindTADs), and matched. This is a runtime sanity check on this benchmark run, not a full correctness
validation; the project's actual correctness validation is the per-tool equivalence-class work recorded
in `cpp/STATUS.md` in the repository, which checks every ported tool's output against the Python
reference implementation on real and synthetic data before it is considered done.
