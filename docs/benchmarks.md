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
| hicPCA | *(pending, orchestrating session will fill in)* | | | | | |

## Accuracy

hicBuildMatrix and hicFindTADs were additionally checked for output size and domain count between the
two runs (matrix dimensions and file size for hicBuildMatrix, boundary and domain count for
hicFindTADs), and matched. This is a runtime sanity check on this benchmark run, not a full correctness
validation; the project's actual correctness validation is the per-tool equivalence-class work recorded
in `cpp/STATUS.md` in the repository, which checks every ported tool's output against the Python
reference implementation on real and synthetic data before it is considered done.
