# HiCExplorer v4 (C++) - per-tool progress ledger

Owner: the orchestrating session; the original supervising agent is no longer
running. Architecture: `cpp/PLAN.md`. Rules: `cpp/AGENTS_CONTRACT.md`.
Optimization rules: `cpp/OPTIMIZATION.md`. Last updated 2026-09-13.

**Current state: 20 of 46 tools ported and committed** on `version4-cpp`. The
harness holds 263 cases for those tools plus the tier 0 round trips; 260 pass
with the memory and CPU-time gates live. All 3 failures are `hicPCA`, see its
row. In progress, each in its own git worktree so no two agents share a source
tree: `hicDifferentialTAD` and `hicInterIntraTAD` on branch `v4-tad-stats`;
`chicQualityControl`, `chicViewpointBackgroundModel` and `chicViewpoint` on
branch `v4-chic-foundation`; a `hicPCA` fix on `version4-cpp`.

**Provenance of the facts in this file.** "(reproduced)" marks a result the
orchestrating session reproduced itself. Unmarked results come from the
implementing agent's report and the tests and case files it committed. The
distinction matters: a monitoring filter that matched only lowercase `fail`
once hid three `FAIL` verdicts, and commit 882164e5 claimed a pass count that
only held because most `hicPCA` cases had not been run.

## Legend

- **Tier**: porting tier from `PLAN.md` section 6. Tier 0, the core library,
  has no tools and is tracked in its own table below.
- **Class**: equivalence class from `PLAN.md` 5.1. **Every tool must clear ED:
  each item within three significant digits, relative `1e-3`, set by the project
  owner 2026-09-01 (`PLAN.md` 5.0). Byte identity is not required.** A stricter
  class here means the tool reaches it and is checked at that level. EN is
  reserved for outputs whose Python reference is not reproducible against
  itself, where the tolerance is measured from repeated reference runs rather
  than chosen. E0 byte identical, E1 HDF5 structural, E2 value exact, E3
  `1e-12` relative, E4 `1e-6` relative, E5 set agreement (Jaccard `>= 0.99`), E6
  image, E7 not equivalent by design.
- **Mem**: peak-RSS budget `alpha*W + beta*D + C`, `W` the stored
  upper-triangle CSR working set, `D` the largest dense per-chromosome block,
  `C = 64 MB`, all in SI MB (`PLAN.md` 4.5). A hard gate in the harness.
- **Py test**: state of the Python test **before** the port: `none`, `weak`
  (assertions that cannot constrain the numbers), `partial`, `good`.
- **Char. test**: `written` once a characterization test pinning the Python
  behaviour is committed; `yes` where one is still required; `no` where the
  existing test is adequate.
- **Port**: `not started`, `in progress`, `done`.
- **Equiv**: harness result as passing cases over total, with the gates live.

## Tools

| Tool | Tier | Class | Mem | Py test | Char. test | Port | Equiv | Notes |
|---|---|---|---|---|---|---|---|---|
| hicInfo | 1 | E0 | 0W cool / 1.05W h5 | none | written | done | 25/25 | Byte identical on all 184 cool and h5 matrices in `test_data`; reproduced on six of them, including the `--no_metadata` float path. `--no_metadata` on `gm12878_chr1.cool` peaks at 829.2 MB against an 842.9 MB budget, 98.4 percent and the tightest case in the corpus, and against 2,854 MB for Python. An earlier figure of 809 MB mixed KiB and SI MB. Pinned quirks: `--outFileName` is reopened per matrix, so with several `-m` only the last block survives; a `bin-size` attribute holding the string `"null"` prints no `Bin_length`. Open: no mcool input case. |
| hicConvertFormat | 1 | E1 cool / E2 h5 / E0 text | 1.3W | partial | written | done | 17/17 | Inputs h5, cool, hicpro and 2D text; outputs cool, h5, homer, ginteractions, hicpro and mcool. Text is byte identical after gunzip. Peak RSS 5x to 10x below Python. `.hic` input and `--chromosome` refuse with an explicit error (see deviations). Reproduces F15. Open: the h5 writer emits `/intervals/extra_list` as float64 where Python infers int64 for hicpro and homer sources; homer input has no harness case because `test_data` holds no homer file. |
| hicSumMatrices | 1 | E1 cool / E2 h5 | **3.0W** | none | written | done | 11/11 | `alpha` raised from 2.2 to 3.0: scipy's binary operations drop exact-zero results, so the output is neither a superset nor a subset of either operand and the peak is `nnz_a + nnz_b + nnz_result`. The largest case measured 100.6 percent of a 2.2W budget and 81.7 percent of 3.0W. 995 MB to 167 MB on the 4.2 M nonzero GSM pair. Reproduces F18. |
| hicCompareMatrices | 1 | E1/E2, ED on `log2ratio` | **3.0W** | partial | written | done | 18/18 | Bit exact except `--operation log2ratio`, where numpy's own float64 `log2` differs from libm by at most `3.2e-16` relative. `alpha` 3.0 for the same reason as hicSumMatrices. |
| hicAdjustMatrix | 1 | E1/E2 | 1.3W | partial | written | done | 21/21 | No longer needs the both-triangles exemption: `select_bins` scatters each entry straight into the upper triangle of `P A P^T`, 140 to 57 MB on the duplicate-bin case. Reproduces F14 and F28. |
| hicMergeMatrixBins | 1 | E1/E2, ED on one case | 1.3W, see note | partial | written | done | 13/13 | Bit identical at `-nb 5` on `Li_et_al_2015.h5` (reproduced). The float `--runningWindow` case is ED because scipy's `csr_sort_indices` leaves the order of duplicates undefined, `1.8e-19` relative. `--runningWindow` produces a denser matrix than it reads, so `1.3W` of the input is the wrong budget shape; it passes only because `C` is large next to a 20 MB input and needs a term in the output nnz. Reproduces F17 and F19. |
| hicFindRestSite | 2 | E0 | n/a | good | written | done | 7/7 | Python 5.0 s CPU against under 0.01 s (reproduced from the verification logs). Sorts in memory under `LC_ALL=C` byte order in place of GNU `sort`. |
| hicMergeLoops | 2 | E0 | n/a | partial | written | done | 5/5 | Python 13.6 s CPU against 0.28 s (reproduced from the verification logs). |
| hicValidateLocations | 2 | E0 | 1.3W | good | written | done | 6/6 | Python 93.2 s CPU against 0.60 s, median peak RSS 177 to 51 MB (reproduced from the verification logs). The 313,762-bin, 93-contig input is the many-contig stress case. |
| hicCreateThresholdFile | 2 | E0 | n/a | partial | written | done | 5/5 | `--resolution` now covered. |
| hicMergeTADbins | 2 | E1/E2 | 1.3W | none | written | done | 4/4 | 404.6 to 34.6 MB. Reproduces F16 and a bare `TypeError` for a domain outside the matrix (`hicMergeTADbins.py:129`). |
| hicAverageRegions | 2 | E2 | 1.3W | partial | written | done | 12/12 | Output through the new `npz_file` core component, which writes scipy `save_npz`. Python 20.5 s CPU against 0.07 s (reproduced from the verification logs). |
| hicNormalize | 2 | E1/E2 | 1.2W | good | written | done | 9/9 | Python 20.2 s CPU against 0.25 s (reproduced from the verification logs). |
| hicTransform | 3 | ED or better | 1.2W obs_exp / 1.1W + 1.15D pearson | weak | written | done | 23/23 | `--method pearson` on `Li_et_al_2015.h5` falls from 5,532 MB to 97 MB and from 64.9 s to 3.9 s CPU (reproduced from the verification logs): the dense block is held once and streamed into the output CSR instead of five or six live dense copies and a LIL accumulator. The original 11 tests asserted to the nearest integer. Reuses the transform components written for hicPCA. |
| hicCorrectMatrix | 3 | ED or better ICE / **EN** KR | 1.2W | weak | written | done | 20/20 | krbalancing is reimplemented, not vendored, so F2, F3, F5 and F6 are not inherited. On `gm12878_chr1.cool` ICE goes from 9,331 to 830 MB and 144.9 to 26.9 s CPU, KR from 9,409 to 829 MB and 65.9 to 18.3 s, both at 87 percent of the 954 MB budget (reproduced). `--compatMode v3` reproduces the float32 behaviour of F1 and F2; on `Li_et_al_2015.h5` it shifts all 3,313,962 values by `1.65e-04` relative against the float64 default, with an identical pattern (reproduced). F4 is reproduced through the call ordering. `diagnostic_plot` refuses explicitly, since it belongs to tier 7. |
| hicPCA | 3 | ED or better | 1.1W + 1.15D | weak | written | done, **3 cases failing** | **9/12** | **Correction to commit 882164e5**, which reported all cases passing when only 2 of 12 had run. Rerun against that revision in a separate worktree, the same 3 fail as at HEAD (reproduced): (a) `dist_norm --ligation_factor`: on chrX the C++ ranks a different eigenvector first, and its vector 2 equals Python's vector 1 exactly, because the Python takes `linalg.eig` columns unsorted (F8); (b) `--ignoreMaskedBins` with `--chromosomes`: the Python silently undoes the masking (F25), the C++ does not, 445 lines against 449; (c) a single round-off value around zero in the whole-genome intermediate matrix, `-9.3e-18` against `-1.2e-19`, which the relative ED rule cannot judge (see known open work). A fix for (a) and (b) is in progress. `mm9_reduced_chr1.cool`: 7,923 s to 369 s CPU and 4,185 to 1,585 MB against a 2,509 MB budget. |
| hicCompartmentalization | 3 | E3 | 1.3W | weak | written | not started | - | Characterization tests committed in af76573c, 8 passed and 1 xpassed (reproduced). No C++ yet. |
| hicInterIntraTAD | 3 | E3 | 1.3W | weak | yes | in progress | - | Branch `v4-tad-stats`. `are_files_equal` at `test_hicInterIntraTAD.py:55` is called without `assert`, so its result is discarded; the only asserted check is an `xfail` image comparison. |
| hicPlotSVL | 3 | E3 data / E6 plot | 1.3W | partial | yes | not started | - | 1 test, `are_files_equal(delta=2)` on the two text outputs; the image comparison is commented out at `:66-67`. Untested: `--distance`, `--chromosomes`, `--threads`, `--colorList`. Plot step follows the tier 7 rule. |
| hicBuildMatrix | 4 | E0 BAM and QC / E1 cool / E2 h5 | `2*nnz_out*12 + threads*64 MB + C` | partial | written | done | 17/17 | New htslib BAM subsystem. The output BAM is byte identical, which also pins the record order `--outBam` produces, and the QC folder is byte identical, which is the contract `hicPrepareQCreport` reads. 21.1 s to 0.5 s CPU and 936 to 142 MB on the small_test BAMs (reproduced). Reproduces F29 and F30. The same build written as cool drops per-bin coverage (F26). Baseline: one failure and one error in the Python `trivial_runs` on master; the diagnosis was not recorded before that agent's session ended. |
| hicBuildMatrixMicroC | 4 | as hicBuildMatrix | as hicBuildMatrix | weak | **yes, not written** | done | 4/4 | Shares the BAM subsystem. **Landed without a characterization test**, on harness cases alone, contrary to contract rule 1: `test_hicBuildMatrixMicroC.py` was not strengthened and still leaves 10 of 16 options untested. `--region` behaves as documented here, because without a cut site file F30 cannot occur. |
| hicQuickQC | 4 | E0 | 1.05W | good | no | not started | - | 1 test with strict line equality on `QC.log`; all options covered. Shares `createMatrix` with hicBuildMatrix, so it is now unblocked. |
| hicFindTADs | 5 | E0 text / E1 cool / E2 h5 | per case, `2.2 * Wz + C` | partial | written | done | 8/8 | The original 4 tests caught 2 of 10 source mutations: their comparator used `zip()`, which stops at the shorter file, and allowed 10 differing characters per line. The replacement kills all 10. Text byte identical, z-score matrix E2 in h5 and E1 in cool. `small_test_matrix.h5 --fdr`: 22.2 to 1.0 s CPU, 890 to 114 MB. The budget is declared per case against the z-score band `Wz` (54.7 MB on a 0.7 MB input), not against the input. An AVX2 pairwise kernel is bit identical to the scalar reference; an AVX-512 kernel was 24 percent slower on Zen 4 and was reverted. F22. |
| hicDetectLoops | 5 | E0 text / EN fitted `size` / E5 calls | 2.2W | partial | written | done | 13/13 | The original tests computed nothing (F23). Across 19 configurations and 2,161 called loops: pooled Jaccard 1.000000, no p-value differences, 18 of 19 byte identical, the exception differing only in line order. Four numpy float32 promotions that decide which loops are called are reproduced behind `NBinomPrecision`. The fitted `prob` clears ED but not the tighter EN envelope. On `gm12878_chr1.cool` peak RSS is **877 MB against Python's 619 MB** (budget 1,696 MB), the one case where the port uses more memory, because the cool reader loads the full pixel table before cutting the band. |
| hicDifferentialTAD | 5 | E5 | 1.3W x 2 | good | no | in progress | - | Branch `v4-tad-stats`. 16 tests, exact, across all `-m` and `-mr` modes and thread counts 1, 4 and 11; `--pValue` never varied. |
| hicMergeDomains | 5 | E0/E5 | n/a | weak | yes | not started | - | **2 of its tests fail on master** in the baseline run (`test_main_two_file_protein`, `test_main_two_file_no_protein`), and every `are_files_equal` call (`:70,84-85,106-107`) lacks `assert`. The reference behaviour has to be established before a port. Needs a `scipy.cluster.hierarchy.linkage` reimplementation and a DOT writer. |
| hicAggregateContacts | 5 | E3 `.tab` / E6 plot | 2.2W | weak | yes | not started | - | All 12 general tests are both `xfail` and `skipif(4 GB > memory)`; the 162 `trivial_runs` items assert nothing. Needs `KMeans(random_state=0)`. |
| chicQualityControl | 6 | E0 text / E6 plot | 1.3W | weak | yes | in progress | - | Branch `v4-chic-foundation`. 1 test, `xfail`; `:66` compares the `_failed_reference_points` output against the `_report` reference, so that output is never checked. |
| chicViewpointBackgroundModel | 6 | **EN** fitted parameters | 1.3W | weak | yes | in progress | - | Branch `v4-chic-foundation`. EN rather than E4 because of F24. The existing test allows 700 and 1000 mismatching values at `eps=0.1`. |
| chicViewpoint | 6 | E1 output / E3 values | 1.3W | weak | yes | in progress | - | Branch `v4-chic-foundation`. 2 tests assert only HDF5 keys and parameter defaults; no numeric comparison. |
| chicSignificantInteractions | 6 | E5 calls | 1.3W | weak | yes | not started | - | Builds on the chicViewpoint core. 3 tests assert structure and attribute echo-back only. Uses `pybedtools` at `:515`. |
| chicAggregateStatistic | 6 | E1 | 1.3W | weak | yes | not started | - | 5 tests assert group names and lengths only. |
| chicDifferentialTest | 6 | E3 p-values / E5 calls | 1.3W | weak | yes | not started | - | No p-value comparison in the existing tests. Needs `fisher_exact`, `chi2_contingency` and `chi2.ppf`. |
| chicExportData | 6 | E0 text / E1 bigwig | 1.3W | partial | yes | not started | - | Text at `delta=1`, bigwig at `decimal=0`. The three branches of `--outputValueBigwig` at `chicExportData.py:158-162` are dead code. |
| hicPlotMatrix | 7 | E6 | 1.1W + 1.15D | weak | yes | not started | - | All 32 tests are `xfail` and `skipif` on memory, six behind a 120 GB gate; nothing in this file has run in CI. |
| hicPlotTADs | 7 | E7 | n/a | none | no | not started | - | A 9-line delegation to `pygenometracks.plotTracks.main`; nothing to port. |
| hicPlotViewpoint | 7 | E0 data / E6 plot | 1.3W | weak | yes | not started | - | All 6 tests `xfail`; `:83-97` asserts against stale file names. |
| hicPlotAverageRegions | 7 | E6 | n/a | weak | yes | not started | - | All 4 tests `xfail`, image only. |
| hicPlotDistVsCounts | 7 | E3 data / E6 plot | 1.3W | weak | yes | not started | - | The only real assertion is a PNG byte-size difference below 2000. |
| hicCorrelate | 7 | E3 matrix / E6 plot | 1.3W x n | weak | yes | not started | - | Both tests correlate a file with itself, so 1.0 is guaranteed. Needs complete-linkage clustering with matching leaf order. |
| hicPrepareQCreport (alias hicQC) | 7 | E0 tables / E6 plots | n/a | none | yes | not started | - | No test file, although `test_data/QC*/` holds six reference outputs. Its input, the hicBuildMatrix QC folder, is now byte identical. |
| chicPlotViewpoint | 7 | E6 | 1.3W | weak | yes | not started | - | All 4 tests `xfail`. |
| hicTADClassifier | 8 | E7 | 1.3W features | weak | yes | not started | - | The shipped `.BIN` models are pickles of `imblearn` and `cleanlab` classes with no ONNX converter. |
| hicTrainTADClassifier | 8 | E7 | 1.3W features | weak | yes | not started | - | 14 of 24 options untested. F13. |
| hicHyperoptDetectLoops | 8 | E7 | inherits hicDetectLoops | partial | yes | not started | - | TPE search is RNG dependent. |
| hicHyperoptDetectLoopsHiCCUPS | 8 | E7 | n/a | weak | yes | not started | - | 1 test, `xfail` and doubly `skipif` on `nvcc` and `juicer.jar`, so it never runs. Its `are_files_equal` (`:29`) is a real line-by-line comparator, but the call at `:64` has no `assert`, so its result is discarded. Shells out to `java -jar juicer.jar hiccups`. |

**Tiers 7 and 8 await a decision by the project owner.** `PLAN.md` recommends a
Python plotting shell over the C++ core for the 8 plotting tools, and C++
feature extraction with inference left in Python for the ML tools. That
recommendation has not been confirmed, and it decides whether the port covers 34
or 46 tools.

## Counts

| tier | tools | ported | in progress | not started |
|---|---|---|---|---|
| 1 | 6 | 6 | 0 | 0 |
| 2 | 7 | 7 | 0 | 0 |
| 3 | 6 | 3 | 1 | 2 (hicCompartmentalization has tests) |
| 4 | 3 | 2 | 0 | 1 |
| 5 | 5 | 2 | 1 | 2 |
| 6 | 7 | 0 | 3 | 4 |
| 7 | 8 | 0 | 0 | 8 |
| 8 | 4 | 0 | 0 | 4 |
| **total** | **46** | **20** | **5** | **21** |

Characterization tests written: 20 tools, 19 of the 20 ported plus hicCompartmentalization; hicBuildMatrixMicroC was ported without one. Originally required before porting:
39 of 46.

**Baseline of the Python suite**, first full run on this branch (2026-09-01,
2 h 11 min): 498 passed, 3 failed, 1 error, 1 skipped, 23 xfailed, 62 xpassed
(reproduced). The failures are `test_hicMergeDomains.py::test_main_two_file_protein`,
`::test_main_two_file_no_protein`, and
`test_hicBuildMatrix_trivial_runs.py::test_build_matrix_restrictionCutFile_two`;
the error is `test_hicBuildMatrix_trivial_runs_2.py::test_build_matrix_restrictionCutFile_six`.

## Memory budgets on the designated large inputs

`budget = alpha*W + beta*D + C`, `C = 64 MB`, SI MB. Python peaks measured on
this machine against the reference oracle.

| tool and input | `W` | `D` | Python peak | budget | measured C++ |
|---|---|---|---|---|---|
| `hicCorrectMatrix` KR, gm12878_chr1.cool | 741.9 | - | 9,409 | 954 | **829 (reproduced)** |
| `hicCorrectMatrix` ICE, gm12878_chr1.cool | 741.9 | - | 9,331 | 954 | **830 (reproduced)** |
| `hicPCA`, mm9_reduced_chr1.cool | 5.6 | 762 | 4,185 | 2,509 | **1,585** |
| `hicTransform --method pearson`, Li_et_al_2015.h5 | 19.9 | 987 | 5,532 | 1,221 | **97 (reproduced)** |
| `hicTransform --method pearson`, gm12878_chr1.cool | 741.9 | 4,971 | > 25,000 (est.) | 6,596 | not measured |
| round trip cool to cool, gm12878_chr1.cool | 741.9 | - | 4,048 | 1,028 | **790 (reproduced)** |
| `hicInfo --no_metadata`, gm12878_chr1.cool | 741.9 | - | 2,854 | 842.9 | **829.2, 98.4 percent (reproduced)** |
| `hicSumMatrices`, GSM pair chr1+chr2 (4.2 M nnz) | 50.5 | - | 995 | 215.6 at alpha 3.0 | **176.2** |
| `hicBuildMatrix`, small_test BAMs | - | - | 936 | 321 | **142 (reproduced)** |
| `hicFindTADs --fdr`, small_test_matrix.h5 | 54.7 (z-score band) | - | 890 | 184 | **114 (reproduced)** |
| `hicDetectLoops`, gm12878_chr1.cool | 741.9 | - | 619 | 1,696 | **877, above Python** |

## Dual-mode tools

| tool | `v3` reproduces | default | measured `v3` against default |
|---|---|---|---|
| `hicCorrectMatrix --correctionMethod KR` | krbalancing's float32 input rounding (F2) and float32 rescale accumulators (F1), summed in a fixed order so the C++ stays deterministic | float64 throughout, `v4` | `Li_et_al_2015.h5`: identical pattern, all 3,313,962 values shifted by `1.65e-04` relative, none by more than `1e-3` (reproduced) |
| `hicPCA` | under revision: matching the Python's eigenvector column order is a correctness requirement the ED gate cannot absorb, and the fix in progress may take LAPACK `dgeev` for ordering | - | pending the fix |

The `--perchr` h5-versus-cool correction factor split (F4) is preserved in
both modes, and neither mode reproduces krbalancing's `exit(0)` (F5).

## Known open work in the landed code

| where | issue | effect |
|---|---|---|
| `cpp/scripts/comparators/base.py` | the relative ED rule has no floor for round-off around zero | a value of `-9.3e-18` against `-1.2e-19` in a correlation matrix fails, although both are zero to machine precision. Fails `hicPCA.intermediateMatrices.h5.wholeGenome`. Needs a floor derived from machine epsilon and the dataset's scale, applied as a named rule, never widened to fit |
| `cpp/tools/hicPCA.cpp` | failures (a) and (b) in the hicPCA row | fix in progress |
| all tools | `PLAN.md` 5.0.1 v4 provenance strings are recorded but not implemented | output still names HiCExplorer 3.7.x in its provenance fields. To land across all tools at once, after the worktree branches are merged |
| HDF5 writers | `H5Pset_obj_track_times` left at its default | two runs a second apart differ in bytes (F27). Same timing as the provenance change |
| `core/include/hicx/sparse_matrix.hpp` | values held as `double` regardless of the stored dtype | an int32 or float32 matrix costs twice what it needs; `W` for `gm12878_chr1.cool` would fall from 741.9 to 494 MB |
| `core/src/cool_file.cpp` | no hyperslab read off `indexes/bin1_offset`; the whole pixel table is read before a band or chromosome is cut | the reason hicDetectLoops uses more memory than Python on `gm12878_chr1.cool`, and why the single-chromosome load reads chromosomes it does not need |
| budget formula | `--runningWindow` and hicFindTADs do not fit `alpha * W_input` | both need a term in the working set they actually build |
| `core/src/lbfgsb.cpp` | a projected L-BFGS-B, not a translation of the Fortran | hicDetectLoops' fitted `prob` clears ED but not its EN envelope; closing that means translating the Fortran line search |
| `cpp/scripts/equiv.py` | no `noise.py` EN driver, no dual-mode runs, no `image`, `bigwig` or `bam` comparator | EN cases are currently measured by per-tool scripts |
| argument parsing | no argparse compatibility layer | unambiguous prefix abbreviations such as `--matr` are rejected; `--help` wrapping is close but not identical |

## Findings against the Python reference

Behaviour found while porting that the C++ reproduces, pins, or deliberately
departs from. Most are defects worth reporting upstream. F4 is intended
behaviour, kept here because the port reproduces it on purpose.

| id | where | finding |
|---|---|---|
| **F1** | `krbalancing.cpp:228-229` | The rescale sums accumulate in float32 inside an `omp critical` within a `parallel for`, which serialises the additions without fixing their order, so **KR is not reproducible run to run**: three runs on `gm12878_chr1.cool` gave normalisation factors 0.00660838, 0.00663618 and 0.00666417, a spread of `8.4e-03` relative (reproduced on `Li_et_al_2015.h5`: 0.0190859, 0.01909, 0.019089). |
| **F2** | `krbalancing.cpp:12,27` | float64 input values are downcast to float32 on load. |
| **F3** | `krbalancing.cpp:11-47` | the matrix is staged through a triplet vector, a reserved destination, `setFromTriplets` and `A = A + I`: four copies against a 742 MB working set, most of the 9.2 GB peak. |
| **F4** | `hicCorrectMatrix.py:715-732` | **Intended, confirmed by the project owner 2026-09-01.** h5 historically could not carry correction factors, so h5 output is written already corrected, while cool output keeps raw counts with the factors as a weight column. On `--perchr` the call ordering makes the factors rescaled for h5 and unrescaled for cool. The port reproduces the call ordering. |
| **F5** | `krbalancing.cpp:115-119` | `exit(0)` after 300 outer iterations: the host process ends with a success status and no output. |
| **F6** | `krbalancing.cpp:212-221`, `krbalancing.hpp:29` | `omp critical` around a body with no shared state; `num_threads` hardcoded to 10. |
| **F7** | `hicCorrectMatrix.py:722,724` | `astype(..., copy=False)` copies anyway when the dtype differs: 989 MB each on `gm12878_chr1.cool`. |
| **F8** | `hicPCA.py:305,316` | `scipy.linalg.eig`, the general solver, on a symmetric matrix, with columns taken by position and never sorted, so which vector counts as first depends on LAPACK's internal ordering. |
| **F9** | `hicmatrix/lib/h5.py:84` | `distance_counts` is read from `/correction_factors`. |
| **F10** | `hicmatrix/HiCMatrix.py:58-59` | `correction_factors` and `distance_counts` are unpacked swapped. |
| **F11** | `hicmatrix/HiCMatrix.py:902-919` | `truncTrans` is a no-op. |
| **F12** | `bin/hicFindEnrichedContacts` | imports a module that does not exist (reproduced). |
| **F13** | `lib/tadClassifier.py:782` | `imblearn.base` used without `import imblearn`. |
| **F14** | `hicAdjustMatrix.py:161-165` | `--maskBadRegions` never opens the BED file it is given: `:162` takes `len()` of `--chromosomes`, which is mutually exclusive with the option and always `None`, so h5 input is written back unchanged at exit 0 and cool input raises `TypeError` (reproduced by reading the source). |
| **F15** | `hicConvertFormat --enforce_integer` | on an already corrected matrix every value rounds to zero: `GSM2644945_Untreated-R1.100000_chr1.h5` to cool writes 1,252,980 pixels, all zero, sum zero, exit 0, no warning (reproduced). |
| **F16** | `hicMergeTADbins`, `reduce_matrix(diagonal=True)` | rebuilds the symmetric matrix as `R + R.T - diag(R)` where the subtracted diagonal is the whole within-TAD block sum, losing 22 percent of the total count (30,482,969.64 to 23,695,524.86). |
| **F17** | `hicMergeMatrixBins.py:248` | groups smaller than `numBins/2` are dropped with only a `log.debug`; at `--numBins 20` two whole chromosomes disappear. |
| **F18** | `hicSumMatrices.py:72` | `maskBins` of the union of both inputs' nan bins discards real counts: 5,142 on the GSM2644945 plus GSM2644947 pair. |
| **F19** | `hicMergeMatrixBins --runningWindow` | the window runs over raw bin indices and crosses chromosome borders. |
| **F20** | `hicmatrix/HiCMatrix.py:84-92` | `save` reuses the handler built during `load`, so the output format follows the input regardless of the name given: `-m a.h5 -o out.cool` writes an h5 file named `out.cool.h5` (reproduced by reading the source). A name ending in neither `cool` nor `h5` writes nothing and exits 0. |
| **F21** | consequence of F10 | for cool input the correction factors land in `distance_counts`, which the cool writer never reads, so a KR-balanced cooler is written back with **no weight column and balanced values as its counts**. |
| **F22** | `hicFindTADs.py:345` | `--numberOfProcessors 12` or more aborts with `ValueError` on `small_test_matrix.h5`, because `np.array_split` can hand a worker a range whose bins are all skipped. Every working process count gives identical output. |
| **F23** | `test_hicDetectLoops.py` | the three cool tests pass `--maxLoopDistance 3000000` against 2.5 Mb bins, so only the main diagonal loads, is deleted, and no loop is called; they pass because the comparator uses `zip()`. The fourth test asserts nothing. Separately, the Python's loop output order depends on `--threads`. |
| **F24** | `fit_nbinom` | not reproducible against itself: perturbing the start of scipy's `fmin_l_bfgs_b` by `1e-9` moves the fitted `size` by a median of 8.7 percent and up to 75 percent, a flat likelihood ridge meeting `pgtol=1e-5`. Fitted negative binomial parameters are therefore class EN everywhere they appear. |
| **F25** | `hicPCA.py:253-258`, `hicmatrix/HiCMatrix.py:603` | `--ignoreMaskedBins` masks and enlarges the bins, then `--chromosomes` calls `keepOnlyTheseChr`, whose first action is `restoreMaskedBins()`. The masking is silently undone whenever both options are given (reproduced). |
| **F26** | h5 against cool | the same `hicBuildMatrix` run written as h5 and as cool gives identical matrices, but the cool file loses the per-bin coverage (18,346 distinct values become a single 1.0) and has no nan-bin list, which hicmatrix rebuilds on load as every all-zero row (0 against 14,845). `hicCorrectMatrix --sequencedCountCutoff` therefore runs on h5 input and **crashes with `AssertionError` on cool input** at `hicCorrectMatrix.py:679` (all reproduced; identical in Python and C++). |
| **F27** | HDF5 output, Python and v4 | every object embeds its modification time, 841 across the corpus, so two runs a second apart differ in bytes. |
| **F28** | `hicAdjustMatrix.py:124,144` | the `--regions --action remove` warning calls `getChrBinRange` with the chromosome of the last readable BED line rather than the region reported. |
| **F29** | `buildMatrixMethods.py:744-762` | `--keepSelfCircles` is a no-op, because the `continue` meant to drop a self circle sits inside the loop over restriction sequences; with two enzymes a self circle is counted once per enzyme. |
| **F30** | `hicBuildMatrix --region` | empties the restriction site list, because `bed2interval_list` keeps a site only when `region_end <= site_end` and `region_end` is the chromosome length; every close inward pair becomes "same fragment" and no self circle is counted. |

## Tier 0 - core library components

| Component | Path | State | Notes |
|---|---|---|---|
| bins, interval lookup | `core/include/hicx/bins.hpp` | done | |
| CSR matrix, upper-triangle storage | `sparse_matrix.hpp`, `hic_matrix.hpp` | done | int32 indices while they fit; values held as double (open work) |
| HDF5 wrapper and blosc filter | `hdf5_util.cpp` | done | compress and decompress; Python reads C++-written h5 back to identical `hicInfo` output (reproduced); `cd_values` byte for byte what PyTables emits |
| cool reader and writer | `cool_file.cpp` | done | writer at E1; `/bins/weight` is float64, deflate 6, no shuffle, unlimited maxshape; no hyperslab band read yet |
| mcool | `cool_file.cpp` | done | `::` group URIs, both group layouts |
| HiCExplorer h5 reader and writer | `h5_file.cpp` | done | writer at E2; metadata writing refactored into `write_h5_metadata` in af76573c, regression 171 of 171 afterwards (reproduced) |
| homer, ginteractions, hicpro, 2D text | text format components | done | |
| npy/npz | `npz_file.cpp` | done | |
| `.hic` reader | - | not started | `hicConvertFormat --inputFormat hic` refuses explicitly |
| numpy-compatible reductions | `numpy_compat.cpp`, `simd_reduce.*` | done | per-8192-buffer pairwise sum verified over 161 sizes; AVX2 kernel bit identical to scalar over 307 lengths |
| statistics | `stats_ops.*`, `lbfgsb.cpp` | done in part | cephes `ndtr`, `erf`, `erfc`, `digamma`, `gammaln`; `rankdata`, `ranksums`, Benjamini-Hochberg, Bonferroni, `betainc`, `nbinom_sf`, `fit_nbinom`, `NBinomPrecision`. Not yet: `fisher_exact`, `chi2_contingency`, `chi2.ppf`, `pearsonr`, `spearmanr`, `anderson_ksamp` |
| KR and ICE | hicCorrectMatrix components | done | KR reimplemented, not vendored |
| obs/exp, z-score | `obsexp_ops.*` | done | the general path without `maxdepth` throws rather than densifying |
| Pearson, covariance, PCA transforms | `transform_ops.*` | done | |
| `reduce_matrix`, bin merging, adjust operations | `reduce_matrix.*`, `adjust_ops.*`, `matrix_ops.*` | done | |
| BED intervals, bedtools replacement | `bedtools_ops.*`, `text_table.*` | done | |
| FASTA reader | `fasta_reader.*` | done | |
| BAM via htslib | hicBuildMatrix components | done | |
| libBigWig | `hicx::bigwig` target | done | used by hicPCA |
| clustering (complete linkage, k-means) | - | not started | for hicCorrelate, hicMergeDomains, hicAggregateContacts |
| deterministic thread pool | `parallel.hpp` | done | fixed partitions, combined in index order |
| peak-RSS self-report | `resource_usage.hpp` | done | |
| argparse compatibility layer | - | not started | |
| equivalence harness | `cpp/scripts/equiv.py`, `comparators/` | done in part | memory gate, CPU-time gate and determinism mode live on every case; comparators `text`, `cool`, `h5`, `interval`, `npz`. Missing: ED round-off floor, `noise.py` EN driver, dual-mode runs, `image`, `bigwig`, `bam` |

## Deliberate deviations from the Python behaviour

| Tool / component | Deviation | Reason | Decided by |
|---|---|---|---|
| acceptance, all tools | three significant digits per item, not byte identity | the owner's gate; SIMD and threading change reduction order | project owner, 2026-09-01 |
| provenance, all tools | output is to name HiCExplorer 4 in its provenance fields (`PLAN.md` 5.0.1) | a file written by v4 must not claim to come from 3.7; recorded, not yet implemented | project owner, 2026-09-02 |
| KR, both modes | `exit(0)` on non-convergence (F5) is not reproduced; the port raises and exits non-zero | a success status with no output breaks every pipeline | orchestrating session |
| KR | deterministic, and class EN rather than a single reference comparison | the reference is not reproducible (F1) | orchestrating session |
| KR | krbalancing reimplemented rather than vendored | vendoring would import F2, F3, F5 and F6 | orchestrating session |
| hicCorrectMatrix, hicPCA | a `--compatMode {v3,v4}` flag that the Python does not have | quantifies what the float32 defects cost without forcing them on users | orchestrating session |
| hicDetectLoops | chromosomes written in a fixed order | the Python's order depends on `--threads` (F23) | implementing agent, reported |
| hicFindTADs | no abort at `--numberOfProcessors 12` or more | F22; an empty partition contributes no rows | implementing agent, reported |
| hicConvertFormat | `.hic` input and `--chromosome` refuse explicitly | `.hic` is tier 4; single-chromosome cooler loading is not in the reader | implementing agent, reported |
| hicCorrectMatrix | `diagnostic_plot` refuses explicitly | matplotlib, tier 7 | orchestrating session |
| all tools | upper-triangle storage with symmetric access | halves the resident matrix; visit order preserved and unit tested | plan |
| tiers 7 and 8 | Python shell over the C++ core; ML inference stays Python | recommended in `PLAN.md`; **not yet confirmed** | pending |
| `hicFindEnrichedContacts` | not ported | F12, dead in the reference | plan |

## Open questions

1. (Resolved.) The blosc compress path produces files PyTables 3.10.1 reads.
2. (Resolved.) `C = 64 MB` holds, although `hicInfo --no_metadata` on
   `gm12878_chr1.cool` is at 98.4 percent of its budget.
3. (In progress, hicPCA fix.) Does LAPACK `dgeev` from the conda prefix
   reproduce `scipy.linalg.eig`'s column order? Failure (a) in the hicPCA row
   turns on it.
4. (Answered in part.) KR `v3` against the default on `Li_et_al_2015.h5`:
   `1.65e-04` relative on every value. Not yet measured on
   `gm12878_raw_values.cool`, where F2 is inert.
5. (Resolved.) Both mcool group layouts are read.
6. (Answered for hicFindTADs.) It needs the z-score band, not both triangles.
   hicDetectLoops and hicAggregateContacts remain open.
7. What is the correct reference behaviour of `hicMergeDomains`, given that two
   of its tests fail on master?
8. Tiers 7 and 8: confirm or replace the recommended strategy (project owner).
