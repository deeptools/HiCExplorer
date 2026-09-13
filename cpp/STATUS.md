# HiCExplorer v4 (C++) - per-tool progress ledger

Owner: the orchestrating session. Architecture: `cpp/PLAN.md`. Rules:
`cpp/AGENTS_CONTRACT.md`. Optimization rules: `cpp/OPTIMIZATION.md`. Last updated
2026-09-13, at commit `3e30a156`.

**Current state: 30 of 46 tools ported and committed** on `version4-cpp`, and no
porting branch is open. The last merged-tree regression, on `63dd2424`, passed
358 cases over 29 tools and the tier 0 round trips (reproduced). hicPCA's 13
cases pass as well: all 9 that use `--chromosomes` reproduced on `63dd2424`, the
rest reproduced on `a3c49fd8`, and mm9_reduced_chr1 in the agent's run only.

Separately, **coolercpp**, an independent C++ library API-compatible with Python
cooler, is being built in its own repository at `~/src/coolercpp`. HiCExplorer v4
is to replace its own cool code with it once coolercpp reaches read and write
parity.

**Provenance of the facts in this file.** "(reproduced)" marks a result the
orchestrating session reproduced itself; unmarked results come from an
implementing agent's report and the tests it committed. The distinction has
mattered repeatedly:
- a filter matching only lowercase `fail` once hid three `FAIL` verdicts;
- commit 882164e5 claimed a pass count that held only because most hicPCA cases
  had not run;
- a verification in the orchestrating session's worktree hid three case inputs
  that were never committed;
- only a build from a clean export exposed a unit test that silently depended on
  a `cpp/build` directory inside the source tree;
- a resumed agent reported reruns it could not have performed in 27 seconds;
  its verdicts reproduced, but the report alone was not evidence.

## Legend

- **Tier**: porting tier from `PLAN.md` section 6.
- **Class**: equivalence class from `PLAN.md` 5.1. **Every tool must clear ED:
  each item within three significant digits (`PLAN.md` 5.0).** Byte identity is
  not required, and how a tool reaches its result is free (`PLAN.md` 5.0 item 3).
  E0 byte identical, E1 HDF5 structural, E2 value exact, E5 set agreement
  (Jaccard `>= 0.99`), EN within the measured oracle-noise envelope, E7 not
  compared by bytes, with the reason recorded.
- **Mem**: peak-RSS budget `alpha*W + beta*D + C`, `C = 64 MB`, SI MB. A hard gate.
- **Py test**: state of the Python test **before** the port.
- **Char. test**: `written`, `yes` (still required) or `no` (adequate).
- **Equiv**: passing harness cases over total, gates live.

## Tools

| Tool | Tier | Class | Mem | Py test | Char. test | Port | Equiv | Notes |
|---|---|---|---|---|---|---|---|---|
| hicInfo | 1 | E0 | 0W cool / 1.05W h5 | none | written | done | 25/25 (reproduced) | Byte identical on all 184 cool and h5 matrices in `test_data`, reproduced on six. `--no_metadata` on `gm12878_chr1.cool`: 829.2 MB against 842.9 MB, the tightest budget in the corpus. |
| hicConvertFormat | 1 | E1 cool / E2 h5 / E0 text | 1.3W | partial | written | done | 17/17 (reproduced) | `.hic` input and `--chromosome` refuse explicitly. Reproduces F15. Open: h5 `extra_list` dtype for hicpro and homer sources. |
| hicSumMatrices | 1 | E1/E2 | 3.0W | none | written | done | 11/11 (reproduced) | Reproduces F18. |
| hicCompareMatrices | 1 | E1/E2, ED on `log2ratio` | 3.0W | partial | written | done | 18/18 (reproduced) | |
| hicAdjustMatrix | 1 | E1/E2 | 1.3W | partial | written | done | 21/21 (reproduced) | Reproduces F14, F28. |
| hicMergeMatrixBins | 1 | E1/E2, ED on one case | 1.3W, see open work | partial | written | done | 13/13 (reproduced) | Reproduces F17, F19. |
| hicFindRestSite | 2 | E0 | n/a | good | written | done | 7/7 (reproduced) | |
| hicMergeLoops | 2 | E0 | n/a | partial | written | done | 5/5 (reproduced) | |
| hicValidateLocations | 2 | E0 | 1.3W | good | written | done | 6/6 (reproduced) | |
| hicCreateThresholdFile | 2 | E0 | n/a | partial | written | done | 5/5 (reproduced) | |
| hicMergeTADbins | 2 | E1/E2 | 1.3W | none | written | done | 4/4 (reproduced) | Reproduces F16. |
| hicAverageRegions | 2 | E2 | 1.3W | partial | written | done | 12/12 (reproduced) | |
| hicNormalize | 2 | E1/E2 | 1.2W | good | written | done | 9/9 (reproduced) | |
| hicTransform | 3 | ED or better | 1.2W / 1.1W + 1.15D | weak | written | done | 23/23 (reproduced) | pearson on `Li_et_al_2015.h5`: 5,532 to 97 MB (reproduced). |
| hicCorrectMatrix | 3 | ED ICE / EN KR | 1.2W | weak | written | done | 20/20 (reproduced) | ICE and KR on `gm12878_chr1.cool` from about 9.3 GB to 830 MB (reproduced). |
| hicPCA | 3 | ED or better | 1.1W + 1.15D | weak | written | done | 13/13 (12 reproduced) | Covariance equal to `np.cov` bit for bit, because eigenvectors are chosen by position from a spectrum with a 169-fold largest eigenvalue (F8, F31). Reproduces F25. |
| hicCompartmentalization | 3 | E0 `_dat` / E2 npz | 1.3W | weak | written | done | 15/15 (reproduced) | Required figure: a plain run refuses, a C++-only `--noPlot` writes the numeric outputs (reproduced by hand). |
| hicInterIntraTAD | 3 | E0 | 1.3W | weak | written | done | 5/5 (reproduced) | Explicitly requested ratio plot refused before writing (reproduced by hand). Reproduces F33 to F35. |
| hicPlotSVL | 3 | E0 text | 1.3W | partial | written | done | 14/14 (reproduced) | `gm12878_chr1.cool`: 3,846 to 829 MB against a 1,028 MB budget, after the single-chromosome cooler load was made to cut in place. The boxplot follows rule 7. Reproduces F37 and F51; F38 is a deviation. |
| hicBuildMatrix | 4 | E0 BAM and QC / E1 / E2 | `2*nnz_out*12 + threads*64 MB + C` | partial | written | done | 17/17 (reproduced) | QC folder lacks the Python's PNGs and HTML report. Reproduces F29, F30. |
| hicBuildMatrixMicroC | 4 | as hicBuildMatrix | as hicBuildMatrix | weak | **yes, not written** | done | 4/4 (reproduced) | Landed without a characterization test, contrary to rule 1. |
| hicQuickQC | 4 | E0 | 192 MB | good | written | done | 5/5 (reproduced) | Reproduces F40. |
| hicFindTADs | 5 | E0 text / E1 / E2 | per case, `2.2 * Wz + C` | partial | written | done | 8/8 (reproduced) | F22. |
| hicDetectLoops | 5 | E0 text / EN fitted size / E5 calls | 2.2W | partial | written | done | 13/13 (reproduced) | Uses Cephes `betainc` in both its float64 and float32 paths, where scipy 1.14 uses Boost.Math for both (F42, open work). `gm12878_chr1.cool` peak 877 MB against Python's 619 MB. |
| hicDifferentialTAD | 5 | E0, E5 calls | 1.3W x 2 | good | written | done | 12/12 (reproduced) | Jaccard 1.0 on every case. |
| hicMergeDomains | 5 | E0 text / E7 PDF | none (reads no matrix) | weak | written | done | 11/11 (reproduced) | Its baseline failures were environmental (F41). DOT source byte identical, rendered by the external graphviz `dot`, found through a per-case `path_prepend`. The PDFs are E7 because `dot` embeds a creation date in compressed streams. Reproduces F39 and F49. |
| hicAggregateContacts | 5 | E0 tables / E2 h5 / E1 cool | 1.3W, z-score 559 MB | weak | written | done | 29/29 (reproduced) | KMeans and ward labels identical in all 108 fits the Python made on real data. z-score case 4.0 GB to 335 MB. Required figure handled with `--noPlot`. Reproduces F52 to F58. |
| chicQualityControl | 6 | E0 text | 1.3W | weak | written | done | 5/5 (reproduced) | Figures follow rule 7. Reproduces F36, F43. |
| chicViewpointBackgroundModel | 6 | exact columns; fitted size and prob by likelihood and downstream E5 | 1.3W | weak | written | done | 5/5 (reproduced) | EN withdrawn for size and prob (F50). Each fit's likelihood must be no worse than the worst of five reference runs by more than `max(spread_i, T_well)`; downstream calls Jaccard 1.000000 on all four fitted cases. Reproduces F45. |
| chicViewpoint | 6 | E1 / bit identical values | 1.3W | weak | written | done | 6/6 (reproduced) | Bit identical through Boost.Math `ibeta`, as scipy uses (F42). Reproduces F44, F46, F47; F48 is a deviation. |
| chicSignificantInteractions | 6 | E5 calls | 1.3W | weak | yes | not started | - | Builds on the chicViewpoint core, now merged. |
| chicAggregateStatistic | 6 | E1 | 1.3W | weak | yes | not started | - | |
| chicDifferentialTest | 6 | E3 p-values / E5 calls | 1.3W | weak | yes | not started | - | Needs `fisher_exact`, `chi2_contingency`, `chi2.ppf`. |
| chicExportData | 6 | E0 text / E1 bigwig | 1.3W | partial | yes | not started | - | |
| hicPlotMatrix | 7 | E6 | 1.1W + 1.15D | weak | yes | not started | - | |
| hicPlotTADs | 7 | E7 | n/a | none | no | not started | - | A delegation to pyGenomeTracks. |
| hicPlotViewpoint | 7 | E0 data / E6 plot | 1.3W | weak | yes | not started | - | |
| hicPlotAverageRegions | 7 | E6 | n/a | weak | yes | not started | - | |
| hicPlotDistVsCounts | 7 | E3 data / E6 plot | 1.3W | weak | yes | not started | - | |
| hicCorrelate | 7 | E3 matrix / E6 plot | 1.3W x n | weak | yes | not started | - | Needs complete-linkage clustering with matching leaf order; `hicx::cluster` has ward only. |
| hicPrepareQCreport (alias hicQC) | 7 | E0 tables / E6 plots | n/a | none | yes | not started | - | |
| chicPlotViewpoint | 7 | E6 | 1.3W | weak | yes | not started | - | |
| hicTADClassifier | 8 | E7 | 1.3W | weak | yes | not started | - | |
| hicTrainTADClassifier | 8 | E7 | 1.3W | weak | yes | not started | - | F13. |
| hicHyperoptDetectLoops | 8 | E7 | inherits | partial | yes | not started | - | |
| hicHyperoptDetectLoopsHiCCUPS | 8 | E7 | n/a | weak | yes | not started | - | |

**Tiers 7 and 8 await the project owner's decision.** Until then no C++ tool draws
figures: an explicitly requested figure is refused before any output is written,
one written only under a default name is skipped with a note, and a figure that
is a required output gets a C++-only `--noPlot` (contract rule 7).

## Counts

| tier | tools | ported | not started |
|---|---|---|---|
| 1 | 6 | 6 | 0 |
| 2 | 7 | 7 | 0 |
| 3 | 6 | 6 | 0 |
| 4 | 3 | 3 | 0 |
| 5 | 5 | 5 | 0 |
| 6 | 7 | 3 | 4 |
| 7 | 8 | 0 | 8 |
| 8 | 4 | 0 | 4 |
| **total** | **46** | **30** | **16** |

## Baseline of the Python suite

First full run on this branch, 2026-09-01: 498 passed, 3 failed, 1 error, 1
skipped, 23 xfailed, 62 xpassed (reproduced). None of the four was a defect in a
HiCExplorer tool (reproduced from the tracebacks): two hicMergeDomains tests could
not find graphviz `dot` on the venv's `PATH`; one hicBuildMatrix trivial run passes
`region="ChrX"`, absent from its input; one takes a `keepSelfLigation` argument no
`parametrize` supplies.

## Memory on the designated large inputs

| tool and input | Python peak | budget | measured C++ |
|---|---|---|---|
| hicCorrectMatrix KR, gm12878_chr1.cool | 9,409 | 954 | 829 (reproduced) |
| hicCorrectMatrix ICE, gm12878_chr1.cool | 9,331 | 954 | 830 (reproduced) |
| hicTransform pearson, Li_et_al_2015.h5 | 5,532 | 1,221 | 97 (reproduced) |
| hicPCA, mm9_reduced_chr1.cool | 4,187 | 2,509 | 1,589 |
| hicAggregateContacts z-score | 4,000 | 559 | 335 |
| round trip cool, gm12878_chr1.cool | 4,048 | 1,028 | 790 (reproduced) |
| hicPlotSVL, gm12878_chr1.cool | 3,846 | 1,028 | 829 |
| hicInfo --no_metadata, gm12878_chr1.cool | 2,854 | 842.9 | 829.2 (reproduced) |
| hicCompartmentalization, small_test_matrix.h5 | 2,534 | - | 26 |
| hicBuildMatrix, small_test BAMs | 936 | 321 | 142 (reproduced) |
| hicQuickQC, median over cases | 921 | 192 | 92 |
| hicDetectLoops, gm12878_chr1.cool | 619 | 1,696 | 877, above Python |

## Known open work

| where | issue | effect |
|---|---|---|
| all tools | `PLAN.md` 5.0.1 v4 provenance and `H5Pset_obj_track_times` are recorded, not implemented | output still names HiCExplorer 3.7.x and embeds object modification times (F27). Now unblocked: every porting branch is merged |
| hicDetectLoops | Cephes `betainc` in both the float64 path (`nbinom_sf`) and the float32 path; the comment at `detect_loops_impl.cpp:423` calls scipy's float32 loop "single precision cephes incbet", which is outdated, since scipy 1.14 uses Boost.Math `ibeta` for both | no case differs on the corpus; a p-value within an ulp of a threshold could. Switch both paths to Boost, which needs a float32 overload in `scipy_special` |
| `cpp/tests/` | three bespoke refusal scripts plus a generic one and two shell scripts | consolidate on `refuses_unavailable_output.cmake` |
| hicAggregateContacts | numpy's argsort fallback for CPUs without AVX-512 was not checked against numpy, because the development machine always takes the AVX-512 path | contact-pair line order could differ on such a CPU |
| hicAggregateContacts | z-score in modes `all` and `inter-chr` has no harness case | the Python needs a dense matrix of several GB |
| hicBuildMatrixMicroC | no characterization test | rule 1 not met |
| hicBuildMatrix, hicQuickQC | the QC folder lacks the Python's PNGs and `hicQC.html` | part of the plotting decision |
| packaging | Boost.Math and x86-simd-sort are fetched at configure time; hicMergeDomains needs graphviz `dot` at run time | an offline package build needs vendored tarballs, and the package needs a graphviz dependency |
| harness | `cpp_args` passes options to the C++ side only; `path_prepend` changes one case's `PATH` | both allowed only where they do not change compared outputs |
| `core/src/cool_file.cpp` | no hyperslab band read | hicDetectLoops uses more memory than Python on `gm12878_chr1.cool`. To be revisited when coolercpp replaces this code |
| `core/include/hicx/sparse_matrix.hpp` | values held as double regardless of dtype | int32 and float32 matrices cost twice what they need |
| budget formula | `--runningWindow` and hicFindTADs do not fit `alpha * W_input` | both need a term in the working set they build |
| argument parsing | no argparse compatibility layer | prefix abbreviations rejected; `--help` close but not identical |

Fixed in this round and recorded for history: `test_transform_ops.cpp` wrote its
scratch file into a `cpp/build` directory inside the source tree and failed on an
out-of-tree build from a fresh checkout (fixed in `3e30a156`). The round-off floor
proposed for class ED was not applied. An orphaned mutation-test run that loaded
the machine for eleven days was stopped on 2026-09-13; CPU-time figures are
unaffected by design.

## Findings against the Python reference

| id | where | finding |
|---|---|---|
| **F1** | `krbalancing.cpp:228-229` | KR is not reproducible run to run (reproduced). |
| **F2** | `krbalancing.cpp:12,27` | float64 input downcast to float32. |
| **F3** | `krbalancing.cpp:11-47` | four full copies of the matrix during construction. |
| **F4** | `hicCorrectMatrix.py:715-732` | **Intended, confirmed by the project owner.** h5 output is written corrected, cool keeps raw counts with weights. |
| **F5** | `krbalancing.cpp:115-119` | `exit(0)` on non-convergence. |
| **F6** | `krbalancing.cpp:212-221` | critical section around a body with no shared state. |
| **F7** | `hicCorrectMatrix.py:722,724` | `astype(copy=False)` copies when the dtype differs. |
| **F8** | `hicPCA.py:305,316` | the general `eig` solver on a symmetric matrix, columns taken by position. |
| **F9** | `hicmatrix/lib/h5.py:84` | `distance_counts` read from `/correction_factors`. |
| **F10** | `hicmatrix/HiCMatrix.py:58-59` | `correction_factors` and `distance_counts` unpacked swapped. |
| **F11** | `hicmatrix/HiCMatrix.py:902-919` | `truncTrans` is a no-op. |
| **F12** | `bin/hicFindEnrichedContacts` | imports a module that does not exist (reproduced). |
| **F13** | `lib/tadClassifier.py:782` | `imblearn.base` used without import. |
| **F14** | `hicAdjustMatrix.py:161-165` | `--maskBadRegions` never opens its BED file (reproduced). |
| **F15** | `hicConvertFormat --enforce_integer` | a corrected matrix becomes all-zero pixels at exit 0 (reproduced). |
| **F16** | `hicMergeTADbins` | loses 22 percent of the total count. |
| **F17** | `hicMergeMatrixBins.py:248` | groups below `numBins/2` silently dropped. |
| **F18** | `hicSumMatrices.py:72` | the nan-bin union discards real counts. |
| **F19** | `--runningWindow` | crosses chromosome borders. |
| **F20** | `hicmatrix/HiCMatrix.py:84-92` | `save` reuses the load handler, so the output format follows the input (reproduced). |
| **F21** | F10 on cool input | a balanced cooler is written back without weights and with balanced counts. |
| **F22** | `hicFindTADs.py:345` | `--numberOfProcessors 12` or more aborts. |
| **F23** | `test_hicDetectLoops.py` | the tests call no loop; the loop order depends on `--threads`. |
| **F24** | `fit_nbinom` | not reproducible against itself on hicDetectLoops' distributions. |
| **F25** | `hicPCA.py:253-258`, `HiCMatrix.py:603` | `--ignoreMaskedBins` is undone by `--chromosomes`, and the restored matrix is float64 (reproduced). |
| **F26** | h5 against cool | cool loses per-bin coverage and the nan-bin list; `--sequencedCountCutoff` crashes on cool (reproduced). |
| **F27** | HDF5 output | object modification times embedded. |
| **F28** | `hicAdjustMatrix.py:124,144` | the remove warning names the last BED line's chromosome. |
| **F29** | `buildMatrixMethods.py:744-762` | `--keepSelfCircles` is a no-op; self circles double counted with two enzymes. |
| **F30** | `hicBuildMatrix --region` | empties the restriction site list. |
| **F31** | numpy and scipy on hicPCA's covariance | `eig` columns beyond the second change with the BLAS thread count. |
| **F32** | original `test_hicPCA.py` | the bigwig comparison was sign-agnostic and integer-rounded. |
| **F33** | hicDifferentialTAD, hicInterIntraTAD | output depends on the process count in three situations. |
| **F34** | hicInterIntraTAD | cannot process any h5 matrix: division by zero, exit 1 without output. |
| **F35** | hicInterIntraTAD, hicDifferentialTAD | the cool right block drops a bin; a two-TAD chromosome's second TAD gets no left test. |
| **F36** | chicQualityControl | accepts a reference point when any matrix passes, contrary to its help text. |
| **F37** | hicPlotSVL | an inf or nan ratio shifts later values onto the wrong chromosome; the checked-in `data.txt` shows 7 values for 12 names. |
| **F38** | `hicPlotSVL.py:183-196` | a worker exception leaves the parent polling forever. |
| **F39** | hicMergeDomains | imports `linkage` and `dendrogram` without calling them; the tree is DOT source rendered by `dot`, and only the first chromosome's graph is `strict`. |
| **F40** | `test_data/hicQuickQC` | the checked-in tables carry an older dangling-end label. |
| **F41** | the baseline suite | its four failures are two environmental and two test defects (reproduced). |
| **F42** | scipy 1.14 `betainc` | resolves to Boost.Math `ibeta`, not Cephes; Cephes differs by up to one ulp. |
| **F43** | chicQualityControl | a blank or short input line shifts later output rows onto the wrong sparsity and ends in `IndexError`. |
| **F44** | chicViewpoint | looks the background model up by bin distance in a model keyed by genomic distance, so every bin but the reference point is tested against the distribution at plus or minus `fixateRange`. |
| **F45** | chicViewpointBackgroundModel | drops the last bin of every viewpoint; a multi-bin reference point shifts its downstream values. |
| **F46** | chicViewpoint | a zero sum over the fixate range gives NaN relative values. |
| **F47** | chicViewpoint | a later reference point can be written into the chromosome group created last. |
| **F48** | chicViewpoint | a gene name used twice on one chromosome hangs the Python forever. |
| **F49** | hicMergeDomains | with the 50 kb file first, one TAD is written twice under the same ID; a start coordinate of 0 in a domain file's first 20 lines makes the Python loop forever when a protein file is given. |
| **F50** | chicViewpointBackgroundModel | its fit is chaotic on the test data: reordering a distribution moves the fitted size by a median of 55 percent, and the reference fails EN against itself. |
| **F51** | hicPlotSVL | float32 values print as the float64 repr of the float32 result. |
| **F52** | hicAggregateContacts | `--spectral` is parsed but never used. |
| **F53** | hicAggregateContacts | with several clusters the contact-pair file puts each cluster's values next to coordinates taken by position from the whole list. |
| **F54** | hicAggregateContacts | a minus-minus strand pair is transposed rather than rotated. |
| **F55** | hicAggregateContacts | the outlier filter divides by the median absolute value instead of the median absolute deviation, so it only removes anything with `--howToCluster center`. |
| **F56** | hicAggregateContacts | row-wise mode writes coordinates unordered. |
| **F57** | hicAggregateContacts | with `--perChr`, k stays 1 for every later chromosome once one has fewer submatrices than clusters. |
| **F58** | hicAggregateContacts | `--chromosomes` on a matrix with NaN bins undoes the masking: 277 submatrices instead of 279 on `Li_et_al_2015`. |

## Deliberate deviations from the Python behaviour

| Tool / component | Deviation | Decided by |
|---|---|---|
| acceptance, all tools | three significant digits per item, not byte identity | project owner, 2026-09-01 |
| implementation, all tools | data structures and algorithms free wherever the result is the same | project owner, 2026-09-13 |
| provenance, all tools | output is to name HiCExplorer 4 (recorded, not implemented) | project owner, 2026-09-02 |
| figures, all tools | explicit requests refused, default-named figures skipped, required figures behind a C++-only `--noPlot` (rule 7) | orchestrating session, pending the plotting decision |
| KR | no `exit(0)`; deterministic; class EN; reimplemented, not vendored | orchestrating session |
| hangs in the reference | hicPlotSVL worker exception (F38), chicViewpoint duplicate gene (F48), hicMergeDomains start coordinate 0 (F49): the port exits 1 with a message | orchestrating session |
| hicMergeDomains | an unknown `-of` format or a missing `dot` is refused before anything is written, where the Python writes text files first and then crashes | implementing agent, reported |
| hicAggregateContacts | an empty cluster exits 1 before writing, where the Python writes earlier clusters and then crashes | implementing agent, reported |
| chicViewpointBackgroundModel | fitted parameters judged by per-distribution likelihood and downstream calls, not EN | orchestrating session |
| hicDetectLoops, hicDifferentialTAD, hicInterIntraTAD | output independent of the thread count (F23, F33) | implementing agents, reported |
| hicFindTADs | no abort at 12 or more processes (F22) | implementing agent, reported |
| hicConvertFormat | `.hic` input and `--chromosome` refuse explicitly | implementing agent, reported |
| hicFindEnrichedContacts | not ported (F12) | plan |

## Open questions

1. Tiers 7 and 8, including the QC figures of hicBuildMatrix and hicQuickQC:
   confirm or replace the recommended strategy (project owner).
2. coolercpp's licence: GPL-3 like HiCExplorer, or BSD-3 like cooler, from whose
   source some of its components derive (project owner, before publication).
