# HiCExplorer v4 (C++) - per-tool progress ledger

Owner: the orchestrating session. Architecture: `cpp/PLAN.md`. Rules:
`cpp/AGENTS_CONTRACT.md`. Optimization rules: `cpp/OPTIMIZATION.md`. Last updated
2026-09-13, at commit `96727a7a`.

**Current state: 24 of 46 tools ported and committed** on `version4-cpp`. The
harness holds 301 cases for them plus the tier 0 round trips, and all 301 pass
with the memory and CPU-time gates live: 288 in a regression over every tool
except hicPCA on the merged tree of `96727a7a` (reproduced), and hicPCA's 13 on
its parent, 12 of them reproduced and mm9_reduced_chr1 in the agent's run only.
In progress, each in its own worktree under `~/src/HiCExplorer-v4-worktrees/`:
chicQualityControl, chicViewpointBackgroundModel and chicViewpoint (branch
`v4-chic-foundation`, complete and awaiting verification and merge);
hicMergeDomains and hicPlotSVL (`v4-mergedomains-svl`); hicAggregateContacts
(`v4-aggregate-contacts`).

**Provenance of the facts in this file.** "(reproduced)" marks a result the
orchestrating session reproduced itself. Unmarked results come from the
implementing agent's report and the tests and case files it committed. The
distinction has mattered repeatedly: a filter matching only lowercase `fail` once
hid three `FAIL` verdicts; commit 882164e5 claimed a pass count that held only
because most hicPCA cases had not run; and a verification that passed in the
orchestrating session's worktree hid three case inputs that were never committed.

## Legend

- **Tier**: porting tier from `PLAN.md` section 6.
- **Class**: equivalence class from `PLAN.md` 5.1. **Every tool must clear ED:
  each item within three significant digits, relative `1e-3` (`PLAN.md` 5.0).**
  Byte identity is not required, and how a tool reaches its result is free
  (`PLAN.md` 5.0 item 3). A stricter class here means the tool reaches it. E0
  byte identical, E1 HDF5 structural, E2 value exact, E5 set agreement (Jaccard
  `>= 0.99`), EN within the measured oracle-noise envelope, E7 not equivalent by
  design.
- **Mem**: peak-RSS budget `alpha*W + beta*D + C`, `C = 64 MB`, SI MB
  (`PLAN.md` 4.5). A hard gate in the harness.
- **Py test**: state of the Python test **before** the port.
- **Char. test**: `written` once a characterization test is committed; `yes`
  where one is still required; `no` where the existing test is adequate.
- **Equiv**: passing harness cases over total, gates live.

## Tools

| Tool | Tier | Class | Mem | Py test | Char. test | Port | Equiv | Notes |
|---|---|---|---|---|---|---|---|---|
| hicInfo | 1 | E0 | 0W cool / 1.05W h5 | none | written | done | 25/25 (reproduced) | Byte identical on all 184 cool and h5 matrices in `test_data`, reproduced on six. `--no_metadata` on `gm12878_chr1.cool`: 829.2 MB against an 842.9 MB budget, the tightest case in the corpus. |
| hicConvertFormat | 1 | E1 cool / E2 h5 / E0 text | 1.3W | partial | written | done | 17/17 (reproduced) | `.hic` input and `--chromosome` refuse explicitly. Reproduces F15. Open: h5 `extra_list` is float64 where Python infers int64 for hicpro and homer sources. |
| hicSumMatrices | 1 | E1/E2 | 3.0W | none | written | done | 11/11 (reproduced) | `alpha` 3.0: scipy binary operations drop exact zeros, so the peak is `nnz_a + nnz_b + nnz_result`. Reproduces F18. |
| hicCompareMatrices | 1 | E1/E2, ED on `log2ratio` | 3.0W | partial | written | done | 18/18 (reproduced) | numpy's own float64 `log2` differs from libm by at most `3.2e-16`. |
| hicAdjustMatrix | 1 | E1/E2 | 1.3W | partial | written | done | 21/21 (reproduced) | Upper triangle of `P A P^T` built directly, no both-triangles exemption. Reproduces F14, F28. |
| hicMergeMatrixBins | 1 | E1/E2, ED on one case | 1.3W, see open work | partial | written | done | 13/13 (reproduced) | Bit identical at `-nb 5` on `Li_et_al_2015.h5` (reproduced). Reproduces F17, F19. |
| hicFindRestSite | 2 | E0 | n/a | good | written | done | 7/7 (reproduced) | |
| hicMergeLoops | 2 | E0 | n/a | partial | written | done | 5/5 (reproduced) | |
| hicValidateLocations | 2 | E0 | 1.3W | good | written | done | 6/6 (reproduced) | |
| hicCreateThresholdFile | 2 | E0 | n/a | partial | written | done | 5/5 (reproduced) | |
| hicMergeTADbins | 2 | E1/E2 | 1.3W | none | written | done | 4/4 (reproduced) | Reproduces F16. |
| hicAverageRegions | 2 | E2 | 1.3W | partial | written | done | 12/12 (reproduced) | Three of its case inputs were missing from git until `b8552cb1`; see the provenance note above. |
| hicNormalize | 2 | E1/E2 | 1.2W | good | written | done | 9/9 (reproduced) | |
| hicTransform | 3 | ED or better | 1.2W / 1.1W + 1.15D | weak | written | done | 23/23 (reproduced) | `--method pearson` on `Li_et_al_2015.h5`: 5,532 to 97 MB, 64.9 to 3.9 s CPU (reproduced). |
| hicCorrectMatrix | 3 | ED ICE / EN KR | 1.2W | weak | written | done | 20/20 (reproduced) | krbalancing reimplemented. ICE and KR on `gm12878_chr1.cool` from about 9.3 GB to 830 MB (reproduced). `--compatMode v3` against the default: `1.65e-04` relative on `Li_et_al_2015.h5` (reproduced). `diagnostic_plot` refuses explicitly. |
| hicPCA | 3 | ED or better | 1.1W + 1.15D | weak | written | done | 13/13 (12 reproduced) | Fixed in `a3c49fd8`. The three failures present since 882164e5 are gone: (a) eigenvectors are chosen by position from `scipy.linalg.eig`, and on chrX the largest eigenvalue occurs 169 times, so the column order depends on the last bits of the covariance, which now equals `np.cov` bit for bit (F8, F31); (b) `--ignoreMaskedBins` with `--chromosomes` reproduces the Python undoing the mask and returning float64 (F25); (c) the whole-genome round-off case passes as a consequence of (a). mm9_reduced_chr1, over two hours of Python CPU, passed in the agent's run only. |
| hicCompartmentalization | 3 | E0 `_dat` / E2 npz | 1.3W | weak | written | done | 15/15 (reproduced) | Largest case 2,534 to 26 MB and 216.7 to 4.4 s CPU. The required `--outputFileName` is the figure, so a plain run refuses; a C++-only `--noPlot` writes the numeric outputs (refusal and `--noPlot` reproduced by hand). Deviation pending the plotting decision. |
| hicInterIntraTAD | 3 | E0 | 1.3W | weak | written | done | 5/5 (reproduced) | An explicitly requested ratio plot exits 1 before writing (reproduced by hand); a default-named one is skipped with a note. Reproduces F33 to F35. |
| hicPlotSVL | 3 | E0 text | 1.3W | partial | yes | in progress | - | Branch `v4-mergedomains-svl`. Reports F37 and F38. The figure follows rule 7. |
| hicBuildMatrix | 4 | E0 BAM and QC / E1 / E2 | `2*nnz_out*12 + threads*64 MB + C` | partial | written | done | 17/17 (reproduced) | Output BAM and QC tables byte identical. The QC folder lacks the Python's PNGs and HTML report (open work). Reproduces F29, F30. |
| hicBuildMatrixMicroC | 4 | as hicBuildMatrix | as hicBuildMatrix | weak | **yes, not written** | done | 4/4 (reproduced) | Landed without a characterization test, contrary to rule 1. |
| hicQuickQC | 4 | E0 | budget 192 MB | good | written | done | 5/5 (reproduced) | Reuses hicBuildMatrix's BAM code. One named harness normalisation for the random temporary matrix name. CPU 47.3 to 0.76 s, median RSS 921 to 92 MB. Reproduces F40. |
| hicFindTADs | 5 | E0 text / E1 / E2 | per case, `2.2 * Wz + C` | partial | written | done | 8/8 (reproduced) | F22. |
| hicDetectLoops | 5 | E0 text / EN fitted size / E5 calls | 2.2W | partial | written | done | 13/13 (reproduced) | F23, F24. Uses cephes `betainc`, one ulp away from scipy 1.14 (F42, open work). Peak RSS on `gm12878_chr1.cool` 877 MB against Python's 619 MB. |
| hicDifferentialTAD | 5 | E0, E5 calls | 1.3W x 2 | good | written | done | 12/12 (reproduced) | Jaccard 1.0 on every case. Largest case at 92.8 percent of budget. Reproduces F33. |
| hicMergeDomains | 5 | E0/E5 | n/a | weak | yes | in progress | - | Branch `v4-mergedomains-svl`. **Its reference is not broken**: the two tests that failed in the baseline are environmental (F41), and the tool writes byte identical output on all three paths. It needs no clustering (F39). Its tree is DOT source rendered by the external graphviz `dot`, which the port reproduces the same way. |
| hicAggregateContacts | 5 | E3 `.tab` / E6 plot | 2.2W | weak | yes | in progress | - | Branch `v4-aggregate-contacts`. All 12 general tests are `xfail` and `skipif(4 GB)`. Which clustering actually executes is to be established before building any. |
| chicQualityControl | 6 | E0 text | 1.3W | weak | written | in progress | - | Branch `v4-chic-foundation`: 5/5 reported, awaiting verification. Figures follow rule 7. |
| chicViewpointBackgroundModel | 6 | exact columns; fitted size and prob by likelihood and downstream E5 | 1.3W | weak | written | in progress | - | Branch `v4-chic-foundation`: 5/5 reported. The EN rule is withdrawn for `size` and `prob`, because the reference fails it against itself at the same rate (F50); they are judged by likelihood against the reference runs and by the downstream calls (Jaccard 1.000000 reported on all four fitted cases). A per-distribution tolerance refinement is in progress. |
| chicViewpoint | 6 | E1 / bit identical values | 1.3W | weak | written | in progress | - | Branch `v4-chic-foundation`: 6/6 reported, bit identical, which required scipy's Boost.Math `betainc` (F42). |
| chicSignificantInteractions | 6 | E5 calls | 1.3W | weak | yes | not started | - | Builds on the chicViewpoint core. |
| chicAggregateStatistic | 6 | E1 | 1.3W | weak | yes | not started | - | |
| chicDifferentialTest | 6 | E3 p-values / E5 calls | 1.3W | weak | yes | not started | - | Needs `fisher_exact`, `chi2_contingency`, `chi2.ppf`. |
| chicExportData | 6 | E0 text / E1 bigwig | 1.3W | partial | yes | not started | - | The three branches of `--outputValueBigwig` are dead code. |
| hicPlotMatrix | 7 | E6 | 1.1W + 1.15D | weak | yes | not started | - | All 32 tests `xfail` and `skipif` on memory. |
| hicPlotTADs | 7 | E7 | n/a | none | no | not started | - | A delegation to pyGenomeTracks. |
| hicPlotViewpoint | 7 | E0 data / E6 plot | 1.3W | weak | yes | not started | - | |
| hicPlotAverageRegions | 7 | E6 | n/a | weak | yes | not started | - | |
| hicPlotDistVsCounts | 7 | E3 data / E6 plot | 1.3W | weak | yes | not started | - | |
| hicCorrelate | 7 | E3 matrix / E6 plot | 1.3W x n | weak | yes | not started | - | Needs complete-linkage clustering with matching leaf order. |
| hicPrepareQCreport (alias hicQC) | 7 | E0 tables / E6 plots | n/a | none | yes | not started | - | Its input, the QC folder, is byte identical from hicBuildMatrix and hicQuickQC. |
| chicPlotViewpoint | 7 | E6 | 1.3W | weak | yes | not started | - | |
| hicTADClassifier | 8 | E7 | 1.3W | weak | yes | not started | - | Shipped models are pickles of imblearn and cleanlab classes. |
| hicTrainTADClassifier | 8 | E7 | 1.3W | weak | yes | not started | - | F13. |
| hicHyperoptDetectLoops | 8 | E7 | inherits | partial | yes | not started | - | |
| hicHyperoptDetectLoopsHiCCUPS | 8 | E7 | n/a | weak | yes | not started | - | Its comparator (`:29`) is sound, but the call at `:64` has no `assert`. |

**Tiers 7 and 8 await the project owner's decision.** Until then no C++ tool draws
figures: an explicitly requested figure is refused before any output is written,
and one the Python would write only under a default name is skipped with a note
(contract rule 7).

## Counts

| tier | tools | ported | in progress | not started |
|---|---|---|---|---|
| 1 | 6 | 6 | 0 | 0 |
| 2 | 7 | 7 | 0 | 0 |
| 3 | 6 | 5 | 1 | 0 |
| 4 | 3 | 3 | 0 | 0 |
| 5 | 5 | 3 | 2 | 0 |
| 6 | 7 | 0 | 3 | 4 |
| 7 | 8 | 0 | 0 | 8 |
| 8 | 4 | 0 | 0 | 4 |
| **total** | **46** | **24** | **6** | **16** |

## Baseline of the Python suite

First full run on this branch, 2026-09-01: 498 passed, 3 failed, 1 error, 1
skipped, 23 xfailed, 62 xpassed (reproduced). **None of the four was a defect in a
HiCExplorer tool** (causes reproduced from the run's tracebacks):

| entry | cause | kind |
|---|---|---|
| `test_hicMergeDomains::test_main_two_file_protein` | graphviz cannot find `dot`, because the venv's `PATH` lacks the conda env's `bin`; the generic exception does not match `xfail(raises=ImageComparisonFailure)` | environment |
| `test_hicMergeDomains::test_main_two_file_no_protein` | the same | environment |
| `test_hicBuildMatrix_trivial_runs::test_build_matrix_restrictionCutFile_two` | passes `region="ChrX"`, absent from the input; `getUserRegion` raises a bare `KeyError: 'ChrX'` | test defect |
| `test_hicBuildMatrix_trivial_runs_2::test_build_matrix_restrictionCutFile_six` (error) | takes `keepSelfLigation` without a matching `parametrize`, so pytest looks for a fixture | test defect |

## Memory on the designated large inputs

| tool and input | Python peak | budget | measured C++ |
|---|---|---|---|
| hicCorrectMatrix KR, gm12878_chr1.cool | 9,409 | 954 | 829 (reproduced) |
| hicCorrectMatrix ICE, gm12878_chr1.cool | 9,331 | 954 | 830 (reproduced) |
| hicPCA, mm9_reduced_chr1.cool | 4,187 | 2,509 | 1,589 |
| hicTransform pearson, Li_et_al_2015.h5 | 5,532 | 1,221 | 97 (reproduced) |
| hicCompartmentalization, small_test_matrix.h5 | 2,534 | - | 26 |
| round trip cool, gm12878_chr1.cool | 4,048 | 1,028 | 790 (reproduced) |
| hicInfo --no_metadata, gm12878_chr1.cool | 2,854 | 842.9 | 829.2 (reproduced) |
| hicBuildMatrix, small_test BAMs | 936 | 321 | 142 (reproduced) |
| hicQuickQC, median over cases | 921 | 192 | 92 |
| hicDetectLoops, gm12878_chr1.cool | 619 | 1,696 | 877, above Python |

## Known open work in the landed code

| where | issue | effect |
|---|---|---|
| all tools | `PLAN.md` 5.0.1 v4 provenance and `H5Pset_obj_track_times` recorded, not implemented | output still names HiCExplorer 3.7.x and embeds object modification times (F27); to land across all tools after the open branches merge |
| hicBuildMatrix, hicQuickQC | the QC folder lacks the Python's PNGs and `hicQC.html` | the plotting decision covers it |
| hicDetectLoops | cephes `betainc`, one ulp from scipy 1.14's Boost.Math `ibeta` (F42) | its 13 cases pass; it should share the Boost implementation once `v4-chic-foundation` merges |
| `cpp/tests/` | two refusal scripts, `refuses_unavailable_output.cmake` (generic) and `check_compartmentalization_refuses_plot.cmake` | consolidate into the generic one |
| harness `cpp_args` | the case field passes options to the C++ side only | allowed only for options that do not change compared outputs; currently used for `--noPlot` alone |
| `core/src/cool_file.cpp` | no hyperslab band read | hicDetectLoops uses more memory than Python on `gm12878_chr1.cool` |
| `core/include/hicx/sparse_matrix.hpp` | values held as double regardless of dtype | int32 and float32 matrices cost twice what they need |
| budget formula | `--runningWindow` and hicFindTADs do not fit `alpha * W_input` | both need a term in the working set they build |
| argument parsing | no argparse compatibility layer | prefix abbreviations rejected; `--help` close but not identical |
| process | an agent cut off by a usage limit left a mutation-test run consuming a CPU core and 4.4 GB for eleven days | killed 2026-09-13; wall-clock figures between 2026-09-02 and 2026-09-13 ran under that extra load, CPU-time figures are unaffected by design |

The round-off floor proposed for class ED is **not applied**: its only
motivating failure disappeared with hicPCA's exact covariance, and the gate is not
widened in advance.

## Findings against the Python reference

| id | where | finding |
|---|---|---|
| **F1** | `krbalancing.cpp:228-229` | KR is not reproducible run to run: float32 rescale sums in an unordered critical section (reproduced). |
| **F2** | `krbalancing.cpp:12,27` | float64 input downcast to float32. |
| **F3** | `krbalancing.cpp:11-47` | four full copies of the matrix during construction. |
| **F4** | `hicCorrectMatrix.py:715-732` | **Intended, confirmed by the project owner.** h5 output is written corrected, cool keeps raw counts with weights; the port reproduces the call ordering. |
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
| **F23** | `test_hicDetectLoops.py` | the tests call no loop; the Python's loop order depends on `--threads`. |
| **F24** | `fit_nbinom` | not reproducible against itself on hicDetectLoops' distributions. |
| **F25** | `hicPCA.py:253-258`, `HiCMatrix.py:603` | `--ignoreMaskedBins` is undone by `--chromosomes`, and the restored matrix is float64, shifting obs/exp values by `3.1e-02` (reproduced). |
| **F26** | h5 against cool | cool loses per-bin coverage and the nan-bin list; `--sequencedCountCutoff` crashes on cool (reproduced). |
| **F27** | HDF5 output | object modification times embedded. |
| **F28** | `hicAdjustMatrix.py:124,144` | the remove warning names the last BED line's chromosome. |
| **F29** | `buildMatrixMethods.py:744-762` | `--keepSelfCircles` is a no-op; self circles double counted with two enzymes. |
| **F30** | `hicBuildMatrix --region` | empties the restriction site list. |
| **F31** | numpy and scipy on hicPCA's covariance | `np.cov`'s bits, and `eig` columns 3 and 4 on chrX, change between one and 32 BLAS threads, so `--whichEigenvectors` beyond 2 depends on the Python's BLAS thread count. |
| **F32** | `hicexplorer/test/general/test_hicPCA.py` (original) | the bigwig comparison was sign-agnostic and integer-rounded. |
| **F33** | hicDifferentialTAD, hicInterIntraTAD | output depends on the process count in three situations (a small chromosome first drops the next chromosome's last TAD; a small one after a large one rewrites stale slots; one TAD per process skips the last left test). The C++ always gives the single-process result. |
| **F34** | hicInterIntraTAD | cannot process any h5 matrix: division by zero on the last TAD's left block, exit 1 without output. |
| **F35** | hicInterIntraTAD, hicDifferentialTAD | the cool right block ends at a literal `-1` and drops a bin; a two-TAD chromosome's second TAD gets no left test. |
| **F37** | hicPlotSVL | a chromosome whose ratio is inf or nan is dropped from the values but not the names, so later values are written against the wrong chromosome and trailing rows are blank; `test_data/hicPlotSVL/data.txt` shows 7 values for 12 names. |
| **F38** | `hicPlotSVL.py:183-196` | an exception inside a worker leaves the parent polling forever. |
| **F39** | hicMergeDomains | imports `linkage` and `dendrogram` but never calls them; the tree is DOT source rendered by the external `dot`, and only the first chromosome's Digraph is `strict`. |
| **F40** | `test_data/hicQuickQC` | the checked-in `QC_table.txt` and `discarded_table.txt` carry an older dangling-end label than the Python now writes. |
| **F41** | the baseline suite | its four failures are two environmental and two test defects; see the baseline table (reproduced). |
| **F42** | scipy 1.14 `betainc` | resolves to Boost.Math `ibeta`, not cephes; cephes differs by up to one ulp, enough to turn a stored p-value of 0.0 into `1.1e-16`. |
| **F50** | chicViewpointBackgroundModel | its fit is chaotic on the cHi-C test data: reordering a distribution moves the fitted size by a median of 55 percent, and the reference fails the EN rule against itself. |

F36 and F43 to F49 are unassigned; the remaining cHi-C findings take them when `v4-chic-foundation` is merged.

## Deliberate deviations from the Python behaviour

| Tool / component | Deviation | Decided by |
|---|---|---|
| acceptance, all tools | three significant digits per item, not byte identity | project owner, 2026-09-01 |
| implementation, all tools | data structures and algorithms free wherever the result is the same | project owner, 2026-09-13 |
| provenance, all tools | output is to name HiCExplorer 4 (recorded, not implemented) | project owner, 2026-09-02 |
| figures, all tools | explicit requests refused, default-named figures skipped (rule 7) | orchestrating session, pending the plotting decision |
| hicCompartmentalization | a C++-only `--noPlot` flag | orchestrating session, pending the plotting decision |
| KR | no `exit(0)`; deterministic; class EN; reimplemented, not vendored | orchestrating session |
| hicCorrectMatrix, hicPCA | `--compatMode`; hicPCA needs no second mode since its covariance is exact | orchestrating session |
| hicDetectLoops, hicDifferentialTAD, hicInterIntraTAD | output independent of the thread count (F23, F33) | implementing agents, reported |
| hicFindTADs | no abort at 12 or more processes (F22) | implementing agent, reported |
| hicPlotSVL | exits non-zero where the Python hangs on a worker exception (F38) | orchestrating session |
| chicViewpointBackgroundModel | fitted parameters judged by likelihood and downstream calls, not EN | orchestrating session, on merge |
| hicConvertFormat | `.hic` input and `--chromosome` refuse explicitly | implementing agent, reported |
| hicFindEnrichedContacts | not ported (F12) | plan |

## Open questions

1. Tiers 7 and 8: confirm or replace the recommended strategy (project owner).
2. hicAggregateContacts: can scikit-learn's `KMeans(random_state=0)` partition
   be reached without reproducing its random number generator?
3. The QC figures and `hicQC.html` of hicBuildMatrix and hicQuickQC: part of the
   plotting decision.
