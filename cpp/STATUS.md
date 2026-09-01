# HiCExplorer v4 (C++) - per-tool progress ledger

Owner: supervising agent. Architecture and rationale: `cpp/PLAN.md`. Rules:
`cpp/AGENTS_CONTRACT.md`. Initialised 2026-09-01 to the true current state:
nothing is ported.

## Legend

- **Tier**: porting tier from `PLAN.md` section 6. Tier 0 (the core library)
  has no tools and is tracked in the second table below.
- **Class**: declared equivalence class from `PLAN.md` section 5.1.
  E0 byte-identical, E1 HDF5-structural, E2 value-exact, E3 tight float
  (1e-12 rel), E4 loose float (1e-6 rel), E5 set agreement (Jaccard >= 0.99),
  E6 image (RMS <= 5), E7 not equivalent by design.
- **Py test**: state of the existing Python test.
  `none` = no test file; `weak` = a test file exists but its assertions do not
  constrain the numbers (structure-only, byte-size, unasserted comparison, or a
  tolerance so loose it cannot fail); `partial` = real value assertions on some
  paths only; `good` = real value assertions on the main paths.
- **Char. test**: is a characterization test required before the C++ port
  (contract rule 1). `yes` for every `none` and `weak` row, and for every
  `partial` row where the port touches an unexercised path.
- **Port**: `not started` / `in progress` / `done`.
- **Equiv**: `-` (not attempted) / `pass` / `fail` / `deviation`.

## Tools

| Tool | Tier | Class | Py test | Char. test | Port | Equiv | Notes |
|---|---|---|---|---|---|---|---|
| hicInfo | 1 | E0 | none | yes | not started | - | No test file. Invoked incidentally at `test_hicBuildMatrix.py:229`. Reports different nnz for the same matrix as cool (1,661,678) vs h5 (3,313,107), and a different set of fields per format; both must be reproduced (`PLAN.md` 2.7 quirk 7). Char. test must cover h5, cool, mcool, `--no_metadata`, `-o`, and multiple `-m`. |
| hicConvertFormat | 1 | E1 cool / E2 h5 / E0 text | partial | yes | not started | - | 9 general tests at `assert_array_almost_equal(decimal=0)` plus 92 `trivial_runs` smoke cases. ginteractions (`:76`), hicpro (`:84`) and mcool (`:94`) conversions assert nothing. `:38-49` asserts an integer dtype without passing `--enforce_integer`. `hic` input is deferred to tier 4 (`PLAN.md` 3.6) and must error explicitly until then. |
| hicSumMatrices | 1 | E1/E2 | none | yes | not started | - | No test file. `Li_et_al_2015.h5` + `Li_et_al_2015_twice.h5` are the obvious char. test inputs. Exercises the `chrBinBoundaries` order check (`hicSumMatrices.py:52`) and `maskBins` union of nan bins. |
| hicCompareMatrices | 1 | E1/E2 | partial | no | not started | - | 2 tests with exact `nt.assert_equal`, but only `--operation diff` and `--noNorm`. `ratio` and `log2ratio` untested. |
| hicAdjustMatrix | 1 | E1/E2 | partial | yes | not started | - | 4 real tests at `decimal=5`; the other 24 collected items are "did not crash". `:157,177` pass the BED file positionally after `--chromosomes`, so `--regions` is not exercised there. `--maskBadRegions` never tested. |
| hicMergeMatrixBins | 1 | E1/E2 | partial | yes | not started | - | 1 test, exact. `--runningWindow` (`hicMergeMatrixBins.py:88`) never tested. Also used internally by `hicConvertFormat` for mcool, so a regression here breaks tier 1. |
| hicFindRestSite | 2 | E0 | good | no | not started | - | 4 tests, `are_files_equal(delta=1)`, all options covered. Shells out to GNU `sort`; the C++ version sorts in memory under `LC_ALL=C` byte order. |
| hicMergeLoops | 2 | E0 | partial | no | not started | - | 1 test, set-based comparison with `delta=2`. All options covered. |
| hicValidateLocations | 2 | E0 | good | no | not started | - | 4 tests, set-based `delta=1`. All options covered. Its `GSM1436265_..._10kb.cool` input (313,762 bins, 93 contigs, 7,987 nnz) is the many-contig stress case for `BinIndex`. |
| hicCreateThresholdFile | 2 | E0 | partial | no | not started | - | 1 test. `--resolution/-r` never tested. Trivial tool. |
| hicMergeTADbins | 2 | E1/E2 | none | yes | not started | - | No test file. Exercises `reduceMatrix.reduce_matrix` with `diagonal=True` and deliberately clears `correction_factors` (`hicMergeTADbins.py:84`). |
| hicAverageRegions | 2 | E2 | partial | yes | not started | - | 8 tests at `decimal=0` on the `.npz` data array. `--considerStrandDirection` never tested. Output is `scipy.sparse.save_npz`, so tier 0 needs the `.npy`/`.npz` writer (`PLAN.md` 2.4). |
| hicNormalize | 2 | E1/E2 | good | no | not started | - | 8 tests with exact `nt.assert_equal` including h5-cool cross-format checks. `--setToZeroThreshold` never tested. The strongest test file in the suite. |
| hicTransform | 3 | E2 obs_exp / E3 pearson, covariance | weak | yes | not started | - | 11 tests but all at `assert_array_almost_equal(decimal=0)`, i.e. agreement to the nearest integer, which cannot detect a real numeric regression. `--chromosomes` never tested. Char. test must re-assert at full float precision. |
| hicCorrectMatrix | 3 | E3 ICE / E2 KR | weak | yes | not started | - | ICE/h5 is exact, but KR/cool is only a range check `3e9 < sum//2 < 3688003604` (`:84`, elementwise comparison commented out at `:85-86`) and KR/partial uses `assert_allclose(rtol=1.0)` (`:106`), a 100 % tolerance. Untested: `--perchr` (a distinct code path in both ICE and KR), `--inflationCutoff`, `--transCutoff`, `--sequencedCountCutoff`, `--skipDiagonal`, `--xMax`, `--verbose`. KR is expected E2 because the krbalancing C++ source is vendored (`PLAN.md` 3.3). |
| hicPCA | 3 | E3 + sign rule | weak | yes | not started | - | 10 tests, but the bigwig comparison is `assert_array_almost_equal(np.absolute(...), decimal=0)` (`:64`): sign-agnostic and integer-rounded. `--histonMarkType` never tested. Highest-risk tool in tier 3: `scipy.linalg.eig` (not `eigh`) with unsorted output (`PLAN.md` 5.4). |
| hicCompartmentalization | 3 | E3 | weak | yes | not started | - | 1 test, image-only at tolerance 60, and it is `xfail`. `--outputMatrix` and `--offset` never tested. Effectively uncovered. |
| hicInterIntraTAD | 3 | E3 | weak | yes | not started | - | `are_files_equal` at `:55` is called without `assert`, so its result is discarded; the only asserted check is an `xfail`-ed image comparison. Effectively "did not crash". |
| hicPlotSVL | 3 | E3 data / E6 plot | partial | yes | not started | - | 1 test, `are_files_equal(delta=2)` on the two text outputs; the image comparison is commented out at `:66-67`. Untested: `--distance`, `--chromosomes`, `--threads`, `--colorList`. Plot step follows the tier 7 rule. |
| hicBuildMatrix | 4 | E1/E2 matrix, E0 QC | partial | yes | not started | - | 10 tests with exact matrix comparison, but the output BAM is only checked by byte size within 80,000 (`:17`). Untested: `--maxDistance`, `--keepSelfLigation`, `--doTestRunLines`. `trivial_runs` parametrizes `region="ChrX"` with the comment `# region does not work!!` and never passes it. Largest single unit of work in the port. |
| hicBuildMatrixMicroC | 4 | E1/E2 matrix, E0 QC | weak | yes | not started | - | 1 test. 10 of 16 options untested: `--maxLibraryInsertSize`, `--genomeAssembly`, `--region`, `--keepSelfCircles`, `--minMappingQuality`, `--inputBufferSize`, `--doTestRun`, `--doTestRunLines`, `--skipDuplicationCheck`, `--chromosomeSizes`. `test_hicBuildMatrixMicroC.py:59` references a `delta` that is not defined in that module. |
| hicQuickQC | 4 | E0 | good | no | not started | - | 1 test with strict line equality on `QC.log`. All options covered. Shares `createMatrix` with hicBuildMatrix, so it comes free once that lands. |
| hicFindTADs | 5 | E3 scores / E5 calls | partial | yes | not started | - | 4 tests: exact on the z-score matrix, but the BED/GFF/BM/bedgraph outputs use `are_files_equal(pDifference=10)`, a per-line budget of 10 differing characters. `--TAD_sep_score_prefix` and `--delta` never tested. Two of the four tests pre-copy the reference z-score matrix, so they only test the downstream calling step. Uses `multiprocessing.Pool` (`hicFindTADs.py:1107`), unlike every other threaded tool. |
| hicDetectLoops | 5 | E3 stats / E5 calls | partial | yes | not started | - | 4 tests, `are_files_equal(delta=0)` on the loop bedgraph, but `test_main_h5` (`:44`) asserts nothing. `--obsExpThreshold` and `--expected` (the entire obs/exp preselection path) never tested. Depends on `fit_nbinom` (E4 on the fitted parameters, `PLAN.md` 3.5). |
| hicDifferentialTAD | 5 | E5 | good | no | not started | - | 16 tests, `are_files_equal(delta=0, skip=4)`, i.e. exact, across all four `-m` modes and both `-mr` modes and thread counts 1/4/11. `--pValue` never varied from its default. Best-covered tool in the suite. |
| hicMergeDomains | 5 | E0/E5 | weak | yes | not started | - | Every `are_files_equal` call (`:70,84-85,106-107`) is missing `assert`; three of four tests are `xfail`. Effectively "did not crash". Untested: `--minimumNumberOfPeaks`, `--value`, `--percent`. Needs a `scipy.cluster.hierarchy.linkage` reimplementation and a DOT writer. |
| hicAggregateContacts | 5 | E3 `.tab` / E6 plot | weak | yes | not started | - | All 12 general tests are both `xfail(ImageComparisonFailure)` and `skipif(4 GB > memory)`, so nothing runs on a small machine and nothing asserts on a large one. The 162 `trivial_runs` items assert nothing, and `_three.py:64-67` hard-codes its arguments so its 72 parametrized cases are identical. Untested: `--considerStrandDirection`, `--largeRegionsOperation`, `--outFileObsExp`, `--spectral`, `--max_deviation`. Needs a `KMeans(random_state=0)` reimplementation (`PLAN.md` 3.5). Plot step follows the tier 7 rule. |
| chicQualityControl | 6 | E0 text / E6 plot | weak | yes | not started | - | 1 test, `xfail`. `:66` compares the `_failed_reference_points` output against the `_report` reference, so that output is never actually checked. Untested: `--fixateRange`, `--dpi`. |
| chicViewpointBackgroundModel | 6 | E4 (NB parameters) | weak | yes | not started | - | `are_files_equal` allows 700 (`:62`) and 1000 (`:74`) mismatching values at `eps=0.1`. Untested: `--averageContactBin`, `--fixateRange`. Depends on `fit_nbinom`; E4 is declared for exactly this reason. |
| chicViewpoint | 6 | E1 hdf5 / E3 values | weak | yes | not started | - | 2 tests asserting only HDF5 keys and the default values of `averageContactBin` and `fixateRange`. No numeric comparison at all. `--averageContactBin` and `--fixateRange` never actually passed. |
| chicSignificantInteractions | 6 | E1 hdf5 / E5 calls | weak | yes | not started | - | 3 tests asserting HDF5 structure and attribute echo-back only. Untested: `--truncateZeroPvalues`, `--fixateRange`, `--peakInteractionsThreshold`. Uses `pybedtools` at `:515`. |
| chicAggregateStatistic | 6 | E1 hdf5 | weak | yes | not started | - | 5 tests asserting group names and `len(...)` only. All 5 options are at least passed. Uses `intervaltree` for target overlap. |
| chicDifferentialTest | 6 | E3 p-values / E5 calls | weak | yes | not started | - | 2 tests asserting HDF5 structure and `attrs['alpha']`/`attrs['test']`. No p-value comparison. Needs `fisher_exact`, `chi2_contingency` and `chi2.ppf` reimplementations. |
| chicExportData | 6 | E0 text / E1 bigwig | partial | yes | not started | - | 16 tests, 5 of them `xfail` by design (bad input). Text at `are_files_equal(delta=1, skip=1)`, bigwig at `decimal=0`. Untested: `--decimalPlaces`, `--oneTargetFile`, `--outputValueBigwig` (all three of its branches at `chicExportData.py:158-162` are dead code), `--threads`. |
| hicPlotMatrix | 7 | E6 | weak | yes | not started | - | 32 tests, **all 32 `xfail(ImageComparisonFailure)` and all 32 `skipif` on memory**, six of them behind a 120 GB gate. Nothing in this file has run in CI. 12 of 33 options untested: `--scoreName`, `--perChromosome`, `--vMin`, `--flipBigwigSign`, `--scaleFactorBigwig`, `--fontsize`, `--rotationX`, `--rotationY`, `--increaseFigureWidth`, `--increaseFigureHeight`, `--loops`, `--loopLargeRegionsOperation`. `:444` passes `--log1`, which argparse prefix-matches to `--log1p`. Python plotting shell over a C++ compute core. |
| hicPlotTADs | 7 | E7 | none | no | not started | - | A 9-line delegation to `pygenometracks.plotTracks.main`. Nothing to port; stays a Python script. Recorded as a deliberate non-port. No characterization test is meaningful because the behaviour is entirely pyGenomeTracks'. |
| hicPlotViewpoint | 7 | E0 data / E6 plot | weak | yes | not started | - | All 6 tests `xfail`. `:83-97` passes `-i viewpoint_interactons` (a literal string) while asserting on the temp file names, so those assertions can only pass against stale files. `--chromosome` never tested. |
| hicPlotAverageRegions | 7 | E6 | weak | yes | not started | - | All 4 tests `xfail`, image-only. `--dpi` never tested. Reads the `.npz` written by hicAverageRegions and uses `scipy.ndimage.rotate`. |
| hicPlotDistVsCounts | 7 | E3 data / E6 plot | weak | yes | not started | - | The only real assertion is a PNG byte-size difference below 2000 (`:19`); the other 32 collected items are "did not crash". `--skipDiagonal` is parametrized and named but never inserted into the argument string (`:41-61`), doubling the case count for nothing. `--domains` never tested. |
| hicCorrelate | 7 | E3 matrix / E6 plot | weak | yes | not started | - | Both tests correlate one file with itself, so a correlation of 1.0 is structurally guaranteed. Untested: `--zMin`, `--zMax`, `--range`, `--threads`, and `--method pearson` entirely. Needs complete-linkage clustering with matching leaf order (`PLAN.md` 3.5). |
| hicPrepareQCreport (alias hicQC) | 7 | E0 tables / E6 plots | none | yes | not started | - | No test file, although `test_data/QC*/` holds six reference output directories produced by it. Table aggregation moves to C++ at E0; the five bar charts and the Jinja2 HTML stay Python. |
| chicPlotViewpoint | 7 | E6 | weak | yes | not started | - | All 4 tests `xfail(ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')`. 9 options untested: `--outputFormat`, `--dpi`, `--colorMapPvalue`, `--maxPValue`, `--minPValue`, `--pValueSignificanceLevels`, `--xFold`, `--truncateZeroPvalues`, `--colorList`. |
| hicTADClassifier | 8 | E7 | weak | yes | not started | - | 6 tests, of which 3 exercise the CLI and assert only `domain_df['Chrom'].iloc[0] == 1`. `--chromosomes` never tested. **Deliberate non-port**: the shipped `.BIN` models are pickles of `imblearn.EasyEnsembleClassifier` and `cleanlab.CleanLearning`, for which no ONNX converter exists (`PLAN.md` tier 8). Feature extraction moves to C++; inference stays Python. |
| hicTrainTADClassifier | 8 | E7 | weak | yes | not started | - | 1 test function, 4 sequential CLI runs, asserting only that the first word of the report is `accuracy`. **14 of 24 options untested**, the worst ratio in the suite: `--threshold`, `--leniency`, `--unselect_border_cases`, `--protein_file`, `--threads`, `--chromosomes`, `--concatenate_before_resample`, `--resampling_method`, `--alternative_resampling_method`, `--distance`, `--impute_value`, `--alternative_classifier`, `--use_cleanlab`, `--chrPrefixProtein`. **Deliberate non-port**, training stays Python. Note a latent `NameError` at `lib/tadClassifier.py:782` (`imblearn.base` referenced without `import imblearn`) that fires whenever `--alternative_resampling_method` is used. |
| hicHyperoptDetectLoops | 8 | E7 | partial | yes | not started | - | 3 tests, `are_files_equal(delta=2)`. `--resolution/-re` and `--threads` never used. **Partial port**: the driver moves to C++, TPE is deferred to v4.1 behind a `--parameterFile` mode; the harness instead checks that a fixed parameter set reproduces the Python result. |
| hicHyperoptDetectLoopsHiCCUPS | 8 | E7 | weak | yes | not started | - | 1 test, `xfail`, and doubly `skipif` on `nvcc` and on `juicer.jar` existing in the CWD, so it never runs. Its `are_files_equal` (`:29`) is a real line-by-line comparator, but its result is discarded because the call at `:64` has no `assert`. Untested: `--chrPrefixLoops`, `--threads`, `--restricted`. **Partial port**: it shells out to `java -jar juicer.jar hiccups`; only the driver is ported. |

Counts: 46 tools. Tier 1: 6. Tier 2: 7. Tier 3: 6. Tier 4: 3. Tier 5: 5.
Tier 6: 7. Tier 7: 8. Tier 8: 4. Plus the `hicQC` alias and the `hicexplorer`
banner script, neither of which is a tool.

Test state: `none` 5, `weak` 22, `partial` 12, `good` 7.
Characterization tests required before porting: **39 of 46**.

## Tier 0 - core library components

| Component | Path | Port | Notes |
|---|---|---|---|
| `CutIntervals`, `BinIndex` | `core/include/hicx/bins.hpp` | not started | replaces `intervaltree`; `PLAN.md` 2.1, 2.2 |
| `CsrMatrix`, `Matrix` | `core/include/hicx/matrix.hpp` | not started | full-symmetric in memory, upper-triangular on disk; `PLAN.md` 2.3 |
| HDF5 C wrapper + blosc filter | `core/src/io/h5c.cpp`, `core/src/io/blosc_filter.c` | not started | vendored from hdf5-blosc; **verify first**, `PLAN.md` risk 4 |
| cool reader/writer | `core/src/io/cool.cpp` | not started | cooler schema v3, gzip-6 + shuffle, ENUM chrom column |
| mcool / scool | `core/src/io/cool.cpp` | not started | `::/resolutions/<r>` and `::/cells/<n>`; note `test_data/matrix.mcool` uses the legacy `/0`../`/4` group naming |
| h5 (PyTables layout) reader/writer | `core/src/io/h5file.cpp` | not started | blosc complevel 5; L3 value-identical only |
| homer, ginteractions, hicpro, 2D-text | `core/src/io/text_formats.cpp` | not started | ginteractions is write-only, as in Python |
| npy/npz | `core/src/io/npz.cpp` | not started | for `hicAverageRegions` / `hicPlotAverageRegions` |
| `.hic` reader | `core/src/io/hic.cpp` | not started | tier 4; until then `hicConvertFormat --inputFormat hic` must error |
| numpy-pairwise reduction | `core/include/hicx/reduce.hpp` | not started | required for E0 on `hicInfo`; `PLAN.md` 5.2 |
| Cephes (`gammaln`, `psi`, `betainc`, `igam`, `igamc`, `igami`) | `core/src/math/cephes/` | not started | vendored; makes NB and chi2 exact |
| krbalancing | `core/src/math/krbalancing/` | not started | vendored C++/Eigen upstream; `PLAN.md` 3.3 |
| ICE | `core/src/math/ice.cpp` | not started | port of `iterativeCorrection.py:10-86` |
| obs/exp, z-score, expected interactions | `core/src/math/obsexp.cpp` | not started | port of `utilities.py:293-604` |
| `reduce_matrix` | `core/src/math/reduce_matrix.cpp` | not started | port of `reduceMatrix.py:12`; the complex-number `np.unique` trick becomes a sort on a packed `(row,col)` key |
| `ranksums`, `anderson_ksamp`, `fisher_exact`, `chi2_contingency`, `pearsonr`, `spearmanr` | `core/src/math/stats.cpp` | not started | |
| complete-linkage clustering, k-means++ with `RandomState(0)` | `core/src/math/cluster.cpp` | not started | for `hicCorrelate`, `hicMergeDomains`, `hicAggregateContacts` |
| L-BFGS-B + NB MLE (`fit_nbinom`) | `core/src/math/nbinom_fit.cpp` | not started | `PLAN.md` risk 5 |
| BED / narrowPeak / broadPeak / bedgraph reader | `core/src/io/bed.cpp` | not started | port of `readBed.py` and `utilities.py:19,38` |
| bedtools replacement (sort, merge, intersect) | `core/src/util/intervals.cpp` | not started | must match bedtools lexicographic chrom order |
| libBigWig | `core/third_party/libBigWig/` | not started | FetchContent; read and write |
| htslib binding | `core/src/io/bam.cpp` | not started | link `$HICX_DEPS/lib/libhts.so.1.21` |
| FASTA reader + IUPAC revcomp | `core/src/io/fasta.cpp` | not started | replaces `Bio.SeqIO` |
| argparse compatibility layer | `core/include/hicx/argparse.hpp` | not started | help text, groups, `choices`, `nargs`, error strings; `PLAN.md` tier 0 |
| Python-repr float formatting | `core/include/hicx/fmt.hpp` | not started | shortest round-trip via `std::to_chars`; `PLAN.md` 5.6 |
| equivalence harness | `cpp/scripts/equiv.py` + `comparators/` | not started | spec in `PLAN.md` section 9 |

## Deliberate deviations from the Python behaviour

None yet. Every entry here must name the tool, what differs, why, and who
decided. Contract rule 4: a tool that cannot be ported faithfully is recorded
here, never quietly dropped and never given a faked equivalence.

Known deviations already planned (they become entries when the code lands):

| Tool / component | Planned deviation | Reason |
|---|---|---|
| h5 writer | Chunk shape and blosc block size will differ from PyTables' | PyTables' chunk heuristic is undocumented and version-dependent; equivalence is declared at L3 value level (`PLAN.md` 2.5) |
| cool writer | `generated-by`, `generated-by-cooler-lib`, `tool-url`, `creation-date` will eventually say `hicx4` | provenance fields; normalised by the comparator. Until the suite is green the writer emits `HiCMatrix-17.2` verbatim so these fields are not a free pass |
| `hicTransform --method pearson`, `hicPCA --pearsonMatrix` | Reduction order changes (single dense buffer instead of five) | fixes a 5.15 GB peak RSS on an 11,104-bin matrix; moves the class from E2 to E3 (`PLAN.md` 4.2) |
| `hicPlotTADs` and the 7 other tier 7 tools | Not a C++ port; a Python plotting shell over the C++ core | matplotlib and pyGenomeTracks have no reproducible C++ equivalent (`PLAN.md` tier 7) |
| `hicTADClassifier`, `hicTrainTADClassifier` | Not a C++ port; feature extraction in C++, model in Python | the shipped `.BIN` models are pickles of imblearn and cleanlab classes with no ONNX converter |
| `hicHyperoptDetectLoops`, `...HiCCUPS` | TPE search deferred; `--parameterFile` mode first | hyperopt's search path is RNG-dependent, so the chosen hyperparameters cannot be equivalent |
| `hicFindEnrichedContacts` | Not ported | `setup.py` installs `bin/hicFindEnrichedContacts`, but the module it imports does not exist. The script is dead in the Python reference |

## Open questions for the implementing agent

1. Does the vendored blosc HDF5 filter, at complevel 5 with shuffle, produce a
   file PyTables 3.10.1 reads without complaint? Verify before anything else in
   tier 0. If not, the whole h5 writer needs a different approach.
2. Does calling LAPACK `dgeev` from `$HICX_DEPS`'s OpenBLAS reproduce
   `scipy.linalg.eig`'s eigenvector column order and signs on
   `hicPCA/mm9_reduced_chr1.cool`? Check this early; it decides whether
   `hicPCA` is portable at all (`PLAN.md` risk 2).
3. Can the krbalancing upstream C++ source be fetched through the proxy? If
   not, KR drops from E2 to E3 and that must be recorded here.
4. (Resolved 2026-09-01.) Both mcool group layouts occur in the corpus and both
   must be readable. `hicConvertFormat --outputFormat mcool` writes
   `/resolutions/<res>`, verified by running it on
   `small_test_matrix_50kb_res.cool` with `-r 100000 200000`;
   `test_data/hicBuildMatrix/multi_small_test_matrix.mcool` is the same layout
   and is opened as `...mcool::/resolutions/5000` at
   `test_hicBuildMatrix.py:249-251`. But `test_data/matrix.mcool` has top-level
   groups `/0`../`/4` and no `/resolutions`. Since hicmatrix opens whatever URI
   it is handed and never enumerates resolutions itself, the C++ reader must do
   the same: treat the part after `::` as an opaque HDF5 group path, and only
   fall back to discovery (`/resolutions/*`, then top-level numeric groups) when
   no `::` suffix is given.
