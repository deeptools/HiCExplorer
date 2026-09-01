# HiCExplorer v4 (C++) - per-tool progress ledger

Owner: supervising agent. Architecture and rationale: `cpp/PLAN.md`. Rules:
`cpp/AGENTS_CONTRACT.md`. Initialised 2026-09-01. Current state: **1 of 46 tools
ported** (`hicInfo`, commit f1ced036), the tier 0 core partially landed, the
equivalence harness running for tier 1. Everything else is not started.

## Legend

- **Tier**: porting tier from `PLAN.md` section 6. Tier 0 (the core library)
  has no tools and is tracked in the second table below.
- **Class**: declared equivalence class from `PLAN.md` section 5.1. **The gate
  every tool must clear is ED: every item within three significant digits,
  relative `1e-3`, set by the project owner 2026-09-01 (`PLAN.md` 5.0).**
  Byte identicality is not required. A stricter class in this column means the
  tool is expected to do better than the gate and is checked at that level,
  because a stricter result is a better regression signal and several tools
  reach E0 for free; a drop from a declared strict class to ED is worth
  investigating even though it still passes. KR is the one tool that cannot be
  held to ED, because the Python reference disagrees with itself by `8.4e-03`
  on `gm12878_chr1.cool`; it stays at EN.
  E0 byte-identical, E1 HDF5-structural, E2 value-exact, E3 tight float
  (1e-12 rel), E4 loose float (1e-6 rel), E5 set agreement (Jaccard >= 0.99),
  E6 image (RMS <= 5), EN within the measured oracle-noise envelope,
  E7 not equivalent by design.
- **Mem**: peak-RSS budget as `alpha`/`beta` in
  `budget = alpha*W + beta*D + C`, with `W` the stored (upper-triangle) CSR
  working set, `D` the largest dense per-chromosome block, `C = 64 MB`
  (`PLAN.md` 4.5). Evaluated numbers for the memory-heavy tools are in the
  budgets table below. The budget is a **hard gate** in the harness
  (`PLAN.md` 8.3 criterion 4), not a report line.
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

| Tool | Tier | Class | Mem | Py test | Char. test | Port | Equiv | Notes |
|---|---|---|---|---|---|---|---|---|
| hicInfo | 1 | E0 | 0W (cool) / 1.05W (h5) | none -> char. test written | **done** | **done** | **pass (E0)** | Ported and byte-identical to the Python tool on all 184 cool and h5 matrices in `test_data`, commit f1ced036. `hicexplorer/test/general/test_hicInfo.py` was written first and passing against Python before any C++ existed (contract rule 1). Peak RSS on `gm12878_chr1.cool --no_metadata`: **809 MB against Python's 2,854 MB**, inside the 843 MB budget. Two further quirks uncovered by the characterization test and now pinned: `--outFileName` is reopened in `w` mode per matrix, so with several `-m` only the last block survives (`test_hicInfo.py:197`); and a cool file whose `bin-size` attribute is the string `"null"` prints no `Bin_length` line (`test_hicInfo.py:128`, which is the case for `Li_et_al_2015.cool`, `bin-type` variable). It also reproduces the cool-versus-h5 nnz split of `PLAN.md` 2.7 quirk 7. Remaining: mcool inputs in both group layouts are not yet in the case file. |
| hicConvertFormat | 1 | E1 cool / E2 h5 / E0 text | 1.3W | partial | yes | not started | - | 9 general tests at `assert_array_almost_equal(decimal=0)` plus 92 `trivial_runs` smoke cases. ginteractions (`:76`), hicpro (`:84`) and mcool (`:94`) conversions assert nothing. `:38-49` asserts an integer dtype without passing `--enforce_integer`. `hic` input is deferred to tier 4 (`PLAN.md` 3.6) and must error explicitly until then. |
| hicSumMatrices | 1 | E1/E2 | 2.2W | none | yes | not started | - | No test file. `Li_et_al_2015.h5` + `Li_et_al_2015_twice.h5` are the obvious char. test inputs. Exercises the `chrBinBoundaries` order check (`hicSumMatrices.py:52`) and `maskBins` union of nan bins. Two matrices live at once, hence alpha 2.2. |
| hicCompareMatrices | 1 | E1/E2 | 2.2W | partial | no | not started | - | 2 tests with exact `nt.assert_equal`, but only `--operation diff` and `--noNorm`. `ratio` and `log2ratio` untested. |
| hicAdjustMatrix | 1 | E1/E2 | 1.3W | partial | yes | not started | - | 4 real tests at `decimal=5`; the other 24 collected items are "did not crash". `:157,177` pass the BED file positionally after `--chromosomes`, so `--regions` is not exercised there. `--maskBadRegions` never tested. Takes the both-triangles exemption (`PLAN.md` 4.4 rule 2) because `reorderBins` permutes rows and columns independently. |
| hicMergeMatrixBins | 1 | E1/E2 | 1.3W | partial | yes | not started | - | 1 test, exact. `--runningWindow` (`hicMergeMatrixBins.py:88`) never tested. Also used internally by `hicConvertFormat` for mcool, so a regression here breaks tier 1. |
| hicFindRestSite | 2 | E0 | n/a (FASTA stream) | good | no | not started | - | 4 tests, `are_files_equal(delta=1)`, all options covered. Shells out to GNU `sort`; the C++ version sorts in memory under `LC_ALL=C` byte order. Budget is `2 * sites * 24 B + C`. |
| hicMergeLoops | 2 | E0 | n/a (BED) | partial | no | not started | - | 1 test, set-based comparison with `delta=2`. All options covered. |
| hicValidateLocations | 2 | E0 | 1.3W | good | no | not started | - | 4 tests, set-based `delta=1`. All options covered. Its `GSM1436265_..._10kb.cool` input (313,762 bins, 93 contigs, 7,987 nnz) is the many-contig stress case for `BinIndex` and the case where a per-bin `std::string` would dominate the footprint. |
| hicCreateThresholdFile | 2 | E0 | n/a (text) | partial | no | not started | - | 1 test. `--resolution/-r` never tested. Trivial tool. |
| hicMergeTADbins | 2 | E1/E2 | 1.3W | none | yes | not started | - | No test file. Exercises `reduceMatrix.reduce_matrix` with `diagonal=True` and deliberately clears `correction_factors` (`hicMergeTADbins.py:84`). |
| hicAverageRegions | 2 | E2 | 1.3W | partial | yes | not started | - | 8 tests at `decimal=0` on the `.npz` data array. `--considerStrandDirection` never tested. Output is `scipy.sparse.save_npz`, so tier 0 needs the `.npy`/`.npz` writer (`PLAN.md` 2.4). |
| hicNormalize | 2 | E1/E2 | 1.2W | good | no | not started | - | 8 tests with exact `nt.assert_equal` including h5-cool cross-format checks. `--setToZeroThreshold` never tested. The strongest test file in the suite. Pure in-place elementwise scaling, hence alpha 1.2. |
| hicTransform | 3 | E2 obs_exp / E3 pearson, covariance | 1.2W (obs_exp) / 1.1W + 1.15D (pearson) | weak | yes | not started | - | 11 tests but all at `assert_array_almost_equal(decimal=0)`, i.e. agreement to the nearest integer, which cannot detect a real numeric regression. `--chromosomes` never tested. Char. test must re-assert at full float precision. The pearson path is the second-largest memory win in the plan: 5,150 MB measured on `Li_et_al_2015.h5`, 259x its 19.9 MB working set. |
| hicCorrectMatrix | 3 | E3 ICE / **EN** KR | 1.2W | weak | yes | not started | - | **Dual-mode** (`--compatMode v3\|v4`), see the modes table. ICE/h5 is exact in the Python test, but KR/cool is only a range check `3e9 < sum//2 < 3688003604` (`:84`, elementwise comparison commented out at `:85-86`) and KR/partial uses `assert_allclose(rtol=1.0)` (`:106`). Given the measured nondeterminism (finding F1) a range test may be the only thing that could have passed reliably. Untested: `--perchr` (a distinct code path in both ICE and KR, and the site of finding F4), `--inflationCutoff`, `--transCutoff`, `--sequencedCountCutoff`, `--skipDiagonal`, `--xMax`, `--verbose`. **The largest memory win in the plan**: 9,199 MB (KR) and 9,064 MB (ICE) measured on `gm12878_chr1.cool`, against a 741.9 MB working set. |
| hicPCA | 3 | E3 + sign rule | 1.1W + 3.2D (v3) / 1.1W + 1.15D (v4) | weak | yes | not started | - | **Dual-mode**, see the modes table. 10 tests, but the bigwig comparison is `assert_array_almost_equal(np.absolute(...), decimal=0)` (`:64`): sign-agnostic and integer-rounded. `--histonMarkType` never tested. Highest-risk tool in tier 3: `scipy.linalg.eig` (not `eigh`) with unsorted output (`PLAN.md` 5.4). Also the worst memory-to-data ratio in the corpus: 4,070 MB and 471 s on a 722 KB input, 727x its 5.6 MB working set. |
| hicCompartmentalization | 3 | E3 | 1.3W | weak | yes | not started | - | 1 test, image-only at tolerance 60, and it is `xfail`. `--outputMatrix` and `--offset` never tested. Effectively uncovered. |
| hicInterIntraTAD | 3 | E3 | 1.3W | weak | yes | not started | - | `are_files_equal` at `:55` is called without `assert`, so its result is discarded; the only asserted check is an `xfail`-ed image comparison. Effectively "did not crash". |
| hicPlotSVL | 3 | E3 data / E6 plot | 1.3W | partial | yes | not started | - | 1 test, `are_files_equal(delta=2)` on the two text outputs; the image comparison is commented out at `:66-67`. Untested: `--distance`, `--chromosomes`, `--threads`, `--colorList`. Plot step follows the tier 7 rule. |
| hicBuildMatrix | 4 | E1/E2 matrix, E0 QC | `2*nnz_out*12 + threads*64 MB + C` | partial | yes | not started | - | 10 tests with exact matrix comparison, but the output BAM is only checked by byte size within 80,000 (`:17`). Untested: `--maxDistance`, `--keepSelfLigation`, `--doTestRunLines`. `trivial_runs` parametrizes `region="ChrX"` with the comment `# region does not work!!` and never passes it. Largest single unit of work in the port. The Python forks workers that each hold a share of the pixel buffers; the C++ pool shares one (`PLAN.md` 4.4 rule 8). |
| hicBuildMatrixMicroC | 4 | E1/E2 matrix, E0 QC | as hicBuildMatrix | weak | yes | not started | - | 1 test. 10 of 16 options untested: `--maxLibraryInsertSize`, `--genomeAssembly`, `--region`, `--keepSelfCircles`, `--minMappingQuality`, `--inputBufferSize`, `--doTestRun`, `--doTestRunLines`, `--skipDuplicationCheck`, `--chromosomeSizes`. `test_hicBuildMatrixMicroC.py:59` references a `delta` that is not defined in that module. |
| hicQuickQC | 4 | E0 | 1.05W | good | no | not started | - | 1 test with strict line equality on `QC.log`. All options covered. Shares `createMatrix` with hicBuildMatrix, so it comes free once that lands. |
| hicFindTADs | 5 | E3 scores / E5 calls | 2.2W | partial | yes | not started | - | 4 tests: exact on the z-score matrix, but the BED/GFF/BM/bedgraph outputs use `are_files_equal(pDifference=10)`, a per-line budget of 10 differing characters. `--TAD_sep_score_prefix` and `--delta` never tested. Two of the four tests pre-copy the reference z-score matrix, so they only test the downstream calling step. Uses `multiprocessing.Pool` (`hicFindTADs.py:1107`), unlike every other threaded tool. Takes the both-triangles exemption; revisit the budget once it works (`PLAN.md` risk 7). |
| hicDetectLoops | 5 | E3 stats / E5 calls | 2.2W | partial | yes | not started | - | 4 tests, `are_files_equal(delta=0)` on the loop bedgraph, but `test_main_h5` (`:44`) asserts nothing. `--obsExpThreshold` and `--expected` (the entire obs/exp preselection path) never tested. Depends on `fit_nbinom` (E4 on the fitted parameters, `PLAN.md` 3.5). Takes the both-triangles exemption. |
| hicDifferentialTAD | 5 | E5 | 1.3W x 2 | good | no | not started | - | 16 tests, `are_files_equal(delta=0, skip=4)`, i.e. exact, across all four `-m` modes and both `-mr` modes and thread counts 1/4/11. `--pValue` never varied from its default. Best-covered tool in the suite. |
| hicMergeDomains | 5 | E0/E5 | n/a (BED) | weak | yes | not started | - | Every `are_files_equal` call (`:70,84-85,106-107`) is missing `assert`; three of four tests are `xfail`. Effectively "did not crash". Untested: `--minimumNumberOfPeaks`, `--value`, `--percent`. Needs a `scipy.cluster.hierarchy.linkage` reimplementation and a DOT writer. |
| hicAggregateContacts | 5 | E3 `.tab` / E6 plot | 2.2W | weak | yes | not started | - | All 12 general tests are both `xfail(ImageComparisonFailure)` and `skipif(4 GB > memory)`, so nothing runs on a small machine and nothing asserts on a large one. The 162 `trivial_runs` items assert nothing, and `_three.py:64-67` hard-codes its arguments so its 72 parametrized cases are identical. Untested: `--considerStrandDirection`, `--largeRegionsOperation`, `--outFileObsExp`, `--spectral`, `--max_deviation`. Needs a `KMeans(random_state=0)` reimplementation (`PLAN.md` 3.5). The 4 GB `skipif` should become unnecessary once the memory rules land. |
| chicQualityControl | 6 | E0 text / E6 plot | 1.3W | weak | yes | not started | - | 1 test, `xfail`. `:66` compares the `_failed_reference_points` output against the `_report` reference, so that output is never actually checked. Untested: `--fixateRange`, `--dpi`. |
| chicViewpointBackgroundModel | 6 | E4 (NB parameters) | 1.3W | weak | yes | not started | - | `are_files_equal` allows 700 (`:62`) and 1000 (`:74`) mismatching values at `eps=0.1`. Untested: `--averageContactBin`, `--fixateRange`. Depends on `fit_nbinom`; E4 is declared for exactly this reason. |
| chicViewpoint | 6 | E1 hdf5 / E3 values | 1.3W | weak | yes | not started | - | 2 tests asserting only HDF5 keys and the default values of `averageContactBin` and `fixateRange`. No numeric comparison at all. `--averageContactBin` and `--fixateRange` never actually passed. |
| chicSignificantInteractions | 6 | E1 hdf5 / E5 calls | 1.3W | weak | yes | not started | - | 3 tests asserting HDF5 structure and attribute echo-back only. Untested: `--truncateZeroPvalues`, `--fixateRange`, `--peakInteractionsThreshold`. Uses `pybedtools` at `:515`. |
| chicAggregateStatistic | 6 | E1 hdf5 | 1.3W | weak | yes | not started | - | 5 tests asserting group names and `len(...)` only. All 5 options are at least passed. Uses `intervaltree` for target overlap. |
| chicDifferentialTest | 6 | E3 p-values / E5 calls | 1.3W | weak | yes | not started | - | 2 tests asserting HDF5 structure and `attrs['alpha']`/`attrs['test']`. No p-value comparison. Needs `fisher_exact`, `chi2_contingency` and `chi2.ppf` reimplementations. |
| chicExportData | 6 | E0 text / E1 bigwig | 1.3W | partial | yes | not started | - | 16 tests, 5 of them `xfail` by design (bad input). Text at `are_files_equal(delta=1, skip=1)`, bigwig at `decimal=0`. Untested: `--decimalPlaces`, `--oneTargetFile`, `--outputValueBigwig` (all three of its branches at `chicExportData.py:158-162` are dead code), `--threads`. |
| hicPlotMatrix | 7 | E6 | 1.1W + 1.15D | weak | yes | not started | - | 32 tests, **all 32 `xfail(ImageComparisonFailure)` and all 32 `skipif` on memory**, six of them behind a 120 GB gate. Nothing in this file has run in CI. 12 of 33 options untested: `--scoreName`, `--perChromosome`, `--vMin`, `--flipBigwigSign`, `--scaleFactorBigwig`, `--fontsize`, `--rotationX`, `--rotationY`, `--increaseFigureWidth`, `--increaseFigureHeight`, `--loops`, `--loopLargeRegionsOperation`. `:444` passes `--log1`, which argparse prefix-matches to `--log1p`. Python plotting shell over a C++ compute core; the 120 GB gate is a direct consequence of the memory blowups and should fall with them. |
| hicPlotTADs | 7 | E7 | n/a | none | no | not started | - | A 9-line delegation to `pygenometracks.plotTracks.main`. Nothing to port; stays a Python script. Recorded as a deliberate non-port. No characterization test is meaningful because the behaviour is entirely pyGenomeTracks'. |
| hicPlotViewpoint | 7 | E0 data / E6 plot | 1.3W | weak | yes | not started | - | All 6 tests `xfail`. `:83-97` passes `-i viewpoint_interactons` (a literal string) while asserting on the temp file names, so those assertions can only pass against stale files. `--chromosome` never tested. |
| hicPlotAverageRegions | 7 | E6 | n/a (npz) | weak | yes | not started | - | All 4 tests `xfail`, image-only. `--dpi` never tested. Reads the `.npz` written by hicAverageRegions and uses `scipy.ndimage.rotate`. |
| hicPlotDistVsCounts | 7 | E3 data / E6 plot | 1.3W | weak | yes | not started | - | The only real assertion is a PNG byte-size difference below 2000 (`:19`); the other 32 collected items are "did not crash". `--skipDiagonal` is parametrized and named but never inserted into the argument string (`:41-61`), doubling the case count for nothing. `--domains` never tested. |
| hicCorrelate | 7 | E3 matrix / E6 plot | 1.3W x n | weak | yes | not started | - | Both tests correlate one file with itself, so a correlation of 1.0 is structurally guaranteed. Untested: `--zMin`, `--zMax`, `--range`, `--threads`, and `--method pearson` entirely. Needs complete-linkage clustering with matching leaf order (`PLAN.md` 3.5). Budget scales with the number of input matrices. |
| hicPrepareQCreport (alias hicQC) | 7 | E0 tables / E6 plots | n/a (text) | none | yes | not started | - | No test file, although `test_data/QC*/` holds six reference output directories produced by it. Table aggregation moves to C++ at E0; the five bar charts and the Jinja2 HTML stay Python. |
| chicPlotViewpoint | 7 | E6 | 1.3W | weak | yes | not started | - | All 4 tests `xfail(ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')`. 9 options untested: `--outputFormat`, `--dpi`, `--colorMapPvalue`, `--maxPValue`, `--minPValue`, `--pValueSignificanceLevels`, `--xFold`, `--truncateZeroPvalues`, `--colorList`. |
| hicTADClassifier | 8 | E7 | 1.3W (features) | weak | yes | not started | - | 6 tests, of which 3 exercise the CLI and assert only `domain_df['Chrom'].iloc[0] == 1`. `--chromosomes` never tested. **Deliberate non-port**: the shipped `.BIN` models are pickles of `imblearn.EasyEnsembleClassifier` and `cleanlab.CleanLearning`, for which no ONNX converter exists (`PLAN.md` tier 8). Feature extraction moves to C++, which is also the memory win here; inference stays Python. |
| hicTrainTADClassifier | 8 | E7 | 1.3W (features) | weak | yes | not started | - | 1 test function, 4 sequential CLI runs, asserting only that the first word of the report is `accuracy`. **14 of 24 options untested**, the worst ratio in the suite: `--threshold`, `--leniency`, `--unselect_border_cases`, `--protein_file`, `--threads`, `--chromosomes`, `--concatenate_before_resample`, `--resampling_method`, `--alternative_resampling_method`, `--distance`, `--impute_value`, `--alternative_classifier`, `--use_cleanlab`, `--chrPrefixProtein`. **Deliberate non-port**, training stays Python. Note a latent `NameError` at `lib/tadClassifier.py:782` (`imblearn.base` referenced without `import imblearn`) that fires whenever `--alternative_resampling_method` is used. |
| hicHyperoptDetectLoops | 8 | E7 | inherits hicDetectLoops | partial | yes | not started | - | 3 tests, `are_files_equal(delta=2)`. `--resolution/-re` and `--threads` never used. **Partial port**: the driver moves to C++, TPE is deferred to v4.1 behind a `--parameterFile` mode; the harness instead checks that a fixed parameter set reproduces the Python result. |
| hicHyperoptDetectLoopsHiCCUPS | 8 | E7 | n/a (shells to java) | weak | yes | not started | - | 1 test, `xfail`, and doubly `skipif` on `nvcc` and on `juicer.jar` existing in the CWD, so it never runs. Its `are_files_equal` is `return True` (`:29`) and is called without `assert` (`:64`). Untested: `--chrPrefixLoops`, `--threads`, `--restricted`. **Partial port**: it shells out to `java -jar juicer.jar hiccups`; only the driver is ported. |

Counts: 46 tools. Tier 1: 6. Tier 2: 7. Tier 3: 6. Tier 4: 3. Tier 5: 5.
Tier 6: 7. Tier 7: 8. Tier 8: 4. Plus the `hicQC` alias and the `hicexplorer`
banner script, neither of which is a tool.

Test state at the start of the port: `none` 5, `weak` 22, `partial` 12,
`good` 7. Characterization tests required before porting: **39 of 46**, of which
**1 is written** (`hicInfo`) and 38 remain.

Ported: **1 of 46** (`hicInfo`, tier 1, equivalence pass at E0 and inside its
memory budget).

## Memory budgets on the designated large inputs

`budget = alpha*W + beta*D + C`, `C = 64 MB` (`PLAN.md` 4.5). Python peaks are
measured on this machine against the reference oracle with
`/usr/bin/time -f "%e %M"`. All figures in SI MB. The harness fails a tool whose
C++ peak RSS exceeds the budget (`PLAN.md` 8.3 criterion 4).

| tool and input | `W` | `D` | Python peak | budget | required reduction | measured C++ |
|---|---|---|---|---|---|---|
| `hicCorrectMatrix --correctionMethod KR`, gm12878_chr1.cool | 741.9 | - | 9,199 | **954** | 9.6x | - |
| `hicCorrectMatrix --correctionMethod ICE`, gm12878_chr1.cool | 741.9 | - | 9,064 | **954** | 9.5x | - |
| `hicPCA` v3 (dgeev), mm9_reduced_chr1.cool | 5.6 | 762 | 4,070 | **2,509** | 1.6x | - |
| `hicPCA` v4 (dsyevr), mm9_reduced_chr1.cool | 5.6 | 762 | 4,070 | **946** | 4.3x | - |
| `hicTransform --method pearson`, Li_et_al_2015.h5 | 19.9 | 987 | 5,150 | **1,221** | 4.2x | - |
| `hicTransform --method pearson`, gm12878_chr1.cool | 741.9 | 4,971 | > 25,000 (est.) | **6,596** | > 3.8x | - |
| `hicCorrectMatrix --correctionMethod KR`, Li_et_al_2015.h5 | 19.9 | - | 486 | **88** | 5.5x | - |
| `hicTransform --method obs_exp`, Li_et_al_2015.h5 | 19.9 | - | 305 | **88** | 3.5x | - |
| `hicConvertFormat` h5 -> cool, Li_et_al_2015.h5 | 19.9 | - | 390 | **90** | 4.3x | - |
| `hicInfo`, Li_et_al_2015.h5 | 19.9 | - | 264 | **85** | 3.1x | to be recorded |
| `hicInfo`, Li_et_al_2015.cool | 19.9 | - | 139 | **64** | 2.2x | to be recorded |
| `hicInfo --no_metadata`, gm12878_chr1.cool | 741.9 | - | 2,854 | **843** | 3.4x | **809 (pass)** |
| `hicBuildMatrix`, small_test_R1/R2_unsorted.bam, 16 threads | - | - | not yet measured | `2*nnz_out*12 + 16*64 + 64` | - | - |
| `hicFindTADs`, gm12878_chr1.cool | 741.9 | - | not yet measured | 1,696 | - | - |
| `hicDetectLoops`, GSE63525_..._2_5mb.cool | 12.5 (est.) | - | not yet measured | 92 | - | - |

Every remaining tool takes `alpha = 1.3, beta = 0` and is expected to land well
inside its budget; those rows are added to this table as they are measured.

The `hicBuildMatrix`, `hicFindTADs` and `hicDetectLoops` Python peaks must be
measured before their tiers start, so that the budget is a target rather than a
guess. That measurement is part of the characterization work for those tools.

## Dual-mode tools

`PLAN.md` 5.8. `v3` is the default and is what the harness compares against
Python; `v4` applies the fix and is validated against `v3` with the difference
quantified here.

| tool | `v3` reproduces | `v4` fixes | measured `v3` - `v4` difference |
|---|---|---|---|
| `hicCorrectMatrix --correctionMethod KR` | float32 downcast of input values (`krbalancing.cpp:27`); float32 rescale accumulators (`:228-229`), but summed in a fixed order so the C++ is deterministic where the Python is not; the `--perchr` `.h5`-vs-`.cool` correction-factor split (finding F4) | float64 throughout; pairwise accumulation. The `--perchr` `.h5`-vs-`.cool` behaviour is **preserved in both modes**: it is intended, not a defect, and unifying it would change what a corrected matrix means per output format | not yet measured; expected to be of the order of the oracle's own noise envelope, 1e-4 to 1e-2 relative depending on matrix size |
| `hicPCA` | `dgeev` on the covariance matrix, all eigenpairs, unsorted, columns taken by index (`hicPCA.py:305,314-322`) | `dsyevr` with `range='I'`, only the requested vectors, sorted by descending eigenvalue, deterministic sign convention | not yet measured; report per chromosome on `mm9_reduced_chr1.cool` |

Neither mode reproduces krbalancing's `exit(0)` on non-convergence; see the
deviations table.

## Known open work in the landed code

| where | issue | effect |
|---|---|---|
| `core/include/hicx/sparse_matrix.hpp` | values are held as `double` regardless of the stored dtype | an int32 or float32 matrix costs twice the memory it needs. On `gm12878_chr1.cool` a dtype-parametric value store would take `W` from 741.9 MB to 494 MB and every dependent budget with it. Recorded by the implementing agent in the header itself |
| `cpp/scripts/equiv.py` | peak RSS is recorded and reported as a Python-to-C++ ratio, but the **budget gate** of `PLAN.md` 8.3 criterion 4 is not implemented | a tool can currently be marked `pass` while over budget. Until the gate lands, budget compliance is checked by hand and written into the budgets table above |
| `cpp/scripts/equiv.py` | the determinism check of `PLAN.md` 8.3 criterion 3 (five repeats, `--threads 1` versus `--threads 16`) is not implemented | needed before any threaded tool can pass |
| `cpp/scripts/comparators/` | only `cool`, `h5` and `text` exist | `image`, `chic_hdf5`, `bigwig`, `npz`, `bam` and the `noise.py` EN driver are needed by tiers 3 and above |
| `cpp/scripts/cases/hicInfo.json` | no mcool cases | both group layouts (`/resolutions/<r>` and the legacy `/0`) need a case, see open question 5 |

## Findings against the Python reference

Behaviour found while planning that the port must account for. Most are defects
worth reporting upstream; F4 is not a defect but intended behaviour, kept here
because the port has to reproduce it deliberately. F1 changes what the port can
promise; F5 is a silent-success failure mode.

| id | where | finding |
|---|---|---|
| **F1** | `krbalancing.cpp:228-229` with `:233-250` | `rescale_norm_vector` accumulates `original_sum` and `norm_vector_sum` in **float32** inside an `omp parallel for` whose body is an `omp critical`. The critical section serialises the additions but does not fix their order, so **KR is not reproducible run to run.** Measured: six runs on `Li_et_al_2015.h5` gave six distinct normalisation factors (0.0190883 to 0.0190899), with a maximum pairwise relative difference of **1.503e-04 on the output matrix values** and **7.513e-05 on the correction factors**; three runs on `gm12878_chr1.cool` gave 0.00660838, 0.00663618 and 0.00666417, a spread of **8.4e-03 relative**. The sparsity pattern is stable. |
| **F2** | `krbalancing.cpp:12,27` | The float64 input is downcast to **float32** (`Eigen::Triplet<float>`, `float(input_values(j_start))`) before being stored into a float64 sparse matrix. Inert for raw integer counts below 2^24, but real for already-float matrices such as `Li_et_al_2015.h5` (values 0.170 to 1914.015). |
| **F3** | `krbalancing.cpp:11-32,35,47` | Construction stages the whole matrix into a `std::vector<Eigen::Triplet<float>>` before `setFromTriplets`, having already reserved the destination; `setFromTriplets` builds a third copy; `A = A + I` builds a fourth; and `triplets.clear()` at `:35` does not release capacity, so the staging buffer is still resident during that copy. On `gm12878_chr1.cool` that is 1,483 + 1,978 + 1,978 + 1,978 MB of avoidable allocation against a 742 MB working set. The input is already CSR and the matrix is symmetric, so CSR equals CSC and an `Eigen::Map` would need none of it. |
| **F4** | `hicCorrectMatrix.py:715-732` | **Not a defect. Intended behaviour, confirmed by the project owner 2026-09-01, and the port must reproduce it rather than unify the two branches.** The output format determines what a corrected matrix means. HiCExplorer's h5 historically could not carry correction factors, so for `.h5` output the correction is applied to the matrix values themselves and the file is written already corrected. cool carries the factors as a weight column, so for `.cool` output the raw matrix is written and the factors are stored beside it to be applied on read. That is why `get_normalised_matrix(True)` is called only for `.h5` (`:726`). Mechanically, that call is also what triggers `rescale_norm_vector()`, so the following `get_normalisation_vector(False)` (`:731`) returns rescaled factors for `.h5` and unrescaled ones for `.cool`. Two consequences for the port: the C++ must reproduce this **call ordering**, not just the arithmetic, because the rescaling is an in-place side effect on `x`; and the split is specific to `--perchr`. On the whole-matrix path (`:754`) `get_normalisation_vector(True)` is called unconditionally, so the factors are rescaled for both output formats there and only the matrix values differ. Note also that the modern h5 schema does carry a `correction_factors` dataset (present in `Li_et_al_2015.h5`, shape 11104, float64) and `ma.setCorrectionFactors()` is called for both formats; the historical limitation explains the pre-corrected matrix, not an absence of the dataset. The Python suite never runs `--perchr`, so the characterization test must cover all four combinations of `--perchr` and output format. |
| **F5** | `krbalancing.cpp:115-119` | `outer_loop` calls `exit(0)` after 300 outer iterations, having printed the entire `x` vector to stdout (also at 100 and 200). A library that terminates the host process with a **success** status on non-convergence: HiCExplorer exits 0 and writes no output file, which is indistinguishable from success to any calling pipeline. Worth reporting upstream as a bug in its own right. |
| **F6** | `krbalancing.cpp:212-221`, `krbalancing.hpp:29` | The loop body of `compute_normalised_matrix` is wrapped in `omp critical` although it only performs an elementwise in-place update with no shared state: pure serialisation with contention and no correctness role. `num_threads` is a hardcoded global of 10 that no caller can set. |
| **F7** | `hicCorrectMatrix.py:722,724` | `.indices.astype(np.int64, copy=False)` and `.data.astype(np.float64, copy=False)` both copy, because the dtypes differ. `copy=False` permits avoiding a copy, it does not achieve one. 989 MB each on `gm12878_chr1.cool`. |
| **F8** | `hicPCA.py:305` | `scipy.linalg.eig`, the general non-symmetric solver, is used on a symmetric covariance matrix and all n eigenpairs are computed and returned complex when two real ones are requested. 471 s and 4,070 MB on a 722 KB input. |
| **F9** | `hicmatrix/lib/h5.py:84` | `distance_counts` is read from `f.root.correction_factors`. |
| **F10** | `hicmatrix/HiCMatrix.py:58-59` | The loader tuple is unpacked with `correction_factors` and `distance_counts` swapped relative to what every loader returns, so an h5-to-h5 round trip moves correction factors into `/distance_counts`. |
| **F11** | `hicmatrix/HiCMatrix.py:902-919` | `truncTrans` is a no-op: a 2-tuple unpack of a 3-value return, and `==` where an assignment was meant. |
| **F12** | `bin/hicFindEnrichedContacts` | Installed by `setup.py` but imports `hicexplorer.hicFindEnrichedContacts`, which does not exist. |
| **F13** | `lib/tadClassifier.py:782` | `imblearn.base` is referenced without `import imblearn`; fires whenever `--alternative_resampling_method` is used. |

## Tier 0 - core library components

| Component | Path | Port | Notes |
|---|---|---|---|
| `CutIntervals`, `BinIndex` | `core/include/hicx/bins.hpp` | **done** | struct-of-arrays with interned chrom ids; replaces `intervaltree`; `PLAN.md` 2.1, 2.2 |
| `CsrMatrix`, `Matrix` | `core/include/hicx/sparse_matrix.hpp`, `hic_matrix.hpp` | **done** | **upper-triangle storage with symmetric access** (`PLAN.md` 4.4 rule 2), int32 indices while `nbins <= INT32_MAX` (rule 3) |
| HDF5 C wrapper + blosc filter | `core/src/hdf5_util.cpp` | **read path done** | the landed filter is **decompress-only**, which is enough to read every PyTables-written `.h5` in the corpus without an external filter plugin. The **compress** path, which is what `PLAN.md` risk 5 is about, is still unwritten and unverified |
| cool reader | `core/src/cool_file.cpp` | **done** | takes CSR row offsets straight from `/indexes/bin1_offset` and never reads `bin1_id`, so ingestion is zero-copy as rule 1 requires |
| cool writer | `core/src/cool_file.cpp` | not started | cooler schema v3, gzip-6 + shuffle, ENUM chrom column |
| mcool / scool | `core/src/io/cool.cpp` | not started | `::` path is opaque; discover only when absent; both `/resolutions/<r>` and the legacy `/0`../`/4` layouts occur in the corpus |
| h5 (PyTables layout) reader | `core/src/h5_file.cpp` | **done** | |
| h5 writer | `core/src/h5_file.cpp` | not started | blosc complevel 5; L3 value-identical only; blocked on the blosc compress path |
| JSON writer for the tier 0 state dump | `core/src/json_lite.cpp` | **done** | used by the harness comparators |
| homer, ginteractions, hicpro, 2D-text | `core/src/io/text_formats.cpp` | not started | ginteractions is write-only, as in Python |
| npy/npz | `core/src/io/npz.cpp` | not started | for `hicAverageRegions` / `hicPlotAverageRegions` |
| `.hic` reader | `core/src/io/hic.cpp` | not started | tier 4; until then `hicConvertFormat --inputFormat hic` must error |
| numpy-pairwise reduction | `core/src/numpy_compat.cpp` | **done** | all three layers of `PLAN.md` 5.2 (8192-element ufunc buffer, `PW_BLOCKSIZE = 128` with eight accumulators, float32 kept in float32), verified against numpy over 161 array sizes |
| Cephes (`gammaln`, `psi`, `betainc`, `igam`, `igamc`, `igami`) | `core/src/math/cephes/` | not started | vendored; makes NB, chi2 and Fisher exact |
| KR balancing | `core/src/math/kr/` | not started | reimplemented from the upstream source with a `v3`/`v4` mode switch; **not** vendored verbatim; `PLAN.md` 3.3, 5.7, 5.8 |
| ICE | `core/src/math/ice.cpp` | not started | port of `iterativeCorrection.py:10-86`; CSR row-order marginals reproduce `coo_matvec` exactly, so no mode switch |
| obs/exp, z-score, expected interactions | `core/src/math/obsexp.cpp` | not started | port of `utilities.py:293-604` |
| dense Pearson / covariance | `core/src/math/dense_corr.cpp` | not started | one dense block, `dsyrk`, streamed output (`PLAN.md` 4.4 rule 7) |
| `reduce_matrix` | `core/src/math/reduce_matrix.cpp` | not started | port of `reduceMatrix.py:12`; the complex-number `np.unique` trick becomes an in-place sort on a packed `(row,col)` key |
| `ranksums`, `anderson_ksamp`, `fisher_exact`, `chi2_contingency`, `pearsonr`, `spearmanr` | `core/src/math/stats.cpp` | not started | |
| complete-linkage clustering, k-means++ with `RandomState(0)` | `core/src/math/cluster.cpp` | not started | for `hicCorrelate`, `hicMergeDomains`, `hicAggregateContacts` |
| L-BFGS-B + NB MLE (`fit_nbinom`) | `core/src/math/nbinom_fit.cpp` | not started | `PLAN.md` risk 6 |
| BED / narrowPeak / broadPeak / bedgraph reader | `core/src/io/bed.cpp` | not started | port of `readBed.py` and `utilities.py:19,38` |
| bedtools replacement (sort, merge, intersect) | `core/src/util/intervals.cpp` | not started | must match bedtools lexicographic chrom order |
| libBigWig | `core/third_party/libBigWig/` | not started | FetchContent; read and write |
| htslib binding | `core/src/io/bam.cpp` | not started | link `$HICX_DEPS/lib/libhts.so.1.21` |
| FASTA reader + IUPAC revcomp | `core/src/io/fasta.cpp` | not started | replaces `Bio.SeqIO` |
| thread pool with deterministic reductions | `core/include/hicx/parallel.hpp` | not started | fixed index partitions, combine in index order, never accumulate into a shared float; `PLAN.md` 4.1 |
| peak-RSS self-report | `core/include/hicx/resource_usage.hpp` | **done** | `/proc/self/status` `VmHWM` at exit, printed on `--verbose`; cross-checked by the harness |
| argparse compatibility layer | `core/include/hicx/argparse.hpp` | not started | help text, groups, `choices`, `nargs`, prefix matching, error strings; `PLAN.md` tier 0 |
| Python-repr float formatting | `core/src/numpy_compat.cpp` | **done** | `float_repr`, `int_with_thousands_separator`, and numpy's `array_str` line wrapping at 75 characters |
| equivalence harness | `cpp/scripts/equiv.py` + `comparators/{base,cool,h5,text}.py` | **partly done** | `run`/`compare`/`report`/`list` work and both runtimes and peak RSS are recorded, with `--data` and `--tmpdir` added beyond the spec. Still missing: the **memory-budget gate** (`PLAN.md` 8.3 criterion 4; RSS is currently reported as a ratio only), the determinism check (criterion 3), `--noise-runs` and the `noise.py` EN driver, dual-mode runs, and the `image`, `chic_hdf5`, `bigwig`, `npz` and `bam` comparators |

## Deliberate deviations from the Python behaviour

None yet, because nothing is ported. Every entry here must name the tool, what
differs, why, and who decided. Contract rule 4: a tool that cannot be ported
faithfully is recorded here, never quietly dropped and never given a faked
equivalence.

Known deviations already planned (they become entries when the code lands):

| Tool / component | Planned deviation | Reason |
|---|---|---|
| KR balancing, both modes | `exit(0)` after 300 outer iterations (finding F5) is **not** reproduced; the port raises an error and exits non-zero | terminating the host process with a success status and no output file is not behaviour any correct pipeline can depend on |
| KR balancing, `v3` mode | the two float32 rescale accumulators are summed in a **fixed** order, so the C++ is deterministic where the Python is not (finding F1) | there is no reproducible target to match; the fixed order is one member of the family of orders the Python produces, and determinism is a hard requirement (`PLAN.md` 4.1) |
| `hicCorrectMatrix --correctionMethod KR` | equivalence class **EN**, not E2 or E3: the tolerance is measured from five Python runs rather than chosen | the oracle disagrees with itself by up to 1.5e-4 on values and 8.4e-3 on the normalisation factor |
| `hicCorrectMatrix`, `hicPCA` | ship a `--compatMode {v3,v4}` flag that does not exist in the Python | the accuracy and algorithm fixes change results; forcing them on users would break equivalence, dropping them would forfeit the fix |
| h5 writer | chunk shape and blosc block size differ from PyTables' | PyTables' chunk heuristic is undocumented and version-dependent; equivalence is declared at L3 value level (`PLAN.md` 2.5) |
| cool writer | `generated-by`, `generated-by-cooler-lib`, `tool-url`, `creation-date` will eventually say `hicx4` | provenance fields, normalised by the comparator. Until the suite is green the writer emits `HiCMatrix-17.2` verbatim so these fields are not a free pass |
| `hicTransform --method pearson`, `hicPCA --pearsonMatrix` | reduction order changes (one dense input block, streamed output, instead of five live copies) | fixes a 5,150 MB peak RSS on an 11,104-bin matrix; moves the class from E2 to E3 (`PLAN.md` 4.4) |
| all tools | upper-triangle storage with symmetric access instead of a materialised symmetric CSR | halves the resident matrix; no arithmetic changes provided the access helper preserves visit order, which has a unit test |
| `hicPlotTADs` and the 7 other tier 7 tools | not a C++ port; a Python plotting shell over the C++ core | matplotlib and pyGenomeTracks have no reproducible C++ equivalent (`PLAN.md` tier 7) |
| `hicTADClassifier`, `hicTrainTADClassifier` | not a C++ port; feature extraction in C++, model in Python | the shipped `.BIN` models are pickles of imblearn and cleanlab classes with no ONNX converter |
| `hicHyperoptDetectLoops`, `...HiCCUPS` | TPE search deferred; `--parameterFile` mode first | hyperopt's search path is RNG-dependent, so the chosen hyperparameters cannot be equivalent |
| `hicFindEnrichedContacts` | not ported | finding F12: the script is dead in the Python reference |

## Open questions for the implementing agent

1. (Partly resolved 2026-09-01.) The **decompress** half of the blosc HDF5
   filter is in and reads every PyTables-written `.h5` in the corpus. The
   question that `PLAN.md` risk 5 actually turns on is still open: does the
   **compress** half, at complevel 5 with shuffle, produce a file PyTables
   3.10.1 reads without complaint? Answer it before the h5 writer is built on
   top of it.
2. (Resolved 2026-09-01.) `C = 64 MB` holds. `hicInfo --no_metadata` on
   `gm12878_chr1.cool` peaks at **809 MB against a budget of 843 MB**
   (`1.05 * 741.9 + 64`), and against Python's 2,854 MB. The budget formula and
   the fixed allowance are therefore validated on the largest matrix in the
   corpus, on the tool with the least headroom to hide in. One known slack
   remains: `sparse_matrix.hpp` holds values as `double` regardless of the stored
   dtype, so an int32 matrix such as this one costs twice what it needs; with a
   dtype-parametric value store, `W` for this input falls from 741.9 MB to
   494 MB and the budget to 583 MB. That optimisation is recorded as open
   work, not as a budget relaxation.
3. Does calling LAPACK `dgeev` from `$HICX_DEPS`'s OpenBLAS reproduce
   `scipy.linalg.eig`'s eigenvector column order and signs on
   `hicPCA/mm9_reduced_chr1.cool`? Check this early; it decides whether
   `hicPCA` `v3` mode is possible at all (`PLAN.md` risk 3). Budget 8 minutes
   per Python reference run.
4. What is the actual `v3`-to-`v4` difference for KR on
   `gm12878_raw_values.cool` (raw int32, where F2 is inert) and on
   `Li_et_al_2015.h5` (float64, where it is not)? Fill it into the dual-mode
   table; it is the number that tells a user whether switching modes matters.
5. (Resolved 2026-09-01.) Both mcool group layouts occur in the corpus and both
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
6. Do `hicFindTADs`, `hicDetectLoops` and `hicAggregateContacts` really need
   both triangles materialised, or only a band around the diagonal? If a band
   suffices, withdraw their rule 2 exemption and tighten `alpha` from 2.2
   towards 1.3 (`PLAN.md` risk 7).
