# HiCExplorer v4 (C++) - per-tool progress ledger

Owner: the orchestrating session. Architecture: `cpp/PLAN.md`. Rules:
`cpp/AGENTS_CONTRACT.md`. Optimization rules: `cpp/OPTIMIZATION.md`. Last updated
2026-09-22, at commit `dab5fe2d`.

**Current state: 42 of 44 tools ported and committed** on `version4-cpp`
(`hicTADClassifier`/`hicTrainTADClassifier` dropped 2026-09-21, see PLAN.md
tier 8: little real-world usage, not worth an ONNX-incompatible Python
shell).
- **Last regression,** merge-style (`--cache refresh`, contract rule 13) on a
  clean export of the capture Hi-C branch `eea2a808` (reproduced); merged in
  `2452a010`.
  - 535 of 535 cases over 43 tools, with no time-gate reruns;
  - determinism for all 7 chic tools;
  - ctest 6, specs 43, argparse 79, gui 138 (1 skipped, snakemake), bindings 99;
  - the Python chic tests, including the new characterization tests: 30 passed, 5 xfailed.
  - The build showed two `-Wformat-truncation` warnings in the new tar writer:
    a member of 8 GiB or more would have got a truncated size. They were fixed
    in `a55555cd`, and a clean export of that commit has 0 warnings, ctest 6,
    and chicExportData 22 of 22 with determinism.
- **Earlier regression,** on clean exports of the plotting branch (`6d27fa8e`,
  then `cedc9f33` for its final fix), built out of tree with the Python module
  and without `HICX_PLOT_PYTHON` set (reproduced). Merged in `b71510ba`.
  - 470 cases over 37 tools and the round trips;
  - hicPCA 13 of 13, now including mm9_reduced_chr1;
  - determinism for all 15 tools the branch changed;
  - ctest 6 of 6;
  - specs 39, argparse 79, gui 117 (1 skipped, snakemake), bindings 99;
  - a refused drawing environment exits 3 and leaves the output directory
    unchanged.
- **Previous regression,** on a clean export of `f4c15dc9`, built out of tree with
  the Python module (reproduced). That commit merges the streamed `.hic` to cool
  and the GUI foundation (PLAN 10.1 to 10.3).
  - 384 cases over 29 tools and the tier 0 round trips;
  - hicPCA 12 of 12, all cases except mm9_reduced_chr1, which passed in an agent's run only;
  - determinism for hicConvertFormat, roundtrip, hicCorrectMatrix, hicDifferentialTAD, hicBuildMatrix, hicPlotSVL and hicAdjustMatrix;
  - tool specifications: 31 of 31 tests against the Python argparse parsers, and 79 of 79 command lines behaving as argparse does;
  - bindings 95 of 95; workflow engine 53 passed, 1 skipped (a Snakefile run, snakemake not installed), including the end-to-end workflow E0 against its logged command lines.
- The earlier merge `bfcf75fa` (cool read speed, hicDifferentialTAD, legacy
  `.hic`) passed the same regression, and `53a9fc91` moved both library pins to
  commits that differ only in the libraries' consumer tests.
- **Merged 2026-09-16, `7769ec55`:** `.pairs` input for hicBuildMatrix (PLAN
  9.2) and hicPCA's `--eigenSolver lanczos` plus a faster exact dense path
  (tier 11), verified merge-style together on a clean export at `9bea3c46`
  (reproduced): full regression 546 of 546 over 43 tools, determinism for
  hicPCA and the three matrix-building tools, ctest 7, all suites, and a
  GM12878 chr1 25 kb lanczos-against-dense spot check (max relative
  difference 0, 8.5 s at 765 MB against 221 s at 1,834 MB).
  - `.pairs`: E2 against `cooler cload pairs` and against the BAM route;
    136 M real ENCODE pairs at 29.5 s CPU / 2.6 GB against cooler's 148 s /
    3.65 GB.
  - hicPCA dense: bit-identical to the pre-change output; chromosomes now run
    concurrently and only the requested eigenvector columns are
    back-transformed. Multithreaded OpenBLAS was tried and dropped: every
    LAPACK stage changes bits with the thread count.
  - hicPCA lanczos: ED against dense wherever LAPACK's column order is
    unambiguous; order-differing cases are documented, not counted as passing.
- **Merged 2026-09-16, `967ff69b`:** the replicate-aware differential engine,
  new tool `hicDifferentialAnalysis {tads,loops,compartments}` (PLAN 9.7, work
  item 2). Verified merge-style from a clean export at `967ff69b` (reproduced):
  full regression 559 of 559 over 44 tools, determinism for hicPCA (15),
  hicDifferentialAnalysis (13) and hicBuildMatrix (26), ctest 7, cache pytest
  17, spec tests 44, argparse 79, gui suite 138 (4 legitimate skips).
  - **Model:** negative binomial GLM per unit, distance-decay and
    library-size offsets, quasi-likelihood dispersion with a robust
    (trimmed) empirical-Bayes prior (Phipson, Lee, Majewski, Alexander and
    Smyth 2016), Benjamini-Hochberg FDR, a TREAT-style minimum-fold-change
    test, a shared bin mask. Refuses single-replicate input unless
    `--exploratory`.
  - **TADs:** Simes-combined per-distance-stratum tests. **Loops:** tested
    against local background at the union of all samples' calls.
    **Compartments:** GC-oriented scores.
  - **Calibrated against the PLAN 9.7 gate** (fixed before implementation) on
    GSE234292 (wt/knockout, 2 replicates each): TADs, boundaries and loops
    meet every criterion (nulls call at most ~5 % at p<=0.05, 0 at FDR 0.05;
    2-fold plant recall >=0.917, observed FDR <=0.083). **Compartments do not
    meet the 0.8 recall gate at 2-fold** (0.771 unpaired, 0.719
    replicate-blocked), even after the robust prior improved it from 0.765
    and 0.655 with a plug-in prior. Not tuned to pass; recorded as an open
    limitation below.
  - The unpaired label-swap null is evaluated with genotype as block, since
    it cannot meet the gate unpaired by construction (a real wt-vs-knockout
    effect inflates its within-group variance).
  - Also fixed: `equiv.py`'s reference cache key now includes the data root
    for tool outputs that embed input paths (hicInfo, hicPlotSVL,
    chicQualityControl, hicValidateLocations), which previously gave false
    mismatches when the cache was built in a different checkout.
- **Merged 2026-09-18, `ab2cff32`:** new tool `hicDetectStripes` (PLAN 9.3),
  a faithful C++ port of Stripenn 1.1.65.22's own detection pipeline.
  Verified merge-style from a clean export at `ab2cff32` (reproduced): full
  regression 562 of 562 over 45 tools, ctest 8, check_case_inputs 0 missing,
  spec tests 45, argparse 79, gui suite 139, cache pytest 17.
  - **Method:** maxpixel-percentile image construction (whole-chromosome
    population, not band-limited), `skimage.feature.canny` reproduced pixel
    for pixel against the installed library (Gaussian smoothing with the real
    bleed-over border correction, separable Sobel gradients matching scipy's
    summation order, bilinear non-maximum suppression, real default
    magnitude thresholds 0.1/0.2, not quantile-based), the `verticalLine`
    Sobel-direction filter, zero-row/column frame compaction before the
    block-run-length scan, up/down column-pairing into stripe boxes,
    `RemoveRedundant` overlap-and-elongation dedup, and Stripenn's real 2D
    rectangular direction-specific background windows with a per-row median
    test. No multiple-testing correction by default, matching Stripenn's own
    real-world practice (raw p-value threshold); `--fdr` is an opt-in
    C++-only option.
  - **Fidelity checks against real Stripenn**, each independently
    established, not assumed from reading the source: per-candidate p-values
    match on identical real GM12878 windows (this port 0.03-0.13, Stripenn
    itself 0.042-0.19); raw candidate generation on the same chromosome
    matches Stripenn's own count (474 vs 456, after fixing two real gaps:
    the frame compaction above, and reading `maxpixel`'s percentiles from the
    whole chromosome instead of the near-diagonal band, which had made them
    2-2.7x too high).
  - **PLAN 9.3's recall (0.8) and precision (0.9) gate on the fixed 60-plant
    real GM12878 test is recorded as measured and not met.** Real,
    unmodified Stripenn on the identical plants recovers 4 of 60 (6.7 %);
    this port recovers 6 of 60 (10 %), matching or slightly exceeding the
    reference. Two further real, correctly-verified, independent paradigms
    were benchmarked against the same plants and both also failed at usable
    precision: Quagga (Gaussian-blur peak detection, Poisson/NB test) 0 of 60
    (0 %); Chromosight (cross-correlation against its own stripe kernels)
    looked like a pass at the gate bucket (16 of 20, 0.80) but was shown to
    be chance-level noise from massive over-calling (79,975 raw calls
    genome-wide against about 2,000 for Stripenn or this port; a control of
    200 random non-plant anchors under the identical matching rule gave a
    41.5 % "recovery" rate by chance alone).
  - **Decided by the project owner (2026-09-18):** ship the port with this
    result documented rather than continue searching for a passing method or
    revise the plant design. See PLAN.md 9.3 for the full evidence trail
    across all four methods.
  - Reported without a gate: GSE234292 rep1/rep2 Jaccard 0.065; agreement
    with Stripenn on the real unplanted genome, Jaccard 0.037-0.06 depending
    on which fix stage was measured.
  - Memory/CPU on 22 real autosomes: hicDetectStripes 383-645 MB, 7-91 s
    wall depending on the FDR/candidate-generation fix stage measured; real
    Stripenn 1099-1101 MB, 732-758 s wall, both the memory and CPU
    reference. hicDetectStripes is lighter and faster throughout.
  - Determinism confirmed at `-t1`/`-t4` (byte-identical, cool and h5).
  - Two further real defects found and fixed while running the multi-method
    comparison: the orientation-label convention was backwards and used
    exact bin equality instead of a distance comparison (measured effect on
    the isolated-geometry recall check: 0/6 to 1/6 on real GM12878
    chr21+chr22); `stripe_calibration.py` itself had a full-genome-table
    pixel load and a `getrusage(RUSAGE_CHILDREN)` memory measurement that
    did not reset between subprocess calls.
- **Merged 2026-09-18, `c8a173f3`:** GUI redesign, multi-project tabs and
  data-driven tool filtering (PLAN tier 10), requested directly by the
  project owner after reviewing the previous single-project layout.
  Verified merge-style from a clean export at `c8a173f3` (reproduced): build
  0 warnings, ctest, gui suite all pass.
  - `MainWindow`'s central widget is now a `QTabWidget` of `ProjectTab`s:
    each open project is its own closable top-level tab with independent
    state (tool browser, run controller and view, workflow editor, matrix
    browser). Closing a project (confirmed if a run is in progress) leaves
    other open projects untouched. The C++ tool directory stays a single
    window-global setting, shared by every open project.
  - `dataformats.py` detects cool/mcool/h5/`.hic`/BAM/`.pairs` by extension
    and, where ambiguous, real content signatures.
  - The per-project tool list filters to whatever accepts the loaded file's
    format, built from each tool's own `--help-json` file `role`/`formats`
    data, not a hand-maintained mapping (verified against the real build:
    loading a cool file filtered 45 of 49 tools down to the 29 that declare
    a cool-accepting input argument). Loading a FASTQ file shows an honest
    message that no tool reads it directly (PLAN tier 13 will close this
    with minibwa) instead of inventing support.
  - Settings is no longer a tab: opened from the File/View menu as a dialog.
  - Verified on the real desktop display, not only offscreen: two
    independently-stated project tabs, the cool-driven filter, and the
    settings dialog, confirmed by screenshot review.
- **Merged 2026-09-18, `7538037c`:** auto-open loaded matrix files in the
  Matrix browser, and an arc-style loop track, both requested directly by
  the project owner after using the redesigned GUI live.
  Verified merge-style from a clean export at `7538037c` (reproduced): build
  0 warnings, ctest 6/8 (2 pre-existing environment-gated skips), gui suite
  155 passed/5 skipped.
  - Loading a cool/mcool/h5/`.hic` file now calls the same `open_matrix()`
    path the manual "Open in Matrix browser" button already used, so it
    appears in the Matrix browser immediately without a second click. The
    manual button is unchanged. Non-matrix formats are untouched.
  - Loops are now drawn as a below-matrix arc track (`loop_arc_path` in
    `browser.py`), stacked with the existing bedgraph/bigwig tracks on the
    same x-axis, replacing the old point-marker overlay drawn on the matrix
    face. Arc apex height scales by the loop's BEDPE score column when
    present, else by span; only loops with both anchors in the current view
    are drawn.
  - Verified on the real desktop display, not only offscreen: auto-open
    after loading a cool file, and three loop arcs of visibly different
    heights matching their scores, confirmed by screenshot review.
- **Corrected 2026-09-21, `e2f0b103`:** two scoping mistakes from the
  orchestrating session, reverted at the project owner's direct instruction.
  - minibwa integration (merged 2026-09-18 in `4c5e79a5` as new C++ tools
    `hicBuildIndex`/`hicAlignReads`, PLAN tier 13) is reverted in full: both
    tools, the `minibwa_bridge` core module, and their catalog/spec entries
    are removed. minibwa was never meant to become a HiCExplorer-branded
    tool; it belongs to `hicexplorer-gui` instead, exposed as itself. See
    PLAN.md tier 13 for the corrected scope.
  - The GUI (`gui/`, all of tier 10's work) is removed from this repository
    entirely. It was split out, with its real commit history (`git subtree
    split`), into the standalone `~/src/hicexplorer-gui`, which depends on
    stock Python HiCExplorer (`>=3.7,<4`) rather than this C++ rewrite. This
    repository is the C++ tool rewrite only from here on.
  - Verified: clean rebuild (0 warnings, 0 errors), ctest 6/8 (2 pre-existing
    environment-gated skips), 45 real tool binaries remain (47 minus the two
    reverted; unaffected by the removal since neither was ever part of the
    original 46 Python-ported tools).
- **Merged 2026-09-22, `73f60939`:** native `.hic` reading wired into the
  shared `ToolMatrix::load`, closing a real gap the project owner caught:
  PLAN.md already claimed "every matrix-reading tool accepts a `.hic`," which
  was not true before this commit. Only `hicConvertFormat` could read a
  `.hic` at all; `hicPlotMatrix`, `hicInfo` and every other generic tool
  failed on one (confirmed: "file signature not found," since the shared
  loader tried to open it as HDF5).
  - Detection by real content signature (the `"HIC"` bytes hicfilecpp's own
    reader already requires), falling back to the `.hic` extension only when
    the signature check finds nothing.
  - Selector syntax follows the existing mcool convention exactly:
    `file.hic::/resolutions/10000`,
    `file.hic::/resolutions/10000/normalizations/KR`.
  - A region or chromosome query reads only the blocks it needs through
    hicfilecpp's own block index (`MatrixZoomData::getRecords`), the same
    genuinely partial read a cool `fetch()` already did, not a whole-file or
    whole-chromosome load. Verified on the real 39.9 GB GSE63525 GM12878
    file: `1:18000000-22000000` at 10 kb loads in well under a second at
    under 100 MB peak RSS (0.63 to 0.69 s, 91 to 96 MB, reproduced
    independently by the orchestrating session on a clean export), against
    several minutes and gigabytes for a whole-genome load, with a
    byte-identical PNG result either way.
  - New unit tests (`cpp/tests/test_hic_adapter.cpp`) cross-check the native
    path against `hic2cool_convert` + `read_cool`, class E2 (bit-identical
    stored values): whole file raw, whole file KR, bare-chromosome region
    raw, sub-region raw, sub-region KR. `hicx_tests` 296/296 (151,046
    assertions), ctest 6/8 (2 pre-existing skips), reproduced on a clean
    export by the orchestrating session.
  - Known follow-up, not done: `hicAdjustMatrix.cpp` has an analogous
    cool-only single-chromosome fast path that would need the same `.hic`
    extension; left for a later pass.
- **Merged 2026-09-22, `bea7c0e9`:** `hicPlotMatrix --matrix2`, a v4-only
  upper/lower triangle comparison heatmap (upper = `--matrix`, lower =
  `--matrix2`, diagonal always `--matrix`, reasoned in the code comment on
  `combine_triangles`), requested directly by the project owner. Both
  matrices load through the same genuinely partial, region-scoped path
  `73f60939` added, confirmed together on a 2 GB cool file and the 39.9 GB
  real `.hic` file: 0.80 s, ~52 MB peak RSS for both matrices combined.
  Element-wise correctness verified via `--plotData` against separate
  single-matrix dumps, pinned by a new test
  (`cpp/tests/cli/hicPlotMatrix_triangle_cli.sh`). `hicx_tests`/`ctest`/the
  existing 23 `hicPlotMatrix` equivalence cases all still pass, reproduced
  on a clean export by the orchestrating session.
- **Merged 2026-09-22, `9696f88a`:** `hicMergeMatrixBins` now requires
  `--chromosomeSizes`/`-cs`, fixing a real correctness bug the project owner
  flagged directly: the merged bin layout used to come from whichever bins
  the input happened to contain, not the true genome length, so two inputs
  of the same genome/resolution with different missing bins produced
  differently-shaped, non-comparable output (the Python original has this
  same bug, reproduced deliberately until now). New
  `merge_bins_genome`/`plan_bin_merge_genome` tile every chromosome from
  position 0 at the matrix's own resolution, independent of input
  presence; scoped to fixed-resolution matrices only (a restriction-
  fragment matrix falls back to the old, Python-faithful path).
  `hicConvertFormat`'s own use of the shared `merge_bins`/`reduce_matrix` is
  untouched (33/33 equivalence cases confirmed unaffected). `hicx_tests`
  299/299 (151,179 assertions), ctest 9/9, `hicMergeMatrixBins` equivalence
  14/14 (fixed-resolution cases reclassified E1/E2 to EX, a deliberate
  v4-only deviation, documented per-case), reproduced on a clean export by
  the orchestrating session.
- **Merged 2026-09-22, `5be73193` then `412c5d45`:** CHiCAGO scoring
  (PLAN.md 9.15), requested directly by the project owner, now a usable
  feature. `5be73193` landed the foundational per-interaction statistics
  (Delaporte p-value, distance-based weighting/score, distance-function
  fit, genome-geometry formulas), independently ED-verified against a real,
  privately provisioned R 4.5.3 + Chicago 1.38.0 + PCHiCdata 1.38.0
  environment (Chicago's real pipeline on PCHiCdata's real GM12878
  chr20/chr21 and mouse ES chr18/chr19 data), including a real precision
  bug found and fixed (naive linear-space convolution loses all precision
  on deep-tail p-values) and one honestly-flagged deviation (9 of ~4,600
  sampled rows where R's own compiled `pdelap_C` itself underflows and
  substitutes an approximation this implementation does not reproduce).
  `412c5d45` finished it: genome-wide parameter estimation (technical
  noise, bait/other-end normalisation, Brownian dispersion via
  `MASS::theta.ml`'s Newton iteration), and three new, standalone,
  C++-only tools, corrected mid-task from an earlier flag-based design to
  separate tools per the project owner's direct naming instruction:
  `chicChicagoBackgroundModel` (mirrors `chicViewpointBackgroundModel`'s
  role), `chicChicagoScores` (mirrors `chicViewpoint`'s role),
  `chicChicagoSignificantInteractions` (mirrors
  `chicSignificantInteractions`'s role, default score threshold 5).
  - Cross-checking against R found and fixed two more real bugs:
    `normaliseOtherEnds` sums the bait-normalised `NNb` count, not raw
    `N`; its per-pool `ntot` sums `nbpb` only over observed pairs, not
    every defined bin.
  - Dispersion 2.550454 exact match to R's own on GM12878; `s_j` max
    relative deviation 4.9e-12 across 648 baits; `s_i` and `Tmean` exact
    matches.
  - End-to-end gate (PLAN.md 9.15): GM12878, score >= 5, Jaccard
    **1.000000** against R Chicago (1,169/1,169 identical calls). Mouse ES
    (over R's default subsample threshold, so R's own dispersion is itself
    not run-to-run deterministic there): Jaccard **0.998665** (4,488/4,494
    shared calls, all 6 disagreements within 1% of the threshold), meeting
    the >=0.99 gate. Determinism: byte-identical repeat runs
    (single-threaded). Defaults regression: all 76 existing chic-tool
    `equiv.py` cases unaffected.
  - **Merged 2026-09-22, `dab5fe2d`:** `chicChicagoBackgroundModel` and
    `chicChicagoScores` also accept `--matrices` (cool/h5/`.hic`), a
    required, mutually-exclusive alternative to `--chinput`, requested
    directly by the project owner to match how `chicViewpoint` already
    takes real matrices rather than CHiCAGO's own bespoke pre-processed
    format. `hicx::chicago::chinput_from_matrices` derives each
    `(baitID, otherEndID)` pair's `N` from the matrix via the existing
    genuinely partial, region-scoped `ToolMatrix::load`, summing every bin
    pair overlapping each `.rmap`/`.baitmap` fragment's span (one cell at
    exact fragment resolution, many at a fixed bin size); `distSign`'s
    round-half-up convention was reverse-engineered and confirmed against
    all 270,441 real cis rows of the GM12878 chinput fixture. Scoped to
    cis-within-`maxLBrownEst`; trans and farther-cis still need `--chinput`.
    Verified independently by the orchestrating session: a fragment-
    resolution cool matrix built from GM12878's own real `N` values
    reproduces byte-identical background-model and scores output to the
    original chinput file, reproduced from scratch on a clean export
    (twice: the isolated branch and the merged tree). `hicx_tests` 314/314
    (796,405 assertions), ctest 9/9.
  - Known, explicitly scoped-out gaps: multi-replicate merging (R's
    `mergeSamples`) is not implemented, so a caller with several replicates
    must sum `N` per bait/other-end pair before calling these tools;
    shrunken normalisation is not implemented (irrelevant to Chicago's own
    defaults); no formal memory/CPU-vs-R benchmark (qualitatively much
    faster: 0.7 s CPU for the full GM12878 fit against R's ~14 s CPU for
    the same input).
  - `hicx_tests` 311/311 (177,671 assertions), ctest 9/9 (2 pre-existing
    environment-gated skips). Independently reproduced by the orchestrating
    session, twice (the isolated branch and the actually-merged tree): a
    clean-export build and full test suite, and the dispersion/score/
    significant-call numbers above reproduced from scratch on real GM12878
    data (265,494 score rows, 1,169 significant calls), not just re-run
    from the agent's own report.

**Harness: Python reference cache and parallel scheduling, merged in
`04f948ce` (contract rule 13).**
Verified on clean exports of `86efecba` and `c34f0d7d` (reproduced).

| Run | Cases | Wall | Python |
|---|---|---|---|
| Cold, into an empty cache | 483 of 483 | 2,472 s | 499 reference processes |
| Warm | 483 of 483 | 385 s | none |

- The warm run's verdicts are identical to the cold run's. A shim around the
  reference interpreter logged only the 8 validator calls that depend on C++
  output.
- `cache verify --sample 20`: 0 failures.
- Invalidation tests pass, but one assertion is flaky (open work).
- No case needed the rerun-alone rule in either run.
- Before the no-history rule, a cold run failed two `hicPlotAverageRegions`
  time gates under parallel load.

**GUI foundation (PLAN 10.1 to 10.3), merged in `f4c15dc9`:**
- `hicx::cli` parses every tool's command line as Python argparse does, and
  emits `--help-json`.
- The `hicx_matrix` module (behind `HICX_BUILD_PYTHON`) reads only the
  requested region.
- `hicexplorer_gui.workflow` runs YAML workflows headless.

**GUI shell and matrix browser (PLAN 10.4 and 10.5), merged in `0b3b56db`.**
Verified on a clean export of `b87172c2` (reproduced): gui suite 108 passed, 1
skipped (snakemake not installed), bindings 99 of 99; C++ outside `cpp/python`
unchanged.
- Generated forms for all 30 tools; each filled form's command line is parsed
  by the Python tool's `parse_arguments()` into the namespace the form values
  mean.
- Runs go through the workflow engine; form runs may write anywhere, and the
  history records absolute paths.
- The matrix browser matches the oracles exactly: cooler, hicmatrix, and
  hicstraw including the 40 GB version 7 file. It peaks at 257 MB over a
  scripted session on gm12878_chr1.cool and that file.
- Screenshots at 1280x720, 1920x1080 and 3840x2160 were reviewed.
  - Side-by-side tracks and axis ranges past the data were fixed before the
    merge.
  - Left open, cosmetic: the track legend overlaps the signal, and the square
    view leaves empty space in wide windows.
- Browser limits: one chromosome at a time, and h5 files with variable bin
  sizes are refused inline.

**GUI analysis views and workflow templates (PLAN 10.6 and 10.7), merged in
`e9c5bbaa`.** Verified on a clean export of `cb8c67e5` (reproduced): build and
ctest 6 of 6, gui suite 134 passed and 1 skipped (snakemake), bindings 99, the
browser tests skipping with a reason without `HICX_LARGE_HIC`, and the cache
test 15 of 15 in five consecutive runs. No C++ changed.
- **Views**, linked to the matrix browser: QC report, distance decay,
  hicPlotViewpoint and chicPlotViewpoint, aggregate contacts, compartment
  saddle, correlation, and hicDifferentialTAD results with a volcano plot.
- **Figure export** equals the CLI figure (byte-identical or E0).
- **Unavailable, with their PLAN sections:** HiCRep, differential loops and
  compartments.
- **Templates:** Hi-C, differential TADs, capture Hi-C, and conversion and QC
  run end to end on committed data. Outputs equal their logged commands, with
  two named normalisations: the QC html table id and hicInfo's `Date:` line.
- **Screenshots** reviewed at 1280x720 (differential, QC report, template
  picker).
- **Session memory:** a session opening every view peaks at 216 MB for the GUI
  process.

**Libraries,** each in its own repository with no remote and no licence yet:
- **coolercpp** (`~/src/coolercpp`): cooler's API in C++. All cool and mcool I/O
  goes through it since `5da5e074`; `core/src/cool_adapter.cpp` keeps what
  hicmatrix adds.
  - Pinned at `fe4d6c4` (0.3).
  - Harness 163 of 163 at `b1791be` (reproduced from a clean export). Its
    consumer test failed there because it asked for version 0.2; `fe4d6c4` fixes
    that.
- **hicfilecpp** (`~/src/hicfilecpp`): Juicer `.hic`, merged in `8a2fa526` and
  `bfcf75fa` (PLAN 9.1).
  - Reads versions 6 to 9 and writes 8 and 9. Writing 6 and 7 is refused,
    because no obtainable Juicer tools release writes them.
  - Pinned at `af08f98` (0.3).
  - Its harness passes 116 of 116 at `99615bd` against hicstraw 1.3.1 and
    Juicer tools 1.22.01 and 2.20.00 (reproduced from a clean export).
  - The 40 GB version 7 GM12878 file converted to cool at 25 kb is E1 against
    the Python (898,978,865 pixels), at 175 s and 5.42 GB against 1,951 s and
    5.25 GB (reproduced).
  - Converting the 423 MB GSM6505198 `.hic` to cool is E1 against the Python,
    at 38 s and 367 MB against 272 s and 462 MB (reproduced).

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
| hicConvertFormat | 1 | E1 cool / E2 h5 / E0 text | 1.3W | partial | written | done | 33/33 (reproduced) | `.hic` versions 6 to 9 read into cool, mcool, h5 and text, and versions 8 and 9 written from h5, cool and mcool, through hicfilecpp. Versions 6 and 7 are refused on output. `.hic` to cool is E1 against hic2cool, the new directions EX (F59). `--chromosome` refuses explicitly. Reproduces F15. Open: h5 `extra_list` dtype for hicpro and homer sources. |
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
| hicCorrectMatrix | 3 | ED ICE / EN KR | 1.2W | weak | written | done | 26/26 (reproduced) | `diagnostic_plot` drawn through `hicexplorer_plot`, E0. ICE and KR on `gm12878_chr1.cool` from about 9.3 GB to 830 MB (reproduced). |
| hicPCA | 3 | ED or better | 1.1W + 1.15D | weak | written | done | 13/13 (reproduced) | Covariance equal to `np.cov` bit for bit, because eigenvectors are chosen by position from a spectrum with a 169-fold largest eigenvalue (F8, F31). Reproduces F25. |
| hicCompartmentalization | 3 | E0 `_dat` / E2 npz | 1.3W | weak | written | done | 15/15 (reproduced) | The figure is drawn through `hicexplorer_plot`, E0; `--noPlot` is kept as an optional C++-only option. |
| hicInterIntraTAD | 3 | E0 | 1.3W | weak | written | done | 5/5 (reproduced) | The ratio plot is drawn through `hicexplorer_plot`, E0. Reproduces F33 to F35. |
| hicPlotSVL | 3 | E0 text | 1.3W | partial | written | done | 14/14 (reproduced) | `gm12878_chr1.cool`: 3,846 to 829 MB against a 1,028 MB budget, after the single-chromosome cooler load was made to cut in place. The boxplot is drawn through `hicexplorer_plot`, E0. Reproduces F37 and F51; F38 is a deviation. |
| hicBuildMatrix | 4 | E0 BAM and QC / E1 / E2 | `2*nnz_out*12 + threads*64 MB + C` | partial | written | done | 17/17 (reproduced) | QC PNGs and `hicQC.html` drawn through `hicexplorer_plot`, E0 (the HTML after a named normalisation). Reproduces F29, F30. |
| hicBuildMatrixMicroC | 4 | as hicBuildMatrix | as hicBuildMatrix | weak | **yes, not written** | done | 4/4 (reproduced) | Landed without a characterization test, contrary to rule 1. |
| hicQuickQC | 4 | E0 | 192 MB | good | written | done | 5/5 (reproduced) | Reproduces F40. |
| hicFindTADs | 5 | E0 text / E1 / E2 | per case, `2.2 * Wz + C` | partial | written | done | 8/8 (reproduced) | F22. |
| hicDetectLoops | 5 | E0 text / EN fitted size / E5 calls | 2.2W | partial | written | done | 13/13 (reproduced) | Uses Cephes `betainc` in both its float64 and float32 paths, where scipy 1.14 uses Boost.Math for both (F42, open work). `gm12878_chr1.cool` peak 256 MB against Python's 619 MB since the band is cut while reading through coolercpp (`5da5e074`, identical output), down from 877 MB. |
| hicDifferentialTAD | 5 | E0, E5 calls | 1.3W x 2 | good | written | done | 22/22 (reproduced) | Jaccard 1.0 on every case. C++-only `--sharedMask` and `--correctForMultipleTesting` (PLAN 9.7 step 1), E0 against a Python reference. On GSE234292 replicates the null drops from 20.5 % of TADs to 0.9 %, and to 0 with FDR (reproduced). |
| hicMergeDomains | 5 | E0 text / E7 PDF | none (reads no matrix) | weak | written | done | 11/11 (reproduced) | Its baseline failures were environmental (F41). DOT source byte identical, rendered by the external graphviz `dot`, found through a per-case `path_prepend`. The PDFs are E7 because `dot` embeds a creation date in compressed streams. Reproduces F39 and F49. |
| hicAggregateContacts | 5 | E0 tables / E2 h5 / E1 cool | 1.3W, z-score 559 MB | weak | written | done | 33/33 (reproduced) | KMeans and ward labels identical in all 108 fits the Python made on real data. z-score case 4.0 GB to 335 MB. Figures, including the diagnostic heatmap and the 3d plot, drawn through `hicexplorer_plot`, E0; `--noPlot` optional. Reproduces F52 to F58. |
| chicQualityControl | 6 | E0 text | 1.3W | weak | written | done | 5/5 (reproduced) | Figures drawn through `hicexplorer_plot`, E0. Reproduces F36, F43. |
| chicViewpointBackgroundModel | 6 | exact columns; fitted size and prob by likelihood and downstream E5 | 1.3W | weak | written | done | 5/5 (reproduced) | EN withdrawn for size and prob (F50). Each fit's likelihood must be no worse than the worst of five reference runs by more than `max(spread_i, T_well)`; downstream calls Jaccard 1.000000 on all four fitted cases. Reproduces F45. |
| chicViewpoint | 6 | E1 / bit identical values | 1.3W | weak | written | done | 6/6 (reproduced) | Bit identical through Boost.Math `ibeta`, as scipy uses (F42). Reproduces F44, F46, F47; F48 is a deviation. |
| chicSignificantInteractions | 6 | E0 HDF5 (implies E5 calls) | 1.3W | weak | written | done | 12/12 (reproduced) | C++-only `--correctForMultipleTesting` (PLAN 9.7), E0 against a Python reference that adjusts the tool's p-values. 15 MB and under 0.01 s against the Python's 218 MB and 1.55 s. Reproduces F62. |
| chicAggregateStatistic | 6 | E1 | 1.3W | weak | written | done | 9/9 (reproduced) | Reproduces F63. |
| chicDifferentialTest | 6 | declared E3, measured bit-identical | 1.3W | weak | written | done | 9/9 (reproduced) | Reproduces scipy 1.14.1 exactly: Boost.Math hypergeometric for `fisher_exact`, and Cephes `chdtrc` and `igami` for `chi2_contingency` and `chi2.ppf`, bit-identical on 45,081 checked values. Also has `--correctForMultipleTesting`. Reproduces F64. |
| chicExportData | 6 | E0 text / E1 bigWig | 1.3W | partial | written | done | 22/22 (reproduced) | bigWig through libBigWig 0.4.8, compared through pyBigWig interval by interval, plus zoom summaries. The archive lists members sorted, with time 0. |
| hicPlotMatrix | 7 | E0 figures | declared 90 and 130 MB | weak | yes | done | 23/23 (reproduced) | C++ compute plus matplotlib drawing (tier 7 option a). The Li et al. 2015 whole-matrix case, behind a 120 GB `skipif` in the Python tests, runs at 14.0 GB against the Python's 14.6 GB, almost all of it matplotlib. |
| hicPlotTADs | 7 | E0 | n/a | none | no | done | 4/4 (reproduced) | Delegation to pyGenomeTracks 3.9, checked at run time. |
| hicPlotViewpoint | 7 | E0 data and figures | 1.3W | weak | yes | done | 9/9 (reproduced) | Exits 1 on a bad region where the Python exits 0. |
| hicPlotAverageRegions | 7 | E0 | n/a | weak | yes | done | 8/8 (reproduced) | |
| hicPlotDistVsCounts | 7 | E0 data and figures | 1.3W | weak | yes | done | 7/7 (reproduced) | gm12878_chr1 per chromosome: 934 MB against the Python's 4,410 MB. |
| hicCorrelate | 7 | E0, one figure E6 (RMS 1.17) | 1.3W x n, three-matrix case alpha 4.5 | weak | yes | done | 11/11 (reproduced) | Complete-linkage clustering with scipy's leaf order. Refuses non-finite values; mixed float32/float64 inputs meet only E3. Exits 1 for labels and range errors where the Python exits 0. |
| hicPrepareQCreport (alias hicQC) | 7 | E0 tables and charts, HTML E0 after `pandas_styler_uuid` | n/a | none | yes | done | 6/6 (reproduced) | Exits 1 for label-row mismatches where the Python exits 0. |
| chicPlotViewpoint | 7 | E0 tar members | 1.3W | weak | yes | done | 8/8 (reproduced) | |
| hicTADClassifier | 8 | E7 | 1.3W | weak | yes | not started | - | |
| hicTrainTADClassifier | 8 | E7 | 1.3W | weak | yes | not started | - | F13. |
| hicHyperoptDetectLoops | 8 | E7 | inherits | partial | yes | not started | - | |
| hicHyperoptDetectLoopsHiCCUPS | 8 | E7 | n/a | weak | yes | not started | - | |

**Tier 7 is done (option a, merged in `b71510ba`).**
- C++ computes each figure's data, and `hicexplorer_plot` draws it with the
  Python's matplotlib calls.
- The drawing interpreter is `HICX_PLOT_PYTHON`, and must have matplotlib 3.8.4
  (pyGenomeTracks 3.9 for hicPlotTADs).
  - A wrong or missing package exits 3 before any input is read or output
    created.
  - `HICX_PLOT_ALLOW_UNPINNED=1` draws anyway, with a warning.
- The figures formerly refused under contract rule 7 are drawn.

**Tier 8 awaits the project owner's decision.**

## Counts

| tier | tools | ported | not started |
|---|---|---|---|
| 1 | 6 | 6 | 0 |
| 2 | 7 | 7 | 0 |
| 3 | 6 | 6 | 0 |
| 4 | 3 | 3 | 0 |
| 5 | 5 | 5 | 0 |
| 6 | 7 | 7 | 0 |
| 7 | 8 | 8 | 0 |
| 8 | 4 | 0 | 4 |
| **total** | **46** | **42** | **4** |

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
| hicPCA | on mm9_reduced_chr1 the C++ takes 312 s of wall time against the Python's 170 s. The Python runs multithreaded OpenBLAS (4,092 s of CPU); the C++ dense eigensolver runs on one thread (314 s of CPU) | the CPU-time gate passes, but wall time is worse. Use a multithreaded LAPACK for the dense path, keeping the covariance bits and the eigenvector choice identical; tier 11 addresses memory |
| hicQuickQC | its charts are drawn but not compared: their legends carry the random temporary file name, so widths vary by up to 21 px between runs of either tool | E7 for those charts |
| hicBuildMatrixMicroC | QC figures drawn, but no case compares them | add cases |
| drawing tools | the C++ saves the matrix and then draws, where the Python renders the QC report before saving | a failed drawing leaves the matrix behind |
| GUI differential view | the volcano plot's y axis reads "-log10 p" while its heading says adjusted p-value | cosmetic |
| GUI chicPlotViewpoint view | its data holds distances to the reference point, not genomic positions | the view cannot move the matrix browser, and says so |
| GUI differential template | the committed data has one replicate per condition | the template sums a single matrix per condition; a replicated real-data run is still open |
| time gate | C++ CPU time grows with parallel load, and drawing tools sit near a ratio of 1.0 because matplotlib dominates both sides | contract rule 13 schedules and reruns such cases alone; the rerun rule is proven only by a fake case so far |
| Python characterization tests of the plotting tools | they call `main()` directly | they cannot run against the C++ entry points |
| hicAggregateContacts | numpy's argsort fallback for CPUs without AVX-512 was not checked against numpy, because the development machine always takes the AVX-512 path | contact-pair line order could differ on such a CPU |
| hicAggregateContacts | z-score in modes `all` and `inter-chr` has no harness case | the Python needs a dense matrix of several GB |
| hicBuildMatrixMicroC | no characterization test | rule 1 not met |
| h5 writer | fixed at the source: attributes on datasets were opened with `H5Oopen` and closed as groups, leaving 27 identifiers open on a small matrix. `H5close()` before the drawing exec remains only as a safety net | a unit test counts open identifiers after a write |
| packaging | Fetched at configure time: Boost.Math, x86-simd-sort, libBigWig, pybind11, and scipy 1.14.1's Cephes headers (per-file SHA256). hicMergeDomains needs graphviz `dot` at run time, and the plotting tools need a Python with matplotlib 3.8.4 | an offline package build needs vendored tarballs, and the package needs graphviz and the drawing environment as dependencies |
| harness | `cpp_args` passes options to the C++ side only; `path_prepend` changes one case's `PATH` | both allowed only where they do not change compared outputs |
| hicPCA | dense per-chromosome matrices, needed for the bit-exact eigenvector choice (PLAN 5.4) | memory grows with the square of the largest chromosome's bin count: 269 MB of the 344 MB peak on `small_test_matrix` (5,801 bins), about 5 GB per matrix for human chr1 at 10 kb. A sparse Lanczos option with implicit centering is proposed, pending the project owner |
| hicDifferentialTAD, chic tools, hicDetectLoops | false positives: per-sample bin filtering and no multiple-testing correction (PLAN 9.7) | 20.5 % of TADs called between replicates. Step 1 is merged in `bfcf75fa` as hicDifferentialTAD options, not defaults; the chic tools and hicDetectLoops still lack a correction |
| `core/include/hicx/sparse_matrix.hpp` | values held as double regardless of dtype | int32 and float32 matrices cost twice what they need |
| budget formula | `--runningWindow` and hicFindTADs do not fit `alpha * W_input` | both need a term in the working set they build |
| usage texts | the stored usage of hicBuildMatrix and hicQuickQC wraps differently from argparse; the usage of hicPCA, hicCompartmentalization and hicDifferentialTAD lists their C++-only options | no case compares usage text |

**Fixed on 2026-09-14, `.hic` to cool memory.** The conversion streams block
column by block column (`0836d72a`). On the 40 GB version 7 file at 25 kb, peak
RSS fell from 5,422 MB to 408 MB and CPU from 176 s to 147 s, with identical
output (reproduced). The Python needs 5.25 GB and 1,951 s.

**Fixed on 2026-09-14, argument parsing.** The shared argument layer
(`f4c15dc9`) replaced the per-tool parsers.
- Behaviour changed only where the old C++ parsers differed from argparse:
  - prefix abbreviations of long options;
  - Python's `int()` and `float()` rules;
  - hicCorrectMatrix's subcommand-only options and its no-subcommand exit;
  - `FileType('r')` inputs checked while parsing.
- No case depended on any of them.
- The only spec difference allowed is `hic` among hicConvertFormat's output
  formats.

Fixed on 2026-09-14: whole-table and chromosome cool reads go through
`bin1_offset` again (`bfcf75fa`). C++ CPU time on gm12878_chr1 is now below the
pre-coolercpp figures, and peak RSS below them too (reproduced):

| case | before coolercpp | after coolercpp | now |
|---|---|---|---|
| ICE | 26.58 s | 32.78 s | 26.06 s, 746 MB |
| KR | 17.21 s | 19.01 s | 15.64 s, 745 MB |
| hicPlotSVL | 1.73 s | 2.52 s | 1.57 s, 735 MB |
| roundtrip | 6.78 s | 8.25 s | 6.93 s, 741 MB |

The roundtrip is the one exception: 6.93 s, 2 % above its pre-coolercpp figure.

Fixed in an earlier round and recorded for history: `test_transform_ops.cpp` wrote its
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
| **F59** | hicConvertFormat (hic2cool) | the installed hic2cool reports `__version__` 0.8.3 while its pip metadata says 1.0.1, so cool files converted from `.hic` say `generated-by: hic2cool-0.8.3`; the port reproduces the string. |
| **F60** | hicAdjustMatrix | `--action mask --regions` zeroed 4,970 rows at 50 kb on GSE234292 for a BED of 4,067 bins. Measured with the C++ port, which matches the Python on all 21 hicAdjustMatrix cases; not rerun on the Python, and the cause is not yet examined. It inflated an early calibration figure (PLAN 9.7). |
| **F61** | cooler (reference library) | `cooler.fileops.is_cooler` returns False for `hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1_chr2.cool`, whose `format` attribute is a fixed-length byte string, while `cooler.Cooler` opens it (3,790 bins, 3,202,457 pixels, reproduced). `hicx_matrix` follows `cooler.Cooler`. `test_data/matrix.mcool` has no `/resolutions/` groups, so it is not an mcool to cooler's layout. |
| **F62** | chicSignificantInteractions | Five defects reproduced by the port. (1) The significant file is written in `np.unique` order but takes reference points in computation order, so Sox17 carries Eya1's reference point. (2) Without preselection, the peak threshold is compared with the x-fold instead of the raw count. (3) `merge_neighbors` drops a final unmerged candidate. (4) Dual mode reads the second sample under the first sample's chromosome and gene names. (5) Threshold files written as attributes crash with a TypeError, and `mode_preselection_calue` is misspelled when there is no preselection. |
| **F63** | chicAggregateStatistic | Interval trees restart whenever the chromosome changes. A four-column target BED works only with matrices named `c_adj_norm` and `t_adj_norm`. A single-mode target file fails with a KeyError. |
| **F64** | chicDifferentialTest | A skipped reference point shifts every later written result and ends in an IndexError. The chi-squared test rejects on the statistic against the critical value, not on the p-value. |

## Deliberate deviations from the Python behaviour

| Tool / component | Deviation | Decided by |
|---|---|---|
| acceptance, all tools | three significant digits per item, not byte identity | project owner, 2026-09-01 |
| implementation, all tools | data structures and algorithms free wherever the result is the same | project owner, 2026-09-13 |
| provenance, all tools | output is to name HiCExplorer 4 (recorded, not implemented) | project owner, 2026-09-02 |
| figures, all tools | drawn by `hicexplorer_plot` with matplotlib 3.8.4 checked at run time. A wrong or missing drawing environment exits 3 before any input is read or output created; `HICX_PLOT_ALLOW_UNPINNED=1` draws anyway with a warning. `--noPlot` remains as a C++-only option | project owner (GUI request, tier 7 option a), 2026-09-15 |
| hicPlotViewpoint, hicCorrelate, hicPrepareQCreport | exit 1 for a bad region, labels and range errors, and label-row mismatches, where the Python exits 0 | implementing agent, reported |
| chic tools | where the Python hangs on an HDF5 group name that is already taken, the C++ exits 1. chicExportData's archive lists members sorted, since the Python's `os.walk` order is not reproducible, and stamps them with time 0 | implementing agent, reported |
| capture Hi-C workflow template | its multiple-testing correction defaults to `fdr`, so the template differs from the Python tools unless it is set to `none` | orchestrating session, 2026-09-15 (false-positive control, PLAN 9.7) |
| KR | no `exit(0)`; deterministic; class EN; reimplemented, not vendored | orchestrating session |
| hangs in the reference | hicPlotSVL worker exception (F38), chicViewpoint duplicate gene (F48), hicMergeDomains start coordinate 0 (F49): the port exits 1 with a message | orchestrating session |
| hicMergeDomains | an unknown `-of` format or a missing `dot` is refused before anything is written, where the Python writes text files first and then crashes | implementing agent, reported |
| hicAggregateContacts | an empty cluster exits 1 before writing, where the Python writes earlier clusters and then crashes | implementing agent, reported |
| chicViewpointBackgroundModel | fitted parameters judged by per-distribution likelihood and downstream calls, not EN | orchestrating session |
| hicDetectLoops, hicDifferentialTAD, hicInterIntraTAD | output independent of the thread count (F23, F33) | implementing agents, reported |
| hicFindTADs | no abort at 12 or more processes (F22) | implementing agent, reported |
| hicConvertFormat | `--chromosome` refuses explicitly | implementing agent, reported |
| hicDifferentialTAD | C++-only `--sharedMask` and `--correctForMultipleTesting` (PLAN 5.8 dual mode); the defaults stay byte-identical to the Python | project owner request (false positives), orchestrating session |
| hicConvertFormat, hicfilecpp | writing `.hic` versions 6 and 7 is refused with the reason; there is no Juicer tools release to validate a writer against | orchestrating session, 2026-09-14 |
| hicConvertFormat, hicfilecpp | `.hic` output (new, EX): bytes differ from Juicer tools' while pixels, expected values and VC, VC_SQRT, KR and SCALE vectors match. No FRAG, GW or INTER norms, and none of `pre`'s filters or statistics. Norms are always computed, where Juicer skips them when its heap looks small. The matrix is written as loaded, correction applied unless `--load_raw_values`. A missing norm, chromosome or resolution is an error on reading | implementing agent, reported; recorded in hicfilecpp `docs/DEVIATIONS.md` |
| hicFindEnrichedContacts | not ported (F12) | plan |

## Open questions

1. Tier 7 is taken as settled in favour of option (a), C++ computing and
   matplotlib drawing, by the GUI request of 2026-09-13 (PLAN tier 10). The
   project owner has not yet confirmed that reading. Tier 8 is still open.
2. Licences of coolercpp and hicfilecpp: GPL-3 like HiCExplorer, or BSD-3 like
   cooler, from whose source some of coolercpp's components derive (project
   owner, before publication).
3. macOS verification for the GUI and the C++ core on arm64 needs a Mac or a CI
   runner. CI means pushing, which needs the project owner's approval.
