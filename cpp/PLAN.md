# HiCExplorer v4 (C++) - architecture and porting plan

Owner: supervising agent. Companion ledger: `cpp/STATUS.md`. Environment facts:
`cpp/AGENTS_CONTRACT.md`. Date of this revision: 2026-09-01.

## 1. What is being ported

The Python reference in this worktree is HiCExplorer 3.7.7-dev: 46 tool modules
under `hicexplorer/` (18,089 lines) plus five shared helpers
(`utilities.py` 729, `parserCommon.py` 43, `readBed.py` 252, `reduceMatrix.py`
232, `iterativeCorrection.py` 86) and four library modules under
`hicexplorer/lib/` (`buildMatrixMethods.py`, `viewpoint.py`, `cnb.py`,
`tadClassifier.py`).

`setup.py:131-142` installs 47 console scripts. The extra script is `bin/hicQC`,
an alias that calls `hicexplorer.hicPrepareQCreport.main`. `bin/hicFindEnrichedContacts`
is dead: it imports `hicexplorer.hicFindEnrichedContacts`, a module that does not
exist. `bin/hicexplorer` is `list_tools.py`, a help banner. So the real tool count
is **46 executables + 1 alias**, and the port must produce 46 binaries plus a
`hicQC` alias and a `hicexplorer` banner.

The matrix layer is **not** in this repository. It is the external package
`hicmatrix` (HiCMatrix 17.2, 2,034 lines) in
`~/miniconda3/envs/__hicexplorer@3.7.6/lib/python3.12/site-packages/hicmatrix/`.
It must be reimplemented in full as part of the core library; it is the single
largest source of behavioural coupling, and it contains defects that HiCExplorer
compensates for (section 2.7).

Reference environment versions that define the oracle: numpy 1.26.4,
scipy 1.14.1, cooler 0.10.2 (format-version 3), PyTables 3.10.1, h5py 3.12.1,
pandas 2.2.3, intervaltree 3.1.0, hic2cool 1.0.1, krbalancing 0.0.5,
pyBigWig 0.3.22, pybedtools 0.10.0, scikit-learn 1.3.x, imbalanced-learn,
cleanlab 2.6.6, hyperopt 0.2.7, pyGenomeTracks 3.9.

## 2. Core library design (`cpp/core`, target `libhicx4`)

### 2.1 Bin and region model

The Python canonical bin table is `cut_intervals`: a `list[(chrom:str,
start:int, end:int, extra:float)]` of length `nbins`, in genome order. `extra`
is a per-bin coverage/interactions value; the cool reader hardcodes it to `1.0`
(`hicmatrix/lib/cool.py:234-235`), the h5 reader restores whatever was written,
and `hicMergeTADbins`/`hicMergeMatrixBins` recompute it as a mean.

C++ model:

```cpp
namespace hicx {
struct CutIntervals {                 // struct-of-arrays, not array-of-structs
  std::vector<uint32_t> chrom_id;     // index into chrom_names
  std::vector<int64_t>  start;        // 0-based, half-open
  std::vector<int64_t>  end;
  std::vector<double>   extra;
  std::vector<std::string> chrom_names;   // genome order of first appearance
  std::vector<int64_t>     chrom_length;  // from cool /chroms/length, or 0
};
}
```

Struct-of-arrays because every consumer works column-wise (`zip(*cut_intervals)`
appears in `hicMergeTADbins.py:56`, `hicPCA.py:307`, `hicCorrectMatrix.py`, and
throughout `hicmatrix`), and because chrom names as an id vector makes the
str/bytes coercion problem (`utilities.py:641 check_chrom_str_bytes`) disappear:
one interning table, one comparison.

`chrom_names` must preserve **file order**, not sorted order. `hicSumMatrices.py:52`
compares `hic.chrBinBoundaries != hic_to_append.chrBinBoundaries` and aborts when
the order differs, so ordering is observable behaviour.

### 2.2 Chromosome boundaries and interval lookup

Python builds, per chromosome, an `intervaltree.IntervalTree` of
`Interval(start, end, bin_id)` plus an `OrderedDict chrBinBoundaries[chrom] =
(first_bin, last_bin_exclusive)` (`hicmatrix/HiCMatrix.py:986-1020`). All queries
are half-open **point** lookups (`tree[pos:pos+1]`, sorted, `[0].data`;
`HiCMatrix.py:261-262`). No general interval-overlap query is ever used.

Therefore an interval tree is the wrong data structure for the port. Bins within a
chromosome are contiguous and sorted, so:

```cpp
struct BinIndex {
  // chrom_id -> [first_bin, last_bin_exclusive)
  std::vector<std::pair<int64_t,int64_t>> chrom_range;
  int64_t bin_at(uint32_t chrom, int64_t pos) const;  // binary search on start[]
};
```

`bin_at` is `std::upper_bound(start.begin()+lo, start.begin()+hi, pos) - 1`,
then a check that `pos < end[bin]`. O(log n) instead of O(log n) with allocation,
and it removes the intervaltree dependency entirely. Two behaviours must be kept:

- A chromosome that appears twice non-contiguously silently overwrites its tree
  and boundary entry in Python (`HiCMatrix.py:986-1020`). The C++ builder must
  reproduce "last occurrence wins" or, better, detect it and fail loudly; the
  choice must be recorded in `STATUS.md` because it changes behaviour on
  malformed input.
- `getRegionBinRange` returns `None` on an out-of-range query (`HiCMatrix.py:223-270`),
  and callers branch on that. Model it as `std::optional<std::pair<int64_t,int64_t>>`,
  with the end bin **inclusive** as Python has it.

### 2.3 Sparse contact matrix

Python state is a `scipy.sparse.csr_matrix`, held **full symmetric in memory**
and **upper-triangular on disk**. `fillLowerTriangle`
(`HiCMatrix.py:106-120`) symmetrizes on load with `m + triu(m,1).T`, and every
writer re-imposes `triu(k=0)` (`cool.py:289`, `h5.py:118`, `ginteractions.py:21`).

```cpp
template <class T>
struct CsrMatrix {          // canonical: full symmetric, sorted indices, no explicit zeros
  int64_t n = 0;            // square
  std::vector<int64_t> indptr;   // n+1
  std::vector<int32_t> indices;  // widen to int64 when n > INT32_MAX
  std::vector<T>       data;
};
using CsrI32 = CsrMatrix<int32_t>;
using CsrF64 = CsrMatrix<double>;
```

Rules the port must honour, because they are observable:

1. **Duplicate pixel accumulation.** `csr_matrix((data,(i,j)))` sums duplicates
   (`cool.py:95`, `hicpro.py:34`). The COO-to-CSR builder must sum, not overwrite.
2. **`eliminate_zeros()` before every write** (`cool.py:265`, `h5.py:121`). An
   explicit stored zero changes `nnz`, and `nnz` is printed by `hicInfo` and
   stored in the cool `nnz` attribute. Keep an explicit `eliminate_zeros()` and
   call it at exactly the same points.
3. **h5 stores the CSR internals verbatim** (`h5.py:32-39`): `data`, `indices`,
   `indptr`, `shape`. No canonicalisation happens on load. Two logically equal
   matrices with differently ordered `indices` within a row produce different
   files. scipy's `triu` output has sorted indices; the C++ builder must too.
4. **Index dtype is chosen by the writer.** cool pixels are forced to `int32`
   for `bin1_id`/`bin2_id` by hicmatrix (`cool.py:300`) regardless of nbins; the
   cool reader picks `int32` unless `nbins > int32max` (`cool.py:68-73`). h5
   inherits scipy's choice (int32 while `nnz` and `n` fit). Reproduce exactly:
   the dtype is visible in the file.
5. **Element type.** Counts are `int32` for raw matrices and `float64` after
   correction/transformation. The type is carried in the file
   (h5 `/matrix/data` atom, cool `/pixels/count` dtype) and `hicConvertFormat
   --enforce_integer` rounds with `np.rint` (round-half-to-even) to int32
   (`cool.py:355-358`). Use `std::nearbyint` with `FE_TONEAREST`, not
   `std::round` (which is half-away-from-zero).

A `Matrix` aggregate mirrors `hiCMatrix`:

```cpp
struct Matrix {
  CsrF64 m;                                  // or a variant over int32/double
  CutIntervals bins;
  std::vector<int64_t> nan_bins;
  std::optional<std::vector<double>> correction_factors;
  std::optional<std::vector<double>> distance_counts;
  BinIndex index;                            // derived, rebuilt by every mutator
  std::optional<int64_t> bin_size;           // derived, median
  bool bin_size_homogeneous = true;
  // masking bookkeeping
  std::vector<int64_t> orig_bin_ids;         // [kept..., masked...]
  CutIntervals orig_cut_intervals;
};
```

The mutators that must rebuild `index`, mirroring
`keepOnlyTheseChr` (603), `reorderBins` (723), `maskBins` (761),
`restoreMaskedBins` (835), `reorderMatrix` (877), `update_matrix` (813),
`setCutIntervals` (121), `setMatrix` (135).

### 2.4 File format layer

One `MatrixFile` interface with a `load()/save()` pair, dispatched by a factory
equivalent to `MatrixFileHandler` (`hicmatrix/lib/matrixFileHandler.py:19`). Format
detection in Python is **purely by filename suffix** (`HiCMatrix.py:49-52`:
`.h5` means h5, everything else means cool). Reproduce that, including the
consequence that `f.mcool::/resolutions/10000` and `f.scool::/cells/x` go
through the cool reader unchanged.

| format | read | write | notes |
|---|---|---|---|
| cool | yes | yes | HDF5, cooler schema v3 |
| mcool | yes | yes (append per resolution) | **two layouts occur in the corpus.** `hicConvertFormat --outputFormat mcool` writes `/resolutions/<res>` (verified), and `test_data/hicBuildMatrix/multi_small_test_matrix.mcool` matches; but `test_data/matrix.mcool` has top-level groups `/0`../`/4` with no `/resolutions`. hicmatrix never enumerates resolutions: it opens whatever comes after `::` as an opaque group path (`HiCMatrix.py:49-52`). The C++ reader must do the same, and discover only when no `::` suffix is given (`/resolutions/*` first, then top-level numeric groups). Root attrs `format=HDF5::MCOOL`, `format-version=2` are **not** written by hicmatrix; the C++ writer must match that omission or the change must be recorded |
| scool | no (per-cell only) | yes | `cooler.create_scool` equivalent; only `hicConvertFormat` reaches it |
| h5 | yes | yes | PyTables CArray + blosc, section 2.5 |
| homer | yes | yes | dense gzipped TSV |
| ginteractions | no (Python `load` is a stub, `ginteractions.py:14-15`) | yes | BEDPE-like 7 columns |
| hicpro | yes | yes | 1-based triplet `.matrix` + `.bed` |
| 2D-text | yes | no | read in `hicConvertFormat.py:155-163` only |
| hic | via hic2cool | no | Juicer binary, section 3.6 |
| npz | yes | yes | not a matrix format in the hicmatrix sense: `scipy.sparse.save_npz`/`load_npz` written by `hicAverageRegions.py:209` and read by `hicPlotAverageRegions.py:15`. It is a ZIP container holding `.npy` arrays (`format`, `shape`, `data`, `indices`, `indptr` for CSR). Implement a minimal `.npy`/`.npz` reader and writer in `core/src/io/npz.cpp`: stored (uncompressed) ZIP entries, `.npy` v1.0 header with a Python dict literal. This is small and self-contained, and it is the only way `hicAverageRegions` output stays readable by the existing Python plotting tool |

#### cool

Groups and dtypes as written by cooler 0.10.2 through hicmatrix:

```
/chroms/name     S<w>  (fixed-width bytes, width = max name length)
/chroms/length   int32
/bins/chrom      HDF5 ENUM over int32   <- categorical, not a string column
/bins/start      int32
/bins/end        int32
/bins/weight     float32   (only when correction factors exist)
/pixels/bin1_id  int32     (hicmatrix override of cooler's int64 default)
/pixels/bin2_id  int32
/pixels/count    int32, or the matrix dtype when non-integer
/indexes/chrom_offset int64, length nchroms+1
/indexes/bin1_offset  int64, length nbins+1
```

Pixel datasets are created resizable: `shape=(min(5*nbins, nnz),)`,
`maxshape=(nnz,)`. Filters: `gzip` level 6 with `shuffle` on
(cooler defaults; hicmatrix passes no `h5opts`).

Root attributes written by cooler, then partially overwritten by hicmatrix
(`cool.py:386-426`), and only when the file is opened in mode `w`:
`format="HDF5::Cooler"`, `format-version=3`,
`format-url="https://github.com/mirnylab/cooler"` (note: hicmatrix overwrites
cooler's `open2c` URL with the old `mirnylab` one),
`bin-type`, `bin-size`, `storage-mode="symmetric-upper"`, `nchroms`, `nbins`,
`nnz`, `sum`, `genome-assembly`, `creation-date`, `generated-by="HiCMatrix-17.2"`,
`generated-by-cooler-lib="cooler-0.10.2"`,
`tool-url="https://github.com/deeptools/HiCMatrix"`, `metadata` (JSON string),
plus optional `matrix-generated-by` / `matrix-generated-by-url`.

The `generated-by` string is version-bearing and is printed by `hicInfo`. The C++
writer must emit `HiCMatrix-17.2` verbatim if `hicInfo` output is to compare
equal, or emit `hicx4-<version>` and have the harness normalise the field. **Decision:
emit `HiCMatrix-17.2`** while porting, and switch to a v4 identity only once the
whole suite is green, at which point the harness normalises `generated-by`,
`generated-by-cooler-lib`, `creation-date` and `tool-url`. Rationale: those four
fields are the only attributes that legitimately differ between implementations,
and normalising them early would hide real regressions in the other twelve.

Correction (`weight`) semantics, from `cool.py:170-230` and `:305-347`:

- On load, `a_ij <- a_ij op (c_i * c_j)` where `op` is `/` when the correction
  column name is in `{KR, VC, SQRT_VC}` and `*` otherwise. Note that cooler's
  own divisive set is `{KR, VC, VC_SQRT}`; hicmatrix spells the third one
  differently. Preserve the hicmatrix spelling.
- Skipped entirely when every weight is NaN (`cool.py:186`).
- On save, factors are inverted (`c <- 1/c`, NaN and Inf to 0) and the operator
  forced to `*` when the source was h5, or hic2cool >= 0.5 (string comparison on
  the version), or the operator was `/` (`cool.py:305-317`).
- With `pApplyCorrection=True` the counts are divided back out so that `count`
  on disk is raw (`cool.py:322-341`).
- `weight` is written as `convertNansToOnes(factors)` (`cool.py:318,345`), so NaN
  weights become 1.0 on disk, not NaN.

nan_bins are not stored in cool; they are **recovered heuristically** on load
(`cool.py:239-255`): after zeroing NaN data and eliminating zeros, a bin whose
index never occurs in `matrix.indices` and whose row is empty is a nan bin. Any
exception in that block yields `nan_bins = None`. This heuristic must be ported
literally, including its failure mode, because it determines the NaN-bin count
that `hicInfo` prints for h5 outputs derived from cool inputs.

#### h5

PyTables layout, all nodes `CArray` with `Filters(complevel=5, complib='blosc')`
(`h5.py:123`), file title `"HiCExplorer matrix"`:

```
/matrix/data /matrix/indices /matrix/indptr /matrix/shape
/intervals/chr_list /intervals/start_list /intervals/end_list /intervals/extra_list
/nan_bins             (optional, omitted when empty)
/correction_factors   (optional, omitted when None/empty)
/distance_counts      (optional)
```

`chr_list` is fixed-width bytes `S<w>`; `start_list`/`end_list` int64;
`extra_list` float64 or int64 depending on what numpy inferred;
`correction_factors` float64 with NaN mapped to 0 before write (`h5.py:157-158`).
An existing output file is `unlink`ed first (`h5.py:105-109`).

Two defects to preserve or consciously fix, recorded either way in `STATUS.md`:

- `h5.py:84` reads `f.root.correction_factors` when asked for `distance_counts`.
- `h5.py:34-37` swallows a failed `/matrix` read and then raises `KeyError`.

### 2.5 On-disk compatibility policy

There are three achievable levels and the plan commits to a different one per
format, because the formats differ in how much of the byte layout is dictated by
a third-party writer.

| level | meaning | applies to |
|---|---|---|
| **L1 byte-identical** | `cmp` of the two files succeeds | homer, ginteractions, hicpro, all bedgraph/bed/tsv/txt outputs, `hicInfo` text |
| **L2 HDF5-logically-identical** | same set of objects; same dtypes, shapes, fill values, chunk shapes and filter pipelines; every dataset decodes to bit-identical bytes; same attributes modulo the four normalised provenance fields | cool, mcool, scool |
| **L3 HDF5-value-identical** | same objects, same dtypes and shapes; every dataset decodes to bit-identical values; chunk shape and filter parameters may differ | h5 |

Why h5 gets L3 and cool gets L2:

- cool is written by `cooler.create_cooler` with `h5opts` left at cooler's
  defaults, which are three explicit, documented, reproducible settings
  (`compression="gzip"`, `compression_opts=6`, `shuffle=True`). Chunk shape is
  cooler's own fixed `CHUNK` logic. All of it is reimplementable against the
  HDF5 C API. Even so, gzip output bytes depend on the zlib build, so L2 stops
  at "decodes bit-identically" rather than "the compressed bytes match".
- h5 is written by PyTables, whose `CArray` chunk shape comes from an internal
  heuristic (`_calc_chunkshape`, expected row count, an L2-cache-size model),
  and whose blosc filter emits `FILTER_ID 32001` with parameters
  `{2, 2, typesize, blocksize, 5, 1}` where `blocksize` is chosen by blosc at
  runtime from the type size and the level. Observed on
  `test_data/small_test_matrix.h5`: chunk shape 7281 for a 33,754-element
  `S9` dataset, blosc params `2 2 9 65529 5 1`. Reimplementing PyTables' chunk
  heuristic to the element is achievable but is pure liability: it is undocumented,
  version-dependent, and buys nothing scientifically. **Decision: C++ writes h5
  with the blosc filter at complevel 5 and shuffle (so PyTables and h5py can
  read it), and the comparator decodes both files and compares values.**

Consequence for the harness: the cool comparator compares the HDF5 structure
including chunking and filters; the h5 comparator compares decoded values only.
Both are implemented as separate comparator plugins (section 9.3).

Round-trip compatibility is tested in both directions for every format in
tier 1: Python writes, C++ reads, C++ writes, Python reads, and the two Python
reads must produce identical in-memory state.

### 2.6 Math kernels

Kernels that appear in more than one tool and therefore belong in `core`:

| kernel | Python source | notes |
|---|---|---|
| ICE iterative correction | `hicexplorer/iterativeCorrection.py:11-86` | COO row/col scaling, 50 passes default, tol 1e-5 |
| KR balancing | `krbalancing` 0.0.5 (C++/Eigen) | see 3.3 |
| obs/exp (three variants) | `utilities.py:488 lieberman`, `:510 non_zero`, `:554 obs_exp_matrix` | per-diagonal means |
| expected interactions | `utilities.py:293,317,341,356` | threaded in Python |
| z-score matrix | `utilities.py:460`, `HiCMatrix.py:346-562` | per-diagonal mean and std |
| Pearson / covariance of a dense per-chromosome matrix | `hicPCA.py:296,301` (`np.corrcoef`, `np.cov`), `hicTransform.py` | the memory blowup, section 4.2 |
| dense eigendecomposition | `hicPCA.py:305` `scipy.linalg.eig` | section 5.4 |
| `reduce_matrix` (bin merging) | `hicexplorer/reduceMatrix.py` | used by `hicMergeMatrixBins`, `hicMergeTADbins`, `hicBuildMatrix` |
| negative-binomial pdf/cdf | `hicexplorer/lib/cnb.py:14,23` | `gammaln`, `betainc` |
| negative-binomial MLE fit | `fit_nbinom` 1.2 | L-BFGS-B on the NB log-likelihood, section 3.5 |
| Wilcoxon rank-sum, Anderson-Darling k-sample | `scipy.stats.ranksums`, `anderson_ksamp` | `hicDetectLoops`, `hicFindTADs`, `hicDifferentialTAD` |
| BED and narrowPeak/broadPeak reading | `hicexplorer/readBed.py`, `utilities.py:19,38` | |

### 2.7 Quirks that must be reproduced, not fixed

These are behaviour, because downstream code depends on them. Each gets a
regression unit test in `cpp/tests` and a note in `STATUS.md`.

1. **Load-time field swap.** All loaders return `(matrix, cut_intervals, nan_bins,
   distance_counts, correction_factors)` (`h5.py:89`, `cool.py:260`), but
   `HiCMatrix.py:58-59` unpacks it as `..., correction_factors, distance_counts`.
   After a file load, `hiCMatrix.correction_factors` holds the loader's
   `distance_counts` and vice versa. `hicConvertFormat.py:222-223` reads them back
   swapped to compensate, while `MatrixFile.set_matrix_variables`
   (`matrixFile.py:29`) uses the unswapped order and `hiCMatrix.save`
   (`HiCMatrix.py:95-96`) feeds it the swapped members. Net effect: an h5-to-h5
   round trip through `hiCMatrix` moves correction factors from
   `/correction_factors` into `/distance_counts`. This is observable in file
   listings and must be reproduced exactly.
2. **`truncTrans` is a no-op** (`HiCMatrix.py:902-919`): the distance list is
   unpacked as a 2-tuple where 3 values are returned, and line 917 uses `==`
   where an assignment was meant.
3. **`reorderBins` does not permute `correction_factors`** while `maskBins` and
   `reorderMatrix` do.
4. **`hicPCA` uses `scipy.linalg.eig`, not `eigh`** (`hicPCA.py:305`), on a
   symmetric covariance matrix, and then takes columns `[k-1:k]` without sorting
   by eigenvalue. See section 5.4.
5. **`hicMergeTADbins` clears `correction_factors` before saving**
   (`hicMergeTADbins.py:84`), by design.
6. **cool `weight` NaN becomes 1.0 on write** (`convertNansToOnes`).
7. **`hicInfo` reports a different nnz for the same matrix in cool and in h5.**
   Verified on this tree: `hicInfo -m Li_et_al_2015.cool` prints
   `Non-zero elements: 1,661,678` while `hicInfo -m Li_et_al_2015.h5` prints
   `3,313,107` for the identical matrix. The cool path reads the stored `nnz`
   attribute, which counts the upper triangle; the h5 path loads the matrix,
   which `fillLowerTriangle` has symmetrized, and reports `matrix.nnz`
   (`hicInfo.py:89-90` vs `:125`). The same split applies to the whole block of
   fields: the cool path prints `Number of chromosomes` and the available bins
   columns but not `Sum of matrix`, `Minimum (non zero)`, `Maximum` or
   `NaN bins`, because those keys are absent from `cooler_file.info`; the h5
   path prints those four and omits the cool-only ones. Any C++ `hicInfo` that
   "fixes" this by reporting one consistent number fails E0 on half the corpus.

## 3. Dependency decisions

All native dependencies come from `$HICX_DEPS = ~/miniconda3/envs/__hicexplorer@3.7.6`
used as a plain prefix, with RPATH to `$HICX_DEPS/lib`. Verified present: HDF5
1.14 (C, C++, HL), Eigen 3, htslib 1.21 (shared and static), zlib-ng, bzip2,
lzma, zstd, lz4, libdeflate, blosc 1.21.6 and blosc2, OpenBLAS 0.3.28 with the
reference `libblas`/`liblapack` shims.

### 3.1 HDF5: C API, not the C++ API

Use the **C API** (`hdf5.h`) behind a thin internal RAII wrapper in
`core/src/io/h5c.cpp`. Reasons, in order of weight:

1. The formats need things the C++ API exposes awkwardly or not at all:
   HDF5 ENUM types for `/bins/chrom`, variable-length strings for cool
   attributes, resizable datasets with explicit `maxshape`, and a
   **user-defined filter** (blosc, id 32001) registered via `H5Zregister`.
   All of that is plain C.
2. `libhdf5_cpp` is a thin wrapper that throws `H5::Exception`, which mixes
   badly with the error model everything else in the port will use.
3. Linking only `libhdf5` keeps the dependency surface minimal.

HighFive was considered and rejected: it would be a `FetchContent` dependency
that still requires dropping to the C API for the enum, the vlen strings and the
filter registration, so it adds a layer without removing one.

Blosc: link `$HICX_DEPS/lib/libblosc.so` and register the HDF5 filter in
`core/src/io/blosc_filter.c` (vendored from the Blosc project's `hdf5-blosc`,
about 200 lines, BSD). This is required both to read existing `.h5` test data
and to write `.h5` that PyTables can read.

### 3.2 Sparse linear algebra: hand-rolled kernels, Eigen only where it already is

Do **not** build the matrix type on `Eigen::SparseMatrix`. Reasons:

1. The hot kernels are not linear algebra. ICE is row/column scaling of a COO
   triple; obs/exp is a per-diagonal reduction; `reduce_matrix` is a grouped
   sum. All are trivially expressed on the CSR arrays directly, all are
   memory-bound, and all need to reproduce numpy's exact reduction order
   (section 5.2). Wrapping them in Eigen expressions makes the reduction order
   opaque, which is exactly the property that must stay explicit.
2. The file layer needs the raw `indptr`/`indices`/`data` arrays anyway, since
   h5 stores them verbatim. An `Eigen::SparseMatrix` would have to be
   round-tripped to those arrays at every boundary.
3. The one place a real sparse solver would help does not exist: no tool solves
   a sparse system.

Eigen **is** used in two places, both mandated by external code:
`krbalancing` (section 3.3) and any dense LAPACK-adjacent work in `hicPCA`
(section 3.4), where Eigen is the least painful way to hold a dense
column-major block. `Eigen::SparseMatrix<double,0,long>` types appear only at
the krbalancing boundary.

### 3.3 KR balancing: vendor the krbalancing C++ source

`krbalancing` 0.0.5 is **already C++**. `nm -DC` on the installed
`krbalancing.cpython-312-x86_64-linux-gnu.so` shows
`kr_balancing::kr_balancing(long const&, long const&, long const&,
Eigen::Ref<Matrix<long,-1,1>>, Eigen::Ref<Matrix<long,-1,1>>,
Eigen::Ref<Matrix<double,-1,1>>)`, `computeKR`, `inner_loop`, `outer_loop`,
`rescale_norm_vector`, `get_normalised_matrix(bool&)`,
`get_normalisation_vector(bool&)`, over `Eigen::SparseMatrix<double,0,long>`.
It is a pybind11 shell over an Eigen implementation.

**Decision: vendor the upstream C++ source (deeptools/Knight-Ruiz-Matrix-balancing-algorithm)
into `cpp/core/src/math/krbalancing/` and call it directly, dropping pybind11.**
This is the single highest-leverage decision in the plan: it makes KR
bit-identical by construction rather than by tolerance, and removes the only
compute kernel that would otherwise have to be reverse-engineered from a binary.
The call site to reproduce is `hicCorrectMatrix.py:718-731` (per chromosome) and
`:744-752` (whole matrix); note the interface takes CSR `indptr`/`indices` as
`int64` and `data` as `float64`, and that `get_normalisation_vector(True)` is used
for the whole-matrix path while `get_normalisation_vector(False)` is used per
chromosome (`hicCorrectMatrix.py:730` vs `:753`).

Fallback if the source cannot be fetched through the proxy: reimplement
Knight-Ruiz from the paper and accept L3 float tolerance for KR only,
recorded in `STATUS.md`.

### 3.4 Dense eigendecomposition and BLAS

Link OpenBLAS from `$HICX_DEPS` (`libopenblasp-r0.3.28.so`, which is what numpy
and scipy in that env are built against) and call LAPACK `dgeev` directly for
`hicPCA`, and `dsyevr`/`dpotrf` nowhere (nothing needs them). Using the *same*
OpenBLAS build as the oracle is what makes the eigenvector column order and
signs match (section 5.4). Do not use Eigen's own `EigenSolver`: it is a
different algorithm and will not reproduce LAPACK's output ordering.

### 3.5 The remaining Python-only dependencies

| Python dependency | used by | C++ decision |
|---|---|---|
| `pysam` | `lib/buildMatrixMethods.py` only | **htslib 1.21 directly.** `pysam` is a Cython wrapper over the same library. Read BAM with `sam_open`/`sam_hdr_read`/`bam_read1`, use `bam_aux_get` for `SA`/`NM`, and the `bam1_core_t` flags. Mate-pair iteration in `readBamFiles` (`buildMatrixMethods.py:458`) assumes name-sorted paired BAMs read in lockstep, which maps directly. |
| `krbalancing` | `hicCorrectMatrix` | vendor the C++ source, section 3.3 |
| `pyBigWig` | `hicPCA` (read `--extraTrack`, write `.bw`), `hicPlotMatrix` (read), `chicExportData` (write) | **vendor libBigWig** (the C library pyBigWig wraps; `FetchContent` from GitHub, MIT). Writing bigWig by hand is not worth it: the zoom-level and R-tree index construction is fiddly and libBigWig is the exact code the oracle runs. |
| `pybedtools` | `hicMergeLoops` (`BedTool.merge`), `hicValidateLocations`, `chicSignificantInteractions:515`, `lib/tadClassifier.py` | **reimplement.** Only `BedTool(...)`, `.sort()`, `.merge()` and `.intersect()` on small in-memory interval sets are used. That is 200 lines of sort-and-sweep and avoids a `bedtools` binary dependency at runtime. Sort order must match bedtools' lexicographic chrom sort, which is *not* the matrix chrom order; this is the likeliest source of ordering diffs and gets a dedicated unit test. |
| `fit_nbinom` | `hicDetectLoops:161`, `chicViewpointBackgroundModel:234` | **reimplement.** It is 60 lines: the NB log-likelihood with `gammaln`, its analytic gradient with `psi`, and `scipy.optimize.fmin_l_bfgs_b` from an initial `(r, p)`. Port the likelihood and gradient verbatim, and vendor a L-BFGS-B implementation (the original Nocedal Fortran translation, or `LBFGSpp` via FetchContent). Different L-BFGS-B stopping behaviour is the main risk; this is why `hicDetectLoops` and `chicViewpointBackgroundModel` get loose tolerances (section 5.5). |
| `scipy.special` (`gammaln`, `psi`, `betainc`) | `lib/cnb.py`, `fit_nbinom` | **vendor Cephes.** scipy's `gammaln`, `psi` and `betainc` are Cephes routines; using the same Cephes source gives bit-identical results, whereas `std::lgamma` does not agree with Cephes in the last ulp. This matters because `cnb.cdf` feeds p-value thresholds that get compared against fixed cutoffs. |
| `scipy.stats.ranksums` | `hicDetectLoops`, `hicFindTADs`, `hicDifferentialTAD` | reimplement: rank transform with tie averaging, normal approximation, `erfc`-based two-sided p-value. Deterministic and easy to match exactly. |
| `scipy.stats.anderson_ksamp` | `hicDetectLoops` | reimplement the Scholz-Stephens k-sample statistic and scipy's interpolation table for the p-value. Table values must be copied from scipy's source. |
| `scipy.stats.fisher_exact`, `scipy.stats.chi2_contingency`, `scipy.stats.chi2.ppf` | `chicDifferentialTest.py:91-130` | reimplement. Fisher on a 2x2 table is the hypergeometric tail sum; `chi2_contingency` is the Pearson statistic with Yates correction plus the chi2 survival function; `chi2.ppf` is the inverse regularised lower incomplete gamma. All three go through Cephes `igam`/`igamc`/`igami`, so vendoring Cephes covers them |
| `scipy.stats.pearsonr`, `scipy.stats.spearmanr` | `hicPCA.py:6`, `hicCorrelate.py` | reimplement; Spearman needs the same tie-averaged ranking as `ranksums` |
| `scipy.cluster.hierarchy.linkage`, `dendrogram` | `hicCorrelate.py:151,164` (`method='complete'`), `hicMergeDomains.py` | reimplement complete-linkage agglomerative clustering and the dendrogram leaf ordering. The leaf order determines the row order of `hicCorrelate`'s heatmap, so it is observable in the image, not just internal |
| `sklearn.cluster` (`KMeans`, hierarchical, spectral) | `hicAggregateContacts.py:16,553` | `KMeans(n_clusters=k, random_state=0)` is sklearn's k-means++ with a fixed seed, 10 restarts, Elkan or Lloyd depending on data. Reproducing it exactly means reproducing sklearn's RNG stream (`check_random_state(0)` -> numpy `RandomState` Mersenne Twister) and its k-means++ candidate selection. **Decision: reimplement k-means++ against numpy's `RandomState(0)` stream**, which is portable (MT19937 is fully specified) and is the only sklearn algorithm the port needs. Hierarchical and spectral clustering in the same tool fall back to the `scipy.cluster.hierarchy` reimplementation and to a dense eigendecomposition |
| `scipy.ndimage.rotate` | `hicPlotAverageRegions.py:15` | spline-interpolated affine rotation; part of the plotting shell (tier 7), so it stays in Python |
| `graphviz.Digraph` | `hicMergeDomains.py` (`create_tree`:267) | emit the DOT source directly (it is a text format) and shell out to `dot` only when a rendered image is requested. `python-graphviz` itself does exactly that |
| `Bio.SeqIO`, `Bio.Seq` | `hicFindRestSite.py`, `buildMatrixMethods.py:27` | reimplement: a streaming FASTA reader (plain and gzip) and reverse-complement with IUPAC codes. `hicFindRestSite` also shells out to the external `sort` binary (`hicFindRestSite.py:115,121`); the C++ version sorts in memory, which changes nothing observable as long as the comparison key matches GNU `sort`'s default byte order under `LC_ALL=C` |
| `unidecode` | `utilities.py:702 remove_non_ascii` | only used for QC report text; a 30-line ASCII-fold table covers the cases that occur |
| `scipy.signal` | none found | not a dependency |
| `matplotlib`, `pygenometracks` | 8 plotting tools plus the plot step of `hicAggregateContacts`, `hicPlotSVL` and `hicCorrectMatrix --diagnostic_plot`; also `utilities.py:7-8` calls `matplotlib.use('Agg')` at package import, so the backend is fixed for every tool | tier 7 of section 6 |
| `imblearn`, `cleanlab`, most of `sklearn` | `lib/tadClassifier.py` | tier 8 of section 6. `sklearn.cluster` in `hicAggregateContacts` is separate and is reimplemented, see the row above |
| `hyperopt` | `hicHyperoptDetectLoops`, `hicHyperoptDetectLoopsHiCCUPS` | tier 8 of section 6 |
| `Bio.Seq` | `buildMatrixMethods.py:27` (reverse complement of restriction sequences) | 20 lines, reimplement |
| `intervaltree` | `hicmatrix`, `buildMatrixMethods.py:25` | replaced by `BinIndex` (section 2.2); `buildMatrixMethods` needs a real interval tree for restriction fragments, use a sorted-vector + binary search since fragments are non-overlapping |
| `pandas` | many tools, mostly for TSV/bedgraph IO and cooler frames | no dependency; the port writes those files directly. Where pandas' float formatting is observable (`to_csv` default `repr`-shortest), the writer must use shortest-round-trip formatting (`std::to_chars`), section 5.6 |
| `hic2cool` | `hicConvertFormat` (`hic` input) | section 3.6 |

### 3.6 `.hic` input

`hic2cool` 1.0.1 is pure Python that parses the Juicer `.hic` binary and writes a
cooler. It is reachable only through `hicConvertFormat --inputFormat hic
--outputFormat cool` (`hicConvertFormat.py:124-138`) and is the only path to
`.hic`. Test data: `test_data/hicHyperoptDetectLoopsHiCCUPS/SRR1791297_30.hic`
(5.4 MB).

**Decision: port the `.hic` reader in C++ (tier 4), not earlier.** The format is a
documented binary layout (magic `HIC\0`, version 8 or 9, master index, block
compression with zlib) and reading it is a few hundred lines. Do not shell out to
Python. Until it is ported, `hicConvertFormat` reports "hic input not yet
supported" and `STATUS.md` records the gap; it must not silently succeed.

## 4. Threading and memory

32 cores are available. The Python code parallelises with `multiprocessing`
(`Process` + `Queue`, and `multiprocessing.sharedctypes.RawArray`), which forces
serialisation of results and, for `hicDetectLoops` and `hicBuildMatrix`, copies
of large arrays per worker.

### 4.1 Threading model

- **One process, `std::jthread` + a small fixed thread pool** in
  `core/include/hicx/parallel.hpp`. No OpenMP: the env's `libgomp` is already
  loaded by krbalancing, and mixing an OpenMP runtime with an explicit pool
  invites oversubscription. Where krbalancing's own OpenMP loops run, set
  `omp_set_num_threads` from the tool's `--threads` value.
- **Determinism is a hard requirement.** Every parallel reduction must produce
  the same result for any thread count. The rule: partition by a *fixed* index
  range (chromosome, diagonal, bin1 chunk), reduce within a partition
  sequentially, and combine partitions in index order. Never accumulate into a
  shared float. This is checkable and is checked: the harness runs every
  threaded tool at `--threads 1` and `--threads 16` and requires byte-identical
  output (contract rule: the falcoAmadeus precedent, `-tN == -t1`).
- Per-tool `--threads` defaults must match the Python defaults exactly, since
  some Python tools change *results* with thread count (a bug the port should
  not inherit; where it exists, record it).

### 4.2 Memory

The known blowups, measured on this machine against the reference env:

| run | wall | peak RSS |
|---|---|---|
| `hicInfo -m Li_et_al_2015.h5` | 0.70 s | 264 MB |
| `hicInfo -m Li_et_al_2015.cool` | 0.43 s | 139 MB |
| `hicConvertFormat Li h5 -> cool` | 1.31 s | 390 MB |
| `hicTransform --method obs_exp Li_et_al_2015.h5` | 1.32 s | 305 MB |
| `hicCorrectMatrix correct --correctionMethod KR` Li | 1.31 s | 486 MB |
| `hicCorrectMatrix correct --correctionMethod ICE` Li | 13.59 s | 396 MB |
| `hicTransform --method pearson Li_et_al_2015.h5` | 26.45 s | **5.15 GB** |

`Li_et_al_2015.h5` is 11,104 bins on one chromosome with 3,313,107 nonzeros
(2.7 % dense). The pearson path allocates the dense per-chromosome matrix
(11,104^2 x 8 B = 987 MB), then `np.corrcoef` of it, then a `lil_matrix`
accumulator, then a CSR copy: five to six live copies, and it writes a 370 MB
`.h5`. Scaled to `test_data/hicTADClassifier/gm12878_chr1.cool` (24,926 bins,
61.8 M nonzeros) the same path needs about 5 GB per dense copy, so over 25 GB.
That is the memory bug to fix, and it is fixable without changing results:

1. Hold the dense per-chromosome block **once**, column-major, in a single
   `std::vector<double>`, and compute Pearson in place (centre columns, scale by
   the column norms, then one `dsyrk`). One dense copy instead of five.
2. Stream the result straight into the output CSR: for a Pearson matrix the
   result is dense, so build `indptr`/`indices` analytically and fill `data`
   row by row rather than going through a LIL accumulator.
3. Never call `todense()` on a whole-genome matrix. Every `todense()` in the
   Python source is per chromosome; keep it that way and assert it.

Note that fixing this changes the reduction order relative to `np.corrcoef`, so
the Pearson and covariance outputs move from bit-identical to float-tolerant.
That is an accepted, recorded trade (section 5.3).

Other memory rules:

- Read cool pixels in chunks, as the Python does (`nbins//32` bin1-chunks,
  `cool.py:65-105`), into arrays preallocated to the `nnz` attribute. Never
  materialise a pandas-like frame.
- `mmap` is not used for HDF5; rely on the HDF5 chunk cache, sized explicitly
  with `H5Pset_chunk_cache` to at least a chunk row.
- The tool binaries must report peak RSS on `--verbose` so the perf harness
  does not have to guess.

## 5. Numeric equivalence and tolerance policy

Every tool is assigned exactly one class. The class is recorded in `STATUS.md`
and enforced by the harness. There is no per-run tolerance tuning.

### 5.1 The classes

| class | criterion | how measured |
|---|---|---|
| **E0 exact** | byte-identical output file | `cmp` |
| **E1 structural** | HDF5 objects, dtypes, shapes, chunking, filters equal; every dataset decodes to bit-identical bytes; attributes equal after normalising `creation-date`, `generated-by`, `generated-by-cooler-lib`, `tool-url` | cool comparator |
| **E2 value-exact** | same sparsity pattern; every stored value bit-identical (`memcmp` of the decoded arrays); integer fields exactly equal | h5 and text comparators |
| **E3 tight float** | same sparsity pattern exactly; `max |a-b| / max(1, |b|) <= 1e-12` over all nonzeros; integer and string fields exactly equal | numeric comparator |
| **E4 loose float** | same sparsity pattern to within 0.1 % of nonzeros; `max |a-b| / max(1e-6, |b|) <= 1e-6` over the common support; Pearson correlation of the two value vectors `>= 1 - 1e-9` | numeric comparator |
| **E5 set agreement** | for tools whose output is a set of called regions: Jaccard index of the called intervals `>= 0.99`, and every disagreeing call has a score within 1 % of its threshold | interval comparator |
| **E6 visual** | image comparison, RMS difference over the pixel array `<= 5` on a 0-255 scale, same dimensions | image comparator |
| **E7 not equivalent** | deliberate deviation, documented | `STATUS.md` note, no automated check |

### 5.2 Why float classes exist at all

IEEE 754 addition is not associative, so a sum's value depends on the order of
accumulation. numpy does not sum naively: `np.add.reduce` over a contiguous
float64 array uses **pairwise summation** with an 8-way unrolled inner block
(`numpy/_core/src/umath/loops_utils.h`), which gives an error bound of
O(log n) eps rather than O(n) eps. A naive C++ `for` loop gives a different,
usually worse, result. Concretely, for the sum of 3.3 M float64 values in
`Li_et_al_2015.h5` (`hicInfo` prints `17548966.536917936`), naive left-to-right
summation and numpy pairwise summation differ in roughly the last 3-4 decimal
digits.

Consequences that shape the policy:

- Any tool whose output is a **sum, mean, variance or correlation over more
  than a few hundred float values** cannot be bit-identical unless the C++
  reproduces numpy's pairwise blocking. For `hicInfo` (a single global sum that
  is printed and compared as text) the port **must** implement numpy's pairwise
  algorithm, because the number is compared as a decimal string. That is about
  40 lines and it is worth it once, in `core/src/math/reduce.hpp`, and reused
  everywhere a whole-array reduction appears.
- Any tool whose float work is **elementwise** (scaling, division, obs/exp with
  a precomputed divisor) is bit-identical for free, because elementwise IEEE
  operations are exactly specified. Most of the matrix arithmetic is in this
  category.
- Any tool that calls **LAPACK or a transcendental** inherits the library's
  results. Using the oracle's own OpenBLAS and Cephes sources moves these back
  into the exact category; using a different implementation does not.

So the policy is: **implement numpy-pairwise reduction and vendor Cephes,
and then most tools become E2 rather than E3.** The float classes are for the
places where an algorithmic change is deliberately made (section 4.2) or where
an iterative solver's stopping point differs.

### 5.3 Assignment rationale, by kind of tool

- **Format conversion and matrix arithmetic on integers** (`hicConvertFormat`,
  `hicSumMatrices`, `hicAdjustMatrix`, `hicMergeMatrixBins`, `hicMergeTADbins`,
  `hicCompareMatrices` in `diff`/`ratio` on integer input): **E1 for cool, E2 for
  h5, E0 for text formats.** These are permutations, selections and integer sums.
  There is no float rounding to hide behind. If they are not exact, something is
  wrong.
- **`hicInfo`, `hicQuickQC`, QC counters**: **E0.** They are text, the integers
  are exact, and the one float (`sum of matrix`) is exact once pairwise reduction
  is in place. Note `hicInfo.py:126` computes
  `((sum - diag.sum())/2) + diag.sum()`, so the reduction order is fixed by that
  expression and must be reproduced literally.
- **ICE correction** (`hicCorrectMatrix --correctionMethod ICE`): **E3.** The
  loop at `iterativeCorrection.py:40-71` is 50 elementwise passes with a
  convergence test on `max|s-1| < 1e-5`. The per-pass marginal `W.sum(axis=1)` is
  a sparse row reduction whose order the port can match exactly (CSR row order),
  so in principle E2; but the loop is chaotic in the sense that a single-ulp
  difference in a marginal can change the iteration count near the tolerance
  boundary, which changes the result by up to one full pass. E3 with a
  `1e-12` relative bound is the honest class. The harness additionally asserts
  that the **iteration count is equal**, which is the real check.
- **KR correction**: **E2**, because the same C++ source runs (section 3.3).
  If the vendoring fallback is used instead, drop to E3.
- **obs/exp, z-score, normalisation, `hicNormalize`**: **E2.** Per-diagonal
  means are reductions over a few thousand values at most; implement them with
  the pairwise reducer and they are exact.
- **Pearson and covariance matrices** (`hicTransform --method pearson|covariance`,
  `hicPCA --pearsonMatrix`): **E3.** The memory rewrite in 4.2 changes the
  reduction order of the centring and the inner products. `1e-12` relative is
  achievable for correlations in `[-1,1]`; anything looser would hide a real
  error.
- **`hicPCA` eigenvectors**: **E3 on |value|, plus an explicit sign and order
  check.** See 5.4.
- **TAD calling, loop calling, differential tests** (`hicFindTADs`,
  `hicDetectLoops`, `hicDifferentialTAD`, `hicMergeDomains`,
  `chicSignificantInteractions`): **E5.** These produce a set of called
  regions after thresholding a p-value. A last-ulp difference in a test statistic
  can flip a call that sits exactly on the threshold, and no float tolerance
  expresses that correctly. The Jaccard floor of 0.99 with the "every
  disagreement is within 1 % of threshold" side condition does. In addition, the
  underlying continuous outputs (`--tadScore`/`.bm` bedgraph matrix for
  `hicFindTADs`, the raw p-value columns for the others) are compared at **E3**,
  which is where a real numeric regression would show up.
- **Tools calling `fit_nbinom`** (`hicDetectLoops`, `chicViewpointBackgroundModel`):
  the fitted `(r, p)` parameters are compared at **E4** (`1e-6` relative), because
  L-BFGS-B stops on its own gradient criterion and a different implementation
  stops at a different point. Downstream calls stay at E5.
- **Plot data files** (the `.tab`/`.bedgraph`/`.txt` a plotting tool writes
  alongside its image): **E0 or E3** depending on whether they carry floats; the
  image itself is **E6**.
- **ML tools**: **E7.** Section 6.7.

### 5.4 The `hicPCA` eigenvector problem

`hicPCA.py:305` calls `scipy.linalg.eig(corrmatrix)`, the **general
non-symmetric** solver (LAPACK `dgeev`), on a symmetric covariance matrix. It
does not sort the result. `--whichEigenvectors 1 2` then takes columns 0 and 1
of whatever order `dgeev` returned (`hicPCA.py:314-322`). Two consequences:

1. The **column order** is LAPACK's, which depends on the LAPACK implementation
   and its blocking. Reproducing it requires calling the same `dgeev` from the
   same OpenBLAS build. This is why section 3.4 pins OpenBLAS from `$HICX_DEPS`.
2. Eigenvector **signs** are arbitrary. The Python partly compensates by
   flipping the sign to correlate positively with a gene-density or histone
   track when `--extraTrack` is given (`hicPCA.py:139-193`), but with no
   `--extraTrack` the sign is whatever LAPACK produced.

The equivalence rule for `hicPCA` is therefore: compare `|eigenvector|`
elementwise at E3; separately assert that the sign pattern is either identical
or globally flipped per chromosome; and when `--extraTrack` is supplied, require
identical signs, since the flip is then deterministic. Any other outcome is a
failure, not a tolerance.

### 5.5 How tolerance is measured, exactly

For a matrix pair `(A_cpp, B_py)`:

1. **Sparsity pattern**: the sets of `(i,j)` with a stored entry must be equal
   for E2 and E3. For E4 the symmetric difference must be at most 0.1 % of
   `max(nnz_A, nnz_B)`, and the excess entries must all satisfy
   `|value| <= 1e-9`.
2. **Max relative error**: over the common support,
   `max_ij |a_ij - b_ij| / max(atol_floor, |b_ij|)`, with `atol_floor = 1`
   for E3 and `1e-6` for E4. Using a floor rather than a pure relative error
   avoids the meaningless blow-up on values near zero.
3. **Correlation floor** (E4 only): Pearson correlation of the two value
   vectors over the common support, computed in float64 with the pairwise
   reducer.
4. **Integer and string fields**: always exactly equal, in every class except
   E7.

For text outputs the comparator parses each line into typed fields using a
per-format schema (bed, bedgraph, bedpe, tsv-with-header) and applies the same
rules field by field, so that a float column at E3 does not force the whole file
to E3. The number of lines and their order must always match exactly.

### 5.6 Float formatting

Several tools print floats into text output. Python's `str(float)` and
`repr(float)` produce the **shortest string that round-trips** (David Gay /
Ryu). `printf("%g")` does not. The port must use `std::to_chars(..., double)`
with no precision argument, which is also shortest-round-trip, and must
replicate Python's exponent formatting (`1e-05` not `1.0e-05`; `inf`, `nan`,
`-inf` spellings) in a single helper `hicx::fmt::py_repr(double)`. Where the
Python code uses an explicit format (`'{:.12f}'` in
`lib/viewpoint.py:285 writeInteractionFile`, `'{:,}'` thousands separators in
`hicInfo.py:147,163`), reproduce the format string semantics instead.

## 6. Porting order

Nine tiers. Tier 0 is the library; tiers 1 to 8 are the 46 tools, each appearing
exactly once. The ordering principle: **the file layer first, because every tool
depends on it and it is the only part where an error is silent and corrupting;
then tools that exercise the file layer without adding math, because they are
the cheapest possible validation of it; then math; then the awkward tools last,
so that their unresolved questions never block anything.**

### Tier 0 - core library (no tools)

`libhicx4`: `CutIntervals`, `BinIndex`, `CsrMatrix`, `Matrix`, the format
readers and writers (cool, mcool, scool-write, h5, homer, ginteractions,
hicpro, 2D-text), the pairwise reducer, Cephes, the BED/narrowPeak reader,
`reduce_matrix`, and the CLI/argparse compatibility layer.

The **argparse compatibility layer** is not optional. Every tool's help text,
argument names, short options, `choices`, defaults, `nargs`, `metavar` and
error messages are part of the interface and are asserted by the Python tests
(and by any Galaxy wrapper). Build one `hicx::ArgParser` that reproduces
argparse's grouping (`Required arguments` / `Optional arguments`), its
`--help` rendering, its `%(prog)s {version}` version action and its error text.
Doing this once in tier 0 costs a day and saves 46 hand-written parsers.

Exit criteria for tier 0: round-trip every matrix in `test_data/` through
C++ read and C++ write in every format the Python supports, and have the Python
loader produce identical in-memory state; and pass the C++ unit tests for the
quirks in section 2.7.

### Tier 1 - file layer exercisers (6 tools)

Rationale: each is a thin shell over the matrix layer, so a failure localises to
the file layer, and together they cover every reader and writer.

| tool | LOC | class | why here |
|---|---|---|---|
| `hicInfo` | 197 | E0 | pure read + text; exercises cool metadata, h5, mcool, the NaN-bin heuristic |
| `hicConvertFormat` | 344 | E1/E2/E0 | the only tool that touches every format; the acceptance test for tier 0 |
| `hicSumMatrices` | 74 | E1/E2 | integer addition, `chrBinBoundaries` ordering check, `maskBins` |
| `hicCompareMatrices` | 96 | E1/E2 | elementwise diff/ratio/log2ratio |
| `hicAdjustMatrix` | 201 | E1/E2 | `keepOnlyTheseChr`, `maskBins`, `reorderBins` |
| `hicMergeMatrixBins` | 282 | E1/E2 | `reduce_matrix`, also used by `hicConvertFormat --outputFormat mcool` |

### Tier 2 - interval and text tools (7 tools)

Rationale: no matrix math at all, so they can be developed in parallel with
tier 3 by a second pair of hands, and they pin down the BED/bedgraph writers and
the bedtools replacement before any tool that depends on them.

`hicFindRestSite` (142), `hicMergeLoops` (173), `hicValidateLocations` (284),
`hicCreateThresholdFile` (54), `hicMergeTADbins` (147), `hicAverageRegions` (209),
`hicNormalize` (153). Classes: E0 for the pure-text ones, E1/E2 for
`hicMergeTADbins` and `hicNormalize` which write matrices.

### Tier 3 - float matrix math (6 tools)

`hicTransform` (260, E2 for obs_exp/norm, E3 for pearson/covariance),
`hicCorrectMatrix` (779, E3 for ICE / E2 for KR),
`hicPCA` (412, E3 plus the sign rule),
`hicCompartmentalization` (223, E3),
`hicInterIntraTAD` (514, E3),
`hicPlotSVL` (261, statistics at E3, plot deferred to tier 6).

Rationale: this is where the numeric policy is proven. `hicTransform` first
because `hicPCA`, `hicDetectLoops` and `hicFindTADs` all reuse its obs/exp and
Pearson kernels; `hicCorrectMatrix` second because it is the most-used tool in
the suite and because KR being bit-identical (section 3.3) is a strong early
signal.

### Tier 4 - alignment and matrix construction (3 tools + the `.hic` reader)

`hicBuildMatrix` (270 + `lib/buildMatrixMethods.py` 900), `hicBuildMatrixMicroC`
(226), `hicQuickQC` (127), plus the `.hic` reader that completes
`hicConvertFormat` (section 3.6).

Rationale: htslib enters here and nothing else depends on it. `hicBuildMatrix`
is the largest single unit of work in the port (BAM pair iteration, supplementary
alignment resolution, dangling-end and self-circle classification, restriction
fragment binning, the QC table). Its output matrix is integer, so **E1/E2**, and
its QC tables are **E0**. `hicQuickQC` shares the same read classification code
and is essentially a free follow-on.

### Tier 5 - TAD, loop and differential calling (5 tools)

`hicFindTADs` (1,368), `hicDetectLoops` (1,093), `hicDifferentialTAD` (518),
`hicMergeDomains` (428), `hicAggregateContacts` (976, data path only; its plot
goes to tier 6).

Rationale: the heaviest algorithms, all depending on tier 3 kernels. Classes E5
for the call sets and E3 for the continuous intermediates (section 5.3).
`hicFindTADs` additionally writes a z-score matrix and a `.bm` bedgraph matrix,
both **E3**, and those are the real regression detectors; the domain BED is E5.

### Tier 6 - cHi-C suite (7 tools)

`chicQualityControl` (307), `chicViewpointBackgroundModel` (257),
`chicViewpoint` (312), `chicSignificantInteractions` (650),
`chicAggregateStatistic` (391), `chicDifferentialTest` (419) and
`chicExportData` (570). `chicPlotViewpoint` belongs to this pipeline but is a
plotting tool and is counted in tier 7.

Rationale: a self-contained pipeline with its own HDF5 schema (written with
h5py directly by `lib/viewpoint.py:301-367`, not through hicmatrix) and its own
shared library (`lib/viewpoint.py` 1,100 lines, `lib/cnb.py`). It depends on
tier 0 and on the NB machinery from tier 5 but on nothing else, so it can be
done by a separate worker in parallel with tier 5. The intermediate `.hdf5`
files are **E1**; the exported text is **E0**/**E3**; `chicExportData`'s bigWig
output is **E1** via libBigWig.

### Tier 7 - plotting (8 tools)

`hicPlotMatrix` (1,114), `hicPlotTADs` (9, a shim over pyGenomeTracks),
`hicPlotViewpoint` (187), `hicPlotAverageRegions` (117),
`hicPlotDistVsCounts` (554), `hicCorrelate` (450), `hicPrepareQCreport`/`hicQC`
(341) and `chicPlotViewpoint` (524). `hicAggregateContacts` is counted in
tier 5; only its plotting step follows the rule below.

There is no honest pure-C++ equivalent. matplotlib's renderer, its default
style sheet, its tick locators, its colormaps and its text layout via FreeType
are a large, undocumented, version-pinned surface. Reproducing a matplotlib PNG
to within the tests' own tolerance (`compare_images` with tolerances up to 40)
would take longer than the rest of the port combined, and reproducing
pyGenomeTracks (which `hicPlotTADs` simply delegates to, `hicPlotTADs.py:4-9`)
is out of scope by any measure.

The four options, and the recommendation per tool:

- **(a) Thin Python shell over the C++ core.** The C++ binary computes and emits
  the plot's data (already a real output for most of these: `.tab` for
  `hicAggregateContacts`, `.bedgraph`/`.txt` for the viewpoint tools,
  `QC_table.txt`/`distance_table.txt` for `hicPrepareQCreport`), and a small
  Python script does only the drawing, using the *same* matplotlib calls as
  today. Cost: a Python runtime dependency for plotting only. Benefit: the
  images stay bit-comparable to the current masters.
- **(b) C++ SVG/PNG backend.** Feasible for the simple line plots
  (`hicPlotViewpoint`, `hicPlotAverageRegions`, `hicPlotDistVsCounts`), not for
  the heatmaps with colorbars and genomic axes, and not at all for
  `hicPlotTADs`.
- **(c) Defer.** Ship v4 without the tool.
- **(d) Drop.**

Recommendation, per tool:

| tool | recommendation | reason |
|---|---|---|
| `hicPlotMatrix` | (a) shell | heatmap + colorbar + optional bigWig track + `--perChromosome` layout; the compute (region extraction, log transform, obs/exp) moves to C++, the draw stays matplotlib |
| `hicPlotTADs` | (a) shell, unchanged | it is already a 9-line delegation to `pygenometracks.plotTracks`; there is nothing to port. Keep it as a Python script that calls pyGenomeTracks |
| `hicPlotViewpoint` | (a) shell, with (b) as a later option | simple line plot; C++ emits the data, matplotlib draws. A native SVG backend is a reasonable v4.1 follow-up |
| `hicPlotAverageRegions` | (a) shell | small imshow |
| `hicPlotDistVsCounts` | (a) shell | the interesting part (distance-vs-count reduction, per-chromosome fits) is compute and moves to C++ |
| `hicCorrelate` | (a) shell | the correlation matrix is C++ at E3; the heatmap and scatter stay matplotlib |
| `hicPrepareQCreport` / `hicQC` | (a) shell | it is pandas table aggregation plus five bar charts and a Jinja2 HTML template; port the table aggregation to C++ at E0 and keep the rendering in Python |
| `chicPlotViewpoint` | (a) shell | as `hicPlotViewpoint` |
| `hicAggregateContacts` | split: data path C++ (tier 5, E3 on the `.tab`), plot (a) shell | the clustering (`sklearn.cluster` k-means and hierarchical, `hicAggregateContacts.py:16`) also has to move; see tier 8 note |

So: **(a) for all nine.** The v4 deliverable is a C++ core plus a thin,
explicitly-declared Python plotting shell, and `STATUS.md` records every one of
these as "not a pure C++ port" so that no one later mistakes them for done.
The equivalence class for the images stays **E6** against the existing masters
in `test_data/`, which is exactly what the current Python tests already do.

### Tier 8 - machine learning and hyperparameter search (4 tools)

`hicTADClassifier` (86 + `lib/tadClassifier.py` 800), `hicTrainTADClassifier`
(218), `hicHyperoptDetectLoops` (152), `hicHyperoptDetectLoopsHiCCUPS` (158).

The situation is different for each:

- **`hicTrainTADClassifier`**: trains an `imblearn.EasyEnsembleClassifier`
  (a bagging ensemble of AdaBoost over randomly undersampled folds), optionally
  wrapped in `cleanlab.classification.CleanLearning`, with resampling chosen
  from `imblearn.under_sampling`. It pickles the result to a `.BIN`
  (`test_data/hicTADClassifier/trained_model.BIN`, 683 KB). There is no C++
  equivalent of that training stack and no serialisation format in common.
  **Recommendation: keep as Python (option (a)), with feature extraction moved
  into the C++ core and exposed through a small binding.** The features are
  matrix rows and obs/exp windows, which is exactly what the core computes.
  Class **E7**.
- **`hicTADClassifier`**: inference from an existing pickled model. ONNX Runtime
  was the obvious candidate and it does not work: `skl2onnx` has no converter
  for `imblearn.ensemble.EasyEnsembleClassifier` or for `cleanlab`'s
  `CleanLearning` wrapper, and the shipped `.BIN` models are pickles of exactly
  those classes. Converting them would mean unwrapping each ensemble into its
  constituent `AdaBoostClassifier`s and re-expressing the resampling-aware
  prediction rule by hand, at which point the "already-trained model" guarantee
  is gone. **Recommendation: keep as Python (option (a)), same shell as
  training.** Revisit ONNX only if the models are ever retrained as plain
  sklearn estimators. Class **E7**.
- **`hicHyperoptDetectLoops`**: a driver that calls `hicDetectLoops.main` and
  `hicValidateLocations.main` in a loop under hyperopt's TPE optimiser
  (`hicHyperoptDetectLoops.py:11-12,89,117`). Once tiers 2 and 5 exist, the
  inner tools are C++ binaries. What remains is TPE. **Recommendation: port the
  driver to C++ and implement TPE.** Hyperopt's `tpe.suggest` is a documented
  algorithm (Bergstra et al.), about 400 lines: an adaptive Parzen estimator
  over the observation history, `n_startup_jobs=20` random trials,
  `gamma=0.25`, 24 candidate draws per iteration. But it is **seeded**, and the
  search path depends on the RNG stream, so results will differ. Class **E7**
  for the chosen hyperparameters; the harness instead checks that the C++ driver,
  given a *fixed* parameter set, produces the same loop file as the Python driver
  with the same set, which is a real check of the plumbing. Defer TPE itself to
  v4.1 and ship a `--parameterFile` mode first.
- **`hicHyperoptDetectLoopsHiCCUPS`**: the same driver, except the inner tool is
  `java -jar juicer.jar hiccups` invoked through `subprocess`
  (`hicHyperoptDetectLoopsHiCCUPS.py:120-127`). Nothing about it is portable and
  nothing about it needs to be: it is process orchestration. **Recommendation:
  port the driver to C++ (it still shells out to java), same TPE caveat.**
  Class **E7**.

### Tier summary

| tier | contents | tools |
|---|---|---|
| 0 | core library | 0 |
| 1 | file layer exercisers | 6 |
| 2 | interval and text tools | 7 |
| 3 | float matrix math | 6 |
| 4 | alignment and matrix construction | 3 (+ the `.hic` reader) |
| 5 | TAD, loop, differential calling | 5 |
| 6 | cHi-C suite | 7 |
| 7 | plotting (Python shell over C++ core) | 8 |
| 8 | ML and hyperparameter search | 4 |
| | **total** | **46 + 1 alias (`hicQC`)** |

Tiers 2, 6 and 7 are independent of tiers 3-5 after tier 0 lands, so up to three
workers can proceed in parallel from that point.

## 7. The state of the Python test suite, honestly

The suite is not a usable oracle as it stands, and the port cannot rely on
"the Python tests pass" as evidence of anything. The numbers, from an audit of
all 41 files in `hicexplorer/test/general/`, all 13 in `trivial_runs/`, and the
two files at the test root:

- **588 items collected** for `hicexplorer/test`
  (`test_pytest_collected_items.py` compares against 605 from
  `test_data/number_of_tests.txt`, which counts the whole package with
  `--doctest-modules`).
- **269 of those 588 (46 %) are `trivial_runs` smoke cases that assert nothing**
  beyond "the tool did not raise". One of them,
  `trivial_runs/test_hicAggregateContacts_trivial_runs_three.py:64-67`, is
  parametrized over four axes into 72 cases whose arguments are hard-coded, so
  it runs the identical command 72 times.
- **About 76 items are `xfail`**, almost all of them
  `xfail(raises=ImageComparisonFailure)`. Every image assertion in the suite is
  xfail-ed except `hicCorrelate.test_correlate`. That means the checked-in
  master PNGs are currently unverified.
- **About 41 items are always skipped**: six `hicPlotMatrix` cases behind
  `skipif(HIGH_MEMORY=120 GB > memory)`, the `hicHyperoptDetectLoopsHiCCUPS`
  case behind `skipif(nvcc)` and `skipif(not isfile('juicer.jar'))`, and the
  `hicAggregateContacts` and `hicPlotMatrix` cases behind 2 GB and 4 GB memory
  gates.
- **Roughly 170 items can actually fail on a regression.**

Assertion defects found, each of which makes a nominally-covered tool
effectively uncovered:

| location | defect |
|---|---|
| `test_compute_function.py:7` | `pTries = 1` overwrites the caller's retry count, so every `compute(main, args, 5)` runs once. Harmless, but the retry the suite thinks it has does not exist |
| `test_hicMergeDomains.py:70,84-85,106-107` | every `are_files_equal(...)` call is missing `assert`; the result is discarded. `hicMergeDomains` is "did not crash" only |
| `test_hicInterIntraTAD.py:55` | same, `are_files_equal` without `assert` |
| `test_hicHyperoptDetectLoopsHiCCUPS.py:64` | `are_files_equal` result discarded, the call has no `assert` (the comparator itself at `:29` is sound) |
| `test_hicCorrectMatrix.py:84` | the KR/cool check is a range test `3e9 < sum//2 < 3688003604`; the elementwise comparison is commented out at :85-86 |
| `test_hicCorrectMatrix.py:106` | `nt.assert_allclose(rtol=1.0)`, a 100 % relative tolerance |
| `test_hicCorrelate.py` | both tests correlate one file with itself, so a correlation of 1.0 is structurally guaranteed regardless of the implementation |
| `test_hicPlotSVL.py:66-67` | the image comparison is commented out |
| `test_hicTransform.py:17`, `test_hicConvertFormat.py:18`, `test_hicPCA.py:64`, `test_hicAverageRegions.py` | `assert_array_almost_equal(decimal=0)`, i.e. agreement to the nearest integer |
| `test_hicPCA.py:64` | compares `np.absolute(...)` at `decimal=0`, so it is both sign-agnostic and integer-rounded |
| `test_chicViewpointBackgroundModel.py:62,74` | `are_files_equal` allows 700 and 1000 mismatching values with `eps=0.1` |
| `test_hicConvertFormat.py:38-49` | asserts an integer dtype without ever passing `--enforce_integer` |
| `test_hicConvertFormat.py:76,84,94` | ginteractions, hicpro and mcool conversions assert nothing |
| `test_hicDetectLoops.py:44` | the h5 path asserts nothing |
| `test_hicPlotDistVsCounts.py:19` | the only real assertion is that the PNG's byte size differs by less than 2000 |
| `test_hicBuildMatrix.py:17` | the output BAM is checked by byte size within 80,000 bytes |
| `test_chicQualityControl.py:66` | the `_failed_reference_points` output is compared against the `_report` reference, so that output is never checked |
| `test_hicPlotViewpoint.py:83-97` | `-i` is given a literal string, not the temp file names the assertions then read |
| `test_hicAdjustMatrix.py:157,177` | the BED path is passed positionally after `--chromosomes`, so `--regions` is not exercised there |
| `test_hicPlotDistVsCounts.py:41-61` | `--skipDiagonal` is parametrized and named but never inserted into the argument string, doubling the case count with no added coverage |
| `test_pytest_collected_items.py:43-44` | rewrites `test_data/number_of_tests.txt` on success: a test with a side effect on the repository, and a ratchet rather than a fixed expectation |

Tools with **no test file at all** (5): `hicInfo` (invoked only incidentally at
`test_hicBuildMatrix.py:229`), `hicSumMatrices`, `hicMergeTADbins`,
`hicPrepareQCreport`/`hicQC`, `hicPlotTADs`. Helper modules with no test:
`utilities.py`, `parserCommon.py`, `readBed.py` (has doctests),
`reduceMatrix.py` (has doctests), `iterativeCorrection.py`.

Tools whose test file exists but whose assertions do not constrain the numbers
(structure-only or "did not crash"), and which therefore also need a
characterization test before porting: all eight `chic*` tools (HDF5 key and
attribute echo-back only), `hicTADClassifier` and `hicTrainTADClassifier` (one
cell value and one header word), `hicMergeDomains`, `hicInterIntraTAD`,
`hicHyperoptDetectLoopsHiCCUPS`, `hicCorrelate`, `hicPlotDistVsCounts`,
`hicCompartmentalization`, and the ginteractions, hicpro and mcool paths of
`hicConvertFormat`.

CLI options never exercised by any test, per tool, are enumerated in
`STATUS.md`. The worst offenders: `hicTrainTADClassifier` (14 of 24 options
untested), `hicPlotMatrix` (12 of 33), `hicBuildMatrixMicroC` (10 of 16),
`chicPlotViewpoint` (9), `hicCorrectMatrix` (7, including `--perchr`, which is
a distinct code path in both ICE and KR).

### 7.1 What this means for the port

1. **Characterization tests come first, and they are broader than the five
   untested tools.** Contract rule 1 says a tool with no test for the code path
   being ported gets a characterization test pinning current Python behaviour on
   real data, committed before the C++ that replaces it. Given the audit, that
   applies to the 5 untested tools, the ~15 tools with non-constraining
   assertions, and every unexercised CLI option listed in `STATUS.md`.
2. **The characterization tests must assert values, not structure.** The pattern
   to use: run the Python tool on a designated real input, and compare the
   output against a checked-in reference with the equivalence-class comparator
   from section 9.3 rather than with a per-file `are_files_equal(delta=...)`.
   Each test file currently defines its own `are_files_equal` with a different
   delta (10 distinct definitions across the suite); the characterization tests
   must import one shared implementation.
3. **The existing masters are the reference, but they are not trusted until
   regenerated.** Because every image assertion is xfail-ed, the PNG masters may
   not match what today's matplotlib produces. Before tier 7 starts, run every
   image test with the xfail removed and record which masters are stale;
   regenerate them from the current Python and commit them with a note, or the
   C++ port will be validated against images no one has verified in years.
4. **`number_of_tests.txt` will go up.** Adding characterization tests raises the
   collected count; `test_pytest_collected_items.py` is a `>=` ratchet, so this
   is benign, but the file will show up in every diff. Consider pinning it once
   at the end of the characterization work rather than letting each run rewrite
   it.

A full coverage run (`pytest hicexplorer/test --cov=hicexplorer
--cov-report=json`) was started for this plan and had completed only 94 of the
~588 items after 25 minutes, so the per-line coverage numbers are not in this
revision. The gap analysis above is derived from reading all 56 test files
directly, which is the stronger evidence anyway: line coverage would count the
269 assertion-free smoke runs as coverage, and they are not.

## 8. Validation protocol on real data

Contract rule 2: hand-written 5x5 matrices are unit tests, never validation.
Validation inputs are the repository's real files. The full corpus is 532 MB in
`hicexplorer/test/test_data/`.

### 8.1 Designated validation inputs

| file | size | shape | what it validates |
|---|---|---|---|
| `Li_et_al_2015.h5` | 14,215,139 B | 11,104 bins, 1,843 bp, chrX only, 1,661,678 stored (upper-triangle) nonzeros and 3,313,107 after symmetrization, 855 NaN bins, float64 data, sum 17548966.536917936 | the h5 reader/writer, ICE, KR, obs/exp, pearson, the pairwise reduction, NaN-bin handling |
| `Li_et_al_2015.cool` | 13,001,960 B | the same matrix in cool; `hicInfo` prints 1,661,678 here and 3,313,107 for the h5, see quirk 7 in section 2.7 | the cool reader/writer and the h5-cool round trip |
| `Li_et_al_2015_twice.h5` | 14,036,532 B | | `hicSumMatrices` |
| `Li_cut.h5` | 365,040 B | | fast smoke variant of the above |
| `small_test_matrix.h5` / `.cool` | 289,027 / 172,846 B | 33,754 bins, 35,857 nnz, 15 chromosomes | multi-chromosome ordering, `chrBinBoundaries`, `keepOnlyTheseChr`, blosc chunk layout reference |
| `small_test_matrix_50kb_res.h5` / `.cool` | 111,138 / 105,170 B | | `hicMergeMatrixBins`, `hicNormalize` |
| `matrix.mcool` | 2,444,203 B | 5 resolution groups named `/0`../`/4` (the legacy layout, **not** `/resolutions/<res>`), format-version 2 | the mcool reader; note the group naming, since a reader that only understands `/resolutions/` will fail here |
| `hicTADClassifier/gm12878_chr1.cool` | 79,370,717 B | 24,926 bins, 10 kb, chr1 only, 61,804,782 nnz | the large-matrix path, memory ceilings, threading determinism |
| `hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1_chr2.cool` | 21,983,415 B | 100 kb | `hicDifferentialTAD`, `hicInterIntraTAD` |
| `hicDifferentialTAD/GSM2644947_Auxin2days-R1.100000_chr1_chr2.cool` | 28,429,077 B | 3,790 bins, 100 kb, chr1+chr2, 4,208,340 nnz | as above; the two-chromosome case |
| `hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1.h5` | 9,654,418 B | | h5 form of the same |
| `hicDetectLoops/GSE63525_GM12878_insitu_primary_2_5mb.cool` | 1,602,554 B | | `hicDetectLoops` |
| `hicCorrectMatrix/gm12878_raw_values.cool` | 1,480,300 B | | ICE and KR on raw integer counts |
| `hicCorrectMatrix/gm12878_KR.cool` | 2,130,904 B | 1,254 bins, has a `/bins/weight` column | the divisive-correction load path (`correctionFactorTable`) |
| `hicPCA/mm9_reduced_chr1.cool` | 722,538 B | 9,760 bins, 20 kb, mm9 chr1, 470,730 nnz | `hicPCA`, with `pca1.bedgraph`/`pca2.bedgraph`/`pca1.bw` as masters |
| `hicPCA/obsexp_norm.h5` | 56,085 B | 4,526 bins, 2,207 nnz | `hicCompartmentalization`, and the obs/exp reference for `hicTransform` |
| `hicAdjustMatrix/gm12878_1_2_3.cool` | 693 bins, 1 Mb, chr1-3, 223,446 nnz | | the multi-chromosome `--interIntraHandling` paths |
| `hicValidateLocations/GSM1436265_RAD21_ENCFF002EMQ_10kb.cool` | 313,762 bins @10 kb, **93 contigs**, only 7,987 nnz | | the many-contig, extremely sparse case: this is the file that will break a `BinIndex` that assumes few chromosomes or contiguous coverage |
| `cHi-C/FL-E13-5_chr1.cool` | 420,484 B | 197,196 bins, 1 kb, chr1, only 83,665 nnz | the very-sparse, very-many-bins case; the whole cHi-C tier |
| `cHi-C/MB-E10-5_chr1.cool` | 461,587 B | | the second sample for every differential cHi-C step |
| `R1_1000.bam` / `R2_1000.bam` | 48,076 / 43,831 B | 1,000 read pairs | `hicBuildMatrix` fast path |
| `small_test_R1_unsorted.bam` / `small_test_R2_unsorted.bam` | 6,135,374 / 6,139,816 B | | `hicBuildMatrix` full path, QC tables, restriction-fragment mode |
| `build_region.bam` | 1,010,780 B | | `hicBuildMatrix --region` |
| `hicBuildMatrix/DpnII.bed` | 13,704,804 B | | restriction-fragment binning |
| `hicFindRestSite/hindIII.bed` | 1,349,619 B | | `hicFindRestSite` master |
| `dekker.txt.gz` | 14,580,081 B | | homer/2D-text input to `hicConvertFormat` |
| `hicConvertFormat/test_matrix.hicpro` | 1,039,906 B | | hicpro reader |
| `hicExport/matrix.GInteractions.tsv` | 1,720,060 B | | ginteractions writer master |
| `hicHyperoptDetectLoopsHiCCUPS/SRR1791297_30.hic` | 5,441,566 B | | the `.hic` reader |
| `master_matrix_plot.png`, `master_TADs_plot*.png`, `hicAggregateContacts/master_aggregate_*.png` | 12 KB - 452 KB | | E6 image masters |

### 8.2 What is compared, per tier

- **Tier 0**: for each format F and each designated matrix M: Python reads M and
  dumps its in-memory state to a canonical JSON+binary sidecar; C++ reads M and
  dumps the same; the two must be identical (section 5.1 class E2 applied to
  each of `matrix.data`, `matrix.indices`, `matrix.indptr`, `cut_intervals`,
  `nan_bins`, `correction_factors`, `distance_counts`). Then C++ writes M' in F
  and Python reads M' and dumps again: identical. Then the file comparator for F
  runs against the Python-written original at that format's class.
- **Tiers 1-6**: run the Python tool and the C++ tool with identical arguments
  on the designated inputs, compare every produced file with the comparator
  chosen by the file's extension at the tool's declared class. All CLI option
  combinations exercised by the Python test suite must be run, plus the
  combinations identified in `STATUS.md` as untested (those need a
  characterization test written first, contract rule 1).
- **Tier 7**: as above, and the image against the checked-in master with the
  same tolerance the Python test already uses.
- **Tier 8**: as declared per tool in tier 8 of section 6.

### 8.3 When a run counts as passed

A tool moves to `equivalence: pass` in `STATUS.md` only when **all** of:

1. every declared invocation produced the declared files, with exit code 0 where
   Python exits 0 and a non-zero exit where Python exits non-zero;
2. every comparator returned pass at the declared class;
3. `--threads 1` and `--threads 16` outputs are byte-identical, for tools that
   take `--threads`;
4. the tool was run at least once under the large input (`gm12878_chr1.cool` or
   `small_test_R*_unsorted.bam`) without exceeding a peak RSS of 1.5x the Python
   peak RSS on the same input;
5. the run is recorded in the harness report with the git commit of the C++
   tree, the input file checksums, and the comparator output.

Any deviation is recorded in `STATUS.md` with its reason; a tool with a recorded
deviation is `equivalence: deviation`, never `pass`.

## 9. The equivalence harness (`cpp/scripts/`, to be implemented by the implementing agent)

### 9.1 CLI

A single Python 3 entry point, run with the contract's venv interpreter, with no
dependencies beyond what that venv already has (numpy, h5py, cooler, PIL via
matplotlib):

```
cpp/scripts/equiv.py run     [--tool NAME]... [--tier N]... [--case ID]...
                             [--cpp-bin DIR] [--py-python PATH] [--jobs N]
                             [--out DIR] [--keep-workdirs] [--update-baseline]
cpp/scripts/equiv.py compare  --format {cool,h5,text,bed,bedgraph,bedpe,tsv,image,hdf5-chic,bigwig}
                              --class {E0,E1,E2,E3,E4,E5,E6} A B
cpp/scripts/equiv.py report   [--out DIR] [--format {md,json}]
cpp/scripts/equiv.py list     [--tool NAME]
```

- `--cpp-bin` defaults to `cpp/build/tools`; `--py-python` defaults to the venv
  in the contract; `PYTHONPATH` is set to the repo root so the *repo* Python
  runs, not the installed 3.7.6 package.
- `--jobs` runs cases in parallel, each in its own temporary work directory
  under `$TMPDIR`, never under the repo.
- Exit code 0 only if every selected case passed.

### 9.2 Case definition

Cases live in `cpp/scripts/cases/<tool>.yaml` (one file per tool, parsed with a
30-line hand-rolled reader so there is no PyYAML dependency, or as JSON if that
is simpler). One case:

```yaml
- id: hicTransform.obs_exp.h5
  tool: hicTransform
  tier: 3
  args: ["-m", "{data}/Li_et_al_2015.h5", "--method", "obs_exp", "-o", "{out}/oe.h5"]
  outputs:
    - path: "{out}/oe.h5"
      format: h5
      class: E2
  threads_arg: null          # or "--threads", triggers the 1-vs-16 determinism check
  expect_exit: 0
  large: false
  notes: ""
```

`{data}` expands to `hicexplorer/test/test_data`, `{out}` to the case work
directory. The runner executes the Python tool with `{out}` = `out_py` and the
C++ tool with `{out}` = `out_cpp`, then compares pairwise.

### 9.3 Comparator plugins

Each is a module in `cpp/scripts/comparators/` exposing
`compare(path_a, path_b, cls, opts) -> Result(passed, class_met, metrics, diffs)`.

- **`cool.py`** - opens both with `h5py`. Walks the object tree; compares the set
  of groups and datasets, then for each dataset: dtype (including the ENUM member
  map for `/bins/chrom`), shape, `maxshape`, chunk shape, the filter pipeline
  (id, name, `cd_values`), the fill value, and finally the decoded bytes. Compares
  root and group attributes, normalising `creation-date`, `generated-by`,
  `generated-by-cooler-lib` and `tool-url`. For E3/E4, decodes
  `bin1_id`/`bin2_id`/`count` into a COO set and applies the section 5.5 rules.
  Handles `::/resolutions/<r>` and `::/cells/<n>` by comparing every group.
  Reports: `nnz_a`, `nnz_b`, `pattern_symdiff`, `max_abs`, `max_rel`, `corr`,
  and the first 20 differing pixels with their coordinates.
- **`h5.py`** - PyTables/h5py reader for the HiCExplorer layout. Compares the
  node set, then decodes `/matrix/{data,indices,indptr,shape}` and
  `/intervals/*` and `/nan_bins` and `/correction_factors` and
  `/distance_counts`. Reconstructs both as CSR and applies the section 5.5
  rules. Does **not** compare chunk shape or filter parameters (section 2.5,
  level L3), but does compare the filter *identity* (blosc must be blosc), so a
  file PyTables cannot read still fails.
- **`text.py`** - schema-driven. A schema names the comment prefix, the delimiter,
  whether there is a header, and per-column type (`int`, `float`, `str`).
  Compares line count and order exactly; `int`/`str` columns exactly; `float`
  columns at the case's class. Built-in schemas: `bed` (3-12 cols),
  `bedgraph` (4), `bedpe` (>=6), `tsv-header`, `plain` (E0 only, byte compare).
  Reports the first 20 differing lines with the column that differs.
- **`image.py`** - loads both PNGs with `matplotlib.image.imread`, requires equal
  dimensions, computes RMS over the 0-255 pixel array exactly as
  `matplotlib.testing.compare.compare_images` does, and passes at
  `rms <= tolerance` from the case (default 5). Writes a side-by-side diff PNG
  into the report directory on failure.
- **`chic_hdf5.py`** - the cHi-C intermediate schema written by
  `lib/viewpoint.py:301-367`. Walks the group hierarchy (matrix name / chromosome
  / reference point / gene), compares the per-viewpoint datasets, floats at the
  case's class.
- **`bigwig.py`** - reads both with `pyBigWig`, compares the chrom list, the
  header (`nBasesCovered`, `minVal`, `maxVal`, `sumData`, `sumSquared`) and all
  intervals; floats at the case's class.
- **`bam.py`** - only for `hicBuildMatrix --outBam`; compares via `pysam` on the
  header (normalising `@PG`) and every record's fields.

### 9.4 Report format

`--out DIR` receives:

- `report.json`: `{harness_version, timestamp, git_commit, host, cases: [{id,
  tool, tier, class_declared, class_met, passed, py_seconds, cpp_seconds,
  py_peak_rss_kb, cpp_peak_rss_kb, outputs: [{path, format, passed, metrics,
  diffs}], stderr_py, stderr_cpp}]}`.
- `report.md`: a per-tier table (tool, cases, passed, class met, speedup, memory
  ratio) followed by a section per failing case with the comparator's diff
  excerpt.
- `workdirs/<case-id>/` when `--keep-workdirs`, holding `out_py/`, `out_cpp/`,
  the exact command lines, and the raw stdout/stderr.

`equiv.py report` regenerates `report.md` from an existing `report.json` so the
comparison run does not have to be repeated.

### 9.5 What the harness must not do

It must not normalise anything not listed in 9.3, must not retry a failing case,
and must not have a per-case tolerance override. Tolerance lives in the class,
the class lives in the case file, and changing a class is a reviewed edit to
`STATUS.md`.

## 10. Performance and memory measurement

Every harness run records, per case and per implementation:

- wall time from `time.perf_counter()` around the subprocess,
- peak RSS from `resource.getrusage(RUSAGE_CHILDREN).ru_maxrss` deltas, or more
  robustly from `/usr/bin/time -f "%e %M"` wrapping each invocation,
- user and system CPU time, so that a "faster" result that merely burns more
  cores is visible,
- the thread count the tool was given.

The report's per-tool row carries `speedup = py_seconds / cpp_seconds` and
`memory_ratio = cpp_peak_rss / py_peak_rss`. Both are informational except for
the tier-gate in 8.3(4), which requires `memory_ratio <= 1.5` on the large
input. Measurements are only comparable when the machine is otherwise idle; the
harness records the 1-minute load average at the start of each case and marks a
case `timing_unreliable` if it exceeded 2.0.

The baseline numbers in section 4.2 were taken while the Python test suite was
running concurrently and are therefore upper bounds on the Python side; they must
be retaken on an idle machine before they are quoted anywhere.

Targets, stated so that a miss is visible rather than rationalised: 5x on the
I/O-bound tools (tier 1, dominated by HDF5 and by not paying 0.4 s of Python
interpreter start-up), 3x on tier 3, 4x on tier 4 and 5 with 16 threads, and a
memory ratio below 0.4 on `hicTransform --method pearson` (section 4.2).

## 11. Principal risks

1. **The plotting and ML tiers are not C++ ports.** Thirteen of the 46 tools
   (tiers 7 and 8) keep a Python shell. That is a third of the tool count,
   although a much smaller share of the compute. This must be stated in the
   README and in `STATUS.md`, not discovered later.
2. **`dgeev` column order in `hicPCA`.** If linking the oracle's OpenBLAS does
   not reproduce the eigenvector order, `hicPCA` cannot be made equivalent
   without changing the Python (sorting by eigenvalue), which changes its output.
   This is the most likely single point of failure in tier 3.
3. **PyTables blosc chunking.** The decision to compare h5 at value level
   (L3) rather than structurally (L2) means a C++-written `.h5` will differ from
   a Python-written one byte-wise. Any downstream consumer that checksums `.h5`
   files will notice.
4. **The h5 blosc filter must actually work.** If `H5Zregister` with the
   vendored blosc filter cannot reproduce what PyTables 3.10.1 writes closely
   enough for PyTables to read it back, the whole h5 writer is blocked. Verify
   this in the first week of tier 0, before anything else is built on it.
5. **`fit_nbinom`'s L-BFGS-B.** A different optimiser stopping point propagates
   into `hicDetectLoops` and the entire cHi-C background model. E4/E5 classes
   absorb it, but if the loop call sets diverge beyond the Jaccard floor the
   only remaining option is to vendor scipy's exact L-BFGS-B Fortran translation.
6. **Characterization tests must exist before the port.** Five tools have no
   Python test at all (`STATUS.md`). Porting one of them without first pinning
   its behaviour means the port defines the behaviour, which is the one outcome
   the contract forbids.
