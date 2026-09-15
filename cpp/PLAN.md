# HiCExplorer v4 (C++) - architecture and porting plan

Owner: supervising agent. Companion ledger: `cpp/STATUS.md`. Environment facts:
`cpp/AGENTS_CONTRACT.md`. Date of this revision: 2026-09-01.

v4 has three goals of equal standing: **numerical equivalence** with the Python
reference, **reduced peak memory**, and speed. Memory is not a side effect of
the rewrite; it is budgeted per tool and enforced by the harness (section 4).

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
one interning table, one comparison. It also removes a per-bin `std::string`,
which on the 313,762-bin `GSM1436265` matrix is the difference between 10 MB of
string headers and 1.2 MB of ids.

`chrom_names` must preserve **file order**, not sorted order. `hicSumMatrices.py:52`
compares `hic.chrBinBoundaries != hic_to_append.chrBinBoundaries` and aborts when
the order differs, so ordering is observable behaviour.

### 2.2 Chromosome boundaries and interval lookup

Python builds, per chromosome, an `intervaltree.IntervalTree` of
`Interval(start, end, bin_id)` plus an `OrderedDict chrBinBoundaries[chrom] =
(first_bin, last_bin_exclusive)` (`hicmatrix/HiCMatrix.py:986-1020`). All queries
are half-open **point** lookups (`tree[pos:pos+1]`, sorted, `[0].data`;
`HiCMatrix.py:261-262`). No general interval-overlap query is ever used.

Therefore an interval tree is the wrong data structure for the port, in time and
in memory: `intervaltree` allocates a Python object per bin. Bins within a
chromosome are contiguous and sorted, so:

```cpp
struct BinIndex {
  // chrom_id -> [first_bin, last_bin_exclusive)
  std::vector<std::pair<int64_t,int64_t>> chrom_range;
  int64_t bin_at(uint32_t chrom, int64_t pos) const;  // binary search on start[]
};
```

`bin_at` is `std::upper_bound(start.begin()+lo, start.begin()+hi, pos) - 1`,
then a check that `pos < end[bin]`. O(log n) with no allocation, and it removes
the intervaltree dependency entirely. Two behaviours must be kept:

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

The C++ port keeps the **upper triangle in storage** and exposes symmetric
access, which halves the resident matrix for every tool that does not genuinely
need both triangles (section 4.4 rule 2).

```cpp
template <class T>
struct CsrMatrix {          // sorted indices, no explicit zeros
  int64_t n = 0;            // square
  bool    upper_only = true;// storage holds triu(k=0); access is symmetric
  std::vector<int64_t> indptr;   // n+1
  std::vector<int32_t> indices;  // widened to int64 only when n > INT32_MAX
  std::vector<T>       data;
};
using CsrI32 = CsrMatrix<int32_t>;
using CsrF64 = CsrMatrix<double>;
```

Rules the port must honour, because they are observable:

1. **Duplicate pixel accumulation.** `csr_matrix((data,(i,j)))` sums duplicates
   (`cool.py:95`, `hicpro.py:34`). The builder must sum, not overwrite, and must
   do so by sorting in place on a packed 64-bit `(row,col)` key and coalescing
   forward, never through a triplet vector or a map (section 4.4 rule 4).
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

Every reader decodes HDF5 chunks straight into the destination CSR arrays,
preallocated from the `nnz` attribute; every writer emits chunk by chunk
(section 4.4 rules 1 and 7).

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
| npz | yes | yes | not a matrix format in the hicmatrix sense: `scipy.sparse.save_npz`/`load_npz` written by `hicAverageRegions.py:209` and read by `hicPlotAverageRegions.py:15`. It is a ZIP container holding `.npy` arrays (`format`, `shape`, `data`, `indices`, `indptr` for CSR). Implement a minimal `.npy`/`.npz` reader and writer in `core/src/io/npz.cpp`: stored (uncompressed) ZIP entries, `.npy` v1.0 header with a Python dict literal. Small, self-contained, and the only way `hicAverageRegions` output stays readable by the existing Python plotting tool |

#### cool

Groups and dtypes as written by cooler 0.10.2 through hicmatrix:

```
/chroms/name     S<w>  (fixed-width bytes, width = max name length)
/chroms/length   int32
/bins/chrom      HDF5 ENUM over int32   <- categorical, not a string column
/bins/start      int32
/bins/end        int32
/bins/weight     float64   (only when correction factors exist; see below)
/pixels/bin1_id  int32     (hicmatrix override of cooler's int64 default)
/pixels/bin2_id  int32
/pixels/count    int32, or the matrix dtype when non-integer
/indexes/chrom_offset int64, length nchroms+1
/indexes/bin1_offset  int64, length nbins+1
```

Pixel datasets are created resizable: `shape=(min(5*nbins, nnz),)`,
`maxshape=(nnz,)`. Filters: `gzip` level 6 with `shuffle` on
(cooler defaults; hicmatrix passes no `h5opts`).

**`/bins/weight` is the exception and does not follow the other bins columns.**
An earlier revision of this plan said float32 with the same filters; that was
wrong on three counts. Verified by writing a fresh cool with
`hicCorrectMatrix correct --correctionMethod KR -m gm12878_raw_values.cool` and
dumping it, and confirmed against the checked-in
`hicCorrectMatrix/gm12878_KR.cool`:

| dataset | dtype | maxshape | chunk | filters |
|---|---|---|---|---|
| `/bins/start`, `/bins/end`, `/bins/chrom` | int32, int32, ENUM | **fixed** `(nbins)` | `nbins` | shuffle + deflate 6 |
| `/bins/weight` | **float64** | **`H5S_UNLIMITED`** | `nbins` | **deflate 6, no shuffle** |
| `/pixels/*` | int32 / count dtype | `(nnz)` | cooler's own | shuffle + deflate 6 |

The cause is that the weight column is appended through `cooler.core.put` rather
than being created with the `h5opts` the bins table is built with, so it inherits
h5py's defaults for a resizable dataset and the pandas column's float64. Note in
particular that `dtype_pixel['weight'] = np.float32` (`cool.py:303` and `:344`)
never reaches the bins table at all: `dtype_pixel` is passed as `dtypes=` to
`cooler.create_cooler`, which applies it to the **pixel** table. The C++ writer
must special-case this column.

**The `sum` root attribute is order-dependent and its order is set by the pixel
chunking.** `cooler/create/_create.py:234` starts `total = 0` (a Python int) and
`:255` does `total += chunk["count"].sum()` for each chunk in sequence. Each
`chunk["count"].sum()` is a numpy pairwise reduction *within that chunk* and in
the count column dtype; the chunk results are then accumulated **sequentially**.
The number of chunks is therefore part of the observable output, and hicmatrix
fixes it: `cool.py:366-368` splits the pixel frame with
`np.array_split(matrix_data_frame, 1e4)` when `len(self.matrix.data) > 1e7`, and
otherwise passes a single `DataFrame`, which `cooler.create.create` wraps as
`iterable = (pixels,)`. So:

- at most 1e7 stored pixels: **one** numpy pairwise reduction over the whole
  count column;
- above 1e7: **exactly 10,000** chunk reductions accumulated sequentially, and
  `np.array_split` with a float `sections` argument does work in numpy 1.26.4
  (with a `DataFrame.swapaxes` FutureWarning), so this path is live, not dead.

Reproducing `sum` bit-for-bit means reproducing that split. The dtype of the
running total follows the count column (float64 stays float64; an int32 count
column promotes to int64 through the Python `0`), so it is not always float64.

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
| **L1 byte-identical** | `cmp` of the two files succeeds | homer, ginteractions, hicpro, npz, all bedgraph/bed/tsv/txt outputs, `hicInfo` text |
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
| Pearson / covariance of a dense per-chromosome matrix | `hicPCA.py:296,301` (`np.corrcoef`, `np.cov`), `hicTransform.py` | the memory blowup, section 4.3 |
| dense eigendecomposition | `hicPCA.py:305` `scipy.linalg.eig` | sections 3.4, 5.4 |
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
   by eigenvalue. See sections 5.4 and 5.8.
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
8. **`--perchr` KR emits different correction factors for `.h5` and `.cool`
   output.** See section 3.3, point 1 of the call-site notes.

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
2. Chunk-level control is the basis of the streaming-IO memory rules
   (section 4.4 rules 1, 7 and 9); `H5Pset_chunk_cache`, partial hyperslab reads
   and appends to a resizable dataset are all C-API calls.
3. `libhdf5_cpp` is a thin wrapper that throws `H5::Exception`, which mixes
   badly with the error model everything else in the port will use.

HighFive was considered and rejected: it would be a `FetchContent` dependency
that still requires dropping to the C API for the enum, the vlen strings and the
filter registration, so it adds a layer without removing one.

Blosc: link `$HICX_DEPS/lib/libblosc.so` and register the HDF5 filter in
`core/src/io/blosc_filter.c` (vendored from the Blosc project's `hdf5-blosc`,
about 200 lines, BSD). This is required both to read existing `.h5` test data
and to write `.h5` that PyTables can read.

### 3.2 Sparse linear algebra: hand-rolled kernels, Eigen only where it already is

Do **not** build the matrix type on `Eigen::SparseMatrix`. Reasons:

1. The hot kernels are not linear algebra. ICE is row/column scaling of the
   stored triple; obs/exp is a per-diagonal reduction; `reduce_matrix` is a
   grouped sum. All are trivially expressed on the CSR arrays directly, all are
   memory-bound, and all need to reproduce numpy's exact reduction order
   (section 5.2). Wrapping them in Eigen expressions makes the reduction order
   opaque, which is exactly the property that must stay explicit.
2. The file layer needs the raw `indptr`/`indices`/`data` arrays anyway, since
   h5 stores them verbatim. An `Eigen::SparseMatrix` would have to be
   round-tripped to those arrays at every boundary, which is a full copy each
   way. The krbalancing memory analysis in section 4.3 is a direct demonstration
   of what that costs.
3. `Eigen::SparseMatrix<double, ColMajor, int64_t>` spends 16 bytes per nonzero
   against 12 for `CsrF64` with int32 indices, a 33 % penalty that applies to
   every matrix in the process.
4. The one place a real sparse solver would help does not exist: no tool solves
   a sparse system.

Eigen **is** used for dense blocks in `hicPCA` and `hicTransform`, where it is
the least painful way to hold a column-major panel and hand it to LAPACK.

### 3.3 KR balancing

`krbalancing` 0.0.5 is **already C++**. `nm -DC` on the installed
`krbalancing.cpython-312-x86_64-linux-gnu.so` shows
`kr_balancing::kr_balancing(long const&, long const&, long const&,
Eigen::Ref<Matrix<long,-1,1>>, Eigen::Ref<Matrix<long,-1,1>>,
Eigen::Ref<Matrix<double,-1,1>>)`, `computeKR`, `inner_loop`, `outer_loop`,
`rescale_norm_vector`, `get_normalised_matrix(bool&)`,
`get_normalisation_vector(bool&)`, over `Eigen::SparseMatrix<double,0,long>`.
It is a pybind11 shell over an Eigen implementation.

**Decision: reimplement the algorithm in `cpp/core/src/math/kr/`, using the
upstream source (deeptools/Knight-Ruiz-Matrix-balancing-algorithm, 363 lines,
fetched and read at `scratchpad/krb/src/krbalancing.{hpp,cpp}`, confirmed to be
the 0.0.5 that matches the installed `krbalancing-0.0.5-py312h28adbb1_9`) as the
specification of the iteration, but not vendoring it verbatim.**

An earlier revision of this plan said vendoring it would make KR "bit-identical
by construction". **That was wrong**, for two independent reasons found by
reading the source and by measurement:

- The upstream code downcasts the float64 input to **float32** on ingestion
  (`krbalancing.cpp:12` `typedef Eigen::Triplet<float> T;`, `:27`
  `float(input_values(j_start))`) and then stores it into a float64 sparse
  matrix. HiCExplorer's KR balances a float32-rounded matrix today.
- `rescale_norm_vector` accumulates `original_sum` and `norm_vector_sum` as
  **`float`** (`krbalancing.cpp:228-229`) over every stored value, inside an
  `omp parallel for` whose body is wrapped in `omp critical` (`:233-250`). The
  critical section serialises the additions but does **not** fix their order, so
  a float32 sum of tens of millions of terms accumulates in a
  thread-scheduling-dependent order. **The result is that HiCExplorer's KR is not
  reproducible run to run.**

Measured, on this machine, against the reference oracle:

| input | runs | normalisation factors observed |
|---|---|---|
| `Li_et_al_2015.h5`, 11,104 bins, 1.66 M stored nonzeros | 6 | 0.0190883, 0.0190885, 0.0190887, 0.0190889, 0.0190894, 0.0190899 |
| `hicTADClassifier/gm12878_chr1.cool`, 24,926 bins, 61.8 M stored nonzeros | 3 | 0.00660838, 0.00663618, 0.00666417 |

Comparing the five `Li` output matrices pairwise: the sparsity pattern is stable
across runs, but the **maximum pairwise relative difference is 1.503e-04 on the
matrix values and 7.513e-05 on the correction factors**, and the matrix sum
ranges over 18,525,884.7 to 18,528,668.8. On the 20x larger `gm12878_chr1`
matrix the normalisation factor alone spans **8.4e-03 relative**, close to one
percent, which is what a float32 accumulator over 61.8 M terms in an arbitrary
order produces.

So there is no bit-identical target to hit. This drives the equivalence class
for KR (section 5.7) and makes KR the one place where the port ships two modes
(section 5.8).

Two further defects in the upstream source that the port must decide about:

- The whole body of the loop in `compute_normalised_matrix` (`:212-221`) is also
  inside `omp critical`, but that loop only does
  `it.valueRef() = it.value() * x.coeff(row) * x.coeff(col)`, a pure elementwise
  in-place update with no shared state. The critical section there is pure
  serialisation with no correctness role, and `num_threads` is a hardcoded
  global of 10 (`krbalancing.hpp:29`) that no caller can change.
- `outer_loop` calls `exit(0)` after 300 outer iterations
  (`krbalancing.cpp:115-119`), after printing the entire `x` vector to stdout
  (it also prints it at 100 and 200). A library that terminates the host process
  with a **success** status on non-convergence: HiCExplorer exits 0 and writes no
  output file. The port raises an error instead; this is a deliberate deviation,
  recorded in `STATUS.md`, and it is worth reporting upstream.

The call site to reproduce is `hicCorrectMatrix.py:715-732` (per chromosome) and
`:743-755` (whole matrix). Note three things about it:

1. `get_normalisation_vector(True)` is used for the whole matrix (`:753`) but
   `get_normalisation_vector(False)` per chromosome (`:731`). Since
   `rescale_norm_vector()` is only triggered by `get_normalised_matrix(True)`,
   and that is only called when the output name ends in `.h5` (`:726`),
   **`--perchr` KR emits rescaled correction factors for `.h5` output and
   unrescaled ones for `.cool` output.** The same command with a different output
   extension produces correction factors differing by the per-chromosome
   normalisation factor. This is a defect, it is observable, and compat mode
   must reproduce it. It is also on a code path the Python suite never runs.
2. `chr_submatrix.count_nonzero()` is passed as `input_nnz` (`:719`) while the
   constructor's loop is driven by `indptr`, so if the CSR holds any explicit
   zeros the `triplets.reserve(input_nnz)` under-reserves and the vector
   reallocates, transiently doubling its footprint.
3. `.indices.astype(np.int64, copy=False)` (`:722`) and
   `.data.astype(np.float64, copy=False)` (`:724`) both **do** copy, because the
   dtypes differ from the CSR's int32 and (for a raw cool) int32. `copy=False`
   only permits avoiding a copy, it does not achieve one. On the 123.6 M-nonzero
   matrix of section 4.3 those two lines cost 989 MB each.

The float32 downcast deserves one nuance: raw Hi-C bin counts are integers well
below 2^24, so `float(x)` is exact for them and the downcast is inert on the
normal input. It only changes results when KR is applied to an already-float
matrix, which the corpus does contain (`Li_et_al_2015.h5` holds float64 values
from 0.170 to 1914.015). The float32 *accumulator*, by contrast, bites on every
input and is what produces the nondeterminism measured above.

### 3.4 Dense eigendecomposition and BLAS

Link OpenBLAS from `$HICX_DEPS` (`libopenblasp-r0.3.28.so`, which is what numpy
and scipy in that env are built against) and call LAPACK directly:
`dgeev` for `hicPCA` compatibility mode, `dsyevr` with `range='I'` for its
corrected mode, `dsyrk` for the Pearson and covariance kernels. Using the *same*
OpenBLAS build as the oracle is what makes the eigenvector column order and
signs match in compatibility mode (section 5.4). Do not use Eigen's own
`EigenSolver`: it is a different algorithm and will not reproduce LAPACK's output
ordering.

### 3.5 The remaining Python-only dependencies

| Python dependency | used by | C++ decision |
|---|---|---|
| `pysam` | `lib/buildMatrixMethods.py` only | **htslib 1.21 directly.** `pysam` is a Cython wrapper over the same library. Read BAM with `sam_open`/`sam_hdr_read`/`bam_read1`, use `bam_aux_get` for `SA`/`NM`, and the `bam1_core_t` flags. Mate-pair iteration in `readBamFiles` (`buildMatrixMethods.py:458`) assumes name-sorted paired BAMs read in lockstep, which maps directly, and lets the port stream rather than buffer whole mate blocks |
| `krbalancing` | `hicCorrectMatrix` | reimplement from the upstream source, section 3.3 |
| `pyBigWig` | `hicPCA` (read `--extraTrack`, write `.bw`), `hicPlotMatrix` (read), `chicExportData` (write) | **vendor libBigWig** (the C library pyBigWig wraps; `FetchContent` from GitHub, MIT). Writing bigWig by hand is not worth it: the zoom-level and R-tree index construction is fiddly and libBigWig is the exact code the oracle runs |
| `pybedtools` | `hicMergeLoops` (`BedTool.merge`), `hicValidateLocations`, `chicSignificantInteractions:515`, `lib/tadClassifier.py` | **reimplement.** Only `BedTool(...)`, `.sort()`, `.merge()` and `.intersect()` on small in-memory interval sets are used. That is 200 lines of sort-and-sweep and avoids a `bedtools` binary dependency at runtime. Sort order must match bedtools' lexicographic chrom sort, which is *not* the matrix chrom order; this is the likeliest source of ordering diffs and gets a dedicated unit test |
| `fit_nbinom` | `hicDetectLoops:161`, `chicViewpointBackgroundModel:234` | **reimplement.** It is 60 lines: the NB log-likelihood with `gammaln`, its analytic gradient with `psi`, and `scipy.optimize.fmin_l_bfgs_b` from an initial `(r, p)`. Port the likelihood and gradient verbatim, and vendor a L-BFGS-B implementation (the original Nocedal Fortran translation, or `LBFGSpp` via FetchContent). Different L-BFGS-B stopping behaviour is the main risk; this is why `hicDetectLoops` and `chicViewpointBackgroundModel` get loose tolerances (section 5.3) |
| `scipy.special` (`gammaln`, `psi`, `betainc`) | `lib/cnb.py`, `fit_nbinom` | **vendor Cephes.** scipy's `gammaln`, `psi` and `betainc` are Cephes routines; using the same Cephes source gives bit-identical results, whereas `std::lgamma` does not agree with Cephes in the last ulp. This matters because `cnb.cdf` feeds p-value thresholds that get compared against fixed cutoffs |
| `scipy.stats.ranksums` | `hicDetectLoops`, `hicFindTADs`, `hicDifferentialTAD` | reimplement: rank transform with tie averaging, normal approximation, `erfc`-based two-sided p-value. Deterministic and easy to match exactly |
| `scipy.stats.anderson_ksamp` | `hicDetectLoops` | reimplement the Scholz-Stephens k-sample statistic and scipy's interpolation table for the p-value. Table values must be copied from scipy's source |
| `scipy.stats.fisher_exact`, `chi2_contingency`, `chi2.ppf` | `chicDifferentialTest.py:91-130` | reimplement. Fisher on a 2x2 table is the hypergeometric tail sum; `chi2_contingency` is the Pearson statistic with Yates correction plus the chi2 survival function; `chi2.ppf` is the inverse regularised lower incomplete gamma. All three go through Cephes `igam`/`igamc`/`igami`, so vendoring Cephes covers them |
| `scipy.stats.pearsonr`, `spearmanr` | `hicPCA.py:6`, `hicCorrelate.py` | reimplement; Spearman needs the same tie-averaged ranking as `ranksums` |
| `scipy.cluster.hierarchy.linkage`, `dendrogram` | `hicCorrelate.py:151,164` (`method='complete'`), `hicMergeDomains.py` | reimplement complete-linkage agglomerative clustering and the dendrogram leaf ordering. The leaf order determines the row order of `hicCorrelate`'s heatmap, so it is observable in the image, not just internal |
| `sklearn.cluster` (`KMeans`, hierarchical, spectral) | `hicAggregateContacts.py:16,553` | `KMeans(n_clusters=k, random_state=0)` is sklearn's k-means++ with a fixed seed. Reproducing it exactly means reproducing sklearn's RNG stream (`check_random_state(0)` -> numpy `RandomState` MT19937) and its k-means++ candidate selection. **Decision: reimplement k-means++ against numpy's `RandomState(0)` stream**, which is portable because MT19937 is fully specified, and is the only sklearn algorithm the port needs. Hierarchical and spectral clustering in the same tool fall back to the `scipy.cluster.hierarchy` reimplementation and to a dense eigendecomposition |
| `scipy.ndimage.rotate` | `hicPlotAverageRegions.py:15` | spline-interpolated affine rotation; part of the plotting shell (tier 7), so it stays in Python |
| `graphviz.Digraph` | `hicMergeDomains.py` (`create_tree`:267) | emit the DOT source directly (it is a text format) and shell out to `dot` only when a rendered image is requested. `python-graphviz` itself does exactly that |
| `Bio.SeqIO`, `Bio.Seq` | `hicFindRestSite.py`, `buildMatrixMethods.py:27` | reimplement: a streaming FASTA reader (plain and gzip) and reverse-complement with IUPAC codes. `hicFindRestSite` also shells out to the external `sort` binary (`hicFindRestSite.py:115,121`); the C++ version sorts in memory, which changes nothing observable as long as the comparison key matches GNU `sort`'s default byte order under `LC_ALL=C` |
| `unidecode` | `utilities.py:702 remove_non_ascii` | only used for QC report text; a 30-line ASCII-fold table covers the cases that occur |
| `intervaltree` | `hicmatrix`, `buildMatrixMethods.py:25` | replaced by `BinIndex` (section 2.2); `buildMatrixMethods` needs a real interval lookup for restriction fragments, use a sorted vector plus binary search since fragments are non-overlapping |
| `pandas` | many tools, mostly for TSV/bedgraph IO and cooler frames | no dependency; the port writes those files directly. Where pandas' float formatting is observable (`to_csv` default `repr`-shortest), the writer must use shortest-round-trip formatting (`std::to_chars`), section 5.6 |
| `matplotlib`, `pygenometracks` | 8 plotting tools plus the plot step of `hicAggregateContacts`, `hicPlotSVL` and `hicCorrectMatrix --diagnostic_plot`; also `utilities.py:7-8` calls `matplotlib.use('Agg')` at package import, so the backend is fixed for every tool | tier 7 of section 6 |
| `imblearn`, `cleanlab`, most of `sklearn` | `lib/tadClassifier.py` | tier 8 of section 6 |
| `hyperopt` | `hicHyperoptDetectLoops`, `hicHyperoptDetectLoopsHiCCUPS` | tier 8 of section 6 |
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

**Superseded 2026-09-13 by tier 9, item 9.1:** the project owner asked for native
`.hic` reading **and writing**, in its own library like coolercpp. The reader is
no longer a tier 4 item.

## 4. Threading and memory (memory is a v4 goal, not a side effect)

Reduced peak memory is a first-class objective of v4, alongside equivalence and
speed, and it is enforced rather than hoped for: every tool carries a numeric
budget in `STATUS.md`, and the equivalence harness fails a tool whose C++ peak
RSS exceeds it (sections 4.5, 8.3, 9.4). The Python reference routinely spends
five to twelve times the size of the data it operates on, and in two places that
is the difference between a run that fits on a workstation and one that does not.

32 cores are available. The Python code parallelises with `multiprocessing`
(`Process` + `Queue`, `multiprocessing.sharedctypes.RawArray`, and a `Pool` in
`hicFindTADs.py:1107`), which forces serialisation of results and, for
`hicDetectLoops` and `hicBuildMatrix`, a copy of large arrays per worker.

### 4.1 Threading model

- **One process, `std::jthread` plus a small fixed thread pool** in
  `core/include/hicx/parallel.hpp`. No OpenMP as the primary mechanism: mixing an
  OpenMP runtime with an explicit pool invites oversubscription. This is also a
  memory decision: threads share the one matrix, whereas the Python's forked
  workers each accumulate their own dirty pages (section 4.4 rule 8).
- **Determinism is a hard requirement.** Every parallel reduction must produce
  the same result for any thread count. The rule: partition by a *fixed* index
  range (chromosome, diagonal, bin1 chunk), reduce within a partition
  sequentially, and combine partitions in index order. Never accumulate into a
  shared float. This is checkable and is checked: the harness runs every
  threaded tool at `--threads 1` and `--threads 16` and requires byte-identical
  output (contract rule, the falcoAmadeus precedent: `-tN == -t1`).
  The KR measurement in section 3.3 is what this rule exists to prevent: the
  Python reference violates it and is not reproducible against itself.
- Per-tool `--threads` defaults must match the Python defaults exactly.

### 4.2 Measured baselines

All measured on this machine against the reference oracle
(`PYTHONPATH=$PWD $SP/hicx-venv/bin/python bin/<tool>`), peak RSS from
`/usr/bin/time -f "%e %M"`.

| run | wall | peak RSS |
|---|---|---|
| `hicInfo -m Li_et_al_2015.h5` | 0.70 s | 264 MB |
| `hicInfo -m Li_et_al_2015.cool` | 0.43 s | 139 MB |
| `hicConvertFormat Li h5 -> cool` | 1.31 s | 390 MB |
| `hicTransform --method obs_exp Li_et_al_2015.h5` | 1.32 s | 305 MB |
| `hicCorrectMatrix correct --correctionMethod KR` Li | 1.31 s | 486 MB |
| `hicCorrectMatrix correct --correctionMethod ICE` Li | 13.59 s | 396 MB |
| `hicTransform --method pearson Li_et_al_2015.h5` | 26.45 s | **5,150 MB** |
| `hicCorrectMatrix correct --correctionMethod KR` gm12878_chr1.cool -> .cool | 91.6 s | **9,280 MB** |
| `hicCorrectMatrix correct --correctionMethod KR` gm12878_chr1.cool -> .h5 | 106.0 s | **9,199 MB** |
| `hicCorrectMatrix correct --correctionMethod ICE` gm12878_chr1.cool -> .cool | 711.2 s | **9,064 MB** |
| `hicPCA --whichEigenvectors 1 2` mm9_reduced_chr1.cool (9,760 bins, 722 KB input) | 471.2 s | **4,070 MB** |

Working-set sizes of the inputs, for reference:

| matrix | bins | stored (upper-triangle) nonzeros | symmetric nonzeros | CSR working set `W` (float64 values, int32 indices) |
|---|---|---|---|---|
| `Li_et_al_2015.h5` | 11,104 | 1,661,678 | 3,313,107 | 19.9 MB stored / 39.8 MB symmetric |
| `hicPCA/mm9_reduced_chr1.cool` | 9,760 | 470,730 | - | 5.6 MB |
| `hicTADClassifier/gm12878_chr1.cool` | 24,926 | 61,804,782 | 123,587,194 | **741.9 MB** stored / 1,483 MB symmetric |

Note the shape of the KR and ICE rows: 9.2 and 9.1 GB for a matrix whose stored
form is 742 MB. Both correction methods cost roughly 12x their data.

### 4.3 Where the Python memory goes

Two blowups dominate, and both are fully accounted for.

**KR on `gm12878_chr1.cool`: 9,199 MB against a 742 MB working set, a factor of
12.4.** The accounting, from `hicCorrectMatrix.py:743-755` and
`krbalancing.cpp:10-49`:

| allocation | site | size |
|---|---|---|
| CSR from cool, int32 counts + int32 indices | `cool.py` load | 988 MB |
| `.data.astype(np.float64)` (the source is int32, so `copy=False` cannot help) | `hicCorrectMatrix.py:747` | 989 MB |
| `.indices.astype(np.int64)` | `hicCorrectMatrix.py:746` | 989 MB |
| `std::vector<Eigen::Triplet<float>>`, 12 B per entry | `krbalancing.cpp:13-14,26` | 1,483 MB |
| `A.reserve(input_nnz)` on `SparseMatrix<double,ColMajor,int64_t>`, 16 B per entry | `krbalancing.cpp:11` | 1,978 MB |
| `setFromTriplets`' internal transposed intermediate | `krbalancing.cpp:32` | 1,978 MB |
| the temporary from `A = A + I` | `krbalancing.cpp:47` | 1,978 MB |
| interpreter, HDF5, buffers | | ~300 MB |

Not all peak simultaneously, but enough of them do. Two details make it worse
than it looks: `triplets.clear()` at `krbalancing.cpp:35` does **not** release the
vector's capacity, so the 1,483 MB is still resident during the `A = A + I` copy
at `:47`; and the whole staging is pure waste, because the caller already holds
the matrix in CSR, the matrix is symmetric so CSR equals CSC, and the constructor
could take an `Eigen::Map` over the caller's arrays and allocate nothing at all.

A separable point about the same source: the working vectors `x`, `v`, `p`, `Z`,
`rho_km1`, `rho_km2` are declared `SparseMatrixCol`, that is n-by-1 sparse
matrices, although they are dense from the first iteration
(`x = e.sparseView()` where `e` is all ones). At 16 bytes per element against 8
that is real but small: six vectors over 24,926 bins is 2.4 MB against 1.2 MB.
The cost is **time**, not memory: `x.coeff(row, 0)` in the hot loop of
`compute_normalised_matrix` is a binary search per access, executed twice per
stored nonzero, and `A * (x.cwiseProduct(p))` in the inner loop is a general
sparse-times-sparse product where a sparse-times-dense matvec would do. The port
holds these as `std::vector<double>` for both reasons, but the memory saving is
not the argument.

**`hicTransform --method pearson` on `Li_et_al_2015.h5`: 5,150 MB against a
19.9 MB working set, a factor of 259.** The per-chromosome path densifies
(11,104^2 x 8 B = 987 MB), then `np.corrcoef` of that, then a `lil_matrix`
accumulator, then a CSR copy: five to six live copies of a 987 MB block. Scaled
to `gm12878_chr1.cool` the dense block alone is 4.97 GB, so the same path needs
over 25 GB.

Smaller but systematic contributors elsewhere:

- `corrected_matrix = lil_matrix(ma.matrix.shape)` at `hicCorrectMatrix.py:700`
  is allocated unconditionally. Empty it costs only two object arrays of `nbins`
  pointers, and on the whole-matrix path it is discarded unfilled, so it does not
  show up in the measurements above. On the `--perchr` path with `.h5` output it
  **is** filled, and LIL costs roughly 60 bytes per nonzero against 12 in CSR.
  That path is untested by the Python suite (`STATUS.md`), which is presumably
  why nobody has hit it.
- `chr_submatrix = ma.matrix[a:b, a:b]` (`hicCorrectMatrix.py:707`) copies each
  chromosome block out of the whole-genome CSR while the whole-genome CSR stays
  live.
- Every `multiprocessing.Process` worker in `hicDetectLoops` and
  `hicBuildMatrix` gets a copy-on-write fork of the parent and then writes to its
  share, so the resident set multiplies by the number of workers.

### 4.4 Cross-cutting memory rules

These apply to every tool, not only the two above. They are design constraints on
`libhicx4`, checked by code review and by the harness budget (4.5).

1. **Zero-copy ingestion.** The file layer decodes HDF5 chunks directly into the
   destination `CsrMatrix`'s `indptr`/`indices`/`data` vectors, preallocated from
   the `nnz` attribute. No triplet vector, no COO staging, no intermediate frame,
   no `astype`-style dtype copy: the reader converts while decoding, in one pass.
   Where an on-disk dtype differs from the in-memory one (int32 counts to
   float64), the conversion happens chunk by chunk into the final buffer.
2. **Upper-triangle storage is the default representation.** `Matrix` holds the
   upper triangle and a symmetry flag, and every kernel that needs a full row
   uses a symmetric access helper rather than materialising the mirror. This
   halves the resident matrix everywhere. `fillLowerTriangle` semantics
   (`HiCMatrix.py:106-120`) are preserved at the API level, not in storage.
   Tools that genuinely need both triangles materialised, and must declare it:
   `hicAdjustMatrix` (`reorderBins` permutes rows and columns independently),
   `hicTransform --method pearson|covariance` and `hicPCA` (both densify per
   chromosome anyway), `hicFindTADs` (its sliding-window cut weights read
   arbitrary off-diagonal blocks), `hicDetectLoops` (neighbourhood windows
   straddle the diagonal), and `hicAggregateContacts` (submatrix extraction
   around arbitrary bed pairs). Everything else stays triangular.
3. **Index width from the bin count.** `indices` is `int32` while
   `nbins <= INT32_MAX` and `int64` above it, decided once at load. Never
   promote to `int64` to satisfy a callee's signature; that single habit costs
   989 MB in the KR path.
4. **No intermediate staging structures anywhere.** No triplet vectors, no
   LIL-equivalent, no `std::map`-keyed accumulator. Where duplicate pixels must
   be summed (section 2.3 rule 1), sort in place on a packed 64-bit
   `(row, col)` key and coalesce forward.
5. **Streaming per chromosome.** Any tool whose Python form loops over
   chromosomes processes one chromosome's submatrix at a time and releases it
   before the next, rather than slicing every block out of a live whole-genome
   matrix. Where the whole-genome matrix is only a source of blocks, read the
   blocks from the file directly using `/indexes/bin1_offset` and never hold the
   whole thing.
6. **In-place transformation where the operation allows it.** Elementwise
   scaling, obs/exp division, correction-factor application and normalisation all
   rewrite `data` in place. Reserve a second buffer only when the sparsity
   pattern changes.
7. **Streaming output.** HDF5 datasets are written chunk by chunk as the result
   is produced. A dense result (Pearson, covariance) is emitted row by row into a
   resizable dataset instead of being assembled in memory first. This is what
   makes rule 2's exemption for the pearson path affordable: the input dense
   block is unavoidable, the output copy is not.
8. **Threads share, they do not fork.** Section 4.1's single-process pool means
   the matrix exists once, not once per worker.
9. **Explicit HDF5 chunk cache.** `H5Pset_chunk_cache` sized to one chunk row and
   no more, so the cache is a bounded, budgeted cost rather than a default that
   scales with the file.

Rule 7 supersedes what an earlier revision of this plan said about the pearson
path: the fix is not only "one dense copy instead of five" but "one dense input
block, zero dense output copies". The consequence for equivalence is unchanged:
the reduction order moves relative to `np.corrcoef`, so
`hicTransform --method pearson|covariance` and `hicPCA --pearsonMatrix` are E3,
not E2 (section 5.3). Note that none of rules 1 to 9 changes a single arithmetic
result on its own; they are representation changes. The memory workstream is
therefore **not** gated on the dual-mode question of section 5.8, and can land
first.

### 4.5 The budget, and how it is enforced

Every tool has a peak-RSS budget expressed as a formula, evaluated against the
input, and recorded per tool in `STATUS.md` together with its value on the
designated large validation input. The harness fails a tool whose C++ peak RSS
exceeds its budget (section 8.3, criterion 4).

```
W      = nnz_stored * (sizeof(value) + sizeof(index)) + (nbins+1) * 8
D      = max_chromosome_bins^2 * 8          # largest dense per-chromosome block
C      = 64 MB                              # process, HDF5, buffers, output staging
budget = alpha * W + beta * D + C
```

`W` uses the **stored** (upper-triangle) nonzero count and `sizeof(index) = 4`
while `nbins <= INT32_MAX`, because that is what rules 2 and 3 make achievable.
`C = 64 MB` is a C++ process with HDF5 loaded and a 32 MB chunk cache; the
Python equivalent is 130 to 270 MB of interpreter and imports before any data is
touched, which is itself a large part of the win on the small tools.

| tool group | alpha | beta | rationale |
|---|---|---|---|
| `hicInfo` on cool with metadata | 0 | 0 | reads attributes only, never loads the matrix |
| `hicInfo` on h5, `hicQuickQC` | 1.05 | 0 | one matrix, read and reduce |
| `hicConvertFormat`, `hicNormalize`, `hicAdjustMatrix`, `hicMergeMatrixBins`, `hicMergeTADbins`, `hicAverageRegions` | 1.3 | 0 | one matrix plus a reshaped output that cannot always be built in place |
| `hicSumMatrices`, `hicCompareMatrices` | 2.2 | 0 | two matrices live, result written into the first |
| `hicCorrectMatrix` ICE and KR | 1.2 | 0 | in-place scaling, dense n-vectors are negligible |
| `hicTransform` obs_exp, norm | 1.2 | 0 | per-diagonal reduction, in place |
| `hicTransform` pearson, covariance | 1.1 | 1.15 | one dense block plus a dsyrk workspace, streamed output |
| `hicPCA` compat mode (`dgeev`) | 1.1 | 3.2 | dgeev workspace plus scipy's complex eigenvector matrix |
| `hicPCA` corrected mode (`dsyevr`, requested vectors only) | 1.1 | 1.15 | |
| `hicFindTADs`, `hicDetectLoops`, `hicAggregateContacts` | 2.2 | 0 | full symmetric materialisation (rule 2 exemption) plus per-thread scratch |
| `hicBuildMatrix`, `hicBuildMatrixMicroC` | n/a | n/a | budget is `2 * nnz_out * 12 + threads * 64 MB + C`; the matrix does not exist on input |
| cHi-C tier 6 tools | 1.3 | 0 | viewpoint extraction touches narrow row ranges |
| everything else | 1.3 | 0 | |

Targets on the designated large inputs, against the measurements in 4.2:

| tool and input | Python peak | budget | reduction |
|---|---|---|---|
| `hicCorrectMatrix --correctionMethod KR`, gm12878_chr1.cool | 9,199 MB | 1.2 x 741.9 + 64 = **954 MB** | **9.6x** |
| `hicCorrectMatrix --correctionMethod ICE`, gm12878_chr1.cool | 9,064 MB | **954 MB** | **9.5x** |
| `hicTransform --method pearson`, Li_et_al_2015.h5 | 5,150 MB | 1.1 x 19.9 + 1.15 x 987 + 64 = **1,221 MB** | **4.2x** |
| `hicTransform --method pearson`, gm12878_chr1.cool | > 25,000 MB (est.) | 1.1 x 741.9 + 1.15 x 4,971 + 64 = **6,596 MB** | > 3.8x |
| `hicPCA` compat mode (`dgeev`), mm9_reduced_chr1.cool | 4,070 MB | 1.1 x 5.6 + 3.2 x 762 + 64 = **2,509 MB** | **1.6x** |
| `hicPCA` corrected mode (`dsyevr`), mm9_reduced_chr1.cool | 4,070 MB | 1.1 x 5.6 + 1.15 x 762 + 64 = **946 MB** | **4.3x** |
| `hicInfo`, Li_et_al_2015.h5 | 264 MB | 1.05 x 19.9 + 64 = **85 MB** | **3.1x** |
| `hicInfo`, Li_et_al_2015.cool | 139 MB | **64 MB** | **2.2x** |
| `hicConvertFormat` h5 -> cool, Li | 390 MB | 1.3 x 19.9 + 64 = **90 MB** | **4.3x** |
| `hicTransform --method obs_exp`, Li | 305 MB | 1.2 x 19.9 + 64 = **88 MB** | **3.5x** |

The KR floor deserves its own line, because it is the case that prompted the
requirement. The upper triangle of `gm12878_chr1.cool` as float64 values with
int32 indices is 741.9 MB, the ten working n-vectors are 2 MB, and the balanced
output is produced by scaling `data` in place. Everything above roughly 0.95 GB
is avoidable, so the 9,199 MB the Python spends is about 9.6 times the floor
rather than an intrinsic cost of the algorithm.

`hicPCA` also has a time problem that the port should fix in the same breath:
`scipy.linalg.eig` is the general non-symmetric solver (`hicPCA.py:305`), it
computes **all** 9,760 eigenvectors when `--whichEigenvectors` asks for 2, and it
takes **471 s and 4,070 MB on a 722 KB input file**: a factor of 727 over the
5.6 MB working set, and the worst memory-to-data ratio anywhere in the corpus.
Most of that is the dense covariance block (762 MB), `dgeev`'s copy of it, and
scipy's complex128 eigenvector matrix (9,760^2 x 16 B = 1,524 MB) holding all
9,760 vectors when two were asked for. `dsyevr` with `range='I'` on the same
symmetric matrix returns the requested two in a small fraction of the time and a
quarter of the memory. That fix changes results, so it is gated behind a mode
flag (section 5.8); compatibility mode still saves 1.6x by avoiding the
intermediate copies.

## 5. Numeric equivalence and tolerance policy

Every tool is assigned exactly one class. The class is recorded in `STATUS.md`
and enforced by the harness. There is no per-run tolerance tuning.

### 5.0 What a port must meet (project acceptance gate)

**Set by the project owner, 2026-09-01: byte identicality is not required. Every
item must agree to three significant digits, that is a relative difference of at
most `1e-3` per element.** That is class **ED** below, and it is the bar every
tool has to clear.

It is expressed as a relative tolerance rather than a number of decimal places
because the corpus spans many orders of magnitude: corrected matrices hold
values from `7.89e-06` to `0.09`, while raw count matrices reach `10^8`. Three
decimal places would pass every corrected matrix trivially and would be
unreachable on raw counts. An exact zero in the reference is matched exactly,
since a zero has no significant digits to agree with and dividing by it would
let a spurious nonzero through.

The stricter classes E0 to E4 are **not retired**. They are recorded when they
are met, for two reasons: a stricter result is a much better regression signal,
and several tools reach E0 for free once the numpy-compatible reduction
machinery is in place, so there is nothing to gain by loosening them. A case
declared stricter than ED asserts that the tool is expected to achieve it, and
a drop from E0 to ED is a real regression worth investigating even though it
still passes the gate. What changes is that no tool is *blocked* on reaching
better than ED.

Two consequences worth stating plainly:

1. **Work already done is not wasted.** `hicInfo` and both writers are already
   byte-identical or structurally identical, which is strictly stronger than
   required, so they stay declared at E0/E1/E2.
2. **ED does not rescue every tool, because the oracle is not always that
   reproducible itself.** Measured: three `hicCorrectMatrix --correctionMethod KR`
   runs on `gm12878_chr1.cool` produced normalisation factors 0.00660838,
   0.00663618 and 0.00666417, a spread of `8.4e-03` relative. **Python against
   Python fails a three-significant-digit test on that input**, so KR keeps
   class EN (section 5.7) rather than ED: the comparison is against the median
   of N reference runs within the measured envelope. On the smaller
   `Li_et_al_2015.h5` the same spread is `1.6e-04` and ED would hold. The rule
   is therefore ED everywhere except where the reference's own noise exceeds
   `1e-3`, which so far means KR only.
3. **Implementation freedom (set by the project owner, 2026-09-13).** Nothing
   has to be implemented one to one with the Python. Data structures and
   algorithms may be redesigned wherever the result is the same, judged by the
   acceptance gate and the harness. The port already relies on this, with
   upper-triangle storage, streamed writes and a single dense block reduced by
   `dsyrk` in place of numpy's copies. What stays fixed is the result, not the
   route. Two consequences follow. Pinned reference defects are still
   reproduced, because they change results. And where a result depends on the
   low-order bits of an intermediate, matching that intermediate exactly is part
   of producing the same result: hicPCA selects eigenvectors positionally from a
   spectrum whose largest eigenvalue occurs 169 times, so its covariance has to
   equal numpy's bit for bit, while nothing else about how it is computed is
   constrained.

### 5.0.1 Provenance: v4 identifies itself as v4

**Set by the project owner, 2026-09-02.** Output is byte-identical to the Python
reference everywhere except the fields that record which program wrote the file,
and those must honestly name HiCExplorer version 4. Reproducing `3.7.x` there
was a side effect of chasing byte identity and is wrong: a file written by v4
that claims to come from 3.7 misleads every later reader about which
implementation produced it, including any bug report filed against it.

The version string is `4.0.0-dev` until a release is tagged, generated once from
`core/include/hicx/version.hpp.in` and never spelled out in a tool.

Fields that carry it, all of which must be written with the v4 identity:

| output | field |
|---|---|
| text reports | header lines such as `# Matrix information file. Created with HiCExplorer's hicInfo version ...` |
| cool | root attributes `generated-by` and `tool-url`; `metadata` JSON keys `matrix-generated-by` and `matrix-generated-by-url` where the tool sets them |
| HiCExplorer h5 | any producer or version attribute the writer emits |
| QC folders | the version line in `QC.log` and the HTML report |
| `--version` | the program's own version output |

The comparators normalise exactly these fields and nothing else, the same way
`creation-date` is already normalised for cool. The normalisation must be a list
of named fields, never a pattern that could hide a real difference, and a case
where a provenance field differs in anything other than the version number still
fails. The characterization tests stay pinned to the Python output, since they
test the reference, not v4.

The same concern covers a second reproducibility defect found by the harness
determinism mode: every HDF5 object v4 writes embeds its modification time,
because `H5Pset_obj_track_times` is left at its default, so two runs a second
apart differ in bytes. v4 clears it on the file, group and dataset creation
property lists so that repeated runs are byte-identical apart from the declared
`creation-date` attribute.

### 5.1 The classes

| class | criterion | how measured |
|---|---|---|
| **ED acceptance gate** | every item agrees to three significant digits, `abs(a-b) / abs(b) <= 1e-3`; an exact zero in the reference must stay an exact zero; integer and string fields exactly equal | all numeric comparators |
| **E0 exact** | byte-identical output file | `cmp` |
| **E1 structural** | HDF5 objects, dtypes, shapes, chunking, filters equal; every dataset decodes to bit-identical bytes; attributes equal after normalising `creation-date`, `generated-by`, `generated-by-cooler-lib`, `tool-url` | cool comparator |
| **E2 value-exact** | same sparsity pattern; every stored value bit-identical (`memcmp` of the decoded arrays); integer fields exactly equal | h5 and text comparators |
| **E3 tight float** | same sparsity pattern exactly; `max abs(a-b) / max(1, abs(b)) <= 1e-12` over all nonzeros; integer and string fields exactly equal | numeric comparator |
| **E4 loose float** | same sparsity pattern to within 0.1 % of nonzeros; `max abs(a-b) / max(1e-6, abs(b)) <= 1e-6` over the common support; Pearson correlation of the two value vectors `>= 1 - 1e-9` | numeric comparator |
| **E5 set agreement** | for tools whose output is a set of called regions: Jaccard index of the called intervals `>= 0.99`, and every disagreeing call has a score within 1 % of its threshold | interval comparator |
| **E6 visual** | image comparison, RMS difference over the pixel array `<= 5` on a 0-255 scale, same dimensions | image comparator |
| **EN within oracle noise** | the Python reference is **not reproducible against itself**; the tolerance is measured, not chosen. See 5.7 | noise-envelope comparator |
| **E7 not equivalent** | deliberate deviation, documented | `STATUS.md` note, no automated check |
| **EX external reference** | only for tier 9 features, which have no Python version: agreement with a named external implementation on the same real input, and, for calling tools, recovery of signal planted into a real matrix. The criterion is fixed in this plan before the feature is implemented, never after seeing results | per-feature comparator named in tier 9 |

### 5.2 Why float classes exist at all

IEEE 754 addition is not associative, so a sum's value depends on the order of
accumulation. numpy does not sum naively, and the exact scheme has now been
reproduced and verified bit-for-bit against numpy over 161 array sizes
(`core/src/numpy_compat.cpp`, `tests/test_numpy_compat.cpp`). It has three
layers, all of which matter:

1. `np.add.reduce` processes the array through the ufunc **buffer**, whose
   default size is 8192 elements, and accumulates the per-buffer results
   **sequentially**.
2. Within a buffer, the float64 loop in
   `numpy/core/src/umath/loops_arithm_fp.dispatch.c.src` uses **pairwise
   summation**: eight accumulators inside blocks of `PW_BLOCKSIZE = 128`
   elements, splitting recursively above that. The error bound is O(log n) eps
   rather than O(n) eps.
3. A **float32** matrix sums in float32 all the way through `matrix.sum()` and
   is only promoted afterwards, so the accumulation has to happen in single
   precision to reproduce the printed value.

A naive C++ `for` loop gives a different, usually worse, result. Concretely, for
the sum of 3.3 M float64 values in `Li_et_al_2015.h5` (`hicInfo` prints
`17548966.536917936`), naive left-to-right summation and numpy's scheme differ
in roughly the last 3-4 decimal digits. Reproducing all three layers is what made
`hicInfo` byte-identical on all 184 cool and h5 matrices in the corpus.

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
places where an algorithmic change is deliberately made (section 4.4) or where
an iterative solver's stopping point differs.

**Revised 2026-09-01 by the acceptance gate in 5.0.** Half of that policy is
now optional. The split is by cost:

- **Keep the numpy-pairwise reduction.** It is already written, unit-tested
  against numpy over 161 array sizes, and used by every reduction in the core.
  It costs nothing further and it is what makes `hicInfo` and the writers land
  at E0 instead of ED. Removing it would only lose signal.
- **Drop the vendoring of Cephes, and drop pinning the oracle's OpenBLAS.**
  Those existed to make transcendentals and LAPACK bit-reproducible, and a
  different `libm` agrees with Cephes to far better than `1e-3`.
- **Except where the result depends on the bits (revised 2026-09-13).** hicPCA
  keeps LAPACK `dgeev` with positional column selection, as the Python does,
  because on real data its spectrum has a largest eigenvalue occurring 169
  times, and the column order `dgeev` returns depends on the last bits of its
  input. Its covariance is therefore computed to equal `np.cov` bit for bit. An
  earlier revision of this paragraph proposed `dsyevr` alone for hicPCA; that
  selects different eigenvectors and is withdrawn. See 5.0 item 3.
- **Use the library scipy uses.** Where scipy itself relies on a vendored
  library, the port uses the same one: scipy 1.14 computes `betainc` with
  Boost.Math `ibeta`, not Cephes, and Cephes differs by up to one ulp, enough to
  turn a stored p-value of 0.0 into `1.1e-16` (`STATUS.md` F42).

The one place ED does not help is a reference that is not reproducible against
itself. See 5.0 consequence 2 and section 5.7: KR stays at EN.

Note that the upper-triangle storage decision (4.4 rule 2) is *not* one of those
places, as long as the symmetric access helper visits elements in the same order
the full CSR would. It must, and there is a unit test for it.

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
  convergence test on `max abs(s-1) < 1e-5`. The per-pass marginal
  `W.sum(axis=1)` is a `coo_matvec` against a vector of ones, which accumulates
  in COO storage order; the COO comes from `.tocoo()` on a CSR, so that order is
  exactly CSR row-major order and the port reproduces it without a compat mode.
  In principle that makes ICE E2; but the loop is chaotic in the sense that a
  single-ulp difference in a marginal can change the iteration count near the
  tolerance boundary, which changes the result by up to one full pass. E3 with a
  `1e-12` relative bound is the honest class. The harness additionally asserts
  that the **iteration count is equal**, which is the real check.
- **KR correction**: **EN**, section 5.7.
- **obs/exp, z-score, normalisation, `hicNormalize`**: **E2.** Per-diagonal
  means are reductions over a few thousand values at most; implement them with
  the pairwise reducer and they are exact.
- **Pearson and covariance matrices** (`hicTransform --method pearson|covariance`,
  `hicPCA --pearsonMatrix`): **E3.** The memory rewrite in 4.4 changes the
  reduction order of the centring and the inner products. `1e-12` relative is
  achievable for correlations in `[-1,1]`; anything looser would hide a real
  error.
- **`hicPCA` eigenvectors**: **E3 on abs(value), plus an explicit sign and order
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
- **ML tools**: **E7.** Tier 8 of section 6.

### 5.4 The `hicPCA` eigenvector problem

`hicPCA.py:305` calls `scipy.linalg.eig(corrmatrix)`, the **general
non-symmetric** solver (LAPACK `dgeev`), on a symmetric covariance matrix. It
does not sort the result. `--whichEigenvectors 1 2` then takes columns 0 and 1
of whatever order `dgeev` returned (`hicPCA.py:314-322`). Three consequences:

1. The **column order** is LAPACK's, which depends on the LAPACK implementation
   and its blocking. Reproducing it requires calling the same `dgeev` from the
   same OpenBLAS build. This is why section 3.4 pins OpenBLAS from `$HICX_DEPS`.
2. Eigenvector **signs** are arbitrary. The Python partly compensates by
   flipping the sign to correlate positively with a gene-density or histone
   track when `--extraTrack` is given (`hicPCA.py:139-193`), but with no
   `--extraTrack` the sign is whatever LAPACK produced.
3. It is also the reason `hicPCA` is slow and memory-hungry: `dgeev` on an
   n-by-n matrix computes all n eigenpairs and returns them complex, when two
   real ones were asked for. Section 4.5.

The equivalence rule for `hicPCA` compatibility mode is therefore: compare
`abs(eigenvector)` elementwise at E3; separately assert that the sign pattern is
either identical or globally flipped per chromosome; and when `--extraTrack` is
supplied, require identical signs, since the flip is then deterministic. Any
other outcome is a failure, not a tolerance. The corrected mode is E7 against
Python and is validated against the compatibility mode instead (section 5.8).

### 5.5 How tolerance is measured, exactly

For a matrix pair `(A_cpp, B_py)`:

1. **Sparsity pattern**: the sets of `(i,j)` with a stored entry must be equal
   for E2 and E3. For E4 the symmetric difference must be at most 0.1 % of
   `max(nnz_A, nnz_B)`, and the excess entries must all satisfy
   `abs(value) <= 1e-9`.
2. **Max relative error**: over the common support,
   `max_ij abs(a_ij - b_ij) / max(atol_floor, abs(b_ij))`, with `atol_floor = 1`
   for E3 and `1e-6` for E4. Using a floor rather than a pure relative error
   avoids the meaningless blow-up on values near zero.
3. **Correlation floor** (E4 only): Pearson correlation of the two value
   vectors over the common support, computed in float64 with the pairwise
   reducer.
4. **Integer and string fields**: always exactly equal, in every class except
   EN and E7.

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

### 5.7 Class EN: comparing against a nondeterministic oracle

`hicCorrectMatrix --correctionMethod KR` is the one tool whose Python reference
does not agree with itself (section 3.3). Fixing a tolerance by judgement would
be arbitrary; the tolerance is therefore **measured from the oracle**:

1. The harness runs the Python tool `N = 5` times on the same input, in the same
   environment, with everything else held fixed.
2. It computes `S`, the maximum pairwise relative difference among those five
   outputs, field by field, using the section 5.5 machinery.
3. It requires the sparsity pattern of the C++ output to equal the sparsity
   pattern shared by all five Python runs, exactly. (Measured: it is stable.)
4. It requires the C++ output to lie within `max(2 * S, 1e-12)` relative of the
   **median** Python run, elementwise.
5. It requires the C++ output to be **deterministic**: five C++ runs, and the
   `--threads 1` versus `--threads 16` runs, must be byte-identical. The C++
   port is held to a standard the reference does not meet, deliberately.
6. `S` is recorded in the report, per input, so that a change in the oracle's
   noise level is visible rather than silently absorbed.

Observed values of `S` for KR, for calibration: 1.5e-04 on matrix values and
7.5e-05 on correction factors for `Li_et_al_2015.h5`; the normalisation factor
alone spans 8.4e-03 on `gm12878_chr1.cool`. A fixed tolerance chosen before
these measurements would have been wrong by two orders of magnitude in one
direction or the other, which is the argument for measuring it.

If a future krbalancing release makes KR deterministic, EN collapses to E3 for
that tool and the class is retired. Nothing else in the corpus currently needs
it; the harness supports it generically so that it can be applied if another
tool turns out to be nondeterministic under the `N`-run check, which every tool
gets as part of criterion 3 in section 8.3.

### 5.8 Dual-mode: where a fix would change results

Rules 1 to 9 of section 4.4 are representation changes and alter no arithmetic,
so the memory workstream needs no mode switch. Two *accuracy and algorithm*
fixes do change results, and those are shipped as explicit modes rather than
being either forced on users or quietly dropped.

**The pattern.** A tool with such a fix gets a `--compatMode {v3,v4}` flag,
defaulting to `v3`. `v3` reproduces the Python reference including its defects
and is what the equivalence harness exercises. `v4` applies the fix. Both modes
are implemented over the same memory-efficient data structures, so choosing `v3`
costs nothing in memory. The difference between the two modes is **measured on
the real validation inputs and quantified in `STATUS.md`**, per tool, so that a
user switching modes knows the size of the change.

**Where it applies.**

| tool | `v3` (default) | `v4` | expected difference |
|---|---|---|---|
| `hicCorrectMatrix --correctionMethod KR` | float32 downcast of input values on ingestion (`krbalancing.cpp:27`); float32 accumulators in `rescale_norm_vector` (`:228-229`), summed in a **fixed** order (ascending outer index, then ascending inner index) so the result is deterministic while staying inside the oracle's noise envelope; `--perchr` reproduces the `.h5`-versus-`.cool` correction-factor discrepancy | float64 throughout; pairwise accumulation for the two rescale sums; `--perchr` rescales consistently regardless of output extension | to be measured; expected of the order of the oracle's own `S` (1e-4 to 1e-2 relative depending on matrix size), which is precisely why `v4` is worth having |
| `hicPCA` | `dgeev` on the covariance matrix, all eigenpairs, unsorted, columns taken by index | `dsyevr` with `range='I'` returning only the requested eigenvectors, sorted by descending eigenvalue, sign fixed by a deterministic convention (largest-magnitude component positive) unless `--extraTrack` decides it | eigenvector selection may differ when `dgeev`'s arbitrary order does not match eigenvalue order; to be measured and reported per chromosome on `mm9_reduced_chr1.cool` |

**Where it deliberately does not apply**, so that the flag does not proliferate:

- The pearson and covariance reduction-order change (4.4 rule 7) is covered by
  class E3. No flag.
- ICE marginals: the CSR row-order reduction reproduces `coo_matvec` exactly
  (section 5.3). No flag.
- `--enforce_integer` rounding: `np.rint` and `std::nearbyint` under
  `FE_TONEAREST` are both round-half-to-even. No flag.
- The `hiCMatrix` load-time field swap and the other quirks of section 2.7 are
  reproduced unconditionally in both modes. They are format semantics, not
  accuracy defects, and a `v4` that wrote `correction_factors` into the correct
  node would produce files the Python reads wrongly.
- krbalancing's `exit(0)` after 300 outer iterations
  (`krbalancing.cpp:115-119`) is **not** reproduced in either mode. Both modes
  raise an error and exit non-zero. Terminating the process with a success status
  and no output is not behaviour worth preserving, and no correct pipeline can
  depend on it. Recorded as a deviation in `STATUS.md`.

Every mode pair is validated three ways: `v3` against Python at the tool's
declared class; `v4` against `v3` with the difference quantified and recorded;
and both modes against themselves for determinism.

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
hicpro, 2D-text, npz), the pairwise reducer, Cephes, the BED/narrowPeak reader,
`reduce_matrix`, and the CLI/argparse compatibility layer. The memory rules of
section 4.4 are properties of this tier: if zero-copy ingestion, triangular
storage and streaming output are not in place here, no later tier can meet its
budget.

The **argparse compatibility layer** is not optional. Every tool's help text,
argument names, short options, `choices`, defaults, `nargs`, `metavar` and
error messages are part of the interface and are asserted by the Python tests
(and by any Galaxy wrapper). Build one `hicx::ArgParser` that reproduces
argparse's grouping (`Required arguments` / `Optional arguments`), its
`--help` rendering, its `%(prog)s {version}` version action, its prefix matching
(`test_hicPlotMatrix.py:444` passes `--log1` and relies on it resolving to
`--log1p`) and its error text. Doing this once in tier 0 costs a day and saves
46 hand-written parsers.

Exit criteria for tier 0: round-trip every matrix in `test_data/` through
C++ read and C++ write in every format the Python supports, and have the Python
loader produce identical in-memory state; pass the C++ unit tests for the
quirks in section 2.7; and demonstrate on `gm12878_chr1.cool` that loading the
matrix and writing it back costs no more than `1.3 * W + C`.

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
`hicMergeTADbins`, `hicAverageRegions` and `hicNormalize` which write matrices.

### Tier 3 - float matrix math (6 tools)

`hicTransform` (260, E2 for obs_exp/norm, E3 for pearson/covariance),
`hicCorrectMatrix` (779, E3 for ICE / EN for KR),
`hicPCA` (412, E3 plus the sign rule, dual-mode),
`hicCompartmentalization` (223, E3),
`hicInterIntraTAD` (514, E3),
`hicPlotSVL` (261, statistics at E3, plot deferred to tier 6).

Rationale: this is where the numeric policy and the memory budget are both
proven. `hicTransform` first because `hicPCA`, `hicDetectLoops` and
`hicFindTADs` all reuse its obs/exp and Pearson kernels, and because the pearson
path is the second-largest memory win in the plan; `hicCorrectMatrix` second
because it is the most-used tool in the suite and holds the largest memory win.

### Tier 4 - alignment and matrix construction (3 tools + the `.hic` reader)

`hicBuildMatrix` (270 + `lib/buildMatrixMethods.py` 900), `hicBuildMatrixMicroC`
(226), `hicQuickQC` (127), plus the `.hic` reader that completes
`hicConvertFormat` (section 3.6).

Rationale: htslib enters here and nothing else depends on it. `hicBuildMatrix`
is the largest single unit of work in the port (BAM pair iteration, supplementary
alignment resolution, dangling-end and self-circle classification, restriction
fragment binning, the QC table). Its output matrix is integer, so **E1/E2**, and
its QC tables are **E0**. `hicQuickQC` shares the same read classification code
and is essentially a free follow-on. This tier is also where the
threads-share-not-fork rule pays: the Python spawns workers that each hold a
share of the pixel buffers.

### Tier 5 - TAD, loop and differential calling (5 tools)

`hicFindTADs` (1,368), `hicDetectLoops` (1,093), `hicDifferentialTAD` (518),
`hicMergeDomains` (428), `hicAggregateContacts` (976, data path only; its plot
follows the tier 7 rule).

Rationale: the heaviest algorithms, all depending on tier 3 kernels. Classes E5
for the call sets and E3 for the continuous intermediates (section 5.3).
`hicFindTADs` additionally writes a z-score matrix and a `.bm` bedgraph matrix,
both **E3**, and those are the real regression detectors; the domain BED is E5.
All three of these tools take the rule 2 exemption and materialise both
triangles, so their budgets are the loosest in the plan and should be revisited
once they are working.

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

**Done, merged 2026-09-15 in `2452a010`.** All seven chic tools are ported.
chicDifferentialTest reproduces scipy 1.14.1's statistics bit for bit.

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
  images stay bit-comparable to the current masters, and the memory-heavy step
  (extracting and transforming the region) moves into the budgeted C++ core.
- **(b) C++ SVG/PNG backend.** Feasible for the simple line plots
  (`hicPlotViewpoint`, `hicPlotAverageRegions`, `hicPlotDistVsCounts`), not for
  the heatmaps with colorbars and genomic axes, and not at all for
  `hicPlotTADs`.
- **(c) Defer.** Ship v4 without the tool.
- **(d) Drop.**

Recommendation, per tool:

| tool | recommendation | reason |
|---|---|---|
| `hicPlotMatrix` | (a) shell | heatmap + colorbar + optional bigWig track + `--perChromosome` layout; the compute (region extraction, log transform, obs/exp) moves to C++, the draw stays matplotlib. Six of its tests are gated behind a 120 GB memory `skipif`, which the C++ compute path should make unnecessary |
| `hicPlotTADs` | (a) shell, unchanged | already a 9-line delegation to `pygenometracks.plotTracks`; there is nothing to port |
| `hicPlotViewpoint` | (a) shell, with (b) as a later option | simple line plot; C++ emits the data, matplotlib draws |
| `hicPlotAverageRegions` | (a) shell | small imshow; reads the `.npz` the C++ `hicAverageRegions` writes |
| `hicPlotDistVsCounts` | (a) shell | the interesting part (distance-vs-count reduction, per-chromosome fits) is compute and moves to C++ |
| `hicCorrelate` | (a) shell | the correlation matrix and the complete-linkage clustering are C++ at E3; the heatmap and scatter stay matplotlib |
| `hicPrepareQCreport` / `hicQC` | (a) shell | pandas table aggregation plus five bar charts and a Jinja2 HTML template; port the table aggregation to C++ at E0 and keep the rendering in Python |
| `chicPlotViewpoint` | (a) shell | as `hicPlotViewpoint` |

**Done, merged 2026-09-15 in `b71510ba`.**
- All eight tools, plus the figures formerly refused under contract rule 7.
- E0 against the Python-drawn figures, with one E6 at RMS 1.17.
- matplotlib 3.8.4 is pinned and checked at run time.

So: **(a) for all eight.** The v4 deliverable is a C++ core plus a thin,
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
  matrix rows and obs/exp windows, which is exactly what the core computes, and
  moving them is also the memory win for this tool. Class **E7**.
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
| 9 | beyond the Python: new features and a differential redesign | 14 items (9.1 to 9.14) |
| 10 | Python GUI (PySide6) with workflows and visualisation, Linux and macOS | 8 items (10.1 to 10.8) |
| 11 | sparse Lanczos eigensolver for hicPCA, C++-only option | 1 item |
| 12 | hicmatrix 18: hicmatrix with a C++ core and its unchanged Python API, then v4 on it | 2 phases |

Tiers 2, 6 and 7 are independent of tiers 3-5 after tier 0 lands, so up to three
workers can proceed in parallel from that point.

### Tier 9 - beyond the Python (added 2026-09-13)

Features HiCExplorer 3.7 does not have, requested by the project owner. None has a
Python oracle, so each is validated under class **EX** with the criterion stated
here. Where part of a feature does have a Python path (for example `.hic` to
cool through hic2cool), that part keeps its ordinary class. Memory, CPU-time and
determinism gates apply as for every port, measured against the external
reference implementation where no Python exists.

**9.1 Native `.hic` reading and writing.** Started first, at the owner's request.
- *Where:* its own library, `~/src/hicfilecpp`, a C++ implementation of the Juicer
  `.hic` format with no Python or Java dependency, consumed by v4 the way
  coolercpp is. HiCExplorer behaviour stays in v4.
- *Scope:* read versions 6 to 9 (header, chromosomes, attributes, BP and FRAG
  resolutions, blocks, normalisation vectors NONE, VC, VC_SQRT, KR, SCALE,
  expected-value vectors, observed and observed/expected). Write versions 8 and
  9 from pixels at one or many resolutions, with expected values and the
  normalisation vectors Juicer's `pre` and `addNorm` compute. Versions 6 and 7
  (older deposits) are read, including version 6 block records. Writing them is
  refused with the reason, because no obtainable Juicer tools release writes
  them to validate against: `Preprocessor.java` has written version 8 since its
  first Juicebox commit (June 2015), and the oldest downloadable jars (1.6.2,
  1.7.5, 1.7.6) write version 8. Versions below 6 and above 9 are refused
  explicitly, as straw does. (Extended 2026-09-14: the project owner asked for
  both the old and the new format; 8 and 9 landed in `8a2fa526`, reading 6 and 7
  on `v4-hic-legacy`.)
- *v4:* `hicConvertFormat` reads `.hic` into cool, mcool, h5 and the text formats,
  and writes `.hic` from h5, cool and mcool. Afterwards every matrix-reading
  tool accepts a `.hic` with a resolution and normalisation selector.
- *Validation:* `.hic` to cool against the Python `hic2cool` path, **E1/E2**.
  Reading against `hicstraw` record by record at every resolution and
  normalisation: pixels exact, vectors **ED**. Writing: files read back through
  `hicstraw`, `hic2cool` and Juicer tools `dump` must give the source pixels
  exactly, and expected and normalisation vectors **ED** against Juicer tools
  `pre` and `addNorm` on the same input. Output bytes need not match Juicer's.
- *Data:* `SRR1791297_30.hic` (in the repository) for cases; the 423 MB
  `GSM6505198` and 5.3 GB `GSE63525_HMEC` files (read-only mounts) for memory
  and time at scale. Those files are version 8. Version 7: the 40 GB
  `GSE63525_GM12878_insitu_primary+replicate_combined_30.hic` (three copies on
  disk) and `GSE63525_IMR90_combined_30.hic`. There is no real version 6 file, so
  version 6 is validated on real GM12878 pixels re-encoded into version 6 block
  records.

**9.2 `.pairs` input for hicBuildMatrix.** 4DN and pairtools `.pairs`, plain or
bgzipped. *Validation:* **E2** against `cooler cload pairs` on the same file, and
**E2** against the BAM route when the pairs file carries exactly the pairs
hicBuildMatrix keeps from that BAM.
*Real data (2026-09-15):* GSE234292's `pairs_sample.zip` holds lists of pixel
paths, not pairs, and GEO has only mcool files for that series. The scale
check uses ENCODE ENCFF849QYT (mm10 in situ Hi-C, the first 135,889,417
sorted upper-triangle pairs).

**9.3 Stripe detection (`hicDetectStripes`).** *Validation:* planted stripes in
real matrices recovered at a recall and precision fixed in this plan before the
method is chosen, and agreement with Stripenn on GM12878 reported as a Jaccard
index with the disagreeing calls examined.

**9.4 SCALE normalisation in hicCorrectMatrix.** *Validation:* **ED** against the
SCALE vector Juicer tools `addNorm` computes on the same matrix.

**9.5 Multi-resolution files (coarsen, zoomify).** Delivered with coolercpp
milestone 3. *Validation:* **E1/E2** against `cooler coarsen` and `cooler
zoomify`.

**9.6 Multiscale loop calling.** Loops called across resolutions and merged.
*Validation:* planted loops as in 9.3, and agreement with HiCCUPS and Mustache on
GM12878, with the criterion fixed here before implementation.

**9.7 Differential analysis that controls false positives.** The project owner
reports too many false positives in differential TADs, loops and A/B
compartments. This was measured on 2026-09-13 on GSE234292 (mouse, wild type and
BPTF knockout, two biological replicates each, `~/data/Hi-cGAN/biological_replicates/GSE234292`).
The runs used the C++ ports, which match the Python on every harness case,
with ICE on chr1-19 and X (`--filterThreshold -1.5 5`). A replicate against
its replicate is a null comparison. Scripts and outputs are in the session
scratchpad, under `fp_calibration.sh` and `fpcal/`.

| workflow | null result (replicate vs replicate) |
|---|---|
| hicDifferentialTAD, defaults (`-m all -mr one`, p 0.05), 50 kb, 2,730 TADs | 561 TADs (20.5 %) wt, 589 (21.6 %) knockout; p <= 0.05 in 11 % of tests instead of 5 % |
| the same with Benjamini-Hochberg on each test | still 290 TADs (10.6 %) |
| the same after masking the union of both samples' filtered bins in both matrices | 25 TADs (0.9 %); p <= 0.05 in 0.3 % of tests |
| the C++ option `--sharedMask` (`74ec19cc`), without and with `--correctForMultipleTesting fdr` | 25 TADs (0.9 %), the same 25; with FDR 0. Knockout replicates: 24 (0.9 %), with FDR 0 |
| hicDetectLoops, defaults, 10 kb, uncorrected | 79 % of rep1 loops and 83 % of rep2 loops not found in the other replicate (1-bin tolerance) |
| hicDetectLoops on the ICE matrices (defaults, and `-pit 1`) | 0 loops |
| hicPCA E1, 100 kb, orientation aligned per chromosome | 1,148 of 24,468 bins (4.7 %) switch A/B; 2.3 % outside the weakest quartile of abs(E1) |

Findings:
- The dominant false-positive source in hicDifferentialTAD is **asymmetric bin filtering**. ICE filters bins per sample, and a bin removed in one matrix is a zero row against real values in the other. TADs within 500 kb of such a bin were called 63.4 % of the time (n = 816), all others 2.3 % (n = 1,914). TAD size, coverage and local depth ratio barely matter.
- Wild type against knockout calls 59.0 % of TADs. With `--sharedMask` it calls 1,508 (55.2 %), and with FDR as well 955 (35.0 %), with one replicate per side. That is not an interpretable result without replicates and an effect size. An earlier figure of 51.4 % came from masking through `hicAdjustMatrix --regions`, which zeroed 4,970 rows for a BED of 4,067 bins (STATUS F60).
- Loops: comparing call lists is dominated by threshold instability, since four fifths of calls do not replicate.
- Compartments: a naive sign switch calls about 115 Mb between replicates of one condition.
- No differential tool applies a multiple-testing correction. That includes chicDifferentialTest and chicSignificantInteractions; only hicFindTADs has one. hicDetectLoops' docstring promises FDR correction, which the code does not apply.

Work:
1. **hicDifferentialTAD, dual mode (section 5.8).** Done, merged 2026-09-14 (`74ec19cc`).
   - Default output stays Python-equivalent.
   - `--sharedMask` masks the union of both matrices' invalid bins in both.
   - `--correctForMultipleTesting {none,fdr,bonferroni}` adjusts each test across TADs.
   - Both options are E0 against a Python reference that masks through hicmatrix and adjusts the Python tool's p-values.
2. **One count-based differential engine** for TADs, loops and compartments.
   - Replicates per condition, with a negative binomial model: dispersion estimated from replicates, offsets for library size and distance decay, a shared bin mask, and Benjamini-Hochberg FDR with a minimum fold change.
   - The tested units: per TAD, aggregated counts per distance stratum plus boundary insulation; per loop, the union of loop positions from all samples, tested against local background; per compartment bin, a GC-oriented compartment score.
   - With a single sample per condition the tool refuses, or runs only with an explicit option that labels its output exploratory.
3. **Multiple-testing correction** for the chic tools and hicDetectLoops, dual mode.
   - Done for chicSignificantInteractions and chicDifferentialTest (`2452a010`):
     `--correctForMultipleTesting`, with `none` equal to the Python.
   - hicDetectLoops is still open.

Gate (EX, fixed before implementation):
- **Null comparisons** at FDR 0.05 call at most 1 % of tested TADs, loops or compartment bins. The nulls: wt rep1 vs rep2, knockout rep1 vs rep2, and the label swap {wt rep1, knockout rep2} vs {wt rep2, knockout rep1}. In each null the fraction of p <= 0.05 lies between 3 % and 7 %.
- **The label-swap null** is evaluated with the genotype as block (decided
  2026-09-15). An unpaired swap cannot meet the 3 % p-value bound by
  construction, because its within-group variance contains the wt against
  knockout effect. The unpaired swap is still reported, without a gate.
- **Planted differences** in real wt matrices are recovered with recall >= 0.8 at 2-fold and an observed FDR <= 0.10. The plants: contacts scaled inside chosen TADs, at chosen loop pixels, and across chosen compartment bins, at 1.5-fold and 2-fold.

**9.8 Structural variant and translocation detection from Hi-C.** Inter- and
intra-chromosomal breakpoints and balanced translocations, relevant for
leukaemia samples. *Validation:* the known rearrangements of K562 (including
BCR-ABL1) recovered from `GSE63525_K562_combined_30.hic`, plus agreement with
hic_breakfinder on the same data.

**9.9 Copy-number-aware normalisation.** Copy number removed before or
alongside balancing, so gains in aneuploid samples do not appear as
differential contacts. *Validation:* ED against OneD (or HiCnv) on K562, and
the 9.7 null gate on a sample pair with differing copy number.

**9.10 Reproducibility score (HiCRep SCC).** *Validation:* ED against the
reference hicrep implementation on the GSE234292 replicate and condition pairs.

**9.11 Count subsampling to equal depth.** Binomial thinning with a seed,
replacing rescaling for comparisons. *Validation:* exact total, per-pixel
expectation checked statistically, byte-identical output for a given seed at
any thread count.

**9.12 Allele-specific matrices** from phased BAMs (haplotype tags). *Validation:*
E2 that the haplotype matrices plus the unassigned matrix sum to the unsplit
matrix, and agreement with an established allele-specific pipeline on the same
reads.

**9.13 Effective resolution estimate** (Rao et al. 2014 criterion: the finest bin
size at which 80 % of bins have at least 1,000 contacts). *Validation:* exact
against a direct computation of the definition.

**9.14 Micro-C patterns: fountains and jets.** *Validation:* planted patterns as
in 9.3, and agreement with a published fountain caller, such as fontanka.

Order: 9.1 (with reading of versions 6 and 7) and 9.7 step 1 are done (2026-09-14); then
9.2, 9.3, 9.4, the rest of 9.7, 9.6 and 9.8 to 9.14, interleaved with the
remaining tier 6 tools. 9.5 follows coolercpp milestone 3.

### Tier 10 - HiCExplorer GUI (added 2026-09-13)

The project owner asked for a Python GUI that runs on Linux (including Ubuntu)
and macOS, runs end-to-end workflows and offers visualisations. It is built on
tier 7 option (a): C++ computes, Python draws. The GUI request therefore settles
tier 7 in favour of (a) for every plotting tool.

**Stack.**
- **Toolkit:** Qt 6 through PySide6 (LGPL). Native on Linux and macOS, arm64
  and x86_64; Windows works with the same code if wanted.
- **Interactive views:** Qt with an OpenGL-backed image widget (pyqtgraph or
  equivalent).
- **Saved figures:** go through the tier 7 matplotlib shells, so a figure
  exported from the GUI equals the CLI tool's figure.
- **Tool runs:** the GUI runs the C++ tools as subprocesses. A tool crash cannot
  take the GUI down, and every run is reproducible from its logged command line.
- **Interactive data access:** Python bindings (pybind11) over libhicx4,
  coolercpp and hicfilecpp. They query only the region and resolution on
  screen, so browsing stays within the memory budget; a matrix is never loaded
  whole.

**Components** (top-level `gui/`, Python package `hicexplorer_gui`):
1. **Tool specifications.** Each C++ tool emits a machine-readable description
   of its arguments (`--help-json`), which is the single source for GUI forms,
   validation and file pickers. This is the argparse layer listed as open work
   in `STATUS.md`.
2. **Workflow engine, usable without the GUI** (`hicexplorer-workflow run
   file.yaml`).
   - A workflow is a DAG of tool steps in a versioned YAML file.
   - It resumes by hashing inputs and parameters, and records per-step logs,
     peak RSS and CPU time, with cancellation and a thread budget.
   - It exports to a shell script and a Snakemake file; a Galaxy workflow export
     can come later.
3. **Visualisation.**
   - *Matrix browser:* cool, mcool, h5 and `.hic`. Raw, log, obs/exp and
     Pearson views, and resolution switching on zoom for mcool and `.hic`.
     Side-by-side and difference views, and coordinate navigation.
   - *Tracks under the matrix:* TADs, loops, BED, bigWig, bedGraph and
     eigenvectors.
   - *Analysis views:* QC report, distance-decay curves, viewpoint (virtual 4C)
     plots, aggregate contacts and saddle plots, and correlation or HiCRep
     heatmaps. Differential results appear as tables and volcano plots linked
     to the matrix browser.
4. **Projects and history:** inputs, runs, parameters and provenance per
   project.

**Workflow templates.**
- **Hi-C:** BAM or `.pairs`, then hicBuildMatrix with QC, merging replicates,
  binning or zoomify, correction, TADs, loops, compartments and a report.
- **Differential:** two conditions with replicates, a shared bin mask, then the
  9.7 engine for TADs, loops and compartments.
- **Capture Hi-C:** the chic tools, from background model to differential test
  and viewpoint plots.
- **Conversion and QC:** format conversions, hicInfo, correlation and HiCRep,
  distance decay.

A step whose tool is not yet ported shows as unavailable with the reason (rule
4); nothing is silently substituted.

**Validation.**
- *Workflow runner:* every step's outputs are E0 against invoking the same tool
  with the logged command line. The Hi-C template runs end to end on real
  BAMs, and the differential template on GSE234292.
- *Exported figures:* E6 against the tier 7 CLI plotting with the same
  parameters.
- *Interactive views:* the displayed values for a region are E2 against the
  tools' region extraction. Peak RSS stays within budget while browsing
  `gm12878_chr1.cool` and the 5.3 GB `.hic`.
- *GUI tests:* pytest-qt with `QT_QPA_PLATFORM=offscreen`, plus screenshot
  review on a real display. Layouts must stay usable from 1280x720 to 4K,
  with no clipped controls.
- *macOS:* the C++ core must build with Apple clang on arm64 and pass the
  harness there. Its x86 SIMD dispatch needs NEON or scalar paths with
  byte-identical output (`OPTIMIZATION.md`), and the non-AVX-512 fallback is
  currently untested (`STATUS.md`). **This cannot be verified on the build
  machine; it needs a Mac or a CI runner (open question for the owner).**

**Packaging:** a conda recipe for the C++ tools and the GUI, a pip wheel for the
GUI, and a macOS app bundle later.

**Order:**
- 10.1 tool specifications;
- 10.2 Python bindings;
- 10.3 headless workflow engine;
- 10.4 GUI shell with tool forms and the run view;
- 10.5 matrix browser;
- 10.6 tier 7 shells and the analysis views;
- 10.7 workflow templates;
- 10.8 macOS.

10.1 to 10.3 are done, merged 2026-09-14 in `f4c15dc9`.
- 10.1: 30 tools, spec equal to argparse.
- 10.2: bindings E2 against cooler, hicmatrix and hicstraw, 32 MB for a 2 Mb
  window on the 40 GB file.
- 10.3: workflow engine with an end-to-end E0 run.

10.4 and 10.5 are done, merged 2026-09-14 in `0b3b56db`: the PySide6 shell with
forms for all 30 tools and the matrix browser, reviewed on screenshots at three
sizes.
- 10.6: done. The tier 7 shells landed in `b71510ba`; the GUI analysis views
  and their figure export landed in `e9c5bbaa`.
- 10.7: done, `e9c5bbaa`. Four templates run end to end with outputs equal to
  their logged commands.
- 10.8 needs a Mac or a CI runner.

### Tier 11 - sparse eigensolver for hicPCA (added 2026-09-14)

Requested by the project owner, to follow the GUI (tier 10).

**Problem.** hicPCA's port holds dense per-chromosome matrices, because the
Python chooses eigenvectors by their position in LAPACK `dgeev`'s unsorted
output and only a covariance bit-identical to `np.cov` reproduces that choice
(section 5.4). Memory therefore grows with the square of the largest
chromosome's bin count:
- 269 MB of the 344 MB peak on `small_test_matrix` (5,801 bins);
- about 0.8 GB per dense matrix for human chr1 at 25 kb, and about 5 GB at 10 kb.

**Method.** A C++-only option `--eigenSolver {dense,lanczos}`. The default,
`dense`, stays Python-equivalent.
- **`lanczos`** computes only the requested eigenvectors with a Krylov solver
  (Spectra, header-only on Eigen, MPL-2.0) and never builds the dense matrix.
- **Implicit centering.** Every product with the Pearson matrix, and with its
  covariance as hicPCA takes it, is evaluated from the sparse obs/exp matrix,
  a rank-one centering term and a diagonal scaling.
- **Memory** is then the sparse input plus a few vectors per chromosome.
- **Eigenvectors** are returned ordered by eigenvalue. Where LAPACK's order is
  not by eigenvalue, `lanczos` picks a different first eigenvector than the
  Python. That is a documented deviation, reported per case, never presented as
  equal.
- **Sign** follows the same gene or histone track rule as the default.

**Validation** (class ED against `dense`, fixed before implementation):
- On every hicPCA case whose eigenvalue order is unambiguous, each requested
  eigenvector agrees with `dense` to three significant digits after sign
  alignment. Cases where the order differs are listed individually with both
  spectra.
- **Real large data:** GSE234292 mouse chromosomes at 25 kb and 10 kb, and
  human GM12878 chr1 at 25 kb. Report peak RSS and CPU time for the Python,
  C++ `dense` and C++ `lanczos`, plus the agreement of `lanczos` with `dense`
  wherever `dense` fits in memory.
- **Gates:** `lanczos` peak RSS at most the sparse input plus a stated
  per-chromosome vector bound. Determinism at `-t1` and `-tN`, with a fixed
  starting vector.

### Tier 12 - hicmatrix 18: hicmatrix with a C++ core and its Python API (added 2026-09-15)

The project owner asked that hicmatrix become its own independent C++ library.
It is not a separately named library but **a new version of hicmatrix itself**,
which keeps offering the Python API. v4's core already reimplements much of
hicmatrix internally: `cool_adapter.cpp` (what hicmatrix adds on top of
cooler), `h5_file.cpp`, `text_formats.cpp`, `tool_matrix.cpp`, `hic_matrix.cpp`,
and parts of `hdf5_util.cpp` and `matrix_ops.cpp`, used by 18 tools. This tier
moves that into hicmatrix 18 and makes v4 consume it.

**Repository.** The upstream clone `~/src/HiCMatrix` (deeptools/HiCMatrix,
GPL-3, so no separate licence decision), on a local branch created from the
`17.2` tag (the oracle version), never pushed. Version 18.0.
- The C++ core has a C++ API that v4 links against.
- The `hicmatrix` Python package keeps 17.2's public API through bindings:
  `hicmatrix.HiCMatrix.hiCMatrix` with its public attributes (`matrix` as a
  `scipy.sparse.csr_matrix`, `cut_intervals`, `nan_bins`,
  `correction_factors`, `distance_counts`, `bin_size` and so on),
  `hicmatrix.lib.MatrixFileHandler`, the format classes and
  `hicmatrix.utilities`. HiCExplorer's 56 importing Python files must work
  unchanged.
- Heavy work (file I/O, obs/exp, z-score, masking, reordering) runs in C++ on
  the CSR buffers, copying only where the Python API requires an owned object.
- Files it writes name HiCMatrix 18 as their producer (section 5.0.1).

The public API to keep, from hicmatrix 17.2:

**Target:** hicmatrix 17.2 (the reference installed in the oracle
environment), about 2,000 lines.
- **`hiCMatrix`:** every public method, including:
  - loading with region cuts, `restoreMaskedBins`, `pUpperTriangleOnly` and
    the matrix format;
  - `save`, `fillLowerTriangle`, `setMatrix`, `getBinSize`,
    `getChrBinRange`, `getRegionBinRange`, `getDistList`;
  - `convert_to_obs_exp_matrix`, `convert_to_zscore_matrix`;
  - `keepOnlyTheseChr`, `maskBins`, `maskChromosomes`, `restoreMaskedBins`,
    `reorderChromosomes`, `reorderBins`;
  - `setCorrectionFactors`, `truncTrans`, `filterOutInterChrCounts`,
    `get_chromosome_sizes`.
- **`lib/`:** `MatrixFileHandler`, and the formats h5 (PyTables layout), cool
  and mcool (correction-factor handling, the hic2cool and hicmatrix version
  logic, metadata), homer, hicpro, ginteractions and scool.
- **`utilities`.**

**Layout:** inside the HiCMatrix repository.
- A C++ core with a CMake package for C++ consumers.
- A Python build (pybind11, a wheel or conda recipe) that installs the
  `hicmatrix` package.
- Harness, `docs/API_MAPPING.md`, `DEVIATIONS.md` and `PROVENANCE.md`.
- Depends on HDF5 and coolercpp, at a pinned commit. PyTables, pandas and
  intervaltree stay only where the Python API hands out their objects.
- Every hicmatrix quirk v4 reproduces today (for example the correction-factor
  swap, NaN bins rebuilt from empty rows, the part split and sum) belongs in
  hicmatrix 18, as the Python behaviour.

**Phase 1, hicmatrix 18,** independent of v4.
- **hicmatrix's own test suite** (`hicmatrix/test`, with its 32 MB of real
  test data) passes unchanged against version 18.
- **HiCExplorer's Python test suite** (the 3.7.7-dev tree) runs with
  hicmatrix 18 installed in place of 17.2. It must give the same pass, fail,
  xfail and xpass results as the recorded baseline, in a separate venv, never
  the oracle venv.
- **Harness against hicmatrix 17.2,** method by method and format by format,
  on real matrices (the h5, cool, mcool, scool, homer, hicpro and ginteractions
  test data, plus gm12878_chr1.cool). Classes are E0 for written files (after
  normalising the producer name), E2 for loaded arrays and intervals, and ED
  only where 17.2 itself computes in floating point (obs/exp, z-score).
- **Mutation checks,** as in coolercpp.
- **Memory and CPU** against 17.2 for load, save and the transforms.
- **Build:** clean-export build of the C++ core and the Python package, ctest
  including a consumer test that reads the project version, and determinism.

**Phase 2, v4 on hicmatrix 18,** after the open branches merge.
- v4's core replaces its internal implementations with hicmatrix 18's C++ API.
  Only v4-specific code stays: streaming paths, memory-budgeted readers and
  tool glue.
- **Pinning** is like `HicxCoolercpp.cmake`.
- **Gates:** a merge-style full regression (contract rule 13), determinism,
  memory and CPU no worse than before, and the gui and bindings suites. The templates grow as tier 6 and
tier 9 tools land.

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
  gates. Those memory gates exist because of the blowups of section 4.3, so the
  memory workstream also buys back test coverage.
- **Roughly 170 items can actually fail on a regression.**

Assertion defects found, each of which makes a nominally-covered tool
effectively uncovered:

| location | defect |
|---|---|
| `test_compute_function.py:7` | `pTries = 1` overwrites the caller's retry count, so every `compute(main, args, 5)` runs once. Harmless, but the retry the suite thinks it has does not exist |
| `test_hicMergeDomains.py:70,84-85,106-107` | every `are_files_equal(...)` call is missing `assert`; the result is discarded. `hicMergeDomains` is "did not crash" only |
| `test_hicInterIntraTAD.py:55` | same, `are_files_equal` without `assert` |
| `test_hicHyperoptDetectLoopsHiCCUPS.py:29,64` | `are_files_equal` is `return True`, and it is called without `assert` |
| `test_hicCorrectMatrix.py:84` | the KR/cool check is a range test `3e9 < sum//2 < 3688003604`; the elementwise comparison is commented out at :85-86. Given the measured nondeterminism of KR (section 3.3) a range test may in fact be the only thing that could have passed reliably |
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
a distinct code path in both ICE and KR and which carries the defect of
section 2.7 quirk 8).

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
3. **A characterization test for a nondeterministic tool must record the noise,
   not a value.** For KR specifically, the characterization test runs the Python
   five times and records the envelope `S` (section 5.7), not a single output.
   Pinning one run's output would produce a test that fails against its own
   reference implementation.
4. **The existing masters are the reference, but they are not trusted until
   regenerated.** Because every image assertion is xfail-ed, the PNG masters may
   not match what today's matplotlib produces. Before tier 7 starts, run every
   image test with the xfail removed and record which masters are stale;
   regenerate them from the current Python and commit them with a note, or the
   C++ port will be validated against images no one has verified in years.
5. **`number_of_tests.txt` will go up.** Adding characterization tests raises the
   collected count; `test_pytest_collected_items.py` is a `>=` ratchet, so this
   is benign, but the file will show up in every diff. Consider pinning it once
   at the end of the characterization work rather than letting each run rewrite
   it.

A full coverage run (`pytest hicexplorer/test --cov=hicexplorer
--cov-report=json`) was started for this plan and had completed only 99 of the
~588 items after 50 minutes, so the per-line coverage numbers are not in this
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
| `Li_et_al_2015.h5` | 14,215,139 B | 11,104 bins, 1,843 bp, chrX only, 1,661,678 stored (upper-triangle) nonzeros and 3,313,107 after symmetrization, 855 NaN bins, float64 data, sum 17548966.536917936 | the h5 reader/writer, ICE, KR, obs/exp, pearson, the pairwise reduction, NaN-bin handling, the KR noise envelope |
| `Li_et_al_2015.cool` | 13,001,960 B | the same matrix in cool; `hicInfo` prints 1,661,678 here and 3,313,107 for the h5, see quirk 7 in section 2.7 | the cool reader/writer and the h5-cool round trip |
| `Li_et_al_2015_twice.h5` | 14,036,532 B | | `hicSumMatrices` |
| `Li_cut.h5` | 365,040 B | | fast smoke variant of the above |
| `small_test_matrix.h5` / `.cool` | 289,027 / 172,846 B | 33,754 bins, 35,857 nnz, 15 chromosomes | multi-chromosome ordering, `chrBinBoundaries`, `keepOnlyTheseChr`, blosc chunk layout reference |
| `small_test_matrix_50kb_res.h5` / `.cool` | 111,138 / 105,170 B | | `hicMergeMatrixBins`, `hicNormalize` |
| `matrix.mcool` | 2,444,203 B | 5 resolution groups named `/0`../`/4` (the legacy layout, **not** `/resolutions/<res>`), format-version 2 | the mcool reader; a reader that only understands `/resolutions/` will fail here |
| `hicBuildMatrix/multi_small_test_matrix.mcool` | 441,383 B | `/resolutions/{5000,10000,20000}` | the other mcool layout |
| `hicTADClassifier/gm12878_chr1.cool` | 79,370,717 B | 24,926 bins, 10 kb, chr1 only, 61,804,782 stored nonzeros (123,587,194 symmetric), `W = 741.9 MB` | **the memory benchmark**: KR, ICE, threading determinism, the budget gate |
| `hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1_chr2.cool` | 21,983,415 B | 100 kb | `hicDifferentialTAD`, `hicInterIntraTAD` |
| `hicDifferentialTAD/GSM2644947_Auxin2days-R1.100000_chr1_chr2.cool` | 28,429,077 B | 3,790 bins, 100 kb, chr1+chr2, 4,208,340 nnz | as above; the two-chromosome case |
| `hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1.h5` | 9,654,418 B | | h5 form of the same |
| `hicDetectLoops/GSE63525_GM12878_insitu_primary_2_5mb.cool` | 1,602,554 B | | `hicDetectLoops` |
| `hicCorrectMatrix/gm12878_raw_values.cool` | 1,480,300 B | int32 counts | ICE and KR on raw integer counts, where the float32 downcast of section 3.3 is inert |
| `hicCorrectMatrix/gm12878_KR.cool` | 2,130,904 B | 1,254 bins, has a `/bins/weight` column | the divisive-correction load path (`correctionFactorTable`) |
| `hicPCA/mm9_reduced_chr1.cool` | 722,538 B | 9,760 bins, 20 kb, mm9 chr1, 470,730 nnz, `W = 5.6 MB`, dense block `D = 762 MB` | `hicPCA` both modes, with `pca1.bedgraph`/`pca2.bedgraph`/`pca1.bw` as masters; the worst memory-to-data ratio in the corpus at 727x |
| `hicPCA/obsexp_norm.h5` | 56,085 B | 4,526 bins, 2,207 nnz | `hicCompartmentalization`, and the obs/exp reference for `hicTransform` |
| `hicAdjustMatrix/gm12878_1_2_3.cool` | | 693 bins, 1 Mb, chr1-3, 223,446 nnz | the multi-chromosome `--interIntraHandling` paths |
| `hicValidateLocations/GSM1436265_RAD21_ENCFF002EMQ_10kb.cool` | | 313,762 bins @10 kb, **93 contigs**, only 7,987 nnz | the many-contig, extremely sparse case: this is the file that will break a `BinIndex` that assumes few chromosomes or contiguous coverage, and the case where a per-bin `std::string` would dominate the footprint |
| `cHi-C/FL-E13-5_chr1.cool` | 420,484 B | 197,196 bins, 1 kb, chr1, only 83,665 nnz | the very-sparse, very-many-bins case; the whole cHi-C tier |
| `cHi-C/MB-E10-5_chr1.cool` | 461,587 B | | the second sample for every differential cHi-C step |
| `R1_1000.bam` / `R2_1000.bam` | 48,076 / 43,831 B | 1,000 read pairs | `hicBuildMatrix` fast path |
| `small_test_R1_unsorted.bam` / `small_test_R2_unsorted.bam` | 6,135,374 / 6,139,816 B | | `hicBuildMatrix` full path, QC tables, restriction-fragment mode, per-thread memory |
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
  dumps the same; the two must be identical (class E2 applied to each of
  `matrix.data`, `matrix.indices`, `matrix.indptr`, `cut_intervals`,
  `nan_bins`, `correction_factors`, `distance_counts`). Then C++ writes M' in F
  and Python reads M' and dumps again: identical. Then the file comparator for F
  runs against the Python-written original at that format's class. The
  `gm12878_chr1.cool` load-and-write round trip additionally has to meet
  `1.3 * W + C` (tier 0 exit criteria, section 6).
- **Tiers 1-6**: run the Python tool and the C++ tool with identical arguments
  on the designated inputs, compare every produced file with the comparator
  chosen by the file's extension at the tool's declared class, and record peak
  RSS for both. All CLI option combinations exercised by the Python test suite
  must be run, plus the combinations identified in `STATUS.md` as untested
  (those need a characterization test written first, contract rule 1).
- **Tier 7**: as above, and the image against the checked-in master with the
  same tolerance the Python test already uses.
- **Tier 8**: as declared per tool in tier 8 of section 6.
- **Dual-mode tools**: three runs, per section 5.8.

### 8.3 When a run counts as passed

A tool moves to `equivalence: pass` in `STATUS.md` only when **all** of:

1. every declared invocation produced the declared files, with exit code 0 where
   Python exits 0 and a non-zero exit where Python exits non-zero;
2. every comparator returned pass at the declared class;
3. the C++ tool is **deterministic**: five repeat runs are byte-identical, and
   for a tool that takes `--threads`, the `--threads 1` and `--threads 16`
   outputs are byte-identical too;
4. the tool's **peak RSS is within its budget** (section 4.5) on the designated
   large input, measured as described in section 10. This is a hard gate, not a
   report line: exceeding the budget fails the tool exactly as a comparator
   mismatch does;
5. the run is recorded in the harness report with the git commit of the C++
   tree, the input file checksums, the comparator output, and the measured peak
   RSS and budget.

Any deviation is recorded in `STATUS.md` with its reason; a tool with a recorded
deviation is `equivalence: deviation`, never `pass`.

Criterion 4 replaces the softer "no more than 1.5x the Python peak RSS" of an
earlier revision. A ratio against the Python is the wrong gate when the Python
number is itself twelve times the working set: it would let the port ship at
13.8 GB for KR and call it a pass.

## 9. The equivalence harness (`cpp/scripts/`, to be implemented by the implementing agent)

### 9.1 CLI

A single Python 3 entry point, run with the contract's venv interpreter, with no
dependencies beyond what that venv already has (numpy, h5py, cooler, PIL via
matplotlib):

```
cpp/scripts/equiv.py run     [--tool NAME]... [--tier N]... [--case ID]...
                             [--cpp-bin DIR] [--py-python PATH] [--jobs N]
                             [--out DIR] [--keep-workdirs] [--update-baseline]
                             [--noise-runs N] [--skip-memory-gate]
cpp/scripts/equiv.py compare  --format {cool,h5,text,bed,bedgraph,bedpe,tsv,image,hdf5-chic,bigwig,npz}
                              --class {E0,E1,E2,E3,E4,E5,E6,EN} A B
cpp/scripts/equiv.py report   [--out DIR] [--format {md,json}]
cpp/scripts/equiv.py list     [--tool NAME]
```

- `--cpp-bin` defaults to `cpp/build/tools`; `--py-python` defaults to the venv
  in the contract; `PYTHONPATH` is set to the repo root so the *repo* Python
  runs, not the installed 3.7.6 package.
- `--jobs` runs cases in parallel, each in its own temporary work directory
  under `$TMPDIR`, never under the repo. (Revised 2026-09-15, contract rule 13.)
  - **Default:** `--jobs auto`, half the physical cores.
  - **Memory:** cases run in parallel only while their expected peaks fit under
    a fraction of available memory.
  - **CPU time grows with parallel load,** so these run alone:
    - cases with a recorded time-gate ratio of at least 0.5;
    - cases without history that draw figures or belong to hicMergeDomains or
      hicTransform;
    - large unmeasured cases.
  - **Reruns:** a time-gate failure beside other cases is rerun alone with both
    sides fresh.
- `--cache {use,refresh,off}` reuses Python reference results keyed by case,
  arguments, input contents, reference code, package versions and harness
  environment. Merge verification uses `refresh` or `off`.
- `--noise-runs` (default 5) sets `N` for class EN and for the determinism check.
- `--skip-memory-gate` is for development only and marks the whole report
  `memory_gate: skipped`; a report with that flag can never record a `pass`.
- Exit code 0 only if every selected case passed.

### 9.2 Case definition

Cases live in `cpp/scripts/cases/<tool>.yaml` (one file per tool, parsed with a
30-line hand-rolled reader so there is no PyYAML dependency, or as JSON if that
is simpler). One case:

```yaml
- id: hicCorrectMatrix.KR.gm12878_chr1
  tool: hicCorrectMatrix
  tier: 3
  args: ["correct", "-m", "{data}/hicTADClassifier/gm12878_chr1.cool",
         "--correctionMethod", "KR", "--filterThreshold", "-1.5", "5",
         "-o", "{out}/kr.cool"]
  outputs:
    - path: "{out}/kr.cool"
      format: cool
      class: EN
  threads_arg: null          # or "--threads", triggers the 1-vs-16 determinism check
  modes: ["v3", "v4"]        # emits --compatMode v3 / v4 runs; v4 is diffed against v3
  expect_exit: 0
  memory:
    nnz_stored: 61804782
    nbins: 24926
    alpha: 1.2
    beta: 0.0
    budget_mb: 954
  large: true
  notes: "KR is nondeterministic in the Python reference; see PLAN.md 3.3, 5.7"
```

`{data}` expands to `hicexplorer/test/test_data`, `{out}` to the case work
directory. The runner executes the Python tool with `{out}` = `out_py` and the
C++ tool with `{out}` = `out_cpp`, then compares pairwise. `memory` may give
`budget_mb` directly or give `nnz_stored`, `nbins`, `alpha` and `beta` and let
the runner evaluate the formula of section 4.5; giving both is an error, so the
budget always has exactly one source of truth.

### 9.3 Comparator plugins

Each is a module in `cpp/scripts/comparators/` exposing
`compare(path_a, path_b, cls, opts) -> Result(passed, class_met, metrics, diffs)`.

- **`cool.py`** - opens both with `h5py`. Walks the object tree; compares the set
  of groups and datasets, then for each dataset: dtype (including the ENUM member
  map for `/bins/chrom`), shape, `maxshape`, chunk shape, the filter pipeline
  (id, name, `cd_values`), the fill value, and finally the decoded bytes. Compares
  root and group attributes, normalising `creation-date`, `generated-by`,
  `generated-by-cooler-lib` and `tool-url`. For E3/E4/EN, decodes
  `bin1_id`/`bin2_id`/`count` into a COO set and applies the section 5.5 rules.
  Handles `::/resolutions/<r>`, `::/cells/<n>` and the legacy `::/0` layout by
  comparing every group. Reports: `nnz_a`, `nnz_b`, `pattern_symdiff`,
  `max_abs`, `max_rel`, `corr`, and the first 20 differing pixels with their
  coordinates.
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
- **`npz.py`** - `scipy.sparse.load_npz` on both, then the section 5.5 rules.
- **`bam.py`** - only for `hicBuildMatrix --outBam`; compares via `pysam` on the
  header (normalising `@PG`) and every record's fields.
- **`noise.py`** - the class EN driver, and not a file-format comparator. Given
  the case, it runs the Python tool `N` times into separate directories, calls
  the format comparator for every pair to obtain the envelope `S` per numeric
  field, then applies the section 5.7 rules to the C++ output. Reports `S`, the
  median run's identity, the C++ deviation from it, and the ratio of the two.

### 9.4 Report format

`--out DIR` receives:

- `report.json`: `{harness_version, timestamp, git_commit, host, cases: [{id,
  tool, tier, class_declared, class_met, passed, py_seconds, cpp_seconds,
  py_peak_rss_kb, cpp_peak_rss_kb, budget_kb, memory_gate: {passed, ratio_to_budget,
  ratio_to_python}, determinism: {repeats, threads_1_vs_16, passed},
  noise_envelope (EN only), mode_delta (dual-mode only),
  outputs: [{path, format, passed, metrics, diffs}], stderr_py, stderr_cpp}]}`.
- `report.md`: a per-tier table (tool, cases, passed, class met, speedup, peak
  RSS, budget, headroom) followed by a section per failing case with the
  comparator's diff excerpt. A **memory summary table** at the top lists every
  tool sorted by `cpp_peak_rss / budget`, so the tools closest to their gate are
  visible without reading the whole report.
- `workdirs/<case-id>/` when `--keep-workdirs`, holding `out_py/`, `out_cpp/`,
  the exact command lines, and the raw stdout/stderr.

`equiv.py report` regenerates `report.md` from an existing `report.json` so the
comparison run does not have to be repeated.

### 9.5 What the harness must not do

It must not normalise anything not listed in 9.3, must not retry a failing case,
and must not have a per-case tolerance override. Tolerance lives in the class,
the class lives in the case file, and changing a class is a reviewed edit to
`STATUS.md`. The same applies to the memory budget: a case may not carry an
inline exemption, and raising a budget is a reviewed edit to `STATUS.md` with a
recorded reason.

## 10. Performance and memory measurement

Every harness run records, per case and per implementation:

- wall time from `time.perf_counter()` around the subprocess,
- **peak RSS from `/usr/bin/time -f "%e %M"` wrapping each invocation.** This is
  the authoritative figure and the one the section 8.3 gate uses.
  `resource.getrusage(RUSAGE_CHILDREN).ru_maxrss` is recorded as a cross-check
  but is not the gate, because it reports the maximum over *all* children since
  process start, which is wrong as soon as `--jobs > 1`,
- an RSS **trace**, sampled from `/proc/<pid>/status` `VmHWM` every 100 ms, for
  cases marked `large: true`. A single peak number says a tool exceeded its
  budget; the trace says where, which is what the implementing agent needs. The
  trace is written to `workdirs/<case-id>/rss_<impl>.tsv`,
- user and system CPU time, so that a "faster" result that merely burns more
  cores is visible,
- the thread count the tool was given.

Every C++ tool also reports its own peak RSS on `--verbose`, read from
`/proc/self/status` `VmHWM` at exit, so a developer does not need the harness to
see the number. The harness cross-checks the two and flags a discrepancy above
5 %, which usually means an allocation outside the tool's own accounting.

The report's per-tool row carries `speedup = py_seconds / cpp_seconds`,
`memory_ratio_python = cpp_peak_rss / py_peak_rss` (informational) and
`memory_headroom = 1 - cpp_peak_rss / budget` (the gate). Measurements are only
comparable when the machine is otherwise idle; the harness records the 1-minute
load average at the start of each case and marks a case `timing_unreliable` if it
exceeded 2.0. Gated cases run in parallel under the scheduling rules of section 9.1 (revised
2026-09-15).

The baseline numbers in section 4.2 were taken while the Python test suite was
running concurrently, so the wall times are upper bounds; the peak RSS figures
are unaffected by CPU contention and stand. They must be retaken on an idle
machine before any of them is quoted outside this document.

Targets, stated so that a miss is visible rather than rationalised:

| dimension | target |
|---|---|
| runtime, tier 1 (I/O bound) | 5x, largely from not paying 0.4 s of interpreter start-up |
| runtime, tier 3 | 3x |
| runtime, tiers 4 and 5 at 16 threads | 4x |
| runtime, `hicPCA` corrected mode | the Python takes 471 s on `mm9_reduced_chr1.cool`; `dsyevr` for 2 of 9,760 eigenpairs should be under 10 s, a 45x reduction |
| peak RSS, `hicCorrectMatrix` KR and ICE on `gm12878_chr1.cool` | 954 MB, a 9.6x and 9.5x reduction |
| peak RSS, `hicTransform --method pearson` on `Li_et_al_2015.h5` | 1,221 MB, a 4.2x reduction |
| peak RSS, `hicPCA` on `mm9_reduced_chr1.cool` | 2,509 MB compat, 946 MB corrected, against 4,070 MB |
| peak RSS, everything else | within the section 4.5 budget, no exceptions |

## 11. Principal risks

1. **KR has no deterministic oracle.** The reference disagrees with itself by up
   to 1.5e-4 relative on an 11 k-bin matrix and by 8.4e-3 on the normalisation
   factor of a 25 k-bin one (section 3.3). Class EN handles this, but it means
   KR can never be validated as tightly as anything else, and a real regression
   of the same magnitude as the oracle's noise would be invisible. Mitigation:
   the `v4` mode is deterministic and float64, so `v4`-against-`v4` regression
   testing across commits is exact even though `v3`-against-Python is not.
2. **The memory budgets are derived, not yet demonstrated.** Every figure in
   section 4.5 is a floor computed from the data layout plus an allowance. The
   first tool to implement, `hicInfo`, will show whether `C = 64 MB` is
   realistic for a C++ process with HDF5, blosc and OpenBLAS linked in. If it is
   not, every budget in the plan shifts and the table must be recomputed rather
   than individually relaxed.
3. **`dgeev` column order in `hicPCA`.** If linking the oracle's OpenBLAS does
   not reproduce the eigenvector order, `hicPCA` compatibility mode cannot be
   made equivalent, and the tool ships as `v4`-only with a recorded deviation.
4. **PyTables blosc chunking.** Comparing h5 at value level (L3) rather than
   structurally (L2) means a C++-written `.h5` will differ from a Python-written
   one byte-wise. Any downstream consumer that checksums `.h5` files will notice.
5. ~~**The h5 blosc filter must actually work.**~~ **Retired 2026-09-01,
   commit bd3b4bff.** The C++ blosc compress path produces files PyTables 3.10.1
   opens as a `CArray` reporting
   `Filters(complevel=5, complib='blosc', shuffle=True)`, with `cd_values`
   byte-for-byte what PyTables emits, and Python reads a C++-written h5 back to
   identical `hicInfo` output including the exact float sum. This was the risk
   most likely to block tier 0 and it is gone.
6. **`fit_nbinom`'s L-BFGS-B.** A different optimiser stopping point propagates
   into `hicDetectLoops` and the entire cHi-C background model. E4/E5 classes
   absorb it, but if the loop call sets diverge beyond the Jaccard floor the
   only remaining option is to vendor scipy's exact L-BFGS-B Fortran translation.
7. **The rule 2 exemption list may be too long.** Five tools materialise both
   triangles, and their budgets are correspondingly loose (`alpha = 2.2`). If
   `hicFindTADs` or `hicDetectLoops` turns out to need only banded access, the
   exemption should be withdrawn and the budget tightened; leaving it unexamined
   would quietly forfeit half the memory win on the heaviest tools.
8. **Characterization tests must exist before the port.** Five tools have no
   Python test at all and 22 more have tests that cannot fail (`STATUS.md`).
   Porting one of them without first pinning its behaviour means the port defines
   the behaviour, which is the one outcome the contract forbids.
