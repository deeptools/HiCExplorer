# v4 optimization contract: threading, SIMD, loop structure

Requested by the project owner 2026-09-01: multithread where possible, use SSE
and AVX wherever possible, flatten loops. This file is the binding design
constraint for that work. It sits alongside `PLAN.md`, which stays the
architecture, and `AGENTS_CONTRACT.md`, which stays the rules.

## 1. Why this is possible now

Under the original policy every tool had to be bit-identical to Python, which
forbade almost all of this: SIMD and multithreading both change floating-point
reduction order, and reduction order changes the last bits. The acceptance gate
set on 2026-09-01 (`PLAN.md` 5.0, three significant digits per item, relative
`1e-3`) removes that obstacle. Vectorised and threaded reductions are now
allowed to differ from numpy, because `1e-3` is roughly twelve orders of
magnitude looser than the difference a reordered float64 reduction produces.

This does **not** mean results may drift. See section 3.

## 2. Hardware target and dispatch

Development machine: AMD Ryzen 9 7950X, Zen 4, 16 cores and 32 threads, with
`avx`, `avx2`, `fma`, `bmi2`, `sse4_1`, `sse4_2`, `ssse3` and AVX-512
(`avx512f`, `avx512bw`, `avx512cd`, `avx512dq`, `avx512vl`, `avx512ifma`,
`avx512vbmi`).

**Do not build with `-march=native`.** This project is meant to ship through
bioconda, where the binary must run on any x86-64 machine. The rule:

- Compile the library at a **baseline of x86-64-v2** (SSE4.2, popcnt), which is
  safe on anything from the last fifteen years.
- Provide **AVX2 + FMA** and **AVX-512** variants of the hot kernels only,
  selected at **runtime** with `__builtin_cpu_supports` (GCC function
  multiversioning, or an explicit function-pointer table per kernel).
- Every kernel must have a scalar reference implementation that is always
  compiled and always correct. The SIMD paths are optimisations of it, and there
  is a unit test asserting each SIMD path agrees with the scalar one.
- Note on this specific chip: Zen 4 implements AVX-512 on a 256-bit datapath, so
  the gain over AVX2 comes from masking, wider registers for gather and better
  instruction selection, not from doubled throughput. Measure before assuming
  AVX-512 is worth a separate code path; if it is not faster than AVX2 for a
  given kernel, do not ship a third path for it.

## 3. Determinism is not negotiable

Loosening the comparison against Python does **not** loosen our own
reproducibility. Two hard requirements, both gated by the harness:

1. **Run-to-run identical.** Two runs of the same C++ binary on the same input
   must produce byte-identical output. No reduction may depend on thread
   scheduling, completion order or work-stealing.
2. **Thread-count invariant.** Output at `--threads 1` and at `--threads 16`
   must be byte-identical. This is the `-tN == -t1` rule the falcoAmadeus work
   established, and it is the difference between a fast tool and a trustworthy
   one.

The way to satisfy both: partition by a **fixed** index range (chromosome,
diagonal, bin1 chunk), reduce **sequentially within** a partition, and combine
partitions in **index order**. Never accumulate into a shared float, never use
an atomic add on a double, never let the number of threads decide how the sum is
grouped. A parallel reduction whose result depends on the thread count is a bug
even when it passes the ED gate.

The same applies to SIMD: a vectorised reduction must use a fixed number of
accumulator lanes combined in a fixed order, so that its result depends only on
the input length, never on the runtime dispatch choice. Since the scalar, AVX2
and AVX-512 paths will generally produce different last bits from each other, a
kernel's dispatch choice must not vary within a run, and the unit test asserting
they agree uses the ED tolerance, not equality.

## 4. What is worth optimising, in order

Ranked by measured cost in the Python reference, not by guesswork. The figures
are from the runs recorded in `PLAN.md` 4.2 and `STATUS.md`.

1. **`hicCorrectMatrix` ICE**, 711 s and 9,064 MB on `gm12878_chr1.cool`. The
   iteration is 50 passes of a sparse row reduction and an elementwise scale.
   Both vectorise, and the pass over rows partitions cleanly. Largest single
   win available anywhere in the project.
2. **`hicCorrectMatrix` KR**, 9,199 MB. The upstream krbalancing wraps its whole
   OpenMP loop body in `#pragma omp critical`, so it is serial with added
   contention, and it does two binary searches per stored nonzero through
   `x.coeff(row, 0)`. Replace both.
3. **`hicTransform --method pearson`**, 26 s and 5,150 MB on an 11k-bin matrix.
   The dense centring and inner products are exactly what `dsyrk` is for; use
   the BLAS in the conda prefix rather than a hand-rolled loop, and thread over
   chromosome blocks.
4. **`hicPCA`**, 471 s on a 722 KB input. Dominated by the eigensolver;
   `dsyevr` replaces `dgeev` and is both faster and lower memory.
5. **The CSR kernels already written** in `matrix_ops` and `sparse_matrix`:
   elementwise add, subtract, multiply, divide, `log2`, the marginal row sums,
   the symmetric streaming sum. These run inside every tool, so a gain here is
   a gain everywhere. They are the right first target because they are already
   covered by 84 harness cases and 53 unit tests, so a regression is caught
   immediately.
6. **The file layer**: HDF5 chunk decompression is single-threaded per dataset
   today and is the bulk of the time in the read-dominated tools. Decompressing
   independent chunks in parallel is straightforward and order-independent.

## 5. Loop structure

- **Flatten the two-level CSR iteration** (`for row { for k in row }`) into a
  single pass over the value array wherever the row index is not needed, which
  is the case for every elementwise operation. This is the single most effective
  change for auto-vectorisation, because it removes the inner loop's variable
  trip count.
- Hoist the row-boundary lookup out of inner loops; do not call `at()` or any
  binary search inside a loop body.
- Prefer `std::span` over pointer plus length pairs, and mark the hot kernels'
  arrays `__restrict` so the compiler may assume no aliasing.
- Structure of arrays, not array of structures, for anything the SIMD paths
  touch.
- Check what the compiler already did before writing an intrinsic:
  `-fopt-info-vec-missed` on the kernel in question. A hand-written intrinsic
  that matches what GCC already emits is a maintenance cost with no benefit.

## 6. Measurement is mandatory

No optimisation is accepted on the argument that it should be faster.

- Report **CPU time (user plus sys)**, not wall clock. This machine runs
  unrelated training jobs continuously and the load average has been between 13
  and 57 throughout the project, so wall clock is not reproducible. The harness
  gates on CPU time.
- Report the before and after for every kernel changed, on **real** matrices.
- Report peak RSS as well, since several of the SIMD-friendly layouts cost
  memory and the budgets in `STATUS.md` are a hard gate.
- If a change does not measurably help, revert it. Complexity that buys nothing
  is a defect.
