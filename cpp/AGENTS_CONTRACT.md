# HiCExplorer v4 (C++) - working contract

This file is the single source of truth for environment facts, layout and rules.
Both the supervising agent and the implementing agent must read it before acting
and must keep it accurate.

## 1. Repository and branch

- Upstream repo: `~/src/HiCExplorer` (branch `z-score` is the user's working
  checkout, **do not touch it**).
- v4 work happens in the git worktree `~/src/HiCExplorer-v4`, branch
  `version4-cpp`, forked from `master` (`c2ad8630`).
- Python reference version in this tree: 3.7.7-dev (47 tools, ~18k lines).

## 2. Environments

### Python reference oracle (never modify)

Dependencies live in the shared conda env
`~/miniconda3/envs/__hicexplorer@3.7.6`. Conda/mamba installs into it are
forbidden (it is used by Galaxy jobs and the local conda solver is broken).

A throwaway venv layered on top of it provides pytest without polluting it:

```
SP=/tmp/claude-1000167729/-home-mh-hannover-local-wolffjoa/3de46a4d-41bb-4fb9-aea9-26cceef9ad01/scratchpad
$SP/hicx-venv/bin/python        # python 3.12, all HiCExplorer deps + pytest 9
```

Run a reference tool on the v4 tree like this (the repo code, not the installed
3.7.6 package):

```
cd ~/src/HiCExplorer-v4
PYTHONPATH=$PWD $SP/hicx-venv/bin/python bin/hicInfo -m <matrix>
```

Run the Python test suite the same way:

```
cd ~/src/HiCExplorer-v4
PYTHONPATH=$PWD $SP/hicx-venv/bin/python -m pytest hicexplorer/test/general/test_hicInfo.py -q
```

### C++ build dependencies

No system dev packages exist. All native dependencies are taken from the conda
env used as a plain prefix:

```
HICX_DEPS=~/miniconda3/envs/__hicexplorer@3.7.6
```

It provides: HDF5 (C and C++ API), zlib, bzip2, htslib, Eigen 3, BLAS/LAPACK.
Configure with `-DCMAKE_PREFIX_PATH=$HICX_DEPS` and make sure the produced
binaries find the runtime libs (RPATH to `$HICX_DEPS/lib`).

Toolchain: g++ 13.3 (C++20 available), cmake 4.4, 32 cores. Network reaches
GitHub and PyPI through the corporate proxy, so CMake `FetchContent` works.

## 3. Directory layout (fixed - do not invent alternatives)

```
cpp/
  CMakeLists.txt        top-level build
  cmake/                helper modules
  core/                 libhicx4: matrix model, file formats, math kernels
    include/hicx/       public headers
    src/
  tools/                one executable per HiCExplorer tool
  tests/                C++ unit tests
  scripts/              python-vs-C++ equivalence harness
  PLAN.md               architecture + porting roadmap (supervisor owns)
  STATUS.md             per-tool progress ledger (supervisor owns)
  AGENTS_CONTRACT.md    this file
```

`.gitignore` at repo root already ignores `build`, `bin`, `lib`, `var`, `dist`,
`*.so`. Therefore: build inside `cpp/build` (ignored, good) and **never name a
source directory** `bin`, `lib`, `var`, `dist`, `eggs` or `parts` - such
directories would silently not be committed.

## 4. Rules

1. **Test first where coverage is missing.** Before porting a tool that has no
   Python test (or no test for the code path being ported), add a
   characterization test to `hicexplorer/test/` that pins the current Python
   behaviour on the real test data. That test is committed on `version4-cpp`
   before the C++ code that replaces it.
2. **Validate on real data, not toys.** Equivalence runs use the repository's
   real matrices (`hicexplorer/test/test_data/`, e.g. `Li_et_al_2015.cool`,
   `hicTADClassifier/gm12878_chr1.cool` at 79 MB, the 100 kb
   GSM2644945/GSM2644947 cool and h5 matrices) and real BAM input where the
   tool consumes alignments. Hand-written 5x5 matrices are for unit tests only
   and never count as validation.
3. **Equivalence is the acceptance criterion.** A ported tool is done when its
   output matches the Python tool on real input. Numeric tolerance policy is
   defined in `PLAN.md`; it must be justified per tool, not chosen ad hoc, and
   any deliberate deviation must be recorded in `STATUS.md`.
4. **No silent scope cuts.** If a tool cannot be ported faithfully (for example
   because it depends on scikit-learn, hyperopt or matplotlib), record that in
   `STATUS.md` with the reason and the proposed strategy. Do not quietly drop
   it and do not fake equivalence.
5. **Git.** Agents never run git write commands (`add`, `commit`, `checkout`,
   `stash`, `push`). Several agents share this one worktree, so a concurrent
   index write would corrupt another agent's work. Leave finished work in the
   working tree and report it; the orchestrating session commits, one commit per
   completed task, only after it is verified working. Nothing is ever pushed.
6. **Style.** English only, ISO dates, no em-dashes, no emojis in code,
   comments, docs or commit messages.
