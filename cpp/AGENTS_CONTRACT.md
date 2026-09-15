# HiCExplorer v4 (C++) - working contract

This file is the single source of truth for environment facts, layout and rules.
Both the supervising agent and the implementing agent must read it before acting
and must keep it accurate.

## 1. Repository and branch

- Upstream repo: `~/src/HiCExplorer` (branch `z-score` is the user's working
  checkout, **do not touch it**).
- v4 work happens in the git worktree `~/src/HiCExplorer-v4`, branch
  `version4-cpp`, forked from `master` (`c2ad8630`).
- Parallel ports run in separate worktrees under
  `~/src/HiCExplorer-v4-worktrees/`, each on its own branch created from
  `version4-cpp`, and the orchestrating session verifies and merges them. Never
  create a worktree or any other directory inside `~/src/HiCExplorer`.
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

Two own libraries are C++ dependencies, each pinned to a commit in
`cpp/cmake/HicxCoolercpp.cmake` and `cpp/cmake/HicxHicfilecpp.cmake`:
- coolercpp, found through `find_package`, or fetched with
  `-DHICX_COOLERCPP_GIT_REPOSITORY=$HOME/src/coolercpp`;
- hicfilecpp, found the same way, or fetched with
  `-DHICX_HICFILECPP_GIT_REPOSITORY=$HOME/src/hicfilecpp`.

`-DHICX_COOLERCPP_FORCE_FETCH=ON` and `-DHICX_HICFILECPP_FORCE_FETCH=ON` skip an
installed package. Moving a pin is a reviewed commit.

The `.hic` oracles are kept out of the reference venv:
- hicstraw 1.3.1 and hictkpy live in `$SP/hicfile-oracle-venv`;
- Juicer tools 1.22.01 and 2.20.00 jars live in `$SP/juicer/` and run on
  `~/miniconda3/bin/java` (OpenJDK 11).

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
  PLAN.md               architecture + porting roadmap (orchestrating session owns)
  STATUS.md             per-tool progress ledger (orchestrating session owns)
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
   any deliberate deviation must be recorded in `STATUS.md`. How a tool reaches
   its result is free: data structures and algorithms may differ from the
   Python wherever the harness shows the same result (project owner,
   2026-09-13; `PLAN.md` 5.0 item 3). Reproduce the Python's route only where
   the result depends on it.
4. **No silent scope cuts.** If a tool cannot be ported faithfully (for example
   because it depends on scikit-learn, hyperopt or matplotlib), record that in
   `STATUS.md` with the reason and the proposed strategy. Do not quietly drop
   it and do not fake equivalence.
5. **Git.** It depends on where the agent works.
   - **Shared worktree** `~/src/HiCExplorer-v4`: the agent runs no git write
     command (`add`, `commit`, `checkout`, `stash`, `push`), because a
     concurrent index write would corrupt another agent's work. Finished work
     stays in the working tree; the orchestrating session commits it, one
     commit per completed task, only after verifying it.
   - **Own worktree** under `~/src/HiCExplorer-v4-worktrees/`: the agent commits
     its verified work on that worktree's branch, one commit per tool, and
     reports the hashes; the orchestrating session verifies and merges.
   - **Always:** never push, never touch another branch or worktree, never
     reset or rewrite history. (Revised 2026-09-13 for the worktree model.)
6. **Style.** English only, ISO dates, no em-dashes, no emojis in code,
   comments, docs or commit messages.
7. **Requested output files.** A tool never exits 0 without writing every file
   the user explicitly asked for.
   - **Figures:** C++ computes the data and `hicexplorer_plot` draws it with
     the Python's matplotlib calls (tier 7 option a). The interpreter is
     `HICX_PLOT_PYTHON`, with matplotlib 3.8.4.
   - **Wrong or missing environment:** the tool checks it right after parsing
     and exits 3 before reading any input or creating any output, so the output
     directory is left unchanged.
   - **Parse-time outputs:** outputs that argparse's `FileType('w')` would
     create while parsing are created only after that check.
   - **Tier 8:** tools (ML) still await the project owner's decision; report
     any figure or model they require before implementing.
   - (Added 2026-09-13; revised 2026-09-15 when tier 7 landed.)
8. **Waiting on jobs.** An agent does not end its turn while a job it started
   is still running. A detached job does not wake a stopped agent. Run long
   commands in the foreground, or block on a wait loop within the same turn.
   (Added 2026-09-13, after an agent stopped mid-run.)
9. **Counting verdicts.** `cpp/scripts/equiv.py` prints `pass` in lowercase and
   `FAIL` in uppercase. Always count verdicts case-insensitively. A filter that
   matched only lowercase `fail` once hid three failures. (Added 2026-09-13.)
10. **Committed inputs.** Before committing, every file a test or harness case
    reads must be tracked. Run `cpp/scripts/check_case_inputs.py <repo> <rev>`
    on the revision, which fails when a case expecting success reads an
    untracked input. An earlier commit shipped three cases whose BED inputs
    existed only untracked in one worktree, so its verification passed there
    and every fresh checkout failed. (Added 2026-09-13.)
11. **Clean-export build.** Verify by building a clean export of the revision
    outside its source tree and running ctest, not only in the worktree that
    produced it. A unit test once wrote into a `cpp/build` directory inside the
    source tree and passed only where that directory happened to exist.
    (Added 2026-09-13.)
12. **No orphaned processes.** Before reporting, check that no process started
    by the task is still running, including mutation-test and noise runs. An
    agent cut off by a usage limit once left a mutation-test run that consumed a
    CPU core and 4.4 GB for eleven days. (Added 2026-09-13.)
13. **Reference cache and parallel runs.**
    - **Day to day:** regressions run `equiv.py run` with the defaults
      `--cache use` and `--jobs auto`.
    - **Merge verification:** runs `--cache refresh` (or `--cache off`) on a
      clean export, with `--jobs auto`.
    - **Runs alone:**
      - cases with a recorded time-gate ratio of at least 0.5;
      - cases without a C++ measurement in the cache that draw figures or
        belong to hicMergeDomains or hicTransform;
      - large cases without any measurement.
    - **Time-gate reruns:** a case that fails only its time gate while other
      cases ran beside it is rerun alone with both sides fresh, and that verdict
      counts; `report.json` records both attempts.
    - **Failures:** a case that fails alone, or fails any other gate, fails.
    - **A clean machine:** a verdict counts only from a run with nothing else,
      such as a tracer or another heavy job, loading the machine.
    - (Added 2026-09-15: C++ CPU time on the development machine grows with
      parallel load, by a median of 1.13 to 1.44 at 8 to 16 jobs on a
      thin-margin sample.)
14. **Parallel agents.** Up to three implementing agents may run at once, each a
    fresh agent with a self-contained brief.
    - **Isolation:** each works in its own worktree under
      `~/src/HiCExplorer-v4-worktrees/`, uses its own harness cache directory
      and at most `--jobs 4`, and starts no subagents.
    - **Load:** CPU time grows with the other agents' load, so a time-gate
      failure seen under load is rechecked on a quiet machine before a merge.
    - **Disk:** keep scratch in a per-task directory, and delete large
      intermediates before reporting.
    - (Approved by the project owner 2026-09-15.)
