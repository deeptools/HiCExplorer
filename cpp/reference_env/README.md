# Pinned reference environment for the HiCExplorer 3.7.6 equivalence tests

The equivalence tests compare the C++ port with the upstream HiCExplorer 3.7.6
Python implementation (extracted by `cpp/scripts/reference_tree.py`). Outputs
such as reprs, pandas frames and matplotlib PNG bytes depend on the versions of
numpy, pandas, matplotlib, freetype and the rest of the dependency stack, so the
reference has to run on exactly the versions it was validated with
(python 3.12.7, numpy 1.26.4, pandas 2.2.3, matplotlib-base 3.8.4, freetype
2.12.1, ...). With numpy 2, pandas 3 or a newer matplotlib, 63 API cases fail
on representation and figure-byte differences alone.

## Contents

| File | Purpose |
|---|---|
| `pins.txt` | name, version and build string of all 280 packages (linux-64) |
| `make_recipe.py` | writes `recipe/recipe.yaml` from `pins.txt` |
| `recipe/recipe.yaml` | rattler-build recipe of the noarch metapackage `hicexplorer-equivalence-reference` 3.7.6.1 |
| `hicx-ref-3.7.6-linux-64.lock.txt` | explicit lock file (`conda list --explicit --md5`) of an environment created from the package |

The package contains no files. Its run requirements pin every package of the
original `__hicexplorer@3.7.6` environment (including bioconda `hicexplorer`
3.7.6 itself, which the environment held; the reference imports the extracted
tree, not that installation) to exact version and build string, plus
`pytest 8.3.4`, `pytest-xdist 3.6.1` and their dependencies (`execnet`,
`iniconfig`, `pluggy`, `tomli`). The suite imports pysam, pybedtools, pybigwig,
tables, cooler, h5py and psutil; all are in the pin list.

## Create the environment

From the lock file, which needs no solving (fastest, and the fallback if the
solver is unavailable):

    mamba create -n hicx-ref-3.7.6 --file cpp/reference_env/hicx-ref-3.7.6-linux-64.lock.txt

From the package (the solver can only return the pinned set; the package must
be available in a channel, e.g. the local one built below):

    mamba create -n hicx-ref-3.7.6 --override-channels \
        -c /path/to/local-channel -c conda-forge -c bioconda \
        hicexplorer-equivalence-reference

The local conda solver of the lab workstation is broken (libmamba plugin);
use `mamba` or `micromamba`, never the base environment's `conda create`.

## Build the package

rattler-build (0.72.2 here) is used, because conda-build is not usable in the
base installation:

    rattler-build build --recipe cpp/reference_env/recipe/recipe.yaml \
        --output-dir /path/to/local-channel -c conda-forge -c bioconda

`--output-dir` becomes a channel (it contains `noarch/repodata.json`). To
change a pin, edit `pins.txt`, run `python3 cpp/reference_env/make_recipe.py`,
bump `VERSION` in `make_recipe.py`, rebuild, and regenerate the lock file from a
fresh environment:

    mamba create -n hicx-ref-new ...            # as above
    conda list -n hicx-ref-new --explicit --md5 | grep -v '^file:' > cpp/reference_env/hicx-ref-3.7.6-linux-64.lock.txt

Nothing here is submitted to bioconda or any remote channel.

## Point the tests at the environment

`HICX_REFERENCE_ENV` is the conda prefix. Both consumers read it:

* `cpp/scripts/equiv.py`: interpreter of the reference processes becomes
  `$HICX_REFERENCE_ENV/bin/python3`, and `HICX_DEPS` (graphviz `dot`, etc.)
  defaults to the same prefix. `HICX_PY_PYTHON` and `HICX_DEPS` still override
  each part individually, and `--py-python` overrides on the command line.
* `hicexplorer/test/api/reference.py` (Python API equivalence cases):
  interpreter becomes `$HICX_REFERENCE_ENV/bin/python3`;
  `HICX_REFERENCE_PYTHON` overrides.

Without any of these variables the old path
`~/miniconda3/envs/__hicexplorer@3.7.6` is used.

Run the tests with the interpreter of the same environment, so the comparators
and pytest see the pinned h5py, PyTables and numpy:

    export HICX_REFERENCE_ENV=$HOME/miniconda3/envs/hicx-ref-3.7.6
    export HICX_BIN_DIR=$PWD/cpp/build/tools PYTHONPATH=$PWD
    $HICX_REFERENCE_ENV/bin/python -m pytest hicexplorer/test/general hicexplorer/test/api
    $HICX_REFERENCE_ENV/bin/python cpp/scripts/equiv.py run --tool hicInfo

The API reference results are cached by interpreter path, so a new environment
starts with an empty cache (the first full API run takes longer). The equiv.py
cache is keyed on the interpreter and package versions, so an environment with
the same versions reuses entries; add `--cache off` to force fresh reference
runs. `hicMergeDomains` needs graphviz `dot` on PATH when the tests start.
Without it the package side raises `ExecutableNotFound`, which fails against a
cached reference result and, with an empty cache, passes vacuously because the
reference raises the same error. Run
`export PATH=$HICX_REFERENCE_ENV/bin:$PATH` or activate the environment.

## Pinning limits

Everything is pinned to exact build strings from conda-forge and bioconda. If
a channel later removes one of these builds, the lock file fails at that URL;
keep a copy of the downloaded packages (`$CONDA_PKGS_DIRS`) to rebuild
offline. Some packages (libarrow and pyarrow) are CUDA variants because that
is what the original environment held.

## Build environment of the C++ tree

The C++ tree is built against a separate environment, `hicx-build`, not against
this reference environment. Its package list (the host and run requirements for
a bioconda recipe) is in `cpp/BUILDING.md`.
