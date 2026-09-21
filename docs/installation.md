# Installation

!!! warning "This is the C++ rewrite, not the Python package"
    The instructions below build **this** repository: a from-scratch C++ rewrite of HiCExplorer. They
    replace the Python installation methods documented for the original tool (`conda install
    hicexplorer`, `pip install hicexplorer`, the Galaxy Tool Shed wrapper, the
    `quay.io/bgruening/galaxy-hicexplorer` Docker image). Those install the *Python* HiCExplorer, which
    is a separate, unrelated codebase from this one; running them does not give you the tools described
    in this site. As of now this C++ rewrite has no conda (bioconda), Galaxy or Docker packaging of its
    own yet: building from source is the only supported way to install it.

## Requirements

- A C++17 (or newer) compiler. The reference build uses g++ 13.3.
- [CMake](https://cmake.org/) 3.x, and Ninja or Make.
- Native dependencies, most easily obtained as a conda environment used as a plain prefix (not
  activated, just pointed at with `-DCMAKE_PREFIX_PATH`): HDF5 (C and C++ API), zlib, bzip2, htslib,
  Eigen 3, and a BLAS/LAPACK implementation.
- Two sibling libraries built specifically for this project, fetched automatically by CMake or pointed
  at local checkouts: [coolercpp](https://github.com/) (cool/mcool file I/O) and
  [hicfilecpp](https://github.com/) (`.hic` file I/O).
- For the handful of plotting tools that still delegate figure drawing to the original Python code (see
  [Tools > Visualization](tools/visualization.md)): a Python interpreter with matplotlib and
  pyGenomeTracks, referenced through the `HICX_PLOT_PYTHON` environment variable. Every other tool has
  no Python runtime dependency.
- `hicMergeDomains` additionally needs the Graphviz `dot` binary on `PATH` at run time, to render the
  TAD hierarchy plots it produces.

## Build from source

```bash
git clone <repository-url> HiCExplorer-v4
cd HiCExplorer-v4

cmake -B cpp/build -S cpp \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=/path/to/conda/env/used/as/prefix \
  -DHICX_COOLERCPP_GIT_REPOSITORY=/path/to/coolercpp \
  -DHICX_HICFILECPP_GIT_REPOSITORY=/path/to/hicfilecpp

cmake --build cpp/build -j
```

- `-DCMAKE_PREFIX_PATH` should point at a conda environment (or any prefix) providing HDF5, zlib,
  bzip2, htslib, Eigen 3 and BLAS/LAPACK; it is used as a plain library search path, not activated.
- `-DHICX_COOLERCPP_GIT_REPOSITORY` / `-DHICX_HICFILECPP_GIT_REPOSITORY` point CMake at local checkouts
  of the two sibling libraries so they do not need to be fetched over the network; each is otherwise
  resolved through `find_package` if already installed, or fetched and built from a pinned commit.
  `-DHICX_COOLERCPP_FORCE_FETCH=ON` / `-DHICX_HICFILECPP_FORCE_FETCH=ON` skip an already-installed copy
  and force the fetch instead.

Each tool builds as its own executable under `cpp/build/tools/`, for example `cpp/build/tools/hicBuildMatrix`,
`cpp/build/tools/hicFindTADs`, and so on; there is no single combined binary or Python entry-point script.
Add `cpp/build/tools` to `PATH`, or invoke the binaries by their full path.

## Verify the installation

```bash
cpp/build/tools/hicInfo --version
```

Every tool accepts `--help` (a few require their required arguments to be present before `--help` is
honored; see the individual tool's page if `--help` alone reports missing arguments instead of usage).

## Running the test suites

This is a development repository; the equivalence-test harness that validates the C++ tools against the
Python reference implementation is documented in `cpp/AGENTS_CONTRACT.md` and `cpp/STATUS.md` in the
repository, not repeated here since it is aimed at contributors to this rewrite rather than at users of
the built tools.
