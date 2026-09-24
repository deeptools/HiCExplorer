# Installation

## conda

```bash
conda install -c conda-forge -c bioconda hicexplorer
```

The package contains the command-line tools and the Python API. Every tool is started with its own
name (`hicBuildMatrix`, `hicFindTADs`, and so on).

## pip

```bash
pip install hicexplorer
```

The wheels contain the tools and the libraries they need for Linux (x86_64 and arm64) and macOS (arm64).

## Requirements for building from source

- A C++20 compiler and [CMake](https://cmake.org/) 3.24 or newer, and Ninja or Make.
- Native dependencies, most easily obtained as a conda environment used as a plain prefix (not
  activated, just pointed at with `-DCMAKE_PREFIX_PATH`): HDF5, blosc, zlib, htslib and BLAS/LAPACK.
- The libraries [coolercpp](https://github.com/) (cool/mcool file I/O) and
  [hicfilecpp](https://github.com/) (`.hic` file I/O), fetched automatically by CMake or pointed at
  local checkouts.
- The plotting tools draw figures with Python. The Python environment needs matplotlib 3.8.4 and, for
  `hicPlotTADs`, pyGenomeTracks. The environment is selected with the `HICX_PLOT_PYTHON` environment
  variable and defaults to the Python that runs the tool.
- `hicMergeDomains` additionally needs the Graphviz `dot` binary on `PATH` at run time, to render the
  TAD hierarchy plots.

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

- `-DCMAKE_PREFIX_PATH` should point at a conda environment (or any prefix) providing HDF5, blosc, zlib,
  htslib and BLAS/LAPACK; it is used as a plain library search path, not activated.
- `-DHICX_COOLERCPP_GIT_REPOSITORY` / `-DHICX_HICFILECPP_GIT_REPOSITORY` point CMake at local checkouts
  of the two sibling libraries so they do not need to be fetched over the network; each is otherwise
  resolved through `find_package` if already installed, or fetched and built from a pinned commit.
  `-DHICX_COOLERCPP_FORCE_FETCH=ON` / `-DHICX_HICFILECPP_FORCE_FETCH=ON` skip an already-installed copy
  and force the fetch instead.

Each tool builds as its own executable under `cpp/build/tools/`, for example `cpp/build/tools/hicBuildMatrix`,
`cpp/build/tools/hicFindTADs`, and so on. Add `cpp/build/tools` to `PATH`, or invoke the binaries by their
full path. To install the Python package together with the executables, run `pip install .` in the
repository root.

## Verify the installation

```bash
cpp/build/tools/hicInfo --version
```

Every tool accepts `--help` (a few require their required arguments to be present before `--help` is
honored; see the individual tool's page if `--help` alone reports missing arguments instead of usage).

## Running the test suites

The equivalence-test harness that checks every tool against HiCExplorer 3.7.6 is documented in
`cpp/AGENTS_CONTRACT.md` and `cpp/STATUS.md` in the repository. It is aimed at contributors.
