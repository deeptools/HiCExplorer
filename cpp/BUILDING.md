# Building the C++ tree

`cmake` configure and build need no network access. Everything the tree
compiles is either in the tree (`core/io/` for coolercpp, hicfilecpp and
hic2cool; `third_party/` for x86-simd-sort; see `THIRD_PARTY_NOTICES.md`) or
comes from the dependency prefix given with `-DCMAKE_PREFIX_PATH`. The one
exception is scipy's Cephes headers, see "Remaining download" below.

## Build requirements

The prefix is a conda environment used as a plain prefix (never activated).
The packages and version ranges below are the ones the output of the tools was
verified with; they can be copied into the `host` and `run` sections of a
bioconda recipe.

host (build time):

| Package | Version | Used for |
|---|---|---|
| `cmake` | >=3.24 | build system |
| `gxx_linux-64` (or the platform compiler) | C++20, g++ 13 or newer | compiler |
| `hdf5` | >=1.14.3,<1.15 | cool, mcool, h5 matrices |
| `blosc` | >=1.21.6,<1.22 | blosc-compressed PyTables datasets |
| `libzlib` | >=1.3.1,<2 | gzip, .hic blocks, bigWig |
| `htslib` | >=1.21,<1.22 | BAM reading and writing |
| `libopenblas` | >=0.3.28,<0.4 | LAPACK (the library numpy and scipy use) |
| `libblas`, `liblapack` | 3.9.0 (openblas build) | LAPACK symbols |
| `eigen` | >=3.4.0,<3.5 | Spectra (hicPCA --eigenSolver lanczos) |
| `spectra-cpp` | 1.2.0 | Spectra headers |
| `libboost-headers` | >=1.85,<1.86 | Boost.Math `ibeta` (scipy.special.betainc) |
| `libbigwig` (bioconda) | 0.4.8 | bigWig reading and writing |
| `pybind11` | 2.13.6 | the Python extensions (`-DHICX_BUILD_PYTHON=ON`) |
| `python` | >=3.10 | the Python extensions |
| `doctest` | 2.4.11 | unit tests (`-DHICX_BUILD_TESTS=ON`, the default) |

run:

| Package | Version | Used for |
|---|---|---|
| `hdf5`, `blosc`, `libzlib`, `htslib`, `libopenblas`, `libbigwig` | as above, pinned by the conda run exports | shared libraries the executables link |
| `libcurl` | pulled in by `libbigwig` | |
| `graphviz` | any | `dot` for `hicMergeDomains` |
| `python`, `matplotlib-base`, `numpy` | see `requirements.txt` | plotting tools and the Python package |

`libboost-headers` is pinned to 1.85 because the Boost.Math headers `betainc`
includes are, for 1.85.0, identical to the ones scipy 1.14.1 vendors (one
`noexcept` specification of an unrelated `long double` specialisation
differs); the CMake configure step refuses another Boost version. Spectra,
pybind11, doctest and the libBigWig headers were compared with the sources the
tree was verified with and are byte-identical.

A build environment with exactly these packages:

    mamba create -n hicx-build -c conda-forge -c bioconda --override-channels \
        python=3.12 pybind11=2.13.6 libboost-headers=1.85 libbigwig=0.4.8 \
        spectra-cpp=1.2.0 eigen=3.4.0 doctest=2.4.11 hdf5=1.14.3 blosc=1.21.6 \
        libzlib=1.3.1 htslib=1.21 libopenblas=0.3.28 \
        libblas=3.9.0=25_linux64_openblas liblapack=3.9.0=25_linux64_openblas \
        graphviz

## Configure and build

    cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_PREFIX_PATH=$HOME/miniconda3/envs/hicx-build \
        -DHICX_BUILD_PYTHON=ON -DPython_EXECUTABLE=$HOME/miniconda3/envs/hicx-build/bin/python
    cmake --build cpp/build -j
    ctest --test-dir cpp/build

## Remaining download

`core/CMakeLists.txt` still downloads scipy's Cephes headers (scipy v1.14.1,
`scipy/special/special/cephes`, 14 files, each checked against a SHA256) into
`<build>/_deps/scipy_special` when they are not there already. They are used by
`hicx/scipy_stats.hpp`. With the files in that directory, or without network
access after the first configure of a build directory, nothing is downloaded.
