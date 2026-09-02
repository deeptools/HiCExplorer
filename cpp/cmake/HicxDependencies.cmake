# Native dependencies of HiCExplorer v4.
#
# There are no system development packages on the target machine. All native
# dependencies are taken from the conda environment that also provides the
# Python reference implementation, used here as a plain prefix:
#
#   cmake -S cpp -B cpp/build -DCMAKE_PREFIX_PATH=$HOME/miniconda3/envs/__hicexplorer@3.7.6
#
# Nothing is ever installed into that environment; it is read only.
#
# Provided targets:
#   hicx::hdf5   HDF5 C library
#   hicx::blosc  c-blosc 1.x, needed to decode PyTables/blosc compressed
#                datasets in the HiCExplorer h5 format
#   hicx::zlib   zlib
#   hicx::lapack BLAS and LAPACK. This is the *same* libopenblasp-r0.3.28.so
#                that numpy and scipy in that prefix are built against, which
#                is what makes the dgeev column order and eigenvector signs of
#                hicPCA compatibility mode reproducible (cpp/PLAN.md 3.4, 5.4).
#   hicx::hts    htslib, the BAM reader and writer behind hicBuildMatrix

include_guard(GLOBAL)

# Derive the dependency root from CMAKE_PREFIX_PATH, so that a single
# -DCMAKE_PREFIX_PATH=<conda env> is enough to configure the project.
if(NOT DEFINED HICX_DEPS_ROOT)
    foreach(prefix IN LISTS CMAKE_PREFIX_PATH)
        if(EXISTS "${prefix}/include/hdf5.h")
            set(HICX_DEPS_ROOT "${prefix}")
            break()
        endif()
    endforeach()
endif()
if(NOT DEFINED HICX_DEPS_ROOT AND DEFINED ENV{CONDA_PREFIX})
    if(EXISTS "$ENV{CONDA_PREFIX}/include/hdf5.h")
        set(HICX_DEPS_ROOT "$ENV{CONDA_PREFIX}")
    endif()
endif()
set(HICX_DEPS_ROOT "${HICX_DEPS_ROOT}" CACHE PATH "Prefix holding HDF5, blosc and zlib")

if(NOT HICX_DEPS_ROOT)
    message(FATAL_ERROR
        "No dependency prefix found. Configure with\n"
        "  -DCMAKE_PREFIX_PATH=<conda env with hdf5, blosc, zlib>")
endif()

find_path(HICX_HDF5_INCLUDE_DIR
    NAMES hdf5.h
    HINTS "${HICX_DEPS_ROOT}/include"
    NO_DEFAULT_PATH)
find_library(HICX_HDF5_LIBRARY
    NAMES hdf5
    HINTS "${HICX_DEPS_ROOT}/lib"
    NO_DEFAULT_PATH)

find_path(HICX_BLOSC_INCLUDE_DIR
    NAMES blosc.h
    HINTS "${HICX_DEPS_ROOT}/include"
    NO_DEFAULT_PATH)
find_library(HICX_BLOSC_LIBRARY
    NAMES blosc
    HINTS "${HICX_DEPS_ROOT}/lib"
    NO_DEFAULT_PATH)

find_path(HICX_ZLIB_INCLUDE_DIR
    NAMES zlib.h
    HINTS "${HICX_DEPS_ROOT}/include"
    NO_DEFAULT_PATH)
find_library(HICX_ZLIB_LIBRARY
    NAMES z
    HINTS "${HICX_DEPS_ROOT}/lib"
    NO_DEFAULT_PATH)

find_path(HICX_HTS_INCLUDE_DIR
    NAMES htslib/sam.h
    HINTS "${HICX_DEPS_ROOT}/include"
    NO_DEFAULT_PATH)
find_library(HICX_HTS_LIBRARY
    NAMES hts
    HINTS "${HICX_DEPS_ROOT}/lib"
    NO_DEFAULT_PATH)

foreach(var HICX_HDF5_INCLUDE_DIR HICX_HDF5_LIBRARY
            HICX_BLOSC_INCLUDE_DIR HICX_BLOSC_LIBRARY
            HICX_ZLIB_INCLUDE_DIR HICX_ZLIB_LIBRARY
            HICX_HTS_INCLUDE_DIR HICX_HTS_LIBRARY)
    if(NOT ${var})
        message(FATAL_ERROR "${var} not found under ${HICX_DEPS_ROOT}")
    endif()
endforeach()

add_library(hicx_hdf5 INTERFACE)
target_include_directories(hicx_hdf5 SYSTEM INTERFACE "${HICX_HDF5_INCLUDE_DIR}")
target_link_libraries(hicx_hdf5 INTERFACE "${HICX_HDF5_LIBRARY}")
add_library(hicx::hdf5 ALIAS hicx_hdf5)

add_library(hicx_blosc INTERFACE)
target_include_directories(hicx_blosc SYSTEM INTERFACE "${HICX_BLOSC_INCLUDE_DIR}")
target_link_libraries(hicx_blosc INTERFACE "${HICX_BLOSC_LIBRARY}")
add_library(hicx::blosc ALIAS hicx_blosc)

add_library(hicx_zlib INTERFACE)
target_include_directories(hicx_zlib SYSTEM INTERFACE "${HICX_ZLIB_INCLUDE_DIR}")
target_link_libraries(hicx_zlib INTERFACE "${HICX_ZLIB_LIBRARY}")
add_library(hicx::zlib ALIAS hicx_zlib)

# BLAS and LAPACK. In this prefix libblas.so and liblapack.so are both symlinks
# to libopenblasp-r0.3.28.so, the build scipy links against, so one library
# provides dgeev, dsyevr, dsyrk and openblas_set_num_threads. There is no
# header: the Fortran symbols are declared in core/include/hicx/lapack_shim.hpp.
find_library(HICX_LAPACK_LIBRARY
    NAMES openblas lapack
    HINTS "${HICX_DEPS_ROOT}/lib"
    NO_DEFAULT_PATH)
if(NOT HICX_LAPACK_LIBRARY)
    message(FATAL_ERROR "No BLAS/LAPACK found under ${HICX_DEPS_ROOT}/lib")
endif()

add_library(hicx_lapack INTERFACE)
target_link_libraries(hicx_lapack INTERFACE "${HICX_LAPACK_LIBRARY}")
add_library(hicx::lapack ALIAS hicx_lapack)

add_library(hicx_hts INTERFACE)
target_include_directories(hicx_hts SYSTEM INTERFACE "${HICX_HTS_INCLUDE_DIR}")
target_link_libraries(hicx_hts INTERFACE "${HICX_HTS_LIBRARY}")
add_library(hicx::hts ALIAS hicx_hts)

# libBigWig, the C library pyBigWig wraps. hicPCA reads a histone mark track
# and writes its eigenvectors in bigWig, and cpp/PLAN.md 3.5 decides to vendor
# the library rather than reimplement the zoom levels and the R-tree index.
#
# Its own CMakeLists declares cmake_minimum_required(VERSION 3.8), which CMake 4
# rejects, so only the sources are fetched (SOURCE_SUBDIR does-not-exist, the
# same trick tests/CMakeLists.txt uses for doctest) and the five translation
# units are compiled here. NOCURL removes the libcurl dependency; the tools
# only ever open local files.
include(FetchContent)
FetchContent_Declare(libbigwig
    GIT_REPOSITORY https://github.com/dpryan79/libBigWig.git
    GIT_TAG 43c294ef1721a73b760803ca5e9410d581b98f17
    GIT_SHALLOW FALSE
    SOURCE_SUBDIR does-not-exist)
FetchContent_MakeAvailable(libbigwig)

add_library(hicx_bigwig STATIC
    "${libbigwig_SOURCE_DIR}/bwRead.c"
    "${libbigwig_SOURCE_DIR}/bwStats.c"
    "${libbigwig_SOURCE_DIR}/bwValues.c"
    "${libbigwig_SOURCE_DIR}/bwWrite.c"
    "${libbigwig_SOURCE_DIR}/io.c")
target_include_directories(hicx_bigwig SYSTEM PUBLIC "${libbigwig_SOURCE_DIR}"
                                                     "${HICX_ZLIB_INCLUDE_DIR}")
target_compile_definitions(hicx_bigwig PUBLIC NOCURL)
target_link_libraries(hicx_bigwig PUBLIC "${HICX_ZLIB_LIBRARY}" m)
set_target_properties(hicx_bigwig PROPERTIES POSITION_INDEPENDENT_CODE ON C_STANDARD 11)
add_library(hicx::bigwig ALIAS hicx_bigwig)

# The conda prefix is not on the loader path, so bake it into the binaries.
get_filename_component(_hicx_libdir "${HICX_HDF5_LIBRARY}" DIRECTORY)
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)
set(CMAKE_BUILD_RPATH "${_hicx_libdir}")
set(CMAKE_INSTALL_RPATH "${_hicx_libdir}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
