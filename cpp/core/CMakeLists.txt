# libhicx4core: matrix model, file formats and the numpy/scipy compatible
# numerics that the textual output of the tools depends on.

# Take the version from the Python package so that both implementations report
# the same string.
set(_version_file "${CMAKE_CURRENT_SOURCE_DIR}/../../hicexplorer/_version.py")
set(HICX_PYTHON_VERSION "unknown")
if(EXISTS "${_version_file}")
    file(READ "${_version_file}" _version_content)
    string(REGEX MATCH "__version__[ ]*=[ ]*'([^']+)'" _match "${_version_content}")
    if(CMAKE_MATCH_1)
        set(HICX_PYTHON_VERSION "${CMAKE_MATCH_1}")
    endif()
endif()
message(STATUS "hicexplorer4: reference version ${HICX_PYTHON_VERSION}")
configure_file(include/hicx/version.hpp.in
               "${CMAKE_CURRENT_BINARY_DIR}/generated/hicx/version.hpp" @ONLY)

add_library(hicx4core
    src/adjust_ops.cpp
    src/bam_file.cpp
    src/bedtools_ops.cpp
    src/bins.cpp
    src/clustering.cpp
    src/correct_ops.cpp
    src/build_matrix.cpp
    src/chic_hdf5.cpp
    src/chic_viewpoint.cpp
    src/cool_file.cpp
    src/fasta_reader.cpp
    src/h5_file.cpp
    src/hdf5_util.cpp
    src/hic_matrix.cpp
    src/json_lite.cpp
    src/lbfgsb.cpp
    src/lbfgsb_scipy.cpp
    src/math/ice.cpp
    src/math/kr_balancing.cpp
    src/math/sparse_kernels.cpp
    src/matrix_ops.cpp
    src/npz_file.cpp
    src/numpy_compat.cpp
    src/numpy_random.cpp
    src/numpy_sort.cpp
    src/numpy_sort_avx512.cpp
    src/obsexp_ops.cpp
    src/reduce_matrix.cpp
    src/resource_usage.cpp
    src/simd_reduce.cpp
    src/scipy_special.cpp
    src/sparse_matrix.cpp
    src/stats_ops.cpp
    src/text_formats.cpp
    src/text_table.cpp
    src/tool_matrix.cpp
    src/transform_ops.cpp)

target_include_directories(hicx4core PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_BINARY_DIR}/generated")

target_link_libraries(hicx4core
    PUBLIC hicx::hdf5
    PRIVATE hicx::blosc hicx::zlib hicx::hts hicx::lapack coolercpp::coolercpp
            hicx_compile_options)

# np.argsort's tie order on a CPU with AVX512_SKX comes from x86-simd-sort, at
# the commit numpy 1.26.4 vendors as a submodule (numpy_sort.hpp). The header is
# compiled into one translation unit with numpy's SKX target flags and is only
# called after a runtime CPU check; everything else stays at the baseline.
include(FetchContent)
FetchContent_Declare(x86simdsort
    GIT_REPOSITORY https://github.com/intel/x86-simd-sort.git
    GIT_TAG 0631a88763a4a0a4c9e84d5eeb0ec5d36053730b
    GIT_SHALLOW FALSE
    SOURCE_SUBDIR does-not-exist)
FetchContent_MakeAvailable(x86simdsort)
set_source_files_properties(src/numpy_sort_avx512.cpp PROPERTIES
    COMPILE_OPTIONS "-mavx512f;-mavx512cd;-mavx512bw;-mavx512dq;-mavx512vl;-isystem${x86simdsort_SOURCE_DIR}/src")

# clustering.cpp loads OpenBLAS with dlopen on first use instead of linking it,
# so that a tool which may cluster does not start the BLAS threads at load.
get_filename_component(_hicx_openblas_real "${HICX_LAPACK_LIBRARY}" REALPATH)
set_source_files_properties(src/clustering.cpp PROPERTIES
    COMPILE_DEFINITIONS "HICX_OPENBLAS_LIBRARY=\"${_hicx_openblas_real}\"")
target_link_libraries(hicx4core PRIVATE ${CMAKE_DL_LIBS})

# transform_ops runs the covariance rows on std::thread workers.
find_package(Threads REQUIRED)
target_link_libraries(hicx4core PRIVATE Threads::Threads)

# Boost.Math, header only, for hicx::scipy::betainc (hicx/scipy_special.hpp).
# scipy.special.betainc is Boost's ibeta since scipy 1.12, so the exact
# commit scipy v1.14.1 vendors as scipy/_lib/boost_math is used, fetched as a
# tarball rather than a clone (the repository history is large), and in
# standalone mode so that no other part of Boost is needed.
include(FetchContent)
FetchContent_Declare(boost_math
    URL https://codeload.github.com/boostorg/math/tar.gz/a53b013c735caa98179532a32ad24d34569b9710
    URL_HASH SHA256=8fe7ea31430e0a8652acf2659f78b2235e4947f0dca7e50a22257908294b35ed
    SOURCE_SUBDIR does-not-exist)
FetchContent_MakeAvailable(boost_math)
add_library(hicx_boost_math INTERFACE)
target_include_directories(hicx_boost_math SYSTEM INTERFACE "${boost_math_SOURCE_DIR}/include")
target_compile_definitions(hicx_boost_math INTERFACE BOOST_MATH_STANDALONE=1)
target_link_libraries(hicx4core PRIVATE hicx_boost_math)

add_library(hicx::core ALIAS hicx4core)
