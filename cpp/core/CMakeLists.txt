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
    src/argparse.cpp
    src/bam_file.cpp
    src/bedtools_ops.cpp
    src/bins.cpp
    src/clustering.cpp
    src/correct_ops.cpp
    src/build_matrix.cpp
    src/chic_hdf5.cpp
    src/chic_viewpoint.cpp
    src/cool_adapter.cpp
    src/fasta_reader.cpp
    src/h5_file.cpp
    src/hdf5_util.cpp
    src/hic_adapter.cpp
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
    src/plot_bridge.cpp
    src/reduce_matrix.cpp
    src/resource_usage.cpp
    src/simd_reduce.cpp
    src/scipy_special.cpp
    src/scipy_stats.cpp
    src/tar_gz.cpp
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
    PRIVATE hicx::blosc hicx::zlib hicx::hts hicx::lapack coolercpp::coolercpp hicfilecpp::hicfilecpp
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

# scipy's C++ translation of Cephes, for hicx/scipy_stats.hpp. In scipy 1.14.1
# special.gammaincinv and special.chdtrc, under stats.chi2.ppf and the p-value
# of stats.chi2_contingency, map to cephes_igami and cephes_chdtrc
# (scipy/special/functions.json), which are header only in
# scipy/special/special/cephes. The headers those two need are downloaded one
# by one at the commit of the v1.14.1 tag and checked against their SHA256.
set(HICX_SCIPY_COMMIT 92d2a8592782ee19a1161d0bf3fc2241ba78bb63)
set(HICX_SCIPY_SPECIAL_DIR "${CMAKE_BINARY_DIR}/_deps/scipy_special")
set(HICX_SCIPY_SPECIAL_HEADERS
    "special/config.h=1d96ebe9556cb3eb8119f5ef340f8c301e399c09810ccea2cbe21c4cf3f5d833"
    "special/error.h=0bb3e91080f4052fa263d214b2eb1e4a62fe6021009bd51a3a0aa56df01e8f4a"
    "special/cephes/chdtr.h=27a20cabeb467a83a2b2af1ee65e5ba5b7652c900b4275ac0d6d462ca595d039"
    "special/cephes/const.h=94a81da5e3f6dd3f7f3cc03a55628c509a57b021e051a0d88f7e448cdf0181b8"
    "special/cephes/gamma.h=e7f1c1b70a58e480bfdb4f038cb61a4ff2d026a7e15220f1a35a921eb29dd619"
    "special/cephes/igam_asymp_coeff.h=21b23621cd3feaa909da3cf1f36de34fc1cfd2cf92368ee3e83eacac06eb9e70"
    "special/cephes/igam.h=ea06624ca774c21a50875e339e51b6942405b1a32162e19086155291a17c1548"
    "special/cephes/igami.h=68757c6422982eb08acf017734fadf9a6d6c951ba71ddb2f907e7eb4e891a4bc"
    "special/cephes/lanczos.h=203895d8c085e89de24158f5ca8b63835d38922a8e33843fad8f22aff44efc35"
    "special/cephes/ndtr.h=3c5b035184d383facd5885266ed1104f5252b1c284381334ab7472c6e3d42158"
    "special/cephes/polevl.h=a333233a8ddd758f7dda84f781187260a4ada3132ae55d047e3e06ba026e3059"
    "special/cephes/trig.h=30e7d4d067e88912832d58d23b0e8384c5b310826ae9bca9347e3ac8a6a4a714"
    "special/cephes/unity.h=aec1b0010fd0dd17508796077e395350606b7b3e7bddad27613147d6891084ac"
    "special/cephes/zeta.h=6f6b8625171b34602a72e981a2573ad34c5083782fc63fbc8bc1e0d9abb348ac")
foreach(entry IN LISTS HICX_SCIPY_SPECIAL_HEADERS)
    string(REPLACE "=" ";" parts "${entry}")
    list(GET parts 0 relative)
    list(GET parts 1 expected)
    set(header "${HICX_SCIPY_SPECIAL_DIR}/${relative}")
    set(actual "")
    if(EXISTS "${header}")
        file(SHA256 "${header}" actual)
    endif()
    if(NOT actual STREQUAL expected)
        file(DOWNLOAD
            "https://raw.githubusercontent.com/scipy/scipy/${HICX_SCIPY_COMMIT}/scipy/special/${relative}"
            "${header}" EXPECTED_HASH SHA256=${expected} STATUS status)
        list(GET status 0 code)
        if(NOT code EQUAL 0)
            message(FATAL_ERROR "cannot download scipy/special/${relative}: ${status}")
        endif()
    endif()
endforeach()
add_library(hicx_scipy_special INTERFACE)
target_include_directories(hicx_scipy_special SYSTEM INTERFACE "${HICX_SCIPY_SPECIAL_DIR}")
target_link_libraries(hicx4core PRIVATE hicx_scipy_special)

add_library(hicx::core ALIAS hicx4core)
