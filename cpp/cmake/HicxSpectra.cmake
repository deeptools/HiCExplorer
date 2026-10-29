# Spectra, the header-only Krylov eigensolver behind hicPCA --eigenSolver
# lanczos. MPL-2.0. The headers of release v1.2.0 (commit 6841bcb) are shipped
# unchanged in third_party/spectra, because conda-forge has no spectra-cpp
# package for linux-aarch64.
#
# Spectra is written on Eigen. The dependency prefix ships Eigen 3.4.0 under
# include/eigen3, which is what this module points at.
#
# Provided target:
#   hicx::spectra   include directories for Spectra and Eigen

include_guard(GLOBAL)

find_path(HICX_SPECTRA_INCLUDE_DIR
    NAMES Spectra/Util/Version.h
    HINTS "${CMAKE_CURRENT_LIST_DIR}/../third_party/spectra/include"
    NO_DEFAULT_PATH)
if(NOT HICX_SPECTRA_INCLUDE_DIR)
    message(FATAL_ERROR "Spectra headers missing from cpp/third_party/spectra/include")
endif()
file(STRINGS "${HICX_SPECTRA_INCLUDE_DIR}/Spectra/Util/Version.h" _hicx_spectra_version
     REGEX "^#define SPECTRA_(MAJOR|MINOR|PATCH)_VERSION [0-9]+$")
string(REGEX REPLACE "[^0-9;]" "" _hicx_spectra_digits "${_hicx_spectra_version}")
if(NOT _hicx_spectra_digits STREQUAL "1;2;0")
    message(FATAL_ERROR "Spectra 1.2.0 is required, found: ${_hicx_spectra_version}")
endif()

find_path(HICX_EIGEN3_INCLUDE_DIR
    NAMES Eigen/Core
    HINTS "${HICX_DEPS_ROOT}/include/eigen3"
    NO_DEFAULT_PATH)
if(NOT HICX_EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Eigen 3 not found under ${HICX_DEPS_ROOT}/include/eigen3")
endif()

add_library(hicx_spectra INTERFACE)
target_include_directories(hicx_spectra SYSTEM INTERFACE
    "${HICX_SPECTRA_INCLUDE_DIR}"
    "${HICX_EIGEN3_INCLUDE_DIR}")
add_library(hicx::spectra ALIAS hicx_spectra)
