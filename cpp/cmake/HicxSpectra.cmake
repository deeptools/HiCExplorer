# Spectra, the header-only Krylov eigensolver behind hicPCA --eigenSolver
# lanczos (cpp/PLAN.md tier 11). MPL-2.0. Pinned to the commit of release
# v1.2.0. Its own CMake project is not used: only the headers are fetched,
# the same way doctest and x86-simd-sort are.
#
# Spectra is written on Eigen. The dependency prefix ships Eigen 3.4.0 under
# include/eigen3, which is what this module points at.
#
# Provided target:
#   hicx::spectra   include directories for Spectra and Eigen

include_guard(GLOBAL)
include(FetchContent)

FetchContent_Declare(spectra
    GIT_REPOSITORY https://github.com/yixuan/spectra.git
    GIT_TAG 6841bcbacaa0f0a8446210314e682057a084be4e
    GIT_SHALLOW FALSE
    SOURCE_SUBDIR does-not-exist)
FetchContent_MakeAvailable(spectra)

find_path(HICX_EIGEN3_INCLUDE_DIR
    NAMES Eigen/Core
    HINTS "${HICX_DEPS_ROOT}/include/eigen3"
    NO_DEFAULT_PATH)
if(NOT HICX_EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Eigen 3 not found under ${HICX_DEPS_ROOT}/include/eigen3")
endif()

add_library(hicx_spectra INTERFACE)
target_include_directories(hicx_spectra SYSTEM INTERFACE
    "${spectra_SOURCE_DIR}/include"
    "${HICX_EIGEN3_INCLUDE_DIR}")
add_library(hicx::spectra ALIAS hicx_spectra)
