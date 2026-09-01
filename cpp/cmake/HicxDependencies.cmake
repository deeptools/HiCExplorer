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

foreach(var HICX_HDF5_INCLUDE_DIR HICX_HDF5_LIBRARY
            HICX_BLOSC_INCLUDE_DIR HICX_BLOSC_LIBRARY
            HICX_ZLIB_INCLUDE_DIR HICX_ZLIB_LIBRARY)
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

# The conda prefix is not on the loader path, so bake it into the binaries.
get_filename_component(_hicx_libdir "${HICX_HDF5_LIBRARY}" DIRECTORY)
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)
set(CMAKE_BUILD_RPATH "${_hicx_libdir}")
set(CMAKE_INSTALL_RPATH "${_hicx_libdir}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
