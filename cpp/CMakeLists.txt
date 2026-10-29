cmake_minimum_required(VERSION 3.24)

project(hicexplorer4
        VERSION 4.0.0
        DESCRIPTION "HiCExplorer version 4, C++ implementation"
        # C is needed for the vendored libBigWig sources (cmake/HicxDependencies.cmake).
        LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Everything produced by the build stays inside cpp/build, which .gitignore
# already excludes. The tool binaries land in cpp/build/tools, which is what
# the equivalence harness expects as --cpp-bin.

option(HICX_BUILD_TESTS "Build the C++ unit tests" ON)
option(HICX_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
# Install the tools, the drawing package and the hicexplorer._hicx extension
# into a Python wheel layout.
option(HICX_INSTALL_PYTHON_PACKAGE "Install the executables for the Python package" OFF)
# The Python extensions (python/, pybind11): hicexplorer._hicx, which the
# hicexplorer package needs, and hicx_matrix. Off by default because they
# need a Python with its development headers; the wheel build turns it on.
option(HICX_BUILD_PYTHON "Build the Python extensions hicexplorer._hicx and hicx_matrix" OFF)
if(HICX_INSTALL_PYTHON_PACKAGE AND NOT HICX_BUILD_PYTHON)
    message(STATUS "hicexplorer4: HICX_INSTALL_PYTHON_PACKAGE needs the extension, enabling HICX_BUILD_PYTHON")
    set(HICX_BUILD_PYTHON ON CACHE BOOL "Build the Python extensions hicexplorer._hicx and hicx_matrix" FORCE)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(HicxDependencies)
# coolercpp is part of libhicx4core (core/io/coolercpp). Until hicfilecpp and
# hic2cool are vendored as well, the fetched ones link this target.
add_library(hicx_coolercpp_compat INTERFACE)
add_library(coolercpp::coolercpp ALIAS hicx_coolercpp_compat)
include(HicxHicfilecpp)
include(HicxHic2cool)

add_library(hicx_compile_options INTERFACE)
target_compile_options(hicx_compile_options INTERFACE
    -Wall -Wextra -Wpedantic
    $<$<BOOL:${HICX_WARNINGS_AS_ERRORS}>:-Werror>)

add_subdirectory(core)
add_subdirectory(tools)

if(HICX_BUILD_PYTHON)
    add_subdirectory(python)
endif()

if(HICX_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(HICX_INSTALL_PYTHON_PACKAGE)
    include(HicxInstall)
endif()

message(STATUS "hicexplorer4: build type      ${CMAKE_BUILD_TYPE}")
message(STATUS "hicexplorer4: dependency root ${HICX_DEPS_ROOT}")
