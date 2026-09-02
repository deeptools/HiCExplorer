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

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(HicxDependencies)

add_library(hicx_compile_options INTERFACE)
target_compile_options(hicx_compile_options INTERFACE
    -Wall -Wextra -Wpedantic
    $<$<BOOL:${HICX_WARNINGS_AS_ERRORS}>:-Werror>)

add_subdirectory(core)
add_subdirectory(tools)

if(HICX_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

message(STATUS "hicexplorer4: build type      ${CMAKE_BUILD_TYPE}")
message(STATUS "hicexplorer4: dependency root ${HICX_DEPS_ROOT}")
