# hicfilecpp: the Juicer .hic reader and writer every .hic read and write of
# HiCExplorer v4 goes through (core/src/hic_adapter.cpp keeps what hic2cool
# and hicConvertFormat add on top of the format).
#
# How hicfilecpp is obtained, in this order:
#
#  1. An installed hicfilecpp, found through its CMake package:
#       find_package(hicfilecpp 0.3 CONFIG)
#     Point CMAKE_PREFIX_PATH (or hicfilecpp_DIR) at the install prefix. The
#     package must be built from the commit pinned below; hicfilecpp raises its
#     minor version whenever its API grows, so 0.3 is the API this tree needs.
#
#  2. Otherwise FetchContent, from a git repository given at configure time
#     and checked out at the pinned commit, so that every build of this
#     revision compiles the same hicfilecpp:
#       -DHICX_HICFILECPP_GIT_REPOSITORY=<URL or local path of a hicfilecpp clone>
#     The repository location is a cache variable because it is a property of
#     the machine, not of this source tree; the commit is not.
#
#  -DHICX_HICFILECPP_FORCE_FETCH=ON skips step 1.
#
# Updating hicfilecpp means changing HICX_HICFILECPP_COMMIT (and the version in
# step 1 when the API changed) in a reviewed commit of this tree.

include_guard(GLOBAL)

set(HICX_HICFILECPP_COMMIT 99615bd06e949b7a9d26398965b8247de93aba07)
set(HICX_HICFILECPP_VERSION 0.3)

set(HICX_HICFILECPP_GIT_REPOSITORY "" CACHE STRING
    "Git repository (URL or local path) hicfilecpp is fetched from when no installed hicfilecpp is found")
option(HICX_HICFILECPP_FORCE_FETCH "Fetch hicfilecpp even if an installed package exists" OFF)

if(NOT HICX_HICFILECPP_FORCE_FETCH)
    find_package(hicfilecpp ${HICX_HICFILECPP_VERSION} CONFIG QUIET)
endif()

if(hicfilecpp_FOUND AND NOT HICX_HICFILECPP_FORCE_FETCH)
    message(STATUS "hicexplorer4: hicfilecpp ${hicfilecpp_VERSION} from ${hicfilecpp_DIR}")
else()
    if(NOT HICX_HICFILECPP_GIT_REPOSITORY)
        message(FATAL_ERROR
            "hicfilecpp ${HICX_HICFILECPP_VERSION} was not found as an installed package.\n"
            "Either install hicfilecpp at commit ${HICX_HICFILECPP_COMMIT} and add its prefix to "
            "CMAKE_PREFIX_PATH, or let the build fetch it:\n"
            "  -DHICX_HICFILECPP_GIT_REPOSITORY=<URL or path of a hicfilecpp git repository>")
    endif()
    include(FetchContent)
    FetchContent_Declare(hicfilecpp
        GIT_REPOSITORY "${HICX_HICFILECPP_GIT_REPOSITORY}"
        GIT_TAG ${HICX_HICFILECPP_COMMIT})
    FetchContent_MakeAvailable(hicfilecpp)
    message(STATUS "hicexplorer4: hicfilecpp fetched at ${HICX_HICFILECPP_COMMIT}")
endif()

if(NOT TARGET hicfilecpp::hicfilecpp)
    message(FATAL_ERROR "hicfilecpp did not provide the hicfilecpp::hicfilecpp target")
endif()
