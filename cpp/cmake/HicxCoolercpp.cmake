# coolercpp: the cooler-compatible library every cool and mcool read and write
# of HiCExplorer v4 goes through (core/src/cool_adapter.cpp keeps only what
# hicmatrix adds on top of cooler).
#
# How coolercpp is obtained, in this order:
#
#  1. An installed coolercpp, found through its CMake package:
#       find_package(coolercpp 0.2 CONFIG)
#     Point CMAKE_PREFIX_PATH (or coolercpp_DIR) at the install prefix. The
#     package must be built from the commit pinned below; coolercpp raises its
#     minor version whenever its API grows, so 0.2 is the API this tree needs.
#
#  2. Otherwise FetchContent, from a git repository given at configure time
#     and checked out at the pinned commit, so that every build of this
#     revision compiles the same coolercpp:
#       -DHICX_COOLERCPP_GIT_REPOSITORY=<URL or local path of a coolercpp clone>
#     The repository location is a cache variable because it is a property of
#     the machine, not of this source tree; the commit is not.
#
#  -DHICX_COOLERCPP_FORCE_FETCH=ON skips step 1.
#
# Updating coolercpp means changing HICX_COOLERCPP_COMMIT (and the version in
# step 1 when the API changed) in a reviewed commit of this tree.

include_guard(GLOBAL)

set(HICX_COOLERCPP_COMMIT a4683d4c0dc0ea90338332a622540d40242c8429)
set(HICX_COOLERCPP_VERSION 0.2)

set(HICX_COOLERCPP_GIT_REPOSITORY "" CACHE STRING
    "Git repository (URL or local path) coolercpp is fetched from when no installed coolercpp is found")
option(HICX_COOLERCPP_FORCE_FETCH "Fetch coolercpp even if an installed package exists" OFF)

if(NOT HICX_COOLERCPP_FORCE_FETCH)
    find_package(coolercpp ${HICX_COOLERCPP_VERSION} CONFIG QUIET)
endif()

if(coolercpp_FOUND AND NOT HICX_COOLERCPP_FORCE_FETCH)
    message(STATUS "hicexplorer4: coolercpp ${coolercpp_VERSION} from ${coolercpp_DIR}")
else()
    if(NOT HICX_COOLERCPP_GIT_REPOSITORY)
        message(FATAL_ERROR
            "coolercpp ${HICX_COOLERCPP_VERSION} was not found as an installed package.\n"
            "Either install coolercpp at commit ${HICX_COOLERCPP_COMMIT} and add its prefix to "
            "CMAKE_PREFIX_PATH, or let the build fetch it:\n"
            "  -DHICX_COOLERCPP_GIT_REPOSITORY=<URL or path of a coolercpp git repository>")
    endif()
    include(FetchContent)
    FetchContent_Declare(coolercpp
        GIT_REPOSITORY "${HICX_COOLERCPP_GIT_REPOSITORY}"
        GIT_TAG ${HICX_COOLERCPP_COMMIT})
    FetchContent_MakeAvailable(coolercpp)
    message(STATUS "hicexplorer4: coolercpp fetched at ${HICX_COOLERCPP_COMMIT}")
endif()

if(NOT TARGET coolercpp::coolercpp)
    message(FATAL_ERROR "coolercpp did not provide the coolercpp::coolercpp target")
endif()
