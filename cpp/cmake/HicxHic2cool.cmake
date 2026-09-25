# hic2cool: the .hic to cool and cool to .hic converter, used as a library
# (hic2cool::core) by hicConvertFormat. It builds on coolercpp and hicfilecpp,
# which are already defined by HicxCoolercpp and HicxHicfilecpp when this file
# is included, so it reuses their targets.
#
# hic2cool is fetched from a git repository given at configure time and
# checked out at the pinned commit, so that every build of this revision
# compiles the same hic2cool:
#
#   -DHICX_HIC2COOL_GIT_REPOSITORY=<URL or local path of a hic2cool clone>
#
# Updating hic2cool means changing HICX_HIC2COOL_COMMIT in a reviewed commit
# of this tree.

include_guard(GLOBAL)

set(HICX_HIC2COOL_COMMIT 4c550d1d7388f793fbcd3b175d6949c54058e603)

set(HICX_HIC2COOL_GIT_REPOSITORY "" CACHE STRING
    "Git repository (URL or local path) hic2cool is fetched from")

if(NOT HICX_HIC2COOL_GIT_REPOSITORY)
    message(FATAL_ERROR
        "hic2cool was not given. Configure with\n"
        "  -DHICX_HIC2COOL_GIT_REPOSITORY=<URL or path of a hic2cool git repository>\n"
        "(the commit ${HICX_HIC2COOL_COMMIT} is checked out)")
endif()

include(FetchContent)
FetchContent_Declare(hic2cool
    GIT_REPOSITORY "${HICX_HIC2COOL_GIT_REPOSITORY}"
    GIT_TAG ${HICX_HIC2COOL_COMMIT})
FetchContent_MakeAvailable(hic2cool)
message(STATUS "hicexplorer4: hic2cool fetched at ${HICX_HIC2COOL_COMMIT}")

if(NOT TARGET hic2cool::core)
    message(FATAL_ERROR "hic2cool did not provide the hic2cool::core target")
endif()
