# hicx_matrix: Python bindings that read one region of a cool, mcool,
# HiCExplorer h5 or .hic matrix (cpp/PLAN.md tier 10, item 10.2).
#
#   cmake -S cpp -B <build> -DHICX_BUILD_PYTHON=ON -DPython_EXECUTABLE=<python3.12>
#   cmake --build <build> --target hicx_matrix
#
# The module lands in <build>/python; put that directory on PYTHONPATH.
# pybind11 is fetched at a pinned release tag and never taken from the
# environment, so every build of this revision compiles the same bindings.

find_package(Python 3.12 COMPONENTS Interpreter Development REQUIRED)

include(FetchContent)
set(HICX_PYBIND11_TAG v2.13.6)
# pybind11 then uses the Python found above instead of its own search.
set(PYBIND11_FINDPYTHON ON)
FetchContent_Declare(pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG ${HICX_PYBIND11_TAG}
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(pybind11)
message(STATUS "hicexplorer4: pybind11 ${HICX_PYBIND11_TAG} for Python ${Python_VERSION}")

# A shared module can only take position independent objects, and libhicx4 is
# a static library.
set_target_properties(hicx4core PROPERTIES POSITION_INDEPENDENT_CODE ON)

# NO_EXTRAS: no link time optimisation and no stripping, so the module is built
# with the same flags as the rest of the tree.
pybind11_add_module(hicx_matrix MODULE NO_EXTRAS
    src/h5_rows.cpp
    src/hicx_matrix.cpp)

target_link_libraries(hicx_matrix PRIVATE
    hicx::core
    hicx::hdf5
    coolercpp::coolercpp
    hicfilecpp::hicfilecpp
    hicx_compile_options)

set_target_properties(hicx_matrix PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
