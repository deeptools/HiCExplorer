# The Python extensions (pybind11):
#
#   hicexplorer._hicx  the C++ implementation of the public Python API
#                      (src/api, cpp/python/API_BINDINGS.md), shipped in
#                      the wheel next to the hicexplorer package
#   hicx_matrix        reads one region of a cool, mcool, HiCExplorer h5 or
#                      .hic matrix (cpp/PLAN.md tier 10, item 10.2)
#
#   cmake -S cpp -B <build> -DHICX_BUILD_PYTHON=ON -DPython_EXECUTABLE=<python3>
#
# hicx_matrix lands in <build>/python, _hicx in <build>/python/hicexplorer,
# which a source checkout's hicexplorer package finds on its own
# (hicexplorer/__init__.py). pybind11 comes from the dependency prefix.

find_package(Python 3.10 COMPONENTS Interpreter Development.Module REQUIRED)

# pybind11 2.13.6 from the dependency prefix (conda package pybind11), the
# release the bindings were verified with.
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 2.13.6 EXACT CONFIG REQUIRED)
message(STATUS "hicexplorer4: pybind11 ${pybind11_VERSION} for Python ${Python_VERSION}")

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
    hicx_compile_options)

set_target_properties(hicx_matrix PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

# ---------------------------------------------------------------------------
# hicexplorer._hicx

include("${CMAKE_CURRENT_SOURCE_DIR}/src/api/modules.cmake")

set(_hicx_registry "")
set(_hicx_module_sources "")
foreach(_module IN LISTS HICX_API_MODULES)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/api/${_module}.cpp")
        message(FATAL_ERROR "src/api/modules.cmake lists ${_module}, but src/api/${_module}.cpp does not exist")
    endif()
    string(APPEND _hicx_registry "HICX_API_MODULE(${_module})\n")
    list(APPEND _hicx_module_sources "src/api/${_module}.cpp")
endforeach()
set(_hicx_generated "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(CONFIGURE OUTPUT "${_hicx_generated}/api_registry.inc"
     CONTENT "// Generated from src/api/modules.cmake.\n${_hicx_registry}")

pybind11_add_module(_hicx MODULE NO_EXTRAS
    src/api/module.cpp
    src/api/support/convert.cpp
    src/api/support/errors.cpp
    src/api/support/runtime.cpp
    ${_hicx_module_sources})

target_include_directories(_hicx PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src/api"
    "${_hicx_generated}")

target_link_libraries(_hicx PRIVATE
    hicx::core
    hicx_compile_options)

set_target_properties(_hicx PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/hicexplorer")
