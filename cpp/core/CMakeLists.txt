# libhicx4core: matrix model, file formats and the numpy/scipy compatible
# numerics that the textual output of the tools depends on.

# Take the version from the Python package so that both implementations report
# the same string.
set(_version_file "${CMAKE_CURRENT_SOURCE_DIR}/../../hicexplorer/_version.py")
set(HICX_PYTHON_VERSION "unknown")
if(EXISTS "${_version_file}")
    file(READ "${_version_file}" _version_content)
    string(REGEX MATCH "__version__[ ]*=[ ]*'([^']+)'" _match "${_version_content}")
    if(CMAKE_MATCH_1)
        set(HICX_PYTHON_VERSION "${CMAKE_MATCH_1}")
    endif()
endif()
message(STATUS "hicexplorer4: reference version ${HICX_PYTHON_VERSION}")
configure_file(include/hicx/version.hpp.in
               "${CMAKE_CURRENT_BINARY_DIR}/generated/hicx/version.hpp" @ONLY)

add_library(hicx4core
    src/bins.cpp
    src/cool_file.cpp
    src/h5_file.cpp
    src/hdf5_util.cpp
    src/hic_matrix.cpp
    src/json_lite.cpp
    src/numpy_compat.cpp
    src/resource_usage.cpp
    src/sparse_matrix.cpp)

target_include_directories(hicx4core PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_BINARY_DIR}/generated")

target_link_libraries(hicx4core
    PUBLIC hicx::hdf5
    PRIVATE hicx::blosc hicx::zlib hicx_compile_options)

add_library(hicx::core ALIAS hicx4core)
