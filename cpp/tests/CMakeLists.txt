# Unit tests. doctest is header only; its own CMake project still declares a
# pre 3.5 minimum, which CMake 4 rejects, so only the sources are fetched and
# the header is exposed through an interface target.

include(FetchContent)

FetchContent_Declare(doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.11
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR does-not-exist)
FetchContent_MakeAvailable(doctest)

add_library(doctest_header INTERFACE)
target_include_directories(doctest_header SYSTEM INTERFACE "${doctest_SOURCE_DIR}")

add_executable(hicx_tests
    test_main.cpp
    test_numpy_compat.cpp
    test_bins.cpp
    test_sparse_matrix.cpp
    test_matrix_ops.cpp
    test_json_lite.cpp
    test_readers.cpp
    test_reduce_matrix.cpp
    test_text_formats.cpp
    test_writers.cpp)

target_link_libraries(hicx_tests PRIVATE hicx::core doctest_header hicx_compile_options)

# Real matrices from the Python test suite, used by the reader tests.
target_compile_definitions(hicx_tests PRIVATE
    HICX_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../../hicexplorer/test/test_data")

add_test(NAME hicx_unit_tests COMMAND hicx_tests)

# Not a HiCExplorer tool: the smallest driver that exercises the tier 0 file
# layer end to end, used by cpp/scripts/equiv.py against its Python twin
# cpp/scripts/py_roundtrip.py.
add_executable(hicx_roundtrip hicx_roundtrip.cpp)
target_link_libraries(hicx_roundtrip PRIVATE hicx::core hicx_compile_options)
