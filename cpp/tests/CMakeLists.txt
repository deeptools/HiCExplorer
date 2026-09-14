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
    test_adjust_ops.cpp
    test_stats_ops.cpp
    test_text_formats.cpp
    test_tier2_support.cpp
    test_transform_ops.cpp
    test_writers.cpp
    test_build_matrix.cpp
    test_correct_math.cpp
    test_detect_loops.cpp
    test_tad_contacts.cpp
    test_compartmentalization.cpp
    test_quick_qc.cpp
    test_chic.cpp
    test_merge_domains.cpp
    test_cool_chromosome_cut.cpp
    test_aggregate_support.cpp
    # The computational cores of hicDetectLoops, of the two TAD statistics tools,
    # of hicCompartmentalization and of hicMergeDomains, so that the tests link
    # against exactly the translation units the tools run.
    ../tools/detect_loops_impl.cpp
    ../tools/tad_contacts_impl.cpp
    ../tools/compartmentalization_impl.cpp
    ../tools/merge_domains_impl.cpp)

target_link_libraries(hicx_tests PRIVATE hicx::core doctest_header hicx_compile_options)

# Real matrices from the Python test suite, used by the reader tests.
target_compile_definitions(hicx_tests PRIVATE
    HICX_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../../hicexplorer/test/test_data"
    # Files the tests write go into the build tree, never into the source tree:
    # an out-of-tree build from a fresh checkout has no cpp/build directory.
    HICX_TEST_OUTPUT_DIR="${CMAKE_CURRENT_BINARY_DIR}")

add_test(NAME hicx_unit_tests COMMAND hicx_tests)

# hicCompartmentalization draws a matplotlib figure the C++ port cannot draw
# yet. A run that asks for it has to fail visibly and write nothing, while
# --noPlot has to write the numeric outputs. The equivalence harness cannot
# express a case whose exit status differs between Python and C++, so this is
# checked here.
add_test(NAME hicCompartmentalization_refuses_plot
         COMMAND ${CMAKE_COMMAND}
                 -DTOOL=$<TARGET_FILE:hicCompartmentalization>
                 -DDATA=${CMAKE_CURRENT_SOURCE_DIR}/../../hicexplorer/test/test_data
                 -DOUT=${CMAKE_CURRENT_BINARY_DIR}/hicCompartmentalization_refuses_plot
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/check_compartmentalization_refuses_plot.cmake)

# Command line behaviour that only the finished executables can show: refusals
# that must leave no output behind, and hicMergeDomains' graphviz rendering.
find_program(HICX_BASH bash REQUIRED)
set(HICX_TEST_DATA "${CMAKE_CURRENT_SOURCE_DIR}/../../hicexplorer/test/test_data")
add_test(NAME hicPlotSVL_cli
         COMMAND ${HICX_BASH} "${CMAKE_CURRENT_SOURCE_DIR}/cli/hicPlotSVL_cli.sh"
                 $<TARGET_FILE:hicPlotSVL> "${HICX_TEST_DATA}")
add_test(NAME hicMergeDomains_cli
         COMMAND ${HICX_BASH} "${CMAKE_CURRENT_SOURCE_DIR}/cli/hicMergeDomains_cli.sh"
                 $<TARGET_FILE:hicMergeDomains> "${HICX_TEST_DATA}" "${HICX_DEPS_ROOT}")

# hicAggregateContacts draws a matplotlib figure the C++ port cannot draw yet,
# the same situation as hicCompartmentalization: a run that asks for a figure
# (a plain run, or --diagnosticHeatmapFile) has to fail visibly and write
# nothing, while --noPlot has to write the numeric outputs.
add_test(NAME hicAggregateContacts_refuses_plot
         COMMAND ${CMAKE_COMMAND}
                 -DTOOL=$<TARGET_FILE:hicAggregateContacts>
                 -DDATA=${CMAKE_CURRENT_SOURCE_DIR}/../../hicexplorer/test/test_data
                 -DOUT=${CMAKE_CURRENT_BINARY_DIR}/hicAggregateContacts_refuses_plot
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/check_aggregate_contacts_refuses_plot.cmake)

# Not a HiCExplorer tool: the smallest driver that exercises the tier 0 file
# layer end to end, used by cpp/scripts/equiv.py against its Python twin
# cpp/scripts/py_roundtrip.py.
add_executable(hicx_roundtrip hicx_roundtrip.cpp)
target_link_libraries(hicx_roundtrip PRIVATE hicx::core hicx_compile_options)

# Declares the parsers of cpp/scripts/test_argparse_layer.py with hicx::cli,
# which compares its behaviour with Python's argparse.
add_executable(hicx_argparse_driver argparse_driver.cpp)
target_link_libraries(hicx_argparse_driver PRIVATE hicx::core hicx_compile_options)
