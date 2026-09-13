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
    # The computational cores of hicDetectLoops, of the two TAD statistics tools
    # and of hicCompartmentalization, so that the tests link against exactly the
    # translation units the tools run.
    ../tools/detect_loops_impl.cpp
    ../tools/tad_contacts_impl.cpp
    ../tools/compartmentalization_impl.cpp)

target_link_libraries(hicx_tests PRIVATE hicx::core doctest_header hicx_compile_options)

# Real matrices from the Python test suite, used by the reader tests.
target_compile_definitions(hicx_tests PRIVATE
    HICX_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../../hicexplorer/test/test_data")

add_test(NAME hicx_unit_tests COMMAND hicx_tests)

# hicInterIntraTAD cannot draw its matplotlib ratio plot. When the user names
# the plot file explicitly it must exit non-zero before writing anything,
# neither the plot nor the table, and leave its working directory empty.
set(_refusal_dir "${CMAKE_CURRENT_BINARY_DIR}/hicInterIntraTAD_refuses_ratio_plot")
set(_refusal_data "${CMAKE_CURRENT_SOURCE_DIR}/../../hicexplorer/test/test_data")
foreach(_option IN ITEMS --outFileNameRatioPlot -op)
    add_test(NAME hicInterIntraTAD_refuses_ratio_plot${_option}
        COMMAND ${CMAKE_COMMAND}
            "-DTOOL=$<TARGET_FILE:hicInterIntraTAD>"
            "-DWORKDIR=${_refusal_dir}${_option}"
            "-DARGS=-m;${_refusal_data}/hicInterIntraTAD/GSM2644947_Auxin2days-R1.100000_chr1_chr2.cool;--tadDomains;${_refusal_data}/hicInterIntraTAD/untreated_R1_domains_chr1_chr2.bed;${_option};${_refusal_dir}${_option}/ratio.png;-o;${_refusal_dir}${_option}/table.txt"
            "-DMUST_NOT_EXIST=${_refusal_dir}${_option}/ratio.png;${_refusal_dir}${_option}/table.txt"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/refuses_unavailable_output.cmake")
endforeach()

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

# Not a HiCExplorer tool: the smallest driver that exercises the tier 0 file
# layer end to end, used by cpp/scripts/equiv.py against its Python twin
# cpp/scripts/py_roundtrip.py.
add_executable(hicx_roundtrip hicx_roundtrip.cpp)
target_link_libraries(hicx_roundtrip PRIVATE hicx::core hicx_compile_options)
