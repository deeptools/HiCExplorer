# One executable per HiCExplorer tool. Add new tools to the list below; each
# name must match a <name>.cpp in this directory and the Python tool it ports.

set(HICX_TOOLS
    hicAggregateContacts
    hicAdjustMatrix
    hicBuildMatrix
    hicBuildMatrixMicroC
    chicQualityControl
    chicViewpoint
    chicViewpointBackgroundModel
    chicSignificantInteractions
    chicAggregateStatistic
    chicDifferentialTest
    # The CHiCAGO pipeline (Cairns et al. 2016), exposed as its own tools
    # rather than as mode flags on the chic* tools above (cpp/PLAN.md 9.15;
    # no Python HiCExplorer counterpart for any of the three).
    chicChicagoBackgroundModel
    chicChicagoScores
    chicChicagoSignificantInteractions
    chicChicagoPlotViewpoint
    hicCompareMatrices
    hicConvertFormat
    hicFindTADs
    hicCorrectMatrix
    hicInfo
    hicMergeMatrixBins
    hicMergeTADbins
    hicAverageRegions
    hicCreateThresholdFile
    hicFindRestSite
    hicMergeLoops
    hicNormalize
    hicQuickQC
    hicValidateLocations
    hicSumMatrices
    hicTransform
    hicPlotSVL
    hicPlotDistVsCounts
    hicPlotViewpoint
    hicPlotAverageRegions
    chicPlotViewpoint
    hicCorrelate
    hicPrepareQCreport
    hicPlotMatrix)

# hicQC is hicPrepareQCreport under the name bin/hicQC gives the same main();
# only the program name in its usage and messages differs.
add_executable(hicQC hicPrepareQCreport.cpp)
target_link_libraries(hicQC PRIVATE hicx::core hicx_compile_options)
target_compile_definitions(hicQC PRIVATE HICX_QC_PROG="hicQC")

foreach(tool IN LISTS HICX_TOOLS)
    add_executable(${tool} "${tool}.cpp")
    target_link_libraries(${tool} PRIVATE hicx::core hicx_compile_options)
endforeach()


# The drawing layer of the plotting tools (cpp/PLAN.md tier 7, option (a)),
# copied next to the tools directory, where hicx::plot::draw finds it as
# <tools>/../plot without an installed package.
add_custom_target(hicx_plot_package ALL
    COMMAND "${CMAKE_COMMAND}" -E copy_directory
            "${PROJECT_SOURCE_DIR}/../plot/hicexplorer_plot"
            "${CMAKE_BINARY_DIR}/plot/hicexplorer_plot"
    COMMENT "Copying the hicexplorer_plot drawing package")

# hicDetectLoops keeps its computational core in a second translation unit so
# that cpp/tests/test_detect_loops.cpp can link against exactly the code the
# tool runs.
add_executable(hicDetectLoops hicDetectLoops.cpp detect_loops_impl.cpp)
target_link_libraries(hicDetectLoops PRIVATE hicx::core hicx_compile_options)

# hicDetectStripes has no Python counterpart (cpp/PLAN.md tier 9 section 9.3).
# Its statistical core is a second translation unit for the same reason as
# hicDetectLoops': cpp/tests/test_detect_stripes.cpp links against exactly the
# code the tool runs.
add_executable(hicDetectStripes hicDetectStripes.cpp stripes_impl.cpp)
target_link_libraries(hicDetectStripes PRIVATE hicx::core hicx_compile_options)

# hicArrowhead has no Python counterpart (cpp/PLAN.md tier 9): a faithful
# reimplementation of Juicer/Juicebox's own Arrowhead contact-domain caller.
# Its algorithmic core is a second translation unit for the same reason as
# hicDetectLoops'/hicDetectStripes': cpp/tests/test_arrowhead.cpp links
# against exactly the code the tool runs.
add_executable(hicArrowhead hicArrowhead.cpp arrowhead_impl.cpp)
target_link_libraries(hicArrowhead PRIVATE hicx::core hicx_compile_options)

# hicDifferentialTAD and hicInterIntraTAD cut their TAD blocks with the same
# code, which is also what cpp/tests/test_tad_contacts.cpp links against.
foreach(tool IN ITEMS hicDifferentialTAD hicInterIntraTAD)
    add_executable(${tool} "${tool}.cpp" tad_contacts_impl.cpp)
    target_link_libraries(${tool} PRIVATE hicx::core hicx_compile_options)
endforeach()

# hicDifferentialAnalysis has no Python counterpart (cpp/PLAN.md 9.7, work item
# 2). Its count model and its contact input are separate translation units, so
# that cpp/tests/test_diff_engine.cpp links against exactly the code it runs.
add_executable(hicDifferentialAnalysis hicDifferentialAnalysis.cpp diff_engine_impl.cpp
                                       diff_contacts_impl.cpp)
target_link_libraries(hicDifferentialAnalysis PRIVATE hicx::core hicx_compile_options)

# hicMergeDomains keeps its algorithm in a second translation unit for the same
# reason: cpp/tests/test_merge_domains.cpp links against it.
add_executable(hicMergeDomains hicMergeDomains.cpp merge_domains_impl.cpp)
target_link_libraries(hicMergeDomains PRIVATE hicx::core hicx_compile_options)

# hicPCA additionally reads and writes bigWig, through the same libBigWig that
# pyBigWig wraps, and uses zlib directly for a gzipped BED extraTrack.
add_executable(hicPCA hicPCA.cpp)
target_link_libraries(hicPCA PRIVATE hicx::core hicx::bigwig hicx::zlib
                                     hicx_compile_options)

# chicExportData writes its bigWig files through the same libBigWig.
add_executable(chicExportData chicExportData.cpp)
target_link_libraries(chicExportData PRIVATE hicx::core hicx::bigwig hicx_compile_options)

# hicCompartmentalization's computational core is in cpp/core
# (hicx/compartmentalization.hpp), shared with the Python API.
add_executable(hicCompartmentalization hicCompartmentalization.cpp)
target_link_libraries(hicCompartmentalization PRIVATE hicx::core hicx_compile_options)

# The `hicexplorer` command (the tool list of list_tools.py). Its table of
# tools and descriptions is generated at build time by HicxToolList.cmake from
# every executable above (its own --help-json) and from the Python drivers of
# pyproject.toml's console scripts (their argparse description), so the list
# is the installed tool set. This block stays last: it lists the executable
# targets defined before it.
set(HICX_PYTHON_DRIVERS hicHyperoptDetectLoops hicHyperoptDetectLoopsHiCCUPS)
get_directory_property(_hicx_listed_targets BUILDSYSTEM_TARGETS)
set(_hicx_tool_list_lines "")
set(_hicx_tool_list_depends "")
foreach(_target IN LISTS _hicx_listed_targets)
    get_target_property(_type ${_target} TYPE)
    if(_type STREQUAL "EXECUTABLE")
        string(APPEND _hicx_tool_list_lines "${_target}|exe|$<TARGET_FILE:${_target}>\n")
        list(APPEND _hicx_tool_list_depends ${_target})
    endif()
endforeach()
foreach(_driver IN LISTS HICX_PYTHON_DRIVERS)
    set(_source "${PROJECT_SOURCE_DIR}/../hicexplorer/${_driver}.py")
    string(APPEND _hicx_tool_list_lines "${_driver}|python|${_source}\n")
    list(APPEND _hicx_tool_list_depends "${_source}")
endforeach()
file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/hicexplorer_tools.txt"
     CONTENT "${_hicx_tool_list_lines}")
add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/hicexplorer_tool_list.inc"
    COMMAND "${CMAKE_COMMAND}" "-DLIST=${CMAKE_CURRENT_BINARY_DIR}/hicexplorer_tools.txt"
            "-DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/hicexplorer_tool_list.inc"
            -P "${PROJECT_SOURCE_DIR}/cmake/HicxToolList.cmake"
    DEPENDS ${_hicx_tool_list_depends} "${CMAKE_CURRENT_BINARY_DIR}/hicexplorer_tools.txt"
            "${PROJECT_SOURCE_DIR}/cmake/HicxToolList.cmake"
    COMMENT "Collecting the tool descriptions for the hicexplorer command"
    VERBATIM)
add_executable(hicexplorer hicexplorer.cpp "${CMAKE_CURRENT_BINARY_DIR}/hicexplorer_tool_list.inc")
target_include_directories(hicexplorer PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(hicexplorer PRIVATE hicx::core hicx_compile_options)
