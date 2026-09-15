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
    hicPlotTADs
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

# hicDifferentialTAD and hicInterIntraTAD cut their TAD blocks with the same
# code, which is also what cpp/tests/test_tad_contacts.cpp links against.
foreach(tool IN ITEMS hicDifferentialTAD hicInterIntraTAD)
    add_executable(${tool} "${tool}.cpp" tad_contacts_impl.cpp)
    target_link_libraries(${tool} PRIVATE hicx::core hicx_compile_options)
endforeach()

# hicMergeDomains keeps its algorithm in a second translation unit for the same
# reason: cpp/tests/test_merge_domains.cpp links against it.
add_executable(hicMergeDomains hicMergeDomains.cpp merge_domains_impl.cpp)
target_link_libraries(hicMergeDomains PRIVATE hicx::core hicx_compile_options)

# hicPCA additionally reads and writes bigWig, through the same libBigWig that
# pyBigWig wraps, and uses zlib directly for a gzipped BED extraTrack.
add_executable(hicPCA hicPCA.cpp)
target_link_libraries(hicPCA PRIVATE hicx::core hicx::bigwig hicx::zlib
                                     hicx_compile_options)

# hicCompartmentalization keeps its computational core in a second translation
# unit, like hicDetectLoops, so that cpp/tests/test_compartmentalization.cpp
# links against exactly the code the tool runs.
add_executable(hicCompartmentalization hicCompartmentalization.cpp
                                       compartmentalization_impl.cpp)
target_link_libraries(hicCompartmentalization PRIVATE hicx::core hicx_compile_options)
