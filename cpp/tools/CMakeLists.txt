# One executable per HiCExplorer tool. Add new tools to the list below; each
# name must match a <name>.cpp in this directory and the Python tool it ports.

set(HICX_TOOLS
    hicAggregateContacts
    hicAdjustMatrix
    hicBuildMatrix
    hicBuildMatrixMicroC
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
    hicTransform)

foreach(tool IN LISTS HICX_TOOLS)
    add_executable(${tool} "${tool}.cpp")
    target_link_libraries(${tool} PRIVATE hicx::core hicx_compile_options)
endforeach()

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
