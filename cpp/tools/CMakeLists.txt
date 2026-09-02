# One executable per HiCExplorer tool. Add new tools to the list below; each
# name must match a <name>.cpp in this directory and the Python tool it ports.

set(HICX_TOOLS
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
    hicSumMatrices)

foreach(tool IN LISTS HICX_TOOLS)
    add_executable(${tool} "${tool}.cpp")
    target_link_libraries(${tool} PRIVATE hicx::core hicx_compile_options)
endforeach()

# hicPCA additionally reads and writes bigWig, through the same libBigWig that
# pyBigWig wraps, and uses zlib directly for a gzipped BED extraTrack.
add_executable(hicPCA hicPCA.cpp)
target_link_libraries(hicPCA PRIVATE hicx::core hicx::bigwig hicx::zlib
                                     hicx_compile_options)
