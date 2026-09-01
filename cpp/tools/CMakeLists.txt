# One executable per HiCExplorer tool. Add new tools to the list below; each
# name must match a <name>.cpp in this directory and the Python tool it ports.

set(HICX_TOOLS
    hicAdjustMatrix
    hicCompareMatrices
    hicConvertFormat
    hicInfo
    hicMergeMatrixBins
    hicMergeTADbins
    hicSumMatrices)

foreach(tool IN LISTS HICX_TOOLS)
    add_executable(${tool} "${tool}.cpp")
    target_link_libraries(${tool} PRIVATE hicx::core hicx_compile_options)
endforeach()
