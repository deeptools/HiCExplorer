# Run by ctest as hicAggregateContacts_refuses_plot; see tests/CMakeLists.txt.
#
#   1. Without --noPlot the tool must exit non-zero, say that plotting is not
#      available in the C++ port, and write nothing: no figure, no tables.
#   2. --diagnosticHeatmapFile names a second figure and is refused the same
#      way even with --noPlot.
#   3. With --noPlot the same command must exit 0, write the two tables, and
#      still no figure.

file(REMOVE_RECURSE "${OUT}")
file(MAKE_DIRECTORY "${OUT}")

set(common --matrix "${DATA}/Li_et_al_2015.h5"
           --BED "${DATA}/hicAggregateContacts/test_regions.bed"
           --mode intra-chr --range 50000:900000 --numberOfBins 30
           --outFileName "${OUT}/aggregate.png"
           --outFilePrefixMatrix "${OUT}/m" --outFileContactPairs "${OUT}/p")

function(expect_refusal)
    execute_process(COMMAND "${TOOL}" ${common} ${ARGN}
                    RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(status EQUAL 0)
        message(FATAL_ERROR "a figure was requested (${ARGN}) and the tool exited 0\n${err}")
    endif()
    if(NOT err MATCHES "plotting is not yet available in the C\\+\\+ port")
        message(FATAL_ERROR "the refusal does not name the missing plotting support:\n${err}")
    endif()
    file(GLOB leftovers "${OUT}/*")
    if(leftovers)
        message(FATAL_ERROR "the refused run (${ARGN}) wrote: ${leftovers}")
    endif()
endfunction()

expect_refusal()
expect_refusal(--noPlot --diagnosticHeatmapFile "${OUT}/heatmap.png")

execute_process(COMMAND "${TOOL}" ${common} --noPlot
                RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "--noPlot exited ${status}\n${err}")
endif()
foreach(name m_genome.tab p_genome.tab)
    if(NOT EXISTS "${OUT}/${name}")
        message(FATAL_ERROR "--noPlot did not write ${name}")
    endif()
endforeach()
if(EXISTS "${OUT}/aggregate.png")
    message(FATAL_ERROR "--noPlot wrote a figure")
endif()
file(REMOVE_RECURSE "${OUT}")
