# Run by ctest as hicAggregateContacts_no_plot; see tests/CMakeLists.txt.
#
# --noPlot is a C++-only option, so the equivalence harness cannot compare it
# against the Python: with it the tool must exit 0, write the two tables, and
# draw neither the aggregate figure nor a requested diagnostic heatmap.

file(REMOVE_RECURSE "${OUT}")
file(MAKE_DIRECTORY "${OUT}")

execute_process(COMMAND "${TOOL}"
                        --matrix "${DATA}/Li_et_al_2015.h5"
                        --BED "${DATA}/hicAggregateContacts/test_regions.bed"
                        --mode intra-chr --range 50000:900000 --numberOfBins 30
                        --outFileName "${OUT}/aggregate.png"
                        --outFilePrefixMatrix "${OUT}/m" --outFileContactPairs "${OUT}/p"
                        --diagnosticHeatmapFile "${OUT}/heatmap.png" --noPlot
                RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "--noPlot exited ${status}\n${err}")
endif()
foreach(name m_genome.tab p_genome.tab)
    if(NOT EXISTS "${OUT}/${name}")
        message(FATAL_ERROR "--noPlot did not write ${name}")
    endif()
endforeach()
foreach(name aggregate.png heatmap.png)
    if(EXISTS "${OUT}/${name}")
        message(FATAL_ERROR "--noPlot wrote ${name}")
    endif()
endforeach()
file(REMOVE_RECURSE "${OUT}")
