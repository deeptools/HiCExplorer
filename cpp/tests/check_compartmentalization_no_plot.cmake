# Run by ctest as hicCompartmentalization_no_plot; see tests/CMakeLists.txt.
#
# --noPlot is a C++-only option, so the equivalence harness cannot compare it
# against the Python: with it the tool must exit 0, write the _dat file and
# the --outputMatrix archive, and draw no figure.

file(REMOVE_RECURSE "${OUT}")
file(MAKE_DIRECTORY "${OUT}")

set(common -m "${DATA}/hicPCA/obsexp_norm.h5"
           --pca "${DATA}/hicCompartmentalization/pca1.bedgraph"
           -o "${OUT}/plot.png" --outputMatrix "${OUT}/matrix.npz")

execute_process(COMMAND "${TOOL}" ${common} --noPlot
                RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "--noPlot exited ${status}\n${err}")
endif()
foreach(name plot.png_dat matrix.npz)
    if(NOT EXISTS "${OUT}/${name}")
        message(FATAL_ERROR "--noPlot did not write ${name}")
    endif()
endforeach()
if(EXISTS "${OUT}/plot.png")
    message(FATAL_ERROR "--noPlot wrote a figure")
endif()
file(REMOVE_RECURSE "${OUT}")
