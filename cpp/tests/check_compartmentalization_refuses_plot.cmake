# Run by ctest as hicCompartmentalization_refuses_plot; see tests/CMakeLists.txt.
#
#   1. Without --noPlot the tool must exit non-zero, say that plotting is not
#      available in the C++ port, and write none of its three outputs.
#   2. With --noPlot the same command must exit 0 and write the _dat file and
#      the --outputMatrix archive, and still no figure.

file(REMOVE_RECURSE "${OUT}")
file(MAKE_DIRECTORY "${OUT}")

set(common -m "${DATA}/hicPCA/obsexp_norm.h5"
           --pca "${DATA}/hicCompartmentalization/pca1.bedgraph"
           -o "${OUT}/plot.png" --outputMatrix "${OUT}/matrix.npz")

execute_process(COMMAND "${TOOL}" ${common}
                RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(status EQUAL 0)
    message(FATAL_ERROR "the plot was requested and the tool exited 0\n${err}")
endif()
if(NOT err MATCHES "plotting is not yet available in the C\\+\\+ port")
    message(FATAL_ERROR "the refusal does not name the missing plotting support:\n${err}")
endif()
foreach(name plot.png plot.png_dat matrix.npz)
    if(EXISTS "${OUT}/${name}")
        message(FATAL_ERROR "the refused run wrote ${name}")
    endif()
endforeach()

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
