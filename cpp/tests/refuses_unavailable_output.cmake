# Runs a tool with an option that requests a file the C++ port cannot produce
# and asserts the interim policy of 2026-09-13: a non-zero exit status, and
# neither the requested file nor the tool's main output written.
#
#   cmake -DTOOL=<binary> -DWORKDIR=<scratch dir> -DARGS=<;-list>
#         -DMUST_NOT_EXIST=<;-list of paths> -P refuses_unavailable_output.cmake

foreach(variable TOOL WORKDIR ARGS MUST_NOT_EXIST)
    if(NOT DEFINED ${variable})
        message(FATAL_ERROR "${variable} is not set")
    endif()
endforeach()

file(REMOVE_RECURSE "${WORKDIR}")
file(MAKE_DIRECTORY "${WORKDIR}")

execute_process(
    COMMAND "${TOOL}" ${ARGS}
    WORKING_DIRECTORY "${WORKDIR}"
    RESULT_VARIABLE status
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err)

if(status EQUAL 0)
    message(FATAL_ERROR "expected a non-zero exit status, got 0\nstderr: ${err}")
endif()
foreach(path IN LISTS MUST_NOT_EXIST)
    if(EXISTS "${path}")
        message(FATAL_ERROR "the refused run wrote '${path}'")
    endif()
endforeach()
file(GLOB leftovers "${WORKDIR}/*")
if(leftovers)
    message(FATAL_ERROR "the refused run left files in its working directory: ${leftovers}")
endif()
message(STATUS "exit status ${status}, nothing written; stderr: ${err}")
