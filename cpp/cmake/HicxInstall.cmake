# Install layout of the Python package (pip wheel). The executables go to
# hicexplorer/_bin, which hicexplorer/_cpp.py searches first, and the drawing
# package to hicexplorer/plot, which hicx::plot::draw finds as <tools>/../plot.
#
#   -DHICX_INSTALL_PYTHON_PACKAGE=ON
#
# The install prefix is the root of the wheel (scikit-build-core sets it).

include_guard(GLOBAL)

get_directory_property(_hicx_tool_targets DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/../tools" BUILDSYSTEM_TARGETS)
foreach(_target IN LISTS _hicx_tool_targets)
    get_target_property(_type ${_target} TYPE)
    if(_type STREQUAL "EXECUTABLE")
        install(TARGETS ${_target} RUNTIME DESTINATION hicexplorer/_bin)
    endif()
endforeach()

install(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/../../plot/hicexplorer_plot"
        DESTINATION hicexplorer/plot
        PATTERN "__pycache__" EXCLUDE)
