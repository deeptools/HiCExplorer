# Install layout of the Python package (pip wheel).
#
#   -DHICX_INSTALL_PYTHON_PACKAGE=ON
#
# The executables are the command-line tools themselves: they go to the
# wheel's scripts directory (SKBUILD_SCRIPTS_DIR, installed as <prefix>/bin),
# so running a tool starts the C++ executable with no Python in between. The
# drawing package goes to the data directory as share/hicexplorer/plot
# (<prefix>/share/hicexplorer/plot), which hicx::plot finds from the binary as
# <bin>/../share/hicexplorer/plot. Outside scikit-build-core the same layout
# is bin/ and share/hicexplorer/plot under CMAKE_INSTALL_PREFIX. The extension
# hicexplorer._hicx goes into the hicexplorer package of the wheel's platlib
# (SKBUILD_PLATLIB_DIR), outside scikit-build-core into
# lib/python<X.Y>/site-packages/hicexplorer.

include_guard(GLOBAL)

if(DEFINED SKBUILD_SCRIPTS_DIR)
    set(_hicx_bin_destination "${SKBUILD_SCRIPTS_DIR}")
else()
    set(_hicx_bin_destination bin)
endif()
if(DEFINED SKBUILD_PLATLIB_DIR)
    set(_hicx_python_destination "${SKBUILD_PLATLIB_DIR}/hicexplorer")
else()
    set(_hicx_python_destination
        "lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/hicexplorer")
endif()
if(DEFINED SKBUILD_DATA_DIR)
    set(_hicx_plot_destination "${SKBUILD_DATA_DIR}/share/hicexplorer/plot")
else()
    set(_hicx_plot_destination share/hicexplorer/plot)
endif()

get_directory_property(_hicx_tool_targets DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/../tools" BUILDSYSTEM_TARGETS)
foreach(_target IN LISTS _hicx_tool_targets)
    get_target_property(_type ${_target} TYPE)
    if(_type STREQUAL "EXECUTABLE")
        install(TARGETS ${_target} RUNTIME DESTINATION "${_hicx_bin_destination}")
    endif()
endforeach()

install(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/../../plot/hicexplorer_plot"
        DESTINATION "${_hicx_plot_destination}"
        PATTERN "__pycache__" EXCLUDE)

install(TARGETS _hicx LIBRARY DESTINATION "${_hicx_python_destination}")
