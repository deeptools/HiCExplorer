# libhicx4core: matrix model, file formats and the numpy/scipy compatible
# numerics that the textual output of the tools depends on.

# Take the version from the Python package so that both implementations report
# the same string.
set(_version_file "${CMAKE_CURRENT_SOURCE_DIR}/../../hicexplorer/_version.py")
set(HICX_PYTHON_VERSION "unknown")
if(EXISTS "${_version_file}")
    file(READ "${_version_file}" _version_content)
    string(REGEX MATCH "__version__[ ]*=[ ]*'([^']+)'" _match "${_version_content}")
    if(CMAKE_MATCH_1)
        set(HICX_PYTHON_VERSION "${CMAKE_MATCH_1}")
    endif()
endif()
message(STATUS "hicexplorer4: reference version ${HICX_PYTHON_VERSION}")
configure_file(include/hicx/version.hpp.in
               "${CMAKE_CURRENT_BINARY_DIR}/generated/hicx/version.hpp" @ONLY)

# hicQC.html's template, compiled in from the Python package's copy
# (hicexplorer/qc_template.html) so that the C++ writes the report itself;
# editing the template reconfigures.
set(_qc_template "${CMAKE_CURRENT_SOURCE_DIR}/../../hicexplorer/qc_template.html")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_qc_template}")
file(READ "${_qc_template}" _qc_template_hex HEX)
string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," _qc_template_bytes "${_qc_template_hex}")
file(CONFIGURE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/qc_template_html.inc"
     CONTENT "// Generated from hicexplorer/qc_template.html by cpp/core/CMakeLists.txt.\nconst unsigned char kQcTemplateHtml[] = {${_qc_template_bytes}};\n")

# coolercpp (BSD-3-Clause, see cpp/THIRD_PARTY_NOTICES.md): cooler compatible
# cool and mcool access, vendored from the commit named there. Its public
# headers keep the include path coolercpp/<header>.
set(HICX_COOLERCPP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/io/coolercpp")
set(HICX_COOLERCPP_SOURCES
    io/coolercpp/src/api.cpp
    io/coolercpp/src/create.cpp
    io/coolercpp/src/dtype.cpp
    io/coolercpp/src/fastio.cpp
    io/coolercpp/src/h5.cpp
    io/coolercpp/src/json.cpp
    io/coolercpp/src/numpy_compat.cpp
    io/coolercpp/src/parallel.cpp
    io/coolercpp/src/rangequery.cpp
    io/coolercpp/src/rangequery_api.cpp
    io/coolercpp/src/region.cpp
    io/coolercpp/src/table.cpp
    io/coolercpp/src/util.cpp)

# hicfilecpp (MIT, see cpp/THIRD_PARTY_NOTICES.md): the Juicer .hic reader and
# writer, vendored from the commit named there. Public headers keep the
# include path hicfilecpp/<header>.
set(HICX_HICFILECPP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/io/hicfilecpp")
set(HICX_HICFILECPP_SOURCES
    io/hicfilecpp/src/codec.cpp
    io/hicfilecpp/src/io.cpp
    io/hicfilecpp/src/norms.cpp
    io/hicfilecpp/src/output.cpp
    io/hicfilecpp/src/reader.cpp
    io/hicfilecpp/src/writer.cpp)

# hic2cool (MIT, see cpp/THIRD_PARTY_NOTICES.md): .hic to cool and cool to .hic
# conversion, called from hic_adapter.cpp for hicConvertFormat. It is a library
# of the core only; no converter executable is built from it.
set(HICX_HIC2COOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/io/hic2cool")
set(HICX_HIC2COOL_SOURCES
    io/hic2cool/src/convert.cpp
    io/hic2cool/src/cool2hic.cpp
    io/hic2cool/src/extract_norms.cpp
    io/hic2cool/src/hic_header.cpp
    io/hic2cool/src/pyformat.cpp
    io/hic2cool/src/runtime.cpp
    io/hic2cool/src/update.cpp)

add_library(hicx4core
    src/adjust_ops.cpp
    src/aggregate_contacts.cpp
    src/average_regions.cpp
    src/argparse.cpp
    src/bam_file.cpp
    src/bigwig_track.cpp
    src/bedtools_ops.cpp
    src/bins.cpp
    src/clustering.cpp
    src/cnb.cpp
    src/compartmentalization.cpp
    src/correct_ops.cpp
    src/detect_loops.cpp
    src/dist_vs_counts.cpp
    src/build_matrix.cpp
    src/build_matrix_api.cpp
    src/build_matrix_run.cpp
    src/chic_hdf5.cpp
    src/chic_viewpoint.cpp
    src/chicago.cpp
    src/cool_adapter.cpp
    src/correlate_api.cpp
    src/fasta_reader.cpp
    src/find_tads.cpp
    src/h5_file.cpp
    src/hdf5_util.cpp
    src/hic_adapter.cpp
    src/hic_matrix.cpp
    src/iterative_correction.cpp
    src/json_lite.cpp
    src/lbfgsb.cpp
    src/lbfgsb_scipy.cpp
    src/math/ice.cpp
    src/math/kr_balancing.cpp
    src/math/sparse_kernels.cpp
    src/matrix_ops.cpp
    src/merge_domains.cpp
    src/merge_loops.cpp
    src/npz_file.cpp
    src/numpy_compat.cpp
    src/numpy_random.cpp
    src/numpy_sort.cpp
    src/obsexp_ops.cpp
    src/pca_tracks.cpp
    src/pairs_file.cpp
    src/qc_html.cpp
    src/matrix_qc.cpp
    src/matrix_qc_report.cpp
    src/qc_report_api.cpp
    src/plot_bridge.cpp
    src/plot_tiles.cpp
    src/plot_matrix_data.cpp
    src/py_text.cpp
    src/read_bed.cpp
    src/genomedisco.cpp
    src/hicspector.cpp
    src/quasar.cpp
    src/reduce_matrix.cpp
    src/reduce_matrix_entries.cpp
    src/resource_usage.cpp
    src/scc.cpp
    src/simd_reduce.cpp
    src/scipy_stats.cpp
    src/tad_contacts.cpp
    src/tar_gz.cpp
    src/tool_list.cpp
    src/sparse_matrix.cpp
    src/stats_ops.cpp
    src/text_formats.cpp
    src/text_table.cpp
    src/tool_matrix.cpp
    src/transform_ops.cpp
    src/validate_locations.cpp
    src/unidecode_table.cpp
    src/utilities.cpp
    src/dgeev_selected.cpp
    src/sparse_pca.cpp
    src/svl.cpp
    ${HICX_COOLERCPP_SOURCES}
    ${HICX_HICFILECPP_SOURCES}
    ${HICX_HIC2COOL_SOURCES})

target_include_directories(hicx4core PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${HICX_COOLERCPP_DIR}/include"
    "${HICX_HICFILECPP_DIR}/include"
    "${HICX_HIC2COOL_DIR}/include"
    "${CMAKE_CURRENT_BINARY_DIR}/generated")

target_link_libraries(hicx4core
    PUBLIC hicx::hdf5
    PRIVATE hicx::blosc hicx::zlib hicx::hts hicx::lapack
            hicx_compile_options)

# clustering.cpp loads OpenBLAS with dlopen on first use instead of linking it,
# so that a tool which may cluster does not start the BLAS threads at load.
get_filename_component(_hicx_openblas_real "${HICX_LAPACK_LIBRARY}" REALPATH)
set_source_files_properties(src/clustering.cpp PROPERTIES
    COMPILE_DEFINITIONS "HICX_OPENBLAS_LIBRARY=\"${_hicx_openblas_real}\"")
target_link_libraries(hicx4core PRIVATE ${CMAKE_DL_LIBS})

# libBigWig for the bigwig tracks of the Python API's hicPlotMatrix.plotBigwig
# (bigwig_track.cpp), the library pyBigWig wraps.
target_link_libraries(hicx4core PRIVATE hicx::bigwig)

# transform_ops runs the covariance rows on std::thread workers.
find_package(Threads REQUIRED)
target_link_libraries(hicx4core PUBLIC Threads::Threads)

# Spectra and Eigen for hicPCA --eigenSolver lanczos (sparse_pca.cpp only).
include("${CMAKE_CURRENT_SOURCE_DIR}/../cmake/HicxSpectra.cmake")
get_target_property(_hicx_spectra_includes hicx_spectra INTERFACE_INCLUDE_DIRECTORIES)
set(_hicx_spectra_options "")
foreach(_dir IN LISTS _hicx_spectra_includes)
    list(APPEND _hicx_spectra_options "-isystem${_dir}")
endforeach()
set_source_files_properties(src/sparse_pca.cpp PROPERTIES COMPILE_OPTIONS "${_hicx_spectra_options}")

# Boost.Math, header only, for the special functions of hicx/stats_ops.hpp
# (erfc, lgamma, digamma, ibeta, gamma_q, gamma_p_inv). The last digits of
# these functions differ between Boost versions, so the version the outputs
# were verified with is required: the headers come from the conda package
# libboost-headers 1.85 and other versions are refused. Standalone mode keeps
# the rest of Boost out of the include path.
find_path(HICX_BOOST_INCLUDE_DIR
    NAMES boost/math/special_functions/beta.hpp
    HINTS "${HICX_DEPS_ROOT}/include"
    NO_DEFAULT_PATH)
if(NOT HICX_BOOST_INCLUDE_DIR)
    message(FATAL_ERROR "Boost.Math headers not found under ${HICX_DEPS_ROOT}/include "
                        "(conda package libboost-headers 1.85)")
endif()
file(STRINGS "${HICX_BOOST_INCLUDE_DIR}/boost/version.hpp" _hicx_boost_version
     REGEX "^#define BOOST_VERSION [0-9]+$")
if(NOT _hicx_boost_version MATCHES "BOOST_VERSION 1085[0-9][0-9]$")
    message(FATAL_ERROR "Boost 1.85 is required for hicx/stats_ops.hpp, found: ${_hicx_boost_version}")
endif()
add_library(hicx_boost_math INTERFACE)
target_include_directories(hicx_boost_math SYSTEM INTERFACE "${HICX_BOOST_INCLUDE_DIR}")
target_compile_definitions(hicx_boost_math INTERFACE BOOST_MATH_STANDALONE=1)
target_link_libraries(hicx4core PRIVATE hicx_boost_math)

add_library(hicx::core ALIAS hicx4core)
