// Shared implementation of hicBuildMatrix and hicBuildMatrixMicroC.
//
// hicexplorer/hicBuildMatrixMicroC.py is hicBuildMatrix.py with the
// restriction enzyme arguments removed and --binSize made required: its main()
// calls the same lib/buildMatrixMethods.createMatrix with pRestrictionCutFile,
// pRestrictionSequence, pDanglingSequence, pMinDistance, pMaxDistance and
// pKeepSelfLigation all None. So the two tools differ in their argument parser
// and in nothing else, and this header is the "nothing else".
//
// The run itself, createMatrix, is hicx::build_matrix::run
// (core/include/hicx/build_matrix_run.hpp), which the public Python function
// buildMatrixMethods.createMatrix calls too; this header adds the tools'
// error printing and the drawing of the QC report. The design notes and the
// list of reproduced defects are in hicBuildMatrix.cpp.

#ifndef HICX_TOOLS_BUILD_MATRIX_IMPL_HPP
#define HICX_TOOLS_BUILD_MATRIX_IMPL_HPP

#include <atomic>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#include "hicx/bam_file.hpp"
#include "hicx/build_matrix.hpp"
#include "hicx/build_matrix_run.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/reduce_matrix.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/version.hpp"

namespace {

// Set by main before anything else, so that the shared code can name the tool
// the user invoked in its usage and error messages.
const char* g_tool = "hicBuildMatrix";

using Arguments = hicx::build_matrix::Arguments;
using hicx::build_matrix::ends_with;
using hicx::build_matrix::PixelAccumulator;
using hicx::build_matrix::read_chromosome_sizes;
using hicx::build_matrix::to_upper;
using hicx::build_matrix::write_matrix_outputs;

// hicexplorer.utilities.genomicRegion as an argparse type: a value that is
// only whitespace becomes None, one that is only separators is an error.
// hicQuickQC includes this file and has no --region.
[[maybe_unused]] std::optional<std::string> genomic_region_check(const std::string& text) {
    bool blank = true;
    for (const char c : text) {
        if (std::isspace(static_cast<unsigned char>(c)) == 0) {
            blank = false;
            break;
        }
    }
    if (!blank && hicx::normalise_region(text).empty()) {
        return text + " is not a valid region";
    }
    return std::nullopt;
}

// createMatrix with the tools' error reporting: exit status 1 and the
// message on stderr.
int run_build_matrix(const Arguments& args) {
    try {
        hicx::build_matrix::run(args);
    } catch (const hicx::build_matrix::Exit& exit) {
        if (exit.logged()) {
            std::fprintf(stderr, "ERROR:hicexplorer.lib.buildMatrixMethods:%s\n", exit.what());
        } else {
            std::fprintf(stderr, "%s\n", exit.what());
        }
        return 1;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s: %s\n", g_tool, error.what());
        return 1;
    }
    hicx::report_resource_usage(g_tool);
    return 0;
}

// buildMatrixMethods.py:1362, QC.main(["-l", QC.log, "-o", QCfolder]): the
// tables are written by run_build_matrix, and the five PNGs and hicQC.html are
// drawn here with hicPrepareQCreport's default --dpi 200. Called by the tools
// after everything else, because on success it does not return.
int draw_qc_report(const std::string& qc_folder) {
    hicx::plot::JsonObject data;
    data.add("outputFolder", hicx::plot::json_string(qc_folder));
    data.add("dpi", hicx::plot::json_int(200));
    return hicx::plot::draw("hicPrepareQCreport", data.str(), std::nullopt);
}

}  // namespace

#endif  // HICX_TOOLS_BUILD_MATRIX_IMPL_HPP
