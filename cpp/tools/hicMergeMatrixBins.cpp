// Port of hicexplorer/hicMergeMatrixBins.py.
//
// Two ways of coarsening a matrix, sharing one preliminary step:
//
//  * merge_bins, which cuts a new bin every --numBins bins and at every
//    chromosome change and sums each block. This is hicx::merge_bins, the same
//    code hicConvertFormat calls for every resolution of an mcool, so the two
//    tools cannot drift apart.
//  * --runningWindow, which keeps the resolution and spreads every count over
//    the --numBins x --numBins window around it. This is hicx::running_window.
//
// Both are preceded by remove_nans_if_needed, which masks the NaN bins away
// for good and drops the correction factors with them.
//
// Four behaviours are reproduced rather than fixed, all pinned by
// hicexplorer/test/general/test_hicMergeMatrixBins.py:
//
//  * hicMergeMatrixBins.py:248 drops a group holding fewer than numBins/2 bins
//    together with its counts. At --numBins 20 on small_test_matrix_50kb_res.h5
//    that silently removes chr2LHet (8 bins) and chrXHet (5 bins) from the
//    output entirely, plus four trailing part groups, and the total count falls
//    from 59,360 to 52,061. No warning above debug level is printed.
//  * --runningWindow ignores chromosome borders: the window is applied to raw
//    bin indices, so the last bins of one chromosome pick up counts from the
//    first bins of the next.
//  * --runningWindow with an even --numBins ends the run with an
//    AssertionError. Reproduced as an error message and exit status 1.
//  * --numBins 1 with --runningWindow returns before the window is applied, so
//    the only thing that happens is that the NaN bins are dropped.
//
// The output format follows the *input*, not the name given to --outFileName,
// because hiCMatrix.save reuses the handler built during the load. An h5 input
// named out.cool therefore produces out.cool.h5.

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "hicx/adjust_ops.hpp"
#include "hicx/argparse.hpp"
#include "hicx/reduce_matrix.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicMergeMatrixBins --matrix matrix.h5 --outFileName OUTFILENAME\n"
    "                          --numBins int [--runningWindow] [-h] [--version]\n";

const char* const kHelp =
    "\n"
    "Merges bins from a Hi-C matrix. For example, using a matrix containing 5kb\n"
    "bins, a matrix of 50kb bins can be derived using --numBins 10. From one type\n"
    "of downstream analysis to another, different bin sizes are used. For example\n"
    "to call TADs, unmerged matrices are recommended while to display Hi-C\n"
    "matrices, bins of approximately 2000bp usually yield the best representations\n"
    "with `hicPlotMatrix` for small regions, and even larger bins (50kb) are\n"
    "recommended for whole chromosome representations or for\n"
    "`hicPlotDistVsCounts`.\n"
    "\n"
    "Required arguments:\n"
    "  --matrix matrix.h5, -m matrix.h5\n"
    "                        Matrix to reduce in h5 format. (default: None)\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        File name to save the resulting matrix. The output is\n"
    "                        also a .h5 file. But don't add the suffix. (default:\n"
    "                        None)\n"
    "  --numBins int, -nb int\n"
    "                        Number of bins to merge. (default: None)\n"
    "\n"
    "Optional arguments:\n"
    "  --runningWindow       Set to merge for using a running window of length\n"
    "                        --numBins. (default: False)\n"
    "  -h, --help            Show this help message and exit.\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::string matrix;
    std::string out_file_name;
    std::int64_t num_bins = 0;
    bool running_window = false;
};

// hicMergeMatrixBins.py parse_arguments.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicMergeMatrixBins",
                       "Merges bins from a Hi-C matrix. For example, using a matrix containing "
                       "5kb bins, a matrix of 50kb bins can be derived using --numBins 10.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrix", "-m"})
        .metavar("matrix.h5")
        .required()
        .input({"h5", "cool", "mcool"})
        .help("Matrix to reduce in h5 format.");
    required.add({"--outFileName", "-o"})
        .required()
        .output({"h5", "cool"})
        .help("File name to save the resulting matrix. The output format follows the input.");
    required.add({"--numBins", "-nb"})
        .metavar("int")
        .type("int")
        .required()
        .help("Number of bins to merge.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--runningWindow"})
        .action(cli::Action::StoreTrue)
        .help("Set to merge for using a running window of length --numBins.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("Show this help message and exit.");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.matrix = ns.str("matrix");
    args.out_file_name = ns.str("outFileName");
    args.num_bins = ns.integer("numBins");
    args.running_window = ns.flag("runningWindow");
    return args;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    try {
        hicx::ToolMatrix hic = hicx::ToolMatrix::load(args.matrix);

        // remove_nans_if_needed: the NaN bins are deleted, not masked, and the
        // correction factors go with them because they no longer line up.
        if (!hic.nan_bins().empty()) {
            const std::vector<std::int64_t> nan_bins = hic.nan_bins();
            hicx::delete_bins(hic.data(), nan_bins);
            hic.data().correction_factors.reset();
            hic.data().correction_factors_are_column = false;
            hic.refresh_boundaries();
            std::fputs("WARNING:hicexplorer.hicMergeMatrixBins:*WARNING*: The matrix "
                       "is probably a corrected matrix that contains NaN bins. This "
                       "bins can not be merged and are removed. It is preferable to "
                       "first merge bins in a uncorrected  matrix and then correct "
                       "the matrix. Correction factors, if present, are removed as "
                       "well.\n",
                       stderr);
        }

        if (args.running_window) {
            if (args.num_bins != 1) {
                if (args.num_bins % 2 == 0) {
                    // The Python assert at hicMergeMatrixBins.py:139.
                    std::fputs("hicMergeMatrixBins: AssertionError: num_bins has to "
                               "be an odd number\n",
                               stderr);
                    return 1;
                }
                hic.matrix() = hicx::running_window(hic.matrix(), args.num_bins);
                hic.data().nan_bins = hicx::empty_column_bins(hic.matrix());
            }
            // num_bins == 1 returns before the window is applied, so nan_bins
            // stays whatever remove_nans_if_needed left, which is empty.
        } else {
            hicx::MatrixData merged = hicx::merge_bins(hic.data(), args.num_bins);
            hic.data() = std::move(merged);
            hic.refresh_boundaries();
        }

        hic.save(args.out_file_name);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicMergeMatrixBins: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicMergeMatrixBins");
    return 0;
}
