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
//
// Deliberate v4-only deviation from the Python (project owner, 2026-09-22),
// NOT a faithfulness fix: the bug above has a real consequence beyond
// dropping some bins. hicMergeMatrixBins.py derives the merged bin layout
// purely from the bins the *input* matrix happens to have, never from the
// true chromosome length. Two inputs of the same genome at the same
// resolution and the same --numBins produce differently shaped output
// (different total bin count, different per-chromosome bin counts) whenever
// one of them is missing a few bins the other has, for example because a
// trailing region has no observed reads in one sample and does in the other.
// That makes the outputs not directly comparable for a reason that has
// nothing to do with the genome being represented.
//
// v4 fixes this by requiring --chromosomeSizes / -cs and using it, together
// with the matrix's own bin size, to build the merge groups from the true,
// fixed genome layout (hicx::plan_bin_merge_genome / merge_bins_genome in
// reduce_matrix.hpp) instead of from whichever bins the input happens to
// hold. A bin the genome layout expects but the input does not have simply
// contributes nothing to its group's sum, exactly like a present but
// all-zero bin already would, rather than shrinking the group count. Because
// the layout no longer depends on what the input happens to contain, the
// numBins/2 drop quirk above is not reproduced in this path either: every
// group the genome layout produces, however short, is kept, so the shape is
// determined only by the genome, the resolution and --numBins.
//
// This only applies to a matrix with a single, fixed bin size. A
// restriction-fragment matrix has no "bin a genomic position belongs to"
// without the restriction cut positions, which a chromosome sizes file does
// not carry, so for one of those (BinTable::bin_size_homogeneous() false)
// this tool falls back to the Python's own input-derived plan_bin_merge,
// numBins/2 drop and all: --chromosomeSizes is still required (uniformly,
// for every matrix this tool merges) but its content plays no role for that
// matrix. A user comparing this tool's output to hicMergeMatrixBins.py on the
// exact reproduced-bug input (fixed resolution, some bins missing) will see
// a different, genome-length-consistent shape here; on a restriction-fragment
// input the two still agree, quirk included.

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
#include "hicx/bins.hpp"
#include "hicx/reduce_matrix.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/text_formats.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicMergeMatrixBins --matrix matrix.h5 --outFileName OUTFILENAME\n"
    "                          --numBins int --chromosomeSizes txt file\n"
    "                          [--runningWindow] [-h] [--version]\n";

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
    "  --chromosomeSizes txt file, -cs txt file\n"
    "                        File with the chromosome sizes for your genome, a\n"
    "                        plain name<TAB>length text file. C++ port only,\n"
    "                        required: unlike hicMergeMatrixBins.py, this tool\n"
    "                        builds the merged bin layout from the true\n"
    "                        chromosome lengths, not from the bins the input\n"
    "                        matrix happens to have, so that two inputs of the\n"
    "                        same genome and resolution always merge into the\n"
    "                        same output shape however much data either of them\n"
    "                        is missing. Not used for a matrix with irregular\n"
    "                        (for example restriction-fragment) bin widths,\n"
    "                        which falls back to the Python's own layout, but\n"
    "                        still required for every matrix this tool merges.\n"
    "                        (default: None)\n"
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
    std::string chromosome_sizes;
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
    required.add({"--chromosomeSizes", "-cs"})
        .file_type("r")
        .metavar("txt file")
        .input({"txt"})
        .required()
        .cpp_only(
            "The C++ port builds the merged bin layout from the true chromosome lengths, "
            "not from the bins the input matrix happens to have, so that two inputs of the "
            "same genome and resolution always merge into the same output shape however "
            "much data either of them is missing (project owner, 2026-09-22). Required for "
            "every run, including --runningWindow, which does not need it, and a "
            "restriction-fragment matrix, which cannot use it (see the top-of-file "
            "comment) and falls back to the Python's own layout.")
        .help("File with the chromosome sizes for your genome, a plain "
              "name<TAB>length text file.");
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
    args.chromosome_sizes = ns.str("chromosomeSizes");
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
            // --chromosomeSizes is required by the CLI for every run (see the
            // top-of-file comment), but it only replaces the layout for a
            // matrix with a single, fixed bin size. A restriction-fragment
            // matrix has no bin-from-position rule a chromosome-lengths file
            // could supply, so it keeps the Python's own input-derived
            // grouping.
            const hicx::BinTable bin_table(hic.data().cut_intervals);
            if (bin_table.bin_size_homogeneous()) {
                const std::vector<std::pair<std::string, std::int64_t>> chromosome_sizes =
                    hicx::read_chromosome_sizes(args.chromosome_sizes);
                hicx::MatrixData merged = hicx::merge_bins_genome(
                    hic.data(), args.num_bins, chromosome_sizes, bin_table.bin_size());
                hic.data() = std::move(merged);
            } else {
                hicx::MatrixData merged = hicx::merge_bins(hic.data(), args.num_bins);
                hic.data() = std::move(merged);
            }
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
