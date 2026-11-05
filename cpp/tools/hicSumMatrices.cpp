// Port of hicexplorer/hicSumMatrices.py.
//
// Adds two or more contact matrices of the same shape and writes the result
// through the file handler of the first input.
//
// Two pieces of the Python behaviour are defects that are reproduced on
// purpose, both pinned by hicexplorer/test/general/test_hicSumMatrices.py:
//
//  * hicSumMatrices.py:72 masks the union of the two inputs' NaN bins after
//    the addition. Masking deletes rows and columns, so every summed count
//    that lives in a NaN row or column of *either* input is discarded. On the
//    GSM2644945 plus GSM2644947 pair the plain sparse sum holds 3,157,763
//    entries and the file that is written holds 3,152,621, so 5,142 real
//    counts are lost. See test_sum_drops_every_entry_that_lives_in_a_nan_bin.
//  * the masking step also turns an integer matrix into a float64 one, because
//    restoreMaskedBins pads with an empty float64 block. See
//    test_sum_of_an_integer_matrix_is_written_as_float64.
//
// hicSumMatrices.py:60-67 wraps the addition in a try/except that reports a
// shape mismatch. That branch cannot fire: the chrBinBoundaries comparison
// above it rejects every pair whose shapes differ, because the boundaries
// encode the shape (test_the_shape_mismatch_branch_is_unreachable_for_valid_matrices).
// It is therefore not reproduced; hicx::add raises on a shape mismatch and the
// message would be a different one, which no caller can observe.

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/matrix_ops.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicSumMatrices --matrices .h5 or cooler file format\n"
    "                      [.h5 or cooler file format ...] --outFileName\n"
    "                      OUTFILENAME [-h] [--version]\n";

const char* const kHelp =
    "\n"
    "Adds Hi-C matrices of the same size. Format has to be hdf5 (.h5) or npz. In\n"
    "order to minimize the the loss of information, it is recommended to to sum\n"
    "uncorrected matrices (before hicCorrectMatrix).\n"
    "\n"
    "Required arguments:\n"
    "  --matrices .h5 or cooler file format [.h5 or cooler file format ...], "
    "-m .h5 or cooler file format [.h5 or cooler file format ...]\n"
    "                        Space-delimited names of the matrices to add. The\n"
    "                        matrices must have the same shape/size. You can verify\n"
    "                        their size by using `hicInfo`. (default: None)\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        File name to save the resulting matrix. The output is\n"
    "                        from the same file type as the input. Please add the\n"
    "                        file ending suffix (either .h5 or .cool), if it is not\n"
    "                        given, there will be no output. (default: None)\n"
    "\n"
    "Optional arguments:\n"
    "  -h, --help            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::vector<std::string> matrices;
    std::string out_file_name;
};

// hicSumMatrices.py parse_arguments.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicSumMatrices",
                       "Adds Hi-C matrices of the same size. In order to minimize the loss of "
                       "information, it is recommended to sum uncorrected matrices (before "
                       "hicCorrectMatrix).");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrices", "-m"})
        .metavar(".h5 or cooler file format")
        .nargs("+")
        .required()
        .input({"h5", "cool", "mcool"})
        .help("Space-delimited names of the matrices to add. The matrices must have the same "
              "shape/size.");
    required.add({"--outFileName", "-o"})
        .required()
        .output({"h5", "cool"})
        .help("File name to save the resulting matrix. The output is from the same file type "
              "as the input.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"-h", "--help"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.matrices = ns.strs("matrices");
    args.out_file_name = ns.str("outFileName");
    return args;
}

// repr() of the list of chromosome names, which is how the Python error
// message renders list(hic.chrBinBoundaries).
std::string chromosome_list_repr(
    const std::vector<std::pair<std::string, hicx::BinRange>>& boundaries) {
    std::string text = "[";
    for (std::size_t i = 0; i < boundaries.size(); ++i) {
        if (i != 0) {
            text += ", ";
        }
        text += "'" + boundaries[i].first + "'";
    }
    text += "]";
    return text;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    try {
        hicx::ToolMatrix hic = hicx::ToolMatrix::load(args.matrices[0]);
        std::set<std::int64_t> nan_bins(hic.nan_bins().begin(), hic.nan_bins().end());

        for (std::size_t i = 1; i < args.matrices.size(); ++i) {
            const hicx::ToolMatrix other = hicx::ToolMatrix::load(args.matrices[i]);
            if (hic.boundaries() != other.boundaries()) {
                // hicexplorer/__init__.py:3 calls logging.basicConfig, so
                // log.error prints with the default levelname:name: prefix.
                std::fprintf(
                    stderr,
                    "ERROR:hicexplorer.hicSumMatrices:The two matrices have different "
                    "chromosome order. Use the tool `hicConvertFormat` to change the "
                    "order.\n%s: %s\n%s: %s\n",
                    args.matrices[0].c_str(), chromosome_list_repr(hic.boundaries()).c_str(),
                    args.matrices[i].c_str(),
                    chromosome_list_repr(other.boundaries()).c_str());
                return 1;
            }
            hic.matrix() = hicx::add(hic.matrix(), other.matrix());
            if (!other.nan_bins().empty()) {
                nan_bins.insert(other.nan_bins().begin(), other.nan_bins().end());
            }
        }

        // maskBins(sorted(nan_bins)) and the restoreMaskedBins that save()
        // performs, as one step. This is where the counts are lost.
        const std::vector<std::int64_t> mask(nan_bins.begin(), nan_bins.end());
        hicx::mask_and_restore_bins(hic.data(), mask);
        hic.save(args.out_file_name);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicSumMatrices: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicSumMatrices");
    return 0;
}
