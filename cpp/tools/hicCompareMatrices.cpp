// Port of hicexplorer/hicCompareMatrices.py.
//
// Takes exactly two matrices, optionally normalises each by the sum of its own
// stored values, and writes their difference, ratio or log2 ratio.
//
// Details of the Python that are behaviour rather than incidental, all pinned
// by hicexplorer/test/general/test_hicCompareMatrices.py:
//
//  * the normaliser is matrix.data.sum(), the sum over the data array of the
//    *symmetric* matrix, not over the stored upper triangle and not
//    matrix.sum(). hicx::data_sum reproduces its accumulation order exactly.
//  * the ratio is computed as a * (1 / b), not as a / b. The two differ by an
//    ulp on a fraction of the entries, which the Python test pins to the entry
//    (test_ratio_with_normalisation_is_one_up_to_one_ulp), so the order is
//    reproduced literally.
//  * inverting b and multiplying means an entry survives only where *both*
//    matrices are non zero. On the GSM pair that is 2,277,549 of the 3,152,621
//    entries the difference keeps.
//  * log2ratio runs eliminate_zeros twice, so a ratio of exactly one, whose
//    log2 is exactly zero, is dropped from the output.
//  * the mask step at the end costs the same counts and the same integer dtype
//    as it does in hicSumMatrices; see matrix_ops.hpp.
//
// np.log2 is numpy's own float64 implementation, not libm's. Measured on this
// corpus the two disagree on about a fifth of the values by up to one ulp, so
// the log2ratio outputs are equal to the reference to roughly 1e-16 relative
// rather than bit for bit. That is thirteen orders of magnitude inside the ED
// gate of cpp/PLAN.md 5.0; the diff and ratio outputs are bit identical.

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "hicx/matrix_ops.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicCompareMatrices --matrices matrix.h5 matrix.h5 --outFileName\n"
    "                          OUTFILENAME [--operation {diff,ratio,log2ratio}]\n"
    "                          [--noNorm] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Takes two matrices as input, normalizes them and applies the given operation.\n"
    "To normalize the matrices each element is divided by the sum of the matrix.\n"
    "\n"
    "Required arguments:\n"
    "  --matrices matrix.h5 matrix.h5, -m matrix.h5 matrix.h5\n"
    "                        Name of the matrices in .h5 format to use, separated\n"
    "                        by a space. (default: None)\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        File name to save the resulting matrix. The output is\n"
    "                        also a .h5 file. (default: None)\n"
    "\n"
    "Optional arguments:\n"
    "  --operation {diff,ratio,log2ratio}\n"
    "                        Operation to apply to the matrices (Default:\n"
    "                        log2ratio).\n"
    "  --noNorm              Do not apply normalisation before computing the\n"
    "                        operation (Default: False). (default: False)\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::vector<std::string> matrices;
    std::string out_file_name;
    std::string operation = "log2ratio";
    bool no_norm = false;
};

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicCompareMatrices: error: %s\n", message.c_str());
    std::exit(2);
}

Arguments parse_arguments(int argc, char** argv) {
    Arguments args;
    bool matrices_seen = false;
    bool out_seen = false;
    std::string* pending_value = nullptr;
    int matrices_remaining = 0;

    for (int i = 1; i < argc; ++i) {
        const std::string token(argv[i]);

        if (pending_value != nullptr) {
            *pending_value = token;
            pending_value = nullptr;
            continue;
        }

        const bool is_option =
            token.size() > 1 && token[0] == '-' &&
            std::isdigit(static_cast<unsigned char>(token[1])) == 0;
        if (!is_option) {
            if (matrices_remaining > 0) {
                args.matrices.push_back(token);
                --matrices_remaining;
                continue;
            }
            fail("unrecognized arguments: " + token);
        }

        matrices_remaining = 0;
        std::string name = token;
        std::optional<std::string> inline_value;
        const std::size_t equals = token.find('=');
        if (equals != std::string::npos && token.rfind("--", 0) == 0) {
            name = token.substr(0, equals);
            inline_value = token.substr(equals + 1);
        }

        if (name == "-h" || name == "--help") {
            std::fputs(kUsage, stdout);
            std::fputs(kHelp, stdout);
            std::exit(0);
        }
        if (name == "--version") {
            std::printf("hicCompareMatrices %s\n", hicx::kVersion);
            std::exit(0);
        }
        if (name == "--noNorm") {
            args.no_norm = true;
            continue;
        }
        if (name == "-m" || name == "--matrices") {
            matrices_seen = true;
            args.matrices.clear();
            matrices_remaining = 2;
            if (inline_value.has_value()) {
                fail("argument --matrices/-m: expected 2 arguments");
            }
            continue;
        }
        if (name == "-o" || name == "--outFileName") {
            out_seen = true;
            if (inline_value.has_value()) {
                args.out_file_name = *inline_value;
            } else {
                pending_value = &args.out_file_name;
            }
            continue;
        }
        if (name == "--operation") {
            if (inline_value.has_value()) {
                args.operation = *inline_value;
            } else {
                pending_value = &args.operation;
            }
            continue;
        }
        fail("unrecognized arguments: " + token);
    }

    if (pending_value != nullptr) {
        fail("expected one argument");
    }
    if (!matrices_seen && !out_seen) {
        fail("the following arguments are required: --matrices/-m, --outFileName/-o");
    }
    if (!matrices_seen) {
        fail("the following arguments are required: --matrices/-m");
    }
    if (!out_seen) {
        fail("the following arguments are required: --outFileName/-o");
    }
    if (args.matrices.size() != 2) {
        fail("argument --matrices/-m: expected 2 arguments");
    }
    if (args.operation != "diff" && args.operation != "ratio" &&
        args.operation != "log2ratio") {
        fail("argument --operation: invalid choice: '" + args.operation +
             "' (choose from 'diff', 'ratio', 'log2ratio')");
    }
    return args;
}

// str() of OrderedDict.keys(), which is what the Python error message
// interpolates for the two chromosome orders.
std::string odict_keys_repr(
    const std::vector<std::pair<std::string, hicx::BinRange>>& boundaries) {
    std::string text = "odict_keys([";
    for (std::size_t i = 0; i < boundaries.size(); ++i) {
        if (i != 0) {
            text += ", ";
        }
        text += "'" + boundaries[i].first + "'";
    }
    text += "])";
    return text;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    try {
        hicx::ToolMatrix hic1 = hicx::ToolMatrix::load(args.matrices[0]);
        hicx::ToolMatrix hic2 = hicx::ToolMatrix::load(args.matrices[1]);

        if (hic1.matrix().rows() != hic2.matrix().rows() ||
            hic1.matrix().cols() != hic2.matrix().cols()) {
            std::fprintf(stderr,
                         "The two matrices have different size. Use matrices having "
                         "the same resolution and created usingthe same parameters. "
                         "Check the matrix values using the tool `hicInfo`.\n");
            return 1;
        }
        if (hic1.boundaries() != hic2.boundaries()) {
            std::fprintf(stderr,
                         "The two matrices have different chromosome order. Use the "
                         "tool `hicAdjustMatrix` to change the order.\n%s: %s\n%s: %s\n",
                         args.matrices[0].c_str(),
                         odict_keys_repr(hic1.boundaries()).c_str(),
                         args.matrices[1].c_str(),
                         odict_keys_repr(hic2.boundaries()).c_str());
            return 1;
        }

        if (!args.no_norm) {
            // matrix.data = matrix.data.astype(float) / matrix.data.sum().
            // The sum is taken in the dtype the matrix still has, so an
            // integer matrix is normalised by an exact integer total.
            const hicx::Scalar total1 = hicx::data_sum(hic1.matrix());
            const hicx::Scalar total2 = hicx::data_sum(hic2.matrix());
            hicx::divide_data_in_place(hic1.matrix(), total1);
            hicx::divide_data_in_place(hic2.matrix(), total2);
        }

        std::set<std::int64_t> nan_bins(hic1.nan_bins().begin(), hic1.nan_bins().end());
        nan_bins.insert(hic2.nan_bins().begin(), hic2.nan_bins().end());

        hicx::CsrMatrix result;
        if (args.operation == "diff") {
            result = hicx::subtract(hic1.matrix(), hic2.matrix());
        } else {
            hicx::reciprocal_data_in_place(hic2.matrix());
            result = hicx::multiply_elementwise(hic1.matrix(), hic2.matrix());
            result.eliminate_zeros();
            if (args.operation == "log2ratio") {
                hicx::log2_data_in_place(result);
                result.eliminate_zeros();
            }
        }

        hic1.matrix() = std::move(result);
        const std::vector<std::int64_t> mask(nan_bins.begin(), nan_bins.end());
        hicx::mask_and_restore_bins(hic1.data(), mask);
        hic1.save(args.out_file_name);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicCompareMatrices: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicCompareMatrices");
    return 0;
}
