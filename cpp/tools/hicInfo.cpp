// Port of hicexplorer/hicInfo.py.
//
// Prints information about one or several matrices. Two paths exist:
//
//   * cool file and metadata allowed: everything comes from the cooler file
//     attributes, the matrix itself is never read
//   * h5 file, or cool with --no_metadata: the matrix is loaded through the
//     HiCMatrix layer and size, sum, extrema and NaN bins are computed
//
// The output is byte identical to the Python implementation; see
// hicexplorer/test/general/test_hicInfo.py for the pinned reference blocks.

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "hicx/cool_adapter.hpp"
#include "hicx/hic_matrix.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicInfo --matrices MATRICES [MATRICES ...] [--outFileName OUTFILENAME]\n"
    "               [--no_metadata] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Prints information about a matrix or matrices including matrix size,\n"
    "number of elements, sum of elements, etc.\n"
    "An example usage is:\n"
    "$ hicInfo -m matrix1.h5 matrix2.h5 matrix3.h5\n"
    "\n"
    "Required arguments:\n"
    "  --matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]\n"
    "                        The matrix (or multiple matrices) to get information\n"
    "                        about. HiCExplorer supports the following file\n"
    "                        formats: h5 (native HiCExplorer format) and cool.\n"
    "\n"
    "Optional arguments:\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        File name to save information of the matrix instead of\n"
    "                        writing it to the bash.\n"
    "  --no_metadata, -nm    Do not use meta data from cooler file to display\n"
    "                        information. This method is slower and was the default\n"
    "                        until version 2.2 of HiCExplorer. H5 files always use\n"
    "                        this parameter.\n"
    "  --help, -h            Show this help message and exit.\n"
    "  --version, -v         show program's version number and exit\n";

struct Arguments {
    std::vector<std::string> matrices;
    std::string out_file_name;
    // argparse stores --no_metadata with action='store_false', so the flag is
    // true unless the user asks for the slow path.
    bool use_metadata = true;
};

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicInfo: error: %s\n", message.c_str());
    std::exit(2);
}

Arguments parse_arguments(int argc, char** argv) {
    Arguments args;
    bool matrices_seen = false;
    std::string* pending_value = nullptr;
    bool collecting_matrices = false;

    for (int i = 1; i < argc; ++i) {
        std::string token(argv[i]);

        if (pending_value != nullptr) {
            *pending_value = token;
            pending_value = nullptr;
            continue;
        }

        const bool is_option = token.size() > 1 && token[0] == '-' &&
                               !(token.size() > 1 && std::isdigit(static_cast<unsigned char>(token[1])) != 0);
        if (!is_option) {
            if (collecting_matrices) {
                args.matrices.push_back(token);
                continue;
            }
            fail("unrecognized arguments: " + token);
        }

        collecting_matrices = false;
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
        if (name == "-v" || name == "--version") {
            std::printf("hicInfo %s\n", hicx::kVersion);
            std::exit(0);
        }
        if (name == "-nm" || name == "--no_metadata") {
            args.use_metadata = false;
            continue;
        }
        if (name == "-m" || name == "--matrices") {
            matrices_seen = true;
            collecting_matrices = true;
            if (inline_value.has_value()) {
                args.matrices.push_back(*inline_value);
                collecting_matrices = false;
            }
            continue;
        }
        if (name == "-o" || name == "--outFileName") {
            if (inline_value.has_value()) {
                args.out_file_name = *inline_value;
            } else {
                pending_value = &args.out_file_name;
            }
            continue;
        }
        fail("unrecognized arguments: " + token);
    }

    if (pending_value != nullptr) {
        fail("expected one argument");
    }
    if (!matrices_seen) {
        fail("the following arguments are required: --matrices/-m");
    }
    if (args.matrices.empty()) {
        fail("argument --matrices/-m: expected at least one argument");
    }
    return args;
}

// "{}".format(x) on a numpy scalar goes through float(), so a float32 value is
// printed with the repr of the double it widens to.
std::string format_scalar(const hicx::Scalar& value) {
    return value.is_integer() ? std::to_string(value.integer_value)
                              : hicx::npy::float_repr(value.float_value);
}

std::string format_json_number(const hicx::json::Value& value) {
    if (value.type() == hicx::json::Type::Int) {
        return std::to_string(value.as_int());
    }
    return hicx::npy::float_repr(value.as_double());
}

// A cooler attribute that hicInfo pushes through hicexplorer.utilities.toString.
std::optional<std::string> string_info(const hicx::CoolFile& cool, const std::string& key) {
    const hicx::json::Value* value = cool.info_value(key);
    if (value == nullptr || value->is_null()) {
        return std::nullopt;
    }
    return value->is_string() ? value->as_string() : value->to_python_string();
}

std::string describe(const std::string& matrix_path, const Arguments& args) {
    std::optional<std::string> generated_by;
    std::optional<std::string> genome_assembly;
    std::optional<std::string> statistics;
    std::optional<std::string> generated_by_cooler_lib;
    std::optional<std::string> tool_url;
    std::optional<std::string> matrix_generated_by;
    std::optional<std::string> matrix_generated_by_url;
    std::optional<std::string> creation_date;
    std::optional<std::string> bin_length;
    std::optional<std::int64_t> size;
    std::optional<std::int64_t> nchroms;
    std::optional<std::int64_t> num_non_zero;
    std::optional<std::string> min_non_zero;
    std::optional<std::string> max_non_zero;
    std::optional<std::string> sum_elements;
    std::optional<std::int64_t> num_nan_bins;
    std::vector<std::pair<std::string, std::int64_t>> chromosome_sizes;

    const bool cooler_input = hicx::check_cooler(matrix_path);

    if (cooler_input && args.use_metadata) {
        const hicx::CoolFile cool(matrix_path);
        if (const hicx::json::Value* value = cool.info_value("bin-size");
            value != nullptr && !value->is_null()) {
            bin_length = format_json_number(*value);
        }
        if (const hicx::json::Value* value = cool.info_value("nbins"); value != nullptr) {
            size = value->as_int();
        }
        if (const hicx::json::Value* value = cool.info_value("nchroms"); value != nullptr) {
            nchroms = value->as_int();
        }
        if (const hicx::json::Value* value = cool.info_value("nnz"); value != nullptr) {
            num_non_zero = value->as_int();
        }
        if (const hicx::json::Value* value = cool.info_value("min-value");
            value != nullptr && !value->is_null()) {
            min_non_zero = format_json_number(*value);
        }
        if (const hicx::json::Value* value = cool.info_value("max-value");
            value != nullptr && !value->is_null()) {
            max_non_zero = format_json_number(*value);
        }
        generated_by = string_info(cool, "generated-by");
        genome_assembly = string_info(cool, "genome-assembly");
        if (const hicx::json::Value* metadata = cool.info_value("metadata");
            metadata != nullptr && metadata->is_object()) {
            if (const hicx::json::Value* entry = metadata->find("statistics");
                entry != nullptr) {
                statistics = entry->is_string() ? entry->as_string()
                                                : entry->to_python_string();
            }
        }
        generated_by_cooler_lib = string_info(cool, "generated-by-cooler-lib");
        tool_url = string_info(cool, "tool-url");
        matrix_generated_by = string_info(cool, "matrix-generated-by");
        matrix_generated_by_url = string_info(cool, "matrix-generated-by-url");
        creation_date = string_info(cool, "creation-date");
        if (const hicx::json::Value* value = cool.info_value("sum-elements");
            value != nullptr && !value->is_null()) {
            sum_elements = format_json_number(*value);
        }

        const std::vector<std::string>& names = cool.chrom_names();
        const std::vector<std::int64_t>& lengths = cool.chrom_lengths();
        for (std::size_t i = 0; i < names.size() && i < lengths.size(); ++i) {
            chromosome_sizes.emplace_back(names[i], lengths[i]);
        }
    } else {
        const hicx::HiCMatrix hic = hicx::HiCMatrix::load(matrix_path);
        size = hic.matrix().rows();
        num_non_zero = static_cast<std::int64_t>(hic.matrix().nnz());
        const hicx::Scalar total = hic.matrix().sum();
        const hicx::Scalar diagonal = hic.matrix().diagonal_sum();
        // ((matrix.sum() - matrix.diagonal().sum()) / 2) + matrix.diagonal().sum()
        // The subtraction happens in the dtype of the matrix, the division by
        // the Python int 2 then promotes to float64.
        double half = 0.0;
        switch (total.kind) {
            case hicx::DType::Integer:
                half = static_cast<double>(total.integer_value - diagonal.integer_value) /
                       2.0;
                break;
            case hicx::DType::Float32:
                half = static_cast<double>(static_cast<float>(total.float_value) -
                                           static_cast<float>(diagonal.float_value)) /
                       2.0;
                break;
            case hicx::DType::Float64:
            default:
                half = (total.as_double() - diagonal.as_double()) / 2.0;
                break;
        }
        sum_elements = hicx::npy::float_repr(half + diagonal.as_double());
        bin_length = std::to_string(hic.bin_size());
        num_nan_bins = static_cast<std::int64_t>(hic.nan_bins().size());
        min_non_zero = format_scalar(hic.matrix().data_min());
        max_non_zero = format_scalar(hic.matrix().data_max());
        chromosome_sizes = hic.chromosome_sizes();
    }

    std::string out;
    out += "# Matrix information file. Created with HiCExplorer's hicInfo version ";
    out += hicx::kVersion;
    out += "\n";
    out += "File:\t" + matrix_path + "\n";
    if (creation_date.has_value()) {
        out += "Date:\t" + *creation_date + "\n";
    }
    if (genome_assembly.has_value()) {
        out += "Genome assembly:\t" + *genome_assembly + "\n";
    }
    if (size.has_value()) {
        out += "Size:\t" + hicx::npy::int_with_thousands_separator(*size) + "\n";
    }
    if (bin_length.has_value()) {
        out += "Bin_length:\t" + *bin_length + "\n";
    }
    if (sum_elements.has_value()) {
        out += "Sum of matrix:\t" + *sum_elements + "\n";
    }
    out += "Chromosomes:length: ";
    for (const auto& [name, length] : chromosome_sizes) {
        out += name + ": " + std::to_string(length) + " bp; ";
    }
    out += "\n";
    if (nchroms.has_value()) {
        out += "Number of chromosomes:\t" + std::to_string(*nchroms) + "\n";
    }
    if (num_non_zero.has_value()) {
        out += "Non-zero elements:\t" +
               hicx::npy::int_with_thousands_separator(*num_non_zero) + "\n";
    }
    if (min_non_zero.has_value()) {
        out += "Minimum (non zero):\t" + *min_non_zero + "\n";
    }
    if (max_non_zero.has_value()) {
        out += "Maximum:\t" + *max_non_zero + "\n";
    }
    if (num_nan_bins.has_value()) {
        out += "NaN bins:\t" + std::to_string(*num_nan_bins) + "\n";
    }
    if (cooler_input) {
        const hicx::CoolFile cool(matrix_path);
        out += "The following columns are available: " +
               hicx::npy::array_str(cool.bin_columns()) + "\n";
    }
    if (generated_by.has_value()) {
        out += "\n\nGenerated by:\t" + *generated_by + "\n";
    }
    if (generated_by_cooler_lib.has_value()) {
        out += "Cooler library version:\t" + *generated_by_cooler_lib + "\n";
    }
    if (tool_url.has_value()) {
        out += "HiCMatrix url:\t" + *tool_url + "\n";
    }
    if (matrix_generated_by.has_value()) {
        out += "Interaction matrix created with:\t" + *matrix_generated_by + "\n";
    }
    if (matrix_generated_by_url.has_value()) {
        out += "URL:\t" + *matrix_generated_by_url + "\n";
    }
    if (statistics.has_value()) {
        out += "\n\nBuild statistics:\n" + *statistics + "\n";
    }
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    try {
        for (const std::string& matrix_path : args.matrices) {
            const std::string information = describe(matrix_path, args);
            if (!args.out_file_name.empty()) {
                // The Python implementation reopens the file for every matrix,
                // so only the last block survives.
                std::ofstream out(args.out_file_name, std::ios::binary | std::ios::trunc);
                if (!out) {
                    std::fprintf(stderr, "hicInfo: cannot write %s\n",
                                 args.out_file_name.c_str());
                    return 1;
                }
                out << information;
            } else {
                std::cout << information << "\n";
            }
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicInfo: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicInfo");
    return 0;
}
