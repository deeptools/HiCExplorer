// Port of hicexplorer/hicNormalize.py.
//
// Three normalization modes, all of them pure elementwise scaling of the value
// array followed by eliminate_zeros and a save. No reduction is involved apart
// from matrix.sum() in the 'smallest' mode, and that goes through
// CsrMatrix::sum(), which already reproduces scipy's evaluation order.
//
// The one thing that is easy to get wrong and is fully observable in the
// output file is the precision. hicNormalize.py:76, :105 and :131 all start by
// casting the value array to float32, and every arithmetic step afterwards
// happens in single precision:
//
//   * numpy 1.26 uses value based casting for a scalar against an array, so
//     `np.divide(float32_array, np.float64_scalar)` selects the 'ff->f' loop:
//     the scalar is rounded to float32 first and the division is a float32
//     division. The same holds for `data *= python_float` and for
//     `data -= np.float32 min_value`.
//   * min_max_difference is computed as a float32 subtraction and only then
//     widened to float64 (hicNormalize.py:84), which the widening cannot undo.
//
// The port therefore holds the values as double, as the whole core does, but
// rounds through `float` at exactly the points numpy does. `to_float32` is
// that rounding, and every arithmetic step below is written in terms of it.
//
// Two further details of the reference that are reproduced rather than fixed:
//
//   * the NaN and Inf masks are applied twice in the 'smallest' and
//     'multiplicative' branches (hicNormalize.py:107-111 then :116-120, and
//     :133-137 then :142-146). The second application is redundant. It is kept
//     because the first one is inside `if i != argmin` in the smallest branch,
//     so for the argmin matrix only the second one runs.
//   * `--setToZeroThreshold` compares with `<`, so a value sitting exactly on
//     the threshold survives. The Python suite never ran this option; the
//     characterization test added in
//     hicexplorer/test/general/test_hicNormalize.py pins both the threshold
//     and the operator, on a master where 1,706 entries sit exactly on it.
//
// Threading and SIMD: deliberately none. The whole tool is one pass over the
// value array between an HDF5 read and an HDF5 write, and on the designated
// input the C++ run is already 15 times cheaper in CPU time than the Python
// one. cpp/OPTIMIZATION.md 6 requires a measurement before an optimisation,
// and the measurement says the work is in the file layer.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <string>
#include <vector>

#include "hicx/resource_usage.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicNormalize --matrices MATRICES [MATRICES ...] --normalize\n"
    "                    {norm_range,smallest,multiplicative} --outFileName\n"
    "                    FILENAME [FILENAME ...]\n"
    "                    [--multiplicativeValue MULTIPLICATIVEVALUE]\n"
    "                    [--setToZeroThreshold SETTOZEROTHRESHOLD] [--help]\n"
    "                    [--version]\n";

const char* const kHelp =
    "\n"
    "Normalizes given matrices either to the smallest given read number of all "
    "matrices or to 0 - 1 range. However, it does NOT compute the contact "
    "probability.\n"
    "\n"
    "We recommend to compute first the normalization (with hicNormalize) and "
    "correct the data (with hicCorrectMatrix) in a second step.\n"
    "\n"
    "Required arguments:\n"
    "  --matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]\n"
    "                        The matrix (or multiple matrices) to get information\n"
    "                        about. HiCExplorer supports the following file\n"
    "                        formats: h5 (native HiCExplorer format) and cool.\n"
    "  --normalize {norm_range,smallest,multiplicative}, -n "
    "{norm_range,smallest,multiplicative}\n"
    "                        Normalize to a) 0 to 1 range, b) all matrices to the\n"
    "                        lowest read count of the given matrices (Default:\n"
    "                        smallest).\n"
    "  --outFileName FILENAME [FILENAME ...], -o FILENAME [FILENAME ...]\n"
    "                        Output file name for the Hi-C matrix.\n"
    "\n"
    "Optional arguments:\n"
    "  --multiplicativeValue MULTIPLICATIVEVALUE, -mv MULTIPLICATIVEVALUE\n"
    "                        Value to multiply if --normalize is set to\n"
    "                        multiplicative. (Default: 1).\n"
    "  --setToZeroThreshold SETTOZEROTHRESHOLD, -sz SETTOZEROTHRESHOLD\n"
    "                        A threshold to set all values after normalization to 0\n"
    "                        if smaller this threshold. Default value is 0 i.e.\n"
    "                        there is no effect.It is recommended to set it for the\n"
    "                        normalize mode \"smallest\" to 1.0. This parameter will\n"
    "                        influence the sparsity of the matrix by removing many\n"
    "                        values close to 0 in smallest normalization mode.\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::vector<std::string> matrices;
    std::vector<std::string> out_file_names;
    std::string normalize;
    double multiplicative_value = 1.0;
    double set_to_zero_threshold = 0.0;
};

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicNormalize: error: %s\n", message.c_str());
    std::exit(2);
}

double parse_double(const std::string& text, const std::string& option) {
    try {
        std::size_t used = 0;
        const double value = std::stod(text, &used);
        if (used != text.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return value;
    } catch (const std::exception&) {
        fail("argument " + option + ": invalid float value: '" + text + "'");
    }
}

Arguments parse_arguments(int argc, char** argv) {
    Arguments args;
    bool matrices_seen = false;
    bool out_seen = false;
    bool normalize_seen = false;

    enum class Collect { None, Matrices, OutFiles };
    Collect collecting = Collect::None;
    std::string* pending_string = nullptr;
    double* pending_double = nullptr;
    std::string pending_option;

    for (int i = 1; i < argc; ++i) {
        const std::string token(argv[i]);

        if (pending_string != nullptr) {
            *pending_string = token;
            pending_string = nullptr;
            continue;
        }
        if (pending_double != nullptr) {
            *pending_double = parse_double(token, pending_option);
            pending_double = nullptr;
            continue;
        }

        // A negative number is a value, not an option. argparse decides this
        // the same way when the parser holds no numeric-looking option string.
        const bool is_option =
            token.size() > 1 && token[0] == '-' &&
            std::isdigit(static_cast<unsigned char>(token[1])) == 0 && token[1] != '.';
        if (!is_option) {
            if (collecting == Collect::Matrices) {
                args.matrices.push_back(token);
                continue;
            }
            if (collecting == Collect::OutFiles) {
                args.out_file_names.push_back(token);
                continue;
            }
            fail("unrecognized arguments: " + token);
        }

        collecting = Collect::None;
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
            std::printf("hicNormalize %s\n", hicx::kVersion);
            std::exit(0);
        }
        if (name == "-m" || name == "--matrices") {
            matrices_seen = true;
            if (inline_value.has_value()) {
                args.matrices.push_back(*inline_value);
            } else {
                collecting = Collect::Matrices;
            }
            continue;
        }
        if (name == "-o" || name == "--outFileName") {
            out_seen = true;
            if (inline_value.has_value()) {
                args.out_file_names.push_back(*inline_value);
            } else {
                collecting = Collect::OutFiles;
            }
            continue;
        }
        if (name == "-n" || name == "--normalize") {
            normalize_seen = true;
            if (inline_value.has_value()) {
                args.normalize = *inline_value;
            } else {
                pending_string = &args.normalize;
            }
            continue;
        }
        if (name == "-mv" || name == "--multiplicativeValue") {
            if (inline_value.has_value()) {
                args.multiplicative_value =
                    parse_double(*inline_value, "--multiplicativeValue/-mv");
            } else {
                pending_double = &args.multiplicative_value;
                pending_option = "--multiplicativeValue/-mv";
            }
            continue;
        }
        if (name == "-sz" || name == "--setToZeroThreshold") {
            if (inline_value.has_value()) {
                args.set_to_zero_threshold =
                    parse_double(*inline_value, "--setToZeroThreshold/-sz");
            } else {
                pending_double = &args.set_to_zero_threshold;
                pending_option = "--setToZeroThreshold/-sz";
            }
            continue;
        }
        fail("unrecognized arguments: " + token);
    }

    if (pending_string != nullptr || pending_double != nullptr) {
        fail("expected one argument");
    }

    std::string missing;
    const auto add_missing = [&missing](const char* name) {
        missing += missing.empty() ? name : std::string(", ") + name;
    };
    if (!matrices_seen || args.matrices.empty()) {
        add_missing("--matrices/-m");
    }
    if (!normalize_seen) {
        add_missing("--normalize/-n");
    }
    if (!out_seen || args.out_file_names.empty()) {
        add_missing("--outFileName/-o");
    }
    if (!missing.empty()) {
        fail("the following arguments are required: " + missing);
    }
    if (args.normalize != "norm_range" && args.normalize != "smallest" &&
        args.normalize != "multiplicative") {
        fail("argument --normalize/-n: invalid choice: '" + args.normalize +
             "' (choose from 'norm_range', 'smallest', 'multiplicative')");
    }
    return args;
}

// The rounding numpy performs whenever a float32 array takes part in an
// operation. Every arithmetic step of the tool goes through it.
inline double to_float32(double value) {
    return static_cast<double>(static_cast<float>(value));
}

// data = data.astype(np.float32)
void cast_to_float32(hicx::CsrMatrix& matrix) {
    std::vector<double>& data = matrix.mutable_data();
    for (double& value : data) {
        value = to_float32(value);
    }
    matrix.set_dtype("float32");
}

// data[np.isnan(data)] = 0; data[np.isinf(data)] = 0
void clear_non_finite(hicx::CsrMatrix& matrix) {
    for (double& value : matrix.mutable_data()) {
        if (std::isnan(value) || std::isinf(value)) {
            value = 0.0;
        }
    }
}

// The tail every branch shares: eliminate_zeros, then the threshold, then
// eliminate_zeros again (hicNormalize.py:94-98, :121-125, :147-151).
void finish(hicx::CsrMatrix& matrix, double threshold) {
    matrix.eliminate_zeros();
    const float threshold32 = static_cast<float>(threshold);
    for (double& value : matrix.mutable_data()) {
        if (static_cast<float>(value) < threshold32) {
            value = 0.0;
        }
    }
    matrix.eliminate_zeros();
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    try {
        std::vector<hicx::ToolMatrix> matrices;
        matrices.reserve(args.matrices.size());
        std::vector<hicx::Scalar> sums;

        for (const std::string& path : args.matrices) {
            matrices.push_back(hicx::ToolMatrix::load(path));
            if (args.normalize == "smallest") {
                // matrix.sum() of the symmetric matrix, in the dtype the file
                // carries, which is what np.argmin then compares.
                sums.push_back(matrices.back().matrix().sum());
            }
        }

        if (args.matrices.size() > args.out_file_names.size()) {
            // hicNormalize.py indexes args.outFileName[i] without a check and
            // dies with an IndexError. Same exit status, with a message.
            //
            // One deliberate difference: the Python raises inside the loop, so
            // it has already written the outputs it had names for. The check
            // here is made before anything is written, so a run that cannot
            // finish leaves no half-normalized set behind. Recorded because it
            // is a behaviour difference, not a bug fix on the numbers.
            std::fprintf(stderr,
                         "hicNormalize: error: %zu matrices were given but only %zu "
                         "output names; hicNormalize.py:100 indexes "
                         "args.outFileName[i] and raises IndexError here\n",
                         args.matrices.size(), args.out_file_names.size());
            return 1;
        }

        std::size_t argmin = 0;
        if (args.normalize == "smallest") {
            // np.argmin: the first index holding the minimum.
            for (std::size_t i = 1; i < sums.size(); ++i) {
                if (sums[i].as_double() < sums[argmin].as_double()) {
                    argmin = i;
                }
            }
        }

        for (std::size_t i = 0; i < matrices.size(); ++i) {
            hicx::CsrMatrix& matrix = matrices[i].matrix();
            cast_to_float32(matrix);

            if (args.normalize == "norm_range") {
                clear_non_finite(matrix);
                const std::vector<double>& data = matrix.data();
                if (!data.empty()) {
                    // np.min / np.max over the float32 value array.
                    const float min_value = static_cast<float>(
                        *std::min_element(data.begin(), data.end()));
                    const float max_value = static_cast<float>(
                        *std::max_element(data.begin(), data.end()));
                    // hicNormalize.py:84 subtracts in float32 and widens the
                    // result to float64; the divide then rounds it straight
                    // back to float32, so the widening changes nothing.
                    const float min_max_difference = max_value - min_value;
                    for (double& value : matrix.mutable_data()) {
                        const float shifted = static_cast<float>(value) - min_value;
                        value = static_cast<double>(shifted / min_max_difference);
                    }
                }
                clear_non_finite(matrix);
            } else if (args.normalize == "smallest") {
                if (i != argmin) {
                    clear_non_finite(matrix);
                    // sum_list[i] / sum_list[argmin] is a float64 scalar
                    // division; np.divide then rounds it to float32 because
                    // the array is float32.
                    const double adjust_factor =
                        sums[i].as_double() / sums[argmin].as_double();
                    const float factor32 = static_cast<float>(adjust_factor);
                    for (double& value : matrix.mutable_data()) {
                        value = static_cast<double>(static_cast<float>(value) /
                                                    factor32);
                    }
                }
                clear_non_finite(matrix);
            } else {
                clear_non_finite(matrix);
                const float factor32 = static_cast<float>(args.multiplicative_value);
                for (double& value : matrix.mutable_data()) {
                    value =
                        static_cast<double>(static_cast<float>(value) * factor32);
                }
                clear_non_finite(matrix);
            }

            finish(matrix, args.set_to_zero_threshold);
            matrices[i].save(args.out_file_names[i]);
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicNormalize: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicNormalize");
    return 0;
}
