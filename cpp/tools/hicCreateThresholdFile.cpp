// Port of hicexplorer/hicCreateThresholdFile.py.
//
// Writes a two line header and one row per bin over the requested range. The
// whole tool is 15 lines of Python and there is nothing numeric in it, but two
// details are worth naming because a port gets them wrong by default:
//
//   * the loop runs to `range[1] + resolution`, exclusive, so the last row is
//     at `range[1]` and is included. Stopping at `range[1]` drops it.
//   * the threshold is written with `'{}'.format(float)`, which is
//     `repr(float)`: the shortest string that round trips. `-tv 1` prints as
//     `1.0` and `-tv 0.00001` as `1e-05`, neither of which printf's %g gives.
//     hicx::npy::float_repr is that function.
//
// The parser has no --help and no --version: hicCreateThresholdFile.py builds
// its ArgumentParser with add_help=False and never adds one back, so `-h` and
// `--version` are unrecognised arguments and the tool exits 2. Reproduced,
// including the order in which argparse complains: it validates what it parsed
// before it mentions the leftovers, so a bare `--version` reports the missing
// required arguments and not the unknown option.
//
// No threading, no SIMD: the tool writes at most a few hundred lines.

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

#include "hicx/numpy_compat.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicCreateThresholdFile --thresholdValue THRESHOLDVALUE --range RANGE\n"
    "                              RANGE [--resolution RESOLUTION] --outFileName\n"
    "                              OUTFILENAME\n";

struct Arguments {
    double threshold_value = 0.0;
    std::int64_t range_upstream = 0;
    std::int64_t range_downstream = 0;
    std::int64_t resolution = 1000;
    std::string out_file_name;
};

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicCreateThresholdFile: error: %s\n", message.c_str());
    std::exit(2);
}

std::int64_t parse_int(const std::string& text, const std::string& option) {
    try {
        std::size_t used = 0;
        const long long value = std::stoll(text, &used);
        if (used != text.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return value;
    } catch (const std::exception&) {
        fail("argument " + option + ": invalid int value: '" + text + "'");
    }
}

Arguments parse_arguments(int argc, char** argv) {
    Arguments args;
    bool threshold_seen = false;
    bool range_seen = false;
    bool out_seen = false;

    std::vector<std::string> tokens;
    for (int i = 1; i < argc; ++i) {
        const std::string token(argv[i]);
        const std::size_t equals = token.find('=');
        if (equals != std::string::npos && token.rfind("--", 0) == 0) {
            tokens.push_back(token.substr(0, equals));
            tokens.push_back(token.substr(equals + 1));
        } else {
            tokens.push_back(token);
        }
    }

    // argparse reports missing required arguments before it reports an
    // unrecognised one, because parse_args validates what it parsed and only
    // then complains about the leftovers. `hicCreateThresholdFile --version`
    // therefore prints "the following arguments are required", not
    // "unrecognized arguments: --version". Reproduced by collecting the
    // unknown tokens and reporting them last.
    std::vector<std::string> unrecognized;

    for (std::size_t i = 0; i < tokens.size(); ++i) {
        const std::string& token = tokens[i];
        const auto next = [&]() -> std::string {
            if (i + 1 >= tokens.size()) {
                fail("argument " + token + ": expected one argument");
            }
            return tokens[++i];
        };
        if (token == "--thresholdValue" || token == "-tv") {
            const std::string text = next();
            try {
                std::size_t used = 0;
                args.threshold_value = std::stod(text, &used);
                if (used != text.size()) {
                    throw std::invalid_argument("trailing characters");
                }
            } catch (const std::exception&) {
                fail("argument --thresholdValue/-tv: invalid float value: '" + text + "'");
            }
            threshold_seen = true;
        } else if (token == "--range") {
            args.range_upstream = parse_int(next(), "--range");
            args.range_downstream = parse_int(next(), "--range");
            range_seen = true;
        } else if (token == "--resolution" || token == "-r") {
            args.resolution = parse_int(next(), "--resolution/-r");
        } else if (token == "--outFileName" || token == "-o") {
            args.out_file_name = next();
            out_seen = true;
        } else {
            unrecognized.push_back(token);
        }
    }

    std::string missing;
    const auto add_missing = [&missing](const char* name) {
        missing += missing.empty() ? name : std::string(", ") + name;
    };
    if (!threshold_seen) {
        add_missing("--thresholdValue/-tv");
    }
    if (!range_seen) {
        add_missing("--range");
    }
    if (!out_seen) {
        add_missing("--outFileName/-o");
    }
    if (!missing.empty()) {
        fail("the following arguments are required: " + missing);
    }
    if (!unrecognized.empty()) {
        std::string joined;
        for (const std::string& token : unrecognized) {
            joined += joined.empty() ? token : " " + token;
        }
        fail("unrecognized arguments: " + joined);
    }
    return args;
}

}  // namespace

int main(int argc, char** argv) {
    Arguments args = parse_arguments(argc, argv);
    try {
        std::ofstream out(args.out_file_name, std::ios::binary | std::ios::trunc);
        if (!out) {
            std::fprintf(stderr, "hicCreateThresholdFile: cannot write '%s'\n",
                         args.out_file_name.c_str());
            return 1;
        }
        const std::string threshold = hicx::npy::float_repr(args.threshold_value);
        out << "# Threshold file of HiCExplorer's hicCreateThresholdFile version "
            << hicx::kVersion << '\n';
        out << "# Standard threshold " << threshold << '\n';

        if (args.range_upstream > 0) {
            args.range_upstream = -args.range_upstream;
        }
        // range(start, stop, step) with stop = downstream + resolution, so the
        // downstream end is included. A non positive step would make the
        // Python loop empty; range() raises for a step of zero.
        if (args.resolution == 0) {
            std::fprintf(stderr,
                         "hicCreateThresholdFile: ValueError: range() arg 3 must not be "
                         "zero\n");
            return 1;
        }
        const std::int64_t stop = args.range_downstream + args.resolution;
        if (args.resolution > 0) {
            for (std::int64_t i = args.range_upstream; i < stop; i += args.resolution) {
                out << i << '\t' << threshold << '\n';
            }
        } else {
            for (std::int64_t i = args.range_upstream; i > stop; i += args.resolution) {
                out << i << '\t' << threshold << '\n';
            }
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicCreateThresholdFile: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicCreateThresholdFile");
    return 0;
}
