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

#include "hicx/argparse.hpp"
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

// hicCreateThresholdFile.py parse_arguments. argparse reports missing required
// arguments before unrecognised ones, so a bare `--version` prints "the
// following arguments are required"; the argument layer does the same.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicCreateThresholdFile", "");
    parser.set_usage(kUsage).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--thresholdValue", "-tv"})
        .type("float")
        .required()
        .help("Standard threshold value for all relative distances.");
    required.add({"--range"})
        .type("int")
        .nargs(2)
        .required()
        .help("Defines the region upstream and downstream of a reference point which should be "
              "included.");
    required.add({"--resolution", "-r"})
        .type("int")
        .default_value(1000)
        .help("Resolution of the bin in genomic units.");
    required.add({"--outFileName", "-o"})
        .required()
        .output({"txt"})
        .help("The name and path of the created threshold file.");

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.threshold_value = ns.real("thresholdValue");
    const std::vector<std::int64_t> range = ns.integers("range");
    args.range_upstream = range[0];
    args.range_downstream = range[1];
    args.resolution = ns.integer("resolution");
    args.out_file_name = ns.str("outFileName");
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
