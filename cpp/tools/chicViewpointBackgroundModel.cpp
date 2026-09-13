// Port of hicexplorer/chicViewpointBackgroundModel.py.
//
// For every relative position to a reference point, over all reference points
// of all matrices: the smoothed viewpoint values are collected into one
// distribution, a negative binomial is fitted to it, and its maximum and its
// mean, the mean normalised by the sum of all positions' means, are written
// next to the fitted size and prob.
//
// The fit is the one decision that needs explaining. fit_nbinom.fit runs
// scipy's L-BFGS-B on a likelihood that is nearly flat along `size` for these
// distributions, and the size it returns is sensitive to the last bits of its
// sums: refitting the same values in another order moves it by a median of 55
// percent. The port therefore runs a statement by statement translation of the
// same L-BFGS-B (hicx/lbfgsb_scipy.hpp) rather than the projected BFGS the
// other tools use, and the fitted size and prob are compared at class EN, the
// envelope of the reference against itself (cpp/PLAN.md 5.7). The position,
// maximum and mean columns are exact.
//
// Behaviour reproduced as the Python has it:
//
//   * compute_background zips range(view_point_range_start,
//     view_point_range_end) with the data. The range end is the bin of the
//     region's end and is exclusive here, while computeViewpoint included it,
//     so the last element of every viewpoint is dropped. For a reference point
//     spanning several bins the data is also shorter than the range by the
//     collapsed bins, which shifts every downstream value one relative
//     position per collapsed bin towards the reference point.
//   * The positions are relative bin distances, multiplied by the bin size of
//     the last matrix when written.
//   * --truncateZeros can empty a distribution; the Python then writes size 10,
//     prob nan, max 0 and mean 0 for it.
//   * A failure in any reference point is reported only after every matrix
//     has been processed, and then the tool exits 1 without an output file.
//
// Threading: the fits, one per position, are independent and run on --threads
// threads into one slot each; that is where the time goes. Measured on the
// test data (43 reference points, two matrices, 1,001 fits), three runs each:
// 0.247 s of CPU and 0.257 s of wall time at --threads 1, 0.267 s of CPU and
// 0.107 s of wall time at --threads 16, with the peak RSS at 45 MB either way.
// Running the reference points on threads as well saved another 0.017 s of
// wall time and cost 0.010 s of CPU, which is not worth a second parallel
// section, so they run sequentially and are merged in file order. That is the
// order the Python merges them in at --threads 1. At a higher thread count the
// Python merges worker results in the order they arrive, which changes the
// order of the values in each distribution and with it the last bits of the
// fit's sums, and so the fitted size (measured: --threads 8 writes different
// digits from --threads 1 on 987 of 1,001 lines). The C++ always uses the
// --threads 1 order and is byte-identical at every thread count.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "chic_arguments.hpp"
#include "hicx/chic_viewpoint.hpp"
#include "hicx/lbfgsb_scipy.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/parallel.hpp"
#include "hicx/version.hpp"

namespace {

using hicx::chic::ReferencePoint;

const char* const kUsage =
    "usage: chicViewpointBackgroundModel --matrices MATRICES [MATRICES ...]\n"
    "                                    --referencePoints REFERENCEPOINTS\n"
    "                                    [--averageContactBin AVERAGECONTACTBIN]\n"
    "                                    [--truncateZeros] [--outFileName OUTFILENAME]\n"
    "                                    [--threads THREADS]\n"
    "                                    [--fixateRange FIXATERANGE] [--help]\n"
    "                                    [--version]\n";

struct ReferencePointSlot {
    std::vector<std::pair<std::int64_t, double>> values;
    bool failed = false;
    std::string message;
};

struct PositionResult {
    hicx::stats::NBinomFit fit;
    double max = 0.0;
    double average = 0.0;
    bool empty = false;
};

}  // namespace

int main(int argc, char** argv) {
    using hicx::chic_cli::Arity;
    using hicx::chic_cli::Option;
    using hicx::chic_cli::Type;

    hicx::chic_cli::Parser parser(
        "chicViewpointBackgroundModel", kUsage,
        "chicViewpointBackgroundModel computes a background model for all given samples with "
        "all reference points.",
        hicx::kVersion);
    parser.add(Option{{"--matrices", "-m"}, "matrices", Arity::OneOrMore, Type::String, true, {},
                      "The input matrices (samples) to build the background model on."});
    parser.add(Option{{"--referencePoints", "-rp"}, "referencePoints", Arity::One, Type::String,
                      true, {}, "Bed file with the reference points."});
    parser.add(Option{{"--averageContactBin"}, "averageContactBin", Arity::One, Type::Int, false,
                      {"5"}, "Average the contacts of n bins via a sliding window approach."});
    parser.add(Option{{"--truncateZeros", "-tz"}, "truncateZeros", Arity::Flag, Type::String,
                      false, {}, "Truncates the zeros before the distributions are fitted."});
    parser.add(Option{{"--outFileName", "-o"}, "outFileName", Arity::One, Type::String, false,
                      {"background_model.txt"}, "The name of the background model file."});
    parser.add(Option{{"--threads", "-t"}, "threads", Arity::One, Type::Int, false, {"4"},
                      "Number of threads."});
    parser.add(Option{{"--fixateRange", "-fs"}, "fixateRange", Arity::One, Type::Int, false,
                      {"500000"}, "Fixate the background model starting at distance x."});
    const hicx::chic_cli::Parsed args = parser.parse(argc, argv);

    const std::vector<std::string>& matrices = args.list("matrices");
    const std::int64_t average_contact_bin = args.integer("averageContactBin");
    const bool truncate_zeros = args.flag("truncateZeros");
    const std::string out_file_name = args.str("outFileName");
    const std::int64_t threads = args.integer("threads");
    const std::int64_t fixate_range = args.integer("fixateRange");

    try {
        const hicx::chic::ReferencePoints reference_points =
            hicx::chic::read_reference_points(args.str("referencePoints"));
        const std::size_t count = reference_points.points.size();
        if (threads == 0) {
            throw std::runtime_error("ZeroDivisionError: integer division or modulo by zero");
        }
        const unsigned int workers = threads < 1 ? 1U : static_cast<unsigned int>(threads);

        std::map<std::int64_t, std::vector<double>> distributions;
        std::int64_t bin_size = 0;
        bool failed = false;
        std::string failure;
        for (const std::string& matrix_path : matrices) {
            const hicx::chic::ViewpointMatrix matrix =
                hicx::chic::ViewpointMatrix::load(matrix_path);
            bin_size = matrix.bin_size();
            std::vector<ReferencePointSlot> slots(count);
            const auto compute_slot = [&](std::size_t index) {
                ReferencePointSlot& slot = slots[index];
                const ReferencePoint& point = reference_points.points[index];
                try {
                    const hicx::chic::ViewpointRange range = hicx::chic::calculate_viewpoint_range(
                        matrix, point, fixate_range, fixate_range);
                    hicx::chic::ComputedViewpoint viewpoint = hicx::chic::compute_viewpoint(
                        matrix, point, point.chromosome, range.region_start, range.region_end);
                    const std::int64_t view_point_start =
                        matrix.reference_point_indices(point).first;
                    const auto bin_range = matrix.region_bin_range(
                        point.chromosome, range.region_start, range.region_end);
                    if (!bin_range.has_value()) {
                        throw hicx::chic::TypeError("cannot unpack non-iterable NoneType object");
                    }
                    std::vector<double> data = std::move(viewpoint.data);
                    if (average_contact_bin > 0) {
                        data = hicx::chic::smooth_interaction_values(
                            std::span<const double>(data), average_contact_bin);
                    }
                    const std::int64_t span = std::max<std::int64_t>(
                        0, bin_range->second - bin_range->first);
                    const std::size_t pairs =
                        std::min(static_cast<std::size_t>(span), data.size());
                    slot.values.reserve(pairs);
                    for (std::size_t k = 0; k < pairs; ++k) {
                        slot.values.emplace_back(
                            bin_range->first + static_cast<std::int64_t>(k) - view_point_start,
                            data[k]);
                    }
                } catch (const std::exception& error) {
                    slot.failed = true;
                    slot.message = error.what();
                }
            };
            for (std::size_t index = 0; index < count; ++index) {
                compute_slot(index);
            }
            for (ReferencePointSlot& slot : slots) {
                if (slot.failed) {
                    failed = true;
                    failure = slot.message;
                    continue;
                }
                for (const auto& [position, value] : slot.values) {
                    distributions[position].push_back(value);
                }
            }
        }
        if (failed) {
            std::fprintf(stderr,
                         "chicViewpointBackgroundModel: An error occurred caused by one or many "
                         "faulty reference points.\n"
                         "chicViewpointBackgroundModel: Please run chicQualityControl to remove "
                         "these from your reference point file: %s\n"
                         "chicViewpointBackgroundModel: %s\n",
                         args.str("referencePoints").c_str(), failure.c_str());
            return 1;
        }

        std::vector<std::int64_t> positions;
        positions.reserve(distributions.size());
        for (const auto& entry : distributions) {
            positions.push_back(entry.first);
        }
        std::vector<PositionResult> results(positions.size());
        hicx::parallel_for(positions.size(), workers, [&](std::size_t index) {
            const std::vector<double>& all = distributions.at(positions[index]);
            std::vector<double> values;
            if (truncate_zeros) {
                for (const double value : all) {
                    if (value > 0.0) {
                        values.push_back(value);
                    }
                }
            } else {
                values = all;
            }
            PositionResult& result = results[index];
            result.fit = hicx::stats::fit_nbinom_scipy(values);
            if (values.empty()) {
                result.empty = true;
                return;
            }
            // np.max propagates NaN.
            double maximum = values[0];
            for (const double value : values) {
                if (std::isnan(value) || std::isnan(maximum)) {
                    maximum = std::numeric_limits<double>::quiet_NaN();
                } else if (value > maximum) {
                    maximum = value;
                }
            }
            result.max = maximum;
            result.average =
                hicx::npy::pairwise_sum(values.data(), values.size()) /
                static_cast<double>(values.size());
        });

        double sum_all_values = 0.0;
        bool any_numpy_average = false;
        for (const PositionResult& result : results) {
            sum_all_values += result.average;
            any_numpy_average = any_numpy_average || !result.empty;
        }
        if (!positions.empty() && !any_numpy_average) {
            // Every mean is the Python float 0.0, so the normalisation is a
            // Python float division by zero rather than a numpy one.
            throw std::runtime_error("ZeroDivisionError: float division by zero");
        }

        std::ofstream out(out_file_name, std::ios::binary);
        if (!out) {
            throw std::runtime_error("[Errno 2] No such file or directory: '" + out_file_name +
                                     "'");
        }
        out << "Relative position\tsize nbinom\tprob nbinom\tmax value\tmean value\n";
        for (std::size_t index = 0; index < positions.size(); ++index) {
            const PositionResult& result = results[index];
            out << positions[index] * bin_size << "\t"
                << hicx::chic::format_fixed(result.fit.size, 12) << "\t"
                << hicx::chic::format_fixed(result.fit.prob, 12) << "\t"
                << hicx::chic::format_fixed(result.max, 12) << "\t"
                << hicx::chic::format_fixed(result.average / sum_all_values, 12) << "\n";
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chicViewpointBackgroundModel: %s\n", error.what());
        return 1;
    }
    return 0;
}
