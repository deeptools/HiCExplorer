// Port of hicexplorer/chicViewpoint.py.
//
// For every matrix and every reference point: the viewpoint over --range, its
// values relative to the sum over --fixateRange, the x-fold over the
// background model's mean, and a negative binomial p-value per bin, written to
// one HDF5 file in the layout hicx/chic_hdf5.hpp describes.
//
// Behaviour reproduced as the Python has it:
//
//   * The p-value of every bin except the reference point itself is taken
//     from the background distribution at -fixateRange or +fixateRange,
//     because Viewpoint.pvalues looks a genomic-distance keyed model up by bin
//     distance. See hicx::chic::p_values. Reported as a defect.
//   * Near a chromosome end the background can be one element longer than the
//     data. adjustViewpointData then rebuilds both from range(start, end) with
//     the end exclusive, which drops the last bin, and returns the data as
//     float32 (chicViewpoint.py:104). The float32 values reach the smoothing
//     and, when a viewpoint is shorter than --averageContactBin, the p-value.
//   * A sum over --fixateRange of 0 is treated as absent by
//     computeRelativeValues, which then divides by the sum of the data, 0 as
//     well, so every relative value and x-fold is NaN.
//   * The chromosome group is created for the first reference point on a
//     chromosome; a later reference point whose chromosome group already
//     exists is written into the group created last, whichever chromosome
//     that is (chicViewpoint.py:301-305).
//   * reference_point_start and reference_point_end are taken from the
//     reference point file by the viewpoint's position in the result list,
//     which is its position in the file (chicViewpoint.py:305).
//
// Where the Python does not terminate, the port stops with an error instead:
// Viewpoint.createUniqueHDFGroup retries a taken gene name forever, because it
// renames pAdditionalGroupName but keeps creating file_name (viewpoint.py:
// 355-364), so a reference point file with a gene name twice on one
// chromosome hangs the Python. The C++ exits 1 with a message. A deviation,
// recorded as such.
//
// Threading: measured, and not used. With the 43 reference points of
// referencePoints.bed on both test matrices, worker threads leave the wall
// time where it is (0.217 s at --threads 1, 0.207 s at 16) and raise the CPU
// time from 0.210 s to 0.383 s and the peak RSS from 47 MB to 121 MB, above
// the memory budget, because the work per reference point is small and every
// worker thread brings its own malloc arena. Per cpp/OPTIMIZATION.md 6 the
// loop is sequential; --threads is accepted, and 0 still fails as the
// Python's integer division by it does. Viewpoints are written in file order.
// The Python collects its worker chunks by index, which is the same order, so
// its output does not depend on --threads either.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/chic_hdf5.hpp"
#include "hicx/chic_viewpoint.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/version.hpp"

namespace {

using hicx::chic::BackgroundModel;
using hicx::chic::InteractionFileData;
using hicx::chic::ReferencePoint;
using hicx::chic::ViewpointMatrix;

const char* const kUsage =
    "usage: chicViewpoint --matrices MATRICES [MATRICES ...] --range RANGE RANGE\n"
    "                     --referencePoints REFERENCEPOINTS --backgroundModelFile\n"
    "                     BACKGROUNDMODELFILE [--outFileName OUTFILENAME]\n"
    "                     [--threads THREADS]\n"
    "                     [--averageContactBin AVERAGECONTACTBIN]\n"
    "                     [--fixateRange FIXATERANGE] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Computes per input matrix all viewpoints which are defined in the reference points "
    "file.\n"
    "\n"
    "Required arguments:\n"
    "  --matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]\n"
    "                        Path to the Hi-C matrices which store the captured\n"
    "                        Hi-C data per sample.\n"
    "  --range RANGE RANGE   Defines the region upstream and downstream of a\n"
    "                        reference point which should be considered in the\n"
    "                        analysis. Please remember to use the same fixate range\n"
    "                        setting as for the background model computation and\n"
    "                        that distances of the range larger than the fixate\n"
    "                        range use the background model of those.Format is\n"
    "                        --region upstream downstream\n"
    "  --referencePoints REFERENCEPOINTS, -rp REFERENCEPOINTS\n"
    "                        Reference point file. Needs to be in the format: 'chr\n"
    "                        100' for a single reference point or 'chr 100 200' for\n"
    "                        a reference region and with a single reference point\n"
    "                        per line\n"
    "  --backgroundModelFile BACKGROUNDMODELFILE, -bmf BACKGROUNDMODELFILE\n"
    "                        path to the background file computed by\n"
    "                        chicViewpointBackgroundModel\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        This hdf5 file contains all created viewpoint files.\n"
    "\n"
    "Optional arguments:\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads (uses the python multiprocessing\n"
    "                        module) (Default: 4).\n"
    "  --averageContactBin AVERAGECONTACTBIN\n"
    "                        Average the contacts of n bins via a sliding window\n"
    "                        approach to smooth the values and be less sensitive\n"
    "                        for outliers (Default: 5).\n"
    "  --fixateRange FIXATERANGE, -fs FIXATERANGE\n"
    "                        Fixate range of background model starting at distance\n"
    "                        x. E.g. all values greater 500kb are set to the value\n"
    "                        of the 500kb bin (Default: 500000).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

std::string basename(const std::string& path) {
    const std::size_t slash = path.find_last_of('/');
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

// pDataList / pBackgroundList with numpy broadcasting of one dimension.
std::vector<double> divide(const std::vector<double>& numerator,
                           const std::vector<double>& denominator) {
    if (numerator.size() == denominator.size()) {
        std::vector<double> result(numerator.size());
        for (std::size_t i = 0; i < numerator.size(); ++i) {
            result[i] = numerator[i] / denominator[i];
        }
        return result;
    }
    if (denominator.size() == 1) {
        std::vector<double> result(numerator.size());
        for (std::size_t i = 0; i < numerator.size(); ++i) {
            result[i] = numerator[i] / denominator[0];
        }
        return result;
    }
    if (numerator.size() == 1) {
        std::vector<double> result(denominator.size());
        for (std::size_t i = 0; i < denominator.size(); ++i) {
            result[i] = numerator[0] / denominator[i];
        }
        return result;
    }
    throw hicx::chic::ValueError("operands could not be broadcast together with shapes (" +
                                 std::to_string(numerator.size()) + ",) (" +
                                 std::to_string(denominator.size()) + ",)");
}

struct Settings {
    std::int64_t range_upstream = 0;
    std::int64_t range_downstream = 0;
    std::int64_t fixate_range = 0;
    std::int64_t average_contact_bin = 0;
};

// chicViewpoint.compute_viewpoint for one reference point.
InteractionFileData compute_entry(const ViewpointMatrix& matrix, const ReferencePoint& point,
                                  const std::string& gene, const Settings& settings,
                                  const BackgroundModel& model, const BackgroundModel& means) {
    const hicx::chic::ViewpointRange fixed = hicx::chic::calculate_viewpoint_range(
        matrix, point, settings.fixate_range, settings.fixate_range);
    const hicx::chic::ComputedViewpoint intermediate = hicx::chic::compute_viewpoint(
        matrix, point, point.chromosome, fixed.region_start, fixed.region_end);
    const double denominator =
        hicx::npy::pairwise_sum(intermediate.data.data(), intermediate.data.size());

    const hicx::chic::ViewpointRange range = hicx::chic::calculate_viewpoint_range(
        matrix, point, settings.range_upstream, settings.range_downstream);
    hicx::chic::ComputedViewpoint viewpoint = hicx::chic::compute_viewpoint(
        matrix, point, point.chromosome, range.region_start, range.region_end);

    std::vector<double> background =
        hicx::chic::interaction_background_data(means, range.upstream, range.downstream);
    std::vector<double> data = std::move(viewpoint.data);
    bool data_is_float32 = false;

    if (data.size() != background.size()) {
        // adjustViewpointData
        const std::int64_t view_point_start = matrix.reference_point_indices(point).first;
        const auto bin_range =
            matrix.region_bin_range(point.chromosome, range.region_start, range.region_end);
        if (!bin_range.has_value()) {
            throw hicx::chic::TypeError("cannot unpack non-iterable NoneType object");
        }
        const std::size_t span =
            static_cast<std::size_t>(std::max<std::int64_t>(0, bin_range->second - bin_range->first));
        std::vector<std::int64_t> data_order;
        std::unordered_map<std::int64_t, double> data_values;
        for (std::size_t k = 0; k < std::min(span, data.size()); ++k) {
            const std::int64_t key =
                bin_range->first + static_cast<std::int64_t>(k) - view_point_start;
            if (data_values.emplace(key, data[k]).second) {
                data_order.push_back(key);
            } else {
                data_values[key] = data[k];
            }
        }
        std::vector<std::int64_t> background_order;
        std::unordered_map<std::int64_t, double> background_values;
        for (std::size_t k = 0; k < std::min(span, background.size()); ++k) {
            const std::int64_t key =
                bin_range->first + static_cast<std::int64_t>(k) - view_point_start;
            if (background_values.emplace(key, background[k]).second) {
                background_order.push_back(key);
            } else {
                background_values[key] = background[k];
            }
        }
        for (const std::int64_t key : background_order) {
            if (data_values.emplace(key, 0.0).second) {
                data_order.push_back(key);
            }
        }
        data.clear();
        for (const std::int64_t key : data_order) {
            // np.fromiter(..., dtype=np.float32)
            data.push_back(static_cast<double>(static_cast<float>(data_values[key])));
        }
        background.clear();
        for (const std::int64_t key : background_order) {
            background.push_back(background_values[key]);
        }
        data_is_float32 = true;
    }

    if (settings.average_contact_bin > 0 &&
        static_cast<std::int64_t>(data.size()) >= settings.average_contact_bin) {
        if (data_is_float32) {
            std::vector<float> single(data.size());
            for (std::size_t i = 0; i < data.size(); ++i) {
                single[i] = static_cast<float>(data[i]);
            }
            data = hicx::chic::smooth_interaction_values(std::span<const float>(single),
                                                         settings.average_contact_bin);
        } else {
            data = hicx::chic::smooth_interaction_values(std::span<const double>(data),
                                                         settings.average_contact_bin);
        }
        data_is_float32 = false;
    }

    const std::vector<double> raw = data;
    const std::vector<double> relative = hicx::chic::compute_relative_values(raw, denominator);
    std::vector<double> xfold = divide(relative, background);
    std::vector<double> pvalues =
        hicx::chic::p_values(model, raw, viewpoint.index_before_viewpoint, data_is_float32);

    const hicx::chic::ViewpointRange range_again = hicx::chic::calculate_viewpoint_range(
        matrix, point, settings.range_upstream, settings.range_downstream);
    return hicx::chic::create_interaction_file_data(
        matrix, point, point.chromosome, range_again.region_start, range_again.region_end,
        relative, raw, gene, denominator, std::move(pvalues), std::move(xfold));
}

struct Slot {
    InteractionFileData data;
    bool failed = false;
    std::string message;
};

}  // namespace

int main(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser(
        "chicViewpoint",
        "Computes per input matrix all viewpoints which are defined in the reference points file.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrices", "-m"})
        .required()
        .nargs("+")
        .input({"cool", "h5"})
        .help("Path to the Hi-C matrices which store the captured Hi-C data per sample.");
    required.add({"--range"})
        .required()
        .type("int")
        .nargs(2)
        .help("The region upstream and downstream of a reference point which should be "
              "considered in the analysis. Format is --range upstream downstream");
    required.add({"--referencePoints", "-rp"})
        .required()
        .input({"bed"})
        .help("Reference point file, 'chr 100' for a single reference point or 'chr 100 200' "
              "for a reference region, one per line.");
    required.add({"--backgroundModelFile", "-bmf"})
        .required()
        .input({"txt"})
        .help("path to the background file computed by chicViewpointBackgroundModel");
    required.add({"--outFileName", "-o"})
        .default_value("chic_files.hdf5")
        .output({"hdf5"})
        .help("This hdf5 file contains all created viewpoint files.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--threads", "-t"})
        .type("int")
        .default_value(4)
        .help("Number of threads.");
    optional.add({"--averageContactBin"})
        .type("int")
        .default_value(5)
        .help("Average the contacts of n bins via a sliding window approach to smooth the values "
              "and be less sensitive for outliers.");
    optional.add({"--fixateRange", "-fs"})
        .type("int")
        .default_value(500000)
        .help("Fixate range of background model starting at distance x.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    const cli::Namespace args = parser.parse(argc, argv);

    Settings settings;
    settings.range_upstream = args.integers("range").at(0);
    settings.range_downstream = args.integers("range").at(1);
    settings.fixate_range = args.integer("fixateRange");
    settings.average_contact_bin = args.integer("averageContactBin");
    const std::vector<std::string> matrices = args.strs("matrices");
    const std::int64_t threads = args.integer("threads");

    try {
        const hicx::chic::ReferencePoints reference_points =
            hicx::chic::read_reference_points(args.str("referencePoints"));
        const std::size_t count = reference_points.points.size();
        if (threads == 0) {
            throw std::runtime_error("ZeroDivisionError: integer division or modulo by zero");
        }

        const BackgroundModel model = hicx::chic::read_background_model(
            args.str("backgroundModelFile"), settings.range_upstream, settings.range_downstream,
            settings.fixate_range, false);
        const BackgroundModel means = hicx::chic::read_background_model(
            args.str("backgroundModelFile"), settings.range_upstream, settings.range_downstream,
            settings.fixate_range, true);

        hicx::chic::Hdf5Writer writer(args.str("outFileName"));
        writer.set_attribute("/", "type", std::string("interactions"));
        writer.set_attribute("/", "version", std::string(hicx::kVersion));
        const std::vector<std::int64_t> range{settings.range_upstream, settings.range_downstream};
        writer.set_attribute("/", "range", std::span<const std::int64_t>(range));
        writer.set_attribute("/", "averageContactBin", settings.average_contact_bin);
        writer.set_attribute("/", "fixateRange", settings.fixate_range);

        // matrix_collection: a dict keyed by the matrix path, so a repeated
        // path keeps its first position and its last result.
        std::vector<std::pair<std::string, std::vector<InteractionFileData>>> collection;
        std::int64_t resolution = 0;
        for (const std::string& matrix_path : matrices) {
            const ViewpointMatrix matrix = ViewpointMatrix::load(matrix_path);
            if (resolution == 0) {
                resolution = matrix.bin_size();
                writer.set_attribute("/", "resolution", resolution);
            }
            std::vector<Slot> slots(count);
            for (std::size_t index = 0; index < count; ++index) {
                Slot& slot = slots[index];
                try {
                    slot.data = compute_entry(matrix, reference_points.points[index],
                                              reference_points.genes[index], settings, model,
                                              means);
                } catch (const std::exception& error) {
                    slot.failed = true;
                    slot.message = error.what();
                }
            }
            std::vector<InteractionFileData> results;
            results.reserve(count);
            for (Slot& slot : slots) {
                if (slot.failed) {
                    std::fprintf(stderr, "chicViewpoint: %s\n", slot.message.c_str());
                    return 1;
                }
                results.push_back(std::move(slot.data));
            }
            bool replaced = false;
            for (auto& entry : collection) {
                if (entry.first == matrix_path) {
                    entry.second = std::move(results);
                    replaced = true;
                    break;
                }
            }
            if (!replaced) {
                collection.emplace_back(matrix_path, std::move(results));
            }
        }

        std::string chromosome_group;
        for (const auto& [matrix_path, results] : collection) {
            const std::string name = basename(matrix_path);
            const std::string matrix_group = name.substr(0, name.find('.'));
            writer.create_group(matrix_group);
            writer.create_group(matrix_group + "/genes");
            for (std::size_t i = 0; i < results.size(); ++i) {
                const InteractionFileData& data = results[i];
                const std::string candidate = matrix_group + "/" + data.chromosome;
                if (!writer.exists(candidate)) {
                    writer.create_group(candidate);
                    chromosome_group = candidate;
                } else if (chromosome_group.empty()) {
                    throw std::runtime_error(
                        "NameError: cannot access local variable 'chromosomeObject'");
                }
                const std::string group = chromosome_group + "/" + data.gene;
                if (writer.exists(group)) {
                    throw std::runtime_error(
                        "the gene name '" + data.gene + "' occurs twice in " + chromosome_group +
                        ". The Python reference does not terminate on this input "
                        "(Viewpoint.createUniqueHDFGroup never changes the name it retries); "
                        "the C++ port stops instead. Make the gene names in the reference "
                        "point file unique.");
                }
                writer.create_group(group);
                hicx::chic::write_interaction_datasets(
                    writer, group, data, hicx::chic::python_int(reference_points.points[i].start),
                    hicx::chic::python_int(reference_points.points[i].end));
                (void)writer.hard_link(group, matrix_group + "/genes/" + data.gene);
            }
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chicViewpoint: %s\n", error.what());
        return 1;
    }
    return 0;
}
