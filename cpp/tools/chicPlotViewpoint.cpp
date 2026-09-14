// Port of hicexplorer/chicPlotViewpoint.py (cpp/PLAN.md tier 7, option (a)).
//
// Every plot's data is read and computed here: the plot groups of the
// combination mode, the interactions, p-values and background model within
// the range (Viewpoint.getDataForPlotting), the differential regions
// (readRejectedFile) and the significant regions (readSignificantRegionsFile).
// plot/hicexplorer_plot/chicPlotViewpoint.py draws the figures with the
// reference's calls and writes the tar.gz.
//
// What is reproduced, pinned by the harness cases:
//
//  1. The plot groups (main, :333-422). dual pairs every sample with every
//     later sample, but lists the chromosomes and genes of the first sample
//     for both (matrix_obj2 is interactionFileHDF5Object[sample]); allGenes
//     follows the HDF5 order of the first sample's genes group; oneGene adds
//     every matrix's viewpoint on its own before the group of all matrices.
//  2. readBackgroundDataFile with pMean: the last column per relative
//     position, both ends clamped with the downstream range, and the model
//     extended by repeating its outermost values.
//  3. getDataForPlotting: the range cut; the viewpoint entry repeated
//     peak_width = |end - start| // resolution times when the reference point
//     is wider than a bin, and dropped when it is not; the viewpoint index
//     from the background keys when a model is given, else from the
//     interaction keys.
//  4. readInteractionFile: the arrays present are collected in a fixed order
//     and indexed by position, so a file lacking one of them, or any other
//     failure while building the records, yields no data.
//  5. readRejectedFile raises for a viewpoint without a reference point, which
//     fails the run; readSignificantRegionsFile is called for the i-th
//     plotted viewpoint with the i-th available entry of the group, inside a
//     try block that ignores its failures.
//  6. The worker chunks: len(groups) // --threads groups per worker, the rest
//     in the last. The archive pairs the i-th image with the i-th file name
//     after the chunks are flattened, and a group without an image still
//     adds its name, so later names shift. The drawing module reproduces that
//     from the chunks written here.
//
// A Python exception inside a worker makes the reference exit 1 without an
// archive. The port computes every plot before drawing, and exits 1 in the
// same situations.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/hdf5_util.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/version.hpp"

namespace {

namespace cli = hicx::cli;
namespace plot = hicx::plot;
namespace h5 = hicx::h5;

const char* const kDescription =
    "\n"
    "chicPlotViewpoint plots one or many viewpoints with the average background model and the "
    "computed p-value per sample. In addition, it can highlight differential interactions of "
    "two samples and/or significant regions.\n"
    "\n";

const char* const kUsage =
    "usage: chicPlotViewpoint --interactionFile INTERACTIONFILE --range RANGE RANGE\n"
    "                         [--backgroundModelFile BACKGROUNDMODELFILE]\n"
    "                         [--differentialTestResult DIFFERENTIALTESTRESULT]\n"
    "                         [--significantInteractions SIGNIFICANTINTERACTIONS]\n"
    "                         [--plotSignificantInteractions]\n"
    "                         [--outFileName OUTFILENAME]\n"
    "                         [--outputFormat OUTPUTFORMAT] [--dpi DPI]\n"
    "                         [--combinationMode {dual,single,allGenes,oneGene}]\n"
    "                         [--combinationName COMBINATIONNAME]\n"
    "                         [--colorMapPvalue COLORMAPPVALUE]\n"
    "                         [--maxPValue MAXPVALUE] [--minPValue MINPVALUE]\n"
    "                         [--pValue]\n"
    "                         [--pValueSignificanceLevels PVALUESIGNIFICANCELEVELS [PVALUESIGNIFICANCELEVELS ...]]\n"
    "                         [--xFold XFOLD] [--truncateZeroPvalues]\n"
    "                         [--colorList COLORLIST [COLORLIST ...]]\n"
    "                         [--threads THREADS] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "chicPlotViewpoint plots one or many viewpoints with the average background model and the "
    "computed p-value per sample. In addition, it can highlight differential interactions of "
    "two samples and/or significant regions.\n"
    "\n"
    "Required arguments:\n"
    "  --interactionFile INTERACTIONFILE, -if INTERACTIONFILE\n"
    "                        path to the interaction files which should be used for\n"
    "                        plotting\n"
    "  --range RANGE RANGE   Defines the region upstream and downstream of a\n"
    "                        reference point which should be included. Format is\n"
    "                        --region upstream downstream, e.g.: --region 500000\n"
    "                        500000 plots 500kb up- and 500kb downstream. This\n"
    "                        value should not exceed the range used in the other\n"
    "                        chic-tools.\n"
    "\n"
    "Optional arguments:\n"
    "  --backgroundModelFile BACKGROUNDMODELFILE, -bmf BACKGROUNDMODELFILE\n"
    "                        path to the background file which should be used for\n"
    "                        plotting\n"
    "  --differentialTestResult DIFFERENTIALTESTRESULT, -dif DIFFERENTIALTESTRESULT\n"
    "                        Path to the H0 rejected files to highlight the regions\n"
    "                        in the plot.\n"
    "  --significantInteractions SIGNIFICANTINTERACTIONS, -si SIGNIFICANTINTERACTIONS\n"
    "                        Path to the files with detected significant\n"
    "                        interactions to highlight the regions in the plot.\n"
    "  --plotSignificantInteractions, -psi\n"
    "                        Highlights the significant interactions in the plot\n"
    "                        itself. If not set, only the p-values are updated\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        Output tar.gz of the files (Default: plots.tar.gz).\n"
    "  --outputFormat OUTPUTFORMAT, -format OUTPUTFORMAT\n"
    "                        Output format of the plot (Default: png).\n"
    "  --dpi DPI             Optional parameter: Resolution for the image, ifoutput\n"
    "                        is a raster graphics image (e.g png, jpg) (Default:\n"
    "                        300).\n"
    "  --combinationMode {dual,single,allGenes,oneGene}, -cm {dual,single,allGenes,oneGene}\n"
    "                        This option defines how the interaction data should be\n"
    "                        computed and combined: dual: Combines as follows:\n"
    "                        [[matrix1_gene1, matrix2_gene1], [matrix2_gene1,\n"
    "                        matrix3_gene1],[matrix1_gene2, matrix2_gene2],\n"
    "                        ...]single: Combines as follows: [matrix1_gene1,\n"
    "                        matrix1_gene2, matrix2_gene1, ...], allGenes: Combines\n"
    "                        as follows: [[matrix1_gene1, matrix2_gene1,\n"
    "                        matrix2_gene1], [matrix1_gene2, matrix2_gene2,\n"
    "                        matrix3_gene2], ...]oneGene: Computes all data of one\n"
    "                        gene, please specify '--'. If a gene is not unique,\n"
    "                        each viewpoint is treated independently. (Default:\n"
    "                        dual).\n"
    "  --combinationName COMBINATIONNAME, -cn COMBINATIONNAME\n"
    "                        Gene name or file name for modes 'oneGene' or 'file'\n"
    "                        of parameter '--combinationMode' (Default: None).\n"
    "  --colorMapPvalue COLORMAPPVALUE\n"
    "                        Color map to use for the p-value. Available values can\n"
    "                        be seen here: http://matplotlib.org/examples/color/col\n"
    "                        ormaps_reference.html (Default: RdYlBu).\n"
    "  --maxPValue MAXPVALUE, -map MAXPVALUE\n"
    "                        Maximal value for p-value. Values above this threshold\n"
    "                        are set to this value (Default: 0.1).\n"
    "  --minPValue MINPVALUE, -mp MINPVALUE\n"
    "                        Minimal value for p-value. Values below this threshold\n"
    "                        are set to this value (Default: 0.0).\n"
    "  --pValue, -p          Plot p-values as a colorbar\n"
    "  --pValueSignificanceLevels PVALUESIGNIFICANCELEVELS [PVALUESIGNIFICANCELEVELS ...], -psl PVALUESIGNIFICANCELEVELS [PVALUESIGNIFICANCELEVELS ...]\n"
    "                        Highlight the p-values by the defined significance\n"
    "                        levels.\n"
    "  --xFold XFOLD, -xf XFOLD\n"
    "                        Plot x-fold region for the mean background.\n"
    "  --truncateZeroPvalues, -tzpv\n"
    "                        Sets all p-values which are equal to zero to one.\n"
    "  --colorList COLORLIST [COLORLIST ...], -cl COLORLIST [COLORLIST ...]\n"
    "                        Colorlist for the viewpoint lines (Default g b c m y\n"
    "                        k).\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads (uses the python multiprocessing\n"
    "                        module) (Default: 4).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n"
    "\n"
    "C++ port: the plot data is read and computed in C++, and the figures and the\n"
    "archive are written by the hicexplorer_plot drawing layer with the calls of the\n"
    "Python tool (HICX_PLOT_PYTHON names the interpreter). The C++-only option\n"
    "--plotData FILE writes the data of the figures as JSON to FILE instead of\n"
    "drawing them.\n";

// A Python exception that ends the reference with exit status 1.
class PythonError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

using Triplet = std::vector<std::string>;

std::string join(const std::vector<std::string>& parts, const std::string& separator) {
    std::string out;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        out += (i > 0 ? separator : "") + parts[i];
    }
    return out;
}

// h5py's `path in file`: every component must exist; a missing intermediate
// group is simply absent.
bool contains(const h5::File& file, const std::string& path) {
    std::string prefix;
    std::size_t begin = 0;
    while (begin <= path.size()) {
        std::size_t end = path.find('/', begin);
        if (end == std::string::npos) {
            end = path.size();
        }
        const std::string component = path.substr(begin, end - begin);
        if (!component.empty()) {
            prefix += "/" + component;
            if (!file.exists(prefix)) {
                return false;
            }
        }
        begin = end + 1;
    }
    return !prefix.empty();
}

std::vector<std::string> sorted_children(const h5::File& file, const std::string& group) {
    std::vector<std::string> names = file.children(group);
    std::sort(names.begin(), names.end());
    return names;
}

std::string string_attribute(const h5::File& file, const std::string& name) {
    const auto attributes = file.attributes("/");
    const auto it = attributes.find(name);
    if (it == attributes.end()) {
        throw PythonError("KeyError: \"Can't open attribute (can't locate attribute: '" + name +
                          "')\"");
    }
    if (const auto* text = std::get_if<std::string>(&it->second)) {
        return *text;
    }
    return "";
}

// One record of readInteractionFile's interaction_file_data.
struct Record {
    std::string chromosome;
    std::int64_t start = 0;
    std::int64_t end = 0;
    double relative = 0.0;
    double interaction = 0.0;
    double pvalue = 0.0;
};

// readInteractionFile: records in first insertion order of their relative
// position, a later duplicate replacing the value; the reference point.
struct InteractionData {
    std::vector<double> order;
    std::map<double, Record> records;
    std::vector<std::optional<double>> reference_point;  // empty when not found
};

InteractionData read_interaction_file(const h5::File& file, const Triplet& triplet) {
    InteractionData out;
    const std::string internal = join(triplet, "/");
    if (!contains(file, internal)) {
        return out;
    }
    const std::string base = "/" + internal + "/";
    std::vector<std::vector<double>> data;
    for (const char* name : {"relative_position_list", "interaction_data_list", "pvalue", "raw",
                             "xfold", "start_list", "end_list"}) {
        if (contains(file, internal + "/" + name)) {
            data.push_back(file.read_doubles(base + name));
        }
    }
    std::optional<std::string> chromosome;
    if (contains(file, internal + "/chromosome")) {
        chromosome = file.read_strings(base + "chromosome").at(0);
    }
    const bool has_gene = contains(file, internal + "/gene");
    const bool has_sum = contains(file, internal + "/sum_of_interactions");
    std::optional<double> reference_start;
    std::optional<double> reference_end;
    if (file.exists(base + "reference_point_start")) {
        reference_start = file.read_doubles(base + "reference_point_start").at(0);
    }
    if (file.exists(base + "reference_point_end")) {
        reference_end = file.read_doubles(base + "reference_point_end").at(0);
    }
    // The try block around `for i in range(len(data[0]))`: without any array
    // data[0] raises; with an empty first array the loop does not run and the
    // reference point is returned; otherwise a missing array, chromosome, gene
    // or sum, or a shorter array, raises in the loop, and every partial record
    // is discarded.
    if (data.empty()) {
        return InteractionData{};
    }
    if (data[0].empty()) {
        out.reference_point = {reference_start, reference_end};
        return out;
    }
    if (data.size() < 7 || !chromosome.has_value() || !has_gene || !has_sum) {
        return InteractionData{};
    }
    for (std::size_t i = 0; i < data[0].size(); ++i) {
        for (std::size_t column = 1; column < 7; ++column) {
            if (i >= data[column].size()) {
                return InteractionData{};
            }
        }
        Record record;
        record.chromosome = *chromosome;
        record.start = static_cast<std::int64_t>(data[5][i]);
        record.end = static_cast<std::int64_t>(data[6][i]);
        record.relative = data[0][i];
        record.interaction = data[1][i];
        record.pvalue = data[2][i];
        const double key = data[0][i];
        if (out.records.find(key) == out.records.end()) {
            out.order.push_back(key);
        }
        out.records[key] = record;
    }
    out.reference_point = {reference_start, reference_end};
    return out;
}

// readBackgroundDataFile(path, range, range[1], pMean=True), in insertion
// order.
std::vector<std::pair<std::int64_t, double>> read_background(const std::string& path,
                                                             std::int64_t upstream,
                                                             std::int64_t downstream) {
    std::ifstream in(path);
    if (!in) {
        throw PythonError("FileNotFoundError: [Errno 2] No such file or directory: '" + path + "'");
    }
    std::vector<std::pair<std::int64_t, double>> distance;
    auto assign = [&](std::int64_t key, double value) {
        for (auto& entry : distance) {
            if (entry.first == key) {
                entry.second = value;
                return;
            }
        }
        distance.emplace_back(key, value);
    };
    auto find = [&](std::int64_t key) -> double {
        for (const auto& entry : distance) {
            if (entry.first == key) {
                return entry.second;
            }
        }
        throw PythonError("KeyError: " + std::to_string(key) + " (readBackgroundDataFile)");
    };
    std::string line;
    std::getline(in, line);  // header
    while (std::getline(in, line)) {
        std::vector<std::string> fields;
        std::size_t begin = 0;
        while (true) {
            const std::size_t tab = line.find('\t', begin);
            fields.push_back(line.substr(begin, tab == std::string::npos ? std::string::npos
                                                                           : tab - begin));
            if (tab == std::string::npos) {
                break;
            }
            begin = tab + 1;
        }
        std::int64_t key = 0;
        double value = 0.0;
        if (!cli::python_int(fields.front(), &key)) {
            throw PythonError("ValueError: invalid literal for int() with base 10: '" +
                              fields.front() + "'");
        }
        if (!cli::python_float(fields.back(), &value)) {
            throw PythonError("ValueError: could not convert string to float: '" + fields.back() +
                              "'");
        }
        assign(key, value);
    }
    if (distance.empty()) {
        throw PythonError("ValueError: max() arg is an empty sequence (readBackgroundDataFile)");
    }
    std::int64_t max_key = distance.front().first;
    std::int64_t min_key = distance.front().first;
    for (const auto& entry : distance) {
        max_key = std::max(max_key, entry.first);
        min_key = std::min(min_key, entry.first);
    }
    if (max_key > downstream) {
        max_key = downstream;
    }
    if (min_key < -downstream) {
        min_key = -downstream;
    }
    if (distance.size() < 2) {
        throw PythonError("IndexError: list index out of range (readBackgroundDataFile)");
    }
    const std::int64_t inc =
        std::llabs(std::llabs(distance[0].first) - std::llabs(distance[1].first));
    if (max_key < downstream) {
        if (inc == 0) {
            throw PythonError("readBackgroundDataFile would not terminate: the first two "
                              "positions have the same distance to the viewpoint");
        }
        const double value = find(max_key);
        for (std::int64_t i = max_key; i < downstream;) {
            i += inc;
            assign(i, value);
        }
    }
    if (min_key > -upstream) {
        if (inc == 0) {
            throw PythonError("readBackgroundDataFile would not terminate: the first two "
                              "positions have the same distance to the viewpoint");
        }
        const double value = find(min_key);
        for (std::int64_t i = min_key; i > -upstream;) {
            i -= inc;
            assign(i, value);
        }
    }
    return distance;
}

struct PlotData {
    bool skip = false;
    std::vector<double> data;
    std::optional<std::vector<double>> background;
    std::vector<double> p_values;
    std::int64_t viewpoint_index_start = 0;
    std::int64_t viewpoint_index_end = 0;
    std::vector<std::optional<double>> viewpoint;
};

std::int64_t peak_width(const Record& record, std::int64_t resolution, bool* wide) {
    const std::int64_t width = std::llabs(record.start - record.end);
    *wide = width > resolution;
    return width / resolution;
}

// Viewpoint.getDataForPlotting.
PlotData data_for_plotting(const InteractionData& interactions, std::int64_t upstream,
                           std::int64_t downstream,
                           const std::optional<std::vector<std::pair<std::int64_t, double>>>& model,
                           std::int64_t resolution) {
    PlotData out;
    out.viewpoint = interactions.reference_point;
    auto in_range = [&](double key) { return key >= -upstream && key <= downstream; };
    std::map<double, const Record*> kept;
    for (const auto& [key, record] : interactions.records) {
        if (in_range(key)) {
            kept[key] = &record;
        }
    }
    std::vector<std::pair<std::int64_t, double>> background;
    if (model.has_value()) {
        for (const auto& entry : *model) {
            if (in_range(static_cast<double>(entry.first))) {
                background.push_back(entry);
            }
        }
        std::sort(background.begin(), background.end());
    }
    std::optional<std::int64_t> index_end;
    auto append_key = [&](double key, const Record& record, bool with_background, double bg) {
        if (key == 0) {
            bool wide = false;
            const std::int64_t width = peak_width(record, resolution, &wide);
            if (wide) {
                index_end = width;
                for (std::int64_t i = 0; i < width; ++i) {
                    out.data.push_back(record.interaction);
                }
                for (std::int64_t i = 0; i < width; ++i) {
                    out.p_values.push_back(record.pvalue);
                }
                if (with_background) {
                    for (std::int64_t i = 0; i < width; ++i) {
                        out.background->push_back(bg);
                    }
                }
            }
        } else {
            out.data.push_back(record.interaction);
            out.p_values.push_back(record.pvalue);
            if (with_background) {
                out.background->push_back(bg);
            }
        }
    };
    if (!background.empty()) {
        const auto zero = std::find_if(background.begin(), background.end(),
                                       [](const auto& entry) { return entry.first == 0; });
        if (zero == background.end()) {
            throw PythonError("ValueError: 0 is not in list (the background model has no "
                              "viewpoint position)");
        }
        out.viewpoint_index_start = zero - background.begin();
        out.background.emplace();
        for (const auto& [key, value] : background) {
            const auto it = kept.find(static_cast<double>(key));
            if (it != kept.end()) {
                append_key(static_cast<double>(key), *it->second, true, value);
            }
        }
    } else {
        for (const auto& [key, record] : kept) {
            append_key(key, *record, false, 0.0);
        }
        const auto zero = kept.find(0.0);
        if (zero == kept.end()) {
            throw PythonError("ValueError: 0 is not in list (no interaction at the viewpoint "
                              "within the range)");
        }
        out.viewpoint_index_start = std::distance(kept.begin(), zero);
    }
    out.viewpoint_index_end = index_end.has_value() ? *index_end + out.viewpoint_index_start
                                                    : out.viewpoint_index_start;
    return out;
}

std::int64_t python_int_of(const std::optional<double>& value) {
    if (!value.has_value()) {
        throw PythonError("TypeError: int() argument must be a string, a bytes-like object or a "
                          "real number, not 'NoneType'");
    }
    return static_cast<std::int64_t>(*value);
}

using Regions = std::vector<std::pair<double, double>>;

// Viewpoint.readRejectedFile.
std::optional<Regions> rejected_regions(const h5::File& file, const std::vector<std::string>& group,
                                        const PlotData& plot, std::int64_t resolution,
                                        std::int64_t upstream, std::int64_t downstream) {
    if (group.empty()) {
        return std::nullopt;
    }
    const std::string internal = join(group, "/");
    if (!contains(file, internal)) {
        return std::nullopt;
    }
    const std::string base = "/" + internal + "/";
    std::vector<double> starts;
    std::vector<double> ends;
    std::vector<double> relative;
    try {
        starts = file.read_doubles(base + "start_list");
        ends = file.read_doubles(base + "end_list");
        relative = file.read_doubles(base + "relative_distance_list");
    } catch (const std::exception&) {
        return std::nullopt;
    }
    if (starts.empty()) {
        return std::nullopt;
    }
    if (plot.viewpoint.size() != 2) {
        throw PythonError("ValueError: not enough values to unpack (expected 2, got " +
                          std::to_string(plot.viewpoint.size()) + ") (readRejectedFile)");
    }
    Regions areas;
    const std::size_t n = std::min({starts.size(), ends.size(), relative.size()});
    for (std::size_t i = 0; i < n; ++i) {
        const auto rel = static_cast<std::int64_t>(relative[i]);
        if (rel >= -upstream && rel <= downstream) {
            const double width = (ends[i] - starts[i]) / static_cast<double>(resolution);
            double genomic = 0.0;
            std::int64_t index = 0;
            if (rel < 0) {
                genomic = starts[i] - static_cast<double>(python_int_of(plot.viewpoint[0]));
                index = plot.viewpoint_index_start;
            } else {
                genomic = starts[i] - static_cast<double>(python_int_of(plot.viewpoint[1]));
                index = plot.viewpoint_index_end;
            }
            const double position = static_cast<double>(index) + genomic / static_cast<double>(resolution);
            areas.emplace_back(position, position + width);
        }
    }
    if (areas.empty()) {
        return std::nullopt;
    }
    return areas;
}

struct Significant {
    Regions regions;
    std::vector<std::tuple<std::int64_t, std::int64_t, double>> p_values;
};

// Viewpoint.readSignificantRegionsFile; nullopt for (None, None).
std::optional<Significant> significant_regions(const h5::File& file, const Triplet& triplet,
                                               const PlotData& plot, std::int64_t resolution,
                                               std::int64_t upstream, std::int64_t downstream) {
    const InteractionData data = read_interaction_file(file, triplet);
    if (plot.viewpoint.size() != 2) {
        return std::nullopt;
    }
    Significant out;
    for (const double key : data.order) {
        const Record& line = data.records.at(key);
        const auto rel = static_cast<std::int64_t>(line.relative);
        if (rel >= -upstream && rel <= downstream) {
            const double width = static_cast<double>(line.end - line.start) / static_cast<double>(resolution);
            double genomic = 0.0;
            std::int64_t index = 0;
            if (rel < 0) {
                genomic = static_cast<double>(line.start - python_int_of(plot.viewpoint[0]));
                index = plot.viewpoint_index_start;
            } else {
                genomic = static_cast<double>(line.start - python_int_of(plot.viewpoint[1]));
                index = plot.viewpoint_index_end;
            }
            const double position = static_cast<double>(index) + genomic / static_cast<double>(resolution);
            out.regions.emplace_back(position, position + width);
            out.p_values.emplace_back(static_cast<std::int64_t>(position),
                                      static_cast<std::int64_t>(position + width), line.pvalue);
        }
    }
    if (out.regions.empty()) {
        return std::nullopt;
    }
    return out;
}

std::string regions_json(const std::optional<Regions>& regions) {
    if (!regions.has_value()) {
        return "null";
    }
    std::vector<std::string> items;
    for (const auto& [first, second] : *regions) {
        items.push_back(plot::json_numbers({first, second}));
    }
    return plot::json_list(items);
}

}  // namespace

int main(int argc, char** argv) {
    cli::Parser parser("chicPlotViewpoint", kDescription);
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--interactionFile", "-if"})
        .required()
        .input({"hdf5"})
        .help("path to the interaction files which should be used for plotting");
    required.add({"--range"})
        .required()
        .type("int")
        .nargs(2)
        .help("Defines the region upstream and downstream of a reference point which should be "
              "included. Format is --region upstream downstream, e.g.: --region 500000 500000 "
              "plots 500kb up- and 500kb downstream. This value should not exceed the range used "
              "in the other chic-tools.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--backgroundModelFile", "-bmf"})
        .required(false)
        .input({"txt"})
        .help("path to the background file which should be used for plotting");
    optional.add({"--differentialTestResult", "-dif"})
        .required(false)
        .input({"hdf5"})
        .help("Path to the H0 rejected files to highlight the regions in the plot.");
    optional.add({"--significantInteractions", "-si"})
        .required(false)
        .input({"hdf5"})
        .help("Path to the files with detected significant interactions to highlight the regions "
              "in the plot.");
    optional.add({"--plotSignificantInteractions", "-psi"})
        .required(false)
        .action(cli::Action::StoreTrue)
        .help("Highlights the significant interactions in the plot itself. If not set, only the "
              "p-values are updated");
    optional.add({"--outFileName", "-o"})
        .required(false)
        .default_value("plots.tar.gz")
        .output({"tar.gz"})
        .help("Output tar.gz of the files (Default: %(default)s).");
    optional.add({"--outputFormat", "-format"})
        .required(false)
        .default_value("png")
        .help("Output format of the plot (Default: %(default)s).");
    optional.add({"--dpi"})
        .type("int")
        .default_value(std::int64_t{300})
        .required(false)
        .help("Optional parameter: Resolution for the image, ifoutput is a raster graphics image "
              "(e.g png, jpg) (Default: %(default)s).");
    optional.add({"--combinationMode", "-cm"})
        .default_value("dual")
        .choices({"dual", "single", "allGenes", "oneGene"})
        .help("This option defines how the interaction data should be computed and combined.");
    optional.add({"--combinationName", "-cn"})
        .required(false)
        .help("Gene name or file name for modes 'oneGene' or 'file' of parameter "
              "'--combinationMode' (Default: %(default)s).");
    optional.add({"--colorMapPvalue"})
        .default_value("RdYlBu")
        .help("Color map to use for the p-value (Default: %(default)s).");
    optional.add({"--maxPValue", "-map"})
        .type("float")
        .default_value(0.1)
        .help("Maximal value for p-value. Values above this threshold are set to this value "
              "(Default: %(default)s).");
    optional.add({"--minPValue", "-mp"})
        .type("float")
        .default_value(0.0)
        .help("Minimal value for p-value. Values below this threshold are set to this value "
              "(Default: %(default)s).");
    optional.add({"--pValue", "-p"}).action(cli::Action::StoreTrue).help("Plot p-values as a colorbar");
    optional.add({"--pValueSignificanceLevels", "-psl"})
        .type("float")
        .nargs("+")
        .help("Highlight the p-values by the defined significance levels.");
    optional.add({"--xFold", "-xf"}).type("float").help("Plot x-fold region for the mean background.");
    optional.add({"--truncateZeroPvalues", "-tzpv"})
        .required(false)
        .action(cli::Action::StoreTrue)
        .help("Sets all p-values which are equal to zero to one.");
    optional.add({"--colorList", "-cl"})
        .required(false)
        .default_value(hicx::json::Value::array({hicx::json::Value::string("g"),
                                                 hicx::json::Value::string("b"),
                                                 hicx::json::Value::string("c"),
                                                 hicx::json::Value::string("m"),
                                                 hicx::json::Value::string("y"),
                                                 hicx::json::Value::string("k")}))
        .type("str")
        .nargs("+")
        .help("Colorlist for the viewpoint lines (Default g b c m y k).");
    optional.add({"--threads", "-t"})
        .required(false)
        .default_value(std::int64_t{4})
        .type("int")
        .help("Number of threads (uses the python multiprocessing module) (Default: "
              "%(default)s).");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    optional.add({"--plotData"})
        .metavar("FILE")
        .output({"json"})
        .cpp_only("The data of the figures as JSON, without drawing them (cpp/PLAN.md tier 7).")
        .help("Write the data the figures are drawn from as JSON to this file and do not draw "
              "them.");

    const cli::Namespace ns = parser.parse(argc, argv);
    if (const int refused = hicx::plot::preflight("chicPlotViewpoint", !ns.given("plotData")); refused != 0) {
        return refused;
    }
    const std::string interaction_path = ns.str("interactionFile");
    const std::vector<std::int64_t> range = ns.integers("range");
    const std::int64_t upstream = range[0];
    const std::int64_t downstream = range[1];
    const std::optional<std::string> background_path = ns.opt_str("backgroundModelFile");
    const std::optional<std::string> differential_path = ns.opt_str("differentialTestResult");
    const std::optional<std::string> significant_path = ns.opt_str("significantInteractions");
    const bool plot_significant = ns.flag("plotSignificantInteractions");
    const std::string mode = ns.str("combinationMode");
    const std::optional<std::string> combination_name = ns.opt_str("combinationName");
    const std::int64_t threads = ns.integer("threads");
    std::optional<std::vector<double>> levels;
    if (ns.given("pValueSignificanceLevels")) {
        levels = ns.reals("pValueSignificanceLevels");
    }

    std::string chunks_json;
    try {
        if (levels.has_value()) {
            double old = -100;
            for (const double element : *levels) {
                if (old < element) {
                    old = element;
                    continue;
                }
                std::fputs("--pValueSignificanceLevels levels need to increase\n", stderr);
                return 1;
            }
        }
        std::optional<std::vector<std::pair<std::int64_t, double>>> background;
        if (background_path.has_value() && !background_path->empty()) {
            background = read_background(*background_path, upstream, downstream);
        }

        const h5::File interactions(interaction_path);
        if (string_attribute(interactions, "type") != "interactions") {
            std::fputs("Please provide a file created by chicViewpoint for the parameter "
                       "--interactionFile.\n",
                       stderr);
            return 1;
        }
        const auto root_attributes = interactions.attributes("/");
        const auto resolution_it = root_attributes.find("resolution");
        if (resolution_it == root_attributes.end()) {
            throw PythonError("KeyError: \"Can't open attribute (can't locate attribute: "
                              "'resolution')\"");
        }
        const std::int64_t resolution =
            std::holds_alternative<std::int64_t>(resolution_it->second)
                ? std::get<std::int64_t>(resolution_it->second)
                : static_cast<std::int64_t>(std::get<double>(resolution_it->second));
        const std::vector<std::string> keys = interactions.children("/");

        if (differential_path.has_value() && !differential_path->empty() && mode != "dual") {
            std::fputs("Cannot use differential data, only possible for two samples in one plot.\n",
                       stderr);
            return 1;
        }

        auto without_genes = [&](const std::string& sample) {
            std::vector<std::string> names = sorted_children(interactions, "/" + sample);
            const auto it = std::find(names.begin(), names.end(), "genes");
            if (it == names.end()) {
                throw PythonError("ValueError: list.remove(x): x not in list ('genes')");
            }
            names.erase(it);
            return names;
        };

        std::vector<std::vector<Triplet>> groups;
        std::vector<std::vector<std::string>> differential_groups;
        if (mode == "dual") {
            if (keys.size() > 1) {
                for (std::size_t i = 0; i < keys.size(); ++i) {
                    for (std::size_t k = i + 1; k < keys.size(); ++k) {
                        const std::vector<std::string> chromosomes = without_genes(keys[i]);
                        for (const std::string& chromosome : chromosomes) {
                            const std::vector<std::string> genes =
                                sorted_children(interactions, "/" + keys[i] + "/" + chromosome);
                            for (const std::string& gene : genes) {
                                groups.push_back({{keys[i], chromosome, gene},
                                                  {keys[k], chromosome, gene}});
                            }
                        }
                    }
                }
                if (differential_path.has_value() && !differential_path->empty()) {
                    const h5::File differential(*differential_path);
                    if (string_attribute(differential, "type") != "differential") {
                        std::fputs("Please provide a file created by chicDifferentialTest for the "
                                   "parameter --differentialTestResult.\n",
                                   stderr);
                        return 1;
                    }
                    const std::vector<std::string> differential_keys = differential.children("/");
                    for (const auto& group : groups) {
                        std::vector<std::string> entry;
                        if (std::find(differential_keys.begin(), differential_keys.end(),
                                      group[0][0]) != differential_keys.end() &&
                            contains(differential, group[0][0] + "/" + group[1][0]) &&
                            contains(differential,
                                     group[0][0] + "/" + group[1][0] + "/" + group[1][1]) &&
                            contains(differential, group[0][0] + "/" + group[1][0] + "/" +
                                                       group[1][1] + "/" + group[1][2])) {
                            entry = {group[0][0], group[1][0], group[1][1], group[1][2],
                                     "rejected"};
                        }
                        differential_groups.push_back(entry);
                    }
                }
            } else {
                std::fputs("Dual mode selected but only one matrix is stored\n", stderr);
            }
        } else if (mode == "allGenes") {
            if (!keys.empty()) {
                if (!contains(interactions, keys[0] + "/genes")) {
                    throw PythonError("KeyError: \"Unable to open object (object 'genes' doesn't "
                                      "exist)\"");
                }
                for (const std::string& gene : interactions.children("/" + keys[0] + "/genes")) {
                    std::vector<Triplet> group;
                    for (const std::string& matrix : keys) {
                        group.push_back({matrix, "genes", gene});
                    }
                    groups.push_back(group);
                }
            }
        } else if (mode == "single") {
            for (const std::string& sample : keys) {
                for (const std::string& chromosome : without_genes(sample)) {
                    for (const std::string& gene :
                         sorted_children(interactions, "/" + sample + "/" + chromosome)) {
                        groups.push_back({{sample, chromosome, gene}});
                    }
                }
            }
        } else if (mode == "oneGene") {
            if (!keys.empty()) {
                if (!combination_name.has_value()) {
                    throw PythonError("TypeError: a bytes-like object or str is required, not "
                                      "'NoneType' (oneGene without --combinationName)");
                }
                if (!contains(interactions, keys[0] + "/genes")) {
                    throw PythonError("KeyError: \"Unable to open object (object 'genes' doesn't "
                                      "exist)\"");
                }
                std::vector<std::string> genes;
                for (std::int64_t counter = 0;; ++counter) {
                    const std::string name = counter == 0
                                                 ? *combination_name
                                                 : *combination_name + "_" + std::to_string(counter);
                    if (!contains(interactions, keys[0] + "/genes/" + name)) {
                        break;
                    }
                    genes.push_back(name);
                }
                for (const std::string& gene : genes) {
                    std::vector<Triplet> group;
                    for (const std::string& matrix : keys) {
                        group.push_back({matrix, "genes", gene});
                        groups.push_back({{matrix, "genes", gene}});
                    }
                    groups.push_back(group);
                }
            }
        }

        std::vector<std::vector<Triplet>> significant_groups;
        std::optional<h5::File> significant;
        if (significant_path.has_value() && !significant_path->empty()) {
            significant.emplace(*significant_path);
            const std::string type = string_attribute(*significant, "type");
            if (type != "significant") {
                std::fprintf(stderr,
                             "Please provide a file created by chicSignificantInteractions for "
                             "the parameter --significantInteractions. File type is %s\n",
                             type.c_str());
                return 1;
            }
            const std::vector<std::string> significant_keys = significant->children("/");
            for (const auto& group : groups) {
                std::vector<Triplet> available;
                for (const Triplet& item : group) {
                    if (std::find(significant_keys.begin(), significant_keys.end(), item[0]) !=
                            significant_keys.end() &&
                        contains(*significant, item[0] + "/" + item[1]) &&
                        contains(*significant, item[0] + "/" + item[1] + "/" + item[2])) {
                        available.push_back(item);
                    }
                }
                significant_groups.push_back(available);
            }
        }

        if (threads == 0) {
            throw PythonError("ZeroDivisionError: integer division or modulo by zero "
                              "(len(interactionFileList) // args.threads)");
        }
        const auto count = static_cast<std::int64_t>(groups.size());
        auto floor_div = [](std::int64_t a, std::int64_t b) {
            std::int64_t q = a / b;
            if ((a % b != 0) && ((a < 0) != (b < 0))) {
                --q;
            }
            return q;
        };
        const std::int64_t per_thread = floor_div(count, threads);
        auto slice = [&](std::int64_t first, std::optional<std::int64_t> last) {
            std::int64_t begin = std::clamp<std::int64_t>(first, 0, count);
            std::int64_t end = last.has_value() ? std::clamp<std::int64_t>(*last, 0, count) : count;
            return std::make_pair(begin, std::max(begin, end));
        };

        std::vector<std::string> chunk_texts;
        for (std::int64_t t = 0; t < threads; ++t) {
            const auto [begin, end] = t < threads - 1
                                          ? slice(t * per_thread, (t + 1) * per_thread)
                                          : slice(t * per_thread, std::nullopt);
            std::vector<std::string> group_texts;
            for (std::int64_t j = begin; j < end; ++j) {
                const auto& group = groups[static_cast<std::size_t>(j)];
                std::vector<std::string> item_texts;
                std::vector<std::string> file_name;
                for (std::size_t i = 0; i < group.size(); ++i) {
                    const Triplet& item = group[i];
                    file_name.push_back(item[0]);
                    const InteractionData data = read_interaction_file(interactions, item);
                    PlotData plot_data =
                        data_for_plotting(data, upstream, downstream, background, resolution);
                    plot::JsonObject item_json;
                    item_json.add("label", plot::json_string(join(item, ":")));
                    if (plot_data.data.size() <= 1 || plot_data.p_values.size() <= 1) {
                        std::fprintf(stderr,
                                     "Only one data point in given range, no plot is created! "
                                     "Interaction file %s Range [%lld, %lld]\n",
                                     join(item, ":").c_str(), static_cast<long long>(upstream),
                                     static_cast<long long>(downstream));
                        item_json.add("skip", "true");
                        item_texts.push_back(item_json.str());
                        continue;
                    }
                    std::optional<Regions> highlight;
                    if (differential_path.has_value() && !differential_path->empty()) {
                        const h5::File differential(*differential_path);
                        const std::vector<std::string>& entry =
                            differential_groups.at(static_cast<std::size_t>(j));
                        highlight = rejected_regions(differential, entry, plot_data, resolution,
                                                     upstream, downstream);
                    }
                    std::optional<Significant> found;
                    if (significant.has_value()) {
                        const auto& available = significant_groups.at(static_cast<std::size_t>(j));
                        if (i < available.size()) {
                            try {
                                found = significant_regions(*significant, available[i], plot_data,
                                                            resolution, upstream, downstream);
                            } catch (const std::exception&) {
                                found.reset();
                            }
                        }
                    }
                    std::optional<Regions> regions;
                    if (found.has_value() && plot_significant) {
                        regions = found->regions;
                    }
                    std::string p_texts = "null";
                    if (found.has_value()) {
                        std::vector<std::string> items;
                        for (const auto& [first, last, value] : found->p_values) {
                            items.push_back(plot::json_list(
                                {plot::json_int(first), plot::json_int(last), plot::json_number(value)}));
                        }
                        p_texts = plot::json_list(items);
                    }
                    item_json.add("skip", "false");
                    item_json.add("data", plot::json_numbers(plot_data.data));
                    item_json.add("background", plot_data.background.has_value()
                                                    ? plot::json_numbers(*plot_data.background)
                                                    : "null");
                    item_json.add("p_values", plot::json_numbers(plot_data.p_values));
                    item_json.add("highlight", regions_json(highlight));
                    item_json.add("significant_regions", regions_json(regions));
                    item_json.add("significant_p_values", p_texts);
                    item_texts.push_back(item_json.str());
                }
                file_name.push_back(group[0][2]);
                plot::JsonObject group_json;
                group_json.add("file_name", plot::json_string(join(file_name, "_")));
                group_json.add("items", plot::json_list(item_texts));
                group_texts.push_back(group_json.str());
            }
            chunk_texts.push_back(plot::json_list(group_texts));
        }
        chunks_json = plot::json_list(chunk_texts);

        plot::JsonObject data;
        data.add("outFileName", plot::json_string(ns.str("outFileName")));
        data.add("outputFormat", plot::json_string(ns.str("outputFormat")));
        data.add("dpi", plot::json_int(ns.integer("dpi")));
        data.add("range", plot::json_ints(range));
        data.add("resolution", plot::json_int(resolution));
        data.add("colorList", plot::json_strings(ns.strs("colorList")));
        const std::optional<double> x_fold = ns.opt_real("xFold");
        data.add("xFold", x_fold.has_value() ? plot::json_number(*x_fold) : "null");
        data.add("pValue", plot::json_bool(ns.flag("pValue")));
        data.add("truncateZeroPvalues", plot::json_bool(ns.flag("truncateZeroPvalues")));
        data.add("minPValue", plot::json_number(ns.real("minPValue")));
        data.add("maxPValue", plot::json_number(ns.real("maxPValue")));
        data.add("colorMapPvalue", plot::json_string(ns.str("colorMapPvalue")));
        data.add("pValueSignificanceLevels",
                 levels.has_value() ? plot::json_numbers(*levels) : "null");
        data.add("chunks", chunks_json);

        hicx::report_resource_usage("chicPlotViewpoint");
        return plot::draw("chicPlotViewpoint", data.str(), ns.opt_str("plotData"));
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chicPlotViewpoint: %s\n", error.what());
        return 1;
    }
}
