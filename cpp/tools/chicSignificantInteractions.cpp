// Port of hicexplorer/chicSignificantInteractions.py.
//
// Per sample viewpoint of a chicViewpoint interaction file, the candidate
// positions are preselected by x-fold over the background (--xFoldBackground)
// or by a loose p-value (--loosePValue), neighbouring candidates within one bin
// are merged, and a merged candidate is accepted when the negative binomial
// p-value of its summed raw interactions is at most --pValue and its raw sum
// reaches --peakInteractionsThreshold. Without a preselection every position is
// accepted on its stored p-value. The accepted candidates go to the significant
// file; per combination of samples their intervals, sorted and merged within
// one bin, go to the target file.
//
// Behaviour reproduced as the Python has it:
//
//   * Dual mode pairs each sample with every later one but takes the second
//     sample's chromosome and gene names from the first sample's tree
//     (chicSignificantInteractions.py:616-617 open the same group twice).
//   * The significant file is written in the order of the unique sorted
//     sample triplets (np.unique), while the reference point of the i-th
//     written viewpoint is taken from the i-th viewpoint in computation order,
//     so reference_point_start and reference_point_end can belong to another
//     viewpoint (on the test data Sox17 carries Eya1's reference point).
//   * Without a preselection, --peakInteractionsThreshold is compared with
//     the x-fold, not with the raw interactions (filter_by_pvalue reads the
//     last element of interaction_data, which is the x-fold).
//   * Viewpoint.merge_neighbors never keeps the last candidate when it is not
//     merged with its predecessor, and keeps nothing of a single candidate.
//   * A threshold file for --pValue, --loosePValue or --xFoldBackground is
//     used for the computation, and writing it as an HDF5 attribute then
//     raises TypeError: the run exits 1 with a significant file holding the
//     attributes written before it and no target file.
//   * --threads 0 raises ZeroDivisionError; a negative --threads starts no
//     worker, and both files are written with their attributes only.
//   * The target file's root attribute for the preselection value is spelled
//     mode_preselection_calue when there is no preselection.
//   * errorLog.txt in the working directory gains a line per sample viewpoint
//     without an accepted candidate.
//
// Where the Python does not terminate, the port stops with an error instead:
// Viewpoint.createUniqueHDFGroup retries a taken group name forever. The names
// are made unique per matrix before it is called, so this needs a gene group
// name that collides with a renamed one (for example genes "Eya1" and
// "Eya1_1" on one chromosome of the same matrix). A deviation, recorded as such.
//
// --correctForMultipleTesting {none,fdr,bonferroni} (C++ only, PLAN.md 9.7
// work item 3). The family is every p-value the tool compares against
// --pValue: the new p-value of each merged candidate whose position is in the
// background model, or, without a preselection, the stored p-value of each
// position, over every distinct sample viewpoint of the run in computation
// order. They are adjusted together (NaN left NaN and not counted), and a
// candidate is accepted when its adjusted p-value is at most the threshold
// with the same peak condition. The stored p-value stays the unadjusted one;
// each viewpoint group of the significant file gains pvalue_adjusted, and
// both files the root attribute correctForMultipleTesting. `none` writes the
// Python output. The reference is cpp/scripts/py_chicSignificantInteractions_calibrated.py.
//
// Threading: the work per viewpoint is a few hundred values, and the Python's
// output does not depend on --threads, because it collects its workers in
// order. As for chicViewpoint (cpp/OPTIMIZATION.md 6) the loop is sequential;
// --threads is accepted with the Python's semantics for 0 and negative values.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/chic_hdf5.hpp"
#include "hicx/chic_viewpoint.hpp"
#include "hicx/hdf5_util.hpp"
#include "hicx/scipy_special.hpp"
#include "hicx/stats_ops.hpp"
#include "hicx/version.hpp"

namespace {

namespace chic = hicx::chic;
namespace h5 = hicx::h5;
namespace cli = hicx::cli;

using chic::InteractionRecord;
using chic::InteractionTable;
using Triplet = std::vector<std::string>;

const char* const kUsage =
    "usage: chicSignificantInteractions --interactionFile INTERACTIONFILE --pValue\n"
    "                                   PVALUE\n"
    "                                   [--xFoldBackground XFOLDBACKGROUND | --loosePValue LOOSEPVALUE]\n"
    "                                   --backgroundModelFile BACKGROUNDMODELFILE\n"
    "                                   --range RANGE RANGE\n"
    "                                   [--outFileNameSignificant OUTFILENAMESIGNIFICANT]\n"
    "                                   [--outFileNameTarget OUTFILENAMETARGET]\n"
    "                                   [--combinationMode {dual,single}]\n"
    "                                   [--threads THREADS] [--truncateZeroPvalues]\n"
    "                                   [--fixateRange FIXATERANGE]\n"
    "                                   [--peakInteractionsThreshold PEAKINTERACTIONSTHRESHOLD]\n"
    "                                   [--correctForMultipleTesting {none,fdr,bonferroni}]\n"
    "                                   [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Per viewpoint the significant interactions are detected based on the background model. "
    "For each viewpoint file, an output file is created with all recorded significant "
    "interactions and\n"
    "a target file. The target file is especially useful in the batch mode context; for two "
    "consecutive listed control and treatment viewpoints it merges the significant "
    "interactions which can then be used\n"
    "to test for a differential interaction scheme.\n"
    "\n"
    "chicSignificantInteractions supports two modes to detect significant interactions, "
    "either by an x-fold over the average background or a loose p-value. In both cases "
    "neighboring significant peaks are merged together and an additional\n"
    "p-value is computed based on the sum of interactions for this neighborhood. Only "
    "interactions with a higher p-value (as specified by the threshold `--pValue`) are "
    "accepted as a significant interaction.\n"
    "\n"
    "options:\n"
    "  --xFoldBackground XFOLDBACKGROUND, -xf XFOLDBACKGROUND\n"
    "                        Filter x-fold over background. Used to merge\n"
    "                        neighboring bins with a broader peak but less\n"
    "                        significant interactions to a single peak with high\n"
    "                        significance. Used only for pValue option.\n"
    "  --loosePValue LOOSEPVALUE, -lp LOOSEPVALUE\n"
    "                        loose p-value threshold to filter target regions in a\n"
    "                        first round. Used to merge neighboring bins with a\n"
    "                        broader peak but less significant interactions to a\n"
    "                        single peak with high significance. Used only for\n"
    "                        pValue option.\n"
    "\n"
    "Required arguments:\n"
    "  --interactionFile INTERACTIONFILE, -if INTERACTIONFILE\n"
    "                        path to the interaction file (HDF5) which should be\n"
    "                        used for aggregation of the statistics.\n"
    "  --pValue PVALUE, -p PVALUE\n"
    "                        p-value threshold to filter target regions for\n"
    "                        inclusion in differential analysis.\n"
    "  --backgroundModelFile BACKGROUNDMODELFILE, -bmf BACKGROUNDMODELFILE\n"
    "                        path to the background file.\n"
    "  --range RANGE RANGE   Defines the region upstream and downstream of a\n"
    "                        reference point which should be included. Format is\n"
    "                        --region upstream downstream, e.g. --region 500000\n"
    "                        500000 plots 500kb up- and 500kb downstream. This\n"
    "                        value should not exceed the range used in the other\n"
    "                        chic-tools.\n"
    "\n"
    "Optional arguments:\n"
    "  --outFileNameSignificant OUTFILENAMESIGNIFICANT, -os OUTFILENAMESIGNIFICANT\n"
    "                        File name suffix to save the results; prefix is the\n"
    "                        input file name (Default:\n"
    "                        significantInteractions.hdf5).\n"
    "  --outFileNameTarget OUTFILENAMETARGET, -ot OUTFILENAMETARGET\n"
    "                        The file to store the target data (Default:\n"
    "                        targetFile.hdf5).\n"
    "  --combinationMode {dual,single}, -cm {dual,single}\n"
    "                        This option defines how the interaction data should be\n"
    "                        computed and combined: dual: Combines as follows:\n"
    "                        [[matrix1_gene1, matrix2_gene1], [matrix2_gene1,\n"
    "                        matrix3_gene1],[matrix1_gene2, matrix2_gene2],\n"
    "                        ...]single: Combines as follows: [matrix1_gene1,\n"
    "                        matrix1_gene2, matrix2_gene1, ...], (Default: dual).\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads (uses the python multiprocessing\n"
    "                        module) (Default: 4).\n"
    "  --truncateZeroPvalues, -tzpv\n"
    "                        Sets all p-values which are equal to zero to one. This\n"
    "                        has the effect that the associated positions are not\n"
    "                        part of the significance decision.\n"
    "  --fixateRange FIXATERANGE, -fs FIXATERANGE\n"
    "                        Fixate range of backgroundmodel starting at distance\n"
    "                        x. E.g. all values greater than 500kb are set to the\n"
    "                        value of the 500kb bin (Default: 500000).\n"
    "  --peakInteractionsThreshold PEAKINTERACTIONSTHRESHOLD, -pit PEAKINTERACTIONSTHRESHOLD\n"
    "                        The minimum number of interactions a detected peak\n"
    "                        needs to have to be considered (Default: 5).\n"
    "  --correctForMultipleTesting {none,fdr,bonferroni}\n"
    "                        Adjust the tested p-values across all viewpoints of\n"
    "                        all samples, Benjamini-Hochberg (fdr) or Bonferroni,\n"
    "                        and apply --pValue to the adjusted values; the\n"
    "                        significant file gains pvalue_adjusted. Not in the\n"
    "                        Python tool; none gives its output (Default: none).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

class PythonError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

// The Python float key of a numpy value as an int dict key: a float equal to
// an integer finds the int key.
std::optional<std::int64_t> integral(double key) {
    if (!std::isfinite(key) || std::trunc(key) != key || std::fabs(key) > 9.0e18) {
        return std::nullopt;
    }
    return static_cast<std::int64_t>(key);
}

std::string python_key_repr(double key) {
    if (const auto value = integral(key)) {
        return std::to_string(*value);
    }
    char buffer[64];
    std::snprintf(buffer, sizeof buffer, "%.17g", key);
    return buffer;
}

// The value of --pValue, --loosePValue or --xFoldBackground after main's
// conversion: a float, read_threshold_file's dict, or the string itself when
// it was empty (`if args.pValue:` leaves an empty string unconverted).
struct Threshold {
    enum class Kind { Float, Dict, Text } kind = Kind::Float;
    double value = 0.0;
    std::map<std::int64_t, double> dict;
    std::string text;

    // pPValue[key]
    [[nodiscard]] double at(double key) const {
        const auto k = integral(key);
        if (k.has_value()) {
            const auto it = dict.find(*k);
            if (it != dict.end()) {
                return it->second;
            }
        }
        throw PythonError("KeyError: " + python_key_repr(key));
    }
};

std::vector<std::string> split(const std::string& text, char separator) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true) {
        const std::size_t end = text.find(separator, begin);
        if (end == std::string::npos) {
            fields.push_back(text.substr(begin));
            return fields;
        }
        fields.push_back(text.substr(begin, end - begin));
        begin = end + 1;
    }
}

// read_threshold_file(pFile)
std::map<std::int64_t, double> read_threshold_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw PythonError("FileNotFoundError: [Errno 2] No such file or directory: '" + path + "'");
    }
    std::map<std::int64_t, double> out;
    for (const std::string& raw : chic::read_lines(path)) {
        const std::string line(chic::strip(raw));
        if (line.rfind('#', 0) == 0) {
            continue;
        }
        if (line.empty()) {
            break;
        }
        const std::vector<std::string> fields = split(line, '\t');
        if (fields.size() != 2) {
            throw PythonError(fields.size() < 2
                                  ? "ValueError: not enough values to unpack (expected 2, got 1)"
                                  : "ValueError: too many values to unpack (expected 2)");
        }
        out[chic::python_int(fields[0])] = chic::python_float(fields[1]);
    }
    return out;
}

Threshold convert_threshold(const std::string& text) {
    Threshold threshold;
    if (text.empty()) {
        threshold.kind = Threshold::Kind::Text;
        return threshold;
    }
    try {
        threshold.value = chic::python_float(text);
        threshold.kind = Threshold::Kind::Float;
    } catch (const chic::ValueError&) {
        threshold.kind = Threshold::Kind::Dict;
        threshold.dict = read_threshold_file(text);
    }
    threshold.text = text;
    return threshold;
}

struct Options {
    Threshold p_value;
    std::optional<Threshold> x_fold;
    std::optional<Threshold> loose;
    bool truncate_zero = false;
    double peak_threshold = 0.0;
    std::int64_t resolution = 0;
    std::string correction = "none";
};

using Entry = std::pair<double, InteractionRecord>;

// Viewpoint.merge_neighbors(pScoresDictionary, pMergeThreshold=resolution):
// the merged scores and, per merged key, the (first start, last end) the
// target line takes, both in dict insertion order.
struct Merged {
    std::vector<Entry> scores;
    std::map<double, std::pair<std::int64_t, std::int64_t>> spans;
    std::map<double, std::string> chromosomes;
};

Merged merge_neighbors(const std::vector<Entry>& candidates, std::int64_t threshold) {
    std::vector<std::vector<std::size_t>> merge_ids;
    std::vector<std::size_t> non_merge;
    for (std::size_t i = 0; i + 1 < candidates.size(); ++i) {
        const std::int64_t pre = candidates[i].second.relative_position;
        const std::int64_t suc = candidates[i + 1].second.relative_position;
        if (std::llabs(pre - suc) <= threshold) {
            if (!merge_ids.empty() && merge_ids.back().back() == i) {
                merge_ids.back().push_back(i + 1);
            } else {
                merge_ids.push_back({i, i + 1});
            }
        } else {
            // `if i == len(key_list) - 1` never holds, so an unmerged last
            // candidate is never kept.
            if (!merge_ids.empty() && merge_ids.back().back() != i) {
                non_merge.push_back(i);
            } else if (merge_ids.empty()) {
                non_merge.push_back(i);
            }
        }
    }
    Merged merged;
    for (const auto& element : merge_ids) {
        const InteractionRecord& first = candidates[element.front()].second;
        double sums[4] = {first.interaction, first.pvalue, first.xfold, first.raw};
        double max_value = first.raw;
        std::size_t index_maximum = 0;
        for (std::size_t k = 1; k < element.size(); ++k) {
            const InteractionRecord& record = candidates[element[k]].second;
            sums[0] += record.interaction;
            sums[1] += record.pvalue;
            sums[2] += record.xfold;
            sums[3] += record.raw;
            if (max_value < record.raw) {
                max_value = record.raw;
                index_maximum = k;
            }
        }
        const Entry& maximum = candidates[element[index_maximum]];
        InteractionRecord result = maximum.second;
        // base_element[-6] of the first element is overwritten with each
        // successor's sum of interactions; it is the result only when the
        // first element is the maximum.
        if (index_maximum == 0) {
            result.sum_of_interactions = candidates[element.back()].second.sum_of_interactions;
        }
        result.interaction = sums[0];
        result.pvalue = sums[1];
        result.xfold = sums[2];
        result.raw = sums[3];
        result.end = candidates[element.back()].second.end;
        result.start = first.start;
        merged.scores.emplace_back(maximum.first, result);
        merged.spans[maximum.first] = {first.start, candidates[element.back()].second.end};
        merged.chromosomes[maximum.first] = first.chromosome;
    }
    for (const std::size_t index : non_merge) {
        const Entry& entry = candidates[index];
        merged.scores.push_back(entry);
        merged.spans[entry.first] = {entry.second.start, entry.second.end};
        merged.chromosomes[entry.first] = entry.second.chromosome;
    }
    return merged;
}

struct TargetLine {
    std::string chromosome;
    std::int64_t start = 0;
    std::int64_t end = 0;
};

// Adjusted p-values of the family, keyed by (sample viewpoint, position).
using Adjusted = std::map<std::pair<std::size_t, double>, double>;

struct SampleResult {
    std::vector<Entry> accepted;
    std::vector<TargetLine> target;
};

// The tested candidates of one sample viewpoint, in the order the tool tests
// them: (position, p-value).
using Tested = std::vector<std::pair<double, double>>;

// compute_interaction_file's body for one sample. `tested`, when given,
// receives the p-values the tool compares against --pValue; `adjusted`, when
// given, replaces them in the decision.
SampleResult process_sample(const InteractionTable& table, const Options& options,
                            const chic::BackgroundModel& background, std::size_t sample_id,
                            Tested* tested, const Adjusted* adjusted) {
    SampleResult result;
    const auto decide = [&](double key, double pvalue) {
        if (adjusted == nullptr) {
            return pvalue;
        }
        const auto it = adjusted->find({sample_id, key});
        return it == adjusted->end() ? std::numeric_limits<double>::quiet_NaN() : it->second;
    };

    std::vector<Entry> data;
    data.reserve(table.keys.size());
    for (const double key : table.keys) {
        data.emplace_back(key, table.records.at(key));
    }

    if (options.x_fold.has_value() || options.loose.has_value()) {
        std::vector<Entry> candidates;
        if (options.x_fold.has_value()) {
            const Threshold& threshold = *options.x_fold;
            if (threshold.kind != Threshold::Kind::Text) {
                for (const Entry& entry : data) {
                    const double limit = threshold.kind == Threshold::Kind::Float
                                             ? threshold.value
                                             : threshold.at(entry.first);
                    if (entry.second.xfold < limit) {
                        continue;
                    }
                    candidates.push_back(entry);
                }
            }
        } else {
            const Threshold& threshold = *options.loose;
            if (threshold.kind != Threshold::Kind::Text) {
                for (const Entry& entry : data) {
                    const double limit = threshold.kind == Threshold::Kind::Float
                                             ? threshold.value
                                             : threshold.at(entry.first);
                    const double pvalue = entry.second.pvalue;
                    if (options.truncate_zero ? (pvalue == 0.0 || pvalue > limit)
                                              : pvalue > limit) {
                        continue;
                    }
                    candidates.push_back(entry);
                }
            }
        }
        if (candidates.empty()) {
            return result;
        }
        Merged merged = merge_neighbors(candidates, options.resolution);
        // compute_new_p_values
        const Threshold& threshold = options.p_value;
        if (threshold.kind == Threshold::Kind::Text) {
            return result;
        }
        for (Entry& entry : merged.scores) {
            const auto key = integral(entry.first);
            if (!key.has_value() || !background.contains(*key)) {
                continue;
            }
            const std::vector<double>& model = background.at(*key);
            InteractionRecord& record = entry.second;
            record.pvalue = 1.0 - hicx::scipy::betainc(model.at(0), record.raw + 1.0, model.at(1));
            if (tested != nullptr) {
                tested->emplace_back(entry.first, record.pvalue);
                continue;
            }
            const double limit = threshold.kind == Threshold::Kind::Float
                                     ? threshold.value
                                     : threshold.at(entry.first);
            if (decide(entry.first, record.pvalue) <= limit) {
                if (record.raw >= options.peak_threshold) {
                    result.accepted.push_back(entry);
                    const auto span = merged.spans.at(entry.first);
                    result.target.push_back(
                        {merged.chromosomes.at(entry.first), span.first, span.second});
                }
            }
        }
        return result;
    }

    // filter_by_pvalue(data[0], pValue, data[1], peakInteractionsThreshold)
    const Threshold& threshold = options.p_value;
    if (threshold.kind == Threshold::Kind::Text) {
        return result;
    }
    for (const Entry& entry : data) {
        if (tested != nullptr) {
            tested->emplace_back(entry.first, entry.second.pvalue);
            continue;
        }
        const double limit = threshold.kind == Threshold::Kind::Float
                                 ? threshold.value
                                 : threshold.at(entry.first);
        if (decide(entry.first, entry.second.pvalue) <= limit) {
            if (entry.second.xfold >= options.peak_threshold) {
                result.accepted.push_back(entry);
                result.target.push_back(
                    {entry.second.chromosome, entry.second.start, entry.second.end});
            }
        }
    }
    return result;
}

std::string python_str_repr(const std::string& text) {
    const bool single = text.find('\'') != std::string::npos;
    const bool dbl = text.find('"') != std::string::npos;
    const char quote = (single && !dbl) ? '"' : '\'';
    std::string out(1, quote);
    for (const unsigned char c : text) {
        if (c == '\\') {
            out += "\\\\";
        } else if (c == static_cast<unsigned char>(quote)) {
            out += '\\';
            out += static_cast<char>(c);
        } else if (c == '\n') {
            out += "\\n";
        } else if (c == '\t') {
            out += "\\t";
        } else if (c == '\r') {
            out += "\\r";
        } else if (c < 0x20 || c == 0x7f) {
            char buffer[8];
            std::snprintf(buffer, sizeof buffer, "\\x%02x", c);
            out += buffer;
        } else {
            out += static_cast<char>(c);
        }
    }
    out += quote;
    return out;
}

std::string python_repr(const std::vector<Triplet>& combination) {
    std::string out = "[";
    for (std::size_t i = 0; i < combination.size(); ++i) {
        out += i > 0 ? ", [" : "[";
        for (std::size_t j = 0; j < combination[i].size(); ++j) {
            out += (j > 0 ? ", " : "") + python_str_repr(combination[i][j]);
        }
        out += "]";
    }
    return out + "]";
}

using ReferencePoint = std::vector<std::optional<std::int64_t>>;

// np.array(list of reference point lists): a ragged list raises.
void check_homogeneous(const std::vector<ReferencePoint>& points) {
    for (const ReferencePoint& point : points) {
        if (point.size() != points.front().size()) {
            throw PythonError(
                "ValueError: setting an array element with a sequence. The requested array has "
                "an inhomogeneous shape after 1 dimensions. The detected shape was (" +
                std::to_string(points.size()) + ",) + inhomogeneous part.");
        }
    }
}

// int(pReferencePointsList[i][k])
std::int64_t reference_value(const std::vector<ReferencePoint>& points, std::size_t i,
                             std::size_t k) {
    if (i >= points.size() || k >= points[i].size()) {
        throw PythonError("IndexError: index " + std::to_string(k) +
                          " is out of bounds for axis 1 with size " +
                          std::to_string(i < points.size() ? points[i].size() : 0));
    }
    if (!points[i][k].has_value()) {
        throw PythonError("TypeError: int() argument must be a string, a bytes-like object or a "
                          "real number, not 'NoneType'");
    }
    return *points[i][k];
}

void set_threshold_attribute(chic::Hdf5Writer& writer, const std::string& name,
                             const Threshold& threshold) {
    switch (threshold.kind) {
        case Threshold::Kind::Float:
            writer.set_attribute("/", name, threshold.value);
            return;
        case Threshold::Kind::Text:
            writer.set_attribute("/", name, threshold.text);
            return;
        case Threshold::Kind::Dict:
            throw PythonError("TypeError: Object dtype dtype('O') has no native HDF5 equivalent");
    }
}

void write_root_attributes(chic::Hdf5Writer& writer, const std::string& type,
                           const Options& options, const cli::Namespace& args) {
    writer.set_attribute("/", "type", type);
    writer.set_attribute("/", "version", std::string(hicx::kVersion));
    set_threshold_attribute(writer, "pvalue", options.p_value);
    if (options.x_fold.has_value()) {
        writer.set_attribute("/", "mode_preselection", std::string("xfold"));
        set_threshold_attribute(writer, "mode_preselection_value", *options.x_fold);
    } else if (options.loose.has_value()) {
        writer.set_attribute("/", "mode_preselection", std::string("loosePValue"));
        set_threshold_attribute(writer, "mode_preselection_value", *options.loose);
    } else {
        writer.set_attribute("/", "mode_preselection", std::string("None"));
        writer.set_attribute("/",
                             type == "target" ? "mode_preselection_calue"
                                              : "mode_preselection_value",
                             std::string("None"));
    }
    const std::vector<std::int64_t> range = args.integers("range");
    writer.set_attribute("/", "range", std::span<const std::int64_t>(range));
    writer.set_attribute("/", "combinationMode", args.str("combinationMode"));
    writer.set_bool_attribute("/", "truncateZeroPvalues", options.truncate_zero);
    writer.set_attribute("/", "fixateRange", args.integer("fixateRange"));
    writer.set_attribute("/", "peakInteractionsThreshold",
                         args.integer("peakInteractionsThreshold"));
    if (options.correction != "none") {
        writer.set_attribute("/", "correctForMultipleTesting", options.correction);
    }
}

// The gene group name made unique per outer matrix (keys_seen[key[0]]).
std::string unique_gene_name(std::set<std::string>& seen, const std::string& gene) {
    for (std::size_t counter = 0;; ++counter) {
        const std::string name = counter == 0 ? gene : gene + "_" + std::to_string(counter);
        if (seen.insert(name).second) {
            return name;
        }
    }
}

void create_unique_group(chic::Hdf5Writer& writer, const std::string& path) {
    if (writer.exists(path)) {
        throw std::runtime_error(
            "the group " + path +
            " exists already. The Python reference does not terminate on this input "
            "(Viewpoint.createUniqueHDFGroup never changes the name it retries); the C++ port "
            "stops instead.");
    }
    writer.create_group(path);
}

// pybedtools BedTool(lines).sort().merge(d=resolution), bedtools 2.31.
std::vector<TargetLine> sort_merge(std::vector<TargetLine> lines, std::int64_t distance) {
    std::stable_sort(lines.begin(), lines.end(), [](const TargetLine& a, const TargetLine& b) {
        if (a.chromosome != b.chromosome) {
            return a.chromosome < b.chromosome;
        }
        return a.start < b.start;
    });
    std::vector<TargetLine> merged;
    for (const TargetLine& line : lines) {
        if (!merged.empty() && merged.back().chromosome == line.chromosome &&
            line.start - merged.back().end <= distance) {
            merged.back().end = std::max(merged.back().end, line.end);
        } else {
            merged.push_back(line);
        }
    }
    return merged;
}

}  // namespace

int main(int argc, char** argv) {
    cli::Parser parser("chicSignificantInteractions",
                       "Per viewpoint the significant interactions are detected based on the "
                       "background model.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--interactionFile", "-if"})
        .required()
        .input({"hdf5"})
        .help("path to the interaction file (HDF5) which should be used for aggregation of the "
              "statistics.");
    required.add({"--pValue", "-p"})
        .required()
        .help("p-value threshold to filter target regions for inclusion in differential "
              "analysis.");
    cli::ArgumentGroup& options_group = parser.group("options");
    cli::MutuallyExclusiveGroup& preselection = parser.mutually_exclusive(options_group);
    preselection.add({"--xFoldBackground", "-xf"})
        .help("Filter x-fold over background. Used to merge neighboring bins with a broader peak "
              "but less significant interactions to a single peak with high significance. Used "
              "only for pValue option.");
    preselection.add({"--loosePValue", "-lp"})
        .help("loose p-value threshold to filter target regions in a first round. Used to merge "
              "neighboring bins with a broader peak but less significant interactions to a "
              "single peak with high significance. Used only for pValue option.");
    required.add({"--backgroundModelFile", "-bmf"})
        .required()
        .input({"txt"})
        .help("path to the background file.");
    required.add({"--range"})
        .required()
        .type("int")
        .nargs(2)
        .help("Defines the region upstream and downstream of a reference point which should be "
              "included. Format is --region upstream downstream.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--outFileNameSignificant", "-os"})
        .default_value("significantInteractions.hdf5")
        .output({"hdf5"})
        .help("File name suffix to save the results; prefix is the input file name.");
    optional.add({"--outFileNameTarget", "-ot"})
        .default_value("targetFile.hdf5")
        .output({"hdf5"})
        .help("The file to store the target data.");
    optional.add({"--combinationMode", "-cm"})
        .default_value("dual")
        .choices({"dual", "single"})
        .help("This option defines how the interaction data should be computed and combined.");
    optional.add({"--threads", "-t"})
        .type("int")
        .default_value(4)
        .help("Number of threads.");
    optional.add({"--truncateZeroPvalues", "-tzpv"})
        .action(cli::Action::StoreTrue)
        .help("Sets all p-values which are equal to zero to one.");
    optional.add({"--fixateRange", "-fs"})
        .type("int")
        .default_value(500000)
        .help("Fixate range of backgroundmodel starting at distance x.");
    optional.add({"--peakInteractionsThreshold", "-pit"})
        .type("int")
        .default_value(5)
        .help("The minimum number of interactions a detected peak needs to have to be "
              "considered.");
    optional.add({"--correctForMultipleTesting"})
        .choices({"none", "fdr", "bonferroni"})
        .default_value("none")
        .cpp_only("Multiple testing correction of the tested p-values across all viewpoints "
                  "(PLAN.md 9.7); none gives the Python output.")
        .help("Adjust the tested p-values across all viewpoints of all samples, "
              "Benjamini-Hochberg (fdr) or Bonferroni, and apply --pValue to the adjusted "
              "values.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    const cli::Namespace args = parser.parse(argc, argv);

    try {
        Options options;
        options.p_value = convert_threshold(args.str("pValue"));
        if (const auto loose = args.opt_str("loosePValue")) {
            options.loose = convert_threshold(*loose);
        }
        if (const auto x_fold = args.opt_str("xFoldBackground")) {
            options.x_fold = convert_threshold(*x_fold);
        }
        options.truncate_zero = args.flag("truncateZeroPvalues");
        options.peak_threshold = static_cast<double>(args.integer("peakInteractionsThreshold"));
        options.correction = args.str("correctForMultipleTesting");
        const std::vector<std::int64_t> range = args.integers("range");

        const chic::BackgroundModel background = chic::read_background_model(
            args.str("backgroundModelFile"), range.at(0), range.at(1), args.integer("fixateRange"),
            false);

        const h5::File file(args.str("interactionFile"));
        const std::vector<std::string> samples = file.children("/");
        const auto root = file.attributes("/");
        const auto type = root.find("type");
        if (type == root.end()) {
            throw PythonError("KeyError: \"Can't open attribute (can't locate attribute: 'type')\"");
        }
        const auto* type_text = std::get_if<std::string>(&type->second);
        if (type_text == nullptr || *type_text != "interactions") {
            std::fprintf(stderr, "Please provide a file created by chicViewpoint for the "
                                 "parameter --interactionFile.\n");
            return 1;
        }

        const auto chromosomes_of = [&](const std::string& sample) {
            std::vector<std::string> names = file.children("/" + sample);
            std::sort(names.begin(), names.end());
            const auto genes = std::find(names.begin(), names.end(), "genes");
            if (genes == names.end()) {
                throw PythonError("ValueError: list.remove(x): x not in list");
            }
            names.erase(genes);
            return names;
        };
        const auto genes_of = [&](const std::string& sample, const std::string& chromosome) {
            std::vector<std::string> names = file.children("/" + sample + "/" + chromosome);
            std::sort(names.begin(), names.end());
            return names;
        };

        std::vector<std::vector<Triplet>> combinations;
        if (args.str("combinationMode") == "dual") {
            if (samples.size() > 1) {
                for (std::size_t i = 0; i < samples.size(); ++i) {
                    for (std::size_t j = i + 1; j < samples.size(); ++j) {
                        for (const std::string& chromosome : chromosomes_of(samples[i])) {
                            for (const std::string& gene : genes_of(samples[i], chromosome)) {
                                combinations.push_back(
                                    {{samples[i], chromosome, gene}, {samples[j], chromosome, gene}});
                            }
                        }
                    }
                }
            } else {
                std::fprintf(stderr, "Dual mode selected but only one matrix is stored\n");
            }
        } else {
            for (const std::string& sample : samples) {
                for (const std::string& chromosome : chromosomes_of(sample)) {
                    for (const std::string& gene : genes_of(sample, chromosome)) {
                        combinations.push_back({{sample, chromosome, gene}});
                    }
                }
            }
        }

        const auto resolution = root.find("resolution");
        if (resolution == root.end()) {
            throw PythonError(
                "KeyError: \"Can't open attribute (can't locate attribute: 'resolution')\"");
        }
        if (const auto* integer = std::get_if<std::int64_t>(&resolution->second)) {
            options.resolution = *integer;
        } else if (const auto* real = std::get_if<double>(&resolution->second)) {
            options.resolution = static_cast<std::int64_t>(*real);
        }

        const std::int64_t threads = args.integer("threads");
        if (threads == 0) {
            throw PythonError("ZeroDivisionError: integer division or modulo by zero");
        }
        if (threads < 0) {
            // No worker is started, so nothing is computed.
            combinations.clear();
        }

        // The sample viewpoints in computation order, each read once.
        std::map<Triplet, std::size_t> sample_ids;
        std::vector<InteractionTable> tables;
        for (const auto& combination : combinations) {
            for (const Triplet& sample : combination) {
                if (sample_ids.emplace(sample, tables.size()).second) {
                    tables.push_back(chic::read_interaction_table(file, sample));
                }
            }
        }

        std::optional<Adjusted> adjusted;
        if (options.correction != "none") {
            std::vector<std::pair<std::size_t, double>> keys;
            std::vector<double> values;
            std::vector<bool> done(tables.size(), false);
            for (const auto& combination : combinations) {
                for (const Triplet& sample : combination) {
                    const std::size_t id = sample_ids.at(sample);
                    if (done[id]) {
                        continue;
                    }
                    done[id] = true;
                    Tested tested;
                    (void)process_sample(tables[id], options, background, id, &tested, nullptr);
                    for (const auto& [key, pvalue] : tested) {
                        keys.emplace_back(id, key);
                        values.push_back(pvalue);
                    }
                }
            }
            const std::vector<double> corrected =
                options.correction == "fdr"
                    ? hicx::stats::benjamini_hochberg_adjusted(values)
                    : hicx::stats::bonferroni_adjusted(values);
            adjusted.emplace();
            for (std::size_t k = 0; k < keys.size(); ++k) {
                adjusted->emplace(keys[k], corrected[k]);
            }
        }

        std::vector<std::vector<Entry>> significant_data;
        std::vector<Triplet> significant_keys;
        std::vector<ReferencePoint> significant_points;
        std::vector<std::vector<TargetLine>> target_data;
        std::vector<Triplet> target_keys;
        std::vector<ReferencePoint> target_points;
        std::string error_log;
        for (const auto& combination : combinations) {
            std::vector<TargetLine> target;
            Triplet prefix;
            ReferencePoint last_point;
            for (const Triplet& sample : combination) {
                const std::size_t id = sample_ids.at(sample);
                SampleResult result = process_sample(tables[id], options, background, id, nullptr,
                                                     adjusted ? &*adjusted : nullptr);
                prefix.push_back(sample.at(0));
                if (result.accepted.empty()) {
                    error_log += "Failed for: " + python_repr(combination) + ".\n";
                }
                target.insert(target.end(), result.target.begin(), result.target.end());
                significant_data.push_back(std::move(result.accepted));
                significant_keys.push_back(sample);
                significant_points.push_back(tables[id].reference_point);
                last_point = tables[id].reference_point;
            }
            prefix.push_back("::");
            prefix.push_back(combination.front().at(1));
            prefix.push_back(combination.front().at(2));
            target_data.push_back(std::move(target));
            target_keys.push_back(std::move(prefix));
            target_points.push_back(std::move(last_point));
        }
        if (!error_log.empty()) {
            std::ofstream log("errorLog.txt", std::ios::app | std::ios::binary);
            log << error_log;
        }
        check_homogeneous(target_points);
        check_homogeneous(significant_points);

        // np.unique(significant_key_list, axis=0, return_index=True): the
        // sorted distinct triplets, each with its first occurrence.
        std::vector<std::size_t> order(significant_keys.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return significant_keys[a] < significant_keys[b];
        });
        std::vector<std::size_t> unique;
        for (const std::size_t index : order) {
            if (unique.empty() || significant_keys[unique.back()] != significant_keys[index]) {
                unique.push_back(index);
            }
        }

        // writeSignificantHDF
        {
            chic::Hdf5Writer writer(args.str("outFileNameSignificant"));
            write_root_attributes(writer, "significant", options, args);
            std::map<std::string, std::set<std::string>> keys_seen;
            for (std::size_t i = 0; i < unique.size(); ++i) {
                const Triplet& key = significant_keys[unique[i]];
                const std::vector<Entry>& data = significant_data[unique[i]];
                if (data.empty()) {
                    continue;
                }
                const InteractionRecord& last = data.back().second;
                chic::InteractionFileData content;
                content.chromosome = last.chromosome;
                content.gene = last.gene;
                content.sum_of_interactions = last.sum_of_interactions;
                std::vector<double> adjusted_values;
                for (const Entry& entry : data) {
                    content.starts.push_back(entry.second.start);
                    content.ends.push_back(entry.second.end);
                    content.relative_positions.push_back(entry.second.relative_position);
                    content.interaction_data.push_back(entry.second.interaction);
                    content.pvalues.push_back(entry.second.pvalue);
                    content.xfold.push_back(entry.second.xfold);
                    content.raw.push_back(entry.second.raw);
                    if (adjusted.has_value()) {
                        const auto it = adjusted->find(
                            {sample_ids.at(key),
                             static_cast<double>(entry.second.relative_position)});
                        adjusted_values.push_back(it == adjusted->end()
                                                      ? std::numeric_limits<double>::quiet_NaN()
                                                      : it->second);
                    }
                }
                const std::string matrix_group = key.at(0);
                if (!writer.exists(matrix_group)) {
                    writer.create_group(matrix_group);
                    keys_seen[matrix_group].clear();
                }
                const std::string chromosome_group = matrix_group + "/" + content.chromosome;
                if (!writer.exists(chromosome_group)) {
                    writer.create_group(chromosome_group);
                }
                if (!writer.exists(matrix_group + "/genes")) {
                    writer.create_group(matrix_group + "/genes");
                }
                const std::string name = unique_gene_name(keys_seen[matrix_group], key.at(2));
                const std::string group = chromosome_group + "/" + name;
                create_unique_group(writer, group);
                chic::write_interaction_datasets(writer, group, content,
                                                 reference_value(significant_points, i, 0),
                                                 reference_value(significant_points, i, 1));
                if (adjusted.has_value()) {
                    writer.write_array(group + "/pvalue_adjusted",
                                       std::span<const double>(adjusted_values), 9);
                }
                (void)writer.hard_link(group, matrix_group + "/genes/" + name);
            }
        }

        // writeTargetHDF
        {
            chic::Hdf5Writer writer(args.str("outFileNameTarget"));
            write_root_attributes(writer, "target", options, args);
            std::map<std::string, std::set<std::string>> keys_seen;
            for (std::size_t i = 0; i < target_keys.size(); ++i) {
                const Triplet& key = target_keys[i];
                if (target_data[i].empty()) {
                    continue;
                }
                const std::vector<TargetLine> merged =
                    sort_merge(target_data[i], options.resolution);
                std::vector<std::string> starts;
                std::vector<std::string> ends;
                for (const TargetLine& line : merged) {
                    starts.push_back(std::to_string(line.start));
                    ends.push_back(std::to_string(line.end));
                }
                const std::string chromosome = merged.back().chromosome;
                std::string matrix_group = key.at(0);
                if (!writer.exists(matrix_group)) {
                    writer.create_group(matrix_group);
                    keys_seen[key.at(0)].clear();
                }
                for (std::size_t k = 1; k < key.size(); ++k) {
                    if (key[k] == "::") {
                        break;
                    }
                    const std::string inner = matrix_group + "/" + key[k];
                    if (!writer.exists(inner)) {
                        writer.create_group(inner);
                        keys_seen[key[k]].clear();
                    }
                    matrix_group = inner;
                }
                if (!writer.exists(matrix_group + "/genes")) {
                    writer.create_group(matrix_group + "/genes");
                }
                const std::string chromosome_group = matrix_group + "/" + chromosome;
                if (!writer.exists(chromosome_group)) {
                    writer.create_group(chromosome_group);
                }
                const std::string name = unique_gene_name(keys_seen[key.at(0)], key.back());
                const std::string group = chromosome_group + "/" + name;
                create_unique_group(writer, group);
                writer.write_string(group + "/chromosome", chromosome);
                writer.write_strings(group + "/start_list", starts);
                writer.write_strings(group + "/end_list", ends);
                writer.write_scalar(group + "/reference_point_start",
                                    reference_value(target_points, i, 0));
                writer.write_scalar(group + "/reference_point_end",
                                    reference_value(target_points, i, 1));
                (void)writer.hard_link(group, matrix_group + "/genes/" + name);
            }
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chicSignificantInteractions: %s\n", error.what());
        return 1;
    }
    return 0;
}
