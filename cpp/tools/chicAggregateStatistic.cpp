// Port of hicexplorer/chicAggregateStatistic.py.
//
// For every pair of samples of a chicViewpoint interaction file and every
// reference point present in both, the positions of each sample are grouped by
// the target region they overlap (from a chicSignificantInteractions target
// file, a BED of three columns, or a BED of four columns naming the gene), and
// the p-values, x-folds and raw interactions of each group are summed into one
// line. The result is the aggregate file chicDifferentialTest reads.
//
// Behaviour reproduced as the Python has it:
//
//   * A position belongs to the first target region, in (start, end, file
//     order), that overlaps it; its line keeps the first position's start,
//     gene and sum of interactions, and takes the last position's end and
//     relative distance.
//   * intervalListToIntervalTree starts a new tree whenever the chromosome
//     changes, so a chromosome that reappears later in a BED file keeps only
//     its last run of regions.
//   * A BED of four columns is always read on the whitespace split (the tab
//     split unpacks four fields into three names and raises), and its targets
//     are filed under the fixed matrix names c_adj_norm and t_adj_norm
//     (utilities.readTargetBed); any other interaction file fails with
//     KeyError. A line of exactly three tab separated fields reuses the gene of
//     the line before.
//   * A target file written in single mode has no genes group below its
//     matrices and fails with KeyError.
//   * With exactly one target for all reference points, the target is not
//     passed to the worker and a target file of any kind but a three column
//     BED fails with TypeError.
//   * --threads is lowered to the number of reference points; nothing to do,
//     or --threads 0, raises ZeroDivisionError; a negative --threads starts no
//     worker and writes the file with its attributes only.
//   * A reference point listed twice (a gene name in two target matrices of an
//     HDF5 target file) fails when its group is created a second time.
//
// Threading: as for chicViewpoint (cpp/OPTIMIZATION.md 6) the loop is
// sequential. The Python collects its workers in order, so its output does not
// depend on --threads either.

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/chic_hdf5.hpp"
#include "hicx/chic_viewpoint.hpp"
#include "hicx/hdf5_util.hpp"
#include "hicx/version.hpp"

namespace {

namespace chic = hicx::chic;
namespace h5 = hicx::h5;
namespace cli = hicx::cli;

using chic::InteractionRecord;
using chic::InteractionTable;
using Triplet = std::vector<std::string>;

const char* const kUsage =
    "usage: chicAggregateStatistic --interactionFile INTERACTIONFILE\n"
    "                              [--targetFile TARGETFILE]\n"
    "                              [--outFileName OUTFILENAME] [--threads THREADS]\n"
    "                              [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "chicAggregateStatistic is a preprocessing tool for chicDifferentialTest. It takes two "
    "consecutive viewpoint files and one target file and creates one\n"
    "file containing all locations which should be tested for differential interactions. "
    "Either one target file for two consecutive viewpoint files or one\n"
    "target file for all viewpoints is accepted.\n"
    "\n"
    "Required arguments:\n"
    "  --interactionFile INTERACTIONFILE, -if INTERACTIONFILE\n"
    "                        path to the interaction files which should be used for\n"
    "                        aggregation of the statistics.\n"
    "  --targetFile TARGETFILE, -tf TARGETFILE\n"
    "                        path to the target files which contains the target\n"
    "                        regions to prepare data for differential analysis.\n"
    "                        This is either the target file in the hdf format\n"
    "                        created by chicSignificantInteractions or a regular,\n"
    "                        three column bed file.\n"
    "\n"
    "Optional arguments:\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        File name to save the result (Default:\n"
    "                        aggregate_target.hdf).\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads (uses the python multiprocessing\n"
    "                        module)ist (Default: 4).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

class PythonError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

std::vector<std::string> split_tab(const std::string& text) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true) {
        const std::size_t end = text.find('\t', begin);
        if (end == std::string::npos) {
            fields.push_back(text.substr(begin));
            return fields;
        }
        fields.push_back(text.substr(begin, end - begin));
        begin = end + 1;
    }
}

// str.split() without arguments, for ASCII whitespace.
std::vector<std::string> split_whitespace(const std::string& text) {
    std::vector<std::string> fields;
    std::size_t i = 0;
    while (i < text.size()) {
        while (i < text.size() && std::isspace(static_cast<unsigned char>(text[i]))) {
            ++i;
        }
        std::size_t j = i;
        while (j < text.size() && !std::isspace(static_cast<unsigned char>(text[j]))) {
            ++j;
        }
        if (j > i) {
            fields.push_back(text.substr(i, j - i));
        }
        i = j;
    }
    return fields;
}

std::string unpack_error(std::size_t expected, std::size_t got) {
    return got < expected ? "ValueError: not enough values to unpack (expected " +
                                std::to_string(expected) + ", got " + std::to_string(got) + ")"
                          : "ValueError: too many values to unpack (expected " +
                                std::to_string(expected) + ")";
}

// A target region as intervalListToIntervalTree stores it.
struct Region {
    std::int64_t begin = 0;
    std::int64_t end = 0;
    std::int64_t id = 0;
    bool operator<(const Region& other) const {
        return std::tie(begin, end, id) < std::tie(other.begin, other.end, other.id);
    }
    bool operator==(const Region& other) const {
        return begin == other.begin && end == other.end && id == other.id;
    }
};

using Forest = std::map<std::string, std::vector<Region>>;

struct TextRegion {
    std::string chromosome;
    std::string start;
    std::string end;
};

// hiCMatrix.intervalListToIntervalTree(interval_list)[0]
Forest interval_forest(const std::vector<TextRegion>& regions) {
    Forest forest;
    std::optional<std::string> previous;
    std::int64_t id = 0;
    for (const TextRegion& region : regions) {
        const std::int64_t start = chic::python_int(region.start);
        const std::int64_t end = chic::python_int(region.end);
        if (!previous.has_value() || *previous != region.chromosome) {
            forest[region.chromosome].clear();
            previous = region.chromosome;
        }
        if (start >= end) {
            throw PythonError("ValueError: IntervalTree: Null Interval objects not allowed in "
                              "IntervalTree: Interval(" +
                              std::to_string(start) + ", " + std::to_string(end) + ", " +
                              std::to_string(id) + ")");
        }
        forest[region.chromosome].push_back({start, end, id});
        ++id;
    }
    return forest;
}

// utilities.readBed(pBedFile)
std::vector<TextRegion> read_bed(const std::string& path) {
    std::vector<TextRegion> regions;
    for (const std::string& line : chic::read_lines(path)) {
        if (line.rfind('#', 0) == 0) {
            continue;
        }
        const std::string stripped(chic::strip(line));
        std::vector<std::string> fields = split_tab(stripped);
        if (fields.size() < 3) {
            fields = split_whitespace(stripped);
            if (fields.size() < 3) {
                throw PythonError(unpack_error(3, fields.size()));
            }
        }
        regions.push_back({fields[0], fields[1], fields[2]});
    }
    return regions;
}

struct TargetPositions {
    std::string chromosome;
    std::vector<std::string> starts;
    std::vector<std::string> ends;
};

// present_genes: {outer: {inner: [gene, ...]}} with insertion order kept.
struct PresentGenes {
    std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::vector<std::string>>>>>
        outer;

    std::vector<std::string>& at(const std::string& outer_name, const std::string& inner_name) {
        for (auto& [name, inner] : outer) {
            if (name == outer_name) {
                for (auto& [inner_key, genes] : inner) {
                    if (inner_key == inner_name) {
                        return genes;
                    }
                }
                throw PythonError("KeyError: '" + inner_name + "'");
            }
        }
        throw PythonError("KeyError: '" + outer_name + "'");
    }
};

struct Targets {
    enum class Kind { None, Hdf5, Bed3, Bed4 } kind = Kind::None;
    PresentGenes present;
    std::map<std::string, std::vector<std::string>> target_dict;  // gene -> path
    std::map<std::string, TargetPositions> positions;             // bed4
};

// utilities.readTargetBed(pBedFile)
void read_target_bed(const std::string& path, Targets& targets) {
    targets.present.outer.push_back({"c_adj_norm", {{"t_adj_norm", {}}}});
    std::vector<std::string>& genes = targets.present.outer.back().second.back().second;
    std::optional<std::string> gene;
    for (const std::string& line : chic::read_lines(path)) {
        if (line.rfind('#', 0) == 0) {
            continue;
        }
        const std::string stripped(chic::strip(line));
        const std::vector<std::string> tab_fields = split_tab(stripped);
        std::string chromosome;
        std::string start;
        std::string end;
        if (tab_fields.size() == 3) {
            chromosome = tab_fields[0];
            start = tab_fields[1];
            end = tab_fields[2];
            if (!gene.has_value()) {
                throw PythonError("UnboundLocalError: cannot access local variable 'gene' where it "
                                  "is not associated with a value");
            }
        } else {
            const std::vector<std::string> fields = split_whitespace(stripped);
            if (fields.size() < 4) {
                throw PythonError(unpack_error(4, fields.size()));
            }
            chromosome = fields[0];
            start = fields[1];
            end = fields[2];
            gene = fields[3];
        }
        if (std::find(genes.begin(), genes.end(), *gene) == genes.end()) {
            genes.push_back(*gene);
        }
        targets.target_dict[*gene] = {"c_adj_norm", "t_adj_norm", "genes", *gene};
        TargetPositions& position = targets.positions[*gene];
        position.chromosome = chromosome;
        position.starts.push_back(start);
        position.ends.push_back(end);
    }
}

// Viewpoint.readTargetHDFFile(pFile)
void read_target_hdf(const std::string& path, Targets& targets) {
    const h5::File file(path);
    for (const std::string& outer : file.children("/")) {
        auto found = std::find_if(targets.present.outer.begin(), targets.present.outer.end(),
                                  [&](const auto& entry) { return entry.first == outer; });
        if (found == targets.present.outer.end()) {
            targets.present.outer.push_back({outer, {}});
            found = std::prev(targets.present.outer.end());
        }
        for (const std::string& inner : file.children("/" + outer)) {
            auto inner_found = std::find_if(found->second.begin(), found->second.end(),
                                            [&](const auto& entry) { return entry.first == inner; });
            if (inner_found == found->second.end()) {
                found->second.push_back({inner, {}});
                inner_found = std::prev(found->second.end());
            }
            const std::string genes_path = "/" + outer + "/" + inner + "/genes";
            if (!chic::contains(file, genes_path.substr(1))) {
                throw PythonError("KeyError: \"Unable to synchronously open object (object 'genes' "
                                  "doesn't exist)\"");
            }
            for (const std::string& gene : file.children(genes_path)) {
                targets.target_dict[gene] = {outer, inner, "genes", gene};
                inner_found->second.push_back(gene);
            }
        }
    }
}

// The start or end values of a target file dataset: bytes are decoded by
// int(), numbers truncated by int().
std::vector<std::string> read_coordinates(const h5::File& file, const std::string& path) {
    std::vector<std::string> out;
    const std::string dtype = file.dataset_dtype(path);
    if (dtype.find("int") != std::string::npos || dtype.find("float") != std::string::npos) {
        for (const double value : file.read_doubles(path)) {
            out.push_back(std::to_string(static_cast<std::int64_t>(value)));
        }
        return out;
    }
    return file.read_strings(path);
}

using Entry = std::pair<double, InteractionRecord>;

// filter_scores_target_list
std::vector<Entry> filter_scores(const InteractionTable& table, const Targets& targets,
                                 const std::vector<std::string>* target_path,
                                 const Forest* bed3_forest, const std::string& target_file) {
    Forest forest;
    const Forest* regions = nullptr;
    if (targets.kind == Targets::Kind::Hdf5 || targets.kind == Targets::Kind::Bed4) {
        if (target_path == nullptr) {
            throw PythonError(targets.kind == Targets::Kind::Hdf5
                                  ? "TypeError: can only join an iterable"
                                  : "TypeError: 'NoneType' object is not subscriptable");
        }
        std::vector<TextRegion> list;
        if (targets.kind == Targets::Kind::Hdf5) {
            const h5::File file(target_file);
            std::string joined;
            for (std::size_t i = 0; i < target_path->size(); ++i) {
                joined += (i > 0 ? "/" : "") + (*target_path)[i];
            }
            if (!chic::contains(file, joined)) {
                throw PythonError("KeyError: \"Unable to synchronously open object (object '" +
                                  target_path->back() + "' doesn't exist)\"");
            }
            const std::string base = "/" + joined + "/";
            if (!file.exists(base + "chromosome")) {
                throw PythonError("TypeError: 'NoneType' object is not subscriptable");
            }
            const std::string chromosome = file.read_strings(base + "chromosome").at(0);
            const std::vector<std::string> starts = read_coordinates(file, base + "start_list");
            const std::vector<std::string> ends = read_coordinates(file, base + "end_list");
            for (std::size_t i = 0; i < std::min(starts.size(), ends.size()); ++i) {
                list.push_back({chromosome, starts[i], ends[i]});
            }
        } else {
            const auto it = targets.positions.find(target_path->back());
            if (it == targets.positions.end()) {
                throw PythonError("KeyError: '" + target_path->back() + "'");
            }
            for (std::size_t i = 0; i < std::min(it->second.starts.size(), it->second.ends.size());
                 ++i) {
                list.push_back({it->second.chromosome, it->second.starts[i], it->second.ends[i]});
            }
        }
        if (list.empty()) {
            return {};
        }
        forest = interval_forest(list);
        regions = &forest;
    } else {
        regions = bed3_forest;
    }

    // same_target_dict, in first insertion order of its targets
    std::vector<std::pair<Region, std::vector<double>>> same_target;
    for (const double key : table.keys) {
        const InteractionRecord& record = table.records.at(key);
        const auto tree = regions->find(record.chromosome);
        if (tree == regions->end()) {
            continue;
        }
        std::optional<Region> first;
        if (record.start < record.end) {
            for (const Region& region : tree->second) {
                if (region.begin < record.end && region.end > record.start) {
                    if (!first.has_value() || region < *first) {
                        first = region;
                    }
                }
            }
        }
        if (!first.has_value()) {
            continue;
        }
        auto it = std::find_if(same_target.begin(), same_target.end(),
                               [&](const auto& entry) { return entry.first == *first; });
        if (it == same_target.end()) {
            same_target.push_back({*first, {key}});
        } else {
            it->second.push_back(key);
        }
    }

    std::vector<Entry> accepted;
    for (auto& [region, keys] : same_target) {
        std::sort(keys.begin(), keys.end());
        double sums[3] = {0.0, 0.0, 0.0};
        for (const double key : keys) {
            const InteractionRecord& record = table.records.at(key);
            sums[0] += record.pvalue;
            sums[1] += record.xfold;
            sums[2] += record.raw;
        }
        InteractionRecord line = table.records.at(keys.front());
        const InteractionRecord& last = table.records.at(keys.back());
        line.end = last.end;
        line.relative_position = last.relative_position;
        line.pvalue = sums[0];
        line.xfold = sums[1];
        line.raw = sums[2];
        accepted.emplace_back(keys.front(), line);
    }
    return accepted;
}

}  // namespace

int main(int argc, char** argv) {
    cli::Parser parser("chicAggregateStatistic",
                       "chicAggregateStatistic is a preprocessing tool for chicDifferentialTest.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--interactionFile", "-if"})
        .required()
        .input({"hdf5"})
        .help("path to the interaction files which should be used for aggregation of the "
              "statistics.");
    required.add({"--targetFile", "-tf"})
        .input({"hdf5", "bed"})
        .help("path to the target files which contains the target regions to prepare data for "
              "differential analysis.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--outFileName", "-o"})
        .default_value("aggregate_target.hdf")
        .output({"hdf5"})
        .help("File name to save the result.");
    optional.add({"--threads", "-t"})
        .type("int")
        .default_value(4)
        .help("Number of threads.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    const cli::Namespace args = parser.parse(argc, argv);

    try {
        const h5::File interactions(args.str("interactionFile"));
        const std::vector<std::string> samples = interactions.children("/");

        const std::optional<std::string> target_file = args.opt_str("targetFile");
        if (!target_file.has_value()) {
            throw PythonError("TypeError: expected str, bytes or os.PathLike object, not NoneType");
        }
        Targets targets;
        std::vector<std::vector<std::string>> target_list;  // per reference point
        if (h5::is_hdf5(*target_file)) {
            read_target_hdf(*target_file, targets);
            targets.kind = Targets::Kind::Hdf5;
        } else {
            std::optional<std::vector<std::string>> first_line;
            for (const std::string& line : chic::read_lines(*target_file)) {
                if (line.rfind('#', 0) == 0) {
                    continue;
                }
                first_line = split_tab(std::string(chic::strip(line)));
                break;
            }
            if (!first_line.has_value()) {
                throw PythonError("UnboundLocalError: cannot access local variable '_line' where "
                                  "it is not associated with a value");
            }
            if (first_line->size() == 4) {
                targets.kind = Targets::Kind::Bed4;
                read_target_bed(*target_file, targets);
            } else if (first_line->size() == 3) {
                targets.kind = Targets::Kind::Bed3;
                target_list.push_back({*target_file});
            } else {
                std::fprintf(stderr, "BED of targets list must have 3 or 4 columns\n");
            }
        }

        std::vector<std::vector<Triplet>> interaction_list;
        std::map<std::string, std::vector<Triplet>> interaction_dict;
        const auto sorted_children = [&](const std::string& path) {
            std::vector<std::string> names = interactions.children(path);
            std::sort(names.begin(), names.end());
            return names;
        };
        if (samples.size() > 1) {
            for (std::size_t i = 0; i < samples.size(); ++i) {
                for (std::size_t j = i + 1; j < samples.size(); ++j) {
                    std::vector<std::string> chromosomes1 = sorted_children("/" + samples[i]);
                    std::vector<std::string> chromosomes2 = sorted_children("/" + samples[j]);
                    for (auto* list : {&chromosomes1, &chromosomes2}) {
                        const auto genes = std::find(list->begin(), list->end(), "genes");
                        if (genes == list->end()) {
                            throw PythonError("ValueError: list.remove(x): x not in list");
                        }
                        list->erase(genes);
                    }
                    for (std::size_t c = 0; c < std::min(chromosomes1.size(), chromosomes2.size());
                         ++c) {
                        const std::vector<std::string> genes1 =
                            sorted_children("/" + samples[i] + "/" + chromosomes1[c]);
                        const std::vector<std::string> genes2 =
                            sorted_children("/" + samples[j] + "/" + chromosomes2[c]);
                        for (std::size_t g = 0; g < std::min(genes1.size(), genes2.size()); ++g) {
                            std::vector<Triplet> pair{{samples[i], chromosomes1[c], genes1[g]},
                                                      {samples[j], chromosomes2[c], genes2[g]}};
                            if (targets.kind != Targets::Kind::Bed3) {
                                const std::vector<std::string>& present =
                                    targets.present.at(samples[i], samples[j]);
                                if (std::find(present.begin(), present.end(), genes1[g]) !=
                                    present.end()) {
                                    interaction_dict[genes1[g]] = pair;
                                }
                            } else {
                                interaction_list.push_back(pair);
                            }
                        }
                    }
                }
            }
        } else {
            std::fprintf(stderr,
                         "To aggregate and prepare the data for the differential test, at least two "
                         "matrices need to be stored, but only one is present.\n");
        }

        if (targets.kind != Targets::Kind::Bed3) {
            for (const auto& [outer, inner_list] : targets.present.outer) {
                for (const auto& [inner, genes] : inner_list) {
                    for (const std::string& gene : genes) {
                        const auto pair = interaction_dict.find(gene);
                        if (pair == interaction_dict.end()) {
                            throw PythonError("KeyError: '" + gene + "'");
                        }
                        interaction_list.push_back(pair->second);
                        target_list.push_back(targets.target_dict.at(gene));
                    }
                }
            }
        }

        // call_multi_core
        std::int64_t threads = args.integer("threads");
        if (static_cast<std::int64_t>(interaction_list.size()) < threads) {
            threads = static_cast<std::int64_t>(interaction_list.size());
        }
        if (threads == 0) {
            throw PythonError("ZeroDivisionError: integer division or modulo by zero");
        }
        const bool one_target = target_list.size() == 1;

        std::vector<std::vector<Triplet>> names_list;
        std::vector<std::vector<std::vector<Entry>>> accepted_list;
        if (threads > 0) {
            Forest bed3_forest;
            if (targets.kind == Targets::Kind::Bed3) {
                bed3_forest = interval_forest(read_bed(one_target ? target_list.front().front()
                                                                  : std::string()));
            } else if (targets.kind == Targets::Kind::None) {
                throw PythonError("Exception: No target list given.");
            }
            for (std::size_t i = 0; i < interaction_list.size(); ++i) {
                std::vector<std::vector<Entry>> per_sample;
                for (const Triplet& sample : interaction_list[i]) {
                    const InteractionTable table = chic::read_interaction_table(interactions, sample);
                    const std::vector<std::string>* path =
                        one_target ? nullptr : &target_list.at(i);
                    per_sample.push_back(filter_scores(table, targets, path, &bed3_forest,
                                                       *target_file));
                }
                names_list.push_back(interaction_list[i]);
                accepted_list.push_back(std::move(per_sample));
            }
        }

        // writeAggregateHDF
        chic::Hdf5Writer writer(args.str("outFileName"));
        writer.set_attribute("/", "type", std::string("aggregate"));
        writer.set_attribute("/", "version", std::string(hicx::kVersion));
        for (std::size_t n = 0; n < names_list.size(); ++n) {
            const std::vector<Triplet>& keys = names_list[n];
            const std::string combination = keys.at(0).at(0) + "_" + keys.at(1).at(0);
            if (!writer.exists(combination)) {
                writer.create_group(combination);
            }
            for (std::size_t k = 0; k < keys.size(); ++k) {
                const std::vector<Entry>& data = accepted_list[n][k];
                if (data.empty()) {
                    continue;
                }
                const Triplet& key = keys[k];
                const InteractionRecord& last = data.back().second;
                std::vector<std::int64_t> starts;
                std::vector<std::int64_t> ends;
                std::vector<std::int64_t> relative;
                std::vector<double> raw;
                for (const Entry& entry : data) {
                    starts.push_back(entry.second.start);
                    ends.push_back(entry.second.end);
                    relative.push_back(entry.second.relative_position);
                    raw.push_back(entry.second.raw);
                }
                const std::string matrix_group = combination + "/" + key.at(0);
                if (!writer.exists(matrix_group)) {
                    writer.create_group(matrix_group);
                }
                std::string chromosome_group = matrix_group + "/" + key.at(1);
                if (!writer.exists(chromosome_group)) {
                    writer.create_group(chromosome_group);
                } else {
                    chromosome_group = matrix_group + "/" + last.chromosome;
                    if (!writer.exists(chromosome_group)) {
                        throw PythonError("KeyError: \"Unable to synchronously open object (object '" +
                                          last.chromosome + "' doesn't exist)\"");
                    }
                }
                if (!writer.exists(matrix_group + "/genes")) {
                    writer.create_group(matrix_group + "/genes");
                }
                const std::string group = chromosome_group + "/" + key.at(2);
                if (writer.exists(group)) {
                    throw PythonError("ValueError: Unable to synchronously create group (name "
                                      "already exists)");
                }
                writer.create_group(group);
                writer.write_string(group + "/chromosome", last.chromosome);
                writer.write_array(group + "/start_list", std::span<const std::int64_t>(starts), 9);
                writer.write_array(group + "/end_list", std::span<const std::int64_t>(ends), 9);
                writer.write_string(group + "/gene_name", last.gene);
                writer.write_scalar(group + "/sum_of_interactions", last.sum_of_interactions);
                writer.write_array(group + "/relative_distance_list",
                                   std::span<const std::int64_t>(relative), 9);
                writer.write_array(group + "/raw_target_list", std::span<const double>(raw), 9);
                if (!writer.hard_link(group, matrix_group + "/genes/" + key.at(2))) {
                    throw PythonError("OSError: Unable to create link (name already exists)");
                }
            }
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chicAggregateStatistic: %s\n", error.what());
        return 1;
    }
    return 0;
}
