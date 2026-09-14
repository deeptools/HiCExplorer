// Port of hicexplorer/hicAggregateContacts.py, numeric outputs only.
//
// The tool pools the submatrices around pairs of BED positions, optionally
// clusters them, and writes a matplotlib figure of the aggregate (required,
// --outFileName), optionally a diagnostic heatmap figure, and three numeric
// outputs: --outFilePrefixMatrix (np.savetxt '%0.5f' of each aggregate),
// --outFileContactPairs (the contact positions per cluster) and
// --outFileObsExp (the obs/exp matrix). The figures follow cpp/PLAN.md tier 7,
// option (a): after the tables the aggregate of every cluster, the cluster
// members and, for --diagnosticHeatmapFile, the diagonal of every kept
// submatrix go to plot/hicexplorer_plot/hicAggregateContacts.py, which draws
// them with the reference's plot_aggregated_contacts and
// plot_diagnostic_heatmaps. The C++-only options: --noPlot writes the numeric
// outputs without any figure; --plotData writes the figures' data as JSON
// instead of drawing.
//
// Clustering: main() reaches exactly two scikit-learn estimators, KMeans
// (--kmeans) and AgglomerativeClustering with ward linkage (--hclust). Both are
// hicx::cluster (core/include/hicx/clustering.hpp), which reproduces their
// labels. --spectral is parsed but never read by main(), so it is accepted and
// ignored here too. Measured 2026-09-13: every fit the Python tool made over
// 108 runs on Li_et_al_2015 (h5 and cool; all four transforms; full, center
// and diagonal vectors; k from 2 to 6; genome-wide and per chromosome; 81 to
// 443 samples of 1 to 961 features) was recorded with its input, and
// hicx::cluster returned identical labels for all 108 (54 KMeans, 54 ward),
// not only the same partitions.
//
// ---------------------------------------------------------------------------
// Reproduced reference defects, pinned by
// hicexplorer/test/general/test_hicAggregateContacts.py
// ---------------------------------------------------------------------------
//  1. --spectral is ignored.
//  2. The contact pair file of a cluster writes the centre value of its
//     cl_idx-th member next to the coordinates of the cl_idx-th submatrix of
//     all clusters (`coords[cl_idx]` instead of
//     `coords[cluster_indices[cl_idx]]`).
//  3. A '-'/'-' strand pair is transposed, not rotated, and the orientations
//     stay with their BED rows when two chromosomes are swapped.
//  4. get_outlier_indices scales by the median absolute *value*, not the
//     median absolute deviation, so for vectors that are mostly zero it
//     returns None and removes nothing.
//  5. Row-wise mode writes the coordinates as given, without reordering them
//     when the bins are swapped.
//  6. With --perChr, a chromosome with fewer submatrices than clusters sets k
//     to 1 for itself *and every chromosome after it*, while the file names
//     keep the `_cluster_1` suffix of the requested cluster count.
//  7. With --chromosomes on a matrix that has NaN bins, keepOnlyTheseChr
//     restores the bins maskBins removed, together with the bin table from
//     before enlarge_bins, and the matrix becomes float64. Positions in the
//     gaps between the original bins then resolve to no bin.
//
// Behaviour that differs from the Python, deliberately:
//  * a cluster left empty makes the Python raise IndexError after it has
//    written the files of the clusters before it; the port exits 1 before
//    writing any of them.
//
// Reduction orders that the output depends on, reproduced:
//  * np.matrix.sum() of a window is numpy's pairwise sum, over the logical
//    order for a fliplr or flipud view and over the memory order for a
//    transposed one (measured on 300 random windows each);
//  * np.sum and np.mean along axis 0 accumulate sequentially over the
//    submatrices; np.median is a selection, the mean of the two middle values
//    for an even count, NaN if any value is NaN;
//  * the contact pairs are written in the order of np.argsort, whose tie order
//    is numpy's own (core/include/hicx/numpy_sort.hpp).
//
// Threading and SIMD: none. The work after loading is a few hundred 31 by 31
// windows; the load and, with --transform, the obs/exp pass dominate.

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "hicx/adjust_ops.hpp"
#include "hicx/argparse.hpp"
#include "hicx/bins.hpp"
#include "hicx/clustering.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/matrix_ops.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/numpy_sort.hpp"
#include "hicx/obsexp_ops.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/transform_ops.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicAggregateContacts --matrix MATRIX --outFileName OUTFILENAME --BED BED\n"
    "                            --mode {inter-chr,intra-chr,all} [--range RANGE]\n"
    "                            [--row_wise] [--BED2 BED2] [--numberOfBins NUMBEROFBINS]\n"
    "                            [--transform {total-counts,z-score,obs/exp,none}]\n"
    "                            [--operationType {sum,mean,median}] [--perChr]\n"
    "                            [--considerStrandDirection]\n"
    "                            [--largeRegionsOperation {first,last,center}] [--help]\n"
    "                            [--version] [--dpi DPI]\n"
    "                            [--outFilePrefixMatrix OUTFILEPREFIXMATRIX]\n"
    "                            [--outFileContactPairs OUTFILECONTACTPAIRS]\n"
    "                            [--outFileObsExp OUTFILEOBSEXP]\n"
    "                            [--diagnosticHeatmapFile DIAGNOSTICHEATMAPFILE]\n"
    "                            [--kmeans KMEANS] [--hclust HCLUST] [--spectral SPECTRAL]\n"
    "                            [--howToCluster {full,center,diagonal}] [--keep_outlier]\n"
    "                            [--max_deviation MAX_DEVIATION]\n"
    "                            [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]\n"
    "                            [--colorMap COLORMAP] [--plotType {2d,3d}]\n"
    "                            [--vMin VMIN] [--vMax VMAX] [--noPlot]\n"
    "                            [--plotData FILE]\n";

const char* const kHelp =
    "\n"
    "Takes a list of positions in the Hi-C matrix and makes a pooled image.\n"
    "\n"
    "The options are those of the Python hicAggregateContacts. The tables are\n"
    "computed in C++, and the figures are drawn by the hicexplorer_plot drawing\n"
    "layer with the matplotlib calls of the Python tool (HICX_PLOT_PYTHON names the\n"
    "interpreter). C++ port only:\n"
    "  --noPlot              write the numeric outputs (--outFilePrefixMatrix,\n"
    "                        --outFileContactPairs, --outFileObsExp) without any\n"
    "                        figure.\n"
    "  --plotData FILE       write the data of the figures as JSON to FILE instead\n"
    "                        of drawing them.\n";

// A failure the Python reports through exit(message), log.error + exit(1), an
// assertion or an uncaught exception: all of them end with status 1.
struct ToolError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct Arguments {
    std::string matrix;
    std::string out_file_name;
    std::string bed;
    std::string mode;
    std::optional<std::string> range;
    bool row_wise = false;
    std::optional<std::string> bed2;
    std::int64_t number_of_bins = 51;
    std::string transform = "none";
    std::string operation_type = "median";
    bool per_chr = false;
    bool consider_strand_direction = false;
    std::string large_regions_operation = "first";
    std::optional<std::string> out_file_prefix_matrix;
    std::optional<std::string> out_file_contact_pairs;
    std::optional<std::string> out_file_obs_exp;
    std::optional<std::string> diagnostic_heatmap_file;
    std::optional<std::int64_t> kmeans;
    std::optional<std::int64_t> hclust;
    std::string how_to_cluster = "full";
    bool keep_outlier = false;
    std::int64_t max_deviation = 2;
    std::vector<std::string> chromosomes;
    std::string color_map = "RdYlBu_r";
    std::string plot_type = "2d";
    std::optional<double> v_min;
    std::optional<double> v_max;
    bool disable_bbox_tight = false;
    std::int64_t dpi = 300;
    bool no_plot = false;
    std::optional<std::string> plot_data;
};

std::optional<std::int64_t> python_int(const std::string& text) {
    std::size_t begin = 0;
    std::size_t end = text.size();
    while (begin < end && std::isspace(static_cast<unsigned char>(text[begin]))) {
        ++begin;
    }
    while (end > begin && std::isspace(static_cast<unsigned char>(text[end - 1]))) {
        --end;
    }
    std::string value = text.substr(begin, end - begin);
    std::string digits;
    for (std::size_t i = 0; i < value.size(); ++i) {
        const char c = value[i];
        if (c == '_' && i > 0 && i + 1 < value.size() &&
            std::isdigit(static_cast<unsigned char>(value[i - 1])) &&
            std::isdigit(static_cast<unsigned char>(value[i + 1]))) {
            continue;
        }
        digits.push_back(c);
    }
    if (digits.empty()) {
        return std::nullopt;
    }
    std::size_t start = (digits[0] == '+' || digits[0] == '-') ? 1 : 0;
    if (start == digits.size()) {
        return std::nullopt;
    }
    for (std::size_t i = start; i < digits.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(digits[i]))) {
            return std::nullopt;
        }
    }
    errno = 0;
    const long long parsed = std::strtoll(digits.c_str(), nullptr, 10);
    if (errno != 0) {
        return std::nullopt;
    }
    return static_cast<std::int64_t>(parsed);
}

std::int64_t to_int(const std::string& text) {
    const std::optional<std::int64_t> value = python_int(text);
    if (!value.has_value()) {
        throw ToolError("ValueError: invalid literal for int() with base 10: '" + text + "'");
    }
    return *value;
}

// hicAggregateContacts.py parse_arguments, plus the C++-only --noPlot. The
// plotting options are parsed and validated like the Python's; the numeric
// outputs do not use them. argparse.FileType('r') opens --BED and --BED2
// while parsing; FileType('w') is not opened here (no figure is written).
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicAggregateContacts",
                       "Takes a list of positions in the Hi-C matrix and makes a pooled image.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrix", "-m"})
        .required()
        .input({"h5", "cool", "mcool"})
        .help("Path of the Hi-C matrix to plot.");
    required.add({"--outFileName", "-out"})
        .file_type("w")
        .required()
        .output({"png", "pdf", "svg"})
        .help("File name to save the image.");
    required.add({"--BED"})
        .file_type("r")
        .required()
        .input({"bed"})
        .help("Interactions between regions in this BED file are plotted.");
    required.add({"--mode"}).choices({"inter-chr", "intra-chr", "all"}).required();

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--range"}).help(
        "Range of contacts considered for the aggregate contacts, as low_range:high_range in bp.");
    optional.add({"--row_wise"})
        .action(cli::Action::StoreTrue)
        .help("Compute the interactions between each row of the BED file and the same row of "
              "the BED2 file.");
    optional.add({"--BED2"})
        .file_type("r")
        .input({"bed"})
        .help("Optional second BED file.");
    optional.add({"--numberOfBins"})
        .type("int")
        .default_value("51")
        .help("Number of bins to include in the submatrix.");
    optional.add({"--transform"})
        .choices({"total-counts", "z-score", "obs/exp", "none"})
        .default_value("none")
        .help("Type of transformation for the matrix.");
    optional.add({"--operationType"})
        .choices({"sum", "mean", "median"})
        .default_value("median")
        .help("Operation that summarizes the submatrices into a single matrix.");
    optional.add({"--perChr"})
        .action(cli::Action::StoreTrue)
        .help("Generate a plot per chromosome (intra-chromosomal contacts only).");
    optional.add({"--considerStrandDirection"})
        .action(cli::Action::StoreTrue)
        .help("Take the strand into account: contacts of a reverse strand region are inverted.");
    optional.add({"--largeRegionsOperation"})
        .choices({"first", "last", "center"})
        .default_value("first")
        .help("Which bin of a region larger than a bin is used: first, last or center.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    cli::ArgumentGroup& output = parser.group("Output options");
    output.add({"--outFilePrefixMatrix"})
        .output({"tsv"}, "prefix")
        .help("Prefix of the tab-delimited tables with the values underlying the output matrix.");
    output.add({"--outFileContactPairs"})
        .output({"tsv"}, "prefix")
        .help("Prefix of the files with the contact positions used for the submatrices.");
    output.add({"--outFileObsExp"})
        .output({"h5", "cool"})
        .help("Writes the obs/exp matrix to a file, if --transform=obs/exp.");
    output.add({"--diagnosticHeatmapFile"})
        .file_type("w")
        .output({"png", "pdf", "svg"})
        .help("A heatmap (per chromosome) of the diagonals of the submatrices.");

    cli::ArgumentGroup& clustering = parser.group("Clustering options");
    clustering.add({"--kmeans"}).type("int").help("Number of clusters to compute with k-means.");
    clustering.add({"--hclust"})
        .type("int")
        .help("Number of clusters to compute (per chromosome) with hierarchical clustering.");
    clustering.add({"--spectral"})
        .type("int")
        .help("Number of clusters to compute (per chromosome) with spectral clustering.");
    clustering.add({"--howToCluster"})
        .choices({"full", "center", "diagonal"})
        .default_value("full")
        .help("Which values of each submatrix are clustered: full, center or diagonal.");
    clustering.add({"--keep_outlier"})
        .action(cli::Action::StoreTrue)
        .help("keep outliers before clustering.");
    clustering.add({"--max_deviation"})
        .type("int")
        .default_value(2)
        .help("max deviation from mean to be determined as outlier.");

    cli::ArgumentGroup& plotting = parser.group("Plotting options");
    plotting.add({"--chromosomes", "-C"}).nargs("+").help("List of chromosomes to plot.");
    plotting.add({"--colorMap"})
        .default_value("RdYlBu_r")
        .help("Color map to use for the heatmap.");
    plotting.add({"--plotType"}).choices({"2d", "3d"}).default_value("2d").help("Plot type.");
    plotting.add({"--vMin"}).type("float").help("Minimum value of the plotted score.");
    plotting.add({"--vMax"}).type("float").help("Maximum value of the plotted score.");
    plotting.add({"--disable_bbox_tight"}).action(cli::Action::StoreTrue).help(cli::kSuppress);
    plotting.add({"--noPlot"})
        .action(cli::Action::StoreTrue)
        .cpp_only("Writes the numeric outputs without drawing any figure.")
        .help("Write the numeric outputs without the figure.");
    plotting.add({"--plotData"})
        .metavar("FILE")
        .output({"json"})
        .cpp_only("The data of the figures as JSON, without drawing them (cpp/PLAN.md tier 7).")
        .help("Write the data the figures are drawn from as JSON to this file and do not draw "
              "them.");
    optional.add({"--dpi"})
        .type("int")
        .default_value(300)
        .help("Resolution for the image in case the output is a raster graphics image.");

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.matrix = ns.str("matrix");
    args.out_file_name = ns.str("outFileName");
    args.bed = ns.str("BED");
    args.mode = ns.str("mode");
    args.range = ns.opt_str("range");
    args.row_wise = ns.flag("row_wise");
    args.bed2 = ns.opt_str("BED2");
    args.number_of_bins = ns.integer("numberOfBins");
    args.transform = ns.str("transform");
    args.operation_type = ns.str("operationType");
    args.per_chr = ns.flag("perChr");
    args.consider_strand_direction = ns.flag("considerStrandDirection");
    args.large_regions_operation = ns.str("largeRegionsOperation");
    args.out_file_prefix_matrix = ns.opt_str("outFilePrefixMatrix");
    args.out_file_contact_pairs = ns.opt_str("outFileContactPairs");
    args.out_file_obs_exp = ns.opt_str("outFileObsExp");
    args.diagnostic_heatmap_file = ns.opt_str("diagnosticHeatmapFile");
    args.kmeans = ns.opt_integer("kmeans");
    args.hclust = ns.opt_integer("hclust");
    args.how_to_cluster = ns.str("howToCluster");
    args.keep_outlier = ns.flag("keep_outlier");
    args.max_deviation = ns.integer("max_deviation");
    args.chromosomes = ns.strs("chromosomes");
    args.color_map = ns.str("colorMap");
    args.plot_type = ns.str("plotType");
    args.v_min = ns.opt_real("vMin");
    args.v_max = ns.opt_real("vMax");
    args.disable_bbox_tight = ns.flag("disable_bbox_tight");
    args.dpi = ns.integer("dpi");
    args.no_plot = ns.flag("noPlot");
    args.plot_data = ns.opt_str("plotData");
    return args;
}

std::vector<std::string> split_whitespace(const std::string& line) {
    std::vector<std::string> fields;
    std::istringstream stream(line);
    std::string field;
    while (stream >> field) {
        fields.push_back(field);
    }
    return fields;
}

std::vector<std::string> read_lines(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        throw ToolError("cannot open '" + path + "'");
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line)) {
        lines.push_back(line);
    }
    return lines;
}

// utilities.change_chrom_names
std::string change_chrom_names(const std::string& chrom) {
    if (chrom.rfind("chr", 0) == 0) {
        return chrom.substr(3);
    }
    return "chr" + chrom;
}

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// --------------------------------------------------------------------------
// the matrix, as count_contacts reads it

class MatrixView {
  public:
    explicit MatrixView(const hicx::CsrMatrix& matrix) : matrix_(matrix) {
        const std::vector<std::int64_t>& indptr = matrix.indptr();
        const std::vector<std::int32_t>& indices = matrix.indices();
        for (std::int64_t row = 0; row < matrix.rows() && sorted_; ++row) {
            for (std::int64_t k = indptr[static_cast<std::size_t>(row)] + 1;
                 k < indptr[static_cast<std::size_t>(row) + 1]; ++k) {
                if (indices[static_cast<std::size_t>(k)] <
                    indices[static_cast<std::size_t>(k) - 1]) {
                    sorted_ = false;
                    break;
                }
            }
        }
    }

    // matrix[row, col], summing duplicates as scipy's indexing does.
    [[nodiscard]] double at(std::int64_t row, std::int64_t col) const {
        if (matrix_.symmetry() == hicx::Symmetry::UpperTriangle && row > col) {
            std::swap(row, col);
        }
        const auto begin = matrix_.indices().begin() +
                           matrix_.indptr()[static_cast<std::size_t>(row)];
        const auto end = matrix_.indices().begin() +
                         matrix_.indptr()[static_cast<std::size_t>(row) + 1];
        double value = 0.0;
        if (sorted_) {
            for (auto it = std::lower_bound(begin, end, static_cast<std::int32_t>(col));
                 it != end && *it == col; ++it) {
                value += matrix_.data()[static_cast<std::size_t>(it - matrix_.indices().begin())];
            }
        } else {
            for (auto it = begin; it != end; ++it) {
                if (*it == col) {
                    value += matrix_.data()[static_cast<std::size_t>(it - matrix_.indices().begin())];
                }
            }
        }
        return value;
    }

  private:
    const hicx::CsrMatrix& matrix_;
    bool sorted_ = true;
};

// --------------------------------------------------------------------------
// aggregation state, agg_info in the Python

struct Group {
    std::string chrom2;
    std::vector<std::vector<double>> submatrices;  // oriented, row major
    std::vector<double> centers;
    std::vector<std::array<std::string, 4>> positions;  // start1, end1, start2, end2
};

struct ChromEntry {
    std::string chrom1;
    std::vector<Group> groups;  // in first seen order

    Group* find(const std::string& chrom2) {
        for (Group& group : groups) {
            if (group.chrom2 == chrom2) {
                return &group;
            }
        }
        return nullptr;
    }
};

struct Region {
    std::string chrom;
    std::string start;
    std::string end;
    std::optional<std::string> strand;
};

struct Aggregator {
    const Arguments& args;
    const hicx::BinTable& bins;
    const MatrixView& view;
    std::int64_t bin_size;
    std::int64_t m_half;
    std::string range;
    std::vector<ChromEntry> entries;
    std::unordered_map<std::string, std::size_t> entry_index;
    std::unordered_map<std::string, std::pair<std::int64_t, std::int64_t>> chrom_coord;
    std::set<std::pair<std::int64_t, std::int64_t>> seen;
    std::set<std::pair<std::string, std::string>> seen_chrs;
    std::int64_t counter = 0;
    std::int64_t used_counter = 0;
    std::int64_t empty_mat = 0;

    // count_contacts (hicAggregateContacts.py:348-462)
    void count_contacts(Region first, Region second) {
        if (chrom_coord.count(first.chrom) == 0 || chrom_coord.count(second.chrom) == 0) {
            return;
        }
        const auto& coord1 = chrom_coord.at(first.chrom);
        const auto& coord2 = chrom_coord.at(second.chrom);
        if (to_int(first.end) > coord1.second || to_int(second.end) > coord2.second) {
            return;
        }
        if (to_int(first.start) < coord1.first || to_int(second.start) < coord2.first) {
            return;
        }
        const auto range1 =
            bins.region_bin_range(first.chrom, to_int(first.start), to_int(first.end));
        const auto range2 =
            bins.region_bin_range(second.chrom, to_int(second.start), to_int(second.end));
        if (range1 == range2) {
            return;  // includes None == None
        }
        if (!range1.has_value() || !range2.has_value()) {
            return;
        }
        const auto pick = [this](const std::pair<std::int64_t, std::int64_t>& range) {
            if (args.large_regions_operation == "first") {
                return range.first;
            }
            if (args.large_regions_operation == "last") {
                return range.second;
            }
            // int(np.floor(np.mean(bin_id)))
            const double mean =
                (static_cast<double>(range.first) + static_cast<double>(range.second)) / 2.0;
            return static_cast<std::int64_t>(std::floor(mean));
        };
        std::int64_t bin1 = pick(*range1);
        std::int64_t bin2 = pick(*range2);
        const std::optional<std::string> orientation1 = first.strand;
        const std::optional<std::string> orientation2 = second.strand;
        if (bin1 > bin2) {
            if (first.chrom == second.chrom) {
                std::swap(bin1, bin2);
            } else {
                // chromosome, start, end and bin move; the orientations do not
                std::swap(first.chrom, second.chrom);
                std::swap(first.start, second.start);
                std::swap(first.end, second.end);
                std::swap(bin1, bin2);
            }
        }
        ++counter;
        const std::string& chrom1 = first.chrom;
        const std::string& chrom2 = second.chrom;
        const hicx::BinRange chrom1_range = *bins.chrom_bin_range(chrom1);
        const hicx::BinRange chrom2_range = *bins.chrom_bin_range(chrom2);

        if (seen_chrs.insert({chrom1, chrom2}).second) {
            ChromEntry& entry = entries[entry_index.at(chrom1)];
            if (Group* existing = entry.find(chrom2)) {
                existing->submatrices.clear();
                existing->centers.clear();
                existing->positions.clear();
            } else {
                entry.groups.push_back(Group{chrom2, {}, {}, {}});
            }
        }

        if (args.mode == "intra-chr") {
            const std::size_t colon = range.find(':');
            const std::int64_t min_bins = floor_div(to_int(range.substr(0, colon)), bin_size);
            const std::int64_t max_bins = floor_div(to_int(range.substr(colon + 1)), bin_size);
            const std::int64_t distance = std::llabs(bin2 - bin1);
            if (min_bins > distance || distance > max_bins) {
                return;
            }
        }
        if (!seen.insert({bin1, bin2}).second) {
            return;
        }
        if (bin1 - m_half < chrom1_range.first || bin1 + m_half >= chrom1_range.last) {
            return;
        }
        if (bin2 - m_half < chrom2_range.first || bin2 + m_half >= chrom2_range.last) {
            return;
        }

        const std::int64_t width = 2 * m_half + 1;
        if (width <= 0) {
            ++empty_mat;  // an empty slice sums to zero
            return;
        }
        const auto w = static_cast<std::size_t>(width);
        std::vector<double> window(w * w);
        for (std::int64_t i = 0; i < width; ++i) {
            for (std::int64_t j = 0; j < width; ++j) {
                window[static_cast<std::size_t>(i) * w + static_cast<std::size_t>(j)] =
                    view.at(bin1 - m_half + i, bin2 - m_half + j);
            }
        }
        enum class Orientation { None, FlipLR, FlipUD, Transpose };
        Orientation orientation = Orientation::None;
        if (!(orientation1.has_value() || orientation2.has_value()) ||
            (orientation1 == std::string("+") && orientation2 == std::string("+"))) {
            orientation = Orientation::None;
        } else if (orientation1 == std::string("+") && orientation2 == std::string("-")) {
            orientation = Orientation::FlipLR;
        } else if (orientation1 == std::string("-") && orientation2 == std::string("+")) {
            orientation = Orientation::FlipUD;
        } else if (orientation1 == std::string("-") && orientation2 == std::string("-")) {
            orientation = Orientation::Transpose;
        }
        std::vector<double> oriented(w * w);
        for (std::size_t i = 0; i < w; ++i) {
            for (std::size_t j = 0; j < w; ++j) {
                std::size_t source = i * w + j;
                switch (orientation) {
                    case Orientation::None: break;
                    case Orientation::FlipLR: source = i * w + (w - 1 - j); break;
                    case Orientation::FlipUD: source = (w - 1 - i) * w + j; break;
                    case Orientation::Transpose: source = j * w + i; break;
                }
                oriented[i * w + j] = window[source];
            }
        }
        // np.matrix.sum(): logical order for the flipped views, memory order
        // for the unflipped and the transposed one.
        const double sum =
            (orientation == Orientation::FlipLR || orientation == Orientation::FlipUD)
                ? hicx::npy::pairwise_sum(oriented.data(), oriented.size())
                : hicx::npy::pairwise_sum(window.data(), window.size());
        if (sum == 0.0) {
            ++empty_mat;
            return;
        }
        ++used_counter;
        if (args.transform == "total-counts" && sum > 0.0) {
            for (double& value : oriented) {
                value = value / sum;
            }
        }
        Group& group = *entries[entry_index.at(chrom1)].find(chrom2);
        group.submatrices.push_back(std::move(oriented));
        group.centers.push_back(view.at(bin1, bin2));
        group.positions.push_back({first.start, first.end, second.start, second.end});
    }

    static std::int64_t floor_div(std::int64_t a, std::int64_t b) {
        std::int64_t quotient = a / b;
        if ((a % b != 0) && ((a < 0) != (b < 0))) {
            --quotient;
        }
        return quotient;
    }
};

// read_bed_per_chrom (:239-266)
std::vector<std::pair<std::string, std::vector<Region>>> read_bed_per_chrom(
    const std::string& path, const std::set<std::string>& chrom_list, bool strand) {
    std::vector<std::pair<std::string, std::vector<Region>>> intervals;
    for (const std::string& raw : read_lines(path)) {
        std::string line = raw;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (!line.empty() && line[0] == '#') {
            continue;
        }
        std::vector<std::string> fields = split_whitespace(line);
        if (strand && fields.size() < 6) {
            throw ToolError("Strand information should be considered, but BED file has not at "
                            "least six columns. Exiting!");
        }
        if (fields.empty()) {
            throw ToolError("IndexError: list index out of range (empty line in " + path + ")");
        }
        if (chrom_list.count(fields[0]) == 0) {
            if (chrom_list.count(change_chrom_names(fields[0])) != 0) {
                fields[0] = change_chrom_names(fields[0]);
            } else {
                continue;
            }
        }
        if (fields.size() < 3) {
            throw ToolError("IndexError: list index out of range (" + path + ")");
        }
        Region region;
        region.chrom = fields[0];
        region.start = std::to_string(to_int(fields[1]));
        region.end = std::to_string(to_int(fields[2]));
        if (strand) {
            region.strand = fields[5];
        }
        auto it = std::find_if(intervals.begin(), intervals.end(),
                               [&](const auto& item) { return item.first == region.chrom; });
        if (it == intervals.end()) {
            intervals.emplace_back(region.chrom, std::vector<Region>{});
            it = intervals.end() - 1;
        }
        it->second.push_back(std::move(region));
    }
    return intervals;
}

bool same_coordinates(const Region& a, const Region& b) {
    return a.start == b.start && a.end == b.end && a.strand == b.strand;
}

// --------------------------------------------------------------------------
// clustering and output

struct ClusterInput {
    std::string name;
    std::vector<const std::vector<double>*> submatrices;
    std::vector<double> centers;
    std::vector<std::array<std::string, 6>> coords;
    std::vector<std::vector<std::int64_t>> clusters;
};

// np.median over a copy of the values: NaN if any value is NaN.
double median(std::vector<double> values) {
    for (const double value : values) {
        if (std::isnan(value)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    const std::size_t n = values.size();
    const std::size_t half = n / 2;
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(half),
                     values.end());
    const double upper = values[half];
    if (n % 2 == 1) {
        return upper;
    }
    const double lower = *std::max_element(values.begin(),
                                           values.begin() + static_cast<std::ptrdiff_t>(half));
    return (lower + upper) / 2.0;
}

// compute_clusters (:500-572)
void compute_clusters(ClusterInput& info, std::int64_t k, const std::string& method,
                      const std::string& how, std::int64_t max_deviation, bool keep_outlier) {
    if (info.submatrices.empty()) {
        throw ToolError("IndexError: list index out of range (no submatrix to cluster)");
    }
    const std::size_t cells = info.submatrices[0]->size();
    const auto w = static_cast<std::size_t>(std::llround(std::sqrt(static_cast<double>(cells))));
    std::size_t features = cells;
    if (how == "diagonal") {
        features = w;
    } else if (how == "center") {
        features = 1;
        if (w < 3) {
            throw ToolError("ValueError: cannot reshape the centre of a " + std::to_string(w) +
                            " by " + std::to_string(w) + " submatrix into shape (1,9)");
        }
    }
    const std::size_t n = info.submatrices.size();
    std::vector<double> vectors(n * features);
    for (std::size_t s = 0; s < n; ++s) {
        const std::vector<double>& m = *info.submatrices[s];
        if (how == "diagonal") {
            for (std::size_t i = 0; i < w; ++i) {
                vectors[s * features + i] = m[i * w + i];
            }
        } else if (how == "center") {
            const std::size_t center_bin = (w + 1) / 2;
            std::array<double, 9> block{};
            std::size_t index = 0;
            for (std::size_t i = center_bin - 2; i <= center_bin; ++i) {
                for (std::size_t j = center_bin - 2; j <= center_bin; ++j) {
                    block[index++] = m[i * w + j];
                }
            }
            vectors[s] = hicx::npy::pairwise_sum(block.data(), 9) / 9.0;
        } else {
            std::copy(m.begin(), m.end(), vectors.begin() + static_cast<std::ptrdiff_t>(s * features));
        }
    }

    std::vector<std::size_t> kept(n);
    for (std::size_t s = 0; s < n; ++s) {
        kept[s] = s;
    }
    if (!keep_outlier) {
        // get_outlier_indices, defect 4
        const double data_median = median(vectors);
        std::vector<double> absolute(vectors.size());
        for (std::size_t i = 0; i < vectors.size(); ++i) {
            absolute[i] = std::abs(vectors[i]);
        }
        const double mad = 1.4826 * median(absolute);
        if (!(mad == 0.0)) {
            std::vector<std::size_t> survivors;
            for (std::size_t s = 0; s < n; ++s) {
                bool outlier = false;
                for (std::size_t f = 0; f < features; ++f) {
                    const double deviation =
                        std::abs(vectors[s * features + f] - data_median) / mad;
                    if (deviation > static_cast<double>(max_deviation)) {
                        outlier = true;
                        break;
                    }
                }
                if (!outlier) {
                    survivors.push_back(s);
                }
            }
            if (survivors.size() != n) {
                if (survivors.empty()) {
                    throw ToolError("ERROR: all submatrices have been detected as outliers. "
                                    "You can consider changing the threshold");
                }
                kept = std::move(survivors);
            }
        }
    }

    hicx::cluster::Samples samples;
    samples.samples = static_cast<std::int64_t>(kept.size());
    samples.features = static_cast<std::int64_t>(features);
    samples.values.reserve(kept.size() * features);
    ClusterInput updated;
    updated.name = info.name;
    for (const std::size_t s : kept) {
        for (std::size_t f = 0; f < features; ++f) {
            const double value = vectors[s * features + f];
            samples.values.push_back(std::isnan(value) ? 0.0 : value);
        }
        updated.submatrices.push_back(info.submatrices[s]);
        updated.centers.push_back(info.centers[s]);
        updated.coords.push_back(info.coords[s]);
    }

    std::vector<std::int64_t> labels;
    try {
        if (method == "kmeans") {
            labels = hicx::cluster::kmeans(samples, k, 0);
        } else if (method == "hierarchical") {
            labels = hicx::cluster::ward(samples, k);
        } else {
            labels.assign(kept.size(), 0);
        }
    } catch (const hicx::cluster::ClusteringError& error) {
        throw ToolError(std::string("ValueError: ") + error.what());
    }
    for (std::int64_t cluster = 0; cluster < k; ++cluster) {
        std::vector<std::int64_t> members;
        for (std::size_t i = 0; i < labels.size(); ++i) {
            if (labels[i] == cluster) {
                members.push_back(static_cast<std::int64_t>(i));
            }
        }
        updated.clusters.push_back(std::move(members));
    }
    info = std::move(updated);
}

std::string format_fixed5(double value) {
    if (std::isnan(value)) {
        return "nan";
    }
    if (std::isinf(value)) {
        return value < 0 ? "-inf" : "inf";
    }
    char buffer[512];
    std::snprintf(buffer, sizeof(buffer), "%0.5f", value);
    return buffer;
}

void write_text(const std::string& path, const std::string& content) {
    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    if (!file) {
        throw ToolError("[Errno " + std::to_string(errno) + "] " + std::strerror(errno) + ": '" +
                        path + "'");
    }
    file << content;
    file.close();
    if (!file) {
        throw ToolError("could not write '" + path + "'");
    }
}

void write_obs_exp(const std::string& path, const hicx::MatrixData& data) {
    // hicAggregateContacts.py:857-866 builds its own MatrixFileHandler from the
    // output name, as hicPCA does for its intermediate matrices.
    hicx::MatrixData out;
    out.matrix = data.matrix;
    out.cut_intervals = data.cut_intervals;
    out.nan_bins = data.nan_bins;
    out.correction_factors = data.correction_factors;
    out.distance_counts = data.distance_counts;
    out.correction_factors_are_column = data.correction_factors_are_column;
    if (ends_with(path, ".h5")) {
        hicx::H5SaveOptions options;
        options.symmetric = true;
        hicx::write_hicexplorer_h5(path, out, options);
        return;
    }
    hicx::CoolSaveOptions options;
    options.symmetric = true;
    options.apply_correction = false;
    hicx::write_cool(path, out, options);
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);

    namespace plot = hicx::plot;
    std::string plot_json;
    try {
        hicx::ToolMatrix hic = hicx::ToolMatrix::load(args.matrix);
        hicx::MatrixData& data = hic.data();

        // ma.maskBins(ma.nan_bins); ma.matrix.data[isnan] = 0; the second
        // maskBins is a no operation because the first emptied nan_bins.
        const std::vector<std::int64_t> nan_bins = data.nan_bins;
        const auto zero_nan = [&data]() {
            for (double& value : data.matrix.mutable_data()) {
                if (std::isnan(value)) {
                    value = 0.0;
                }
            }
        };
        if (!args.chromosomes.empty() && !nan_bins.empty()) {
            // Defect 7: keepOnlyTheseChr restores what maskBins removed and the
            // bin table from before enlarge_bins.
            hicx::mask_and_restore_bins(data, nan_bins);
            zero_nan();
        } else {
            hicx::delete_bins(data, nan_bins);
            zero_nan();
            hicx::enlarge_bins(data.cut_intervals);
        }
        if (!args.chromosomes.empty()) {
            std::vector<std::string> names;
            for (const std::string& name : args.chromosomes) {
                names.push_back(name);
            }
            try {
                hicx::keep_only_chromosomes(data, names);
            } catch (const std::runtime_error& error) {
                throw ToolError(std::string("ValueError: ") + error.what());
            }
        }
        hic.refresh_boundaries();
        const hicx::BinTable bins(data.cut_intervals);
        const std::int64_t bin_size = bins.bin_size();

        const std::string default_range = "1000000:20000000";
        const std::string range = args.range.value_or(default_range);
        const std::size_t colon = range.find(':');
        if (colon == std::string::npos || range.find(':', colon + 1) != std::string::npos) {
            throw ToolError("ValueError: --range must have the form low:high");
        }
        const std::string min_dist = range.substr(0, colon);
        const std::string max_dist = range.substr(colon + 1);
        if (args.mode == "intra-chr") {
            if (!(to_int(min_dist) < to_int(max_dist))) {
                throw ToolError("AssertionError: Error lower range is larger than upper range!");
            }
        }

        if (args.transform == "z-score" || args.transform == "obs/exp") {
            hicx::ObsExpOptions options;
            options.zscore = args.transform == "z-score";
            options.perchr = true;
            if (args.mode == "intra-chr") {
                options.max_depth_bp = static_cast<double>(to_int(max_dist)) * 2.5;
                options.unbounded = options.max_depth_bp == 0.0;  // `if maxdepth:`
            } else {
                options.unbounded = true;
            }
            try {
                hicx::convert_to_obs_exp_matrix(data, bin_size, options);
            } catch (const std::runtime_error& error) {
                throw ToolError(std::string("ValueError: ") + error.what());
            }
            if (args.transform == "obs/exp" && args.out_file_obs_exp.has_value()) {
                write_obs_exp(*args.out_file_obs_exp, data);
            }
        }

        const std::int64_t m_bins =
            args.number_of_bins % 2 == 1 || args.number_of_bins % 2 == -1
                ? args.number_of_bins
                : args.number_of_bins + 1;
        // int((M - 1) // 2), floor division
        const std::int64_t m_half = Aggregator::floor_div(m_bins - 1, 2);

        const MatrixView view(data.matrix);
        Aggregator aggregator{args, bins, view, bin_size, m_half, range, {}, {}, {}, {}, {}};
        std::set<std::string> chrom_set;
        for (const auto& [chrom, bin_range] : bins.chrom_bin_boundaries()) {
            chrom_set.insert(chrom);
            aggregator.entry_index.emplace(chrom, aggregator.entries.size());
            aggregator.entries.push_back(ChromEntry{chrom, {}});
            aggregator.chrom_coord[chrom] = {bins.bin_pos(static_cast<std::size_t>(bin_range.first)).start,
                                             bins.bin_pos(static_cast<std::size_t>(bin_range.last - 1)).end};
        }

        if (args.mode == "inter-chr" && aggregator.chrom_coord.size() == 1) {
            throw ToolError("Error: 'inter-chr' mode can not be applied on matrices of only one "
                            "chromosme.");
        }
        if (args.mode == "inter-chr" && args.per_chr) {
            throw ToolError("Error: 'inter-chr' mode can not be used along with --perChr.");
        }
        if (args.mode == "all" && args.per_chr) {
            throw ToolError("Error: 'all' mode can not be used along with --perChr.");
        }

        if (args.row_wise) {
            const std::vector<std::string> lines1 = read_lines(args.bed);
            if (!args.bed2.has_value()) {
                throw ToolError("Error computing row-wise contacts requires two bed files!");
            }
            const std::vector<std::string> lines2 = read_lines(*args.bed2);
            if (lines1.size() != lines2.size()) {
                throw ToolError("Error row_wise only works if both bed files have the same length.");
            }
            for (std::size_t r = 0; r < lines1.size(); ++r) {
                std::vector<std::string> line1 = split_whitespace(lines1[r]);
                std::vector<std::string> line2 = split_whitespace(lines2[r]);
                if (line1.empty() || line2.empty()) {
                    throw ToolError("IndexError: list index out of range (row-wise BED line " +
                                    std::to_string(r + 1) + ")");
                }
                if (chrom_set.count(line1[0]) == 0) {
                    line1[0] = change_chrom_names(line1[0]);
                    if (chrom_set.count(line1[0]) == 0) {
                        continue;
                    }
                }
                if (chrom_set.count(line2[0]) == 0) {
                    line2[0] = change_chrom_names(line2[0]);
                    if (chrom_set.count(line2[0]) == 0) {
                        continue;
                    }
                }
                if (args.mode == "inter-chr" && line1[0] == line2[0]) {
                    continue;
                }
                if (args.mode == "intra-chr" && line1[0] != line2[0]) {
                    continue;
                }
                const std::size_t needed = args.consider_strand_direction ? 6 : 3;
                if (line1.size() < needed || line2.size() < needed) {
                    throw ToolError("IndexError: list index out of range (row-wise BED line " +
                                    std::to_string(r + 1) + ")");
                }
                Region first{line1[0], line1[1], line1[2], std::nullopt};
                Region second{line2[0], line2[1], line2[2], std::nullopt};
                if (args.consider_strand_direction) {
                    first.strand = line1[5];
                    second.strand = line2[5];
                }
                aggregator.count_contacts(first, second);
            }
            // aggregate_contacts_per_row cleanup (:336-345)
        } else {
            const auto bed1 =
                read_bed_per_chrom(args.bed, chrom_set, args.consider_strand_direction);
            const auto bed2 = args.bed2.has_value()
                                  ? read_bed_per_chrom(*args.bed2, chrom_set,
                                                       args.consider_strand_direction)
                                  : bed1;
            for (const auto& [k1, v1] : bed1) {
                for (const auto& [k2, v2] : bed2) {
                    if (args.mode == "inter-chr" && k1 == k2) {
                        if (bed1.size() == 1 && bed2.size() == 1) {
                            throw ToolError("Error: 'inter-chr' mode needs at least a pair of "
                                            "coordinates with different chromoses to be "
                                            "available in the bed files.");
                        }
                        continue;
                    }
                    if (args.mode == "intra-chr" && k1 != k2) {
                        continue;
                    }
                    for (const Region& coord1 : v1) {
                        for (const Region& coord2 : v2) {
                            if (k1 == k2 && same_coordinates(coord1, coord2)) {
                                continue;
                            }
                            aggregator.count_contacts(coord1, coord2);
                        }
                    }
                }
            }
        }

        // The cleanup of both aggregate functions: chromosomes without any
        // submatrix go. Empty pair groups are dropped here as well; the
        // row-wise cleanup keeps them in `all` and `inter-chr` mode, where
        // they only ever feed the genome-wide concatenation and add nothing.
        std::vector<ChromEntry*> kept_entries;
        std::size_t total_submatrices = 0;
        for (ChromEntry& entry : aggregator.entries) {
            bool keep = false;
            bool has_group = !entry.groups.empty();
            for (const Group& group : entry.groups) {
                keep = keep || !group.submatrices.empty();
                total_submatrices += group.submatrices.size();
            }
            if (args.row_wise && has_group && args.mode != "intra-chr") {
                keep = true;
            }
            if (keep) {
                kept_entries.push_back(&entry);
            }
        }
        if (kept_entries.empty()) {
            throw ToolError("No susbmatrix found to be aggregated.");
        }

        // cluster_matrices (:575-625)
        std::string method = "no_clust";
        std::int64_t k = 1;
        std::string how = "full";
        if (args.kmeans.has_value()) {
            if (!(*args.kmeans > 1)) {
                throw ToolError("AssertionError: --kmeans must be larger than 1");
            }
            method = "kmeans";
            k = *args.kmeans;
            how = args.how_to_cluster;
        } else if (args.hclust.has_value()) {
            if (!(*args.hclust > 1)) {
                throw ToolError("AssertionError: --hclust must be larger than 1");
            }
            method = "hierarchical";
            k = *args.hclust;
            how = args.how_to_cluster;
        }
        const std::int64_t num_clusters = k;

        std::vector<ClusterInput> clustered;
        if (args.per_chr) {
            for (ChromEntry* entry : kept_entries) {
                for (const Group& group : entry->groups) {
                    if (group.submatrices.empty()) {
                        continue;
                    }
                    ClusterInput info;
                    info.name = entry->chrom1;
                    for (std::size_t s = 0; s < group.submatrices.size(); ++s) {
                        info.submatrices.push_back(&group.submatrices[s]);
                        info.centers.push_back(group.centers[s]);
                        const auto& p = group.positions[s];
                        info.coords.push_back({entry->chrom1, p[0], p[1], group.chrom2, p[2], p[3]});
                    }
                    if (static_cast<std::int64_t>(group.submatrices.size()) < k) {
                        k = 1;  // defect 6: persists for the chromosomes after it
                    }
                    compute_clusters(info, k, method, how, args.max_deviation, args.keep_outlier);
                    clustered.push_back(std::move(info));
                }
            }
        } else {
            ClusterInput info;
            info.name = "genome";
            for (ChromEntry* entry : kept_entries) {
                for (const Group& group : entry->groups) {
                    for (std::size_t s = 0; s < group.submatrices.size(); ++s) {
                        info.submatrices.push_back(&group.submatrices[s]);
                        info.centers.push_back(group.centers[s]);
                        const auto& p = group.positions[s];
                        info.coords.push_back({entry->chrom1, p[0], p[1], group.chrom2, p[2], p[3]});
                    }
                }
            }
            compute_clusters(info, k, method, how, args.max_deviation, args.keep_outlier);
            clustered.push_back(std::move(info));
        }

        for (const ClusterInput& info : clustered) {
            for (const auto& members : info.clusters) {
                if (members.empty()) {
                    throw ToolError("IndexError: tuple index out of range (a cluster of " +
                                    info.name + " is empty)");
                }
            }
        }

        const bool integral_centers = data.matrix.integral_dtype();
        std::vector<std::string> chroms_json;
        for (const ClusterInput& info : clustered) {
            const std::size_t cells = info.submatrices[0]->size();
            const auto w = static_cast<std::size_t>(
                std::llround(std::sqrt(static_cast<double>(cells))));
            std::vector<std::string> clusters_json;
            for (std::size_t cluster = 0; cluster < info.clusters.size(); ++cluster) {
                const std::vector<std::int64_t>& members = info.clusters[cluster];
                const std::string suffix =
                    num_clusters == 1 ? "_" + info.name + ".tab"
                                      : "_" + info.name + "_cluster_" +
                                            std::to_string(cluster + 1) + ".tab";
                // compute_avg (:628-644), for the table and the figure alike.
                std::vector<double> average(cells, 0.0);
                if (args.operation_type == "median") {
                    std::vector<double> column(members.size());
                    for (std::size_t c = 0; c < cells; ++c) {
                        for (std::size_t m = 0; m < members.size(); ++m) {
                            column[m] = (*info.submatrices[static_cast<std::size_t>(members[m])])[c];
                        }
                        average[c] = median(column);
                    }
                } else {
                    average = *info.submatrices[static_cast<std::size_t>(members[0])];
                    for (std::size_t m = 1; m < members.size(); ++m) {
                        const std::vector<double>& sub =
                            *info.submatrices[static_cast<std::size_t>(members[m])];
                        for (std::size_t c = 0; c < cells; ++c) {
                            average[c] = average[c] + sub[c];
                        }
                    }
                    if (args.operation_type == "mean") {
                        for (double& value : average) {
                            value = value / static_cast<double>(members.size());
                        }
                    }
                }
                if (!args.no_plot) {
                    std::vector<std::string> rows;
                    for (std::size_t i = 0; i < w; ++i) {
                        rows.push_back(plot::json_numbers(std::vector<double>(
                            average.begin() + static_cast<std::ptrdiff_t>(i * w),
                            average.begin() + static_cast<std::ptrdiff_t>((i + 1) * w))));
                    }
                    plot::JsonObject cluster_json;
                    cluster_json.add("average", plot::json_list(rows));
                    cluster_json.add("indices", plot::json_ints(members));
                    clusters_json.push_back(cluster_json.str());
                }
                if (args.out_file_prefix_matrix.has_value()) {
                    std::string text;
                    for (std::size_t i = 0; i < w; ++i) {
                        for (std::size_t j = 0; j < w; ++j) {
                            if (j != 0) {
                                text.push_back('\t');
                            }
                            text += format_fixed5(average[i * w + j]);
                        }
                        text.push_back('\n');
                    }
                    write_text(*args.out_file_prefix_matrix + suffix, text);
                }
                if (args.out_file_contact_pairs.has_value()) {
                    std::string text;
                    std::vector<std::int64_t> order;
                    std::vector<double> values(members.size());
                    for (std::size_t m = 0; m < members.size(); ++m) {
                        values[m] = info.centers[static_cast<std::size_t>(members[m])];
                    }
                    if (integral_centers) {
                        std::vector<std::int64_t> integers(values.size());
                        for (std::size_t m = 0; m < values.size(); ++m) {
                            integers[m] = static_cast<std::int64_t>(values[m]);
                        }
                        order = hicx::npy::argsort(integers);
                    } else {
                        order = hicx::npy::argsort(values);
                    }
                    for (auto it = order.rbegin(); it != order.rend(); ++it) {
                        const auto cl_idx = static_cast<std::size_t>(*it);
                        const double value = values[cl_idx];
                        const auto& coord = info.coords[cl_idx];  // defect 2
                        for (const std::string& field : coord) {
                            text += field;
                            text.push_back('\t');
                        }
                        text += integral_centers
                                    ? std::to_string(static_cast<std::int64_t>(value))
                                    : hicx::npy::float_repr(value);
                        text.push_back('\n');
                    }
                    write_text(*args.out_file_contact_pairs + suffix, text);
                }
            }
            if (!args.no_plot) {
                // clustered_info[chrom]['diagonal']: mat_to_append.diagonal() of
                // every kept submatrix, only read by plot_diagnostic_heatmaps.
                std::vector<std::string> diagonals;
                if (args.diagnostic_heatmap_file.has_value()) {
                    for (const std::vector<double>* submatrix : info.submatrices) {
                        std::vector<double> diagonal(w);
                        for (std::size_t i = 0; i < w; ++i) {
                            diagonal[i] = (*submatrix)[i * w + i];
                        }
                        diagonals.push_back(plot::json_numbers(diagonal));
                    }
                }
                plot::JsonObject chrom_json;
                chrom_json.add("name", plot::json_string(info.name));
                chrom_json.add("clusters", plot::json_list(clusters_json));
                chrom_json.add("diagonal", plot::json_list(diagonals));
                chroms_json.push_back(chrom_json.str());
            }
        }

        if (!args.no_plot) {
            plot::JsonObject figure;
            figure.add("outFileName", plot::json_string(args.out_file_name));
            figure.add("diagnosticHeatmapFile",
                       args.diagnostic_heatmap_file.has_value()
                           ? plot::json_string(*args.diagnostic_heatmap_file)
                           : "null");
            figure.add("dpi", plot::json_int(args.dpi));
            figure.add("vMin", args.v_min.has_value() ? plot::json_number(*args.v_min) : "null");
            figure.add("vMax", args.v_max.has_value() ? plot::json_number(*args.v_max) : "null");
            figure.add("colorMap", plot::json_string(args.color_map));
            figure.add("plotType", plot::json_string(args.plot_type));
            figure.add("disable_bbox_tight", plot::json_bool(args.disable_bbox_tight));
            figure.add("M_half", plot::json_int(m_half));
            figure.add("num_clusters", plot::json_int(num_clusters));
            figure.add("chroms", plot::json_list(chroms_json));
            plot_json = figure.str();
        }
    } catch (const ToolError& error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicAggregateContacts: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicAggregateContacts");
    if (args.no_plot) {
        return 0;
    }
    return plot::draw("hicAggregateContacts", plot_json, args.plot_data);
}
