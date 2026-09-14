// Port of hicexplorer/hicPlotDistVsCounts.py (cpp/PLAN.md tier 7, option (a)).
//
// The distance versus counts reduction runs here: the mean contact per genomic
// distance, per matrix and per chromosome or domain, the scale factors between
// matrices and the --outFileData table. The figure is drawn by
// plot/hicexplorer_plot/hicPlotDistVsCounts.py with the matplotlib calls of
// the Python tool, from the curves this binary computes (hicx::plot::draw).
//
// What is reproduced, pinned by the harness cases:
//
//  1. compute_distance_mean (:117-271). With --maxdepth (any non-zero value)
//     the upper triangle is cut at int(maxdepth * 1.5 / binsize) bins, the
//     bins' starts are snapped by HiCMatrix.fit_cut_intervals when more than
//     1 % of the bins deviate from the median bin size, and every stored
//     value is added to its distance in scipy's row major order, which is the
//     order np.bincount accumulates the float64 weights in. A distance's
//     mean divides by max(the diagonal length from the unit sizes, the
//     number of stored values). The scan stops after more than ten
//     consecutive distances without contacts (:251-258); the distance that
//     triggers the stop is kept.
//  2. --domains (from_bed_to_cut_interval, :274-339). The BED handle is an
//     argparse FileType and is read once, so with more than one matrix the
//     second matrix finds it exhausted and the reference stops on its
//     assertion "No region overlapped with bins."; that is reproduced.
//  3. The labels and the sums are dicts keyed by the matrix path, so a path
//     given twice has one entry, but the drawing loop iterates the command
//     line and draws it twice. --labels shorter than --matrices raises
//     KeyError once a curve of an unlabelled matrix is drawn.
//  4. --outFileData: one DataFrame.to_csv(sep='\t') per curve, with header and
//     index, appended to the same handle (:398-406).
//  5. --skipDiagonal is accepted and has no effect, as in the reference.
//  6. The FileType('w') arguments --plotFile and --outFileData are created at
//     parse time, as argparse opens them.
//
// Threading: none. The time goes into reading the matrix; the reduction is a
// single pass over the stored upper triangle.

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
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

#include "hicx/argparse.hpp"
#include "hicx/bins.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

namespace cli = hicx::cli;
namespace plot = hicx::plot;

const char* const kUsage =
    "usage: hicPlotDistVsCounts --matrices MATRICES [MATRICES ...] --plotFile file\n"
    "                           name [--labels LABELS [LABELS ...]]\n"
    "                           [--skipDiagonal] [--maxdepth INT bp] [--perchr]\n"
    "                           [--chromosomeExclude CHROMOSOMEEXCLUDE [CHROMOSOMEEXCLUDE ...]]\n"
    "                           [--domains DOMAINS] [--outFileData OUTFILEDATA]\n"
    "                           [--plotsize PLOTSIZE PLOTSIZE] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "This program creates distance vs. Hi-C counts plots. It can use several matrix\n"
    "files to compare them at once. If the `--perchr` option is given, each\n"
    "chromosome is plotted independently. When plotting multiple matrices, denser\n"
    "matrices are scaled down to match the sum of the smallest matrix.\n"
    "\n"
    "Required arguments:\n"
    "  --matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]\n"
    "                        Hi-C normalized (corrected) matrices. Each path should\n"
    "                        be separated by a space.\n"
    "  --plotFile file name, -o file name\n"
    "                        File name to save the file. The given file ending will\n"
    "                        be used to determine the image format. The available\n"
    "                        options are: .png, .emf, .eps, .pdf and .svg.\n"
    "\n"
    "Optional arguments:\n"
    "  --labels LABELS [LABELS ...]\n"
    "                        Label to assign to each matrix file. Each label should\n"
    "                        be separated by a space. Quote labels that contain\n"
    "                        spaces: E.g. --labels label1 \"labels 2\". If no labels\n"
    "                        are given then the file name is used.\n"
    "  --skipDiagonal, -s    If set, diagonal counts are not included.\n"
    "  --maxdepth INT bp     Maximum distance from diagonal to use. In other words,\n"
    "                        distances up to maxDepth are computed. Default is 3\n"
    "                        million bp.\n"
    "  --perchr              If given, computes and display distance versus Hi-C\n"
    "                        counts plots for each chromosome stored in the\n"
    "                        matrices passed to --matrices.\n"
    "  --chromosomeExclude CHROMOSOMEEXCLUDE [CHROMOSOMEEXCLUDE ...]\n"
    "                        Exclude the given list of chromosomes. This is useful\n"
    "                        for example to exclude the Y chromosome. The names of\n"
    "                        the chromosomes should be separated by space.\n"
    "  --domains DOMAINS     Bed file with domains coordinates: instead of\n"
    "                        evaluating the distance vs. Hi-C counts for intra\n"
    "                        chromosomal counts, compute it for intra-domains.\n"
    "  --outFileData OUTFILEDATA\n"
    "                        If given, the data underlying the plots is saved on\n"
    "                        this file.\n"
    "  --plotsize PLOTSIZE PLOTSIZE\n"
    "                        Width and height of the plot (in inches). Default is\n"
    "                        6*number of cols, 4 * number of rows. The maximum\n"
    "                        number of rows is 4. Example: --plotsize 6 5\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n"
    "\n"
    "C++ port: the distance means, the scale factors and --outFileData are computed\n"
    "in C++, and the figure is drawn by the hicexplorer_plot drawing layer with the\n"
    "matplotlib calls of the Python tool (HICX_PLOT_PYTHON names the interpreter).\n"
    "The C++-only option --plotData FILE writes the data of the figure as JSON to\n"
    "FILE instead of drawing it.\n";

struct Arguments {
    std::vector<std::string> matrices;
    std::string plot_file;
    std::optional<std::vector<std::string>> labels;
    std::int64_t maxdepth = 3000000;
    bool perchr = false;
    std::vector<std::string> chromosome_exclude;
    std::optional<std::string> domains;
    std::optional<std::string> out_file_data;
    std::optional<std::vector<double>> plotsize;
    std::optional<std::string> plot_data;
};

// argparse's FileType('w') opens the file while parsing; a failure is an
// argument error.
void open_like_filetype_w(const std::string& path, const std::string& option) {
    std::FILE* handle = std::fopen(path.c_str(), "w");
    if (handle == nullptr) {
        const int error = errno;
        std::fputs(kUsage, stderr);
        std::fprintf(stderr,
                     "hicPlotDistVsCounts: error: argument %s: can't open '%s': [Errno %d] %s: "
                     "'%s'\n",
                     option.c_str(), path.c_str(), error, std::strerror(error), path.c_str());
        std::exit(2);
    }
    std::fclose(handle);
}

Arguments parse_arguments(int argc, char** argv) {
    cli::Parser parser(
        "hicPlotDistVsCounts",
        "This program creates distance vs. Hi-C counts plots. It can use several matrix files "
        "to compare them at once. If the `--perchr` option is given, each chromosome is plotted "
        "independently. When plotting multiple matrices, denser matrices are scaled down to "
        "match the sum of the smallest matrix.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrices", "-m"})
        .nargs("+")
        .required()
        .input({"h5", "cool"})
        .help("Hi-C normalized (corrected) matrices. Each path should be separated by a space.");
    required.add({"--plotFile", "-o"})
        .file_type("w")
        .metavar("file name")
        .required()
        .output({"png", "emf", "eps", "pdf", "svg"})
        .help("File name to save the file. The given file ending will be used to determine the "
              "image format. The available options are: .png, .emf, .eps, .pdf and .svg.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--labels"})
        .nargs("+")
        .help("Label to assign to each matrix file. Each label should be separated by a space. "
              "Quote labels that contain spaces: E.g. --labels label1 \"labels 2\". If no labels "
              "are given then the file name is used.");
    optional.add({"--skipDiagonal", "-s"})
        .action(cli::Action::StoreTrue)
        .help("If set, diagonal counts are not included.");
    optional.add({"--maxdepth"})
        .metavar("INT bp")
        .type("int")
        .default_value(std::int64_t{3000000})
        .help("Maximum distance from diagonal to use. In other words, distances up to maxDepth "
              "are computed. Default is 3 million bp.");
    optional.add({"--perchr"})
        .action(cli::Action::StoreTrue)
        .help("If given, computes and display distance versus Hi-C counts plots for each "
              "chromosome stored in the matrices passed to --matrices.");
    optional.add({"--chromosomeExclude"})
        .nargs("+")
        .help("Exclude the given list of chromosomes. This is useful for example to exclude the "
              "Y chromosome. The names of the chromosomes should be separated by space.");
    optional.add({"--domains"})
        .file_type("r")
        .input({"bed"})
        .help("Bed file with domains coordinates: instead of evaluating the distance vs. Hi-C "
              "counts for intra chromosomal counts, compute it for intra-domains.");
    optional.add({"--outFileData"})
        .file_type("w")
        .output({"txt"})
        .help("If given, the data underlying the plots is saved on this file.");
    optional.add({"--plotsize"})
        .nargs(2)
        .type("float")
        .help("Width and height of the plot (in inches). Default is 6*number of cols, 4 * number "
              "of rows. The maximum number of rows is 4. Example: --plotsize 6 5");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    optional.add({"--plotData"})
        .metavar("FILE")
        .output({"json"})
        .cpp_only("The data of the figure as JSON, without drawing it (cpp/PLAN.md tier 7).")
        .help("Write the data the figure is drawn from as JSON to this file and do not draw "
              "the figure.");

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.matrices = ns.strs("matrices");
    args.plot_file = ns.str("plotFile");
    if (ns.given("labels")) {
        args.labels = ns.strs("labels");
    }
    args.maxdepth = ns.integer("maxdepth");
    args.perchr = ns.flag("perchr");
    if (ns.given("chromosomeExclude")) {
        args.chromosome_exclude = ns.strs("chromosomeExclude");
    }
    args.domains = ns.opt_str("domains");
    args.out_file_data = ns.opt_str("outFileData");
    if (ns.given("plotsize")) {
        args.plotsize = ns.reals("plotsize");
    }
    args.plot_data = ns.opt_str("plotData");

    open_like_filetype_w(args.plot_file, "--plotFile/-o");
    if (args.out_file_data.has_value()) {
        open_like_filetype_w(*args.out_file_data, "--outFileData");
    }
    return args;
}

// A Python exception the reference would end with (exit status 1).
class PythonError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

struct Entry {
    std::int32_t row;
    std::int32_t col;
    double value;
};

// One key of compute_distance_mean's result: the distances (k - 1) * binsize
// with their means (NaN where the diagonal was empty) and bin counts.
struct UnitMeans {
    std::string name;
    std::vector<std::int64_t> distance;
    std::vector<double> mean;
    std::vector<std::int64_t> nb;
};

bool starts_with(const std::string& text, const char* prefix) {
    return text.rfind(prefix, 0) == 0;
}

std::int64_t python_mod(std::int64_t a, std::int64_t m) {
    const std::int64_t r = a % m;
    return (r != 0 && ((r < 0) != (m < 0))) ? r + m : r;
}

// hicmatrix HiCMatrix.fit_cut_intervals: the bin starts the distances are
// computed from.
std::vector<std::int64_t> fitted_starts(const std::vector<hicx::CutInterval>& intervals) {
    std::vector<std::int64_t> starts;
    starts.reserve(intervals.size());
    for (const hicx::CutInterval& interval : intervals) {
        starts.push_back(interval.start);
    }
    if (intervals.size() <= 1) {
        return starts;
    }
    std::vector<std::string> order;
    std::unordered_map<std::string, std::vector<std::int64_t>> by_chrom;
    for (const hicx::CutInterval& interval : intervals) {
        auto [it, inserted] = by_chrom.try_emplace(interval.chrom);
        if (inserted) {
            order.push_back(interval.chrom);
        }
        it->second.push_back(interval.start);
    }
    std::vector<std::int64_t> diffs;
    for (const std::string& chrom : order) {
        const std::vector<std::int64_t>& values = by_chrom[chrom];
        for (std::size_t i = 1; i < values.size(); ++i) {
            diffs.push_back(values[i] - values[i - 1]);
        }
    }
    if (diffs.empty()) {
        throw PythonError("ValueError: need at least one array to concatenate "
                          "(HiCMatrix.fit_cut_intervals: every unit holds a single bin)");
    }
    std::sort(diffs.begin(), diffs.end());
    const std::size_t n = diffs.size();
    const double median_value =
        n % 2 == 1 ? static_cast<double>(diffs[n / 2])
                   : (static_cast<double>(diffs[n / 2 - 1]) + static_cast<double>(diffs[n / 2])) / 2.0;
    const auto median = static_cast<std::int64_t>(median_value);
    std::size_t deviating = 0;
    for (const hicx::CutInterval& interval : intervals) {
        if (interval.end - interval.start != median) {
            ++deviating;
        }
    }
    if (static_cast<double>(deviating) > static_cast<double>(intervals.size()) * 0.01) {
        if (median == 0) {
            throw PythonError("ZeroDivisionError: integer modulo by zero (fit_cut_intervals)");
        }
        for (std::int64_t& start : starts) {
            const std::int64_t down = -python_mod(start, median);
            const std::int64_t up = python_mod(-start, median);
            start += std::llabs(down) <= std::llabs(up) ? down : up;
        }
    }
    return starts;
}

// The sizes of intervalListToIntervalTree(intervals)[1] whose names do not
// start with _ignore_: a name seen again later keeps its first position and
// takes the later range.
std::vector<std::int64_t> unit_sizes(const std::vector<hicx::CutInterval>& intervals) {
    std::vector<std::string> order;
    std::unordered_map<std::string, std::pair<std::int64_t, std::int64_t>> bounds;
    if (intervals.empty()) {
        return {};
    }
    std::int64_t chr_start_id = 0;
    std::optional<std::string> previous;
    auto assign = [&](const std::string& name, std::int64_t first, std::int64_t last) {
        if (bounds.find(name) == bounds.end()) {
            order.push_back(name);
        }
        bounds[name] = {first, last};
    };
    std::int64_t id = 0;
    for (const hicx::CutInterval& interval : intervals) {
        if (!previous.has_value() || *previous != interval.chrom) {
            if (!previous.has_value()) {
                previous = interval.chrom;
            }
            assign(*previous, chr_start_id, id);
            chr_start_id = id;
            previous = interval.chrom;
        }
        ++id;
    }
    assign(intervals.back().chrom, chr_start_id, id);
    std::vector<std::int64_t> sizes;
    for (const std::string& name : order) {
        if (!starts_with(name, "_ignore_")) {
            sizes.push_back(bounds[name].second - bounds[name].first);
        }
    }
    return sizes;
}

// compute_distance_mean for one key: the stored upper triangle entries of the
// unit in row major order, with rows and columns relative to the unit.
UnitMeans unit_means(const std::string& name, const std::vector<Entry>& entries, std::size_t begin,
                     std::size_t end, std::int32_t offset, std::int64_t unit_bins,
                     const std::vector<hicx::CutInterval>& intervals, std::int64_t binsize,
                     std::int64_t maxdepth) {
    const std::vector<std::int64_t> starts = fitted_starts(intervals);
    std::unordered_map<std::string, std::int64_t> ids;
    std::vector<std::int64_t> chrom_id(intervals.size());
    std::vector<char> ignored(intervals.size());
    for (std::size_t i = 0; i < intervals.size(); ++i) {
        chrom_id[i] = ids.try_emplace(intervals[i].chrom, static_cast<std::int64_t>(ids.size()))
                          .first->second;
        ignored[i] = starts_with(intervals[i].chrom, "_ignore_") ? 1 : 0;
    }

    std::vector<double> sums;
    std::vector<std::int64_t> counts;
    for (std::size_t k = begin; k < end; ++k) {
        const auto row = static_cast<std::size_t>(entries[k].row - offset);
        const auto col = static_cast<std::size_t>(entries[k].col - offset);
        std::int64_t dist = 0;
        if (chrom_id[row] != chrom_id[col]) {
            dist = -1;  // inter unit; its chromosome name is '' and it is kept
        } else {
            if (ignored[row]) {
                continue;
            }
            dist = starts[col] - starts[row];
        }
        if (dist == -1) {
            dist = -binsize;
        }
        const std::int64_t bin =
            static_cast<std::int64_t>(static_cast<double>(dist) / static_cast<double>(binsize)) + 1;
        if (bin < 0) {
            throw PythonError("ValueError: 'list' argument must have no negative elements "
                              "(np.bincount over the bin distances)");
        }
        const auto index = static_cast<std::size_t>(bin);
        if (index >= sums.size()) {
            sums.resize(index + 1, 0.0);
            counts.resize(index + 1, 0);
        }
        sums[index] += entries[k].value;
        ++counts[index];
    }

    const std::vector<std::int64_t> sizes = unit_sizes(intervals);
    struct Mu {
        std::size_t k;
        double mean;
        std::int64_t nb;
    };
    std::vector<Mu> mu;
    std::vector<std::size_t> zero_value_bins;
    std::size_t consecutive = 0;
    for (std::size_t k = 0; k < sums.size(); ++k) {
        const double sum_value = sums[k];
        if (maxdepth != 0 && k == 0) {
            mu.push_back({k, std::nan(""), 0});
            continue;
        }
        double diagonal_length = 0.0;
        std::int64_t diagonal_int = 0;
        if (k == 0) {
            std::int64_t total_intra = unit_bins * unit_bins;
            for (const std::int64_t size : sizes) {
                total_intra -= size * size;
            }
            diagonal_length = static_cast<double>(total_intra) / 2.0;
        } else {
            const auto offset_k = static_cast<std::int64_t>(k - 1);
            for (const std::int64_t size : sizes) {
                if (size > offset_k) {
                    diagonal_int += size - offset_k;
                }
            }
            diagonal_length = static_cast<double>(diagonal_int);
        }
        if (static_cast<double>(counts[k]) > diagonal_length) {
            diagonal_length = static_cast<double>(counts[k]);
            diagonal_int = counts[k];
        }
        if (diagonal_length == 0.0) {
            mu.push_back({k, std::nan(""), 0});
            continue;
        }
        mu.push_back({k, sum_value / diagonal_length, diagonal_int});
        if (sum_value == 0.0) {
            if (!zero_value_bins.empty() && k - zero_value_bins.back() == 1) {
                ++consecutive;
            }
            zero_value_bins.push_back(k);
        }
        if (zero_value_bins.size() > 10 && consecutive > 10) {
            break;
        }
    }

    UnitMeans result;
    result.name = name;
    for (const Mu& m : mu) {
        if (m.k == 0) {
            continue;
        }
        const std::int64_t distance = static_cast<std::int64_t>(m.k - 1) * binsize;
        if (distance > maxdepth) {
            continue;
        }
        result.distance.push_back(distance);
        result.mean.push_back(m.mean);
        result.nb.push_back(m.nb);
    }
    return result;
}

std::string change_chrom_names(const std::string& chrom) {
    return starts_with(chrom, "chr") ? chrom.substr(3) : "chr" + chrom;
}

// from_bed_to_cut_interval over the kept bins.
std::vector<hicx::CutInterval> domain_intervals(const std::string& path,
                                                const std::vector<hicx::CutInterval>& original,
                                                const std::vector<std::string>& chrom_list) {
    std::vector<std::pair<std::string, hicx::BinRange>> ranges;
    {
        hicx::BinTable table(original);
        ranges = table.chrom_bin_boundaries();
    }
    auto range_of = [&](const std::string& chrom) -> std::optional<hicx::BinRange> {
        for (const auto& entry : ranges) {
            if (entry.first == chrom) {
                return entry.second;
            }
        }
        return std::nullopt;
    };
    auto known = [&](const std::string& chrom) {
        return std::find(chrom_list.begin(), chrom_list.end(), chrom) != chrom_list.end();
    };

    std::vector<std::optional<hicx::CutInterval>> result(original.size());
    std::ifstream in(path);
    std::string line;
    std::int64_t id = 0;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (!line.empty() && line[0] == '#') {
            continue;
        }
        std::istringstream fields_stream(line);
        std::vector<std::string> fields;
        for (std::string field; fields_stream >> field;) {
            fields.push_back(field);
        }
        if (fields.empty()) {
            throw PythonError("IndexError: list index out of range (an empty line in --domains)");
        }
        if (!known(fields[0])) {
            if (known(change_chrom_names(fields[0]))) {
                fields[0] = change_chrom_names(fields[0]);
            } else {
                continue;
            }
        }
        if (fields.size() < 3) {
            throw PythonError("ValueError: not enough values to unpack (a --domains line with "
                              "fewer than three fields)");
        }
        std::int64_t start = 0;
        std::int64_t end = 0;
        if (!cli::python_int(fields[1], &start)) {
            throw PythonError("ValueError: invalid literal for int() with base 10: '" + fields[1] +
                              "'");
        }
        if (!cli::python_int(fields[2], &end)) {
            throw PythonError("ValueError: invalid literal for int() with base 10: '" + fields[2] +
                              "'");
        }
        std::vector<std::size_t> overlapping;
        if (start < end) {
            const std::optional<hicx::BinRange> range = range_of(fields[0]);
            if (range.has_value()) {
                for (std::int64_t i = range->first; i < range->last; ++i) {
                    const hicx::CutInterval& bin = original[static_cast<std::size_t>(i)];
                    if (bin.chrom == fields[0] && bin.start < end && bin.end > start) {
                        overlapping.push_back(static_cast<std::size_t>(i));
                    }
                }
            }
        }
        if (overlapping.empty()) {
            throw PythonError("ValueError: min() arg is an empty sequence (a --domains region "
                              "overlaps no bin)");
        }
        const std::int64_t start_pos = original[overlapping.front()].start;
        for (const std::size_t i : overlapping) {
            if (result[i].has_value()) {
                throw PythonError("Exception: 2 features must not overlap the same bin.(" +
                                  original[i].chrom + ", " + std::to_string(original[i].start) +
                                  ", " + std::to_string(original[i].end) +
                                  ") is overlapped twice");
            }
            hicx::CutInterval interval = original[i];
            interval.chrom = std::to_string(id);
            interval.start = original[i].start - start_pos;
            interval.end = original[i].end - start_pos;
            result[i] = interval;
        }
        ++id;
    }
    if (id == 0) {
        throw PythonError("AssertionError: No region overlapped with bins.");
    }
    std::vector<hicx::CutInterval> intervals;
    intervals.reserve(original.size());
    for (std::size_t i = 0; i < original.size(); ++i) {
        if (result[i].has_value()) {
            intervals.push_back(*result[i]);
        } else {
            hicx::CutInterval interval = original[i];
            interval.chrom = "_ignore_" + std::to_string(i);
            interval.start = 0;
            interval.end = original[i].end - original[i].start;
            intervals.push_back(interval);
        }
    }
    return intervals;
}

struct MatrixMeans {
    hicx::Scalar sum;
    std::vector<UnitMeans> units;
};

MatrixMeans matrix_means(const std::string& path, const Arguments& args, bool domains_readable) {
    hicx::ToolMatrix matrix = hicx::ToolMatrix::load(path);
    MatrixMeans out;
    out.sum = matrix.matrix().sum();

    // keepOnlyTheseChr([x for x in interval_trees if x not in exclude])
    const std::vector<hicx::CutInterval>& all_intervals = matrix.cut_intervals();
    std::vector<std::int32_t> new_index(all_intervals.size(), -1);
    std::vector<char> selected(all_intervals.size(), 0);
    for (const auto& entry : matrix.boundaries()) {
        if (std::find(args.chromosome_exclude.begin(), args.chromosome_exclude.end(),
                      entry.first) != args.chromosome_exclude.end()) {
            continue;
        }
        for (std::int64_t i = entry.second.first; i < entry.second.last; ++i) {
            selected[static_cast<std::size_t>(i)] = 1;
        }
    }
    std::vector<hicx::CutInterval> kept;
    for (std::size_t i = 0; i < all_intervals.size(); ++i) {
        if (selected[i]) {
            new_index[i] = static_cast<std::int32_t>(kept.size());
            kept.push_back(all_intervals[i]);
        }
    }
    if (kept.empty()) {
        throw PythonError("ValueError: not enough values to unpack (no chromosome is left after "
                          "--chromosomeExclude, HiCMatrix.getBinSize)");
    }
    const hicx::BinTable kept_table(kept);
    std::vector<std::string> chrom_names;
    for (const auto& entry : kept_table.chrom_bin_boundaries()) {
        chrom_names.push_back(entry.first);
    }

    std::optional<std::vector<hicx::CutInterval>> custom;
    if (args.domains.has_value()) {
        if (!domains_readable) {
            throw PythonError("AssertionError: No region overlapped with bins. (the --domains "
                              "file handle was consumed by the first matrix)");
        }
        custom = domain_intervals(*args.domains, kept, chrom_names);
    }

    const std::int64_t binsize = kept_table.bin_size();
    std::int64_t max_depth_in_bins = 0;
    if (args.maxdepth != 0) {
        if (args.maxdepth < binsize) {
            std::fprintf(stderr, "Please specify a maxDepth larger than bin size (%lld)\n",
                         static_cast<long long>(binsize));
            std::exit(1);
        }
        max_depth_in_bins = static_cast<std::int64_t>(
            static_cast<double>(args.maxdepth) * 1.5 / static_cast<double>(binsize));
    }

    // triu(k=0) minus triu(k=max_depth_in_bins), zeros eliminated, row major.
    std::vector<Entry> entries;
    {
        const hicx::CsrMatrix& csr = matrix.matrix();
        const auto& indptr = csr.indptr();
        const auto& indices = csr.indices();
        const auto& data = csr.data();
        std::vector<Entry> row_entries;
        for (std::int64_t row = 0; row < csr.rows(); ++row) {
            const std::int32_t new_row = new_index[static_cast<std::size_t>(row)];
            if (new_row < 0) {
                continue;
            }
            row_entries.clear();
            bool sorted = true;
            for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
                 k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
                if (data[k] == 0.0 || indices[k] < row) {
                    continue;
                }
                const std::int32_t new_col = new_index[static_cast<std::size_t>(indices[k])];
                if (new_col < 0) {
                    continue;
                }
                if (args.maxdepth != 0 && new_col - new_row >= max_depth_in_bins) {
                    continue;
                }
                if (!row_entries.empty() && row_entries.back().col > new_col) {
                    sorted = false;
                }
                row_entries.push_back({new_row, new_col, data[k]});
            }
            if (!sorted) {
                std::sort(row_entries.begin(), row_entries.end(),
                          [](const Entry& a, const Entry& b) { return a.col < b.col; });
            }
            entries.insert(entries.end(), row_entries.begin(), row_entries.end());
        }
    }
    matrix = hicx::ToolMatrix();  // the per distance reduction needs only the entries

    const std::vector<hicx::CutInterval>& unit_intervals = custom.has_value() ? *custom : kept;
    if (args.perchr) {
        for (const auto& entry : kept_table.chrom_bin_boundaries()) {
            const auto first = static_cast<std::int32_t>(entry.second.first);
            const auto last = static_cast<std::int32_t>(entry.second.last);
            const auto lower = std::lower_bound(
                entries.begin(), entries.end(), first,
                [](const Entry& e, std::int32_t value) { return e.row < value; });
            const auto upper = std::lower_bound(
                entries.begin(), entries.end(), last,
                [](const Entry& e, std::int32_t value) { return e.row < value; });
            std::vector<Entry> block;
            for (auto it = lower; it != upper; ++it) {
                if (it->col < last) {
                    block.push_back(*it);
                }
            }
            const std::vector<hicx::CutInterval> slice(unit_intervals.begin() + first,
                                                       unit_intervals.begin() + last);
            out.units.push_back(unit_means(entry.first, block, 0, block.size(), first,
                                           last - first, slice, binsize, args.maxdepth));
        }
    } else {
        out.units.push_back(unit_means("all", entries, 0, entries.size(), 0,
                                       static_cast<std::int64_t>(unit_intervals.size()),
                                       unit_intervals, binsize, args.maxdepth));
    }
    return out;
}

std::string basename_of(const std::string& path) {
    const std::string::size_type slash = path.rfind('/');
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

// A to_csv field: quoted when it holds the separator, a quote or a line break.
std::string csv_field(const std::string& text) {
    if (text.find_first_of("\t\"\n\r") == std::string::npos) {
        return text;
    }
    std::string out = "\"";
    for (const char c : text) {
        out += c;
        if (c == '"') {
            out += '"';
        }
    }
    return out + "\"";
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    if (const int refused = hicx::plot::preflight("hicPlotDistVsCounts", !args.plot_data.has_value()); refused != 0) {
        return refused;
    }

    // labels = OrderedDict(...) keyed by matrix path.
    std::map<std::string, std::string> labels;
    if (!args.labels.has_value()) {
        for (const std::string& path : args.matrices) {
            labels[path] = basename_of(path);
        }
    } else {
        for (std::size_t i = 0; i < args.matrices.size() && i < args.labels->size(); ++i) {
            labels[args.matrices[i]] = (*args.labels)[i];
        }
    }

    std::map<std::string, MatrixMeans> means;
    std::set<std::string> chroms;
    try {
        bool domains_readable = true;
        for (const std::string& path : args.matrices) {
            MatrixMeans result = matrix_means(path, args, domains_readable);
            domains_readable = false;
            for (const UnitMeans& unit : result.units) {
                if (unit.distance.size() > 1) {
                    chroms.insert(unit.name);
                }
            }
            means[path] = std::move(result);
        }
    } catch (const PythonError& error) {
        std::fprintf(stderr, "hicPlotDistVsCounts: %s\n", error.what());
        return 1;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicPlotDistVsCounts: %s\n", error.what());
        return 1;
    }

    // scale_factor = float(min_sum) / mat_sum; numpy scalars divide in float64.
    double min_sum = 0.0;
    bool first = true;
    for (const auto& entry : means) {
        const double value = entry.second.sum.as_double();
        if (first || value < min_sum) {
            min_sum = value;
            first = false;
        }
    }

    std::string table;
    std::vector<std::string> matrices_json;
    for (const std::string& path : args.matrices) {
        const MatrixMeans& result = means[path];
        const double scale = min_sum / result.sum.as_double();
        std::vector<std::string> series_json;
        const auto label = labels.find(path);
        for (const UnitMeans& unit : result.units) {
            if (unit.distance.size() <= 1) {
                continue;
            }
            std::vector<std::int64_t> x;
            std::vector<double> y;
            std::vector<std::int64_t> z;
            for (std::size_t i = 0; i < unit.distance.size(); ++i) {
                if (unit.mean[i] > 0.0) {
                    x.push_back(unit.distance[i]);
                    y.push_back(unit.mean[i] * scale);
                    z.push_back(unit.nb[i]);
                }
            }
            if (x.size() <= 1) {
                continue;
            }
            const bool needs_label = !(args.perchr && args.matrices.size() == 1) ||
                                     args.out_file_data.has_value();
            if (needs_label && label == labels.end()) {
                std::fprintf(stderr, "hicPlotDistVsCounts: KeyError: '%s' (--labels has fewer "
                                     "entries than --matrices)\n",
                             path.c_str());
                return 1;
            }
            if (args.out_file_data.has_value()) {
                table += "\tMatrix\tChromosome\tDistance\tContacts\tNumber_bins\tScale_factor\n";
                const std::string prefix = csv_field(label->second) + "\t" + csv_field(unit.name);
                const std::string scale_text = hicx::npy::float_repr(scale);
                for (std::size_t i = 0; i < x.size(); ++i) {
                    table += std::to_string(i) + "\t" + prefix + "\t" + std::to_string(x[i]) +
                             "\t" + hicx::npy::float_repr(y[i]) + "\t" + std::to_string(z[i]) +
                             "\t" + scale_text + "\n";
                }
            }
            plot::JsonObject series;
            series.add("chrom", plot::json_string(unit.name));
            series.add("x", plot::json_ints(x));
            series.add("y", plot::json_numbers(y));
            series_json.push_back(series.str());
        }
        plot::JsonObject matrix;
        matrix.add("label", label == labels.end() ? "null" : plot::json_string(label->second));
        matrix.add("series", plot::json_list(series_json));
        matrices_json.push_back(matrix.str());
    }

    if (args.out_file_data.has_value()) {
        std::ofstream out(*args.out_file_data, std::ios::binary | std::ios::trunc);
        out << table;
        out.close();
        if (!out) {
            std::fprintf(stderr, "hicPlotDistVsCounts: cannot write %s\n",
                         args.out_file_data->c_str());
            return 1;
        }
    }

    plot::JsonObject data;
    data.add("plotFile", plot::json_string(args.plot_file));
    data.add("perchr", plot::json_bool(args.perchr));
    data.add("num_chroms", plot::json_int(static_cast<std::int64_t>(chroms.size())));
    data.add("maxdepth", plot::json_int(args.maxdepth));
    data.add("plotsize", args.plotsize.has_value() ? plot::json_numbers(*args.plotsize) : "null");
    data.add("matrices", plot::json_list(matrices_json));

    hicx::report_resource_usage("hicPlotDistVsCounts");
    return plot::draw("hicPlotDistVsCounts", data.str(), args.plot_data);
}
