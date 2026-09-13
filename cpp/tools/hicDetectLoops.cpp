// Port of hicexplorer/hicDetectLoops.py.
//
// The pipeline, per chromosome:
//
//   1. load that chromosome's matrix, cut to --maxLoopDistance    (this file)
//   2. triu, drop the main diagonal                               (this file)
//   3. observed over expected                       (hicx::loops::obs_exp_*)
//   4. fit a negative binomial per genomic distance and keep the pixels whose
//      tail probability clears --pValuePreselection
//                                          (hicx::loops::preselect_candidates)
//   5. intersect with --peakInteractionsThreshold on the raw counts
//   6. keep a candidate only if it is the maximum of its own window
//                                           (hicx::loops::neighborhood_merge)
//   7. the donut test: three Wilcoxon rank sums against the horizontal, the
//      vertical and the bottom left corner, then one against the whole
//      background                        (hicx::loops::candidate_region_test)
//   8. map to genomic coordinates, drop pairs beyond --maxLoopDistance, write
//
// Everything numeric lives in detect_loops_impl.hpp so that the unit tests can
// reach it. This file is the argument parser, the loader and the writer.
//
// Two deliberate differences from the Python, both forced by
// cpp/OPTIMIZATION.md 3 and both recorded in the report rather than hidden:
//
//   * **Line order.** The Python starts one process per chromosome and appends
//     each result to the output as that process happens to finish, so its line
//     order is not reproducible: measured on small_test_matrix.h5 with
//     --threads 4, three runs gave two different orderings of the same 52
//     loops, and --threads 8 --threadsPerChromosome 4 gave a third. The set is
//     always identical. This port writes the chromosomes in a fixed order
//     instead: the order of --chromosomes when it is given, otherwise the
//     order the chromosomes appear in the file. Where the Python's order is
//     itself well defined, which is --threads 1 or a single chromosome, the
//     two agree line for line.
//   * **The chromosome shuffle.** hicDetectLoops.py:982-1001 reorders the
//     chromosomes of a cool file as "largest, then --threads smallest, then
//     second largest, ..." to cap the peak memory of the process pool. The
//     resulting order depends on --threads, so reproducing it would make the
//     output depend on the thread count. It is a scheduling heuristic over a
//     permutation of the same set, and this port does not need it because it
//     holds one chromosome at a time rather than --threads of them.
//
// Threading. Chromosomes are processed one after another and the work inside a
// chromosome is threaded over individual distances and individual candidates,
// on --threads * --threadsPerChromosome workers, which is the same total the
// Python's help text promises. Doing it this way rather than the other way
// round is what keeps the resident set at one chromosome instead of --threads
// of them, and it makes the result independent of both thread counts because
// every parallel item is independent and is combined in index order.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "detect_loops_impl.hpp"
#include "hicx/adjust_ops.hpp"
#include "hicx/bins.hpp"
#include "hicx/cool_file.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/matrix_ops.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicDetectLoops --matrix MATRIX --outFileName OUTFILENAME\n"
    "                      [--peakWidth PEAKWIDTH] [--windowSize WINDOWSIZE]\n"
    "                      [--pValuePreselection PVALUEPRESELECTION]\n"
    "                      [--peakInteractionsThreshold PEAKINTERACTIONSTHRESHOLD]\n"
    "                      [--obsExpThreshold OBSEXPTHRESHOLD] [--pValue PVALUE]\n"
    "                      [--maxLoopDistance MAXLOOPDISTANCE]\n"
    "                      [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]\n"
    "                      [--threads THREADS]\n"
    "                      [--threadsPerChromosome THREADSPERCHROMOSOME]\n"
    "                      [--expected {mean,mean_nonzero,mean_nonzero_ligation}]\n"
    "                      [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Computes enriched regions (peaks) or long range contacts on the given "
    "contact matrix.\n"
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        The matrix to compute the loop detection on.\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        Outfile name to store the detected loops. The file\n"
    "                        will in bedgraph format.\n"
    "\n"
    "Optional arguments:\n"
    "  --peakWidth PEAKWIDTH, -pw PEAKWIDTH\n"
    "                        The width of the peak region in bins. (Default: 2).\n"
    "  --windowSize WINDOWSIZE, -w WINDOWSIZE\n"
    "                        The window size for the neighborhood region the peak\n"
    "                        is located in. (Default: 5).\n"
    "  --pValuePreselection PVALUEPRESELECTION, -pp PVALUEPRESELECTION\n"
    "                        Only candidates with p-values less the given threshold\n"
    "                        will be considered as candidates. Can a single value or\n"
    "                        a threshold file created by hicCreateThresholdFile\n"
    "                        (Default: 0.1).\n"
    "  --peakInteractionsThreshold PEAKINTERACTIONSTHRESHOLD, -pit PEAKINTERACTIONSTHRESHOLD\n"
    "                        The minimum number of interactions a detected peaks\n"
    "                        needs to have to be considered (Default: 10).\n"
    "  --obsExpThreshold OBSEXPTHRESHOLD, -oet OBSEXPTHRESHOLD\n"
    "                        The minimum number of obs/exp interactions a detected\n"
    "                        peaks needs to have to be considered (Default: 1.5).\n"
    "  --pValue PVALUE, -p PVALUE\n"
    "                        Rejection level for Anderson-Darling or Wilcoxon-rank\n"
    "                        sum test for H0. (Default: 0.025).\n"
    "  --maxLoopDistance MAXLOOPDISTANCE\n"
    "                        Maximum genomic distance of a loop (Default: 2000000).\n"
    "  --chromosomes CHROMOSOMES [CHROMOSOMES ...]\n"
    "                        Chromosomes to include in the analysis. If not set, all\n"
    "                        chromosomes are included.\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads to use (Default: 4).\n"
    "  --threadsPerChromosome THREADSPERCHROMOSOME, -tpc THREADSPERCHROMOSOME\n"
    "                        Number of threads to use per parallel thread processing\n"
    "                        a chromosome (Default: 4).\n"
    "  --expected {mean,mean_nonzero,mean_nonzero_ligation}, -exp {mean,mean_nonzero,mean_nonzero_ligation}\n"
    "                        Method to compute the expected value per distance\n"
    "                        (Default: mean).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::string matrix;
    std::string out_file_name;
    std::int64_t peak_width = 2;
    std::int64_t window_size = 5;
    std::string p_value_preselection = "0.1";
    double peak_interactions_threshold = 10.0;
    double obs_exp_threshold = 1.5;
    double p_value = 0.025;
    std::int64_t max_loop_distance = 2000000;
    std::optional<std::vector<std::string>> chromosomes;
    int threads = 4;
    int threads_per_chromosome = 4;
    std::string expected = "mean";
};

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicDetectLoops: error: %s\n", message.c_str());
    std::exit(2);
}

std::int64_t parse_int(const std::string& name, const std::string& text) {
    try {
        std::size_t consumed = 0;
        const long long value = std::stoll(text, &consumed);
        if (consumed != text.size()) {
            throw std::invalid_argument("trailing");
        }
        return value;
    } catch (const std::exception&) {
        fail("argument " + name + ": invalid int value: '" + text + "'");
    }
}

double parse_double(const std::string& name, const std::string& text) {
    try {
        std::size_t consumed = 0;
        const double value = std::stod(text, &consumed);
        if (consumed != text.size()) {
            throw std::invalid_argument("trailing");
        }
        return value;
    } catch (const std::exception&) {
        fail("argument " + name + ": invalid float value: '" + text + "'");
    }
}

Arguments parse_arguments(int argc, char** argv) {
    Arguments args;
    bool matrix_seen = false;
    bool out_seen = false;

    std::vector<std::string> tokens;
    tokens.reserve(static_cast<std::size_t>(argc));
    for (int i = 1; i < argc; ++i) {
        tokens.emplace_back(argv[i]);
    }

    for (std::size_t i = 0; i < tokens.size(); ++i) {
        std::string name = tokens[i];
        std::optional<std::string> inline_value;
        const std::size_t equals = name.find('=');
        if (equals != std::string::npos && name.rfind("--", 0) == 0) {
            inline_value = name.substr(equals + 1);
            name = name.substr(0, equals);
        }
        const auto next_value = [&](const char* option) -> std::string {
            if (inline_value.has_value()) {
                return *inline_value;
            }
            if (i + 1 >= tokens.size()) {
                fail(std::string("argument ") + option + ": expected one argument");
            }
            return tokens[++i];
        };

        if (name == "-h" || name == "--help") {
            std::fputs(kUsage, stdout);
            std::fputs(kHelp, stdout);
            std::exit(0);
        }
        if (name == "--version") {
            std::printf("hicDetectLoops %s\n", hicx::kVersion);
            std::exit(0);
        }
        if (name == "-m" || name == "--matrix") {
            args.matrix = next_value("--matrix/-m");
            matrix_seen = true;
        } else if (name == "-o" || name == "--outFileName") {
            args.out_file_name = next_value("--outFileName/-o");
            out_seen = true;
        } else if (name == "-pw" || name == "--peakWidth") {
            args.peak_width = parse_int("--peakWidth/-pw", next_value("--peakWidth/-pw"));
        } else if (name == "-w" || name == "--windowSize") {
            args.window_size =
                parse_int("--windowSize/-w", next_value("--windowSize/-w"));
        } else if (name == "-pp" || name == "--pValuePreselection") {
            args.p_value_preselection = next_value("--pValuePreselection/-pp");
        } else if (name == "-pit" || name == "--peakInteractionsThreshold") {
            args.peak_interactions_threshold =
                parse_double("--peakInteractionsThreshold/-pit",
                             next_value("--peakInteractionsThreshold/-pit"));
        } else if (name == "-oet" || name == "--obsExpThreshold") {
            args.obs_exp_threshold = parse_double("--obsExpThreshold/-oet",
                                                  next_value("--obsExpThreshold/-oet"));
        } else if (name == "-p" || name == "--pValue") {
            args.p_value = parse_double("--pValue/-p", next_value("--pValue/-p"));
        } else if (name == "--maxLoopDistance") {
            args.max_loop_distance =
                parse_int("--maxLoopDistance", next_value("--maxLoopDistance"));
        } else if (name == "--chromosomes") {
            std::vector<std::string> names;
            if (inline_value.has_value()) {
                names.push_back(*inline_value);
            }
            while (i + 1 < tokens.size() && tokens[i + 1].rfind("-", 0) != 0) {
                names.push_back(tokens[++i]);
            }
            if (names.empty()) {
                fail("argument --chromosomes: expected at least one argument");
            }
            args.chromosomes = std::move(names);
        } else if (name == "-t" || name == "--threads") {
            args.threads =
                static_cast<int>(parse_int("--threads/-t", next_value("--threads/-t")));
        } else if (name == "-tpc" || name == "--threadsPerChromosome") {
            args.threads_per_chromosome = static_cast<int>(parse_int(
                "--threadsPerChromosome/-tpc", next_value("--threadsPerChromosome/-tpc")));
        } else if (name == "-exp" || name == "--expected") {
            args.expected = next_value("--expected/-exp");
            if (args.expected != "mean" && args.expected != "mean_nonzero" &&
                args.expected != "mean_nonzero_ligation") {
                fail("argument --expected/-exp: invalid choice: '" + args.expected +
                     "' (choose from 'mean', 'mean_nonzero', 'mean_nonzero_ligation')");
            }
        } else {
            fail("unrecognized arguments: " + tokens[i]);
        }
    }

    std::string missing;
    const auto require = [&missing](bool seen, const char* option) {
        if (!seen) {
            missing += missing.empty() ? option : std::string(", ") + option;
        }
    };
    require(matrix_seen, "--matrix/-m");
    require(out_seen, "--outFileName/-o");
    if (!missing.empty()) {
        fail("the following arguments are required: " + missing);
    }
    return args;
}

// hicDetectLoops.read_threshold_file (:940-952).
std::map<std::int64_t, double> read_threshold_file(const std::string& path) {
    std::map<std::int64_t, double> thresholds;
    std::ifstream input(path);
    if (!input) {
        std::fprintf(stderr, "hicDetectLoops: error: cannot open %s\n", path.c_str());
        std::exit(1);
    }
    std::string line;
    while (std::getline(input, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n' ||
                                 line.back() == ' ' || line.back() == '\t')) {
            line.pop_back();
        }
        if (!line.empty() && line[0] == '#') {
            continue;
        }
        if (line.empty()) {
            break;  // the Python stops at the first blank line
        }
        const std::size_t tab = line.find('\t');
        if (tab == std::string::npos) {
            continue;
        }
        thresholds[std::stoll(line.substr(0, tab))] = std::stod(line.substr(tab + 1));
    }
    return thresholds;
}

// The state one chromosome contributes.
struct ChromosomeMatrix {
    hicx::CsrMatrix matrix;  // triu, diagonal removed, Symmetry::Full
    std::vector<hicx::CutInterval> cut_intervals;
    std::int64_t bin_size = 0;
};

// triu(matrix, k=1) of the represented matrix, keeping only entries whose
// distance is within `limit` bins. Both operations are symmetric, so applying
// them to the stored upper triangle is the same as applying them to the full
// matrix and taking triu afterwards, which is what the Python does.
hicx::CsrMatrix upper_band(const hicx::CsrMatrix& matrix, double distance_limit,
                           bool inclusive) {
    hicx::CsrMatrix::Arrays arrays;
    arrays.rows = matrix.rows();
    arrays.cols = matrix.cols();
    arrays.dtype = matrix.dtype();
    arrays.symmetry = hicx::Symmetry::Full;
    arrays.indptr.assign(static_cast<std::size_t>(matrix.rows()) + 1, 0);

    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& values = matrix.data();
    // Whether the matrix is stored as its upper triangle or with both
    // triangles materialised makes no difference here, because only entries
    // with column > row are kept either way.
    for (std::int64_t row = 0; row < matrix.rows(); ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t column = static_cast<std::int64_t>(indices[k]);
            if (column <= row) {
                continue;  // triu(k=0) then the removal of the diagonal
            }
            const double distance = static_cast<double>(column - row);
            const bool inside = inclusive ? distance <= distance_limit
                                          : distance < distance_limit;
            if (!inside) {
                continue;
            }
            if (values[k] == 0.0) {
                continue;  // eliminate_zeros
            }
            arrays.indices.push_back(indices[k]);
            arrays.data.push_back(values[k]);
            ++arrays.indptr[static_cast<std::size_t>(row) + 1];
        }
    }
    for (std::size_t i = 1; i < arrays.indptr.size(); ++i) {
        arrays.indptr[i] += arrays.indptr[i - 1];
    }
    return hicx::CsrMatrix::adopt(std::move(arrays));
}

// hicmatrix.lib.Cool.load for a single chromosome with pDistance set
// (cool.py:120-146). The pixels of the chromosome block are read, the ones
// further from the diagonal than `distance // binsize` bins are dropped, and
// the counts are stored as float32, which is the dtype of the lil_matrix the
// Python stages them in and therefore the dtype obs_exp_matrix casts back to.
//
// The pixel table and the bin table are read once for the whole file and
// passed in, not reopened per chromosome. Reading them inside the loop was
// measured first and is what the timing note in the report compares against:
// on the 25 chromosome GSE63525 cool it cost 0.89 s of CPU against 0.24 s,
// because the whole 705k pixel table was decompressed 25 times over.
ChromosomeMatrix load_cool_chromosome(const hicx::CoolFile& cool,
                                      const hicx::CsrMatrix& whole,
                                      const std::vector<hicx::CutInterval>& bins,
                                      const std::string& chromosome,
                                      std::int64_t max_loop_distance) {
    std::int64_t first = -1;
    std::int64_t last = 0;
    for (std::size_t bin = 0; bin < bins.size(); ++bin) {
        if (bins[bin].chrom != chromosome) {
            continue;
        }
        if (first < 0) {
            first = static_cast<std::int64_t>(bin);
        }
        last = static_cast<std::int64_t>(bin) + 1;
    }
    if (first < 0) {
        throw std::runtime_error("Chromosome name not in matrix. '" + chromosome + "'");
    }

    // cooler.Cooler.binsize: the file's bin-size attribute, absent for a
    // variable bin cooler, in which case the Python fetches the whole
    // chromosome and never applies the distance filter.
    std::optional<std::int64_t> file_bin_size;
    const hicx::json::Value* attribute = cool.info_value("bin-size");
    if (attribute != nullptr && attribute->is_number()) {
        const double value = attribute->as_double();
        if (value > 0) {
            file_bin_size = static_cast<std::int64_t>(value);
        }
    }

    hicx::CsrMatrix::Arrays arrays;
    arrays.rows = last - first;
    arrays.cols = last - first;
    arrays.symmetry = hicx::Symmetry::Full;
    arrays.indptr.assign(static_cast<std::size_t>(arrays.rows) + 1, 0);
    const bool restrict_distance = file_bin_size.has_value();
    const std::int64_t band = restrict_distance ? max_loop_distance / *file_bin_size : 0;
    arrays.dtype = restrict_distance ? "float32" : whole.dtype();

    const std::vector<std::int64_t>& indptr = whole.indptr();
    const std::vector<std::int32_t>& indices = whole.indices();
    const std::vector<double>& values = whole.data();
    for (std::int64_t row = first; row < last; ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t column = static_cast<std::int64_t>(indices[k]);
            if (column < first || column >= last) {
                continue;
            }
            if (restrict_distance && (column - row) >= band) {
                continue;
            }
            double value = values[k];
            if (restrict_distance) {
                value = static_cast<double>(static_cast<float>(value));
            }
            if (value == 0.0) {
                continue;  // lil_matrix does not store a zero
            }
            arrays.indices.push_back(static_cast<std::int32_t>(column - first));
            arrays.data.push_back(value);
            ++arrays.indptr[static_cast<std::size_t>(row - first) + 1];
        }
    }
    for (std::size_t i = 1; i < arrays.indptr.size(); ++i) {
        arrays.indptr[i] += arrays.indptr[i - 1];
    }

    ChromosomeMatrix result;
    result.matrix = hicx::CsrMatrix::adopt(std::move(arrays));
    result.cut_intervals.assign(bins.begin() + static_cast<std::ptrdiff_t>(first),
                                bins.begin() + static_cast<std::ptrdiff_t>(last));
    result.bin_size = hicx::BinTable(result.cut_intervals).bin_size();
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);

    if (args.window_size <= args.peak_width) {
        std::fprintf(stderr,
                     "ERROR:hicexplorer.hicDetectLoops:The window size (%lld) must be "
                     "larger than the peakWidth (%lld)\n",
                     static_cast<long long>(args.window_size),
                     static_cast<long long>(args.peak_width));
        return 1;
    }

    const unsigned int threads_per_chromosome = static_cast<unsigned int>(
        std::max(1, args.threads_per_chromosome));
    const unsigned int workers =
        std::max(1U, static_cast<unsigned int>(std::max(1, args.threads)) *
                         threads_per_chromosome);

    // --pValuePreselection is either a float or the name of a threshold file,
    // decided by whether float() succeeds (hicDetectLoops.py:895-898).
    double preselection_value = 0.1;
    std::map<std::int64_t, double> preselection_table;
    try {
        std::size_t consumed = 0;
        preselection_value = std::stod(args.p_value_preselection, &consumed);
        if (consumed != args.p_value_preselection.size()) {
            throw std::invalid_argument("trailing");
        }
    } catch (const std::exception&) {
        preselection_table = read_threshold_file(args.p_value_preselection);
    }

    const bool is_cooler = hicx::check_cooler(args.matrix);

    // The chromosome list, and the whole matrix when the input is an h5, which
    // has no per chromosome loader.
    std::vector<std::string> chromosomes;
    hicx::MatrixData whole;
    std::optional<hicx::CoolFile> cool;
    hicx::CsrMatrix cool_pixels;
    std::vector<hicx::CutInterval> cool_bins;
    if (!is_cooler) {
        whole = hicx::read_hicexplorer_h5(args.matrix);
        // hiCMatrix.__init__ with pUpperTriangleOnly unset calls
        // fillLowerTriangle, which the port expresses as a relabelling
        // whenever nothing is stored below the diagonal.
        whole.matrix.symmetrize_in_place();
        if (args.chromosomes.has_value()) {
            chromosomes = *args.chromosomes;
        } else {
            for (const auto& [name, range] :
                 hicx::chrom_bin_boundaries(whole.cut_intervals)) {
                (void)range;
                chromosomes.push_back(name);
            }
        }
    } else {
        cool.emplace(args.matrix);
        cool_bins = cool->read_bins();
        cool_pixels = cool->read_matrix();
        if (args.chromosomes.has_value()) {
            chromosomes = *args.chromosomes;
        } else {
            // cooler.Cooler(path).chromsizes, in chroms table order. The
            // Python then permutes it as "largest, then --threads smallest,
            // ..." to cap the peak memory of its process pool; that
            // permutation depends on --threads and is not reproduced, see the
            // header comment.
            chromosomes = cool->chrom_names();
        }
    }

    struct Loop {
        std::string chrom_x;
        std::int64_t start_x = 0;
        std::int64_t end_x = 0;
        std::string chrom_y;
        std::int64_t start_y = 0;
        std::int64_t end_y = 0;
        double pvalue = 0.0;
    };
    std::vector<Loop> mapped_loops;
    const bool single_chromosome = chromosomes.size() == 1;

    for (std::size_t chromosome_index = 0; chromosome_index < chromosomes.size();
         ++chromosome_index) {
        const std::string& chromosome = chromosomes[chromosome_index];
        ChromosomeMatrix block;
        try {
            if (is_cooler) {
                block = load_cool_chromosome(*cool, cool_pixels, cool_bins,
                                             chromosome, args.max_loop_distance);
                if (chromosome_index + 1 == chromosomes.size()) {
                    // The whole pixel table is not needed once the last block
                    // has been cut out of it, and on a single chromosome cool
                    // it is the largest thing in the process by an order of
                    // magnitude: 742 MB against the 51 MB band that is
                    // actually processed on gm12878_chr1.cool. Releasing it
                    // here rather than at the end of main takes the peak from
                    // 962 MB to 877 MB, measured, three runs each.
                    cool_pixels = hicx::CsrMatrix();
                }
            } else {
                // keepOnlyTheseChr: a monotone selection of that chromosome's
                // bins out of the whole matrix.
                std::vector<std::int64_t> selection;
                for (std::size_t bin = 0; bin < whole.cut_intervals.size(); ++bin) {
                    if (whole.cut_intervals[bin].chrom == chromosome) {
                        selection.push_back(static_cast<std::int64_t>(bin));
                    }
                }
                if (selection.empty()) {
                    throw std::runtime_error("Chromosome name not in matrix. '" +
                                             chromosome + "'");
                }
                block.matrix = hicx::select_bins(whole.matrix, selection);
                block.cut_intervals.reserve(selection.size());
                for (const std::int64_t bin : selection) {
                    block.cut_intervals.push_back(
                        whole.cut_intervals[static_cast<std::size_t>(bin)]);
                }
                block.bin_size = hicx::BinTable(block.cut_intervals).bin_size();
            }
        } catch (const std::exception& error) {
            std::fprintf(stderr, "ERROR:hicexplorer.hicDetectLoops:%s\n", error.what());
            return 1;
        }

        if (block.bin_size <= 0) {
            continue;
        }

        // The two paths cut the distance differently and the difference is
        // real, not a rounding choice. The cool loader keeps bin2 - bin1
        // strictly below `maxLoopDistance // binsize`, an integer floor
        // (cool.py:128,137). The h5 path masks `distances > maxLoopDistance /
        // binSize`, a float division with a non strict comparison
        // (hicDetectLoops.py:820-823). With maxLoopDistance 3000000 on a
        // 2500000 bp matrix the first keeps distance 0 only and the second
        // keeps distances 0 and 1.
        hicx::CsrMatrix matrix;
        if (is_cooler) {
            // The band was already applied while loading; only triu and the
            // diagonal removal are left.
            matrix = upper_band(block.matrix, std::numeric_limits<double>::infinity(),
                                true);
        } else {
            const double limit = static_cast<double>(args.max_loop_distance) /
                                 static_cast<double>(block.bin_size);
            matrix = upper_band(block.matrix, limit, true);
        }
        block.matrix = hicx::CsrMatrix();

        if (matrix.stored_nnz() == 0 || matrix.rows() < 5 || matrix.cols() < 5) {
            continue;
        }

        const hicx::loops::DistanceGroups groups =
            hicx::loops::group_by_distance(matrix);

        hicx::CsrMatrix obs_exp;
        if (args.expected == "mean") {
            obs_exp = hicx::loops::obs_exp_matrix(matrix, groups, workers);
        } else {
            obs_exp = hicx::loops::obs_exp_matrix_non_zero(
                matrix, args.expected == "mean_nonzero_ligation");
        }
        obs_exp.eliminate_zeros();
        if (obs_exp.stored_nnz() != matrix.stored_nnz()) {
            // hicDetectLoops.py:888: the raw and the obs/exp matrix must have
            // the same number of stored values, otherwise the mask built from
            // one cannot index the other. An integer obs/exp whose ratios
            // truncate to zero lands here.
            continue;
        }

        const hicx::loops::DistanceGroups obs_exp_groups =
            hicx::loops::group_by_distance(obs_exp);

        hicx::loops::PreselectionResult preselection;
        try {
            preselection = hicx::loops::preselect_candidates(
                obs_exp, obs_exp_groups, preselection_value, preselection_table,
                block.bin_size, args.obs_exp_threshold, workers);
        } catch (const std::exception& error) {
            std::fprintf(stderr, "ERROR:hicexplorer.hicDetectLoops:%s\n", error.what());
            return 1;
        }

        // mask &= matrix.data >= peakInteractionsThreshold, position by
        // position in CSR order; the two matrices share their sparsity
        // pattern, which is what the length check above established.
        std::vector<hicx::loops::Candidate> candidates;
        {
            const std::vector<std::int64_t>& indptr = obs_exp.indptr();
            const std::vector<std::int32_t>& indices = obs_exp.indices();
            const std::vector<double>& raw = matrix.data();
            for (std::int64_t row = 0; row < obs_exp.rows(); ++row) {
                const std::size_t begin =
                    static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
                const std::size_t end =
                    static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
                for (std::size_t k = begin; k < end; ++k) {
                    if (preselection.mask[k] == 0) {
                        continue;
                    }
                    if (!(raw[k] >= args.peak_interactions_threshold)) {
                        continue;
                    }
                    candidates.push_back(hicx::loops::Candidate{
                        row, static_cast<std::int64_t>(indices[k])});
                }
            }
        }
        if (candidates.empty()) {
            continue;
        }

        candidates = hicx::loops::neighborhood_merge(candidates, args.window_size,
                                                     obs_exp, workers);
        if (candidates.empty()) {
            continue;
        }
        const hicx::loops::RegionTestResult tested = hicx::loops::candidate_region_test(
            obs_exp, candidates, args.window_size, args.p_value, args.peak_width,
            workers);
        if (tested.candidates.empty()) {
            continue;
        }

        for (std::size_t i = 0; i < tested.candidates.size(); ++i) {
            const hicx::CutInterval& x =
                block.cut_intervals[static_cast<std::size_t>(tested.candidates[i].row)];
            const hicx::CutInterval& y =
                block.cut_intervals[static_cast<std::size_t>(tested.candidates[i].col)];
            const std::int64_t distance = std::abs(x.start - y.start);
            if (distance > args.max_loop_distance) {
                continue;
            }
            mapped_loops.push_back(Loop{x.chrom, x.start, x.end, y.chrom, y.start, y.end,
                                        tested.pvalues[i]});
        }
    }

    if (mapped_loops.empty() && single_chromosome) {
        // hicDetectLoops.py:1024-1026, the single chromosome branch.
        std::fputs("ERROR:hicexplorer.hicDetectLoops:No loops could be detected. "
                   "Please change your input parameters, use a matrix with a better "
                   "read coverage or contact the develops on "
                   "https://github.com/deeptools/HiCExplorer/issues\n",
                   stderr);
        return 1;
    }

    if (!mapped_loops.empty()) {
        std::ofstream output(args.out_file_name);
        if (!output) {
            std::fprintf(stderr, "hicDetectLoops: error: cannot write %s\n",
                         args.out_file_name.c_str());
            return 1;
        }
        for (const Loop& loop : mapped_loops) {
            output << loop.chrom_x << '\t' << loop.start_x << '\t' << loop.end_x << '\t'
                   << loop.chrom_y << '\t' << loop.start_y << '\t' << loop.end_y << '\t'
                   << hicx::npy::float_repr(loop.pvalue) << '\n';
        }
    }
    std::fprintf(stderr,
                 "INFO:hicexplorer.hicDetectLoops:Number of detected loops for all "
                 "regions: %zu\n",
                 mapped_loops.size());
    hicx::report_resource_usage("hicDetectLoops");
    return 0;
}
