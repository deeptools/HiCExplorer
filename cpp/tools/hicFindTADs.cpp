// Port of hicexplorer/hicFindTADs.py.
//
// The pipeline, and where each piece lives:
//
//   1. load the matrix, optionally restricting and reordering the chromosomes
//   2. drop the diagonal, mask the NaN bins           (this file)
//   3. z-score the matrix per chromosome              (hicx::convert_to_obs_exp_matrix)
//   4. close the gaps the masking left                (hicx::enlarge_bins)
//   5. truncate to 2 * maxDepth around the diagonal   (this file)
//   6. the TAD-separation score per bin and window    (this file, threaded)
//   7. local minima of the mean score                 (peakdetect, this file)
//   8. a Wilcoxon rank sum p-value per minimum        (hicx::stats::ranksums, threaded)
//   9. FDR or Bonferroni, then the call sets          (hicx::stats)
//
// Threading. The Python uses multiprocessing.Pool over a contiguous range of
// bins per worker (hicFindTADs.py:1103-1109). This port uses threads over the
// same fixed partition of the bin index range, writes every result into its
// own preallocated slot and combines the slots in index order, so the output
// does not depend on the thread count (cpp/OPTIMIZATION.md 3). The p-values
// are threaded the same way. Nothing is accumulated into a shared float.
//
// Behaviour of the Python that is reproduced deliberately, not fixed:
//
//   * hicFindTADs.py:902-920 assigns chr_end_idx = chr_start_idx and then
//     mutates it in place, so the chromosome *start* indices are lost and the
//     end indices are added to the boundary list twice.
//   * delta_wrt_window (:612-613) averages matrix_avg[i-w : i+3] together with
//     matrix_avg[i+4 : i+w], which includes the minimum itself and the two
//     bins after it and skips the bin at i+3, although the docstring says the
//     minimum is excluded.
//   * get_incremental_step_size grows as min + int(step * k ** 1.5), which is
//     not the progression the --step help text describes.
//   * hiCMatrix caches its bin size at construction time, so the bin size used
//     after the bins have been enlarged is still the one from the input.
//   * the z-score matrix keeps NaN entries that lie beyond the truncation
//     distance, because NaN - NaN is NaN and scipy stores every result that is
//     not exactly zero.
//
// One behaviour is deliberately *not* reproduced, and it is a defect rather
// than a numeric difference. np.array_split hands every worker a contiguous
// range of bins, and a worker whose bins are all skipped reaches
// `zip(*positions_array)` with an empty list and ends the whole run with
// "ValueError: not enough values to unpack". On small_test_matrix.h5 that
// happens from --numberOfProcessors 12 upwards while 1 to 10 succeed and give
// byte-identical output, so in the Python the process count decides whether
// the tool runs at all. Here an empty partition simply contributes no rows.
// Pinned by
// test_hicFindTADs.py::test_find_TADs_high_numberOfProcessors_is_a_defect.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
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
#include "hicx/bins.hpp"
#include "hicx/matrix_ops.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/obsexp_ops.hpp"
#include "hicx/parallel.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/simd_reduce.hpp"
#include "hicx/stats_ops.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

const char* const kUsage =
    "usage: hicFindTADs --matrix MATRIX --outPrefix OUTPREFIX\n"
    "                   --correctForMultipleTesting {fdr,bonferroni,None}\n"
    "                   [--minDepth INT bp] [--maxDepth INT bp] [--step INT bp]\n"
    "                   [--TAD_sep_score_prefix TAD_SEP_SCORE_PREFIX]\n"
    "                   [--thresholdComparisons THRESHOLDCOMPARISONS]\n"
    "                   [--delta DELTA] [--minBoundaryDistance MINBOUNDARYDISTANCE]\n"
    "                   [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]\n"
    "                   [--numberOfProcessors NUMBEROFPROCESSORS] [--help]\n"
    "                   [--version]\n";

const char* const kHelp =
    "\n"
    "Uses a measure called TAD-separation score to identify the degree of "
    "separation between\n"
    "the left and right regions at each Hi-C matrix bin. This is done for a\n"
    "running window of different sizes. Then, TADs are called as those\n"
    "positions having a local TAD-separation score minimum.\n"
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        Corrected Hi-C matrix to use for the computations.\n"
    "  --outPrefix OUTPREFIX\n"
    "                        File prefix to save the resulting files.\n"
    "  --correctForMultipleTesting {fdr,bonferroni,None}\n"
    "                        Select the bonferroni or false discovery rate for a\n"
    "                        multiple comparison. (Default: fdr).\n"
    "\n"
    "Optional arguments:\n"
    "  --minDepth INT bp     Minimum window length (in bp).\n"
    "  --maxDepth INT bp     Maximum window length (in bp).\n"
    "  --step INT bp         Step size when moving from --minDepth to --maxDepth.\n"
    "  --TAD_sep_score_prefix TAD_SEP_SCORE_PREFIX\n"
    "                        Prefix of an existing TAD-separation score.\n"
    "  --thresholdComparisons THRESHOLDCOMPARISONS\n"
    "                        P-value threshold for the Bonferroni correction /\n"
    "                        q-value for FDR (Default: 0.01).\n"
    "  --delta DELTA         Minimum threshold of the difference between the\n"
    "                        TAD-separation score of a putative boundary and the\n"
    "                        mean of the TAD-sep. score of surrounding bins\n"
    "                        (Default: 0.01).\n"
    "  --minBoundaryDistance MINBOUNDARYDISTANCE\n"
    "                        Minimum distance between boundaries (in bp).\n"
    "  --chromosomes CHROMOSOMES [CHROMOSOMES ...]\n"
    "                        Chromosomes and order in which the chromosomes should\n"
    "                        be plotted.\n"
    "  --numberOfProcessors NUMBEROFPROCESSORS, -p NUMBEROFPROCESSORS\n"
    "                        Number of processors to use (Default: 1).\n"
    "  --help, -h            show this help message and exit.\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::string matrix;
    std::string out_prefix;
    std::string correct_for_multiple_testing;
    std::optional<std::int64_t> min_depth;
    std::optional<std::int64_t> max_depth;
    std::optional<std::int64_t> step;
    std::optional<std::string> tad_sep_score_prefix;
    double threshold_comparisons = 0.01;
    double delta = 0.01;
    std::optional<std::int64_t> min_boundary_distance;
    std::optional<std::vector<std::string>> chromosomes;
    int processors = 1;
};

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicFindTADs: error: %s\n", message.c_str());
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
    bool prefix_seen = false;
    bool correction_seen = false;

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
            std::printf("hicFindTADs %s\n", hicx::kVersion);
            std::exit(0);
        }
        if (name == "-m" || name == "--matrix") {
            args.matrix = next_value("--matrix/-m");
            matrix_seen = true;
        } else if (name == "--outPrefix") {
            args.out_prefix = next_value("--outPrefix");
            prefix_seen = true;
        } else if (name == "--correctForMultipleTesting") {
            args.correct_for_multiple_testing = next_value("--correctForMultipleTesting");
            if (args.correct_for_multiple_testing != "fdr" &&
                args.correct_for_multiple_testing != "bonferroni" &&
                args.correct_for_multiple_testing != "None") {
                fail("argument --correctForMultipleTesting: invalid choice: '" +
                     args.correct_for_multiple_testing +
                     "' (choose from 'fdr', 'bonferroni', 'None')");
            }
            correction_seen = true;
        } else if (name == "--minDepth") {
            args.min_depth = parse_int("--minDepth", next_value("--minDepth"));
        } else if (name == "--maxDepth") {
            args.max_depth = parse_int("--maxDepth", next_value("--maxDepth"));
        } else if (name == "--step") {
            args.step = parse_int("--step", next_value("--step"));
        } else if (name == "--TAD_sep_score_prefix") {
            args.tad_sep_score_prefix = next_value("--TAD_sep_score_prefix");
        } else if (name == "--thresholdComparisons") {
            args.threshold_comparisons =
                parse_double("--thresholdComparisons", next_value("--thresholdComparisons"));
        } else if (name == "--delta") {
            args.delta = parse_double("--delta", next_value("--delta"));
        } else if (name == "--minBoundaryDistance") {
            args.min_boundary_distance =
                parse_int("--minBoundaryDistance", next_value("--minBoundaryDistance"));
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
        } else if (name == "-p" || name == "--numberOfProcessors") {
            args.processors = static_cast<int>(
                parse_int("--numberOfProcessors", next_value("--numberOfProcessors")));
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
    require(prefix_seen, "--outPrefix");
    require(correction_seen, "--correctForMultipleTesting");
    if (!missing.empty()) {
        fail("the following arguments are required: " + missing);
    }
    return args;
}

bool file_exists(const std::string& path) {
    std::ifstream probe(path);
    return probe.good();
}

// Python compares two lists of cut interval tuples with ==, which tests
// identity before value, so a NaN in `extra` that came from the same object
// compares equal. enlarge_bins reuses the extra of the tuple it rewrites, so
// bitwise equality is the right test here.
bool same_intervals(const std::vector<hicx::CutInterval>& a,
                    const std::vector<hicx::CutInterval>& b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].chrom != b[i].chrom || a[i].start != b[i].start ||
            a[i].end != b[i].end || a[i].extra_text != b[i].extra_text) {
            return false;
        }
        const double x = a[i].extra;
        const double y = b[i].extra;
        if (!(x == y) && !(std::isnan(x) && std::isnan(y))) {
            return false;
        }
    }
    return true;
}

// get_incremental_step_size (hicFindTADs.py:283-301).
std::vector<std::int64_t> incremental_step_size(std::int64_t min_window,
                                                std::int64_t max_window,
                                                std::int64_t start_step) {
    std::vector<std::int64_t> steps;
    std::int64_t step = -1;
    while (true) {
        ++step;
        const std::int64_t increment =
            min_window + static_cast<std::int64_t>(static_cast<double>(start_step) *
                                                   std::pow(static_cast<double>(step), 1.5));
        if (step > 1 && !steps.empty() && increment == steps.back()) {
            continue;
        }
        if (increment > max_window) {
            break;
        }
        steps.push_back(increment);
    }
    return steps;
}

// --------------------------------------------------------------------------
// The matrix operations hiCMatrix performs that no core component covers yet.

// hiCMatrix.diagflat(value=0): the diagonal is replaced by zeros, and scipy's
// sparse addition drops results that are exactly zero, so the entries simply
// disappear.
void drop_diagonal(hicx::CsrMatrix& matrix) {
    const std::int64_t n = matrix.rows();
    std::vector<std::int64_t> indptr(static_cast<std::size_t>(n) + 1, 0);
    std::vector<std::int32_t> indices;
    std::vector<double> values;
    indices.reserve(matrix.stored_nnz());
    values.reserve(matrix.stored_nnz());
    for (std::int64_t row = 0; row < n; ++row) {
        const std::size_t begin =
            static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            if (static_cast<std::int64_t>(matrix.indices()[k]) == row) {
                continue;
            }
            indices.push_back(matrix.indices()[k]);
            values.push_back(matrix.data()[k]);
        }
        indptr[static_cast<std::size_t>(row) + 1] =
            static_cast<std::int64_t>(values.size());
    }
    const hicx::Symmetry symmetry = matrix.symmetry();
    matrix = hicx::CsrMatrix(n, n, std::move(indptr), std::move(indices),
                             std::move(values), matrix.dtype());
    matrix.set_symmetry(symmetry);
}

// The bookkeeping half of hiCMatrix.maskBins / restoreMaskedBins. The Python
// keeps the removed bins so that save() can put them back as empty rows, and
// hicFindTADs relies on that: when enlarge_bins does not change the bin table,
// the masked bins reappear in the saved z-score matrix.
struct MaskState {
    std::vector<std::int64_t> orig_bin_ids;
    std::vector<hicx::CutInterval> orig_cut_intervals;
};

MaskState mask_bins(hicx::MatrixData& data, const std::vector<std::int64_t>& bin_ids) {
    MaskState state;
    if (bin_ids.empty()) {
        return state;  // maskBins returns before touching anything
    }
    std::set<std::int64_t> masked(bin_ids.begin(), bin_ids.end());
    // The union with the NaN bins the matrix already carries
    // (HiCMatrix.py:794-799).
    masked.insert(data.nan_bins.begin(), data.nan_bins.end());
    data.nan_bins.clear();

    std::vector<std::int64_t> rows;
    rows.reserve(static_cast<std::size_t>(data.matrix.rows()));
    for (std::int64_t bin = 0; bin < data.matrix.rows(); ++bin) {
        if (masked.count(bin) == 0) {
            rows.push_back(bin);
        }
    }
    data.matrix = hicx::select_bins(data.matrix, rows);

    state.orig_bin_ids = rows;
    state.orig_bin_ids.insert(state.orig_bin_ids.end(), masked.begin(), masked.end());
    state.orig_cut_intervals.reserve(state.orig_bin_ids.size());
    for (const std::int64_t bin : state.orig_bin_ids) {
        state.orig_cut_intervals.push_back(
            data.cut_intervals[static_cast<std::size_t>(bin)]);
    }
    std::vector<hicx::CutInterval> kept;
    kept.reserve(rows.size());
    for (const std::int64_t bin : rows) {
        kept.push_back(data.cut_intervals[static_cast<std::size_t>(bin)]);
    }
    data.cut_intervals = std::move(kept);
    if (data.correction_factors.has_value()) {
        std::vector<double> factors;
        factors.reserve(rows.size());
        for (const std::int64_t bin : rows) {
            factors.push_back((*data.correction_factors)[static_cast<std::size_t>(bin)]);
        }
        data.correction_factors = std::move(factors);
    }
    return state;
}

void restore_masked_bins(hicx::MatrixData& data, MaskState& state) {
    if (state.orig_bin_ids.empty()) {
        return;
    }
    const std::int64_t m = data.matrix.rows();
    const std::int64_t total = static_cast<std::int64_t>(state.orig_bin_ids.size());
    const std::int64_t added = total - m;

    // argsort(orig_bin_ids), then matrix[rows, :][:, rows] on the matrix padded
    // with `added` empty rows and columns.
    std::vector<std::int64_t> order(static_cast<std::size_t>(total));
    for (std::int64_t i = 0; i < total; ++i) {
        order[static_cast<std::size_t>(i)] = i;
    }
    std::stable_sort(order.begin(), order.end(),
                     [&state](std::int64_t a, std::int64_t b) {
                         return state.orig_bin_ids[static_cast<std::size_t>(a)] <
                                state.orig_bin_ids[static_cast<std::size_t>(b)];
                     });

    if (added > 0) {
        // Pad: the extra rows and columns hold nothing, so only the shape and
        // the row offsets change.
        std::vector<std::int64_t> indptr(static_cast<std::size_t>(total) + 1, 0);
        for (std::int64_t row = 0; row <= m; ++row) {
            indptr[static_cast<std::size_t>(row)] =
                data.matrix.indptr()[static_cast<std::size_t>(row)];
        }
        for (std::int64_t row = m + 1; row <= total; ++row) {
            indptr[static_cast<std::size_t>(row)] =
                static_cast<std::int64_t>(data.matrix.stored_nnz());
        }
        std::vector<std::int32_t> indices(data.matrix.indices());
        std::vector<double> values(data.matrix.data());
        const hicx::Symmetry symmetry = data.matrix.symmetry();
        data.matrix = hicx::CsrMatrix(total, total, std::move(indptr), std::move(indices),
                                      std::move(values), "float64");
        data.matrix.set_symmetry(symmetry);
    }
    data.matrix = hicx::select_bins(data.matrix, order);

    std::vector<hicx::CutInterval> restored;
    restored.reserve(order.size());
    for (const std::int64_t index : order) {
        restored.push_back(state.orig_cut_intervals[static_cast<std::size_t>(index)]);
    }
    data.cut_intervals = std::move(restored);
    data.nan_bins.assign(state.orig_bin_ids.begin() + static_cast<std::ptrdiff_t>(m),
                         state.orig_bin_ids.end());
    std::sort(data.nan_bins.begin(), data.nan_bins.end());

    if (data.correction_factors.has_value()) {
        std::vector<double> padded = *data.correction_factors;
        padded.resize(static_cast<std::size_t>(total), kNaN);
        std::vector<double> reordered;
        reordered.reserve(padded.size());
        for (const std::int64_t index : order) {
            reordered.push_back(padded[static_cast<std::size_t>(index)]);
        }
        data.correction_factors = std::move(reordered);
    }
    state.orig_bin_ids.clear();
    state.orig_cut_intervals.clear();
}

// triu(m, 0) - triu(m, limit): entries at or beyond `limit` cancel, except
// NaN, which cannot cancel and stays stored.
void truncate_to_distance(hicx::CsrMatrix& matrix, std::int64_t limit) {
    const std::int64_t n = matrix.rows();
    std::vector<std::int64_t> indptr(static_cast<std::size_t>(n) + 1, 0);
    std::vector<std::int32_t> indices;
    std::vector<double> values;
    indices.reserve(matrix.stored_nnz());
    values.reserve(matrix.stored_nnz());
    for (std::int64_t row = 0; row < n; ++row) {
        const std::size_t begin =
            static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t column = static_cast<std::int64_t>(matrix.indices()[k]);
            if (column < row) {
                continue;
            }
            const double value = matrix.data()[k];
            if (column - row >= limit && !std::isnan(value)) {
                continue;
            }
            if (value == 0.0) {
                continue;  // eliminate_zeros
            }
            indices.push_back(static_cast<std::int32_t>(column));
            values.push_back(value);
        }
        indptr[static_cast<std::size_t>(row) + 1] =
            static_cast<std::int64_t>(values.size());
    }
    matrix = hicx::CsrMatrix(n, n, std::move(indptr), std::move(indices),
                             std::move(values), "float64");
    matrix.set_symmetry(hicx::Symmetry::Full);
}

// --------------------------------------------------------------------------
// The TAD-separation score

// A dense block of the matrix, row major, exactly what
// matrix[r0:r1, c0:c1].todense() produces. `stored` is the .nnz of the slice.
struct DenseBlock {
    std::vector<double> values;
    std::size_t stored = 0;
    std::int64_t rows = 0;
    std::int64_t columns = 0;
};

void extract_block(const hicx::CsrMatrix& matrix, std::int64_t r0, std::int64_t r1,
                   std::int64_t c0, std::int64_t c1, DenseBlock& block) {
    block.rows = std::max<std::int64_t>(0, r1 - r0);
    block.columns = std::max<std::int64_t>(0, c1 - c0);
    const std::size_t count =
        static_cast<std::size_t>(block.rows) * static_cast<std::size_t>(block.columns);
    block.values.assign(count, 0.0);
    block.stored = 0;
    for (std::int64_t row = r0; row < r1; ++row) {
        const std::size_t begin =
            static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t column = static_cast<std::int64_t>(matrix.indices()[k]);
            if (column < c0) {
                continue;
            }
            if (column >= c1) {
                break;  // column indices are sorted
            }
            block.values[static_cast<std::size_t>(row - r0) *
                             static_cast<std::size_t>(block.columns) +
                         static_cast<std::size_t>(column - c0)] = matrix.data()[k];
            ++block.stored;
        }
    }
}

// The state one TAD-score computation needs: the matrix, the bin table and the
// chromosome end positions.
struct ScoreContext {
    const hicx::CsrMatrix* matrix = nullptr;
    const std::vector<hicx::CutInterval>* cut_intervals = nullptr;
    const hicx::BinTable* bins = nullptr;
    std::unordered_map<std::string, std::int64_t> chromosome_sizes;
};

// get_idx_of_bins_at_given_distance. Returns false where the Python's
// getRegionBinRange returns None and the caller ends up with a TypeError.
bool bins_at_distance(const ScoreContext& context, std::int64_t bin,
                      std::int64_t window_length, std::int64_t& left,
                      std::int64_t& right) {
    const hicx::CutInterval& interval =
        (*context.cut_intervals)[static_cast<std::size_t>(bin)];
    const std::int64_t left_start = std::max<std::int64_t>(0, interval.start - window_length);
    const auto left_range =
        context.bins->region_bin_range(interval.chrom, left_start, left_start + 1);
    if (!left_range.has_value()) {
        return false;
    }
    const auto size = context.chromosome_sizes.find(interval.chrom);
    if (size == context.chromosome_sizes.end()) {
        return false;
    }
    const std::int64_t right_end =
        std::min(size->second, interval.end + window_length) - 1;
    const auto right_range =
        context.bins->region_bin_range(interval.chrom, right_end, right_end);
    if (!right_range.has_value()) {
        return false;
    }
    left = left_range->first;
    right = right_range->first;
    return true;
}

// get_cut_weight(..., return_mean=True). `defined` is false where the Python
// returns None.
double cut_weight_mean(const ScoreContext& context, std::int64_t cut,
                       std::int64_t window_length, DenseBlock& scratch, bool& defined) {
    defined = true;
    if (cut < 0 || cut > context.matrix->rows()) {
        defined = false;
        return kNaN;
    }
    std::int64_t left = 0;
    std::int64_t right = 0;
    if (!bins_at_distance(context, cut, window_length, left, right)) {
        defined = false;  // the TypeError branch of get_cut_weight
        return kNaN;
    }
    extract_block(*context.matrix, left, cut, cut, right, scratch);
    if (scratch.stored == 0) {
        return 0.0;  // the Python returns the integer 0 here
    }
    // The hot reduction of the whole tool: one dense block per bin and window
    // size. hicx::simd::pairwise_sum is bit-identical to the numpy compatible
    // scalar reduction, so this stays byte-identical to the Python.
    const std::size_t count = scratch.values.size();
    return hicx::simd::pairwise_sum(scratch.values.data(), count) /
           static_cast<double>(count);
}

// get_cut_weight(..., return_mean=False): the flattened dense block.
bool cut_weight_values(const ScoreContext& context, std::int64_t cut,
                       std::int64_t window_length, DenseBlock& scratch) {
    if (cut < 0 || cut > context.matrix->rows()) {
        return false;
    }
    std::int64_t left = 0;
    std::int64_t right = 0;
    if (!bins_at_distance(context, cut, window_length, left, right)) {
        return false;
    }
    extract_block(*context.matrix, left, cut, cut, right, scratch);
    return true;
}

// --------------------------------------------------------------------------
// peakdetect (hicFindTADs.py:466-572)

struct Peak {
    std::int64_t position = 0;
    double value = 0.0;
};

void peakdetect(const std::vector<double>& y, const std::vector<std::string>& chrom,
                std::int64_t lookahead, std::vector<Peak>& maxima,
                std::vector<Peak>& minima) {
    maxima.clear();
    minima.clear();
    std::vector<bool> dump;
    const std::int64_t n = static_cast<std::int64_t>(y.size());
    const std::int64_t limit = n - lookahead;

    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    std::int64_t max_pos = 0;
    std::int64_t min_pos = 0;
    int search_for = 0;  // 0 none, 1 min, 2 max
    bool have_previous_chrom = false;
    std::string previous_chrom;

    for (std::int64_t index = 0; index < limit; ++index) {
        const double value = y[static_cast<std::size_t>(index)];
        if (!(value > -std::numeric_limits<double>::infinity() &&
              value < std::numeric_limits<double>::infinity())) {
            throw std::runtime_error(
                "Error, infinity value detected for value at position " +
                std::to_string(index));
        }
        if (!have_previous_chrom) {
            previous_chrom = chrom[static_cast<std::size_t>(index)];
            have_previous_chrom = true;
        }
        if (previous_chrom != chrom[static_cast<std::size_t>(index)]) {
            min_y = std::numeric_limits<double>::infinity();
            max_y = -std::numeric_limits<double>::infinity();
            search_for = 0;
        }
        previous_chrom = chrom[static_cast<std::size_t>(index)];

        if (value > max_y) {
            max_y = value;
            max_pos = index;
        }
        if (value < min_y) {
            min_y = value;
            min_pos = index;
        }

        if (value < max_y - 0.0 && max_y != std::numeric_limits<double>::infinity() &&
            search_for != 1) {
            double ahead = -std::numeric_limits<double>::infinity();
            for (std::int64_t k = index; k < index + lookahead && k < n; ++k) {
                ahead = std::max(ahead, y[static_cast<std::size_t>(k)]);
            }
            if (ahead < max_y) {
                maxima.push_back(Peak{max_pos, max_y});
                dump.push_back(true);
                max_y = value;
                min_y = value;
                min_pos = index;
                search_for = 1;
                continue;
            }
        }

        if (value > min_y + 0.0 && min_y != -std::numeric_limits<double>::infinity() &&
            search_for != 2) {
            double ahead = std::numeric_limits<double>::infinity();
            for (std::int64_t k = index; k < index + lookahead && k < n; ++k) {
                ahead = std::min(ahead, y[static_cast<std::size_t>(k)]);
            }
            if (ahead > min_y) {
                minima.push_back(Peak{min_pos, min_y});
                dump.push_back(false);
                min_y = value;
                max_y = value;
                max_pos = index;
                search_for = 2;
            }
        }
    }

    if (!dump.empty()) {
        if (dump[0]) {
            if (!maxima.empty()) {
                maxima.erase(maxima.begin());
            }
        } else if (!minima.empty()) {
            minima.erase(minima.begin());
        }
    }
}

// delta_wrt_window (hicFindTADs.py:574-616).
std::map<std::int64_t, double> delta_wrt_window(const std::vector<std::int64_t>& minima,
                                                const std::vector<double>& scores,
                                                const std::vector<std::string>& chrom,
                                                std::int64_t window_length = 10) {
    // np.unique(chrom, return_index=True) gives the first index of every
    // distinct name; the ranges are then built from the sorted indices plus
    // len(chrom) - 1.
    std::map<std::string, std::int64_t> first_index;
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(chrom.size()); ++i) {
        first_index.emplace(chrom[static_cast<std::size_t>(i)], i);
    }
    std::vector<std::int64_t> boundaries;
    boundaries.reserve(first_index.size() + 1);
    for (const auto& [name, index] : first_index) {
        boundaries.push_back(index);
    }
    boundaries.push_back(static_cast<std::int64_t>(chrom.size()) - 1);
    std::sort(boundaries.begin(), boundaries.end());

    std::map<std::int64_t, double> delta;
    for (const std::int64_t minimum : minima) {
        bool close_to_border = true;
        for (std::size_t r = 0; r + 1 < boundaries.size(); ++r) {
            const std::int64_t start = boundaries[r];
            const std::int64_t end = boundaries[r + 1];
            if (start < minimum && minimum < end) {
                if (minimum - window_length >= start && minimum + window_length < end) {
                    close_to_border = false;
                }
            }
        }
        if (close_to_border) {
            delta[minimum] = kNaN;
            continue;
        }
        // matrix_avg[i - w : i + 3] concatenated with matrix_avg[i + 4 : i + w].
        std::vector<double> window;
        const std::int64_t n = static_cast<std::int64_t>(scores.size());
        const auto append = [&](std::int64_t from, std::int64_t to) {
            from = std::max<std::int64_t>(0, from);
            to = std::min(n, to);
            for (std::int64_t k = from; k < to; ++k) {
                window.push_back(scores[static_cast<std::size_t>(k)]);
            }
        };
        append(minimum - window_length, minimum + 3);
        append(minimum + 4, minimum + window_length);
        const double mean = hicx::npy::pairwise_sum(window.data(), window.size()) /
                            static_cast<double>(window.size());
        delta[minimum] = mean - scores[static_cast<std::size_t>(minimum)];
    }
    return delta;
}

// --------------------------------------------------------------------------
// Output formatting

std::string format_fixed(double value, int decimals) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.*f", decimals, value);
    return buffer;
}

struct BedgraphMatrix {
    std::vector<std::string> chrom;
    std::vector<std::int64_t> start;
    std::vector<std::int64_t> end;
    std::vector<double> values;  // row major, `columns` per row
    std::size_t columns = 0;

    [[nodiscard]] std::size_t rows() const { return chrom.size(); }
    [[nodiscard]] std::vector<double> row_means() const {
        std::vector<double> means(rows(), 0.0);
        for (std::size_t r = 0; r < rows(); ++r) {
            means[r] = hicx::simd::pairwise_sum(values.data() + r * columns, columns) /
                       static_cast<double>(columns);
        }
        return means;
    }
};

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    try {
        // hicFindTADs.py:1319-1325 checks this only after the matrix has been
        // loaded, so the diagnostic never fires; here it fires.
        std::string suffix = args.matrix;
        const std::size_t dot = suffix.rfind('.');
        suffix = dot == std::string::npos ? suffix : suffix.substr(dot + 1);
        if (suffix != "cool" && suffix != "h5" && suffix.find("mcool::") == std::string::npos) {
            std::fprintf(stderr,
                         "ERROR:hicexplorer.hicFindTADs:Could not determine file "
                         "ending. Please use either .h5, .cool or "
                         ".mcool::/path/to/matrix.\n");
            std::fprintf(stderr,
                         "ERROR:hicexplorer.hicFindTADs:used input: %s, ending %s\n",
                         args.matrix.c_str(), suffix.c_str());
            return 1;
        }

        hicx::ToolMatrix hic = hicx::ToolMatrix::load(args.matrix);

        // reorderChromosomes(valid_chromosomes)
        if (args.chromosomes.has_value()) {
            std::vector<std::int64_t> order;
            std::vector<std::string> invalid;
            for (const std::string& chrom : *args.chromosomes) {
                const auto range = std::find_if(
                    hic.boundaries().begin(), hic.boundaries().end(),
                    [&chrom](const auto& entry) { return entry.first == chrom; });
                if (range == hic.boundaries().end()) {
                    invalid.push_back(chrom);
                    continue;
                }
                for (std::int64_t bin = range->second.first; bin < range->second.last;
                     ++bin) {
                    order.push_back(bin);
                }
            }
            if (!invalid.empty()) {
                std::fputs(
                    "WARNING:hicexplorer.hicFindTADs:WARNING: The following "
                    "chromosome/scaffold names were not found. Please check"
                    "the correct spelling of the chromosome names. \n\n",
                    stderr);
                for (const std::string& chrom : invalid) {
                    std::fprintf(stderr, "%s\n", chrom.c_str());
                }
            }
            hicx::reorder_bins(hic.data(), order);
            hic.refresh_boundaries();
        }

        // hiCMatrix caches the bin size at construction, so this value stays in
        // force even after the bins have been enlarged.
        const hicx::BinTable initial_bins(hic.data().cut_intervals);
        const std::int64_t bin_size = initial_bins.bin_size();

        // set_variables (hicFindTADs.py:419-464)
        std::int64_t max_depth = 0;
        std::int64_t min_depth = 0;
        std::int64_t step = 0;
        if (args.max_depth.has_value() && args.min_depth.has_value() &&
            *args.max_depth <= *args.min_depth) {
            std::fputs("ERROR:hicexplorer.hicFindTADs:Please check that maxDepth is "
                       "larger than minDepth.\n",
                       stderr);
            return 0;  // the Python calls exit() with no argument
        }
        if (!args.max_depth.has_value()) {
            max_depth = bin_size < 1000    ? bin_size * 60
                        : bin_size < 20000 ? bin_size * 40
                                           : bin_size * 10;
        } else if (*args.max_depth < bin_size * 5) {
            std::fputs("ERROR:hicexplorer.hicFindTADs:Please specify a --maxDepth "
                       "that is at least 5 times larger than the matrix bin size\n",
                       stderr);
            return 1;
        } else {
            max_depth = *args.max_depth;
        }
        if (!args.min_depth.has_value()) {
            min_depth = bin_size < 1000    ? bin_size * 30
                        : bin_size < 20000 ? bin_size * 10
                                           : bin_size * 5;
        } else if (*args.min_depth < bin_size * 3) {
            std::fputs("ERROR:hicexplorer.hicFindTADs:Please specify a --minDepth "
                       "that is at least 3 times larger than the matrix bin size\n",
                       stderr);
            return 1;
        } else {
            min_depth = *args.min_depth;
        }
        if (!args.step.has_value()) {
            step = bin_size < 1000 ? bin_size * 4 : bin_size * 2;
        } else if (*args.step < bin_size) {
            std::fputs("ERROR:hicexplorer.hicFindTADs:Please specify a --step that "
                       "is at least the size of the matrix bin size\n",
                       stderr);
            return 1;
        } else {
            step = *args.step;
        }

        const std::string score_suffix = "_tad_score.bm";
        std::string tad_score_file = args.out_prefix + score_suffix;
        std::string zscore_matrix_file = args.out_prefix + "_zscore_matrix." + suffix;
        if (args.tad_sep_score_prefix.has_value()) {
            tad_score_file = *args.tad_sep_score_prefix + score_suffix;
            zscore_matrix_file = *args.tad_sep_score_prefix + "_zscore_matrix." + suffix;
            if (!file_exists(tad_score_file)) {
                std::fprintf(stderr,
                             "ERROR:hicexplorer.hicFindTADs:The given "
                             "TAD_sep_score_prefix does not contain a valid "
                             "TAD-separation score. Please check.\nCould not find "
                             "file %s\n",
                             tad_score_file.c_str());
                return 1;
            }
            if (!file_exists(zscore_matrix_file)) {
                std::fprintf(stderr,
                             "ERROR:hicexplorer.hicFindTADs:The given "
                             "TAD_sep_score_prefix does not contain a valid z-score "
                             "matrix. Please check.\nCould not find file %s\n",
                             zscore_matrix_file.c_str());
                return 1;
            }
        }

        BedgraphMatrix bedgraph;
        const bool reuse_existing =
            args.tad_sep_score_prefix.has_value() || file_exists(tad_score_file);

        if (reuse_existing) {
            // set_matrix(zscore_matrix_file) followed by load_bedgraph_matrix.
            hic = hicx::ToolMatrix::load(zscore_matrix_file);
            if (args.chromosomes.has_value()) {
                std::vector<std::int64_t> order;
                for (const std::string& chrom : *args.chromosomes) {
                    const auto range = std::find_if(
                        hic.boundaries().begin(), hic.boundaries().end(),
                        [&chrom](const auto& entry) { return entry.first == chrom; });
                    if (range == hic.boundaries().end()) {
                        continue;
                    }
                    for (std::int64_t bin = range->second.first;
                         bin < range->second.last; ++bin) {
                        order.push_back(bin);
                    }
                }
                hicx::reorder_bins(hic.data(), order);
                hic.refresh_boundaries();
            }

            std::ifstream input(tad_score_file);
            if (!input) {
                std::fprintf(stderr, "hicFindTADs: cannot open %s\n",
                             tad_score_file.c_str());
                return 1;
            }
            std::string line;
            while (std::getline(input, line)) {
                if (!line.empty() && line.back() == '\r') {
                    line.pop_back();
                }
                if (line.empty()) {
                    continue;
                }
                if (line[0] == '#') {
                    // The parameters the spectrum was computed with replace the
                    // command line values (hicFindTADs.py:1157-1160). Only the
                    // four keys the writer emits occur.
                    const auto read = [&line](const char* key) -> std::int64_t {
                        const std::string needle = std::string("\"") + key + "\":";
                        const std::size_t at = line.find(needle);
                        if (at == std::string::npos) {
                            return 0;
                        }
                        return std::strtoll(line.c_str() + at + needle.size(), nullptr, 10);
                    };
                    min_depth = read("minDepth");
                    max_depth = read("maxDepth");
                    step = read("step");
                    continue;
                }
                std::vector<std::string> fields;
                std::size_t start = 0;
                while (true) {
                    const std::size_t tab = line.find('\t', start);
                    if (tab == std::string::npos) {
                        fields.push_back(line.substr(start));
                        break;
                    }
                    fields.push_back(line.substr(start, tab - start));
                    start = tab + 1;
                }
                if (fields.size() < 4) {
                    continue;
                }
                if (args.chromosomes.has_value() &&
                    std::find(args.chromosomes->begin(), args.chromosomes->end(),
                              fields[0]) == args.chromosomes->end()) {
                    continue;
                }
                bedgraph.chrom.push_back(fields[0]);
                bedgraph.start.push_back(
                    static_cast<std::int64_t>(std::strtod(fields[1].c_str(), nullptr)));
                bedgraph.end.push_back(
                    static_cast<std::int64_t>(std::strtod(fields[2].c_str(), nullptr)));
                for (std::size_t c = 3; c < fields.size(); ++c) {
                    bedgraph.values.push_back(std::strtod(fields[c].c_str(), nullptr));
                }
                bedgraph.columns = fields.size() - 3;
            }
        } else {
            // ------------------------------------------------------------------
            // compute_spectra_matrix (hicFindTADs.py:1017-1129)
            std::fputs("INFO:hicexplorer.hicFindTADs:removing diagonal values\n\n",
                       stderr);
            drop_diagonal(hic.matrix());

            MaskState mask = mask_bins(hic.data(), hic.data().nan_bins);
            const std::vector<hicx::CutInterval> original_intervals =
                hic.data().cut_intervals;

            std::fputs("INFO:hicexplorer.hicFindTADs:Computing z-score matrix...\n\n",
                       stderr);
            hicx::ObsExpOptions options;
            options.max_depth_bp = static_cast<double>(max_depth) * 2.5;
            options.zscore = true;
            options.perchr = true;
            hicx::convert_to_obs_exp_matrix(hic.data(), bin_size, options);

            hicx::enlarge_bins(hic.data().cut_intervals);
            if (!same_intervals(hic.data().cut_intervals, original_intervals)) {
                mask.orig_bin_ids.clear();
                mask.orig_cut_intervals.clear();
                hic.data().nan_bins.clear();
            }
            hic.refresh_boundaries();

            std::fputs(
                "INFO:hicexplorer.hicFindTADs:Computing TAD-separation scores...\n\n",
                stderr);
            const std::int64_t min_depth_in_bins = min_depth / bin_size;
            const std::int64_t max_depth_in_bins = max_depth / bin_size;
            if (step / bin_size == 0) {
                std::fprintf(stderr,
                             "ERROR:hicexplorer.hicFindTADs:Please select a step size "
                             "larger than %lld\n",
                             static_cast<long long>(bin_size));
                return 1;
            }
            if (min_depth_in_bins <= 1 || max_depth_in_bins <= 1) {
                std::fputs("ERROR:hicexplorer.hicFindTADs:depth length too small\n",
                           stderr);
                return 0;  // the Python calls exit(0) on both of these
            }
            truncate_to_distance(hic.matrix(), 2 * max_depth_in_bins);

            const std::vector<std::int64_t> windows =
                incremental_step_size(min_depth, max_depth, step);

            const hicx::BinTable bins(hic.data().cut_intervals);
            ScoreContext context;
            context.matrix = &hic.matrix();
            context.cut_intervals = &hic.data().cut_intervals;
            context.bins = &bins;
            for (const auto& [chrom, size] : bins.chromosome_sizes()) {
                context.chromosome_sizes[chrom] = size;
            }

            // bins_to_consider, in chrBinBoundaries order.
            std::vector<std::int64_t> bins_to_consider;
            bins_to_consider.reserve(hic.data().cut_intervals.size());
            for (const auto& [chrom, range] : hic.boundaries()) {
                for (std::int64_t bin = range.first; bin < range.last; ++bin) {
                    bins_to_consider.push_back(bin);
                }
            }

            const std::size_t count = bins_to_consider.size();
            std::vector<char> keep(count, 0);
            std::vector<double> scores(count * windows.size(), 0.0);
            const unsigned int threads =
                static_cast<unsigned int>(std::max(1, args.processors));
            hicx::parallel_for(count, threads, [&](std::size_t index) {
                DenseBlock scratch;
                const std::int64_t bin = bins_to_consider[index];
                bool usable = true;
                for (std::size_t w = 0; w < windows.size(); ++w) {
                    bool defined = true;
                    const double value =
                        cut_weight_mean(context, bin, windows[w], scratch, defined);
                    if (!defined || std::isnan(value)) {
                        usable = false;
                        break;
                    }
                    scores[index * windows.size() + w] = value;
                }
                keep[index] = usable ? 1 : 0;
            });

            bedgraph.columns = windows.size();
            for (std::size_t index = 0; index < count; ++index) {
                if (keep[index] == 0) {
                    continue;
                }
                const hicx::CutInterval& interval =
                    hic.data().cut_intervals[static_cast<std::size_t>(
                        bins_to_consider[index])];
                bedgraph.chrom.push_back(interval.chrom);
                bedgraph.start.push_back(interval.start);
                bedgraph.end.push_back(interval.end);
                bedgraph.values.insert(
                    bedgraph.values.end(), scores.begin() + static_cast<std::ptrdiff_t>(
                                                                index * windows.size()),
                    scores.begin() +
                        static_cast<std::ptrdiff_t>((index + 1) * windows.size()));
            }
            if (bedgraph.rows() == 0) {
                std::fputs("hicFindTADs: no bin produced a TAD-separation score\n",
                           stderr);
                return 1;
            }

            // hiCMatrix.save restores the masked bins first.
            restore_masked_bins(hic.data(), mask);
            hic.refresh_boundaries();
            hic.save(args.out_prefix + "_zscore_matrix." + suffix);

            // save_bedgraph_matrix
            std::ofstream out(tad_score_file);
            out << "#{\"step\":" << step << ",\"minDepth\":" << min_depth
                << ",\"maxDepth\":" << max_depth << ",\"binsize\":" << bin_size << "}\n";
            for (std::size_t r = 0; r < bedgraph.rows(); ++r) {
                out << bedgraph.chrom[r] << '\t' << bedgraph.start[r] << '\t'
                    << bedgraph.end[r];
                for (std::size_t c = 0; c < bedgraph.columns; ++c) {
                    out << '\t'
                        << format_fixed(bedgraph.values[r * bedgraph.columns + c], 6);
                }
                out << '\n';
            }
        }

        // ----------------------------------------------------------------------
        // find_boundaries (hicFindTADs.py:1251-1294)
        std::vector<std::int64_t> widths(bedgraph.rows());
        for (std::size_t r = 0; r < bedgraph.rows(); ++r) {
            widths[r] = bedgraph.end[r] - bedgraph.start[r];
        }
        std::vector<std::int64_t> sorted_widths = widths;
        std::sort(sorted_widths.begin(), sorted_widths.end());
        double average_bin_size = 0.0;
        if (!sorted_widths.empty()) {
            const std::size_t n = sorted_widths.size();
            average_bin_size = n % 2 == 1
                                   ? static_cast<double>(sorted_widths[n / 2])
                                   : 0.5 * (static_cast<double>(sorted_widths[n / 2 - 1]) +
                                            static_cast<double>(sorted_widths[n / 2]));
        }
        const double min_boundary_distance =
            args.min_boundary_distance.has_value()
                ? static_cast<double>(*args.min_boundary_distance)
                : average_bin_size * 4.0;
        const std::int64_t lookahead =
            static_cast<std::int64_t>(min_boundary_distance / average_bin_size);
        if (lookahead < 1) {
            std::fputs("hicFindTADs: minBoundaryDistance must be '1' or above in value\n",
                       stderr);
            return 1;
        }

        const std::vector<double> row_means = bedgraph.row_means();
        std::vector<Peak> maxima;
        std::vector<Peak> minima;
        peakdetect(row_means, bedgraph.chrom, lookahead, maxima, minima);
        std::vector<std::int64_t> minimum_indices;
        minimum_indices.reserve(minima.size());
        for (const Peak& peak : minima) {
            minimum_indices.push_back(peak.position);
        }
        std::map<std::int64_t, double> delta_of_min =
            delta_wrt_window(minimum_indices, row_means, bedgraph.chrom);

        // min_pvalue (hicFindTADs.py:1168-1249)
        const hicx::BinTable bins(hic.data().cut_intervals);
        ScoreContext context;
        context.matrix = &hic.matrix();
        context.cut_intervals = &hic.data().cut_intervals;
        context.bins = &bins;
        for (const auto& [chrom, size] : bins.chromosome_sizes()) {
            context.chromosome_sizes[chrom] = size;
        }

        std::vector<char> usable(minimum_indices.size(), 0);
        std::vector<double> raw_pvalues(minimum_indices.size(), kNaN);
        {
            const unsigned int threads =
                static_cast<unsigned int>(std::max(1, args.processors));
            hicx::parallel_for(
                minimum_indices.size(), threads, [&](std::size_t slot) {
                    const std::int64_t index = minimum_indices[slot];
                    const auto range = bins.region_bin_range(
                        bedgraph.chrom[static_cast<std::size_t>(index)],
                        bedgraph.start[static_cast<std::size_t>(index)],
                        bedgraph.end[static_cast<std::size_t>(index)]);
                    if (!range.has_value()) {
                        return;
                    }
                    usable[slot] = 1;
                    const std::int64_t matrix_index = range->first;
                    std::int64_t left_bin = 0;
                    std::int64_t right_bin = 0;
                    if (!bins_at_distance(context, matrix_index, min_depth, left_bin,
                                          right_bin)) {
                        // The Python unpacks None here and dies with a TypeError.
                        throw std::runtime_error(
                            "bin " + std::to_string(matrix_index) +
                            " has no window at the minimum depth; the Python "
                            "reference raises a bare TypeError here "
                            "(hicFindTADs.py:1201)");
                    }
                    DenseBlock left_block;
                    DenseBlock right_block;
                    DenseBlock boundary_block;
                    const bool has_left =
                        cut_weight_values(context, left_bin, min_depth, left_block);
                    const bool has_right =
                        cut_weight_values(context, right_bin, min_depth, right_block);
                    const bool has_boundary = cut_weight_values(context, matrix_index,
                                                               min_depth, boundary_block);
                    const std::size_t left_size =
                        has_left ? left_block.values.size() : 0;
                    const std::size_t right_size =
                        has_right ? right_block.values.size() : 0;
                    const std::size_t boundary_size =
                        has_boundary ? boundary_block.values.size() : 0;
                    if (left_size == 0 || right_size == 0 || boundary_size == 0) {
                        raw_pvalues[slot] = kNaN;
                        return;
                    }
                    const double p1 =
                        hicx::stats::ranksums(boundary_block.values, left_block.values)
                            .pvalue;
                    const double p2 =
                        hicx::stats::ranksums(boundary_block.values, right_block.values)
                            .pvalue;
                    raw_pvalues[slot] = std::min(p1, p2);
                });
        }

        std::vector<std::int64_t> pvalue_index;
        std::vector<double> pvalues;
        for (std::size_t slot = 0; slot < minimum_indices.size(); ++slot) {
            if (usable[slot] == 0) {
                continue;
            }
            pvalue_index.push_back(minimum_indices[slot]);
            pvalues.push_back(raw_pvalues[slot]);
        }

        double pvalue_fdr = 0.0;
        if (args.correct_for_multiple_testing == "fdr") {
            for (double& p : pvalues) {
                if (std::isnan(p)) {
                    p = 1.0;
                }
            }
            pvalue_fdr = hicx::stats::benjamini_hochberg_cutoff(
                pvalues, args.threshold_comparisons);
        } else if (args.correct_for_multiple_testing == "bonferroni") {
            hicx::stats::bonferroni_in_place(pvalues);
        }
        std::map<std::int64_t, double> pvalue_of_min;
        for (std::size_t i = 0; i < pvalue_index.size(); ++i) {
            pvalue_of_min.emplace(pvalue_index[i], pvalues[i]);
        }

        if (minimum_indices.empty()) {
            std::fputs("ERROR:hicexplorer.hicFindTADs:\n*ERROR*\nNo boundaries were "
                       "found.\n",
                       stderr);
            return 1;
        }

        // ----------------------------------------------------------------------
        // save_domains_and_boundaries (hicFindTADs.py:888-1015)
        const std::size_t rows = bedgraph.rows();
        std::map<std::string, std::int64_t> first_index;
        for (std::size_t i = 0; i < rows; ++i) {
            first_index.emplace(bedgraph.chrom[i], static_cast<std::int64_t>(i));
        }
        std::vector<std::int64_t> chr_end_idx;
        if (first_index.size() == 1) {
            chr_end_idx.push_back(static_cast<std::int64_t>(rows) - 1);
        } else {
            for (const auto& [name, index] : first_index) {
                chr_end_idx.push_back(index == 0 ? static_cast<std::int64_t>(rows) - 1
                                                 : index - 1);
            }
        }
        // chr_end_idx is the same array object as chr_start_idx in the Python,
        // so the start indices are lost and the end indices enter the list
        // twice. Reproduced.
        std::vector<std::int64_t> candidates;
        candidates.insert(candidates.end(), chr_end_idx.begin(), chr_end_idx.end());
        candidates.insert(candidates.end(), chr_end_idx.begin(), chr_end_idx.end());
        candidates.insert(candidates.end(), minimum_indices.begin(),
                          minimum_indices.end());
        std::sort(candidates.begin(), candidates.end());

        std::vector<std::int64_t> filtered;
        for (const std::int64_t index : candidates) {
            const auto delta_entry = delta_of_min.find(index);
            const double delta =
                delta_entry == delta_of_min.end() ? kNaN : delta_entry->second;
            if (delta_entry == delta_of_min.end()) {
                delta_of_min[index] = kNaN;
            }
            const auto pvalue_entry = pvalue_of_min.find(index);
            if (pvalue_entry == pvalue_of_min.end()) {
                continue;
            }
            const double threshold = args.correct_for_multiple_testing == "fdr"
                                         ? pvalue_fdr
                                         : args.threshold_comparisons;
            if (delta >= args.delta && pvalue_entry->second <= threshold) {
                filtered.push_back(index);
            }
        }

        const std::set<std::int64_t> chromosome_ends(chr_end_idx.begin(),
                                                     chr_end_idx.end());
        const std::string delta_text = hicx::npy::float_repr(args.delta);

        std::ofstream boundaries_bed(args.out_prefix + "_boundaries.bed");
        std::ofstream domains_bed(args.out_prefix + "_domains.bed");
        std::ofstream gff(args.out_prefix + "_boundaries.gff");
        int count = 1;
        for (std::size_t position = 0; position < filtered.size(); ++position) {
            const std::int64_t index = filtered[position];
            if (chromosome_ends.count(index) != 0) {
                continue;
            }
            const std::size_t at = static_cast<std::size_t>(index);
            // chr_start[min_bin_id - 1] with min_bin_id == 0 is Python's
            // negative index, that is the last row. It cannot be reached in
            // practice because delta_wrt_window gives index 0 a NaN delta and
            // NaN fails the delta filter, but the wraparound is reproduced
            // rather than left as undefined behaviour.
            const std::size_t before = at == 0 ? rows - 1 : at - 1;
            const std::int64_t right_center =
                bedgraph.start[at] + (bedgraph.end[at] - bedgraph.start[at]) / 2;
            const std::int64_t left_center =
                bedgraph.start[before] +
                (bedgraph.end[before] - bedgraph.start[before]) / 2;
            if (bedgraph.chrom[at] != bedgraph.chrom[before]) {
                continue;
            }
            char identifier[16];
            std::snprintf(identifier, sizeof(identifier), "B%05lld",
                          static_cast<long long>(index));
            boundaries_bed << bedgraph.chrom[at] << '\t' << left_center << '\t'
                           << right_center << '\t' << identifier << '\t'
                           << format_fixed(row_means[at], 12) << "\t.\n";
            gff << bedgraph.chrom[at] << "\tHiCExplorer\tboundary\t" << left_center
                << '\t' << right_center << '\t' << format_fixed(row_means[at], 12)
                << "\t.\t.\tID=" << identifier
                << ";delta=" << format_fixed(delta_of_min[index], 12)
                << ";pvalue=" << format_fixed(pvalue_of_min[index], 12)
                << ";tad_sep=" << format_fixed(row_means[at], 12) << '\n';

            if (position + 1 == filtered.size() ||
                bedgraph.chrom[at] !=
                    bedgraph.chrom[static_cast<std::size_t>(filtered[position + 1])]) {
                continue;
            }
            const std::int64_t start = bedgraph.start[at];
            const std::int64_t end =
                bedgraph.start[static_cast<std::size_t>(filtered[position + 1])];
            const char* rgb = count % 2 == 0 ? "51,160,44" : "31,120,180";
            domains_bed << bedgraph.chrom[at] << '\t' << start << '\t' << end << "\tID_"
                        << delta_text << '_' << count << '\t'
                        << format_fixed(row_means[at], 12) << "\t.\t" << start << '\t'
                        << end << '\t' << rgb << '\n';
            ++count;
        }
        boundaries_bed.close();
        domains_bed.close();
        gff.close();

        std::ofstream score(args.out_prefix + "_score.bedgraph");
        for (std::size_t index = 1; index < rows; ++index) {
            const std::int64_t right_center =
                bedgraph.start[index] +
                (bedgraph.end[index] - bedgraph.start[index]) / 2;
            const std::int64_t left_center =
                bedgraph.start[index - 1] +
                (bedgraph.end[index - 1] - bedgraph.start[index - 1]) / 2;
            if (right_center <= left_center) {
                continue;
            }
            score << bedgraph.chrom[index] << '\t' << left_center << '\t' << right_center
                  << '\t' << format_fixed(row_means[index], 12) << '\n';
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicFindTADs: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicFindTADs");
    return 0;
}
