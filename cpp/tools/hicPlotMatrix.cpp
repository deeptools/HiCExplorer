// Port of hicexplorer/hicPlotMatrix.py (cpp/PLAN.md tier 7, option (a)).
//
// The compute runs here: loading the matrix (a cool region or chromosome
// through cooler's fetch when the reference takes that route), --clearMaskedBins
// (maskBins and enlarge_bins), --chromosomeOrder (reorderChromosomes),
// --region and --region2 (translate_region and getRegion, whose bin
// selection differs between cool and h5 inputs), the dense matrices pcolormesh
// draws (getMatrix fills the NaN bins with NaN on the whole-matrix path), the
// replacement of zeros, NaN and infinite values for --log and --log1p, the
// +1 of --log1p, the bin start positions and the chromosome extents.
// plot/hicexplorer_plot/hicPlotMatrix.py draws the figure with the
// reference's plotHeatmap, bigwig, loop and TAD functions; the bigwig tracks,
// loops and TADs are read there, as the reference reads them while drawing.
//
// What is reproduced, pinned by the harness cases:
//
//  1. A cool input without --region2 and with at most one --chromosomeOrder
//     entry is loaded as the region or chromosome only, and the dense matrix
//     is that whole block even when getRegion selects fewer bins.
//  2. Otherwise the whole matrix is loaded; with --region the dense matrix is
//     matrix[idx1, :][:, idx2], without it getMatrix() with NaN rows and
//     columns for the NaN bins.
//  3. getRegion keeps, for a cool input, bins that overlap the region by any
//     of three conditions, and for an h5 input bins inside it; --region2 uses
//     a strict end.
//  4. --log and --log1p on the whole matrix: zeros take the smallest non-zero
//     value, then NaN and infinite values the smallest value that is not NaN,
//     not infinite. Per chromosome all three masks are taken before any
//     replacement, and a failing replacement stops the rest silently.
//  5. On the whole-matrix path --chromosomeOrder becomes bytes, so loops and
//     TADs are filtered against bytes names (the drawing module reproduces
//     that from a flag).
//
// Memory: the dense matrices are what pcolormesh needs, so they are built
// here and written to .npy files; the sparse matrix is freed before the
// drawing process starts.

#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/bins.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/hic_adapter.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

namespace cli = hicx::cli;
namespace plot = hicx::plot;

const char* const kUsage =
    "usage: hicPlotMatrix --matrix MATRIX --outFileName OUTFILENAME\n"
    "                     [--title TITLE] [--scoreName SCORENAME]\n"
    "                     [--perChromosome] [--clearMaskedBins]\n"
    "                     [--chromosomeOrder CHROMOSOMEORDER [CHROMOSOMEORDER ...]]\n"
    "                     [--region REGION] [--region2 REGION2] [--log1p] [--log]\n"
    "                     [--colorMap COLORMAP] [--vMin VMIN] [--vMax VMAX]\n"
    "                     [--dpi DPI] [--bigwig BIGWIG [BIGWIG ...]]\n"
    "                     [--bigwigAdditionalVerticalAxis]\n"
    "                     [--vMinBigwig VMINBIGWIG] [--vMaxBigwig VMAXBIGWIG]\n"
    "                     [--flipBigwigSign]\n"
    "                     [--scaleFactorBigwig SCALEFACTORBIGWIG]\n"
    "                     [--fontsize FONTSIZE] [--rotationX ROTATIONX]\n"
    "                     [--rotationY ROTATIONY]\n"
    "                     [--increaseFigureWidth INCREASEFIGUREWIDTH]\n"
    "                     [--increaseFigureHeight INCREASEFIGUREHEIGHT]\n"
    "                     [--loops LOOPS]\n"
    "                     [--loopLargeRegionsOperation {first,last,center}]\n"
    "                     [--tads TADS] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Creates a heatmap of a Hi-C matrix.\n"
    "\n"
    "The C++ port computes the matrices, positions and extents of the figure and\n"
    "draws it through the hicexplorer_plot drawing layer with the calls of the\n"
    "Python tool (HICX_PLOT_PYTHON names the interpreter); see the Python tool's\n"
    "help for every option. The C++-only option --plotData FILE writes the data of\n"
    "the figure as JSON to FILE (and its matrices as FILE.<n>.npy) instead of\n"
    "drawing it. The C++-only option --matrix2 MATRIX2 has no equivalent in the\n"
    "Python original (cpp/PLAN.md tier 9): it draws a single heatmap whose upper\n"
    "triangle (column index greater than row index) comes from --matrix and whose\n"
    "lower triangle (column index less than row index) comes from --matrix2, the\n"
    "side-by-side comparison of two matrices used throughout the Hi-C literature.\n"
    "Both matrices are read for the same requested region and must resolve to the\n"
    "same shape there; the diagonal itself is always taken from --matrix. It is\n"
    "incompatible with --perChromosome.\n";

class ExitError : public std::runtime_error {
  public:
    explicit ExitError(const std::string& message, int code = 1)
        : std::runtime_error(message), code_(code) {}
    [[nodiscard]] int code() const noexcept { return code_; }

  private:
    int code_;
};

// The hiCMatrix state hicPlotMatrix works on: the bins and a symmetric sparse
// matrix, as a list of full-matrix entries after any selection.
struct Model {
    std::vector<hicx::CutInterval> intervals;
    hicx::CsrMatrix matrix;  // Full or UpperTriangle, never copied
    // The bins after maskBins and reorderChromosomes, as positions in the
    // loaded matrix: the selections are kept as this view instead of a copy
    // of the matrix, and only the dense extraction reads through it.
    std::vector<std::int64_t> view;
    std::vector<std::int64_t> nan_bins;
};

// Every entry of the symmetric matrix, both triangles.
template <typename Visit>
void for_each_full_entry(const hicx::CsrMatrix& matrix, Visit&& visit) {
    const bool upper = matrix.symmetry() == hicx::Symmetry::UpperTriangle;
    const auto& indptr = matrix.indptr();
    const auto& indices = matrix.indices();
    const auto& data = matrix.data();
    for (std::int64_t row = 0; row < matrix.rows(); ++row) {
        for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
             k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
            const std::int64_t col = indices[k];
            visit(row, col, data[k]);
            if (upper && col != row) {
                visit(col, row, data[k]);
            }
        }
    }
}

// matrix[order, :][:, order] with the cut intervals, as hiCMatrix's
// reorderBins and maskBins select bins; the matrix itself stays as loaded.
void select_bins(Model& model, const std::vector<std::int64_t>& order) {
    std::vector<std::int64_t> view;
    std::vector<hicx::CutInterval> intervals;
    view.reserve(order.size());
    intervals.reserve(order.size());
    for (const std::int64_t bin : order) {
        view.push_back(model.view[static_cast<std::size_t>(bin)]);
        intervals.push_back(model.intervals[static_cast<std::size_t>(bin)]);
    }
    model.view = std::move(view);
    model.intervals = std::move(intervals);
}

// Whether matrix[view, :][:, view] stores any value.
bool stores_any_value(const Model& model) {
    std::vector<char> in_view(static_cast<std::size_t>(model.matrix.rows()), 0);
    for (const std::int64_t bin : model.view) {
        in_view[static_cast<std::size_t>(bin)] = 1;
    }
    bool any = false;
    for_each_full_entry(model.matrix, [&](std::int64_t r, std::int64_t c, double) {
        any = any || (in_view[static_cast<std::size_t>(r)] && in_view[static_cast<std::size_t>(c)]);
    });
    return any;
}

// The chromosome boundaries of the current bins, in order.
std::vector<std::pair<std::string, hicx::BinRange>> boundaries(const Model& model) {
    if (model.intervals.empty()) {
        return {};
    }
    return hicx::BinTable(model.intervals).chrom_bin_boundaries();
}

std::optional<hicx::BinRange> bin_range(const Model& model, const std::string& chrom) {
    for (const auto& entry : boundaries(model)) {
        if (entry.first == chrom) {
            return entry.second;
        }
    }
    return std::nullopt;
}

// hiCMatrix.maskBins(nan_bins): the NaN bins are removed and forgotten.
void mask_nan_bins(Model& model) {
    if (model.nan_bins.empty()) {
        return;
    }
    std::vector<std::int64_t> removed = model.nan_bins;
    std::sort(removed.begin(), removed.end());
    removed.erase(std::unique(removed.begin(), removed.end()), removed.end());
    std::vector<std::int64_t> keep;
    std::size_t r = 0;
    for (std::int64_t bin = 0; bin < static_cast<std::int64_t>(model.intervals.size()); ++bin) {
        while (r < removed.size() && removed[r] < bin) {
            ++r;
        }
        if (r < removed.size() && removed[r] == bin) {
            continue;
        }
        keep.push_back(bin);
    }
    model.nan_bins.clear();
    select_bins(model, keep);
}

// hicexplorer.utilities.enlarge_bins, mutating the list as it walks it.
void enlarge_bins(std::vector<hicx::CutInterval>& bins) {
    if (bins.empty()) {
        return;
    }
    bool chr_start = true;
    for (std::size_t idx = 0; idx + 1 < bins.size(); ++idx) {
        hicx::CutInterval& current = bins[idx];
        hicx::CutInterval& next = bins[idx + 1];
        if (chr_start) {
            current.start = 0;
            chr_start = false;
        }
        if (current.chrom == next.chrom && current.end != next.start) {
            const std::int64_t middle = next.start - static_cast<std::int64_t>(
                                                         static_cast<double>(next.start - current.end) / 2.0);
            current.end = middle;
            next.start = middle;
        }
        if (current.chrom != next.chrom) {
            chr_start = true;
        }
    }
}

// hiCMatrix.reorderChromosomes(names).
void reorder_chromosomes(Model& model, const std::vector<std::string>& names) {
    std::vector<std::int64_t> order;
    const auto ranges = boundaries(model);
    for (const std::string& name : names) {
        for (const auto& entry : ranges) {
            if (entry.first == name) {
                for (std::int64_t bin = entry.second.first; bin < entry.second.last; ++bin) {
                    order.push_back(bin);
                }
            }
        }
    }
    std::vector<char> was_nan(model.intervals.size(), 0);
    for (const std::int64_t bin : model.nan_bins) {
        was_nan[static_cast<std::size_t>(bin)] = 1;
    }
    std::vector<std::int64_t> nan_bins;
    for (std::size_t i = 0; i < order.size(); ++i) {
        if (was_nan[static_cast<std::size_t>(order[i])]) {
            nan_bins.push_back(static_cast<std::int64_t>(i));
        }
    }
    select_bins(model, order);
    model.nan_bins = std::move(nan_bins);
}

std::string change_chrom_names(const std::string& chrom) {
    return chrom.rfind("chr", 0) == 0 ? chrom.substr(3) : "chr" + chrom;
}

struct Region {
    std::string chrom;
    std::int64_t start = 0;
    std::int64_t end = 0;
};

// translate_region(region_string, ma).
Region translate_region(std::string text, const Model& model) {
    std::string cleaned;
    for (const char c : text) {
        if (c == ',' || c == ';' || c == '!') {
            continue;
        }
        cleaned += c == '-' ? ':' : c;
    }
    std::vector<std::string> fields;
    std::string current;
    for (const char c : cleaned) {
        if (c == ':') {
            fields.push_back(current);
            current.clear();
        } else {
            current += c;
        }
    }
    fields.push_back(current);
    Region region;
    region.chrom = fields[0];
    std::optional<hicx::BinRange> range = bin_range(model, region.chrom);
    if (!range.has_value()) {
        region.chrom = change_chrom_names(region.chrom);
        range = bin_range(model, region.chrom);
        if (!range.has_value()) {
            throw ExitError("Chromosome name " + change_chrom_names(region.chrom) +
                            " in --region not in matrix");
        }
    }
    auto parse = [](const std::string& value) {
        std::int64_t number = 0;
        if (!cli::python_int(value, &number)) {
            throw ExitError("ValueError: invalid literal for int() with base 10: '" + value + "'");
        }
        return number;
    };
    region.start = fields.size() > 1 ? parse(fields[1])
                                     : model.intervals[static_cast<std::size_t>(range->first)].start;
    region.end = fields.size() > 2 ? parse(fields[2])
                                   : model.intervals[static_cast<std::size_t>(range->last - 1)].end;
    return region;
}

struct Selection {
    Region region;
    std::vector<std::int64_t> idx1;
    std::vector<std::int64_t> start_pos1;
    std::optional<Region> region2;
    std::vector<std::int64_t> idx2;
    std::vector<std::int64_t> start_pos2;
};

// getRegion(args, ma).
Selection get_region(const std::string& region_text, const std::optional<std::string>& region2_text,
                     const Model& model, bool is_cooler) {
    Selection out;
    out.region = translate_region(region_text, model);
    const Region& r = out.region;
    for (std::size_t idx = 0; idx < model.intervals.size(); ++idx) {
        const hicx::CutInterval& x = model.intervals[idx];
        if (x.chrom != r.chrom) {
            continue;
        }
        const bool keep = is_cooler ? ((x.start >= r.start && x.end <= r.end) ||
                                       (x.start < r.end && x.end < r.end && x.end > r.start) ||
                                       (x.start > r.start && x.start < r.end))
                                    : (x.start >= r.start && x.end <= r.end);
        if (keep) {
            out.idx1.push_back(static_cast<std::int64_t>(idx));
            out.start_pos1.push_back(x.start);
        }
    }
    if (out.idx1.empty()) {
        throw ExitError("ValueError: not enough values to unpack (expected 2, got 0) (no bin of "
                        "--region)");
    }
    if (region2_text.has_value() && !region2_text->empty()) {
        out.region2 = translate_region(*region2_text, model);
        const Region& r2 = *out.region2;
        for (std::size_t idx = 0; idx < model.intervals.size(); ++idx) {
            const hicx::CutInterval& x = model.intervals[idx];
            if (x.chrom != r2.chrom) {
                continue;
            }
            const bool keep = is_cooler ? ((x.start >= r2.start && x.end < r2.end) ||
                                           (x.start < r2.end && x.end < r2.end && x.end > r2.start) ||
                                           (x.start > r2.start && x.start < r2.end))
                                        : (x.start >= r2.start && x.end < r2.end);
            if (keep) {
                out.idx2.push_back(static_cast<std::int64_t>(idx));
                out.start_pos2.push_back(x.start);
            }
        }
        if (out.idx2.empty()) {
            throw ExitError("ValueError: not enough values to unpack (expected 2, got 0) (no bin "
                            "of --region2)");
        }
    } else {
        out.idx2 = out.idx1;
        out.start_pos2 = out.start_pos1;
    }
    return out;
}

// matrix[rows, :][:, cols].todense().astype(float).
std::vector<double> dense(const Model& model, const std::vector<std::int64_t>& rows,
                          const std::vector<std::int64_t>& cols) {
    const auto n = static_cast<std::size_t>(model.matrix.rows());
    std::vector<std::vector<std::int64_t>> row_positions(n);
    std::vector<std::vector<std::int64_t>> col_positions(n);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        row_positions[static_cast<std::size_t>(model.view[static_cast<std::size_t>(rows[i])])]
            .push_back(static_cast<std::int64_t>(i));
    }
    for (std::size_t j = 0; j < cols.size(); ++j) {
        col_positions[static_cast<std::size_t>(model.view[static_cast<std::size_t>(cols[j])])]
            .push_back(static_cast<std::int64_t>(j));
    }
    std::vector<double> out(rows.size() * cols.size(), 0.0);
    for_each_full_entry(model.matrix, [&](std::int64_t r, std::int64_t c, double v) {
        for (const std::int64_t i : row_positions[static_cast<std::size_t>(r)]) {
            for (const std::int64_t j : col_positions[static_cast<std::size_t>(c)]) {
                out[static_cast<std::size_t>(i) * cols.size() + static_cast<std::size_t>(j)] += v;
            }
        }
    });
    return out;
}

std::vector<std::int64_t> all_bins(const Model& model) {
    std::vector<std::int64_t> bins(model.intervals.size());
    for (std::size_t i = 0; i < bins.size(); ++i) {
        bins[i] = static_cast<std::int64_t>(i);
    }
    return bins;
}

// np.nanmin over the selected values; nullopt for an empty selection, which
// raises ValueError, and NaN when every value is NaN.
std::optional<double> nanmin(const std::vector<double>& values, const std::vector<char>& excluded) {
    bool any = false;
    double smallest = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (excluded[i]) {
            continue;
        }
        any = true;
        if (!std::isnan(values[i]) && (std::isnan(smallest) || values[i] < smallest)) {
            smallest = values[i];
        }
    }
    if (!any) {
        return std::nullopt;
    }
    return smallest;
}

void replace_masked(std::vector<double>& values, const std::vector<char>& mask, double value) {
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (mask[i]) {
            values[i] = value;
        }
    }
}

// main() for --log and --log1p on the whole-matrix path.
void mask_whole(std::vector<double>& values) {
    std::vector<char> zero(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        zero[i] = values[i] == 0.0 ? 1 : 0;
    }
    const std::optional<double> smallest = nanmin(values, zero);
    replace_masked(values, zero, smallest.has_value() ? *smallest
                                                      : std::numeric_limits<double>::min());
    bool any_nan = false;
    bool any_inf = false;
    for (const double v : values) {
        any_nan = any_nan || std::isnan(v);
        any_inf = any_inf || std::isinf(v);
    }
    if (any_nan || any_inf) {
        std::vector<char> nan_mask(values.size());
        std::vector<char> inf_mask(values.size());
        for (std::size_t i = 0; i < values.size(); ++i) {
            nan_mask[i] = std::isnan(values[i]) ? 1 : 0;
            inf_mask[i] = std::isinf(values[i]) ? 1 : 0;
        }
        const std::optional<double> without_nan = nanmin(values, nan_mask);
        if (!without_nan.has_value()) {
            throw ExitError("ValueError: zero-size array to reduction operation fmin which has no "
                            "identity");
        }
        replace_masked(values, nan_mask, *without_nan);
        const std::optional<double> without_inf = nanmin(values, inf_mask);
        if (!without_inf.has_value()) {
            throw ExitError("ValueError: zero-size array to reduction operation fmin which has no "
                            "identity");
        }
        replace_masked(values, inf_mask, *without_inf);
    }
}

// plotPerChr() for --log and --log1p: the three masks first, then the
// replacements, the first failure ending them.
void mask_per_chromosome(std::vector<double>& values) {
    std::vector<char> zero(values.size());
    std::vector<char> nan_mask(values.size());
    std::vector<char> inf_mask(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        zero[i] = values[i] == 0.0 ? 1 : 0;
        nan_mask[i] = std::isnan(values[i]) ? 1 : 0;
        inf_mask[i] = std::isinf(values[i]) ? 1 : 0;
    }
    for (const std::vector<char>* mask : {&zero, &nan_mask, &inf_mask}) {
        const std::optional<double> smallest = nanmin(values, *mask);
        if (!smallest.has_value()) {
            return;
        }
        replace_masked(values, *mask, *smallest);
    }
}

// chromosome_start_end(ma).
std::string chromosome_start_end_json(const Model& model) {
    std::vector<std::string> items;
    for (const auto& [name, range] : boundaries(model)) {
        items.push_back(plot::json_list(
            {plot::json_string(name),
             plot::json_int(model.intervals[static_cast<std::size_t>(range.first)].start),
             plot::json_int(model.intervals[static_cast<std::size_t>(range.last - 1)].end)}));
    }
    return plot::json_list(items);
}

// make_start_pos_array(ma).
std::vector<std::int64_t> make_start_pos_array(const Model& model) {
    std::vector<std::pair<std::string, std::int64_t>> sizes;
    for (const auto& [name, range] : boundaries(model)) {
        sizes.emplace_back(name, model.intervals[static_cast<std::size_t>(range.last - 1)].end -
                                     model.intervals[static_cast<std::size_t>(range.first)].start);
    }
    auto size_of = [&](const std::string& name) {
        for (const auto& entry : sizes) {
            if (entry.first == name) {
                return entry.second;
            }
        }
        return std::int64_t{0};
    };
    std::vector<std::int64_t> start_pos;
    if (model.intervals.empty()) {
        return start_pos;
    }
    std::string prev_chrom = model.intervals[0].chrom;
    std::int64_t prev_chroms_sum = 0;
    std::int64_t shift = 0;
    for (std::size_t index = 0; index < model.intervals.size(); ++index) {
        const hicx::CutInterval& bin = model.intervals[index];
        if (index == 0 && bin.start != 0) {
            shift = bin.start;
        }
        if (bin.chrom != prev_chrom) {
            prev_chroms_sum += size_of(prev_chrom);
            prev_chrom = bin.chrom;
            shift = bin.start != 0 ? bin.start : 0;
        }
        start_pos.push_back(bin.start - shift + prev_chroms_sum);
    }
    return start_pos;
}

std::string region_json(const Region& region) {
    return plot::json_list({plot::json_string(region.chrom), plot::json_int(region.start),
                            plot::json_int(region.end)});
}

}  // namespace

namespace {

// Whether fopen(path, "w") would succeed, decided without creating or
// truncating anything: an existing file that is not a directory must be
// writable, otherwise its directory must be.
bool writable_without_creating(const std::string& path) {
    if (path.empty() || path.back() == '/') {
        return false;
    }
    struct stat status {};
    if (::stat(path.c_str(), &status) == 0) {
        return !S_ISDIR(status.st_mode) && ::access(path.c_str(), W_OK) == 0;
    }
    if (errno != ENOENT) {
        return false;
    }
    const std::size_t slash = path.rfind('/');
    const std::string directory =
        slash == std::string::npos ? "." : (slash == 0 ? "/" : path.substr(0, slash));
    return ::access(directory.c_str(), W_OK | X_OK) == 0;
}

// The result of loading one matrix into a dense, non-per-chromosome block:
// exactly what --matrix's own load produces in main() below, factored out so
// --matrix2 can be read through the identical region-scoped or whole-matrix
// path (never a separate, naive whole load).
struct DenseLoad {
    std::vector<double> values;
    std::int64_t rows = 0;
    std::int64_t cols = 0;
};

// Loads path for the given request the same way --matrix's own load does
// (is_partial_source fast path first, whole-matrix path otherwise; see the
// top-of-file comment), producing the dense matrix pcolormesh would draw for
// that request. region/region2 are taken by value: the fast path may consume
// them (as main()'s --matrix load does), and each caller needs its own
// request untouched by the other's.
DenseLoad load_dense_matrix(const std::string& path, std::optional<std::string> region,
                            std::optional<std::string> region2,
                            const std::optional<std::vector<std::string>>& chromosome_order,
                            bool clear_masked_bins) {
    const bool is_cooler = hicx::check_cooler(path);
    const bool is_partial_source = is_cooler || hicx::is_hic_path(path);
    const bool open_cooler_chromosome_order =
        !(chromosome_order.has_value() && chromosome_order->size() > 1);

    Model model;
    auto load = [&](const std::optional<std::string>& chrname) {
        hicx::ToolMatrix loaded = hicx::ToolMatrix::load(path, chrname);
        model.intervals = loaded.cut_intervals();
        model.nan_bins = loaded.nan_bins();
        model.matrix = std::move(loaded.matrix());
        model.view.resize(model.intervals.size());
        for (std::size_t i = 0; i < model.view.size(); ++i) {
            model.view[i] = static_cast<std::int64_t>(i);
        }
    };
    auto clear_masked = [&]() {
        if (clear_masked_bins) {
            mask_nan_bins(model);
            enlarge_bins(model.intervals);
        }
    };

    DenseLoad out;
    if (is_partial_source && !(region2.has_value() && !region2->empty()) && open_cooler_chromosome_order) {
        std::optional<std::string> retrieve;
        if (region.has_value() && !region->empty()) {
            retrieve = *region;
        }
        if (chromosome_order.has_value()) {
            region.reset();
            region2.reset();
            retrieve = chromosome_order->front();
        }
        load(retrieve);
        clear_masked();
        const std::vector<std::int64_t> bins = all_bins(model);
        out.values = dense(model, bins, bins);
        out.rows = out.cols = static_cast<std::int64_t>(bins.size());
    } else {
        load(std::nullopt);
        clear_masked();
        if (chromosome_order.has_value()) {
            region.reset();
            region2.reset();
            std::vector<std::string> valid;
            std::vector<std::string> invalid;
            for (const std::string& chrom : *chromosome_order) {
                if (bin_range(model, chrom).has_value()) {
                    valid.push_back(chrom);
                } else {
                    invalid.push_back(chrom);
                }
            }
            if (!invalid.empty()) {
                std::fputs("WARNING: The following chromosome/scaffold names were not found. "
                           "Please checkthe correct spelling of the chromosome names. \n",
                           stderr);
            }
            reorder_chromosomes(model, valid);
        }
        if (!stores_any_value(model)) {
            throw ExitError("ValueError: zero-size array to reduction operation minimum which "
                            "has no identity");
        }
        if (region.has_value() && !region->empty()) {
            const Selection selection = get_region(*region, region2, model, is_partial_source);
            out.values = dense(model, selection.idx1, selection.idx2);
            out.rows = static_cast<std::int64_t>(selection.idx1.size());
            out.cols = static_cast<std::int64_t>(selection.idx2.size());
        } else {
            const std::vector<std::int64_t> bins = all_bins(model);
            out.values = dense(model, bins, bins);
            out.rows = out.cols = static_cast<std::int64_t>(bins.size());
            for (const std::int64_t bin : model.nan_bins) {
                for (std::int64_t k = 0; k < out.rows; ++k) {
                    out.values[static_cast<std::size_t>(bin * out.cols + k)] =
                        std::numeric_limits<double>::quiet_NaN();
                    out.values[static_cast<std::size_t>(k * out.cols + bin)] =
                        std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
    }
    return out;
}

// Combines two same-shaped dense matrices into one for the --matrix2 upper
// or lower triangle heatmap: column index greater than row index (the upper
// triangle) takes upper's value, column index less than row index (the lower
// triangle) takes lower's value. The diagonal (row == column) is always
// taken from upper (--matrix): it is a single bin's self contact, not a
// choice between two off-diagonal directions, so there is no upper/lower
// distinction to make there, and keeping it filled (rather than NaN) avoids
// a blank diagonal line under --log/--log1p, where the diagonal is usually
// the strongest signal in the plot.
void combine_triangles(std::vector<double>& upper, const std::vector<double>& lower,
                       std::int64_t rows, std::int64_t cols) {
    for (std::int64_t row = 0; row < rows; ++row) {
        for (std::int64_t col = 0; col < cols; ++col) {
            if (col < row) {
                const auto index = static_cast<std::size_t>(row * cols + col);
                upper[index] = lower[index];
            }
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    cli::Parser parser("hicPlotMatrix", "Creates a heatmap of a Hi-C matrix.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrix", "-m"})
        .required()
        .input({"h5", "cool", "mcool"})
        .help("Path of the Hi-C matrix to plot.");
    required.add({"--outFileName", "-out"})
        .type("writableFile")
        .check([](const std::string& value) -> std::optional<std::string> {
            // hicexplorer.utilities.writableFile opens the file with 'w' while
            // parsing. The probe decides the same way without creating or
            // truncating it; the file is opened after the drawing environment
            // check, so a refused environment leaves it as it was.
            if (!writable_without_creating(value)) {
                return value + " file can be opened for writing";
            }
            return std::nullopt;
        })
        .required()
        .output({"png", "pdf", "svg"})
        .help("File name to save the image.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--matrix2"})
        .input({"h5", "cool", "mcool"})
        .cpp_only("Draws one heatmap whose upper triangle (column index greater than row "
                  "index) is --matrix and whose lower triangle (column index less than row "
                  "index) is --matrix2, the two-matrix comparison heatmap common in the Hi-C "
                  "literature; has no equivalent in the Python original (cpp/PLAN.md tier 9). "
                  "Both matrices are read for the same requested region and must resolve to "
                  "the same shape there. The diagonal is always taken from --matrix. "
                  "Incompatible with --perChromosome.")
        .help("Path of a second Hi-C matrix; combined with --matrix into a single upper/lower "
              "triangle comparison heatmap.");
    optional.add({"--title", "-t"}).help("Plot title.");
    optional.add({"--scoreName", "-s"}).help("Score name label for the heatmap legend.");
    optional.add({"--perChromosome"})
        .action(cli::Action::StoreTrue)
        .help("Instead of plotting the whole matrix, each chromosome is plotted next to the "
              "other. This parameter is not compatible with --region.");
    optional.add({"--clearMaskedBins"})
        .action(cli::Action::StoreTrue)
        .help("If set, masked bins are removed from the matrix and the nearest bins are extended "
              "to cover the empty space instead of plotting black lines.");
    optional.add({"--chromosomeOrder"})
        .nargs("+")
        .help("Chromosomes and order in which the chromosomes should be plotted. This option "
              "overrides --region and --region2.");
    optional.add({"--region"}).help("Plot only this region. The format is chr:start-end.");
    optional.add({"--region2"}).help("If given, then only the region defined by --region and "
                                     "--region2 is given.");
    optional.add({"--log1p"}).action(cli::Action::StoreTrue).help("Plot the log1p of the matrix values.");
    optional.add({"--log"}).action(cli::Action::StoreTrue).help("Plot the *MINUS* log of the matrix values.");
    optional.add({"--colorMap"})
        .default_value("RdYlBu_r")
        .help("Color map to use for the heatmap (Default: %(default)s).");
    optional.add({"--vMin"}).type("float").help("Minimum score value.");
    optional.add({"--vMax"}).type("float").help("Maximum score value.");
    optional.add({"--dpi"})
        .type("int")
        .default_value(std::int64_t{72})
        .help("Resolution for the image in case theoutput is a raster graphics image (e.g png, "
              "jpg) (Default: %(default)s).");
    optional.add({"--bigwig"})
        .type("str")
        .nargs("+")
        .input({"bw", "bigwig"})
        .help("Bigwig file to plot below the matrix.");
    optional.add({"--bigwigAdditionalVerticalAxis"})
        .action(cli::Action::StoreTrue)
        .help("Add an additional axis to determine the values of a bigwig file in 2D better.");
    optional.add({"--vMinBigwig"}).type("float").help("Minimum score value for bigwig.");
    optional.add({"--vMaxBigwig"}).type("float").help("Maximum score value for bigwig");
    optional.add({"--flipBigwigSign"})
        .action(cli::Action::StoreTrue)
        .help("The sign of the bigwig values are flipped.");
    optional.add({"--scaleFactorBigwig"})
        .type("float")
        .default_value(1.0)
        .help("Scale the values of a bigwig file by the given factor (Default: %(default)s).");
    optional.add({"--fontsize"})
        .type("float")
        .default_value(std::int64_t{10})
        .help("Fontsize in the plot for x and y axis (Default: %(default)s).");
    optional.add({"--rotationX"})
        .type("float")
        .default_value(std::int64_t{0})
        .help("Rotation in degrees for the labels of x axis (Default: %(default)s).");
    optional.add({"--rotationY"})
        .type("float")
        .default_value(std::int64_t{0})
        .help("Rotation in degrees for the labels of y axis (Default: %(default)s).");
    optional.add({"--increaseFigureWidth"})
        .type("float")
        .default_value(0.5)
        .help("Increase the figure width for additional bigwig tracks (Default: %(default)s).");
    optional.add({"--increaseFigureHeight"})
        .type("float")
        .default_value(0.5)
        .help("Increase the figure height for additional bigwig tracks (Default: %(default)s).");
    optional.add({"--loops"})
        .type("str")
        .input({"bedgraph", "bedpe"})
        .help("Bedgraph file to plot detected long range contacts from hicDetectLoops.");
    optional.add({"--loopLargeRegionsOperation"})
        .choices({"first", "last", "center"})
        .default_value("first")
        .help("Which bin of a loop coordinate larger than a bin is plotted.");
    optional.add({"--tads"}).type("str").input({"bed"}).help("Bedgraph file to plot detected tads");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    optional.add({"--disable_tight_layout"})
        .action(cli::Action::StoreTrue)
        .help(cli::kSuppress);
    optional.add({"--plotData"})
        .metavar("FILE")
        .output({"json"})
        .cpp_only("The data of the figure as JSON (and its matrices as FILE.<n>.npy), without "
                  "drawing it (cpp/PLAN.md tier 7).")
        .help("Write the data the figure is drawn from as JSON to this file and do not draw it.");

    const cli::Namespace ns = parser.parse(argc, argv);
    if (const int refused = hicx::plot::preflight("hicPlotMatrix", !ns.given("plotData")); refused != 0) {
        return refused;
    }
    // writableFile's open(string, 'w').close(), now that drawing was accepted.
    {
        const std::string out_file = ns.str("outFileName");
        std::FILE* handle = std::fopen(out_file.c_str(), "w");
        if (handle == nullptr) {
            // The argparse error the check above reports while parsing.
            std::fputs(kUsage, stderr);
            std::fprintf(stderr,
                         "hicPlotMatrix: error: argument --outFileName/-out: %s file can be "
                         "opened for writing\n",
                         out_file.c_str());
            return 2;
        }
        std::fclose(handle);
    }
    const std::string matrix_path = ns.str("matrix");
    const std::optional<std::string> matrix2_path = ns.opt_str("matrix2");
    std::optional<std::string> region = ns.opt_str("region");
    std::optional<std::string> region2 = ns.opt_str("region2");
    // Untouched copies for --matrix2's own load: the --matrix load below
    // mutates region/region2 (it resets them once it has read --chromosomeOrder
    // or a --region into a Selection), and --matrix2 must be read against the
    // same request the user made, not against what is left of it afterwards.
    const std::optional<std::string> requested_region = region;
    const std::optional<std::string> requested_region2 = region2;
    std::optional<std::vector<std::string>> chromosome_order;
    if (ns.given("chromosomeOrder")) {
        chromosome_order = ns.strs("chromosomeOrder");
    }
    const bool per_chromosome = ns.flag("perChromosome");
    const bool log1p = ns.flag("log1p");
    const bool log = ns.flag("log");
    const std::optional<std::string> plot_data = ns.opt_str("plotData");
    std::vector<std::string> temporary_files;

    try {
        if (per_chromosome && region.has_value() && !region->empty()) {
            throw ExitError("ERROR, choose from the option --perChromosome or --region, the two "
                            "options at the same time are not compatible.");
        }
        if (matrix2_path.has_value() && per_chromosome) {
            throw ExitError("ERROR, --matrix2 and --perChromosome are not compatible: --matrix2 "
                            "draws a single upper/lower triangle heatmap, not one heatmap per "
                            "chromosome.");
        }
        if (ns.given("bigwig") && ns.strs("bigwig").size() > 1 &&
            ns.flag("bigwigAdditionalVerticalAxis")) {
            std::fputs("Either multiple bigwig files on x axis or additional vertical axis are "
                       "supported\n",
                       stderr);
        }
        const bool is_cooler = hicx::check_cooler(matrix_path);
        // v4's own addition: a .hic source is block indexed the same way a
        // cooler is (hic_adapter.hpp read_hic), so the single-region fast path
        // below and getRegion's cool style bin overlap test apply to it too,
        // not only to check_cooler's cool/mcool. hicx::CoolFile itself is
        // never opened in this file, so nothing here needs the file to
        // actually be a cooler.
        const bool is_partial_source = is_cooler || hicx::is_hic_path(matrix_path);
        const bool open_cooler_chromosome_order =
            !(chromosome_order.has_value() && chromosome_order->size() > 1);
        bool chromosome_order_as_bytes = false;

        Model model;
        std::optional<std::string> xlabel;
        std::optional<std::string> ylabel;
        std::optional<Selection> selection;
        std::vector<double> matrix_values;
        std::int64_t matrix_rows = 0;
        std::int64_t matrix_cols = 0;

        auto load = [&](const std::optional<std::string>& chrname) {
            hicx::ToolMatrix loaded = hicx::ToolMatrix::load(matrix_path, chrname);
            model.intervals = loaded.cut_intervals();
            model.nan_bins = loaded.nan_bins();
            model.matrix = std::move(loaded.matrix());
            model.view.resize(model.intervals.size());
            for (std::size_t i = 0; i < model.view.size(); ++i) {
                model.view[i] = static_cast<std::int64_t>(i);
            }
        };
        auto clear_masked = [&]() {
            if (ns.flag("clearMaskedBins")) {
                mask_nan_bins(model);
                enlarge_bins(model.intervals);
            }
        };

        if (is_partial_source && !(region2.has_value() && !region2->empty()) && open_cooler_chromosome_order) {
            std::optional<std::string> retrieve;
            if (region.has_value() && !region->empty()) {
                retrieve = *region;
            }
            if (chromosome_order.has_value()) {
                region.reset();
                region2.reset();
                retrieve = chromosome_order->front();
            }
            load(retrieve);
            clear_masked();
            if (region.has_value() && !region->empty()) {
                selection = get_region(*region, region2, model, is_partial_source);
                xlabel = selection->region.chrom;
                ylabel = selection->region2.has_value() ? selection->region2->chrom
                                                        : selection->region.chrom;
            }
            if (!per_chromosome) {
                const std::vector<std::int64_t> bins = all_bins(model);
                matrix_values = dense(model, bins, bins);
                matrix_rows = matrix_cols = static_cast<std::int64_t>(bins.size());
            }
        } else {
            load(std::nullopt);
            clear_masked();
            if (chromosome_order.has_value()) {
                region.reset();
                region2.reset();
                std::vector<std::string> valid;
                std::vector<std::string> invalid;
                for (const std::string& chrom : *chromosome_order) {
                    if (bin_range(model, chrom).has_value()) {
                        valid.push_back(chrom);
                    } else {
                        invalid.push_back(chrom);
                    }
                }
                if (!invalid.empty()) {
                    std::fputs("WARNING: The following chromosome/scaffold names were not found. "
                               "Please checkthe correct spelling of the chromosome names. \n",
                               stderr);
                }
                reorder_chromosomes(model, valid);
                chromosome_order_as_bytes = true;
            }
            if (!stores_any_value(model)) {
                throw ExitError("ValueError: zero-size array to reduction operation minimum which "
                                "has no identity");
            }
            if (region.has_value() && !region->empty()) {
                selection = get_region(*region, region2, model, is_partial_source);
                xlabel = selection->region.chrom;
                ylabel = selection->region2.has_value() ? selection->region2->chrom
                                                        : selection->region.chrom;
                if (!per_chromosome) {
                    matrix_values = dense(model, selection->idx1, selection->idx2);
                    matrix_rows = static_cast<std::int64_t>(selection->idx1.size());
                    matrix_cols = static_cast<std::int64_t>(selection->idx2.size());
                }
            } else if (!per_chromosome) {
                // getMatrix(): the NaN bins as NaN rows and columns.
                const std::vector<std::int64_t> bins = all_bins(model);
                matrix_values = dense(model, bins, bins);
                matrix_rows = matrix_cols = static_cast<std::int64_t>(bins.size());
                for (const std::int64_t bin : model.nan_bins) {
                    for (std::int64_t k = 0; k < matrix_rows; ++k) {
                        matrix_values[static_cast<std::size_t>(bin * matrix_cols + k)] =
                            std::numeric_limits<double>::quiet_NaN();
                        matrix_values[static_cast<std::size_t>(k * matrix_cols + bin)] =
                            std::numeric_limits<double>::quiet_NaN();
                    }
                }
            }
        }

        if (matrix2_path.has_value()) {
            // --perChromosome was already refused together with --matrix2 above,
            // so matrix_values/matrix_rows/matrix_cols were filled by exactly one
            // of the two dense branches above. --matrix2 is read for the same
            // originally requested region/region2/chromosomeOrder, through the
            // identical region-scoped or whole-matrix path load_dense_matrix
            // factors out of the block above, never a separate whole load.
            DenseLoad second;
            try {
                second = load_dense_matrix(*matrix2_path, requested_region, requested_region2,
                                           chromosome_order, ns.flag("clearMaskedBins"));
            } catch (const ExitError& error) {
                throw ExitError(std::string("--matrix2: ") + error.what(), error.code());
            } catch (const std::exception& error) {
                throw ExitError(std::string("--matrix2: ") + error.what());
            }
            if (second.rows != matrix_rows || second.cols != matrix_cols) {
                throw ExitError(
                    "ERROR: --matrix2 does not have the same shape as --matrix for the "
                    "requested region (--matrix: " +
                    std::to_string(matrix_rows) + "x" + std::to_string(matrix_cols) +
                    ", --matrix2: " + std::to_string(second.rows) + "x" +
                    std::to_string(second.cols) + ").");
            }
            combine_triangles(matrix_values, second.values, matrix_rows, matrix_cols);
        }

        const std::int64_t resolution = hicx::BinTable(model.intervals).bin_size();

        plot::JsonObject options;
        auto optional_number = [&](const char* dest) {
            const std::optional<double> value = ns.opt_real(dest);
            return value.has_value() ? plot::json_number(*value) : std::string("null");
        };
        auto optional_string = [&](const char* dest) {
            const std::optional<std::string> value = ns.opt_str(dest);
            return value.has_value() ? plot::json_string(*value) : std::string("null");
        };
        auto number_or_default = [&](const char* dest) {
            return ns.given(dest) ? plot::json_number(ns.real(dest))
                                  : ns.default_of(dest).to_python_string();
        };
        options.add("outFileName", plot::json_string(ns.str("outFileName")));
        options.add("dpi", plot::json_int(ns.integer("dpi")));
        options.add("title", optional_string("title"));
        options.add("scoreName", optional_string("scoreName"));
        options.add("perChromosome", plot::json_bool(per_chromosome));
        options.add("log", plot::json_bool(log));
        options.add("log1p", plot::json_bool(log1p));
        options.add("vMin", optional_number("vMin"));
        options.add("vMax", optional_number("vMax"));
        options.add("colorMap", plot::json_string(ns.str("colorMap")));
        options.add("fontsize", number_or_default("fontsize"));
        options.add("rotationX", number_or_default("rotationX"));
        options.add("rotationY", number_or_default("rotationY"));
        options.add("bigwig", ns.given("bigwig") ? plot::json_strings(ns.strs("bigwig")) : "null");
        options.add("bigwigAdditionalVerticalAxis", plot::json_bool(ns.flag("bigwigAdditionalVerticalAxis")));
        options.add("vMinBigwig", optional_number("vMinBigwig"));
        options.add("vMaxBigwig", optional_number("vMaxBigwig"));
        options.add("flipBigwigSign", plot::json_bool(ns.flag("flipBigwigSign")));
        options.add("scaleFactorBigwig", plot::json_number(ns.real("scaleFactorBigwig")));
        options.add("increaseFigureWidth", plot::json_number(ns.real("increaseFigureWidth")));
        options.add("increaseFigureHeight", plot::json_number(ns.real("increaseFigureHeight")));
        options.add("loops", optional_string("loops"));
        options.add("loopLargeRegionsOperation", plot::json_string(ns.str("loopLargeRegionsOperation")));
        options.add("tads", optional_string("tads"));
        options.add("disable_tight_layout", plot::json_bool(ns.flag("disable_tight_layout")));
        options.add("chromosomeOrder",
                    chromosome_order.has_value() ? plot::json_strings(*chromosome_order) : "null");
        options.add("region", "null");
        options.add("region2", region2.has_value() ? plot::json_string(*region2) : "null");

        auto matrix_file = [&](std::size_t number) {
            if (plot_data.has_value()) {
                return *plot_data + "." + std::to_string(number) + ".npy";
            }
            std::string path = plot::temporary_file();
            temporary_files.push_back(path);
            return path;
        };

        std::vector<std::string> names;
        for (const auto& entry : boundaries(model)) {
            names.push_back(entry.first);
        }

        plot::JsonObject data;
        data.add("options", options.str());
        data.add("resolution", plot::json_int(resolution));
        data.add("chrom_names", plot::json_strings(names));
        data.add("chromosome_start_end", chromosome_start_end_json(model));
        data.add("chromosome_order_as_bytes", plot::json_bool(chromosome_order_as_bytes));

        if (per_chromosome) {
            std::vector<std::string> entries;
            std::size_t number = 0;
            for (const auto& [name, range] : boundaries(model)) {
                std::vector<std::int64_t> bins;
                for (std::int64_t bin = range.first; bin < range.last; ++bin) {
                    bins.push_back(bin);
                }
                std::vector<double> values = dense(model, bins, bins);
                if (log || log1p) {
                    mask_per_chromosome(values);
                }
                if (log1p) {
                    for (double& value : values) {
                        value += 1.0;
                    }
                }
                const Selection chromosome = get_region(name, region2, model, is_partial_source);
                const std::string path = matrix_file(number++);
                plot::write_npy_float64(path, values, static_cast<std::int64_t>(bins.size()),
                                        static_cast<std::int64_t>(bins.size()));
                plot::JsonObject entry;
                entry.add("name", plot::json_string(name));
                entry.add("matrix", plot::json_string(path));
                entry.add("start_pos", plot::json_ints(chromosome.start_pos1));
                entry.add("start_pos2", plot::json_ints(chromosome.start_pos2));
                entry.add("region", region_json(chromosome.region));
                entries.push_back(entry.str());
            }
            data.add("chromosomes", plot::json_list(entries));
        } else {
            if (log || log1p) {
                mask_whole(matrix_values);
            }
            if (log1p) {
                for (double& value : matrix_values) {
                    value += 1.0;
                }
            }
            const std::string path = matrix_file(0);
            plot::write_npy_float64(path, matrix_values, matrix_rows, matrix_cols);
            std::vector<double>().swap(matrix_values);
            data.add("matrix", plot::json_string(path));
            if (selection.has_value()) {
                data.add("start_pos", plot::json_ints(selection->start_pos1));
                data.add("start_pos2", plot::json_ints(selection->start_pos2));
                data.add("region", region_json(selection->region));
            } else {
                data.add("start_pos", plot::json_ints(make_start_pos_array(model)));
                data.add("start_pos2", "null");
                data.add("region", "null");
            }
            data.add("xlabel", xlabel.has_value() ? plot::json_string(*xlabel) : "null");
            data.add("ylabel", ylabel.has_value() ? plot::json_string(*ylabel) : "null");
        }
        data.add("temporary_files", plot::json_strings(temporary_files));
        model = Model();

        hicx::report_resource_usage("hicPlotMatrix");
        const int status = plot::draw("hicPlotMatrix", data.str(), plot_data);
        for (const std::string& path : temporary_files) {
            std::remove(path.c_str());  // reached only when the drawing could not start
        }
        return status;
    } catch (const ExitError& error) {
        std::fprintf(stderr, "%s\n", error.what());
        for (const std::string& path : temporary_files) {
            std::remove(path.c_str());
        }
        return error.code();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicPlotMatrix: %s\n", error.what());
        for (const std::string& path : temporary_files) {
            std::remove(path.c_str());
        }
        return 1;
    }
}
