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
// The tracks are read here too: plotBigwig's scores per bigwig file (through
// libBigWig's bwStats, which pyBigWig's stats calls; the chromosome name
// fallbacks, np.linspace positions, the sign flip, scaling, clip and axis
// limits) for every heatmap, and the positions plotLongRangeContacts and
// plotTADs read from --loops and --tads. plot/hicexplorer_plot/
// hicPlotMatrix.py draws the figure with the reference's plotHeatmap and
// drawing calls from that data.
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
// Memory: pcolormesh needs a value per cell, which for a whole genome is far
// more than the figure has pixels. The matrices are therefore never built
// whole: they are streamed from the sparse matrix as tiles (see "The dense
// matrices of the figure, as tiles" below), aggregated in blocks where an axis
// has more bins than pixels, and written to .npy files as int32 (counts) or
// float32; the sparse matrix is freed before the drawing process starts.

#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "hicx/bigwig_track.hpp"

#include "hicx/argparse.hpp"
#include "hicx/bins.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/hic_adapter.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/plot_tiles.hpp"
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
    "incompatible with --perChromosome.\n"
    "\n"
    "A matrix with more bins along an axis than the figure has pixels there is drawn\n"
    "as blocks of ceil(bins / pixels) bins (the sum of the counts, or the mean of the\n"
    "other values; blocks stop at chromosome ends; --log and --log1p apply to the\n"
    "blocks). The C++-only option --noDownsample draws single bins instead.\n";

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

std::vector<std::int64_t> all_bins(const Model& model) {
    std::vector<std::int64_t> bins(model.intervals.size());
    for (std::size_t i = 0; i < bins.size(); ++i) {
        bins[i] = static_cast<std::int64_t>(i);
    }
    return bins;
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

// ---------------------------------------------------------------------------
// The tracks plotHeatmap draws over the heatmap: plotBigwig's scores and
// plotLongRangeContacts' and plotTADs' positions, read here.

struct ChromExtent {
    std::string name;
    std::int64_t start = 0;
    std::int64_t end = 0;
};

// chromosome_start_end(ma), in its order.
std::vector<ChromExtent> chromosome_extents(const Model& model) {
    std::vector<ChromExtent> extents;
    for (const auto& [name, range] : boundaries(model)) {
        extents.push_back({name, model.intervals[static_cast<std::size_t>(range.first)].start,
                           model.intervals[static_cast<std::size_t>(range.last - 1)].end});
    }
    return extents;
}

const ChromExtent& extent_of(const std::vector<ChromExtent>& extents, const std::string& name) {
    for (const ChromExtent& extent : extents) {
        if (extent.name == name) {
            return extent;
        }
    }
    throw std::runtime_error("KeyError: '" + name + "'");
}

// np.linspace(start, stop, num) as numpy 1.26 computes it.
std::vector<double> linspace(double start, double stop, std::int64_t num) {
    if (num < 0) {
        throw std::runtime_error("ValueError: Number of samples, " + std::to_string(num) +
                                 ", must be non-negative.");
    }
    std::vector<double> y(static_cast<std::size_t>(num));
    const std::int64_t div = num - 1;
    const double delta = stop - start;
    for (std::int64_t i = 0; i < num; ++i) {
        double value = static_cast<double>(i);
        if (div > 0) {
            const double step = delta / static_cast<double>(div);
            if (step == 0.0) {
                value = value / static_cast<double>(div);
                value = value * delta;
            } else {
                value = value * step;
            }
        } else {
            value = value * delta;
        }
        y[static_cast<std::size_t>(i)] = value + start;
    }
    if (num > 1) {
        y.back() = stop;
    }
    return y;
}

std::int64_t floor_div(std::int64_t a, std::int64_t b) {
    std::int64_t q = a / b;
    if ((a % b != 0) && ((a < 0) != (b < 0))) {
        --q;
    }
    return q;
}

// A --bigwig track read as pyBigWig reads it (hicx::bigwig::Track in
// cpp/core, shared with the Python API's plotBigwig), with this tool's
// messages.
struct BigwigFile {
    std::optional<hicx::bigwig::Track> track;
    explicit BigwigFile(const std::string& path) {
        try {
            track.emplace(path);
        } catch (const std::runtime_error&) {
            throw std::runtime_error("RuntimeError: Received an error during file opening! (" +
                                     path + ")");
        }
    }

    // The chromosome's length, or -1 when the file has no such chromosome.
    [[nodiscard]] std::int64_t length(const std::string& chrom) const {
        return track->length(chrom).value_or(-1);
    }

    // bw.stats(chrom, start, end, nBins=bins) of pyBigWig: the mean per bin
    // (NaN, None there, for a bin without data, which astype(float) makes NaN).
    [[nodiscard]] std::vector<double> stats(const std::string& chrom, std::int64_t start,
                                            std::int64_t end, std::int64_t bins) const {
        if (bins < 0) {
            throw std::runtime_error("RuntimeError: An error was encountered while fetching the "
                                     "statistics.");
        }
        try {
            return track->mean(chrom, start, end, bins);
        } catch (const std::invalid_argument&) {
            throw std::runtime_error("RuntimeError: Invalid interval bounds!");
        } catch (const std::runtime_error&) {
            throw std::runtime_error("RuntimeError: An error was encountered while fetching the "
                                     "statistics.");
        }
    }
};

struct BigwigTrackOptions {
    bool flip = false;
    double scale = 1.0;
    std::optional<double> v_min;
    std::optional<double> v_max;
    std::int64_t resolution = 1;
};

// One plotBigwig call: per file in `files` order, {"missing": false, "x",
// "y", "lim": [low, high]}, up to and including the first file whose
// chromosome is missing ({"missing": true}), where the reference returns.
// The scores are flipped, scaled and clipped as plotBigwig does; the axis
// limits are those it sets.
std::string bigwig_tracks_json(const std::vector<std::string>& files,
                               const std::optional<Region>& region,
                               const std::vector<ChromExtent>& extents,
                               const BigwigTrackOptions& options) {
    std::vector<std::string> tracks;
    for (const std::string& path : files) {
        const BigwigFile bw(path);
        std::vector<double> x_values;
        std::vector<double> scores;
        std::string low = "null";
        std::string high = "null";
        bool missing = false;
        if (region.has_value()) {
            std::string chrom = region->chrom;
            const std::int64_t region_start = region->start;
            const std::int64_t region_end = std::min(region->end, extent_of(extents, chrom).end);
            if (bw.length(chrom) < 0) {
                chrom = chrom.rfind("chr", 0) == 0 ? chrom.substr(3) : "chr" + chrom;
                missing = bw.length(chrom) < 0;
            }
            if (!missing) {
                const std::int64_t bigwig_end = std::min(bw.length(chrom), region_end);
                const std::int64_t num_bins =
                    floor_div(bigwig_end - region_start, options.resolution);
                scores = bw.stats(chrom, region_start, bigwig_end, num_bins);
                x_values = linspace(static_cast<double>(region_start),
                                    static_cast<double>(region_end), num_bins);
                if (x_values.size() != scores.size()) {
                    throw std::runtime_error("AssertionError: the bigwig scores of " + path +
                                             " do not match the bins of the region");
                }
                low = plot::json_int(region_start);
                high = plot::json_int(region_end);
            }
        } else {
            std::int64_t chrom_length_sum = 0;
            std::optional<std::int64_t> min_start;
            for (const ChromExtent& extent : extents) {
                std::string chrom = extent.name;
                if (bw.length(chrom) < 0) {
                    chrom = "chr" + chrom;
                    if (bw.length(chrom) < 0) {
                        missing = true;
                        break;
                    }
                }
                const std::int64_t bigwig_end = std::min(bw.length(chrom), extent.end);
                const std::int64_t num_bins =
                    floor_div(bigwig_end - extent.start, options.resolution);
                const std::vector<double> chrom_scores =
                    bw.stats(chrom, extent.start, bigwig_end, num_bins);
                const std::vector<double> chrom_x = linspace(
                    static_cast<double>(chrom_length_sum + extent.start),
                    static_cast<double>(chrom_length_sum + extent.end - extent.start), num_bins);
                if (chrom_x.size() != chrom_scores.size()) {
                    throw std::runtime_error("AssertionError: the bigwig scores of " + path +
                                             " do not match the bins of " + extent.name);
                }
                x_values.insert(x_values.end(), chrom_x.begin(), chrom_x.end());
                scores.insert(scores.end(), chrom_scores.begin(), chrom_scores.end());
                chrom_length_sum += extent.end - extent.start;
                if (!min_start.has_value() || *min_start > extent.start) {
                    min_start = extent.start;
                }
            }
            if (min_start.has_value()) {
                low = plot::json_int(*min_start);
            }
            high = plot::json_int(chrom_length_sum);
        }
        plot::JsonObject track;
        if (missing) {
            track.add("missing", "true");
            tracks.push_back(track.str());
            break;
        }
        for (double& value : scores) {
            if (options.flip) {
                value *= -1;
            }
            if (options.scale != 1.0) {
                value *= options.scale;
            }
            if (options.v_min.has_value() || options.v_max.has_value()) {
                // np.clip: NaN stays NaN, the upper bound applied last.
                if (!std::isnan(value)) {
                    if (options.v_min.has_value() && value < *options.v_min) {
                        value = *options.v_min;
                    }
                    if (options.v_max.has_value() && value > *options.v_max) {
                        value = *options.v_max;
                    }
                }
            }
        }
        track.add("missing", "false");
        track.add("x", plot::json_numbers(x_values));
        track.add("y", plot::json_numbers(scores));
        track.add("lim", plot::json_list({low, high}));
        tracks.push_back(track.str());
    }
    return plot::json_list(tracks);
}

// The lines of a loops or TADs file as plotLongRangeContacts and plotTADs
// read them: readlines() of the bytes, decoded as ASCII (any other byte
// raises there), str.strip(), split on tabs.
std::vector<std::vector<std::string>> track_file_fields(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("FileNotFoundError: [Errno 2] No such file or directory: '" +
                                 path + "'");
    }
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::vector<std::vector<std::string>> lines;
    std::size_t begin = 0;
    while (begin < text.size()) {
        std::size_t end = text.find('\n', begin);
        end = end == std::string::npos ? text.size() : end + 1;
        std::string line = text.substr(begin, end - begin);
        begin = end;
        for (const char c : line) {
            if (static_cast<unsigned char>(c) >= 0x80) {
                throw std::runtime_error("UnicodeDecodeError: 'ascii' codec can't decode a byte in " +
                                         path);
            }
        }
        auto space = [](char c) {
            return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\x0b' || c == '\x0c' ||
                   (c >= '\x1c' && c <= '\x1f');
        };
        std::size_t first = 0;
        std::size_t last = line.size();
        while (first < last && space(line[first])) {
            ++first;
        }
        while (last > first && space(line[last - 1])) {
            --last;
        }
        line = line.substr(first, last - first);
        std::vector<std::string> fields;
        std::size_t from = 0;
        for (std::size_t tab = line.find('\t'); tab != std::string::npos;
             tab = line.find('\t', from)) {
            fields.push_back(line.substr(from, tab - from));
            from = tab + 1;
        }
        fields.push_back(line.substr(from));
        lines.push_back(std::move(fields));
    }
    return lines;
}

// `chrom not in pChromosomeOrder`: always true when the order became bytes.
bool not_in_order(const std::string& chrom,
                  const std::optional<std::vector<std::string>>& chromosome_order,
                  bool order_as_bytes) {
    if (order_as_bytes) {
        return true;
    }
    return std::find(chromosome_order->begin(), chromosome_order->end(), chrom) ==
           chromosome_order->end();
}

// plotLongRangeContacts' x_list and y_list (pRegion None raises inside the
// try block on every line, so nothing is kept without --region).
std::string loop_positions_json(const std::string& path, const std::optional<Region>& region,
                                const std::optional<std::vector<std::string>>& chromosome_order,
                                bool order_as_bytes, const std::string& operation) {
    std::vector<std::int64_t> xs;
    std::vector<std::int64_t> ys;
    for (const std::vector<std::string>& fields : track_file_fields(path)) {
        if (fields.size() < 6) {
            continue;  // the unpacking of fields[0:3] or fields[3:6] raises
        }
        const std::string& chrom_x = fields[0];
        const std::string& chrom_y = fields[3];
        if (region.has_value() && (chrom_x != region->chrom || chrom_y != region->chrom)) {
            continue;
        }
        if (chromosome_order.has_value() &&
            (not_in_order(chrom_x, chromosome_order, order_as_bytes) ||
             not_in_order(chrom_y, chromosome_order, order_as_bytes))) {
            continue;
        }
        std::int64_t start_x = 0, end_x = 0, start_y = 0, end_y = 0, x = 0, y = 0;
        if (operation == "first") {
            if (!hicx::cli::python_int(fields[1], &start_x) ||
                !hicx::cli::python_int(fields[4], &start_y)) {
                continue;
            }
            x = start_x;
            y = start_y;
        } else if (operation == "last") {
            if (!hicx::cli::python_int(fields[2], &end_x) ||
                !hicx::cli::python_int(fields[5], &end_y)) {
                continue;
            }
            x = end_x;
            y = end_y;
        } else {
            if (!hicx::cli::python_int(fields[1], &start_x) ||
                !hicx::cli::python_int(fields[2], &end_x) ||
                !hicx::cli::python_int(fields[4], &start_y) ||
                !hicx::cli::python_int(fields[5], &end_y)) {
                continue;
            }
            x = floor_div(start_x + end_x, 2);
            y = floor_div(start_y + end_y, 2);
        }
        if (!region.has_value()) {
            continue;  // int(pRegion[1]) raises TypeError
        }
        if (x >= region->start && x <= region->end && y >= region->start && y <= region->end) {
            xs.push_back(x);
            ys.push_back(y);
        }
    }
    plot::JsonObject out;
    out.add("x", plot::json_ints(xs));
    out.add("y", plot::json_ints(ys));
    return out.str();
}

// plotTADs' x_list and y_list: the start and end of every domain inside the
// region.
std::string tad_positions_json(const std::string& path, const std::optional<Region>& region,
                               const std::optional<std::vector<std::string>>& chromosome_order,
                               bool order_as_bytes) {
    std::vector<std::int64_t> xs;
    std::vector<std::int64_t> ys;
    for (const std::vector<std::string>& fields : track_file_fields(path)) {
        if (fields.size() < 3) {
            continue;
        }
        const std::string& chrom_x = fields[0];
        if (region.has_value() && chrom_x != region->chrom) {
            continue;
        }
        if (chromosome_order.has_value() &&
            not_in_order(chrom_x, chromosome_order, order_as_bytes)) {
            continue;
        }
        std::int64_t x = 0;
        std::int64_t y = 0;
        if (!hicx::cli::python_int(fields[1], &x) || !hicx::cli::python_int(fields[2], &y)) {
            continue;
        }
        if (!region.has_value()) {
            continue;
        }
        if (x >= region->start && x <= region->end && y >= region->start && y <= region->end) {
            xs.push_back(x);
            ys.push_back(y);
        }
    }
    plot::JsonObject out;
    out.add("x", plot::json_ints(xs));
    out.add("y", plot::json_ints(ys));
    return out.str();
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

// ---------------------------------------------------------------------------
// The dense matrices of the figure, as tiles.
//
// pcolormesh needs the values of every cell, so the figure's matrix has one
// value per bin pair; for a whole genome that is far more than the figure has
// pixels, and a dense array of it does not fit in memory. The matrix is
// therefore never built whole. It is streamed from the sparse matrix tile by
// tile (a rectangle of at most kDefaultTileBins x kDefaultTileBins cells), the
// statistics of all tiles are gathered in a first pass, and the tiles are
// written as .npy files in a second one; the drawing layer draws one tile at a
// time on the same axes with a shared colour scale.
//
// When an axis has more bins than the figure has pixels along it, blocks of
// k bins are aggregated before the handoff instead of single bins
// (k = ceil(bins / pixels)): the sum of the block for count matrices, the mean
// of its valid cells otherwise. The blocks never cross a chromosome boundary
// (the last block of a chromosome is shorter) unless that would give more than
// twice as many blocks as pixels. --log and --log1p are applied to the
// aggregated values, which is the plot of the coarser matrix: the visual
// equivalent of the full resolution plot, whose bins are smaller than a pixel.
// --noDownsample keeps single bins (the tiles then cover the full resolution).
//
// A block is NaN when it holds no valid cell: a cell is invalid when its row or
// column is a NaN bin (the whole-matrix path of getMatrix) or when its stored
// value is NaN. The valid cells of a block that is not empty are the bin pairs
// of non-NaN bins, with the zero of every unstored pair; the mean divides by
// their number.
//
// All sums are in double (exact for integers below 2^53) and become int32,
// float32 or, where those do not hold the values, int64 or float64 at the
// handoff (hicx/plot_tiles.hpp).

constexpr std::int64_t kDefaultTileBins = 4096;

// One matrix and the bins of its rows and columns; rows and columns of the
// figure are positions in `rows` and `cols`.
struct Source {
    const Model* model = nullptr;
    std::vector<std::int64_t> rows;
    std::vector<std::int64_t> cols;
    bool mark_nan = false;  // getMatrix(): the NaN bins are NaN rows and columns
};

struct AxisMap {
    std::vector<std::int64_t> loaded;      // row (column) of the loaded matrix per position
    std::vector<std::int64_t> inv_off;     // positions of each loaded row, CSR style
    std::vector<std::int64_t> inv_pos;
    std::vector<char> nan;                 // the position is a NaN bin
    std::vector<std::int32_t> chromosome;  // run number of the chromosome
};

AxisMap make_axis(const Model& model, const std::vector<std::int64_t>& bins, bool mark_nan) {
    AxisMap axis;
    const std::size_t n = bins.size();
    const auto loaded_size = static_cast<std::size_t>(model.matrix.rows());
    axis.loaded.resize(n);
    axis.nan.assign(n, 0);
    axis.chromosome.resize(n);
    std::vector<char> flagged;
    if (mark_nan) {
        flagged.assign(model.intervals.size(), 0);
        for (const std::int64_t bin : model.nan_bins) {
            flagged[static_cast<std::size_t>(bin)] = 1;
        }
    }
    std::int32_t run = 0;
    for (std::size_t p = 0; p < n; ++p) {
        const auto bin = static_cast<std::size_t>(bins[p]);
        axis.loaded[p] = model.view[bin];
        if (mark_nan) {
            axis.nan[p] = flagged[bin];
        }
        if (p > 0 && model.intervals[bin].chrom !=
                         model.intervals[static_cast<std::size_t>(bins[p - 1])].chrom) {
            ++run;
        }
        axis.chromosome[p] = run;
    }
    axis.inv_off.assign(loaded_size + 1, 0);
    for (const std::int64_t row : axis.loaded) {
        ++axis.inv_off[static_cast<std::size_t>(row) + 1];
    }
    for (std::size_t i = 0; i < loaded_size; ++i) {
        axis.inv_off[i + 1] += axis.inv_off[i];
    }
    axis.inv_pos.resize(n);
    std::vector<std::int64_t> fill(axis.inv_off.begin(), axis.inv_off.end() - 1);
    for (std::size_t p = 0; p < n; ++p) {
        axis.inv_pos[static_cast<std::size_t>(fill[static_cast<std::size_t>(axis.loaded[p])]++)] =
            static_cast<std::int64_t>(p);
    }
    return axis;
}

// Whether the stored column indices of every row are sorted, and whether the
// stored values are all integers (a count matrix).
struct MatrixTraits {
    bool sorted = true;
    bool integral = true;
};

MatrixTraits scan_traits(const hicx::CsrMatrix& matrix) {
    MatrixTraits traits;
    const auto& indptr = matrix.indptr();
    const auto& indices = matrix.indices();
    const auto& data = matrix.data();
    for (std::int64_t row = 0; row < matrix.rows(); ++row) {
        const auto begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const auto end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            if (k > begin && indices[k] < indices[k - 1]) {
                traits.sorted = false;
            }
            const double v = data[k];
            if (std::isfinite(v) && v != std::floor(v)) {
                traits.integral = false;
            }
        }
    }
    return traits;
}

struct Prepared {
    const Model* model = nullptr;
    AxisMap row;
    AxisMap col;
    MatrixTraits traits;
    std::vector<std::int64_t> valid_row;  // valid (non NaN bin) positions per block
    std::vector<std::int64_t> valid_col;
};

Prepared prepare(const Source& source) {
    Prepared p;
    p.model = source.model;
    p.row = make_axis(*source.model, source.rows, source.mark_nan);
    p.col = make_axis(*source.model, source.cols, source.mark_nan);
    p.traits = scan_traits(source.model->matrix);
    return p;
}

std::vector<std::int64_t> valid_per_block(const AxisMap& axis, const hicx::plot::AxisBlocks& blocks) {
    std::vector<std::int64_t> valid(static_cast<std::size_t>(blocks.count()), 0);
    for (std::size_t p = 0; p < axis.nan.size(); ++p) {
        if (!axis.nan[p]) {
            ++valid[static_cast<std::size_t>(blocks.block_of[p])];
        }
    }
    return valid;
}

struct Geometry {
    hicx::plot::AxisBlocks rows;
    hicx::plot::AxisBlocks cols;
    bool mean = false;        // average the valid cells of a block instead of summing them
    bool transposed = false;  // tile[cols][rows], pcolormesh's orientation
};

// Calls emit(i, j, value) for every stored value of the (symmetric) matrix at
// a pair of positions i in [i0, i1) and j in [j0, j1).
template <typename Emit>
void scan_tile(const Prepared& p, std::int64_t i0, std::int64_t i1, std::int64_t j0,
               std::int64_t j1, Emit&& emit) {
    const hicx::CsrMatrix& matrix = p.model->matrix;
    const bool upper = matrix.symmetry() == hicx::Symmetry::UpperTriangle;
    const auto& indptr = matrix.indptr();
    const auto& indices = matrix.indices();
    const auto& data = matrix.data();
    auto span = [](const AxisMap& axis, std::int64_t from, std::int64_t to) {
        std::int64_t low = std::numeric_limits<std::int64_t>::max();
        std::int64_t high = -1;
        for (std::int64_t q = from; q < to; ++q) {
            low = std::min(low, axis.loaded[static_cast<std::size_t>(q)]);
            high = std::max(high, axis.loaded[static_cast<std::size_t>(q)]);
        }
        return std::pair<std::int64_t, std::int64_t>(low, high);
    };
    // Entries of one stored row with a column in [low, high].
    auto each_in_row = [&](std::int64_t row, std::int64_t low, std::int64_t high, auto&& visit) {
        auto begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        auto end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        if (p.traits.sorted) {
            const auto* base = indices.data();
            end = static_cast<std::size_t>(std::upper_bound(base + begin, base + end,
                                                            static_cast<std::int32_t>(high)) -
                                           base);
            begin = static_cast<std::size_t>(std::lower_bound(base + begin, base + end,
                                                              static_cast<std::int32_t>(low)) -
                                             base);
        }
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t column = indices[k];
            if (column >= low && column <= high) {
                visit(column, data[k]);
            }
        }
    };
    const auto [col_low, col_high] = span(p.col, j0, j1);
    for (std::int64_t i = i0; i < i1; ++i) {
        each_in_row(p.row.loaded[static_cast<std::size_t>(i)], col_low, col_high,
                    [&](std::int64_t column, double value) {
                        for (auto q = p.col.inv_off[static_cast<std::size_t>(column)];
                             q < p.col.inv_off[static_cast<std::size_t>(column) + 1]; ++q) {
                            const std::int64_t j = p.col.inv_pos[static_cast<std::size_t>(q)];
                            if (j >= j0 && j < j1) {
                                emit(i, j, value);
                            }
                        }
                    });
    }
    if (upper) {
        // The mirrored entries: stored (c, r) with r > c is the cell (r, c).
        const auto [row_low, row_high] = span(p.row, i0, i1);
        for (std::int64_t j = j0; j < j1; ++j) {
            const std::int64_t stored_row = p.col.loaded[static_cast<std::size_t>(j)];
            each_in_row(stored_row, row_low, row_high, [&](std::int64_t column, double value) {
                if (column == stored_row) {
                    return;
                }
                for (auto q = p.row.inv_off[static_cast<std::size_t>(column)];
                     q < p.row.inv_off[static_cast<std::size_t>(column) + 1]; ++q) {
                    const std::int64_t i = p.row.inv_pos[static_cast<std::size_t>(q)];
                    if (i >= i0 && i < i1) {
                        emit(i, j, value);
                    }
                }
            });
        }
    }
}

// The cells of the blocks [bi0, bi1) x [bj0, bj1) of one source.
std::vector<double> compute_tile(const Prepared& p, const Geometry& g, std::int64_t bi0,
                                 std::int64_t bi1, std::int64_t bj0, std::int64_t bj1) {
    const std::int64_t nbi = bi1 - bi0;
    const std::int64_t nbj = bj1 - bj0;
    std::vector<double> sums(static_cast<std::size_t>(nbi * nbj), 0.0);
    std::vector<std::uint32_t> nan_count;  // allocated on the first stored NaN
    auto index = [&](std::int64_t bi, std::int64_t bj) {
        return static_cast<std::size_t>(g.transposed ? (bj - bj0) * nbi + (bi - bi0)
                                                     : (bi - bi0) * nbj + (bj - bj0));
    };
    scan_tile(p, g.rows.first[static_cast<std::size_t>(bi0)],
              g.rows.first[static_cast<std::size_t>(bi1)],
              g.cols.first[static_cast<std::size_t>(bj0)],
              g.cols.first[static_cast<std::size_t>(bj1)],
              [&](std::int64_t i, std::int64_t j, double value) {
                  if (p.row.nan[static_cast<std::size_t>(i)] || p.col.nan[static_cast<std::size_t>(j)]) {
                      return;
                  }
                  const std::size_t at = index(g.rows.block_of[static_cast<std::size_t>(i)],
                                               g.cols.block_of[static_cast<std::size_t>(j)]);
                  if (std::isnan(value)) {
                      if (nan_count.empty()) {
                          nan_count.assign(sums.size(), 0);
                      }
                      ++nan_count[at];
                  } else {
                      sums[at] += value;
                  }
              });
    for (std::int64_t bi = bi0; bi < bi1; ++bi) {
        for (std::int64_t bj = bj0; bj < bj1; ++bj) {
            const std::size_t at = index(bi, bj);
            const std::int64_t valid = p.valid_row[static_cast<std::size_t>(bi)] *
                                           p.valid_col[static_cast<std::size_t>(bj)] -
                                       (nan_count.empty() ? 0 : static_cast<std::int64_t>(nan_count[at]));
            if (valid <= 0) {
                sums[at] = std::numeric_limits<double>::quiet_NaN();
            } else if (g.mean) {
                sums[at] /= static_cast<double>(valid);
            }
        }
    }
    return sums;
}

struct TileRange {
    std::int64_t r0, r1, c0, c1;
};

// Everything a heatmap of the figure is drawn from.
struct Heatmap {
    Source first;
    const Source* second = nullptr;  // --matrix2: the lower triangle
    std::vector<std::int64_t> start_pos_rows;
    std::vector<std::int64_t> start_pos_cols;
    bool per_chromosome = false;     // the masking of plotPerChr
    double pixels_x = 0.0;           // upper bounds of the axes in pixels
    double pixels_y = 0.0;
};

struct TileSettings {
    bool log = false;
    bool log1p = false;
    bool downsample = true;
    std::int64_t tile_bins = kDefaultTileBins;
    std::function<std::string()> new_file;  // the path of the next .npy file
};

template <typename T>
void put_values(std::vector<char>& buffer, const std::vector<double>& values, double add,
                const hicx::plot::MaskParams& mask) {
    buffer.resize(values.size() * sizeof(T));
    auto* out = reinterpret_cast<T*>(buffer.data());
    for (std::size_t i = 0; i < values.size(); ++i) {
        out[i] = static_cast<T>(mask.apply(values[i]) + add);
    }
}

// Writes the tiles of `heat` and describes them in `target`: "matrix" is the
// path of the one file when the whole matrix is a single tile at full
// resolution (the figure then takes the path of a single pcolormesh call and
// the function returns true), else null and "tile_spec" lists the tiles.
bool write_heatmap(const Heatmap& heat, const TileSettings& settings,
                   hicx::plot::JsonObject& target) {
    namespace hp = hicx::plot;
    const Source& a = heat.first;
    Prepared prepared = prepare(a);
    std::optional<Prepared> prepared2;
    bool integral = prepared.traits.integral;
    if (heat.second != nullptr) {
        prepared2 = prepare(*heat.second);
        integral = integral && prepared2->traits.integral;
    }

    const auto n_rows = static_cast<std::int64_t>(a.rows.size());
    const auto n_cols = static_cast<std::int64_t>(a.cols.size());
    const auto px = static_cast<std::int64_t>(std::ceil(heat.pixels_x));
    const auto py = static_cast<std::int64_t>(std::ceil(heat.pixels_y));
    Geometry g;
    std::int64_t k_rows = hp::positions_per_block(n_rows, px, settings.downsample);
    std::int64_t k_cols = hp::positions_per_block(n_cols, py, settings.downsample);
    if (n_rows == n_cols || prepared2) {
        // A square matrix keeps square blocks, so its diagonal is the diagonal of
        // the blocks (and the triangles of --matrix2 are triangles of blocks).
        k_rows = k_cols = std::max(k_rows, k_cols);
    }
    g.rows = hp::make_blocks(n_rows, k_rows, prepared.row.chromosome, 2 * px);
    g.cols = hp::make_blocks(n_cols, k_cols, prepared.col.chromosome, 2 * py);
    const bool reduced = g.rows.k > 1 || g.cols.k > 1;
    g.mean = reduced && !integral;
    const std::int64_t nbr = g.rows.count();
    const std::int64_t nbc = g.cols.count();
    const std::int64_t tile = std::max<std::int64_t>(settings.tile_bins, 1);
    const bool single_matrix = !reduced && nbr <= tile && nbc <= tile;
    g.transposed = !single_matrix;
    prepared.valid_row = valid_per_block(prepared.row, g.rows);
    prepared.valid_col = valid_per_block(prepared.col, g.cols);
    if (prepared2) {
        prepared2->valid_row = valid_per_block(prepared2->row, g.rows);
        prepared2->valid_col = valid_per_block(prepared2->col, g.cols);
    }

    std::vector<TileRange> tiles;
    for (std::int64_t r = 0; r < nbr; r += tile) {
        for (std::int64_t c = 0; c < nbc; c += tile) {
            tiles.push_back({r, std::min(r + tile, nbr), c, std::min(c + tile, nbc)});
        }
    }
    auto make_tile = [&](const TileRange& t) {
        std::vector<double> values = compute_tile(prepared, g, t.r0, t.r1, t.c0, t.c1);
        if (prepared2) {
            const std::vector<double> lower = compute_tile(*prepared2, g, t.r0, t.r1, t.c0, t.c1);
            const std::int64_t nbi = t.r1 - t.r0;
            const std::int64_t nbj = t.c1 - t.c0;
            for (std::int64_t bi = t.r0; bi < t.r1; ++bi) {
                for (std::int64_t bj = t.c0; bj < t.c1; ++bj) {
                    // the block of column index below the row index: the lower triangle
                    if (bj < bi) {
                        const auto at = static_cast<std::size_t>(
                            g.transposed ? (bj - t.c0) * nbi + (bi - t.r0)
                                         : (bi - t.r0) * nbj + (bj - t.c0));
                        values[at] = lower[at];
                    }
                }
            }
        }
        return values;
    };

    // First pass: the statistics of all cells. One tile is kept for the second.
    hp::MaskStats stats;
    std::vector<double> held;
    for (const TileRange& t : tiles) {
        std::vector<double> values = make_tile(t);
        for (const double v : values) {
            stats.add(v);
        }
        if (tiles.size() == 1) {
            held = std::move(values);
        }
    }
    const bool masked = settings.log || settings.log1p;
    hp::MaskParams mask;
    if (masked) {
        try {
            mask = heat.per_chromosome ? hp::mask_params_per_chromosome(stats)
                                       : hp::mask_params_whole(stats);
        } catch (const std::runtime_error& error) {
            throw ExitError(error.what());
        }
    }
    const double add = settings.log1p ? 1.0 : 0.0;
    hp::MaskParams::Summary summary;
    if (masked) {
        summary = mask.summarize(stats);
    } else {
        summary.nonfinite = stats.nans > 0 || stats.infs > 0;
        summary.integral = stats.all_integral;
        summary.max_abs = stats.max_abs_finite;
    }
    const hp::HandoffChoice choice =
        hp::choose_handoff(summary.integral, summary.nonfinite, summary.max_abs + add);
    if (choice.widened) {
        std::fprintf(stderr,
                     "hicPlotMatrix: warning: the largest value of this matrix (%.17g) does not "
                     "fit %s; it is handed to the drawing layer as %s.\n",
                     summary.max_abs + add,
                     choice.type == hp::Handoff::Int64 ? "int32" : "a float32 next to NaN",
                     choice.type == hp::Handoff::Int64 ? "int64" : "float64");
    }

    // Second pass: the tiles.
    double value_min = std::numeric_limits<double>::infinity();
    double value_max = -std::numeric_limits<double>::infinity();
    std::vector<std::string> tile_items;
    std::vector<char> buffer;
    std::string single_path;
    for (const TileRange& t : tiles) {
        std::vector<double> values = tiles.size() == 1 ? std::move(held) : make_tile(t);
        for (const double v : values) {
            const double x = mask.apply(v) + add;
            if (std::isfinite(x) && (!masked || x > 0.0)) {
                value_min = std::min(value_min, x);
                value_max = std::max(value_max, x);
            }
        }
        switch (choice.type) {
            case hp::Handoff::Int32:
                put_values<std::int32_t>(buffer, values, add, mask);
                break;
            case hp::Handoff::Int64:
                put_values<std::int64_t>(buffer, values, add, mask);
                break;
            case hp::Handoff::Float32:
                put_values<float>(buffer, values, add, mask);
                break;
            case hp::Handoff::Float64:
                put_values<double>(buffer, values, add, mask);
                break;
        }
        const std::string path = settings.new_file();
        const std::int64_t nbi = t.r1 - t.r0;
        const std::int64_t nbj = t.c1 - t.c0;
        plot::write_npy_raw(path, hp::npy_dtype(choice.type),
                            g.transposed ? std::vector<std::int64_t>{nbj, nbi}
                                         : std::vector<std::int64_t>{nbi, nbj},
                            buffer.data(), buffer.size());
        if (single_matrix) {
            single_path = path;
        }
        plot::JsonObject item;
        item.add("file", plot::json_string(path));
        item.add("r0", plot::json_int(t.r0));
        item.add("r1", plot::json_int(t.r1));
        item.add("c0", plot::json_int(t.c0));
        item.add("c1", plot::json_int(t.c1));
        tile_items.push_back(item.str());
    }

    if (single_matrix) {
        target.add("matrix", plot::json_string(single_path));
        return true;
    }
    target.add("matrix", "null");
    plot::JsonObject spec;
    spec.add("x_edges", plot::json_numbers(hp::block_edges(hp::nearest_edges(heat.start_pos_rows), g.rows)));
    spec.add("y_edges", plot::json_numbers(hp::block_edges(hp::nearest_edges(heat.start_pos_cols), g.cols)));
    spec.add("bins", plot::json_ints({n_rows, n_cols}));
    spec.add("blocks", plot::json_ints({nbr, nbc}));
    spec.add("bins_per_block", plot::json_ints({g.rows.k, g.cols.k}));
    spec.add("blocks_by_chromosome", plot::json_bool(g.rows.by_chromosome || g.cols.by_chromosome));
    spec.add("aggregate", plot::json_string(reduced ? (g.mean ? "mean" : "sum") : "none"));
    spec.add("dtype", plot::json_string(hp::npy_dtype(choice.type)));
    const bool bounded = tiles.size() > 1 && value_min <= value_max;
    spec.add("vmin", bounded ? plot::json_number(value_min) : std::string("null"));
    spec.add("vmax", bounded ? plot::json_number(value_max) : std::string("null"));
    spec.add("tiles", plot::json_list(tile_items));
    target.add("tile_spec", spec.str());
    if (reduced) {
        std::fprintf(stderr,
                     "hicPlotMatrix: %lld x %lld bins are drawn as %lld x %lld blocks of %lld x "
                     "%lld bins (%s of the cells of a block%s); --noDownsample draws single bins.\n",
                     static_cast<long long>(n_rows), static_cast<long long>(n_cols),
                     static_cast<long long>(nbr), static_cast<long long>(nbc),
                     static_cast<long long>(g.rows.k), static_cast<long long>(g.cols.k),
                     g.mean ? "mean" : "sum",
                     (g.rows.by_chromosome || g.cols.by_chromosome) ? ", cut at chromosome ends"
                                                                    : "");
    }
    return false;
}

// A second matrix (--matrix2) read for the same request as the first.
struct SecondMatrix {
    Model model;
    Source source;
};

// Loads path for the given request the same way --matrix's own load does
// (is_partial_source fast path first, whole-matrix path otherwise; see the
// top-of-file comment), and says which bins the heatmap rows and columns are.
// region/region2 are taken by value: the fast path may consume them (as
// main()'s --matrix load does), and each caller needs its own request untouched
// by the other's.
void load_second(SecondMatrix& out, const std::string& path, std::optional<std::string> region,
                 std::optional<std::string> region2,
                 const std::optional<std::vector<std::string>>& chromosome_order,
                 bool clear_masked_bins) {
    const bool is_cooler = hicx::check_cooler(path);
    const bool is_partial_source = is_cooler || hicx::is_hic_path(path);
    const bool open_cooler_chromosome_order =
        !(chromosome_order.has_value() && chromosome_order->size() > 1);

    Model& model = out.model;
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

    out.source.model = &model;
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
        out.source.rows = all_bins(model);
        out.source.cols = out.source.rows;
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
            out.source.rows = selection.idx1;
            out.source.cols = selection.idx2;
        } else {
            out.source.rows = all_bins(model);
            out.source.cols = out.source.rows;
            out.source.mark_nan = true;
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
    optional.add({"--noDownsample"})
        .action(cli::Action::StoreTrue)
        .cpp_only("Draw single bins even where the matrix has more bins along an axis than the "
                  "figure has pixels. By default the bins are aggregated in blocks of "
                  "ceil(bins / pixels) bins (sum of the counts, mean of the other values, blocks "
                  "cut at chromosome ends), --log and --log1p are applied to the blocks, and the "
                  "axes are unchanged. The matrix is handed to the drawing layer in tiles either "
                  "way; the drawing process keeps every tile it has drawn until the figure is "
                  "saved (about 55 bytes per cell).")
        .help("Do not aggregate bins into blocks when the matrix has more bins than the figure "
              "has pixels.");
    optional.add({"--tileBins"})
        .type("int")
        .default_value(kDefaultTileBins)
        .cpp_only("Largest edge of a tile of the matrix handed to the drawing layer, in bins or "
                  "blocks (Default: %(default)s). A matrix that fits one tile is drawn with one "
                  "pcolormesh call.")
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
        Source primary;
        primary.model = &model;

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
                primary.rows = all_bins(model);
                primary.cols = primary.rows;
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
                    primary.rows = selection->idx1;
                    primary.cols = selection->idx2;
                    }
            } else if (!per_chromosome) {
                // getMatrix(): the NaN bins as NaN rows and columns.
                primary.rows = all_bins(model);
                primary.cols = primary.rows;
                primary.mark_nan = true;
            }
        }

        SecondMatrix second;
        if (matrix2_path.has_value()) {
            // --perChromosome was already refused together with --matrix2 above,
            // so `primary` was filled by exactly one of the two branches above.
            // --matrix2 is read for the same originally requested
            // region/region2/chromosomeOrder, through the identical
            // region-scoped or whole-matrix path load_second factors out of the
            // block above, never a separate whole load.
            try {
                load_second(second, *matrix2_path, requested_region, requested_region2,
                            chromosome_order, ns.flag("clearMaskedBins"));
            } catch (const ExitError& error) {
                throw ExitError(std::string("--matrix2: ") + error.what(), error.code());
            } catch (const std::exception& error) {
                throw ExitError(std::string("--matrix2: ") + error.what());
            }
            if (second.source.rows.size() != primary.rows.size() ||
                second.source.cols.size() != primary.cols.size()) {
                throw ExitError(
                    "ERROR: --matrix2 does not have the same shape as --matrix for the "
                    "requested region (--matrix: " +
                    std::to_string(primary.rows.size()) + "x" + std::to_string(primary.cols.size()) +
                    ", --matrix2: " + std::to_string(second.source.rows.size()) + "x" +
                    std::to_string(second.source.cols.size()) + ").");
            }
        }
        const auto matrix_rows = static_cast<std::int64_t>(primary.rows.size());

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

        std::size_t file_number = 0;
        auto matrix_file = [&]() {
            const std::size_t number = file_number++;
            if (plot_data.has_value()) {
                return *plot_data + "." + std::to_string(number) + ".npy";
            }
            std::string path = plot::temporary_file();
            temporary_files.push_back(path);
            return path;
        };
        TileSettings tile_settings;
        tile_settings.log = log;
        tile_settings.log1p = log1p;
        tile_settings.downsample = !ns.flag("noDownsample");
        tile_settings.tile_bins = ns.integer("tileBins");
        tile_settings.new_file = matrix_file;
        if (tile_settings.tile_bins < 1) {
            throw ExitError("ERROR: --tileBins must be at least 1.");
        }
        // The axes of the figure in pixels, at most: the heatmap axes without
        // bigwig tracks are 5.0 x 4.8 inches, the figure is 8 x 7 inches plus
        // the increase per track, and per chromosome the cell of the grid
        // (plotPerChr) bounds each axes.
        const double dpi = static_cast<double>(ns.integer("dpi"));
        const std::size_t bigwig_count = ns.given("bigwig") ? ns.strs("bigwig").size() : 0;
        const double figure_width = 8.0 + static_cast<double>(bigwig_count) * ns.real("increaseFigureWidth");
        const double figure_height = 7.0 + static_cast<double>(bigwig_count) * ns.real("increaseFigureHeight");
        const std::size_t chromosome_count = boundaries(model).size();
        const auto grid_rows = static_cast<double>((chromosome_count + 4) / 5);
        const auto grid_cols = static_cast<double>(std::min<std::size_t>(chromosome_count, 5));

        std::vector<std::string> names;
        for (const auto& entry : boundaries(model)) {
            names.push_back(entry.first);
        }

        // The tracks of every heatmap plotHeatmap draws (a matrix of fewer
        // than five rows is skipped there, tracks included).
        const std::vector<ChromExtent> extents = chromosome_extents(model);
        const bool bigwig = ns.given("bigwig");
        const std::vector<std::string> bigwig_files =
            bigwig ? ns.strs("bigwig") : std::vector<std::string>{};
        const bool vertical = ns.flag("bigwigAdditionalVerticalAxis");
        BigwigTrackOptions track_options;
        track_options.flip = ns.flag("flipBigwigSign");
        track_options.scale = ns.real("scaleFactorBigwig");
        track_options.v_min = ns.opt_real("vMinBigwig");
        track_options.v_max = ns.opt_real("vMaxBigwig");
        track_options.resolution = resolution;
        bool bigwig_checked = false;
        // plotBigwig's horizontal and vertical tracks for one heatmap.
        auto bigwig_json = [&](const std::optional<Region>& track_region) -> std::string {
            if (!bigwig) {
                return "null";
            }
            if (!bigwig_checked) {
                for (const std::string& file : bigwig_files) {
                    const std::string extension = file.substr(file.rfind('.') + 1);
                    if (extension != "bigwig" && extension != "bw") {
                        throw ExitError("ERROR:hicexplorer.hicPlotMatrix:Given files are not bigwig");
                    }
                }
                // hicx::bigwig::Track initialises libBigWig on first use
                bigwig_checked = true;
            }
            plot::JsonObject tracks;
            tracks.add("horizontal",
                       bigwig_tracks_json(bigwig_files, track_region, extents, track_options));
            if (vertical) {
                const std::vector<std::string> reversed(bigwig_files.rbegin(), bigwig_files.rend());
                tracks.add("vertical",
                           bigwig_tracks_json(reversed, track_region, extents, track_options));
            } else {
                tracks.add("vertical", "null");
            }
            return tracks.str();
        };

        plot::JsonObject data;
        data.add("options", options.str());
        data.add("resolution", plot::json_int(resolution));
        data.add("chrom_names", plot::json_strings(names));
        data.add("chromosome_start_end", chromosome_start_end_json(model));
        data.add("chromosome_order_as_bytes", plot::json_bool(chromosome_order_as_bytes));

        if (per_chromosome) {
            std::vector<std::string> entries;
            for (const auto& [name, range] : boundaries(model)) {
                std::vector<std::int64_t> bins;
                for (std::int64_t bin = range.first; bin < range.last; ++bin) {
                    bins.push_back(bin);
                }
                const Selection chromosome = get_region(name, region2, model, is_partial_source);
                Heatmap heat;
                heat.first.model = &model;
                heat.first.rows = bins;
                heat.first.cols = bins;
                heat.start_pos_rows = chromosome.start_pos1;
                heat.start_pos_cols = chromosome.start_pos2;
                heat.per_chromosome = true;
                heat.pixels_x = (6.0 * (grid_cols * 1.05 + 0.10) +
                                 static_cast<double>(bigwig_count) * ns.real("increaseFigureWidth")) /
                                grid_cols * dpi;
                heat.pixels_y = (6.0 * grid_rows +
                                 static_cast<double>(bigwig_count) * ns.real("increaseFigureHeight")) /
                                grid_rows * dpi;
                plot::JsonObject entry;
                entry.add("name", plot::json_string(name));
                const bool single = write_heatmap(heat, tile_settings, entry);
                if (single) {
                    entry.add("start_pos", plot::json_ints(chromosome.start_pos1));
                    entry.add("start_pos2", plot::json_ints(chromosome.start_pos2));
                } else {
                    entry.add("start_pos", "null");
                    entry.add("start_pos2", "null");
                }
                entry.add("region", region_json(chromosome.region));
                entry.add("bigwig_tracks", bins.size() >= 5
                                               ? bigwig_json(chromosome.region)
                                               : std::string("null"));
                entries.push_back(entry.str());
            }
            data.add("chromosomes", plot::json_list(entries));
        } else {
            Heatmap heat;
            heat.first = primary;
            heat.second = matrix2_path.has_value() ? &second.source : nullptr;
            if (selection.has_value()) {
                heat.start_pos_rows = selection->start_pos1;
                heat.start_pos_cols = selection->start_pos2;
            } else {
                heat.start_pos_rows = make_start_pos_array(model);
                heat.start_pos_cols = heat.start_pos_rows;
            }
            heat.pixels_x = bigwig_count == 0 ? 5.0 * dpi : figure_width * dpi;
            heat.pixels_y = bigwig_count == 0 ? 4.8 * dpi : figure_height * dpi;
            const bool single = write_heatmap(heat, tile_settings, data);
            if (single) {
                data.add("start_pos", plot::json_ints(heat.start_pos_rows));
                data.add("start_pos2", selection.has_value() ? plot::json_ints(selection->start_pos2)
                                                             : std::string("null"));
            } else {
                data.add("start_pos", "null");
                data.add("start_pos2", "null");
            }
            if (selection.has_value()) {
                data.add("region", region_json(selection->region));
            } else {
                data.add("region", "null");
            }
            std::optional<Region> whole_region;
            if (selection.has_value()) {
                whole_region = selection->region;
            }
            const bool heatmap_drawn = matrix_rows >= 5;
            data.add("bigwig_tracks", heatmap_drawn ? bigwig_json(whole_region) : "null");
            const std::optional<std::string> loops = ns.opt_str("loops");
            const std::optional<std::string> tads = ns.opt_str("tads");
            data.add("loop_positions",
                     heatmap_drawn && loops.has_value() && !loops->empty()
                         ? loop_positions_json(*loops, whole_region, chromosome_order,
                                               chromosome_order_as_bytes,
                                               ns.str("loopLargeRegionsOperation"))
                         : "null");
            data.add("tad_positions",
                     heatmap_drawn && tads.has_value() && !tads->empty()
                         ? tad_positions_json(*tads, whole_region, chromosome_order,
                                              chromosome_order_as_bytes)
                         : "null");
            data.add("xlabel", xlabel.has_value() ? plot::json_string(*xlabel) : "null");
            data.add("ylabel", ylabel.has_value() ? plot::json_string(*ylabel) : "null");
        }
        data.add("temporary_files", plot::json_strings(temporary_files));
        model = Model();
        second = SecondMatrix();

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
