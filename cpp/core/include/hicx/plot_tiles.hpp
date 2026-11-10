// Helpers of the tiled, reduced dense handoff of hicPlotMatrix (see the file
// comment of cpp/tools/hicPlotMatrix.cpp): the replacement of zeros, NaN and
// infinite values for --log and --log1p from streamed statistics instead of
// a dense array, the choice of the element type of the .npy handed to the
// drawing layer, and the block layout of an axis.
//
// Nothing here allocates per cell; every function works on a few numbers or
// on one vector with an entry per bin.

#ifndef HICX_PLOT_TILES_HPP
#define HICX_PLOT_TILES_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace hicx::plot {

// ---------------------------------------------------------------------------
// Statistics of the cells of one matrix, gathered tile by tile.

struct MaskStats {
    std::int64_t cells = 0;
    std::int64_t zeros = 0;       // cells equal to 0
    std::int64_t nans = 0;
    std::int64_t infs = 0;        // +inf and -inf
    bool has_nonzero_value = false;   // some non-zero, non-NaN cell
    double nonzero_min = 0.0;         // smallest of them (infinities included)
    bool has_finite_nonzero = false;  // some finite non-zero cell
    double finite_nonzero_min = 0.0;
    // For the choice of the element type of the handoff.
    bool all_integral = true;     // every finite cell holds an integer
    double max_abs_finite = 0.0;

    void add(double value) noexcept;
    void merge(const MaskStats& other) noexcept;
};

// The replacement values of --log / --log1p. A cell equal to 0 becomes
// zero_value, then a NaN cell nan_value, then an infinite cell inf_value, each
// only when its flag is set (hicPlotMatrix.py replaces them in this order).
struct MaskParams {
    // Whole matrix: a zero replaced by an infinite or NaN value is replaced
    // again by the next step. Per chromosome: the masks come from the original
    // cells, so it is not.
    bool chain_zero = true;
    bool replace_zero = false;
    double zero_value = 0.0;
    bool replace_nan = false;
    double nan_value = 0.0;
    bool replace_inf = false;
    double inf_value = 0.0;

    [[nodiscard]] double apply(double value) const noexcept;
    // What the cells of `stats` hold after apply(): whether a NaN or an
    // infinite value remains, whether all finite cells are integers, and the
    // largest finite absolute value.
    struct Summary {
        bool nonfinite = false;
        bool integral = true;
        double max_abs = 0.0;
    };
    [[nodiscard]] Summary summarize(const MaskStats& stats) const noexcept;
};

// main() of hicexplorer/hicPlotMatrix.py on the whole matrix: zeros take the
// smallest non-zero value, then NaN and infinite values the smallest value that
// is not NaN, not infinite. Throws std::runtime_error("ValueError: ...") where
// numpy raises on an empty reduction.
[[nodiscard]] MaskParams mask_params_whole(const MaskStats& stats);

// plotPerChr(): the three masks are taken before any replacement, and a
// failing replacement stops the rest silently.
[[nodiscard]] MaskParams mask_params_per_chromosome(const MaskStats& stats);

// ---------------------------------------------------------------------------
// The element type of the handoff.

enum class Handoff { Int32, Int64, Float32, Float64 };

struct HandoffChoice {
    Handoff type = Handoff::Float32;
    // Set when integer counts did not fit the preferred type (int32, or the
    // 24 bit mantissa of float32 next to NaN) and a wider one was chosen.
    bool widened = false;
};

// Integer count matrices go to int32, accumulated in int64 and converted at the
// handoff; when a value would exceed INT32_MAX the matrix is handed over as
// int64 (widened = true). An integer matrix that must hold NaN or infinity
// (which an integer array cannot) goes to float32, which is exact up to 2^24,
// and to float64 beyond that (widened = true). Everything else is float32.
[[nodiscard]] HandoffChoice choose_handoff(bool all_integral, bool has_nonfinite,
                                           double max_abs_value);

[[nodiscard]] const char* npy_dtype(Handoff type) noexcept;
[[nodiscard]] std::size_t element_size(Handoff type) noexcept;

// ---------------------------------------------------------------------------
// Blocks of an axis: the positions 0..n-1 are cut into runs of consecutive
// positions. `labels[p]` identifies the chromosome of position p; a run never
// holds two labels when `by_label` is set.

struct AxisBlocks {
    std::vector<std::int64_t> first;     // first position of each block, plus n
    std::vector<std::int32_t> block_of;  // block of each position
    std::int64_t k = 1;                  // requested positions per block
    bool by_chromosome = false;          // blocks stop at chromosome boundaries
    [[nodiscard]] std::int64_t count() const noexcept {
        return static_cast<std::int64_t>(first.size()) - 1;
    }
};

// Blocks of k positions. With `labels` non-empty the blocks follow the runs of
// equal labels; if that gives more than max_blocks blocks the label-free
// blocking of k positions is used instead (by_chromosome = false).
[[nodiscard]] AxisBlocks make_blocks(std::int64_t n, std::int64_t k,
                                     const std::vector<std::int32_t>& labels,
                                     std::int64_t max_blocks);

// k = ceil(bins / pixels) when bins > pixels, else 1.
[[nodiscard]] std::int64_t positions_per_block(std::int64_t bins, std::int64_t pixels,
                                               bool downsample) noexcept;

// The cell edges pcolormesh(shading="nearest") derives from cell centres:
// midpoints between them, the ends extended by half the last difference.
[[nodiscard]] std::vector<double> nearest_edges(const std::vector<std::int64_t>& centres);

// The edges of the blocks: edges[first[b]] and the last edge.
[[nodiscard]] std::vector<double> block_edges(const std::vector<double>& edges,
                                              const AxisBlocks& blocks);

}  // namespace hicx::plot

#endif
