#include "hicx/plot_tiles.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hicx::plot {

void MaskStats::add(double value) noexcept {
    ++cells;
    if (std::isnan(value)) {
        ++nans;
        return;
    }
    if (value == 0.0) {
        ++zeros;
        return;
    }
    if (!has_nonzero_value || value < nonzero_min) {
        nonzero_min = value;
        has_nonzero_value = true;
    }
    if (std::isinf(value)) {
        ++infs;
        return;
    }
    if (!has_finite_nonzero || value < finite_nonzero_min) {
        finite_nonzero_min = value;
        has_finite_nonzero = true;
    }
    if (value != std::floor(value)) {
        all_integral = false;
    }
    max_abs_finite = std::max(max_abs_finite, std::fabs(value));
}

void MaskStats::merge(const MaskStats& other) noexcept {
    cells += other.cells;
    zeros += other.zeros;
    nans += other.nans;
    infs += other.infs;
    if (other.has_nonzero_value && (!has_nonzero_value || other.nonzero_min < nonzero_min)) {
        nonzero_min = other.nonzero_min;
        has_nonzero_value = true;
    }
    if (other.has_finite_nonzero &&
        (!has_finite_nonzero || other.finite_nonzero_min < finite_nonzero_min)) {
        finite_nonzero_min = other.finite_nonzero_min;
        has_finite_nonzero = true;
    }
    all_integral = all_integral && other.all_integral;
    max_abs_finite = std::max(max_abs_finite, other.max_abs_finite);
}

double MaskParams::apply(double value) const noexcept {
    if (value == 0.0) {
        if (!replace_zero) {
            return value;
        }
        value = zero_value;
        if (!chain_zero) {
            return value;
        }
    }
    if (std::isnan(value)) {
        return replace_nan ? nan_value : value;
    }
    if (std::isinf(value)) {
        return replace_inf ? inf_value : value;
    }
    return value;
}

MaskParams::Summary MaskParams::summarize(const MaskStats& stats) const noexcept {
    Summary out;
    out.integral = stats.all_integral;
    out.max_abs = stats.max_abs_finite;
    auto visit = [&](double value) {
        if (!std::isfinite(value)) {
            out.nonfinite = true;
        } else {
            if (value != std::floor(value)) {
                out.integral = false;
            }
            out.max_abs = std::max(out.max_abs, std::fabs(value));
        }
    };
    if (stats.zeros > 0) {
        visit(apply(0.0));
    }
    if (stats.nans > 0) {
        visit(apply(std::numeric_limits<double>::quiet_NaN()));
    }
    if (stats.infs > 0) {
        visit(apply(std::numeric_limits<double>::infinity()));
    }
    return out;
}

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// np.nanmin over a few candidate values that are present: NaN when there is
// none that is not NaN.
double smallest_of(const double* values, const bool* present, int count) {
    double best = kNaN;
    for (int i = 0; i < count; ++i) {
        if (present[i] && !std::isnan(values[i]) && (std::isnan(best) || values[i] < best)) {
            best = values[i];
        }
    }
    return best;
}

[[noreturn]] void empty_reduction() {
    throw std::runtime_error("ValueError: zero-size array to reduction operation fmin which has "
                             "no identity");
}

}  // namespace

MaskParams mask_params_whole(const MaskStats& s) {
    MaskParams p;
    p.chain_zero = true;
    // Zeros take the smallest non-zero value (the smallest normal double when
    // every cell is zero, NaN when every non-zero cell is NaN).
    const std::int64_t nonzero_cells = s.cells - s.zeros;
    const double s1 = nonzero_cells == 0 ? std::numeric_limits<double>::min()
                                         : (s.has_nonzero_value ? s.nonzero_min : kNaN);
    p.replace_zero = true;
    p.zero_value = s1;
    const bool any_nan = s.nans > 0 || (s.zeros > 0 && std::isnan(s1));
    const bool any_inf = s.infs > 0 || (s.zeros > 0 && std::isinf(s1));
    if (!any_nan && !any_inf) {
        return p;
    }
    // NaN cells take the smallest value that is not NaN.
    const std::int64_t not_nan = s.cells - s.nans - (std::isnan(s1) ? s.zeros : 0);
    if (not_nan == 0) {
        empty_reduction();
    }
    {
        const double values[2] = {s.nonzero_min, s1};
        const bool present[2] = {s.has_nonzero_value, s.zeros > 0};
        p.nan_value = smallest_of(values, present, 2);
        p.replace_nan = true;
    }
    // Infinite cells (as they were before the NaN replacement) take the
    // smallest value among the cells that were not infinite; cells that only
    // became infinite through the NaN replacement count there and are not
    // replaced.
    const std::int64_t finite_nonzero = s.cells - s.zeros - s.nans - s.infs;
    const std::int64_t selected =
        finite_nonzero + (std::isinf(s1) ? 0 : s.zeros) + s.nans;
    if (selected == 0) {
        empty_reduction();
    }
    {
        const double zero_now = std::isnan(s1) ? p.nan_value : s1;
        const double values[3] = {s.finite_nonzero_min, zero_now, p.nan_value};
        const bool present[3] = {s.has_finite_nonzero, s.zeros > 0 && !std::isinf(s1), s.nans > 0};
        p.inf_value = smallest_of(values, present, 3);
        p.replace_inf = true;
    }
    return p;
}

MaskParams mask_params_per_chromosome(const MaskStats& s) {
    MaskParams p;
    p.chain_zero = false;
    if (s.cells - s.zeros == 0) {
        return p;
    }
    const double s1 = s.has_nonzero_value ? s.nonzero_min : kNaN;
    p.replace_zero = true;
    p.zero_value = s1;
    if (s.cells - s.nans == 0) {
        return p;
    }
    {
        const double values[2] = {s.nonzero_min, s1};
        const bool present[2] = {s.has_nonzero_value, s.zeros > 0};
        p.nan_value = smallest_of(values, present, 2);
        p.replace_nan = true;
    }
    if (s.cells - s.infs == 0) {
        return p;
    }
    {
        const double values[3] = {s.finite_nonzero_min, s1, p.nan_value};
        const bool present[3] = {s.has_finite_nonzero, s.zeros > 0, s.nans > 0};
        p.inf_value = smallest_of(values, present, 3);
        p.replace_inf = true;
    }
    return p;
}

HandoffChoice choose_handoff(bool all_integral, bool has_nonfinite, double max_abs_value) {
    HandoffChoice choice;
    if (!all_integral) {
        choice.type = Handoff::Float32;
        return choice;
    }
    if (!has_nonfinite) {
        if (max_abs_value <= 2147483647.0) {
            choice.type = Handoff::Int32;
        } else if (max_abs_value < 9.0e18) {
            choice.type = Handoff::Int64;
            choice.widened = true;
        } else {
            choice.type = Handoff::Float64;
            choice.widened = true;
        }
        return choice;
    }
    if (max_abs_value <= 16777216.0) {
        choice.type = Handoff::Float32;
    } else {
        choice.type = Handoff::Float64;
        choice.widened = true;
    }
    return choice;
}

const char* npy_dtype(Handoff type) noexcept {
    switch (type) {
        case Handoff::Int32:
            return "<i4";
        case Handoff::Int64:
            return "<i8";
        case Handoff::Float32:
            return "<f4";
        case Handoff::Float64:
            return "<f8";
    }
    return "<f8";
}

std::size_t element_size(Handoff type) noexcept {
    return (type == Handoff::Int32 || type == Handoff::Float32) ? 4 : 8;
}

std::int64_t positions_per_block(std::int64_t bins, std::int64_t pixels, bool downsample) noexcept {
    if (!downsample || pixels < 1 || bins <= pixels) {
        return 1;
    }
    return (bins + pixels - 1) / pixels;
}

AxisBlocks make_blocks(std::int64_t n, std::int64_t k, const std::vector<std::int32_t>& labels,
                       std::int64_t max_blocks) {
    AxisBlocks out;
    out.k = std::max<std::int64_t>(k, 1);
    auto build = [&](bool by_label) {
        out.first.clear();
        out.block_of.assign(static_cast<std::size_t>(n), 0);
        out.by_chromosome = by_label && out.k > 1;
        std::int64_t start = 0;
        for (std::int64_t p = 0; p < n; ++p) {
            const bool new_block =
                p == 0 || p - start >= out.k ||
                (by_label && labels[static_cast<std::size_t>(p)] !=
                                 labels[static_cast<std::size_t>(p - 1)]);
            if (new_block) {
                out.first.push_back(p);
                start = p;
            }
            out.block_of[static_cast<std::size_t>(p)] =
                static_cast<std::int32_t>(out.first.size() - 1);
        }
        out.first.push_back(n);
    };
    const bool can_label = out.k > 1 && labels.size() == static_cast<std::size_t>(n);
    build(can_label);
    if (can_label && out.count() > max_blocks) {
        build(false);
    }
    return out;
}

std::vector<double> nearest_edges(const std::vector<std::int64_t>& centres) {
    const std::size_t n = centres.size();
    std::vector<double> edges;
    if (n == 0) {
        return edges;
    }
    if (n == 1) {
        return {static_cast<double>(centres[0]), static_cast<double>(centres[0])};
    }
    edges.resize(n + 1);
    auto half_step = [&](std::size_t i) {
        return static_cast<double>(centres[i + 1] - centres[i]) * 0.5;
    };
    edges[0] = static_cast<double>(centres[0]) - half_step(0);
    for (std::size_t i = 1; i < n; ++i) {
        edges[i] = static_cast<double>(centres[i - 1]) + half_step(i - 1);
    }
    edges[n] = static_cast<double>(centres[n - 1]) + half_step(n - 2);
    return edges;
}

std::vector<double> block_edges(const std::vector<double>& edges, const AxisBlocks& blocks) {
    std::vector<double> out;
    const std::int64_t count = blocks.count();
    out.reserve(static_cast<std::size_t>(count) + 1);
    for (std::int64_t b = 0; b < count; ++b) {
        out.push_back(edges[static_cast<std::size_t>(blocks.first[static_cast<std::size_t>(b)])]);
    }
    out.push_back(edges.back());
    return out;
}

}  // namespace hicx::plot
