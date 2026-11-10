#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <random>
#include <stdexcept>
#include <vector>

#include "hicx/plot_tiles.hpp"

namespace {

using hicx::plot::Handoff;
using hicx::plot::MaskParams;
using hicx::plot::MaskStats;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

// The dense implementations of hicPlotMatrix.cpp before the tiled handoff
// (main() and plotPerChr() of the Python tool), kept here as the reference.
std::optional<double> nanmin(const std::vector<double>& values, const std::vector<char>& excluded) {
    bool any = false;
    double smallest = kNaN;
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

// Returns false where the Python tool raises.
bool dense_whole(std::vector<double>& values) {
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
            return false;
        }
        replace_masked(values, nan_mask, *without_nan);
        const std::optional<double> without_inf = nanmin(values, inf_mask);
        if (!without_inf.has_value()) {
            return false;
        }
        replace_masked(values, inf_mask, *without_inf);
    }
    return true;
}

void dense_per_chromosome(std::vector<double>& values) {
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

bool same(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || a == b;
}

MaskStats stats_of(const std::vector<double>& values) {
    MaskStats s;
    // in two parts, to exercise merge()
    MaskStats t;
    for (std::size_t i = 0; i < values.size(); ++i) {
        (i % 2 == 0 ? s : t).add(values[i]);
    }
    s.merge(t);
    return s;
}

std::vector<double> random_cells(std::mt19937& rng, std::size_t n) {
    static const double pool[] = {0.0, 0.0, 0.0, 1.0, 2.0, 3.0, -1.5, 7.25, kNaN, kInf, -kInf};
    std::uniform_int_distribution<int> kind(0, 9);
    std::uniform_int_distribution<std::size_t> pick(0, sizeof(pool) / sizeof(pool[0]) - 1);
    std::uniform_real_distribution<double> real(0.001, 50.0);
    // each draw picks a restricted pool so the cases without NaN, without
    // infinities, and with only zeros all occur
    const int mode = static_cast<int>(rng() % 6);
    std::vector<double> v(n);
    for (double& x : v) {
        const double candidate = kind(rng) < 3 ? real(rng) : pool[pick(rng)];
        if ((mode == 0 && std::isnan(candidate)) || (mode == 1 && std::isinf(candidate)) ||
            (mode == 2 && (std::isnan(candidate) || std::isinf(candidate))) ||
            (mode == 3 && candidate != 0.0)) {
            x = 0.0;
        } else if (mode == 4 && candidate != 0.0 && !std::isnan(candidate)) {
            x = kNaN;
        } else {
            x = candidate;
        }
    }
    return v;
}

}  // namespace

TEST_CASE("mask parameters from statistics equal the dense replacement, whole matrix") {
    std::mt19937 rng(12345);
    int raised = 0;
    for (int round = 0; round < 4000; ++round) {
        const std::vector<double> cells = random_cells(rng, 1 + rng() % 12);
        std::vector<double> expected = cells;
        const bool ok = dense_whole(expected);
        const MaskStats stats = stats_of(cells);
        if (!ok) {
            ++raised;
            CHECK_THROWS_AS(hicx::plot::mask_params_whole(stats), std::runtime_error);
            continue;
        }
        const MaskParams params = hicx::plot::mask_params_whole(stats);
        bool nonfinite = false;
        for (std::size_t i = 0; i < cells.size(); ++i) {
            const double got = params.apply(cells[i]);
            CHECK(same(got, expected[i]));
            nonfinite = nonfinite || !std::isfinite(expected[i]);
        }
        CHECK(params.summarize(stats).nonfinite == nonfinite);
    }
    CHECK(raised > 0);
}

TEST_CASE("mask parameters from statistics equal the dense replacement, per chromosome") {
    std::mt19937 rng(777);
    for (int round = 0; round < 4000; ++round) {
        const std::vector<double> cells = random_cells(rng, 1 + rng() % 12);
        std::vector<double> expected = cells;
        dense_per_chromosome(expected);
        const MaskParams params = hicx::plot::mask_params_per_chromosome(stats_of(cells));
        for (std::size_t i = 0; i < cells.size(); ++i) {
            CHECK(same(params.apply(cells[i]), expected[i]));
        }
    }
}

TEST_CASE("integer counts go to int32 and fall back to int64 above INT32_MAX") {
    using hicx::plot::choose_handoff;
    CHECK(choose_handoff(true, false, 0.0).type == Handoff::Int32);
    CHECK(choose_handoff(true, false, 2147483647.0).type == Handoff::Int32);
    CHECK_FALSE(choose_handoff(true, false, 2147483647.0).widened);
    const auto wide = choose_handoff(true, false, 2147483648.0);
    CHECK(wide.type == Handoff::Int64);
    CHECK(wide.widened);
    CHECK(choose_handoff(true, false, 1.0e19).type == Handoff::Float64);
    // an integer matrix that holds NaN cannot be an integer array
    CHECK(choose_handoff(true, true, 5.0).type == Handoff::Float32);
    CHECK(choose_handoff(true, true, 16777216.0).type == Handoff::Float32);
    const auto big_nan = choose_handoff(true, true, 16777217.0);
    CHECK(big_nan.type == Handoff::Float64);
    CHECK(big_nan.widened);
    // balanced and mean values
    CHECK(choose_handoff(false, false, 3.0).type == Handoff::Float32);
    CHECK(choose_handoff(false, true, 3.0).type == Handoff::Float32);
    CHECK(std::string(hicx::plot::npy_dtype(Handoff::Int32)) == "<i4");
}

TEST_CASE("statistics find non-integral values and the largest magnitude") {
    MaskStats s;
    for (const double v : {0.0, 3.0, -2147483648.0, kNaN}) {
        s.add(v);
    }
    CHECK(s.all_integral);
    CHECK(s.max_abs_finite == 2147483648.0);
    s.add(0.5);
    CHECK_FALSE(s.all_integral);
}

TEST_CASE("blocks stop at chromosome boundaries and cover every position") {
    const std::vector<std::int32_t> labels = {0, 0, 0, 0, 0, 1, 1, 1, 2};
    const auto blocks = hicx::plot::make_blocks(9, 2, labels, 100);
    CHECK(blocks.by_chromosome);
    const std::vector<std::int64_t> expected = {0, 2, 4, 5, 7, 8, 9};
    CHECK(blocks.first == expected);
    CHECK(blocks.block_of[4] == 2);
    CHECK(blocks.block_of[8] == 5);
    // too many blocks per chromosome: plain blocks
    const auto plain = hicx::plot::make_blocks(9, 2, labels, 5);
    CHECK_FALSE(plain.by_chromosome);
    CHECK(plain.count() == 5);
    CHECK(hicx::plot::positions_per_block(33754, 576, true) == 59);
    CHECK(hicx::plot::positions_per_block(500, 576, true) == 1);
    CHECK(hicx::plot::positions_per_block(33754, 576, false) == 1);
}

TEST_CASE("nearest edges are the midpoints with the ends extended by half a step") {
    const auto e = hicx::plot::nearest_edges({0, 10, 20, 35});
    const std::vector<double> expected = {-5.0, 5.0, 15.0, 27.5, 42.5};
    CHECK(e == expected);
    const auto blocks = hicx::plot::make_blocks(4, 2, {}, 10);
    const auto b = hicx::plot::block_edges(e, blocks);
    const std::vector<double> expected_blocks = {-5.0, 15.0, 42.5};
    CHECK(b == expected_blocks);
}
