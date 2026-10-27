// Unit tests for hicx/find_tads.hpp. The general forms the Python API uses
// must agree with the kernels the hicFindTADs tool runs wherever the tool's
// special case applies; the expected values of the other cases were computed
// with HiCExplorer 3.7.6 (hicexplorer/hicFindTADs.py).

#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include <doctest/doctest.h>

#include "hicx/find_tads.hpp"
#include "hicx/python_error.hpp"
#include "hicx/stats_ops.hpp"

namespace ft = hicx::find_tads;

TEST_CASE("find_tads: incremental step sizes over integers and floats agree") {
    // get_incremental_step_size(8, 50, 2) == [8, 10, 13, 18, 24, 30, 37, 45]
    const std::vector<std::int64_t> ints = ft::incremental_step_size(8, 50, 2);
    CHECK(ints == std::vector<std::int64_t>{8, 10, 13, 18, 24, 30, 37, 45});
    const std::vector<double> reals = ft::incremental_step_size_real(8, 50, 2);
    REQUIRE(reals.size() == ints.size());
    for (std::size_t i = 0; i < ints.size(); ++i) {
        CHECK(reals[i] == static_cast<double>(ints[i]));
    }
    CHECK(ft::incremental_step_size_real(10, 5, 2).empty());
    CHECK_THROWS_AS((void)ft::incremental_step_size_real(1, 5, 0), hicx::PythonError);
    CHECK_THROWS_AS((void)ft::incremental_step_size_real(1, 5, std::nan("")), hicx::PythonError);
}

TEST_CASE("find_tads: diagflat_zero drops the diagonal and stored zeros, keeps NaN") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    // [[1, 0*, 2], [0, nan, 3], [4, 0, 5]] with an explicit zero at (0, 1)
    hicx::CsrMatrix m(3, 3, {0, 3, 5, 7}, {0, 1, 2, 1, 2, 0, 2}, {1.0, 0.0, 2.0, nan, 3.0, 4.0, 5.0},
                      "float32");
    const hicx::CsrMatrix out = ft::diagflat_zero(m);
    CHECK(out.dtype() == "float64");
    CHECK(out.indptr() == std::vector<std::int64_t>{0, 1, 3, 4});
    CHECK(out.indices() == std::vector<std::int32_t>{2, 1, 2, 0});
    CHECK(out.data()[0] == 2.0);
    CHECK(std::isnan(out.data()[1]));
    CHECK(out.data()[2] == 3.0);
    CHECK(out.data()[3] == 4.0);
    hicx::CsrMatrix ints(2, 2, {0, 2, 3}, {0, 1, 1}, {7.0, 1.0, 2.0}, "int32");
    CHECK(ft::diagflat_zero(ints).dtype() == "int64");
}

TEST_CASE("find_tads: the general peakdetect agrees with the tool's") {
    const std::vector<double> y = {0.5, 0.2, -0.3, -0.1, 0.4, 0.6, 0.1, -0.5, -0.2, 0.3,
                                   0.7, 0.2, -0.4, 0.0,  0.5, 0.1, 0.2, 0.3};
    const std::vector<std::string> chrom(y.size(), "chr1");
    for (std::int64_t lookahead = 1; lookahead <= 4; ++lookahead) {
        std::vector<ft::Peak> maxima;
        std::vector<ft::Peak> minima;
        ft::peakdetect(y, chrom, lookahead, maxima, minima);
        const ft::PeakResult general =
            ft::peakdetect_general(y, &chrom, lookahead, 0.0, true, true, "ndarray");
        REQUIRE(general.minima.size() == minima.size());
        REQUIRE(general.maxima.size() == maxima.size());
        for (std::size_t i = 0; i < minima.size(); ++i) {
            CHECK(general.minima[i].position == minima[i].position);
            CHECK(y[static_cast<std::size_t>(general.minima[i].value_index)] == minima[i].value);
        }
        for (std::size_t i = 0; i < maxima.size(); ++i) {
            CHECK(general.maxima[i].position == maxima[i].position);
        }
    }
    // lookahead 0 iterates over y[:-0], which is empty
    CHECK(ft::peakdetect_general(y, &chrom, 0, 0.0, true, true, "ndarray").minima.empty());
    // chrom=None raises as soon as the loop runs
    CHECK_THROWS_AS((void)ft::peakdetect_general(y, nullptr, 3, 0.0, true, true, "ndarray"),
                    hicx::PythonError);
    std::vector<double> with_inf = y;
    with_inf[2] = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS((void)ft::peakdetect_general(with_inf, &chrom, 3, 0.0, true, true, "ndarray"),
                    hicx::PythonError);
}

TEST_CASE("find_tads: delta_wrt_window values agree with the tool's map") {
    std::vector<double> scores;
    std::vector<std::string> chrom;
    for (int i = 0; i < 60; ++i) {
        scores.push_back(std::sin(i * 0.37) + 0.01 * i);
        chrom.push_back(i < 35 ? "chr2" : "chr1");
    }
    const std::vector<std::int64_t> minima = {3, 12, 20, 33, 40, 47, 58};
    const std::map<std::int64_t, double> map = ft::delta_wrt_window(minima, scores, chrom, 10);
    const std::vector<double> values = ft::delta_wrt_window_values(minima, scores, chrom, 10);
    for (std::size_t i = 0; i < minima.size(); ++i) {
        const double expected = map.at(minima[i]);
        if (std::isnan(expected)) {
            CHECK(std::isnan(values[i]));
        } else {
            CHECK(values[i] == expected);
        }
    }
}

TEST_CASE("find_tads: the FDR cutoff agrees with stats::benjamini_hochberg_cutoff") {
    const std::vector<double> p = {0.001, 0.2, 0.004, 0.03, 0.5, 0.0005, 0.011};
    const std::optional<double> cutoff = ft::fdr_cutoff(p, 0.05);
    REQUIRE(cutoff.has_value());
    CHECK(*cutoff == hicx::stats::benjamini_hochberg_cutoff(p, 0.05));
    CHECK_FALSE(ft::fdr_cutoff({0.9, 0.8}, 0.01).has_value());
}

TEST_CASE("find_tads: numpy statistics and np.split") {
    CHECK(ft::numpy_median({3.0, 1.0, 2.0, 10.0}) == 2.5);
    CHECK(ft::numpy_percentile({1.0, 2.0, 3.0, 4.0}, 75) == 3.25);
    CHECK(ft::numpy_percentile({1.0, 2.0, 3.0, 4.0}, 25) == 1.75);
    // np.split(np.arange(6), [4, 2]) gives [0:4], [4:2] (empty) and [2:6]
    const auto pieces = ft::numpy_split_bounds(6, {4, 2});
    REQUIRE(pieces.size() == 3);
    CHECK(pieces[0] == std::pair<std::size_t, std::size_t>{0, 4});
    CHECK(pieces[1].first == pieces[1].second);
    CHECK(pieces[2] == std::pair<std::size_t, std::size_t>{2, 6});
}

TEST_CASE("find_tads: hierarchical clustering of one chromosome") {
    // boundaries at 0, 100, 250, 400, 600 with values 1, 0.5, 0.9, 0.2, 1:
    // the inner values are 0.5, 0.9, 0.2, merged from the highest down. 3.7.6
    // gives z (1, 2, 0.9, 2, 175, 325), (0, 4, 0.5, 3, 50, 250),
    // (5, 3, 0.2, 4, 150, 500) and, at cutoff 0.3, the clusters (0, 100),
    // (100, 400), (400, 600).
    ft::ClusterInput input;
    input.start = {0, 100, 250, 400, 600};
    input.value = {1, 0.5, 0.9, 0.2, 1};
    input.start_integer = true;
    const ft::ClusterResult result = ft::cluster_chromosome(input, {0.3});
    REQUIRE(result.merges.size() == 3);
    CHECK(result.merges[0].left == 1);
    CHECK(result.merges[0].right == 2);
    CHECK(result.merges[0].value_index == 2);
    CHECK(result.merges[0].size == 2);
    CHECK(result.merges[0].x_left == 175);
    CHECK(result.merges[0].x_right == 325);
    CHECK(result.merges[1].left == 0);
    CHECK(result.merges[1].right == 4);
    CHECK(result.merges[1].x_left == 50);
    CHECK(result.merges[1].x_right == 250);
    CHECK(result.merges[2].left == 5);
    CHECK(result.merges[2].right == 3);
    CHECK(result.merges[2].size == 4);
    CHECK(result.merges[2].x_left == 150);
    CHECK(result.merges[2].x_right == 500);
    REQUIRE(result.levels.size() == 1);
    REQUIRE(result.levels[0].size() == 3);
    CHECK(result.levels[0][1].start_domain == 1);
    CHECK(result.levels[0][1].end_domain == 2);
    ft::ClusterInput empty;
    CHECK_THROWS_AS((void)ft::cluster_chromosome(empty, {}), hicx::PythonError);
}

TEST_CASE("find_tads: format_fixed writes NaN as Python does") {
    CHECK(ft::format_fixed(-std::numeric_limits<double>::quiet_NaN(), 12) == "nan");
    CHECK(ft::format_fixed(std::numeric_limits<double>::infinity(), 3) == "inf");
    CHECK(ft::format_fixed(0.5, 3) == "0.500");
}
