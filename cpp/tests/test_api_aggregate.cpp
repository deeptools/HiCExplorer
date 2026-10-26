// Unit tests of the core code behind the Python API of hicAverageRegions,
// hicCompartmentalization, hicPCA, hicPlotDistVsCounts, hicPlotSVL,
// hicPlotMatrix and hicAggregateContacts (cpp/python/API_BINDINGS.md).
// The equivalence with 3.7.6 is tested by hicexplorer/test/api/cases_*.py;
// these cases pin the arithmetic and the clamping rules.

#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>

#include "hicx/average_regions.hpp"

namespace ar = hicx::average_regions;

TEST_CASE("average_regions: centre truncates towards zero") {
    CHECK(ar::center_of(10.0, 21.0) == 15);
    CHECK(ar::center_of(2500001.5, 3500000.0) == 3000000);
    CHECK(ar::center_of(-3.0, 0.0) == -1);
}

TEST_CASE("average_regions: genomic window clamps to 0 and to the chromosome end minus one") {
    ar::Window w = ar::viewpoint_range(100, 30, 50, 1000);
    CHECK(w.start == 70);
    CHECK(w.end == 150);
    CHECK_FALSE(w.start_out_of_range);
    CHECK_FALSE(w.end_out_of_range);
    w = ar::viewpoint_range(10, 30, 5000, 1000);
    CHECK(w.start == 0);
    CHECK(w.end == 999);
    CHECK(w.start_out_of_range);
    CHECK(w.end_out_of_range);
    // reaching exactly the end is not out of range
    w = ar::viewpoint_range(500, 0, 500, 1000);
    CHECK(w.end == 1000);
    CHECK_FALSE(w.end_out_of_range);
}

TEST_CASE("average_regions: bin window clamps to the chromosome's bins, last bin inclusive") {
    ar::Window w = ar::viewpoint_range_bins(105, 10, 20, 100, 120);
    CHECK(w.start == 100);
    CHECK(w.end == 120);
    CHECK(w.start_out_of_range);
    CHECK(w.end_out_of_range);
    w = ar::viewpoint_range_bins(110, 10, 10, 100, 120);
    CHECK(w.start == 100);
    CHECK(w.end == 120);
    CHECK_FALSE(w.start_out_of_range);
    CHECK_FALSE(w.end_out_of_range);
}

#include "hicx/sparse_matrix.hpp"
#include "hicx/svl.hpp"

namespace {

// [[1 2 0], [2 3 4], [0 4 5]] as the full matrix and as its upper triangle.
hicx::CsrMatrix svl_fixture(bool upper, const std::string& dtype) {
    if (upper) {
        hicx::CsrMatrix m(3, 3, {0, 2, 4, 5}, {0, 1, 1, 2, 2}, {1, 2, 3, 4, 5}, dtype);
        m.set_symmetry(hicx::Symmetry::UpperTriangle);
        return m;
    }
    return hicx::CsrMatrix(3, 3, {0, 2, 5, 7}, {0, 1, 0, 1, 2, 1, 2}, {1, 2, 2, 3, 4, 4, 5}, dtype);
}

}  // namespace

TEST_CASE("svl: short and long range sums of a block") {
    for (const bool upper : {false, true}) {
        const hicx::svl::BlockSums ints = hicx::svl::block_sums(svl_fixture(upper, "int32"), 0, 3, 0.0);
        CHECK(ints.kind == hicx::DType::Integer);
        CHECK(ints.smaller_int == 9);
        CHECK(ints.greater_int == 12);
        CHECK(ints.ratio == 9.0 / 12.0);
        CHECK(ints.kept);
        const hicx::svl::BlockSums floats = hicx::svl::block_sums(svl_fixture(upper, "float32"), 0, 3, 1.0);
        CHECK(floats.kind == hicx::DType::Float32);
        CHECK(floats.smaller == 21.0);
        CHECK(floats.greater == 0.0);
        CHECK_FALSE(floats.kept);
        const hicx::svl::BlockSums block = hicx::svl::block_sums(svl_fixture(upper, "float64"), 1, 3, 0.5);
        CHECK(block.smaller == 8.0);
        CHECK(block.greater == 8.0);
        CHECK(block.ratio == 1.0);
    }
}

TEST_CASE("svl: stored zeros of a full block are the Python's IndexError") {
    hicx::CsrMatrix m(2, 2, {0, 2, 3}, {0, 1, 1}, {1, 0, 3}, "float64");
    CHECK_THROWS_AS((void)hicx::svl::block_sums(m, 0, 2, 1.0), hicx::svl::ExplicitZeros);
    // outside the block they do not matter
    CHECK(hicx::svl::block_sums(m, 1, 2, 1.0).smaller == 3.0);
}

#include <vector>

#include "hicx/compartmentalization.hpp"

TEST_CASE("compartments: within_vs_between of a rectangular array stops at its edges") {
    // [[1 2 3], [4 5 6]] with four quantiles: q = 1 sums a[0:1, 0:1] + a[3:4, 3:4]
    // (empty) over a[0:1, 3:4] (empty) + a[3:4, 0:1] (empty)
    const std::vector<double> values{1, 2, 3, 4, 5, 6};
    const std::vector<double> ratios = hicx::compartments::within_vs_between(values, 2, 3, 4);
    REQUIRE(ratios.size() == 3);
    CHECK(std::isinf(ratios[0]));
    // q = 2: a[0:2, 0:2] + a[2:4, 2:4] = 12 + 0 over a[0:2, 2:4] + a[2:4, 0:2] = 9 + 0
    CHECK(ratios[1] == 12.0 / 9.0);
    // q = 3: a[0:3, 0:3] + a[1:4, 1:4] = 21 + 11 over a[0:3, 1:4] + a[1:4, 0:3] = 16 + 15
    CHECK(ratios[2] == 32.0 / 31.0);
    // the square overload is the general one
    const std::vector<double> square{1, 2, 3, 4};
    CHECK(hicx::compartments::within_vs_between(square, 2) ==
          hicx::compartments::within_vs_between(square, 2, 2, 2));
}
