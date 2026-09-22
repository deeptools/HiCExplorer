// Bin merging.
//
// The expected matrices are the ones in the doctests of
// hicexplorer/reduceMatrix.py:48-144 and hicexplorer/hicMergeMatrixBins.py:
// 205-231, which is the only place where the Python states what the reduction
// is supposed to produce.

#include <doctest/doctest.h>

#include <cmath>
#include <string>
#include <vector>

#include "hicx/reduce_matrix.hpp"

using hicx::CsrMatrix;
using hicx::CutInterval;
using hicx::MatrixData;

namespace {

// A dense square matrix as a CSR, zeros dropped, in the given dtype.
CsrMatrix dense(const std::vector<std::vector<double>>& rows,
                const std::string& dtype) {
    const std::int64_t n = static_cast<std::int64_t>(rows.size());
    std::vector<std::int32_t> row_index;
    std::vector<std::int32_t> col_index;
    std::vector<double> values;
    for (std::int64_t r = 0; r < n; ++r) {
        for (std::int64_t c = 0; c < n; ++c) {
            const double value = rows[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
            if (value == 0.0) {
                continue;
            }
            row_index.push_back(static_cast<std::int32_t>(r));
            col_index.push_back(static_cast<std::int32_t>(c));
            values.push_back(value);
        }
    }
    return CsrMatrix::from_coo(n, n, row_index, col_index, std::move(values), dtype);
}

void check_dense(const CsrMatrix& matrix,
                 const std::vector<std::vector<double>>& expected) {
    REQUIRE(matrix.rows() == static_cast<std::int64_t>(expected.size()));
    for (std::int64_t r = 0; r < matrix.rows(); ++r) {
        for (std::int64_t c = 0; c < matrix.cols(); ++c) {
            CHECK(matrix.at(r, c) ==
                  expected[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)]);
        }
    }
}

// The 5x5 upper triangular int32 matrix of reduceMatrix.py:52.
CsrMatrix example_a() {
    return dense({{5, 5, 2, 2, 0},
                  {0, 5, 2, 2, 1},
                  {0, 0, 1, 1, 0},
                  {0, 0, 0, 1, 0},
                  {0, 0, 0, 0, 0}},
                 "int32");
}

}  // namespace

TEST_CASE("reduce_matrix returns the input when nothing is merged") {
    const CsrMatrix matrix = dense({{1, 0}, {0, 1}}, "int32");
    const CsrMatrix reduced = hicx::reduce_matrix(matrix, {{0}, {1}}, true, true);
    // reduceMatrix.py:154-155: len(bins_to_merge) == shape[0] short circuits.
    check_dense(reduced, {{1, 0}, {0, 1}});
}

TEST_CASE("reduce_matrix merges a two by two into a single bin") {
    const CsrMatrix matrix = dense({{1, 0}, {0, 1}}, "int32");
    const CsrMatrix reduced = hicx::reduce_matrix(matrix, {{0, 1}}, true, true);
    check_dense(reduced, {{2}});
}

TEST_CASE("reduce_matrix without the triangle selection") {
    const CsrMatrix matrix = example_a();
    check_dense(hicx::reduce_matrix(matrix, {{0, 1}, {2, 3}, {4}}, false, true),
                {{15, 8, 1}, {0, 3, 0}, {0, 0, 0}});
    check_dense(hicx::reduce_matrix(matrix, {{0, 1, 2}, {3, 4}}, false, true),
                {{20, 6}, {0, 1}});
}

TEST_CASE("reduce_matrix with the triangle selection symmetrises the result") {
    const CsrMatrix matrix = example_a();
    check_dense(hicx::reduce_matrix(matrix, {{0, 1}, {2, 3}, {4}}, true, true),
                {{15, 8, 1}, {8, 3, 0}, {1, 0, 0}});
}

TEST_CASE("reduce_matrix on a symmetric matrix, with and without the triangle") {
    const CsrMatrix matrix = dense({{2, 2, 1, 1, 1, 1},
                                    {2, 2, 1, 1, 1, 1},
                                    {1, 1, 1, 1, 1, 1},
                                    {1, 1, 1, 1, 1, 1},
                                    {1, 1, 1, 1, 1, 1},
                                    {1, 1, 1, 1, 1, 1}},
                                   "int32");
    check_dense(hicx::reduce_matrix(matrix, {{0, 1}, {2, 3}, {4, 5}}, false, true),
                {{8, 4, 4}, {4, 4, 4}, {4, 4, 4}});
    check_dense(hicx::reduce_matrix(matrix, {{0, 1}, {2, 3}, {4, 5}}, true, true),
                {{6, 4, 4}, {4, 3, 4}, {4, 4, 3}});
    check_dense(hicx::reduce_matrix(matrix, {{0, 1, 2}, {3, 4, 5}}, false, true),
                {{13, 9}, {9, 9}});
}

TEST_CASE("reduce_matrix drops the indices that appear in no group") {
    const CsrMatrix matrix = dense({{5, 1, 5, 1, 0},
                                    {1, 1, 2, 1, 1},
                                    {5, 2, 5, 1, 0},
                                    {1, 1, 1, 1, 0},
                                    {0, 1, 0, 0, 0}},
                                   "int32");
    check_dense(hicx::reduce_matrix(matrix, {{0, 2}, {1, 3}, {4}}, false, true),
                {{20, 5, 0}, {5, 4, 1}, {0, 1, 0}});
    // Rows 0 and 4 are not listed, so they and their columns disappear.
    check_dense(hicx::reduce_matrix(matrix, {{1, 2}, {3}}, false, true),
                {{10, 2}, {2, 1}});
}

TEST_CASE("reduce_matrix propagates NaN and keeps float values") {
    const double nan = std::nan("");
    const CsrMatrix matrix = dense({{0.1, 0.1, 0.2, 0.2, nan},
                                    {0.1, 0.1, 0.2, 0.2, 1.1},
                                    {0.2, 0.2, 0.2, 0.2, 0.0},
                                    {0.2, 0.2, 0.2, 0.1, 0.0},
                                    {nan, 1.1, 0.0, 0.0, 0.0}},
                                   "float64");
    const CsrMatrix reduced =
        hicx::reduce_matrix(matrix, {{0, 1}, {2, 3}, {4}}, false, true);
    CHECK(reduced.at(0, 0) == doctest::Approx(0.4));
    CHECK(reduced.at(0, 1) == doctest::Approx(0.8));
    CHECK(std::isnan(reduced.at(0, 2)));
    CHECK(reduced.at(1, 1) == doctest::Approx(0.7));
    CHECK(reduced.at(1, 2) == 0.0);
    CHECK(std::isnan(reduced.at(2, 0)));
}

TEST_CASE("reduce_matrix with pDiagonal false zeroes the main diagonal") {
    const CsrMatrix matrix = example_a();
    const CsrMatrix reduced =
        hicx::reduce_matrix(matrix, {{0, 1}, {2, 3}, {4}}, true, false);
    CHECK(reduced.at(0, 0) == 0.0);
    CHECK(reduced.at(1, 1) == 0.0);
    CHECK(reduced.at(2, 2) == 0.0);
    // The off diagonal is unaffected.
    CHECK(reduced.at(0, 1) == 8.0);
    CHECK(reduced.at(0, 2) == 1.0);
}

TEST_CASE("reduce_matrix casts the sums back to the input dtype") {
    // np.bincount sums in float64 and coo_matrix(dtype=int32) truncates.
    const CsrMatrix matrix = dense({{0.6, 0.6}, {0.0, 0.6}}, "int32");
    const CsrMatrix reduced = hicx::reduce_matrix(matrix, {{0, 1}}, true, true);
    CHECK(reduced.at(0, 0) == 1.0);
}

TEST_CASE("plan_bin_merge cuts on the bin count and on every chromosome change") {
    const std::vector<CutInterval> intervals{
        {"a", 0, 10, 0.5, ""},  {"a", 10, 20, 1.0, ""}, {"a", 20, 30, 1.0, ""},
        {"a", 30, 40, 0.1, ""}, {"b", 40, 50, 1.0, ""}};
    const hicx::BinMergePlan plan = hicx::plan_bin_merge(intervals, 2);
    REQUIRE(plan.intervals.size() == 3);
    CHECK(plan.intervals[0] == CutInterval{"a", 0, 20, 0.75, ""});
    CHECK(plan.intervals[1] == CutInterval{"a", 20, 40, 0.55, ""});
    CHECK(plan.intervals[2] == CutInterval{"b", 40, 50, 1.0, ""});
    CHECK(plan.bins_to_merge ==
          std::vector<std::vector<std::int64_t>>{{0, 1}, {2, 3}, {4}});
}

TEST_CASE("merge_bins reproduces the hicMergeMatrixBins doctest") {
    MatrixData input;
    input.cut_intervals = {{"a", 0, 10, 0.5, ""},  {"a", 10, 20, 1.0, ""},
                           {"a", 20, 30, 1.0, ""}, {"a", 30, 40, 0.1, ""},
                           {"b", 40, 50, 1.0, ""}};
    input.matrix = dense({{50, 10, 5, 3, 0},
                          {10, 60, 15, 5, 1},
                          {5, 15, 80, 7, 3},
                          {3, 5, 7, 90, 1},
                          {0, 1, 3, 1, 100}},
                         "int32");
    const MatrixData merged = hicx::merge_bins(input, 2);
    REQUIRE(merged.cut_intervals.size() == 3);
    CHECK(merged.cut_intervals[0] == CutInterval{"a", 0, 20, 0.75, ""});
    CHECK(merged.cut_intervals[1] == CutInterval{"a", 20, 40, 0.55, ""});
    CHECK(merged.cut_intervals[2] == CutInterval{"b", 40, 50, 1.0, ""});
    check_dense(merged.matrix, {{120, 28, 1}, {28, 177, 4}, {1, 4, 100}});
    CHECK(merged.nan_bins.empty());
}

TEST_CASE("merge_bins reports the empty columns as nan bins") {
    MatrixData input;
    input.cut_intervals = {{"a", 0, 10, 1.0, ""},
                           {"a", 10, 20, 1.0, ""},
                           {"a", 20, 30, 1.0, ""},
                           {"a", 30, 40, 1.0, ""}};
    input.matrix = dense({{1, 1, 0, 0}, {1, 1, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                         "int32");
    const MatrixData merged = hicx::merge_bins(input, 2);
    REQUIRE(merged.matrix.rows() == 2);
    CHECK(merged.nan_bins == std::vector<std::int64_t>{1});
}

TEST_CASE("plan_bin_merge_genome keeps the output shape when the input is missing bins") {
    // Two chromosomes of 45 bp at a 10 bp resolution: 'a' should tile into 5
    // bins (0,10,20,30,40-45) and 'b' into 4 (0,10,20,30-45)... actually a
    // clean 40 bp chromosome tiles into 4 bins and a 45 bp one into 5, with a
    // short last bin. Full input has every bin; the other input is missing
    // the tail of 'a' (as if no reads were observed there).
    const hicx::ChromosomeLengths sizes{{"a", 40}, {"b", 45}};
    const std::vector<CutInterval> full{
        {"a", 0, 10, 1.0, ""},  {"a", 10, 20, 1.0, ""}, {"a", 20, 30, 1.0, ""},
        {"a", 30, 40, 1.0, ""}, {"b", 0, 10, 1.0, ""},  {"b", 10, 20, 1.0, ""},
        {"b", 20, 30, 1.0, ""}, {"b", 30, 40, 1.0, ""}, {"b", 40, 45, 1.0, ""}};
    // Missing the last two bins of chromosome 'a'.
    const std::vector<CutInterval> partial{
        {"a", 0, 10, 1.0, ""},  {"a", 10, 20, 1.0, ""}, {"b", 0, 10, 1.0, ""},
        {"b", 10, 20, 1.0, ""}, {"b", 20, 30, 1.0, ""}, {"b", 30, 40, 1.0, ""},
        {"b", 40, 45, 1.0, ""}};

    const hicx::BinMergePlan plan_full = hicx::plan_bin_merge_genome(full, 2, sizes, 10);
    const hicx::BinMergePlan plan_partial = hicx::plan_bin_merge_genome(partial, 2, sizes, 10);

    // Same output shape from both: 2 groups on 'a' (0-20, 20-40) and 3 on 'b'
    // (0-20, 20-40, 40-45), a total of 5, whether or not the input actually
    // has data for every one of them.
    REQUIRE(plan_full.intervals.size() == 5);
    REQUIRE(plan_partial.intervals.size() == 5);
    for (std::size_t i = 0; i < plan_full.intervals.size(); ++i) {
        CHECK(plan_full.intervals[i].chrom == plan_partial.intervals[i].chrom);
        CHECK(plan_full.intervals[i].start == plan_partial.intervals[i].start);
        CHECK(plan_full.intervals[i].end == plan_partial.intervals[i].end);
    }
    CHECK(plan_full.intervals[4] == CutInterval{"b", 40, 45, 1.0, ""});
    CHECK(plan_partial.intervals[4] == CutInterval{"b", 40, 45, 1.0, ""});

    // The group covering the missing tail of 'a' (bins 20-40) is empty for
    // the partial input, not dropped: it still produces a zero row/column
    // rather than shrinking the layout, unlike plan_bin_merge.
    CHECK(plan_full.bins_to_merge[1] == std::vector<std::int64_t>{2, 3});
    CHECK(plan_partial.bins_to_merge[1].empty());
}

TEST_CASE("merge_bins_genome produces the same shape from two differently incomplete inputs") {
    // The concrete failure the project owner flagged: two matrices of the
    // same genome and resolution, one missing a trailing region the other
    // has data for, must merge into the same shape with the same
    // --chromosomeSizes and --numBins.
    const hicx::ChromosomeLengths sizes{{"a", 40}, {"b", 45}};

    MatrixData full;
    full.cut_intervals = {{"a", 0, 10, 1.0, ""},  {"a", 10, 20, 1.0, ""},
                          {"a", 20, 30, 1.0, ""}, {"a", 30, 40, 1.0, ""},
                          {"b", 0, 10, 1.0, ""},  {"b", 10, 20, 1.0, ""},
                          {"b", 20, 30, 1.0, ""}, {"b", 30, 40, 1.0, ""},
                          {"b", 40, 45, 1.0, ""}};
    full.matrix = dense({{5, 1, 0, 0, 0, 0, 0, 0, 0},
                        {1, 5, 1, 0, 0, 0, 0, 0, 0},
                        {0, 1, 5, 1, 0, 0, 0, 0, 0},
                        {0, 0, 1, 5, 1, 0, 0, 0, 0},
                        {0, 0, 0, 1, 5, 1, 0, 0, 0},
                        {0, 0, 0, 0, 1, 5, 1, 0, 0},
                        {0, 0, 0, 0, 0, 1, 5, 1, 0},
                        {0, 0, 0, 0, 0, 0, 1, 5, 1},
                        {0, 0, 0, 0, 0, 0, 0, 1, 5}},
                       "int32");

    MatrixData partial;
    partial.cut_intervals = {{"a", 0, 10, 1.0, ""}, {"a", 10, 20, 1.0, ""},
                             {"b", 0, 10, 1.0, ""}, {"b", 10, 20, 1.0, ""},
                             {"b", 20, 30, 1.0, ""}, {"b", 30, 40, 1.0, ""},
                             {"b", 40, 45, 1.0, ""}};
    partial.matrix = dense({{5, 1, 0, 0, 0, 0, 0},
                           {1, 5, 0, 0, 0, 0, 0},
                           {0, 0, 5, 1, 0, 0, 0},
                           {0, 0, 1, 5, 1, 0, 0},
                           {0, 0, 0, 1, 5, 1, 0},
                           {0, 0, 0, 0, 1, 5, 1},
                           {0, 0, 0, 0, 0, 1, 5}},
                          "int32");

    const MatrixData merged_full = hicx::merge_bins_genome(full, 2, sizes, 10);
    const MatrixData merged_partial = hicx::merge_bins_genome(partial, 2, sizes, 10);

    REQUIRE(merged_full.matrix.rows() == 5);
    REQUIRE(merged_partial.matrix.rows() == 5);
    REQUIRE(merged_full.cut_intervals.size() == 5);
    REQUIRE(merged_partial.cut_intervals.size() == 5);
    for (std::size_t i = 0; i < merged_full.cut_intervals.size(); ++i) {
        CHECK(merged_full.cut_intervals[i].chrom == merged_partial.cut_intervals[i].chrom);
        CHECK(merged_full.cut_intervals[i].start == merged_partial.cut_intervals[i].start);
        CHECK(merged_full.cut_intervals[i].end == merged_partial.cut_intervals[i].end);
    }

    // Values that both inputs actually have data for still sum correctly and
    // agree between the two runs: group 0 ('a' bins 0-20) is fully present in
    // both, groups on 'b' are fully present in both.
    check_dense(merged_full.matrix, {{11, 1, 0, 0, 0},
                                     {1, 11, 1, 0, 0},
                                     {0, 1, 11, 1, 0},
                                     {0, 0, 1, 11, 1},
                                     {0, 0, 0, 1, 5}});
    // Group 1 ('a' bins 20-40) has no data in the partial input, so it is an
    // all zero row and column rather than absent, and every other group's
    // value is unaffected by that.
    check_dense(merged_partial.matrix, {{11, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0},
                                        {0, 0, 11, 1, 0},
                                        {0, 0, 1, 11, 1},
                                        {0, 0, 0, 1, 5}});
    for (std::int64_t r = 0; r < 5; ++r) {
        for (std::int64_t c = r; c < 5; ++c) {
            if (r == 1 || c == 1) {
                continue;  // The group with no data in the partial input.
            }
            CHECK(merged_full.matrix.at(r, c) == merged_partial.matrix.at(r, c));
        }
    }
}

TEST_CASE("merge_bins_genome with no missing data matches merge_bins exactly") {
    // A regression check: when the input has every bin the genome layout
    // expects, the chromosome-sizes-driven path must be a pure relabelling of
    // the same groups plan_bin_merge would find on this input (its own bin
    // table already covers the genome with no numBins/2 drop in play), so the
    // merged matrix, at least, is byte identical to the pre-fix code path.
    const hicx::ChromosomeLengths sizes{{"a", 50}, {"b", 10}};
    MatrixData input;
    input.cut_intervals = {{"a", 0, 10, 0.5, ""},  {"a", 10, 20, 1.0, ""},
                           {"a", 20, 30, 1.0, ""}, {"a", 30, 40, 0.1, ""},
                           {"a", 40, 50, 1.0, ""}, {"b", 0, 10, 1.0, ""}};
    input.matrix = dense({{50, 10, 5, 3, 0, 0},
                          {10, 60, 15, 5, 1, 0},
                          {5, 15, 80, 7, 3, 0},
                          {3, 5, 7, 90, 1, 0},
                          {0, 1, 3, 1, 40, 2},
                          {0, 0, 0, 0, 2, 100}},
                         "int32");

    const MatrixData by_input = hicx::merge_bins(input, 2);
    const MatrixData by_genome = hicx::merge_bins_genome(input, 2, sizes, 10);

    REQUIRE(by_input.cut_intervals.size() == by_genome.cut_intervals.size());
    for (std::size_t i = 0; i < by_input.cut_intervals.size(); ++i) {
        CHECK(by_input.cut_intervals[i].chrom == by_genome.cut_intervals[i].chrom);
        CHECK(by_input.cut_intervals[i].start == by_genome.cut_intervals[i].start);
        CHECK(by_input.cut_intervals[i].end == by_genome.cut_intervals[i].end);
    }
    REQUIRE(by_input.matrix.rows() == by_genome.matrix.rows());
    for (std::int64_t r = 0; r < by_input.matrix.rows(); ++r) {
        for (std::int64_t c = 0; c < by_input.matrix.cols(); ++c) {
            CHECK(by_input.matrix.at(r, c) == by_genome.matrix.at(r, c));
        }
    }
    CHECK(by_input.nan_bins == by_genome.nan_bins);
}

TEST_CASE("plan_bin_merge drops a short leading chromosome but never the last") {
    const std::vector<CutInterval> intervals{
        {"a", 0, 10, 1.0, ""},  {"b", 0, 10, 1.0, ""},  {"b", 10, 20, 1.0, ""},
        {"b", 20, 30, 1.0, ""}, {"b", 30, 40, 1.0, ""}, {"c", 0, 10, 1.0, ""}};
    const hicx::BinMergePlan plan = hicx::plan_bin_merge(intervals, 4);
    // 'a' has one bin, below 4/2, so it is dropped with its row and column.
    // 'c' also has one bin but is the final group, which is appended outside
    // the loop and therefore never tested.
    REQUIRE(plan.intervals.size() == 2);
    CHECK(plan.intervals[0].chrom == "b");
    CHECK(plan.intervals[1].chrom == "c");
    CHECK(plan.bins_to_merge ==
          std::vector<std::vector<std::int64_t>>{{1, 2, 3, 4}, {5}});
}
