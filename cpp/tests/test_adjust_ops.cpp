// Bin selection, deletion and block zeroing, plus the two additions to the bin
// merging layer that hicMergeMatrixBins and hicMergeTADbins need.
//
// The expected matrices for running_window come from the doctests of
// hicexplorer/hicMergeMatrixBins.py:100-135, which is the only place where the
// Python states what the running window is supposed to produce. Everything
// else is checked against what scipy's fancy indexing and hiCMatrix.maskBins
// do, spelled out per cell.

#include <doctest/doctest.h>

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "hicx/adjust_ops.hpp"
#include "hicx/matrix_ops.hpp"
#include "hicx/reduce_matrix.hpp"

using hicx::CsrMatrix;
using hicx::CutInterval;
using hicx::InterIntra;
using hicx::MatrixData;
using hicx::Symmetry;

namespace {

CsrMatrix dense(const std::vector<std::vector<double>>& rows, const std::string& dtype) {
    const std::int64_t n = static_cast<std::int64_t>(rows.size());
    std::vector<std::int32_t> row_index;
    std::vector<std::int32_t> col_index;
    std::vector<double> values;
    for (std::int64_t r = 0; r < n; ++r) {
        for (std::int64_t c = 0; c < n; ++c) {
            const double value =
                rows[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
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

// The upper triangle of a symmetric matrix, marked as such, which is the state
// every matrix is in after hiCMatrix's fillLowerTriangle.
CsrMatrix upper(const std::vector<std::vector<double>>& rows, const std::string& dtype) {
    const std::int64_t n = static_cast<std::int64_t>(rows.size());
    std::vector<std::vector<double>> triangle(rows);
    for (std::int64_t r = 0; r < n; ++r) {
        for (std::int64_t c = 0; c < r; ++c) {
            triangle[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] = 0.0;
        }
    }
    CsrMatrix matrix = dense(triangle, dtype);
    matrix.set_symmetry(Symmetry::UpperTriangle);
    return matrix;
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

std::vector<CutInterval> four_bins() {
    return {CutInterval{"a", 0, 10, 0.5, ""}, CutInterval{"a", 10, 20, 1.0, ""},
            CutInterval{"a", 20, 30, 1.0, ""}, CutInterval{"b", 40, 50, 1.0, ""}};
}

}  // namespace

TEST_CASE("select_bins keeps the upper triangle for an increasing selection") {
    const CsrMatrix matrix = upper({{1, 2, 3, 4},
                                    {2, 5, 6, 7},
                                    {3, 6, 8, 9},
                                    {4, 7, 9, 10}},
                                   "int32");
    const CsrMatrix selected = hicx::select_bins(matrix, {0, 2, 3});
    CHECK(selected.symmetry() == Symmetry::UpperTriangle);
    check_dense(selected, {{1, 3, 4}, {3, 8, 9}, {4, 9, 10}});
}

TEST_CASE("select_bins reorders rows and columns without materialising the mirror") {
    const CsrMatrix matrix = upper({{1, 2, 3}, {2, 4, 5}, {3, 5, 6}}, "int32");
    const CsrMatrix reordered = hicx::select_bins(matrix, {2, 0, 1});
    // P A P^T is symmetric, so the result is still only half stored.
    CHECK(reordered.symmetry() == Symmetry::UpperTriangle);
    check_dense(reordered, {{6, 3, 5}, {3, 1, 2}, {5, 2, 4}});
}

TEST_CASE("select_bins duplicates a bin that is selected twice") {
    // scipy's matrix[[0, 1, 1], :][:, [0, 1, 1]] repeats the row and the
    // column rather than merging them, which is what an overlapping --regions
    // BED produces in hicAdjustMatrix.
    const CsrMatrix matrix = upper({{1, 2}, {2, 3}}, "int32");
    const CsrMatrix selected = hicx::select_bins(matrix, {0, 1, 1});
    check_dense(selected, {{1, 2, 2}, {2, 3, 3}, {2, 3, 3}});
}

TEST_CASE("select_bins on a matrix that stores both triangles gathers rows") {
    CsrMatrix matrix = dense({{1, 2, 0}, {0, 3, 4}, {5, 0, 6}}, "float64");
    CHECK(matrix.symmetry() == Symmetry::Full);
    const CsrMatrix selected = hicx::select_bins(matrix, {2, 0});
    CHECK(selected.symmetry() == Symmetry::Full);
    check_dense(selected, {{6, 5}, {0, 1}});
}

TEST_CASE("reorder_bins carries the bin table and the nan bins but not the factors") {
    MatrixData data;
    data.matrix = upper({{1, 2, 3, 4},
                         {2, 5, 6, 7},
                         {3, 6, 8, 9},
                         {4, 7, 9, 10}},
                        "int32");
    data.cut_intervals = four_bins();
    data.nan_bins = {1, 3};
    data.correction_factors = std::vector<double>{10.0, 20.0, 30.0, 40.0};

    hicx::reorder_bins(data, {3, 1});
    check_dense(data.matrix, {{10, 7}, {7, 5}});
    REQUIRE(data.cut_intervals.size() == 2);
    CHECK(data.cut_intervals[0].chrom == "b");
    CHECK(data.cut_intervals[1].start == 10);
    CHECK(data.nan_bins == std::vector<std::int64_t>{0, 1});
    // HiCMatrix.py:723-746 never touches correction_factors, cpp/PLAN.md 2.7
    // quirk 3. Reproduced, not fixed: the vector still has four entries.
    REQUIRE(data.correction_factors.has_value());
    CHECK(data.correction_factors->size() == 4);
}

TEST_CASE("delete_bins drops the requested bins together with the nan bins") {
    MatrixData data;
    data.matrix = upper({{1, 2, 3, 4},
                         {2, 5, 6, 7},
                         {3, 6, 8, 9},
                         {4, 7, 9, 10}},
                        "int32");
    data.cut_intervals = four_bins();
    data.nan_bins = {0};
    data.correction_factors = std::vector<double>{10.0, 20.0, 30.0, 40.0};

    hicx::delete_bins(data, {2});
    // bin 0 goes because it is a nan bin, bin 2 because it was asked for.
    check_dense(data.matrix, {{5, 7}, {7, 10}});
    REQUIRE(data.cut_intervals.size() == 2);
    CHECK(data.cut_intervals[0].start == 10);
    CHECK(data.cut_intervals[1].chrom == "b");
    CHECK(data.nan_bins.empty());
    REQUIRE(data.correction_factors.has_value());
    CHECK(*data.correction_factors == std::vector<double>{20.0, 40.0});
}

TEST_CASE("delete_bins with an empty list is a no operation") {
    MatrixData data;
    data.matrix = upper({{1, 2}, {2, 3}}, "int32");
    data.cut_intervals = {CutInterval{"a", 0, 10, 1.0, ""},
                          CutInterval{"a", 10, 20, 1.0, ""}};
    data.nan_bins = {0};
    hicx::delete_bins(data, {});
    // maskBins returns before it looks at the nan bins, so bin 0 survives.
    CHECK(data.matrix.rows() == 2);
    CHECK(data.nan_bins == std::vector<std::int64_t>{0});
}

TEST_CASE("empty_column_bins sees the mirror of an upper triangle") {
    const CsrMatrix matrix = upper({{0, 0, 3}, {0, 0, 0}, {3, 0, 0}}, "int32");
    // Column 0 is empty in the stored triangle but not in the matrix it
    // represents, because the entry at (0, 2) mirrors to (2, 0).
    CHECK(hicx::empty_column_bins(matrix) == std::vector<std::int64_t>{1});
}

TEST_CASE("zero_inter_or_intra blanks the blocks the Python slices assign to") {
    const std::vector<std::pair<std::string, hicx::BinRange>> boundaries = {
        {"a", hicx::BinRange{0, 2}}, {"b", hicx::BinRange{2, 4}}};

    SUBCASE("inter") {
        CsrMatrix matrix = upper({{1, 2, 3, 4},
                                  {2, 5, 6, 7},
                                  {3, 6, 8, 9},
                                  {4, 7, 9, 10}},
                                 "int32");
        hicx::zero_inter_or_intra(matrix, boundaries, InterIntra::Inter);
        check_dense(matrix, {{1, 2, 0, 0}, {2, 5, 0, 0}, {0, 0, 8, 9}, {0, 0, 9, 10}});
    }
    SUBCASE("intra") {
        CsrMatrix matrix = upper({{1, 2, 3, 4},
                                  {2, 5, 6, 7},
                                  {3, 6, 8, 9},
                                  {4, 7, 9, 10}},
                                 "int32");
        hicx::zero_inter_or_intra(matrix, boundaries, InterIntra::Intra);
        check_dense(matrix, {{0, 0, 3, 4}, {0, 0, 6, 7}, {3, 6, 0, 0}, {4, 7, 0, 0}});
    }
}

TEST_CASE("running_window reproduces the hicMergeMatrixBins doctests") {
    SUBCASE("two bins, window of three") {
        // hicMergeMatrixBins.py:112-121: a 2x2 matrix of ones becomes 3s.
        const CsrMatrix matrix = upper({{1, 1}, {1, 1}}, "int64");
        const CsrMatrix windowed = hicx::running_window(matrix, 3);
        check_dense(windowed, {{3, 3}, {3, 3}});
    }
    SUBCASE("four bins, window of three") {
        // hicMergeMatrixBins.py:123-135.
        const CsrMatrix matrix = upper({{1, 1, 1, 1},
                                        {1, 1, 1, 1},
                                        {1, 1, 1, 1},
                                        {1, 1, 1, 1}},
                                       "int64");
        const CsrMatrix windowed = hicx::running_window(matrix, 3);
        check_dense(windowed,
                    {{3, 5, 6, 4}, {5, 6, 8, 6}, {6, 8, 6, 5}, {4, 6, 5, 3}});
    }
}

TEST_CASE("running_window rejects an even window") {
    const CsrMatrix matrix = upper({{1, 1}, {1, 1}}, "int64");
    CHECK_THROWS_AS(hicx::running_window(matrix, 4), std::invalid_argument);
}

TEST_CASE("running_window ignores chromosome borders") {
    // The window is applied to raw bin indices, so a cell that spans the
    // border between two chromosomes picks up counts from both. Pinned by
    // test_hicMergeMatrixBins.py and reproduced rather than fixed.
    const CsrMatrix matrix = upper({{1, 0, 0}, {0, 0, 0}, {0, 0, 1}}, "int64");
    const CsrMatrix windowed = hicx::running_window(matrix, 3);
    CHECK(windowed.at(0, 1) == 1.0);
    CHECK(windowed.at(1, 2) == 1.0);
}

TEST_CASE("plan_tad_merge cuts at every boundary and at every chromosome change") {
    const std::vector<CutInterval> intervals = {
        CutInterval{"a", 0, 10, 1.0, ""},  CutInterval{"a", 10, 20, 2.0, ""},
        CutInterval{"a", 20, 30, 3.0, ""}, CutInterval{"a", 30, 40, 4.0, ""},
        CutInterval{"b", 40, 50, 5.0, ""}};
    // A boundary at bin 2 cuts between bin 1 and bin 2; the chromosome change
    // at bin 4 cuts again. A boundary at bin 0 is ignored, because the count
    // is zero at the start of a group.
    const hicx::BinMergePlan plan = hicx::plan_tad_merge(intervals, {0, 2});
    REQUIRE(plan.bins_to_merge.size() == 3);
    CHECK(plan.bins_to_merge[0] == std::vector<std::int64_t>{0, 1});
    CHECK(plan.bins_to_merge[1] == std::vector<std::int64_t>{2, 3});
    CHECK(plan.bins_to_merge[2] == std::vector<std::int64_t>{4});
    REQUIRE(plan.intervals.size() == 3);
    CHECK(plan.intervals[0].chrom == "a");
    CHECK(plan.intervals[0].start == 0);
    CHECK(plan.intervals[0].end == 20);
    CHECK(plan.intervals[0].extra == doctest::Approx(1.5));
    CHECK(plan.intervals[1].start == 20);
    CHECK(plan.intervals[1].end == 40);
    CHECK(plan.intervals[2].chrom == "b");
    CHECK(plan.intervals[2].extra == doctest::Approx(5.0));
}

TEST_CASE("plan_tad_merge keeps a group shorter than half a bin") {
    // The difference from plan_bin_merge, which drops a trailing group of
    // fewer than num_bins/2 bins. A TAD is a TAD whatever its size.
    const std::vector<CutInterval> intervals = {CutInterval{"a", 0, 10, 1.0, ""},
                                                CutInterval{"a", 10, 20, 1.0, ""},
                                                CutInterval{"a", 20, 30, 1.0, ""}};
    const hicx::BinMergePlan plan = hicx::plan_tad_merge(intervals, {2});
    REQUIRE(plan.bins_to_merge.size() == 2);
    CHECK(plan.bins_to_merge[1] == std::vector<std::int64_t>{2});
}

TEST_CASE("plan_tad_merge on an empty bin table produces nothing to merge") {
    const hicx::BinMergePlan plan = hicx::plan_tad_merge({}, {});
    CHECK(plan.bins_to_merge.empty());
    CHECK(plan.intervals.empty());
}

TEST_CASE("reduce_matrix with diagonal loses the within group sum from the total") {
    // The defect hicMergeTADbins ships: R + R.T - diag(R) subtracts the whole
    // within-TAD block sum rather than the original main diagonal, so the
    // symmetric total is not conserved while the upper triangle is.
    const CsrMatrix matrix = upper({{1, 2, 0}, {2, 4, 0}, {0, 0, 8}}, "int64");
    const CsrMatrix reduced = hicx::reduce_matrix(matrix, {{0, 1}, {2}}, true, true);
    // The upper triangle of the input sums to 1 + 2 + 4 + 8 = 15, and so does
    // the upper triangle of the result.
    CHECK(reduced.at(0, 0) == 7.0);
    CHECK(reduced.at(1, 1) == 8.0);
    CHECK(reduced.at(0, 1) == 0.0);
    // The symmetric input sums to 1 + 2 + 2 + 4 + 8 = 17, the symmetric result
    // to 7 + 8 = 15. Two counts vanish, which is the same arithmetic that
    // costs 22 percent on Li_et_al_2015.h5.
    CHECK(reduced.at(0, 0) + reduced.at(1, 1) == 15.0);
}
