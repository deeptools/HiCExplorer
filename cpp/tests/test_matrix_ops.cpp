// Tests for the matrix arithmetic and the mask semantics that hicSumMatrices
// and hicCompareMatrices are built on.
//
// The small matrices check the scipy and numpy rules in isolation; the cases
// on the real matrices of hicexplorer/test/test_data check the two numbers the
// Python characterization tests pin, so a regression in the reduction order or
// in the mask is caught here and not only by the equivalence harness.

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "hicx/matrix_data.hpp"
#include "hicx/matrix_ops.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/tool_matrix.hpp"

using hicx::CsrMatrix;

namespace {

std::string data_path(const std::string& relative) {
    return std::string(HICX_TEST_DATA_DIR) + "/" + relative;
}

// Upper triangle of the symmetric matrix
//   1 2 0
//   2 0 3
//   0 3 4
CsrMatrix upper_triangle(std::string dtype = "float64") {
    const std::vector<std::int32_t> row{0, 0, 1, 2};
    const std::vector<std::int32_t> col{0, 1, 2, 2};
    std::vector<double> data{1.0, 2.0, 3.0, 4.0};
    CsrMatrix matrix = CsrMatrix::from_coo(3, 3, row, col, std::move(data),
                                           std::move(dtype));
    matrix.symmetrize_in_place();
    return matrix;
}

std::vector<hicx::CutInterval> three_bins(const std::string& first,
                                          const std::string& second,
                                          const std::string& third) {
    return {hicx::CutInterval{first, 0, 10, 1.0, {}},
            hicx::CutInterval{second, 10, 20, 1.0, {}},
            hicx::CutInterval{third, 20, 30, 1.0, {}}};
}

}  // namespace

TEST_CASE("promote_dtype follows numpy.promote_types") {
    CHECK(hicx::promote_dtype("int32", "int32") == "int32");
    CHECK(hicx::promote_dtype("int32", "int64") == "int64");
    CHECK(hicx::promote_dtype("float32", "float32") == "float32");
    CHECK(hicx::promote_dtype("float32", "float64") == "float64");
    // numpy widens int32 or int64 mixed with float32 all the way to float64.
    CHECK(hicx::promote_dtype("int32", "float32") == "float64");
    CHECK(hicx::promote_dtype("int64", "float32") == "float64");
    CHECK(hicx::promote_dtype("int64", "float64") == "float64");
}

TEST_CASE("adding two upper triangles keeps the upper triangle storage") {
    const CsrMatrix a = upper_triangle();
    const CsrMatrix b = upper_triangle();
    const CsrMatrix sum = hicx::add(a, b);
    CHECK(sum.symmetry() == hicx::Symmetry::UpperTriangle);
    CHECK(sum.at(0, 0) == 2.0);
    CHECK(sum.at(0, 1) == 4.0);
    CHECK(sum.at(1, 0) == 4.0);
    CHECK(sum.at(2, 2) == 8.0);
    CHECK(sum.stored_nnz() == 4);
}

TEST_CASE("a binary operation drops an exactly zero result, like csr_binop_csr") {
    const CsrMatrix a = upper_triangle();
    const CsrMatrix b = upper_triangle();
    const CsrMatrix difference = hicx::subtract(a, b);
    CHECK(difference.stored_nnz() == 0);
}

TEST_CASE("elementwise multiply keeps only the common support") {
    const std::vector<std::int32_t> row_a{0, 0};
    const std::vector<std::int32_t> col_a{0, 1};
    const std::vector<std::int32_t> row_b{0, 1};
    const std::vector<std::int32_t> col_b{1, 2};
    CsrMatrix a = CsrMatrix::from_coo(3, 3, row_a, col_a, {2.0, 3.0}, "float64");
    CsrMatrix b = CsrMatrix::from_coo(3, 3, row_b, col_b, {5.0, 7.0}, "float64");
    a.symmetrize_in_place();
    b.symmetrize_in_place();
    const CsrMatrix product = hicx::multiply_elementwise(a, b);
    CHECK(product.stored_nnz() == 1);
    CHECK(product.at(0, 1) == 15.0);
}

TEST_CASE("an integer sum stays integer and a float32 pair stays float32") {
    const CsrMatrix a = upper_triangle("int32");
    const CsrMatrix b = upper_triangle("int32");
    CHECK(hicx::add(a, b).dtype() == "int32");

    const CsrMatrix c = upper_triangle("float32");
    CHECK(hicx::add(c, upper_triangle("float64")).dtype() == "float64");
    CHECK(hicx::add(c, c).dtype() == "float32");
}

TEST_CASE("data_sum walks the symmetric data array, not the stored triangle") {
    const CsrMatrix matrix = upper_triangle();
    // The symmetric matrix holds 1, 2, 2, 3, 3, 4.
    CHECK(hicx::data_sum(matrix).as_double() == 15.0);
    // matrix.sum() is the same total here, but it is a different reduction and
    // a different function; the tools need the data array one.
    CHECK(matrix.sum().as_double() == 15.0);
}

TEST_CASE("data_sum of an integer matrix is an exact integer") {
    const CsrMatrix matrix = upper_triangle("int32");
    const hicx::Scalar total = hicx::data_sum(matrix);
    CHECK(total.is_integer());
    CHECK(total.integer_value == 15);
}

TEST_CASE("reciprocal and log2 follow numpy's dtype rules") {
    CsrMatrix matrix = upper_triangle("int32");
    hicx::reciprocal_data_in_place(matrix);
    CHECK(matrix.dtype() == "float64");
    CHECK(matrix.at(0, 1) == 0.5);

    CsrMatrix single = upper_triangle("float32");
    hicx::reciprocal_data_in_place(single);
    // A Python float is a weak scalar, so float32 stays float32.
    CHECK(single.dtype() == "float32");

    CsrMatrix values = upper_triangle();
    hicx::log2_data_in_place(values);
    CHECK(values.at(0, 1) == 1.0);
    CHECK(values.at(2, 2) == 2.0);
    CHECK(values.at(0, 0) == 0.0);
}

TEST_CASE("chrom_bin_boundaries reproduces the OrderedDict of hicmatrix") {
    const auto boundaries = hicx::chrom_bin_boundaries(three_bins("chr1", "chr1", "chr2"));
    REQUIRE(boundaries.size() == 2);
    CHECK(boundaries[0].first == "chr1");
    CHECK(boundaries[0].second.first == 0);
    CHECK(boundaries[0].second.last == 2);
    CHECK(boundaries[1].first == "chr2");
    CHECK(boundaries[1].second.first == 2);
    CHECK(boundaries[1].second.last == 3);

    // Order is part of the value: the tools compare two of these for equality.
    CHECK(hicx::chrom_bin_boundaries(three_bins("chr1", "chr2", "chr2")) !=
          hicx::chrom_bin_boundaries(three_bins("chr2", "chr2", "chr1")));

    // A chromosome that reappears overwrites its entry and keeps its position,
    // which is what assigning to an existing OrderedDict key does.
    const auto repeated = hicx::chrom_bin_boundaries(three_bins("chrA", "chrB", "chrA"));
    REQUIRE(repeated.size() == 2);
    CHECK(repeated[0].first == "chrA");
    CHECK(repeated[0].second.first == 2);
    CHECK(repeated[0].second.last == 3);
}

TEST_CASE("chrom_bin_boundaries agrees with BinTable on a real bin table") {
    const hicx::ToolMatrix matrix =
        hicx::ToolMatrix::load(data_path("small_test_matrix.cool"));
    const hicx::BinTable table(matrix.cut_intervals());
    CHECK(matrix.boundaries() == table.chrom_bin_boundaries());
    CHECK(matrix.boundaries().size() == 15);
}

TEST_CASE("mask_and_restore_bins deletes rows and columns and upcasts to float64") {
    hicx::MatrixData data;
    data.matrix = upper_triangle("int32");
    data.cut_intervals = three_bins("chr1", "chr1", "chr1");

    hicx::mask_and_restore_bins(data, {1});
    CHECK(data.matrix.rows() == 3);
    CHECK(data.matrix.dtype() == "float64");
    CHECK(data.matrix.symmetry() == hicx::Symmetry::UpperTriangle);
    // (0, 1) and (1, 2) are gone, (0, 0) and (2, 2) survive.
    CHECK(data.matrix.at(0, 0) == 1.0);
    CHECK(data.matrix.at(0, 1) == 0.0);
    CHECK(data.matrix.at(1, 2) == 0.0);
    CHECK(data.matrix.at(2, 2) == 4.0);
    CHECK(data.nan_bins == std::vector<std::int64_t>{1});
}

TEST_CASE("an empty mask changes nothing, not even the dtype") {
    hicx::MatrixData data;
    data.matrix = upper_triangle("int32");
    data.cut_intervals = three_bins("chr1", "chr1", "chr1");
    hicx::mask_and_restore_bins(data, {});
    CHECK(data.matrix.dtype() == "int32");
    CHECK(data.matrix.stored_nnz() == 4);
}

TEST_CASE("mask_and_restore_bins folds in the NaN bins the matrix already has") {
    hicx::MatrixData data;
    data.matrix = upper_triangle();
    data.cut_intervals = three_bins("chr1", "chr1", "chr1");
    data.nan_bins = {0};
    data.correction_factors = std::vector<double>{1.0, 2.0, 3.0};

    hicx::mask_and_restore_bins(data, {2});
    CHECK(data.nan_bins == std::vector<std::int64_t>{0, 2});
    CHECK(data.matrix.stored_nnz() == 0);
    REQUIRE(data.correction_factors.has_value());
    CHECK(std::isnan((*data.correction_factors)[0]));
    CHECK((*data.correction_factors)[1] == 2.0);
    CHECK(std::isnan((*data.correction_factors)[2]));
}

TEST_CASE("data_sum reproduces the normalisers hicCompareMatrices pins") {
    // hicexplorer/test/general/test_hicCompareMatrices.py:245-248 asserts
    // these two values of matrix.data.sum() at full float64 precision.
    const hicx::ToolMatrix untreated = hicx::ToolMatrix::load(
        data_path("hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1.h5"));
    const hicx::ToolMatrix auxin = hicx::ToolMatrix::load(
        data_path("hicDifferentialTAD/GSM2644947_Auxin2days-R1.100000_chr1.h5"));
    CHECK(hicx::data_sum(untreated.matrix()).as_double() == 1514.2970371802683);
    CHECK(hicx::data_sum(auxin.matrix()).as_double() == 1385.1134232101283);
    // The cool twin of the same matrix must give the same number.
    const hicx::ToolMatrix untreated_cool = hicx::ToolMatrix::load(
        data_path("hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1.cool"));
    CHECK(hicx::data_sum(untreated_cool.matrix()).as_double() == 1514.2970371802683);
}

TEST_CASE("the mask step loses the 5,142 counts the Python test pins") {
    // test_hicSumMatrices.py:116-136. The plain sparse sum has 3,157,763
    // entries and the file the tool writes has 3,152,621.
    hicx::ToolMatrix a = hicx::ToolMatrix::load(
        data_path("hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1.h5"));
    const hicx::ToolMatrix b = hicx::ToolMatrix::load(
        data_path("hicDifferentialTAD/GSM2644947_Auxin2days-R1.100000_chr1.h5"));
    CHECK(a.nan_bins().size() == 83);
    CHECK(b.nan_bins().size() == 81);

    a.matrix() = hicx::add(a.matrix(), b.matrix());
    CHECK(a.matrix().nnz() == 3157763);

    std::vector<std::int64_t> mask = a.nan_bins();
    for (const std::int64_t bin : b.nan_bins()) {
        if (std::find(mask.begin(), mask.end(), bin) == mask.end()) {
            mask.push_back(bin);
        }
    }
    std::sort(mask.begin(), mask.end());
    CHECK(mask.size() == 83);
    hicx::mask_and_restore_bins(a.data(), mask);
    CHECK(a.matrix().nnz() == 3152621);
}

TEST_CASE("loading through ToolMatrix swaps the two optional vectors") {
    // cpp/PLAN.md 2.7 quirk 1: hiCMatrix.__init__ unpacks the loader tuple with
    // correction_factors and distance_counts the wrong way round.
    // Li_et_al_2015.h5 carries /correction_factors and no /distance_counts, so
    // after the swap the factors sit in distance_counts.
    const hicx::ToolMatrix matrix = hicx::ToolMatrix::load(data_path("Li_et_al_2015.h5"));
    CHECK_FALSE(matrix.data().correction_factors.has_value());
    REQUIRE(matrix.data().distance_counts.has_value());
    CHECK(matrix.data().distance_counts->size() == 11104);
}
