#include <doctest/doctest.h>

#include <vector>

#include "hicx/numpy_compat.hpp"
#include "hicx/sparse_matrix.hpp"

using hicx::CsrMatrix;

namespace {

// Upper triangle of
//   1 2 0
//   0 0 3
//   0 0 4
CsrMatrix upper_triangle() {
    const std::vector<std::int32_t> row{0, 0, 1, 2};
    const std::vector<std::int32_t> col{0, 1, 2, 2};
    std::vector<double> data{1.0, 2.0, 3.0, 4.0};
    return CsrMatrix::from_coo(3, 3, row, col, std::move(data), "float64");
}

}  // namespace

TEST_CASE("from_coo sorts the columns of every row and sums duplicates") {
    const std::vector<std::int32_t> row{1, 0, 1, 0};
    const std::vector<std::int32_t> col{2, 1, 2, 0};
    std::vector<double> data{1.5, 2.0, 0.5, 7.0};
    const CsrMatrix matrix = CsrMatrix::from_coo(2, 3, row, col, std::move(data),
                                                 "float64");
    CHECK(matrix.nnz() == 3);
    CHECK(matrix.at(0, 0) == 7.0);
    CHECK(matrix.at(0, 1) == 2.0);
    CHECK(matrix.at(1, 2) == 2.0);
    CHECK(matrix.indices() == std::vector<std::int32_t>{0, 1, 2});
}

TEST_CASE("materialize_full mirrors the upper triangle") {
    CsrMatrix matrix = upper_triangle();
    CHECK(matrix.lower_triangle_is_zero());
    CHECK_FALSE(matrix.has_lower_entries());
    matrix.materialize_full();
    CHECK(matrix.symmetry() == hicx::Symmetry::Full);

    CHECK(matrix.nnz() == 6);
    CHECK(matrix.at(0, 1) == 2.0);
    CHECK(matrix.at(1, 0) == 2.0);
    CHECK(matrix.at(1, 2) == 3.0);
    CHECK(matrix.at(2, 1) == 3.0);
    CHECK(matrix.at(0, 0) == 1.0);
    CHECK(matrix.at(2, 2) == 4.0);
    // Column indices stay sorted inside every row.
    for (std::int64_t r = 0; r < matrix.rows(); ++r) {
        for (std::int64_t k = matrix.indptr()[r] + 1; k < matrix.indptr()[r + 1]; ++k) {
            CHECK(matrix.indices()[k - 1] < matrix.indices()[k]);
        }
    }

    // Idempotent: a matrix that already has a filled lower triangle is left
    // alone, because tril(matrix, -1).sum() is no longer zero.
    const std::size_t nnz_before = matrix.nnz();
    matrix.materialize_full();
    CHECK(matrix.nnz() == nnz_before);
}

TEST_CASE("symmetrize_in_place keeps the upper triangle and answers the same") {
    CsrMatrix lazy = upper_triangle();
    lazy.symmetrize_in_place();
    CsrMatrix eager = upper_triangle();
    eager.materialize_full();

    CHECK(lazy.symmetry() == hicx::Symmetry::UpperTriangle);
    // Half the entries in memory, the same matrix.
    CHECK(lazy.stored_nnz() == 4);
    CHECK(eager.stored_nnz() == 6);
    CHECK(lazy.nnz() == eager.nnz());
    CHECK(lazy.sum().as_double() == eager.sum().as_double());
    CHECK(lazy.diagonal_sum().as_double() == eager.diagonal_sum().as_double());
    CHECK(lazy.data_min().as_double() == eager.data_min().as_double());
    CHECK(lazy.data_max().as_double() == eager.data_max().as_double());
    for (std::int64_t r = 0; r < 3; ++r) {
        for (std::int64_t c = 0; c < 3; ++c) {
            CHECK(lazy.at(r, c) == eager.at(r, c));
        }
    }

    // Expanding afterwards gives the same matrix as expanding directly.
    lazy.materialize_full();
    CHECK(lazy.symmetry() == hicx::Symmetry::Full);
    CHECK(lazy.indptr() == eager.indptr());
    CHECK(lazy.indices() == eager.indices());
    CHECK(lazy.data() == eager.data());
}

TEST_CASE("adding the transpose drops entries that cancel to zero") {
    // scipy's csr_binop_csr only stores non zero results, so an explicit zero
    // in the input disappears from the symmetric matrix.
    const std::vector<std::int32_t> row{0, 0, 1};
    const std::vector<std::int32_t> col{0, 1, 1};
    std::vector<double> data{1.0, 0.0, 2.0};
    CsrMatrix matrix = CsrMatrix::from_coo(2, 2, row, col, std::move(data), "float64");
    CHECK(matrix.stored_nnz() == 3);  // the explicit zero is still stored
    CsrMatrix lazy = matrix;
    matrix.materialize_full();
    CHECK(matrix.nnz() == 2);
    lazy.symmetrize_in_place();
    CHECK(lazy.nnz() == 2);
}

TEST_CASE("sum, diagonal and extrema follow scipy") {
    CsrMatrix matrix = upper_triangle();
    matrix.materialize_full();

    const hicx::Scalar total = matrix.sum();
    const hicx::Scalar diagonal = matrix.diagonal_sum();
    CHECK_FALSE(total.is_integer());
    CHECK(total.as_double() == 15.0);
    CHECK(diagonal.as_double() == 5.0);

    // hicInfo's "Sum of matrix"
    const double sum_elements = ((total.as_double() - diagonal.as_double()) / 2.0) +
                                diagonal.as_double();
    CHECK(hicx::npy::float_repr(sum_elements) == "10.0");

    CHECK(matrix.data_min().as_double() == 1.0);
    CHECK(matrix.data_max().as_double() == 4.0);
}

TEST_CASE("integer matrices keep integer semantics") {
    const std::vector<std::int32_t> row{0, 0, 1};
    const std::vector<std::int32_t> col{0, 1, 1};
    std::vector<double> data{3.0, 5.0, 7.0};
    CsrMatrix matrix = CsrMatrix::from_coo(2, 2, row, col, std::move(data), "int32");
    matrix.materialize_full();

    CHECK(matrix.integral_dtype());
    CHECK(matrix.sum().is_integer());
    CHECK(matrix.sum().integer_value == 20);
    CHECK(matrix.diagonal_sum().integer_value == 10);
    CHECK(matrix.data_max().is_integer());
    CHECK(matrix.data_max().integer_value == 7);
}

TEST_CASE("eliminate_zeros removes stored zeros") {
    std::vector<std::int64_t> indptr{0, 2, 3};
    std::vector<std::int32_t> indices{0, 1, 1};
    std::vector<double> data{0.0, 4.0, 0.0};
    CsrMatrix matrix(2, 2, std::move(indptr), std::move(indices), std::move(data),
                     "float64");
    CHECK(matrix.nnz() == 3);
    matrix.eliminate_zeros();
    CHECK(matrix.nnz() == 1);
    CHECK(matrix.at(0, 1) == 4.0);
    CHECK(matrix.indptr() == std::vector<std::int64_t>{0, 1, 1});
}
