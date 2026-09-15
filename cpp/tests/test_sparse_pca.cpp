// hicPCA --eigenSolver lanczos (hicx/sparse_pca.hpp): the implicit covariance
// product against the dense covariance, the Krylov eigenvectors against the
// dense symmetric solver, and byte identity across thread counts.

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

#include "hicx/adjust_ops.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/sparse_pca.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/transform_ops.hpp"

namespace {

std::string data_path(const std::string& relative) {
    return std::string(HICX_TEST_DATA_DIR) + "/" + relative;
}

hicx::CsrMatrix chrx_obs_exp() {
    hicx::ToolMatrix hic = hicx::ToolMatrix::load(data_path("small_test_matrix_50kb_res.h5"));
    hicx::keep_only_chromosomes(hic.data(), {"chrX"});
    hicx::CsrMatrix block = hic.matrix();
    block.materialize_full();
    hicx::obs_exp_lieberman_in_place(block, block.rows(), 1);
    return block;
}

}  // namespace

TEST_CASE("the implicit covariance product equals the dense covariance times a vector") {
    const hicx::CsrMatrix block = chrx_obs_exp();
    const std::int64_t n = block.rows();
    const hicx::DenseSymmetric dense = hicx::covariance_of_symmetric(block, 1);
    const hicx::CovarianceOperator op(block, 1);
    REQUIRE(op.masked_rows() == 0);
    std::vector<double> x(static_cast<std::size_t>(n));
    for (std::int64_t i = 0; i < n; ++i) {
        x[static_cast<std::size_t>(i)] = std::cos(0.013 * static_cast<double>(i * i));
    }
    std::vector<double> y(static_cast<std::size_t>(n));
    op.apply(x.data(), y.data());
    double scale = 0.0;
    double worst = 0.0;
    for (std::int64_t i = 0; i < n; ++i) {
        double want = 0.0;
        for (std::int64_t j = 0; j < n; ++j) {
            want += dense.at(i, j) * x[static_cast<std::size_t>(j)];
        }
        scale = std::max(scale, std::fabs(want));
        worst = std::max(worst, std::fabs(y[static_cast<std::size_t>(i)] - want));
    }
    CHECK(worst <= 1e-11 * scale);
}

TEST_CASE("lanczos eigenvectors agree with dsyevr on a real chromosome block") {
    const hicx::CsrMatrix block = chrx_obs_exp();
    hicx::pin_blas_to_one_thread();
    hicx::DenseSymmetric dense = hicx::covariance_of_symmetric(block, 1);
    hicx::zero_non_finite_in_place(dense, 1);
    const std::vector<int> which{1, 2, 3};
    const hicx::EigenResult reference =
        hicx::leading_eigenvectors(dense, which, hicx::EigenSolver::Dsyevr);
    hicx::LanczosDiagnostics diagnostics;
    const hicx::EigenResult got =
        hicx::covariance_eigenvectors_lanczos(block, which, 1, &diagnostics);
    CHECK_FALSE(diagnostics.dense_fallback);
    CHECK(diagnostics.ncv == hicx::lanczos_subspace_size(3, block.rows()));
    for (std::size_t k = 0; k < which.size(); ++k) {
        INFO("eigenvector " << which[k]);
        CHECK(got.values[k] == doctest::Approx(reference.values[k]).epsilon(1e-10));
        // Both apply the same sign rule, so no alignment is needed.
        double worst = 0.0;
        for (std::size_t i = 0; i < got.vectors[k].size(); ++i) {
            worst = std::max(worst, std::fabs(got.vectors[k][i] - reference.vectors[k][i]));
        }
        CHECK(worst <= 1e-9);
    }
}

TEST_CASE("lanczos eigenvectors are byte-identical for any thread count") {
    // chrX of the 50 kb matrix is below the entry count at which the product
    // is split, so the operator is exercised on a larger block: chr3R of the
    // full-resolution small_test_matrix.h5, 5,582 bins.
    hicx::ToolMatrix hic = hicx::ToolMatrix::load(data_path("small_test_matrix.h5"));
    hicx::keep_only_chromosomes(hic.data(), {"chr3R"});
    hicx::CsrMatrix block = hic.matrix();
    block.materialize_full();
    hicx::obs_exp_lieberman_in_place(block, block.rows(), 1);
    const std::vector<int> which{1, 2};
    const hicx::EigenResult one = hicx::covariance_eigenvectors_lanczos(block, which, 1);
    for (const int threads : {2, 7, 16}) {
        const hicx::EigenResult many = hicx::covariance_eigenvectors_lanczos(block, which, threads);
        for (std::size_t k = 0; k < which.size(); ++k) {
            INFO("threads " << threads << " eigenvector " << which[k]);
            REQUIRE(many.vectors[k].size() == one.vectors[k].size());
            CHECK(std::memcmp(many.vectors[k].data(), one.vectors[k].data(),
                              one.vectors[k].size() * sizeof(double)) == 0);
            CHECK(std::memcmp(&many.values[k], &one.values[k], sizeof(double)) == 0);
        }
    }
}

TEST_CASE("a chromosome at or below the dense limit is solved densely") {
    hicx::ToolMatrix small = hicx::ToolMatrix::load(data_path("small_test_matrix_50kb_res.h5"));
    hicx::keep_only_chromosomes(small.data(), {"chr4"});  // 28 bins
    hicx::CsrMatrix tiny = small.matrix();
    tiny.materialize_full();
    hicx::obs_exp_lieberman_in_place(tiny, tiny.rows(), 1);
    hicx::LanczosDiagnostics diagnostics;
    const hicx::EigenResult got =
        hicx::covariance_eigenvectors_lanczos(tiny, {1, 2, 40}, 1, &diagnostics);
    CHECK(diagnostics.dense_fallback);
    CHECK(got.vectors[0].size() == 28);
    CHECK(got.vectors[2].empty());
}
