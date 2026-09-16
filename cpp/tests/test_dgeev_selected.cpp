// dgeev_selected_eigenvectors must return, bit for bit, the columns dgeev_
// itself returns (hicx/dgeev_selected.hpp). Checked on a real chromosome block
// for the real branch and on a small non-symmetric matrix, a unit-test fixture
// only, for the complex conjugate pair branch that a symmetric covariance does
// not reach.

#include <doctest/doctest.h>

#include <cmath>
#include <cstddef>
#include <cstring>
#include <string>
#include <vector>

#include "hicx/adjust_ops.hpp"
#include "hicx/dgeev_selected.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/transform_ops.hpp"

extern "C" {
void dgeev_(const char*, const char*, const int*, double*, const int*, double*, double*, double*,
            const int*, double*, const int*, double*, const int*, int*, std::size_t, std::size_t);
}

namespace {

std::string data_path(const std::string& relative) {
    return std::string(HICX_TEST_DATA_DIR) + "/" + relative;
}

struct Full {
    std::vector<double> wr, wi, vr;
};

Full full_dgeev(int n, std::vector<double> a) {
    Full out;
    out.wr.assign(static_cast<std::size_t>(n), 0.0);
    out.wi.assign(static_cast<std::size_t>(n), 0.0);
    out.vr.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
    const char cN = 'N';
    const char cV = 'V';
    int info = 0;
    int lwork = -1;
    double query = 0.0;
    double dummy = 0.0;
    dgeev_(&cN, &cV, &n, a.data(), &n, out.wr.data(), out.wi.data(), &dummy, &n, out.vr.data(), &n,
           &query, &lwork, &info, 1, 1);
    lwork = static_cast<int>(query);
    std::vector<double> work(static_cast<std::size_t>(lwork));
    dgeev_(&cN, &cV, &n, a.data(), &n, out.wr.data(), out.wi.data(), &dummy, &n, out.vr.data(), &n,
           work.data(), &lwork, &info, 1, 1);
    REQUIRE(info == 0);
    return out;
}

bool same_bits(const double* a, const double* b, std::size_t count) {
    return std::memcmp(a, b, count * sizeof(double)) == 0;
}

}  // namespace

TEST_CASE("the selected dgeev columns are dgeev_'s own on a real chromosome block") {
    hicx::pin_blas_to_one_thread();
    hicx::ToolMatrix hic = hicx::ToolMatrix::load(data_path("small_test_matrix_50kb_res.h5"));
    hicx::keep_only_chromosomes(hic.data(), {"chrX"});
    hicx::CsrMatrix block = hic.matrix();
    block.materialize_full();
    hicx::obs_exp_lieberman_in_place(block, block.rows(), 1);
    hicx::DenseSymmetric covariance = hicx::covariance_of_symmetric(block, 1);
    const int n = static_cast<int>(covariance.size());
    std::vector<double> copy(covariance.data(), covariance.data() + covariance.element_count());
    const Full reference = full_dgeev(n, copy);

    const std::vector<int> which{1, 2, 3, 5, n, n + 1};
    {
        hicx::DenseSymmetric input(n);
        std::memcpy(input.data(), copy.data(), copy.size() * sizeof(double));
        const hicx::EigenResult got = hicx::dgeev_selected_eigenvectors(input, which);
        REQUIRE(got.vectors.size() == which.size());
        for (std::size_t entry = 0; entry + 1 < which.size(); ++entry) {
            const int column = which[entry] - 1;
            INFO("column " << column);
            REQUIRE(got.vectors[entry].size() == static_cast<std::size_t>(n));
            CHECK(got.values[entry] == reference.wr[static_cast<std::size_t>(column)]);
            CHECK(same_bits(got.vectors[entry].data(),
                            reference.vr.data() + static_cast<std::size_t>(column) * n,
                            static_cast<std::size_t>(n)));
        }
        CHECK(got.vectors.back().empty());
        CHECK(input.size() == 0);  // consumed
    }
}

TEST_CASE("the selected dgeev columns follow dgeev_ through complex conjugate pairs") {
    hicx::pin_blas_to_one_thread();
    // A deterministic non-symmetric 300 by 300 matrix: a rotation-heavy
    // structure gives many complex pairs, large enough for the blocked
    // back-transform (NB = 128) to take several blocks.
    constexpr int n = 300;
    std::vector<double> a(static_cast<std::size_t>(n) * n);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            a[static_cast<std::size_t>(i) + static_cast<std::size_t>(j) * n] =
                std::sin(0.37 * i + 0.91 * j) + (i == j + 1 ? 2.0 : 0.0) - (j == i + 1 ? 2.0 : 0.0);
        }
    }
    const Full reference = full_dgeev(n, a);
    int complex_columns = 0;
    for (double value : reference.wi) {
        complex_columns += value != 0.0 ? 1 : 0;
    }
    REQUIRE(complex_columns > 10);

    std::vector<int> which;
    for (int k = 1; k <= n; k += 7) {
        which.push_back(k);
    }
    {
        std::vector<double> input = a;
        std::vector<double> imaginary;
        const hicx::EigenResult got =
            hicx::dgeev_selected_general(n, input.data(), which, &imaginary);
        for (std::size_t entry = 0; entry < which.size(); ++entry) {
            const int index = which[entry] - 1;
            int source = index;
            if (reference.wi[static_cast<std::size_t>(index)] < 0.0 && index > 0) {
                source = index - 1;
            }
            INFO("index " << index);
            CHECK(got.values[entry] == reference.wr[static_cast<std::size_t>(index)]);
            CHECK(imaginary[entry] == reference.wi[static_cast<std::size_t>(index)]);
            CHECK(same_bits(got.vectors[entry].data(),
                            reference.vr.data() + static_cast<std::size_t>(source) * n,
                            static_cast<std::size_t>(n)));
        }
    }
}
