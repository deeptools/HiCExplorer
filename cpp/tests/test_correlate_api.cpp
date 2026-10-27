// The core code behind hicexplorer.hicCorrelate (hicx/correlate_api.hpp).
// Their equivalence with HiCExplorer 3.7.6 is tested from Python
// (hicexplorer/test/api); these tests pin the C++ contracts.

#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "hicx/correlate_api.hpp"
#include "hicx/python_error.hpp"

namespace corr = hicx::correlate;

namespace {

corr::SparseInput<double> csr(std::int64_t rows, std::int64_t cols, std::vector<std::int64_t> indptr,
                              std::vector<std::int64_t> minor, std::vector<double> data) {
    corr::SparseInput<double> m;
    m.format = corr::SparseFormat::Csr;
    m.rows = rows;
    m.cols = cols;
    m.indptr = std::move(indptr);
    m.minor = std::move(minor);
    m.data = std::move(data);
    return m;
}

}  // namespace

TEST_CASE("get_vectors: the values of the other matrix over the union of the entries") {
    // a = [[1 0 2], [0 3 0]], b = [[0 0 5], [4 3 0]]; the union is (0,0) (0,2)
    // (1,0) (1,1); (a + b + 1) - a is b + 1 there, so values1 is b and values2 is a
    const auto a = csr(2, 3, {0, 2, 3}, {0, 2, 1}, {1, 2, 3});
    const auto b = csr(2, 3, {0, 1, 3}, {2, 0, 1}, {5, 4, 3});
    const auto [first, second] = corr::get_vectors(a, b);
    CHECK(first == std::vector<double>{0, 5, 4, 3});
    CHECK(second == std::vector<double>{1, 2, 0, 3});
}

TEST_CASE("get_vectors: unsorted indices come out in scipy's general order") {
    // a stores columns 2, 0, 1 of its row (not sorted), b column 1; the general
    // scan lists the columns of every result in reverse order of appearance
    const auto a = csr(1, 3, {0, 3}, {2, 0, 1}, {3, 1, 2});
    const auto b = csr(1, 3, {0, 1}, {1}, {10});
    const auto [first, second] = corr::get_vectors(a, b);
    CHECK(first == std::vector<double>{0, 0, 10});
    CHECK(second == std::vector<double>{3, 1, 2});
}

TEST_CASE("get_vectors: integers wrap like numpy's") {
    corr::SparseInput<std::int32_t> a;
    a.format = corr::SparseFormat::Csr;
    a.rows = 1;
    a.cols = 1;
    a.indptr = {0, 1};
    a.minor = {0};
    a.data = {std::numeric_limits<std::int32_t>::max()};
    auto b = a;
    b.data = {1};
    const auto [first, second] = corr::get_vectors(a, b);
    CHECK(first == std::vector<std::int32_t>{1});
    CHECK(second == std::vector<std::int32_t>{std::numeric_limits<std::int32_t>::max()});
}

TEST_CASE("matrix_minimum, looks_like_distance_matrix") {
    CHECK(corr::matrix_minimum({3.0, -1.0, 2.0}) == -1.0);
    CHECK(std::isnan(corr::matrix_minimum({3.0, std::nan(""), 2.0})));
    CHECK_THROWS_AS((void)corr::matrix_minimum({}), hicx::PythonError);
    CHECK(corr::looks_like_distance_matrix({0, 1, 1, 0}, 2, 2));
    CHECK_FALSE(corr::looks_like_distance_matrix({1, 1, 1, 1}, 2, 2));
    CHECK_FALSE(corr::looks_like_distance_matrix({0, 1, 2, 0}, 2, 2));
    CHECK_FALSE(corr::looks_like_distance_matrix({0, -1, -1, 0}, 2, 2));
    CHECK_FALSE(corr::looks_like_distance_matrix({0, 1, 0, 1, 0, 0}, 2, 3));
}

TEST_CASE("complete_linkage_of_rows: scipy's errors and a three row example") {
    CHECK_THROWS_WITH((void)corr::complete_linkage_of_rows({1.0}, 1, 1),
                      "The number of observations cannot be determined on an empty distance matrix.");
    CHECK_THROWS_WITH((void)corr::complete_linkage_of_rows({1.0, std::nan(""), 2.0, 3.0}, 2, 2),
                      "The condensed distance matrix must contain only finite values.");
    CHECK_THROWS_WITH((void)corr::complete_linkage_of_rows({1e308, 0.0, -1e308, 0.0}, 2, 2),
                      "The condensed distance matrix must contain only finite values.");
    const std::vector<double> z = corr::complete_linkage_of_rows({0, 0, 1, 0, 10, 0}, 3, 2);
    REQUIRE(z.size() == 8);
    CHECK(z[0] == 0);
    CHECK(z[1] == 1);
    CHECK(z[2] == 1.0);
    CHECK(z[3] == 2);
    CHECK(z[6] == 10.0);
    CHECK(z[7] == 3);
}

