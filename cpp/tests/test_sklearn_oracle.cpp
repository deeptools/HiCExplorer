// Randomized comparison against the original libraries: k-means labels and
// cluster centres (bit for bit), agglomerative ward labels, and the draws of
// numpy's legacy RandomState. The expected values are written by
// cpp/scripts/sklearn_oracle_fixture.py from scikit-learn 1.3.2 and numpy 1.26.4
// (two OpenMP threads; see the generator for why the centres of datasets above
// 512 samples are compared to rounding error only, the labels always exactly).

#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "doctest/doctest.h"
#include "hicx/clustering.hpp"
#include "hicx/numpy_random.hpp"

namespace {

// numpy's einsum and the BLAS kernels differ between the machines, so each has
// its own fixtures (see the generator).
#if defined(__aarch64__) || defined(_M_ARM64)
const std::string kFixtures = std::string(HICX_SKLEARN_FIXTURE_DIR) + "/arm64";
#else
const std::string kFixtures = std::string(HICX_SKLEARN_FIXTURE_DIR) + "/x86_64";
#endif

// The BLAS the fixtures were written with, "OpenBLAS <version> <kernel>".
std::string fixture_blas() {
    std::ifstream in(kFixtures + "/blas.txt");
    std::string line;
    std::getline(in, line);
    return line;
}

class Reader {
  public:
    explicit Reader(const std::string& name) : in_(kFixtures + "/" + name, std::ios::binary) {
        REQUIRE_MESSAGE(in_.good(), "cannot open fixture " << name);
    }
    template <class T>
    T one() {
        T value;
        in_.read(reinterpret_cast<char*>(&value), sizeof(T));
        REQUIRE(in_.good());
        return value;
    }
    template <class T>
    std::vector<T> many(std::int64_t count) {
        std::vector<T> values(static_cast<std::size_t>(count));
        in_.read(reinterpret_cast<char*>(values.data()),
                 static_cast<std::streamsize>(values.size() * sizeof(T)));
        REQUIRE(in_.good());
        return values;
    }

  private:
    std::ifstream in_;
};

bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() &&
           (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

}  // namespace

TEST_CASE("k-means and ward agree with scikit-learn on random data") {
    // The distances come from BLAS, whose kernels round differently; the bits
    // are only comparable with the BLAS the fixtures were written with.
    const std::string expected_blas = fixture_blas();
    const std::string loaded_blas = hicx::cluster::detail::blas_kernel();
    const bool same_blas = loaded_blas == expected_blas;
    if (!same_blas) {
        MESSAGE("k-means part skipped (ward still checked): fixtures from " << expected_blas
                                                                           << ", running with "
                                                                           << loaded_blas);
    }
    Reader in("kmeans_oracle.bin");
    const int cases = in.one<std::int32_t>();
    int checked = 0;
    for (int c = 0; c < cases; ++c) {
        const std::int64_t n = in.one<std::int32_t>();
        const std::int64_t d = in.one<std::int32_t>();
        const std::int64_t k = in.one<std::int32_t>();
        const auto seed = static_cast<std::uint32_t>(in.one<std::int32_t>());
        hicx::cluster::Samples X;
        X.samples = n;
        X.features = d;
        X.values = in.many<double>(n * d);
        const auto labels = in.many<std::int64_t>(n);
        const auto centers = in.many<double>(k * d);
        const auto ward = in.many<std::int64_t>(n);
        CAPTURE(c);
        CAPTURE(n);
        CAPTURE(d);
        CAPTURE(k);
        if (same_blas) {
            const hicx::cluster::detail::KMeansFit fit = hicx::cluster::detail::kmeans_fit(X, k, seed);
            CHECK(fit.labels == labels);
            if (n <= 512) {
                // up to two blocks of the assignment step the sums are
                // independent of the thread count
                CHECK(same_bits(fit.centers, centers));
            } else {
                REQUIRE(fit.centers.size() == centers.size());
                for (std::size_t i = 0; i < centers.size(); ++i) {
                    CHECK(fit.centers[i] == doctest::Approx(centers[i]).epsilon(1e-12));
                }
            }
            CHECK(hicx::cluster::kmeans(X, k, seed) == labels);
        }
        if (n >= 2) {
            CHECK(hicx::cluster::ward(X, k) == ward);
        }
        ++checked;
    }
    CHECK(checked == cases);
}

TEST_CASE("RandomState draws agree with numpy for many seeds") {
    Reader in("random_state_oracle.bin");
    const std::vector<double> p = {0.05, 0.2, 0.1, 0.25, 0.15, 0.05, 0.2};
    const int seeds = in.one<std::int32_t>();
    for (int s = 0; s < seeds; ++s) {
        const auto seed = in.one<std::uint32_t>();
        const auto samples = in.many<double>(1500);
        const auto uniform = in.many<double>(7);
        const auto choices = in.many<std::int64_t>(40);
        CAPTURE(seed);
        hicx::npy::RandomState state(seed);
        std::vector<double> got(1500);
        for (double& value : got) value = state.random_sample();
        CHECK(same_bits(got, samples));
        CHECK(same_bits(state.uniform(7), uniform));
        hicx::npy::RandomState chooser(seed);
        std::vector<std::int64_t> picked;
        for (int i = 0; i < 40; ++i) picked.push_back(chooser.choice(p));
        CHECK(picked == choices);
    }
}

TEST_CASE("einsum row norm agrees with numpy for every length") {
    Reader in("einsum_oracle.bin");
    const int rows = in.one<std::int32_t>();
    for (int r = 0; r < rows; ++r) {
        const std::int64_t length = in.one<std::int32_t>();
        const auto x = in.many<double>(length);
        const double expected = in.one<double>();
        CAPTURE(length);
        const double got = hicx::cluster::detail::einsum_row_norm_squared(x.data(), length);
        CHECK(std::memcmp(&got, &expected, sizeof(double)) == 0);
    }
}
