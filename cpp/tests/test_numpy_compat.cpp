#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "hicx/numpy_compat.hpp"

using hicx::npy::array_str;
using hicx::npy::float_repr;
using hicx::npy::int_with_thousands_separator;
using hicx::npy::pairwise_sum;

TEST_CASE("pairwise_sum reproduces numpy, not a sequential loop") {
    // The reference values were produced with numpy 1.26.4:
    //   a = np.array([1.0 / (i + 1) for i in range(n)]); repr(a.sum())
    // A sequential loop over the same values gives a different last digit,
    // which is what makes this test meaningful.
    std::vector<double> a;
    for (int i = 0; i < 1000; ++i) {
        a.push_back(1.0 / (i + 1));
    }
    double naive = 0.0;
    for (const double value : a) {
        naive += value;
    }
    CHECK(float_repr(pairwise_sum(a)) == "7.485470860550345");
    CHECK(float_repr(naive) == "7.485470860550343");

    std::vector<double> b;
    for (int i = 0; i < 137; ++i) {
        b.push_back(1.0 / (i + 1));
    }
    CHECK(float_repr(pairwise_sum(b)) == "5.50084178584451");

    // Blocks of eight accumulators keep the small values alive next to a large
    // one, where a sequential sum loses them completely.
    std::vector<double> c(100, 1.0);
    c[0] = 1e16;
    CHECK(float_repr(pairwise_sum(c)) == "1.0000000000000084e+16");

    // Short arrays fall back to the sequential loop.
    const std::vector<double> tiny{1.0, 2.0, 3.0};
    CHECK(pairwise_sum(tiny) == doctest::Approx(6.0));
    CHECK(pairwise_sum(std::vector<double>{}) == 0.0);
}

TEST_CASE("float_repr matches Python's repr") {
    CHECK(float_repr(17548966.536917936) == "17548966.536917936");
    CHECK(float_repr(1.3271409292967762e-05) == "1.3271409292967762e-05");
    CHECK(float_repr(0.12034209360197595) == "0.12034209360197595");
    CHECK(float_repr(1914.0146478655674) == "1914.0146478655674");
    CHECK(float_repr(0.1702102389180371) == "0.1702102389180371");
    CHECK(float_repr(742.330787829914) == "742.330787829914");
    CHECK(float_repr(0.0) == "0.0");
    CHECK(float_repr(-0.0) == "-0.0");
    CHECK(float_repr(5.0) == "5.0");
    CHECK(float_repr(-2.5) == "-2.5");
    CHECK(float_repr(0.0001) == "0.0001");
    CHECK(float_repr(0.00001) == "1e-05");
    CHECK(float_repr(1e16) == "1e+16");
    CHECK(float_repr(1e17) == "1e+17");
    CHECK(float_repr(1.2345678901234567e16) == "1.2345678901234568e+16");
    CHECK(float_repr(1e22) == "1e+22");
    CHECK(float_repr(1234567890123456.0) == "1234567890123456.0");
}

TEST_CASE("int_with_thousands_separator matches Python's {:,} format") {
    CHECK(int_with_thousands_separator(0) == "0");
    CHECK(int_with_thousands_separator(999) == "999");
    CHECK(int_with_thousands_separator(1000) == "1,000");
    CHECK(int_with_thousands_separator(11104) == "11,104");
    CHECK(int_with_thousands_separator(1661678) == "1,661,678");
    CHECK(int_with_thousands_separator(61804782) == "61,804,782");
    CHECK(int_with_thousands_separator(-1234567) == "-1,234,567");
}

TEST_CASE("array_str matches str() of a numpy string array") {
    CHECK(array_str({"chrom", "start", "end"}) == "['chrom' 'start' 'end']");
    CHECK(array_str({"chrom", "start", "end", "weight"}) ==
          "['chrom' 'start' 'end' 'weight']");
    CHECK(array_str({"chrom", "start", "end", "KR", "VC", "VC_SQRT"}) ==
          "['chrom' 'start' 'end' 'KR' 'VC' 'VC_SQRT']");
    CHECK(array_str({}) == "[]");
}

TEST_CASE("PairwiseSumStream is bit identical to pairwise_sum over the same sequence") {
    // Sizes around the 8 accumulator unroll, the 128 element pairwise block
    // and the 8192 element reduction buffer, where a streaming implementation
    // that got the buffering wrong would first diverge.
    std::uint64_t state = 0x9e3779b97f4a7c15ULL;
    const auto next = [&state]() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        return static_cast<double>(state % 2000003) / 7.0 - 142857.0;
    };
    for (const std::size_t n : {0UL, 1UL, 7UL, 8UL, 9UL, 127UL, 128UL, 129UL, 1000UL,
                                8191UL, 8192UL, 8193UL, 16384UL, 20000UL, 65537UL}) {
        std::vector<double> values(n);
        for (double& value : values) {
            value = next() * 1e-3;
        }
        hicx::npy::PairwiseSumStream<double> stream;
        for (const double value : values) {
            stream.add(value);
        }
        CAPTURE(n);
        CHECK(stream.count() == n);
        const double expected = pairwise_sum(values);
        const double got = stream.result();
        CHECK(std::memcmp(&expected, &got, sizeof(double)) == 0);

        std::vector<float> narrow(values.begin(), values.end());
        hicx::npy::PairwiseSumStream<float> stream32;
        for (const float value : narrow) {
            stream32.add(value);
        }
        const float expected32 = pairwise_sum(narrow);
        const float got32 = stream32.result();
        CHECK(std::memcmp(&expected32, &got32, sizeof(float)) == 0);
    }

    // An empty tail is not added: -0.0 stays -0.0, as in buffered_pairwise_sum.
    hicx::npy::PairwiseSumStream<double> negative_zero;
    for (int i = 0; i < 8192; ++i) {
        negative_zero.add(-0.0);
    }
    const std::vector<double> zeros(8192, -0.0);
    CHECK(std::signbit(negative_zero.result()) == std::signbit(pairwise_sum(zeros)));
}
