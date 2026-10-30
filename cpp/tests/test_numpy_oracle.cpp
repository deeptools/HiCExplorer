// Comparison with numpy 1.26.4 and Python on randomised inputs. The expected
// values are in tests/fixtures/numpy_oracle/, written by
// scripts/gen_numpy_oracle_fixture.py; the inputs are regenerated here from
// the seeds on each fixture line by the same splitmix64 stream.

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "hicx/numpy_compat.hpp"
#include "hicx/numpy_sort.hpp"
#include "hicx/simd_reduce.hpp"

namespace {

std::vector<std::string> fixture_lines(const std::string& name) {
    std::ifstream in(std::string(HICX_NUMPY_ORACLE_DIR) + "/" + name);
    REQUIRE_MESSAGE(in.good(), "cannot open fixture " << name);
    std::vector<std::string> lines;
    for (std::string line; std::getline(in, line);) {
        if (!line.empty()) {
            lines.push_back(line);
        }
    }
    return lines;
}

// Outputs 1..n of the splitmix64 sequence started at `seed`.
std::vector<std::uint64_t> stream(std::uint64_t seed, std::size_t n) {
    constexpr std::uint64_t gamma = 0x9E3779B97F4A7C15ULL;
    std::vector<std::uint64_t> out(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::uint64_t z = seed * gamma + (i + 1) * gamma;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        out[i] = z ^ (z >> 31);
    }
    return out;
}

std::vector<double> unit(std::uint64_t seed, std::size_t n) {
    std::vector<double> v(n);
    const auto u = stream(seed, n);
    for (std::size_t i = 0; i < n; ++i) {
        v[i] = static_cast<double>(u[i] >> 11) * 0x1p-53;
    }
    return v;
}

std::vector<double> sort_values(int mode, std::uint64_t levels, std::size_t n, std::uint64_t seed) {
    const auto u = stream(seed, n);
    const auto u2 = stream(seed + 1000003, n);
    std::vector<double> v(n);
    switch (mode) {
        case 0:
            return unit(seed, n);
        case 1:
            for (std::size_t i = 0; i < n; ++i) {
                v[i] = static_cast<double>(u[i] % levels);
            }
            return v;
        case 2:
            for (std::size_t i = 0; i < n; ++i) {
                v[i] = (u2[i] % 7 == 0) ? std::nan("") : static_cast<double>(u[i] % levels);
            }
            return v;
        case 3:
            v = unit(seed, n);
            for (std::size_t i = 0; i < n; ++i) {
                if (u2[i] % 3 == 0) {
                    v[i] = 0.0;
                }
            }
            return v;
        case 4:
            for (std::size_t i = 0; i < n; ++i) {
                v[i] = (u[i] % 10 == 0) ? 0.0 : static_cast<double>(n - i);
            }
            return v;
        default:
            for (std::size_t i = 0; i < n; ++i) {
                v[i] = static_cast<double>(i < n / 2 ? i : n - i);
            }
            return v;
    }
}

std::uint64_t fnv(const std::vector<std::int64_t>& index) {
    std::uint64_t h = 0xCBF29CE484222325ULL;
    for (const std::int64_t x : index) {
        std::uint64_t bits = static_cast<std::uint64_t>(x);
        for (int b = 0; b < 8; ++b) {
            h = (h ^ ((bits >> (8 * b)) & 0xFF)) * 0x100000001B3ULL;
        }
    }
    return h;
}

}  // namespace

TEST_CASE("argsort and argpartition give numpy's permutation on random inputs") {
    int checked = 0;
    for (const std::string& line : fixture_lines("sort.txt")) {
        std::istringstream in(line);
        std::string kind;
        in >> kind;
        if (kind == "adversary") {
            std::uint64_t expected = 0;
            in >> expected;
            std::vector<std::int64_t> keys;
            for (std::int64_t k; in >> k;) {
                keys.push_back(k);
            }
            CHECK(fnv(hicx::npy::argsort(keys)) == expected);
        } else if (kind == "argsort") {
            int mode = 0;
            std::uint64_t levels = 0, seed = 0, expected = 0;
            std::size_t n = 0;
            in >> mode >> levels >> n >> seed >> expected;
            const auto v = sort_values(mode, levels, n, seed);
            CHECK_MESSAGE(fnv(hicx::npy::argsort(v)) == expected,
                          "mode " << mode << " n " << n << " seed " << seed);
        } else {
            int mode = 0;
            std::uint64_t levels = 0, seed = 0, expected = 0;
            std::size_t n = 0;
            std::int64_t k = 0;
            in >> mode >> levels >> n >> seed >> k >> expected;
            const auto v = sort_values(mode, levels, n, seed);
            CHECK_MESSAGE(fnv(hicx::npy::argpartition(v, k)) == expected,
                          "mode " << mode << " n " << n << " seed " << seed << " kth " << k);
        }
        ++checked;
    }
    CHECK(checked > 600);
}

namespace {

std::vector<double> sum_values(int mode, std::size_t n, std::uint64_t seed) {
    std::vector<double> v(n);
    const auto u = stream(seed, n);
    const auto e = stream(seed + 7, n);
    const auto x = unit(seed, n);
    for (std::size_t i = 0; i < n; ++i) {
        switch (mode) {
            case 0:
                v[i] = x[i];
                break;
            case 1:
                v[i] = (x[i] - 0.5) * std::exp2(static_cast<double>(static_cast<std::int64_t>(e[i] % 33) - 16));
                break;
            case 2:
                v[i] = static_cast<double>(static_cast<std::int64_t>(u[i] % 2001) - 1000) / 7;
                break;
            default:
                v[i] = 1e8 + x[i] * 1e-3;
                break;
        }
    }
    return v;
}

}  // namespace

TEST_CASE("pairwise_sum is bit identical to np.sum on random float64 and float32 arrays") {
    int checked = 0;
    for (const std::string& line : fixture_lines("sum.txt")) {
        std::istringstream in(line);
        std::string tag, kind;
        int mode = 0;
        std::size_t n = 0;
        std::uint64_t seed = 0, bits = 0;
        in >> tag >> kind >> mode >> n >> seed >> bits;
        const std::vector<double> v = sum_values(mode, n, seed);
        if (kind == "f64") {
            const double sum = hicx::npy::pairwise_sum(v.data(), n);
            std::uint64_t got = 0;
            std::memcpy(&got, &sum, 8);
            CHECK_MESSAGE(got == bits, "f64 mode " << mode << " n " << n);
            const double dispatched = hicx::simd::pairwise_sum(v.data(), n);
            std::memcpy(&got, &dispatched, 8);
            CHECK_MESSAGE(got == bits, "simd f64 mode " << mode << " n " << n);
        } else {
            std::vector<float> f(v.begin(), v.end());
            const float sum = hicx::npy::pairwise_sum(f.data(), n);
            std::uint32_t got = 0;
            std::memcpy(&got, &sum, 4);
            CHECK_MESSAGE(got == bits, "f32 mode " << mode << " n " << n);
        }
        ++checked;
    }
    CHECK(checked > 400);
}

TEST_CASE("float_repr and float32_repr match Python's repr and numpy's str") {
    int checked = 0;
    for (const std::string& line : fixture_lines("repr.txt")) {
        std::istringstream in(line);
        std::string tag, kind, text;
        std::uint64_t bits = 0;
        in >> tag >> kind >> bits >> text;
        if (kind == "f64") {
            double d = 0;
            std::memcpy(&d, &bits, 8);
            CHECK_MESSAGE(hicx::npy::float_repr(d) == text, "bits " << bits);
        } else {
            const auto b32 = static_cast<std::uint32_t>(bits);
            float f = 0;
            std::memcpy(&f, &b32, 4);
            CHECK_MESSAGE(hicx::npy::float32_repr(f) == text, "bits " << bits);
        }
        ++checked;
    }
    CHECK(checked > 9000);
}

TEST_CASE("array_str matches str() of a numpy string array") {
    for (const std::string& line : fixture_lines("arraystr.txt")) {
        const auto first = line.find('\t');
        const auto second = line.find('\t', first + 1);
        std::istringstream items_in(line.substr(first + 1, second - first - 1));
        std::vector<std::string> items;
        for (std::string w; items_in >> w;) {
            items.push_back(w);
        }
        std::string expected = line.substr(second + 1);
        std::string got = hicx::npy::array_str(items);
        std::string escaped;
        for (const char c : got) {
            if (c == '\n') {
                escaped += "\\n";
            } else {
                escaped += c;
            }
        }
        CHECK(escaped == expected);
    }
}

TEST_CASE("int_with_thousands_separator") {
    using hicx::npy::int_with_thousands_separator;
    CHECK(int_with_thousands_separator(0) == "0");
    CHECK(int_with_thousands_separator(999) == "999");
    CHECK(int_with_thousands_separator(1000) == "1,000");
    CHECK(int_with_thousands_separator(-1234567) == "-1,234,567");
    CHECK(int_with_thousands_separator(INT64_MIN) == "-9,223,372,036,854,775,808");
}

#include "hicx/numpy_stats.hpp"

namespace {

std::uint64_t canonical_bits(double x) {
    if (std::isnan(x)) {
        return 0x7FF8000000000000ULL;
    }
    std::uint64_t bits = 0;
    std::memcpy(&bits, &x, 8);
    return bits;
}

std::uint64_t fnv_bits(const std::vector<double>& values) {
    std::uint64_t h = 0xCBF29CE484222325ULL;
    for (const double v : values) {
        const std::uint64_t bits = canonical_bits(v);
        for (int b = 0; b < 8; ++b) {
            h = (h ^ ((bits >> (8 * b)) & 0xFF)) * 0x100000001B3ULL;
        }
    }
    return h;
}

double from_bits(std::uint64_t bits) {
    double x = 0;
    std::memcpy(&x, &bits, 8);
    return x;
}

std::vector<double> stat_values(int mode, std::size_t n, std::uint64_t seed) {
    const auto u = stream(seed, n);
    const auto u2 = stream(seed + 5, n);
    std::vector<double> v = unit(seed, n);
    for (std::size_t i = 0; i < n; ++i) {
        switch (mode) {
            case 0:
                v[i] = v[i] * 100.0 - 20.0;
                break;
            case 1:
                v[i] = static_cast<double>(u[i] % 6);
                break;
            case 2:
                v[i] = (u2[i] % 5 == 0) ? std::nan("") : v[i] - 0.5;
                break;
            default: {
                const std::uint64_t r = u2[i] % 9;
                v[i] = r == 0 ? std::nan("")
                       : r == 1 ? HUGE_VAL
                       : r == 2 ? -HUGE_VAL
                       : r == 3 ? -0.0
                                : v[i] - 0.5;
            }
        }
    }
    return v;
}

}  // namespace

TEST_CASE("median, quantile, percentile, linspace and searchsorted match numpy") {
    int checked = 0;
    for (const std::string& line : fixture_lines("stats.txt")) {
        std::istringstream in(line);
        std::string kind;
        in >> kind;
        if (kind == "linspace") {
            std::uint64_t a = 0, b = 0, expected = 0;
            std::int64_t num = 0;
            in >> a >> b >> num >> expected;
            CHECK_MESSAGE(fnv_bits(hicx::npy::linspace(from_bits(a), from_bits(b), num)) == expected,
                          line);
        } else if (kind == "median") {
            int mode = 0;
            std::size_t n = 0;
            std::uint64_t seed = 0, expected = 0;
            in >> mode >> n >> seed >> expected;
            CHECK_MESSAGE(canonical_bits(hicx::npy::median(stat_values(mode, n, seed))) == expected,
                          line);
        } else if (kind == "nanquantile") {
            int mode = 0;
            std::size_t n = 0, count = 0;
            std::uint64_t seed = 0, expected = 0;
            in >> mode >> n >> seed >> count >> expected;
            std::vector<double> qs = {0.0, 1.0, 0.5, 0.25};
            for (const double q : unit(seed + 9, count - 4)) {
                qs.push_back(q);
            }
            CHECK_MESSAGE(
                fnv_bits(hicx::npy::nanquantile_linear(stat_values(mode, n, seed), qs)) == expected,
                line);
        } else if (kind == "percentile") {
            int mode = 0;
            std::size_t n = 0;
            std::uint64_t seed = 0, expected = 0;
            double q = 0;
            in >> mode >> n >> seed >> q >> expected;
            CHECK_MESSAGE(canonical_bits(hicx::npy::percentile_linear(stat_values(mode, n, seed), q)) ==
                              expected,
                          line);
        } else {
            int mode = 0, sorted_flag = 0;
            std::size_t n = 0;
            std::uint64_t seed = 0, expected = 0;
            in >> mode >> n >> seed >> sorted_flag >> expected;
            std::vector<double> arr = stat_values(mode, n, seed);
            if (sorted_flag) {
                // np.sort: NaN last.
                std::sort(arr.begin(), arr.end(),
                          [](double x, double y) { return x < y || (y != y && x == x); });
            }
            std::vector<double> keys = stat_values(mode, 40, seed + 3);
            keys.push_back(std::nan(""));
            keys.push_back(HUGE_VAL);
            keys.push_back(-HUGE_VAL);
            CHECK_MESSAGE(fnv(hicx::npy::searchsorted_right(arr, keys)) == expected, line);
        }
        ++checked;
    }
    CHECK(checked > 900);
}
