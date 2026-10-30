// Randomised comparisons of the numerical building blocks of hicDetectStripes
// (tools/stripes_numpy.hpp) with the libraries they stand in for. The expected
// values are what numpy, Python's random module, scikit-image, OpenCV and SciPy
// returned for random inputs with ties, extremes and degenerate shapes; they
// are written by cpp/scripts/stripenn_oracle/make_primitive_fixtures.py into
// tests/fixtures/stripenn (files prim_*.npy). Every comparison is exact (bit
// equality for floating point values). The larger randomised comparison against
// a freshly built probe is cpp/scripts/stripenn_oracle/check_primitives.py.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include <doctest/doctest.h>

#include "../tools/stripes_npy.hpp"
#include "../tools/stripes_numpy.hpp"

namespace np = hicx::stripes::np;
namespace npy = hicx::stripes::npy;

namespace {

const std::string kFixtures = HICX_STRIPES_FIXTURE_DIR;

npy::Array load(const std::string& name) { return npy::load(kFixtures + "/prim_" + name + ".npy"); }

template <class T>
std::vector<T> load_as(const std::string& name) {
    return npy::as<T>(load(name));
}

bool same_bits(double a, double b) {
    return std::memcmp(&a, &b, sizeof(double)) == 0 || (std::isnan(a) && std::isnan(b));
}

}  // namespace

TEST_CASE("the Mersenne Twister seeded as Python seeds it gives Python's random() and choices()") {
    const auto seeds = load_as<std::int64_t>("seeds");
    const auto draws = load_as<double>("draws");
    const auto picks = load_as<std::int64_t>("picks");
    const auto cum = load_as<double>("cumweights");
    const auto wpicks = load_as<std::int64_t>("wpicks");
    const std::size_t per_seed = draws.size() / seeds.size();
    REQUIRE(per_seed == 1300);  // more than 2 * 624 words: the state is regenerated twice
    for (std::size_t i = 0; i < seeds.size(); ++i) {
        np::PyRandom rng(static_cast<std::uint64_t>(seeds[i]));
        for (std::size_t k = 0; k < per_seed; ++k) {
            CHECK(rng.random() == draws[i * per_seed + k]);
        }
        std::size_t at = i * 3 * 200;
        for (const std::size_t size : {std::size_t{1}, std::size_t{17}, std::size_t{99991}}) {
            np::PyRandom again(static_cast<std::uint64_t>(seeds[i]));
            for (int k = 0; k < 200; ++k) {
                CHECK(static_cast<std::int64_t>(again.choice_index(size)) == picks[at++]);
            }
        }
        np::PyRandom weighted(static_cast<std::uint64_t>(seeds[i]));
        for (int k = 0; k < 200; ++k) {
            CHECK(static_cast<std::int64_t>(weighted.weighted_choice_index(cum)) ==
                  wpicks[i * 200 + static_cast<std::size_t>(k)]);
        }
    }
}

TEST_CASE("median, nanmedian, quantile and round follow numpy on random inputs") {
    const auto values = load_as<double>("ragged_values");
    const auto lengths = load_as<std::int64_t>("ragged_lengths");
    const auto med = load_as<double>("ragged_median");
    const auto nanmed = load_as<double>("ragged_nanmedian");
    std::size_t at = 0;
    for (std::size_t i = 0; i < lengths.size(); ++i) {
        const std::vector<double> x(values.begin() + static_cast<std::ptrdiff_t>(at),
                                    values.begin() + static_cast<std::ptrdiff_t>(at) + lengths[i]);
        at += static_cast<std::size_t>(lengths[i]);
        CHECK(same_bits(np::median(x), med[i]));
        CHECK(same_bits(np::nanmedian(x), nanmed[i]));
    }
    for (const std::string tag : {"quantile", "quantile_one"}) {
        auto sorted = load_as<double>(tag == "quantile" ? "quantile_in" : "quantile_one_in");
        std::sort(sorted.begin(), sorted.end());
        const auto qs = load_as<double>("quantile_q");
        const auto expected = load_as<double>(tag == "quantile" ? "quantile_out" : "quantile_one_out");
        for (std::size_t i = 0; i < qs.size(); ++i) {
            const np::QuantileRanks r = np::quantile_ranks(static_cast<std::int64_t>(sorted.size()), qs[i]);
            CHECK(np::quantile_lerp(sorted[static_cast<std::size_t>(r.previous)],
                                    sorted[static_cast<std::size_t>(r.next)], r.gamma) == expected[i]);
        }
    }
    const auto in = load_as<double>("round_in");
    const auto out = load_as<double>("round_out");
    for (std::size_t i = 0; i < in.size(); ++i) {
        CHECK(same_bits(np::round_half_even(in[i]), out[i]));
    }
}

TEST_CASE("Canny stages and edges match scikit-image for random images, sigmas and thresholds") {
    for (const std::string tag : {"big", "tiny"}) {
        const npy::Array a = load("canny_" + tag + "_in");
        const auto images = npy::as<float>(a);
        const auto params = load_as<double>("canny_" + tag + "_params");
        const auto edges = load_as<std::uint8_t>("canny_" + tag + "_edges");
        const auto suppressed = load_as<float>("canny_" + tag + "_suppressed");
        const auto smoothed = load_as<float>("canny_" + tag + "_smoothed");
        const int count = static_cast<int>(a.shape[0]);
        const int rows = static_cast<int>(a.shape[1]);
        const int cols = static_cast<int>(a.shape[2]);
        const std::size_t plane = static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);
        std::size_t edge_pixels = 0;
        for (int i = 0; i < count; ++i) {
            np::ImageF image(rows, cols);
            image.v.assign(images.begin() + static_cast<std::ptrdiff_t>(i * plane),
                           images.begin() + static_cast<std::ptrdiff_t>((i + 1) * plane));
            const np::CannyStages st = np::canny_stages(image, params[3 * i], params[3 * i + 1], params[3 * i + 2]);
            CHECK(st.smoothed.v == std::vector<float>(smoothed.begin() + static_cast<std::ptrdiff_t>(i * plane),
                                                      smoothed.begin() + static_cast<std::ptrdiff_t>((i + 1) * plane)));
            CHECK(st.suppressed.v == std::vector<float>(suppressed.begin() + static_cast<std::ptrdiff_t>(i * plane),
                                                        suppressed.begin() + static_cast<std::ptrdiff_t>((i + 1) * plane)));
            const std::vector<std::uint8_t> expected(edges.begin() + static_cast<std::ptrdiff_t>(i * plane),
                                                     edges.begin() + static_cast<std::ptrdiff_t>((i + 1) * plane));
            CHECK(st.edges == expected);
            CHECK(np::canny(image, params[3 * i], params[3 * i + 1], params[3 * i + 2]) == expected);
            for (const auto e : expected) {
                edge_pixels += e;
            }
        }
        if (tag == "big") {
            CHECK(edge_pixels > 100);  // the comparison is not one of empty edge maps
        }
    }
}

TEST_CASE("the box filter matches OpenCV for kernels of 1, 3, 5 and 7 and for planes smaller than the kernel") {
    const npy::Array a = load("box_in");
    const auto planes = npy::as<double>(a);
    const auto ks = load_as<std::int64_t>("box_k");
    const auto expected = load_as<double>("box_out");
    const int rows = static_cast<int>(a.shape[1]);
    const int cols = static_cast<int>(a.shape[2]);
    const std::size_t plane = static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);
    for (std::size_t i = 0; i < ks.size(); ++i) {
        const std::vector<double> in(planes.begin() + static_cast<std::ptrdiff_t>(i * plane),
                                     planes.begin() + static_cast<std::ptrdiff_t>((i + 1) * plane));
        const std::vector<double> want(expected.begin() + static_cast<std::ptrdiff_t>(i * plane),
                                       expected.begin() + static_cast<std::ptrdiff_t>((i + 1) * plane));
        CHECK(np::box_filter(in, rows, cols, static_cast<int>(ks[i])) == want);
    }
    for (const std::string name : {"small", "row"}) {
        const npy::Array p = load("box_" + name + "_in");
        const auto in = npy::as<double>(p);
        const auto want = load_as<double>("box_" + name + "_out");
        const std::size_t size = in.size();
        std::size_t at = 0;
        for (const int k : {1, 3, 5}) {
            const std::vector<double> got = np::box_filter(in, static_cast<int>(p.shape[0]), static_cast<int>(p.shape[1]), k);
            CHECK(got == std::vector<double>(want.begin() + static_cast<std::ptrdiff_t>(at),
                                             want.begin() + static_cast<std::ptrdiff_t>(at + size)));
            at += size;
        }
    }
}

TEST_CASE("RGB to grey matches OpenCV for float32 and for 8 bit pixels") {
    const auto px = load_as<float>("gray_px");
    const auto gray = load_as<float>("gray_out");
    for (std::size_t i = 0; i < gray.size(); ++i) {
        CHECK(np::rgb_to_gray(px[3 * i], px[3 * i + 1], px[3 * i + 2]) == gray[i]);
    }
    const auto px8 = load_as<std::uint8_t>("gray8_px");
    const auto gray8 = load_as<std::uint8_t>("gray8_out");
    for (std::size_t i = 0; i < gray8.size(); ++i) {
        CHECK(np::rgb_to_gray_u8(px8[3 * i], px8[3 * i + 1], px8[3 * i + 2]) == gray8[i]);
    }
}

TEST_CASE("convolve2d_same_3x3 matches scipy.signal.convolve2d for random matrices and kernels") {
    const npy::Array a = load("conv_in");
    const auto ms = npy::as<std::uint8_t>(a);
    const auto kernels = load_as<std::int32_t>("conv_kernel");
    const auto expected = load_as<std::int32_t>("conv_out");
    const int rows = static_cast<int>(a.shape[1]);
    const int cols = static_cast<int>(a.shape[2]);
    const std::size_t plane = static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);
    for (std::size_t i = 0; i < a.shape[0]; ++i) {
        const std::vector<std::uint8_t> m(ms.begin() + static_cast<std::ptrdiff_t>(i * plane),
                                          ms.begin() + static_cast<std::ptrdiff_t>((i + 1) * plane));
        const std::int32_t* k = kernels.data() + 9 * i;
        const int kernel[3][3] = {{k[0], k[1], k[2]}, {k[3], k[4], k[5]}, {k[6], k[7], k[8]}};
        const std::vector<int> got = np::convolve2d_same_3x3(m, rows, cols, kernel);
        CHECK(got == std::vector<int>(expected.begin() + static_cast<std::ptrdiff_t>(i * plane),
                                      expected.begin() + static_cast<std::ptrdiff_t>((i + 1) * plane)));
    }
}
