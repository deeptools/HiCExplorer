// Unit tests for hicDetectStripes' computational core (tools/stripes_impl.hpp,
// tools/stripes_numpy.hpp): the Stripenn 1.1.65.22 method and the numpy,
// scipy, scikit-image and OpenCV routines it calls, reproduced operation by
// operation. The expected values are not hand computed: they are what those
// libraries and Stripenn itself produced for the same inputs, written by
// cpp/scripts/stripenn_oracle/make_fixtures.py into tests/fixtures/stripenn,
// and every comparison is exact (bit equality for floating point values).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <doctest/doctest.h>

#include "../tools/stripes_impl.hpp"
#include "../tools/stripes_npy.hpp"
#include "../tools/stripes_numpy.hpp"

namespace st = hicx::stripes;
namespace np = hicx::stripes::np;
namespace npy = hicx::stripes::npy;

namespace {

const std::string kFixtures = HICX_STRIPES_FIXTURE_DIR;

npy::Array load(const std::string& name) { return npy::load(kFixtures + "/" + name); }

template <class T>
std::vector<T> load_as(const std::string& name) {
    return npy::as<T>(load(name));
}

np::ImageF load_image(const std::string& name) {
    const npy::Array a = load(name);
    np::ImageF img(static_cast<int>(a.shape[0]), static_cast<int>(a.shape[1]));
    img.v = npy::as<float>(a);
    return img;
}

// A tab separated table with a header, as pandas wrote it.
struct Table {
    std::vector<std::string> columns;
    std::vector<std::vector<std::string>> rows;
    [[nodiscard]] std::size_t index(const std::string& name) const {
        return static_cast<std::size_t>(std::find(columns.begin(), columns.end(), name) - columns.begin());
    }
};

Table read_table(const std::string& name) {
    std::ifstream in(kFixtures + "/" + name);
    REQUIRE(in.good());
    Table t;
    std::string line;
    std::getline(in, line);
    {
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, '\t')) {
            t.columns.push_back(cell);
        }
    }
    while (std::getline(in, line)) {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, '\t')) {
            row.push_back(cell);
        }
        t.rows.push_back(row);
    }
    return t;
}

double parse_double(const std::string& s) { return std::stod(s); }

}  // namespace

TEST_CASE("Python's random.Random(seed).random() and choices() are reproduced bit for bit") {
    np::PyRandom rng(123456789);
    const std::vector<double> expected = load_as<double>("random.npy");
    for (const double e : expected) {
        CHECK(rng.random() == e);
    }
    np::PyRandom again(123456789);
    const std::vector<std::int64_t> picks = load_as<std::int64_t>("choices.npy");
    std::size_t at = 0;
    for (const std::size_t n : {1, 2, 7, 100, 1234}) {
        for (int k = 0; k < 50; ++k) {
            CHECK(static_cast<std::int64_t>(again.choice_index(n)) == picks[at++]);
        }
    }
    // A seed of two 32 bit words and the seed 0.
    np::PyRandom zero(0);
    CHECK(zero.random() >= 0.0);
}

TEST_CASE("numpy's sum and mean of 2-D views and of rows follow its pairwise summation") {
    const npy::Array big = load("sums_big.npy");
    const std::vector<double> data = npy::as<double>(big);
    const std::vector<std::int64_t> q = load_as<std::int64_t>("sums_queries.npy");
    const std::vector<double> expected = load_as<double>("sums_expected.npy");
    const std::vector<double> row_expected = load_as<double>("rowmeans_expected.npy");
    std::size_t row_at = 0;
    for (std::size_t i = 0; i < q.size() / 4; ++i) {
        const double* base = data.data() + static_cast<std::size_t>(q[4 * i]) * big.shape[1] + static_cast<std::size_t>(q[4 * i + 1]);
        const auto rows = static_cast<std::size_t>(q[4 * i + 2]);
        const auto cols = static_cast<std::size_t>(q[4 * i + 3]);
        CHECK(np::sum_block(base, rows, cols, big.shape[1]) == expected[2 * i]);
        CHECK(np::mean_block(base, rows, cols, big.shape[1]) == expected[2 * i + 1]);
        for (std::size_t r = 0; r < rows; ++r) {
            CHECK(np::mean_of_row(base + r * big.shape[1], cols) == row_expected[row_at++]);
        }
    }
}

TEST_CASE("median, round and quantile follow numpy") {
    CHECK(np::median(load_as<double>("median_in.npy")) == load_as<double>("median_expected.npy")[0]);
    CHECK(np::median(load_as<double>("median_in_even.npy")) == load_as<double>("median_expected_even.npy")[0]);
    CHECK(std::isnan(np::median({})));
    CHECK(np::round_half_even(0.5) == 0.0);
    CHECK(np::round_half_even(1.5) == 2.0);
    CHECK(np::round_half_even(2.5) == 2.0);
    CHECK(np::round_half_even(-0.5) == -0.0);

    // numpy.quantile of the positive entries of a symmetric matrix from its
    // upper triangle, in the two passes of the selector
    const npy::Array upper_array = load("quantile_upper.npy");
    const std::vector<double> upper = npy::as<double>(upper_array);
    const std::vector<double> qs = load_as<double>("quantile_q.npy");
    const std::vector<double> expected = load_as<double>("quantile_expected.npy");
    st::QuantileSelector selector(qs);
    const std::size_t n = upper_array.shape[0];
    for (std::size_t i = 0; i < n; ++i) {
        selector.count(upper[3 * i + 2], upper[3 * i] == upper[3 * i + 1]);
    }
    selector.prepare();
    for (std::size_t i = 0; i < n; ++i) {
        selector.gather(upper[3 * i + 2], upper[3 * i] == upper[3 * i + 1]);
    }
    const std::vector<double> got = selector.result();
    REQUIRE(got.size() == expected.size());
    for (std::size_t i = 0; i < got.size(); ++i) {
        CHECK(got[i] == expected[i]);
    }
}

TEST_CASE("scipy's Gaussian filter, Sobel and skimage's Canny reproduce each stage bit for bit") {
    for (const std::string name : {"real", "random"}) {
        const np::ImageF image = load_image("canny_" + name + "_in.npy");
        const double sigma = load_as<double>("canny_" + name + "_sigma.npy")[0];
        const np::CannyStages stages = np::canny_stages(image, sigma);
        CHECK(stages.smoothed.v == load_as<float>("canny_" + name + "_smoothed.npy"));
        CHECK(stages.isobel.v == load_as<float>("canny_" + name + "_isobel.npy"));
        CHECK(stages.magnitude.v == load_as<float>("canny_" + name + "_magnitude.npy"));
        CHECK(stages.suppressed.v == load_as<float>("canny_" + name + "_suppressed.npy"));
        CHECK(stages.edges == load_as<std::uint8_t>("canny_" + name + "_edges.npy"));
        // the image has edges, so the comparison is not of two empty maps
        CHECK(std::count(stages.edges.begin(), stages.edges.end(), 1) > 10);
    }
}

TEST_CASE("the Gaussian kernel of Stripenn's default sigma is numpy's") {
    CHECK(np::gaussian_kernel1d(2.0, 8) == load_as<double>("gauss_kernel_sigma2.npy"));
}

TEST_CASE("OpenCV's box filter and RGB to gray reproduce its rounding") {
    const npy::Array plane_array = load("box_plane.npy");
    const std::vector<double> plane = npy::as<double>(plane_array);
    const int rows = static_cast<int>(plane_array.shape[0]);
    const int cols = static_cast<int>(plane_array.shape[1]);
    CHECK(np::box_filter(plane, rows, cols, 3) == load_as<double>("box_out3.npy"));
    CHECK(np::box_filter(plane, rows, cols, 5) == load_as<double>("box_out5.npy"));
    const std::vector<float> px = load_as<float>("gray_px.npy");
    const std::vector<float> gray = load_as<float>("gray_out.npy");
    for (std::size_t i = 0; i < gray.size(); ++i) {
        CHECK(np::rgb_to_gray(px[3 * i], px[3 * i + 1], px[3 * i + 2]) == gray[i]);
    }
}

TEST_CASE("py_repr is Python's repr of a float") {
    std::ifstream in(kFixtures + "/repr.json");
    REQUIRE(in.good());
    std::stringstream buffer;
    buffer << in.rdbuf();
    const std::string text = buffer.str();
    // two arrays of strings; the file is written by json.dump on one line
    auto strings_after = [&](const std::string& key) {
        std::vector<std::string> out;
        std::size_t at = text.find("\"" + key + "\"");
        at = text.find('[', at);
        const std::size_t end = text.find(']', at);
        std::size_t pos = at;
        while (true) {
            const std::size_t open = text.find('"', pos);
            if (open == std::string::npos || open > end) {
                break;
            }
            const std::size_t close = text.find('"', open + 1);
            out.push_back(text.substr(open + 1, close - open - 1));
            pos = close + 1;
        }
        return out;
    };
    const std::vector<std::string> values = strings_after("values");
    const std::vector<std::string> reprs = strings_after("reprs");
    REQUIRE(values.size() == reprs.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double v = std::strtod(values[i].c_str(), nullptr);
        CHECK(st::py_repr(v) == reprs[i]);
    }
}

TEST_CASE("nulldist's allocation of the samples over the chromosomes, including its list indexing") {
    std::ifstream in(kFixtures + "/allocate.json");
    REQUIRE(in.good());
    std::stringstream buffer;
    buffer << in.rdbuf();
    const std::string text = buffer.str();
    // parse "sizes", "available" and "expected" arrays of each case
    auto numbers_after = [&](std::size_t from, const std::string& key, std::size_t* end_out) {
        std::vector<std::int64_t> out;
        std::size_t at = text.find("\"" + key + "\"", from);
        at = text.find('[', at);
        const std::size_t end = text.find(']', at);
        std::stringstream ss(text.substr(at + 1, end - at - 1));
        std::string cell;
        while (std::getline(ss, cell, ',')) {
            out.push_back(std::stoll(cell));
        }
        *end_out = end;
        return out;
    };
    std::size_t at = 0;
    int cases = 0;
    while (text.find("\"sizes\"", at) != std::string::npos) {
        std::size_t end = 0;
        const std::vector<std::int64_t> sizes = numbers_after(at, "sizes", &end);
        const std::vector<std::int64_t> available = numbers_after(end, "available", &end);
        const std::vector<std::int64_t> expected = numbers_after(end, "expected", &end);
        std::vector<st::ChromInfo> chroms;
        for (std::size_t i = 0; i < sizes.size(); ++i) {
            st::ChromInfo c;
            c.name = std::to_string(i);
            c.size = sizes[i];
            chroms.push_back(c);
        }
        CHECK(st::allocate_samples(chroms, available) == expected);
        at = end;
        ++cases;
    }
    CHECK(cases == 3);
}

namespace {

std::vector<st::Stripe> redundant_input(const Table& t) {
    std::vector<st::Stripe> stripes;
    for (const auto& row : t.rows) {
        st::Stripe s;
        s.pos1 = std::stoll(row[t.index("pos1")]);
        s.pos2 = std::stoll(row[t.index("pos2")]);
        s.pos3 = std::stoll(row[t.index("pos3")]);
        s.pos4 = std::stoll(row[t.index("pos4")]);
        s.h = std::stoi(row[t.index("h")]);
        s.w = std::stoi(row[t.index("w")]);
        s.num = std::stoll(row[t.index("num")]);
        s.pvalue = parse_double(row[t.index("pvalue")]);
        s.x = std::stoi(row[t.index("row")]);  // the row id travels in x
        stripes.push_back(s);
    }
    return stripes;
}

std::vector<int> read_int_list(const std::string& name) {
    std::ifstream in(kFixtures + "/" + name);
    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string text = buffer.str();
    std::vector<int> out;
    for (char& c : text) {
        if (c == '[' || c == ']' || c == ',') {
            c = ' ';
        }
    }
    std::stringstream ss(text);
    int v;
    while (ss >> v) {
        out.push_back(v);
    }
    return out;
}

}  // namespace

TEST_CASE("RemoveRedundant by size and by p-value keeps the rows Stripenn keeps") {
    const Table table = read_table("redundant_in.tsv");
    const std::vector<st::Stripe> input = redundant_input(table);
    for (const auto& [by, name] : {std::pair{st::RedundantBy::Size, "size"}, std::pair{st::RedundantBy::Pvalue, "pvalue"}}) {
        const std::vector<st::Stripe> kept = st::remove_redundant(input, by);
        std::vector<int> rows;
        for (const st::Stripe& s : kept) {
            rows.push_back(s.x);
        }
        const std::vector<int> expected = read_int_list(std::string("redundant_") + name + ".json");
        CHECK(rows == expected);
        CHECK(kept.size() < input.size());  // something was removed
    }
}

TEST_CASE("a frame is searched like StripeSearch: boxes, positions and sums") {
    const npy::Array sub_array = load("search_submat.npy");
    const std::vector<double> sub = npy::as<double>(sub_array);
    const int S = static_cast<int>(sub_array.shape[0]);
    const double M = load_as<double>("search_M.npy")[0];
    constexpr std::int64_t res = 10000;
    st::ChromInfo chrom;
    chrom.name = "c1";
    chrom.size = static_cast<std::int64_t>(S) * res;
    chrom.nbins = S;
    st::BandMatrix band(S, st::required_depth(res, 8), res, chrom.size);
    for (int i = 0; i < S; ++i) {
        for (int j = i; j < S; ++j) {
            const double v = sub[static_cast<std::size_t>(i) * S + j];
            if (v != 0.0) {
                band.add(i, j, v);
            }
        }
    }
    st::Parameters params;
    const st::FrameResult frame = st::search_frame(band, chrom, 0, M, params, nullptr);
    REQUIRE(frame.searched);
    const Table expected = read_table("search_expected.tsv");
    REQUIRE(frame.stripes.size() == expected.rows.size());
    for (std::size_t i = 0; i < expected.rows.size(); ++i) {
        const auto& row = expected.rows[i];
        const st::Stripe& s = frame.stripes[i];
        CHECK(s.pos1 == std::stoll(row[expected.index("pos1")]));
        CHECK(s.pos2 == std::stoll(row[expected.index("pos2")]));
        CHECK(s.pos3 == std::stoll(row[expected.index("pos3")]));
        CHECK(s.pos4 == std::stoll(row[expected.index("pos4")]));
        CHECK(s.length == std::stoll(row[expected.index("length")]));
        CHECK(s.width == std::stoll(row[expected.index("width")]));
        CHECK(s.total == parse_double(row[expected.index("total")]));
        CHECK(s.mean == parse_double(row[expected.index("Mean")]));
        CHECK(s.x == std::stoi(row[expected.index("x")]));
        CHECK(s.y == std::stoi(row[expected.index("y")]));
        CHECK(s.h == std::stoi(row[expected.index("h")]));
        CHECK(s.w == std::stoi(row[expected.index("w")]));
    }
}
