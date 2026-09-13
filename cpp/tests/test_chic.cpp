// Unit tests for the cHi-C core components: hicx/chic_viewpoint.hpp,
// hicx/chic_hdf5.hpp and hicx/lbfgsb_scipy.hpp.
//
// Expected values were produced by the reference environment (numpy 1.26.4,
// scipy 1.14.1, the repository's lib/viewpoint.py) on the cHi-C test data.

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include <hdf5.h>

#include "doctest/doctest.h"
#include "hicx/chic_hdf5.hpp"
#include "hicx/chic_viewpoint.hpp"
#include "hicx/lbfgsb_scipy.hpp"
#include "hicx/scipy_special.hpp"
#include "hicx/stats_ops.hpp"

namespace {

const std::string kChic = std::string(HICX_TEST_DATA_DIR) + "/cHi-C/";

std::string temporary_path(const std::string& name) {
    return (std::filesystem::temp_directory_path() / ("hicx_test_chic_" + name)).string();
}

}  // namespace

TEST_CASE("chic: Python text semantics") {
    using namespace hicx::chic;
    CHECK(python_int(" 42\n") == 42);
    CHECK(python_int("-1_000") == -1000);
    CHECK_THROWS_AS((void)python_int("4.5"), ValueError);
    CHECK_THROWS_AS((void)python_int(""), ValueError);
    CHECK(python_float("0.000105073709\n") == 0.000105073709);
    CHECK(std::isnan(python_float("nan")));
    CHECK_THROWS_AS((void)python_float("x"), ValueError);
    CHECK(format_fixed(58.254577634279, 12) == "58.254577634279");
    CHECK(format_fixed(std::numeric_limits<double>::quiet_NaN(), 12) == "nan");
    CHECK(format_fixed(-std::numeric_limits<double>::quiet_NaN(), 12) == "nan");
    CHECK(format_fixed(2.5e-13, 12) == "0.000000000000");
    CHECK(strip(" \tchr1\t1\t2\r\n") == "chr1\t1\t2");
}

TEST_CASE("chic: reference point files") {
    const auto points = hicx::chic::read_reference_points(kChic + "referencePoints_edge_cases.bed");
    REQUIRE(points.points.size() == 7);
    CHECK(points.genes[4] == "ThreeColumns");
    // A three field line is (chrom, start, start).
    CHECK(points.points[4].start == "197195000");
    CHECK(points.points[4].end == "197195000");
    CHECK(hicx::chic::reference_point_string(points.points[1]) == "chr1_24202042_24205042");
    // A blank last line is skipped by the reader.
    const auto three = hicx::chic::read_reference_points(kChic + "referencePoints_chicViewpoint.bed");
    CHECK(three.points.size() == 3);
    CHECK(hicx::chic::read_lines(kChic + "referencePoints_chicViewpoint.bed").size() == 4);
}

TEST_CASE("chic: smoothing keeps the right border asymmetry of even windows") {
    std::vector<double> data;
    for (int i = 1; i <= 21; ++i) {
        data.push_back(std::pow(static_cast<double>(i), 1.5));
    }
    const auto mean = [&](std::size_t begin, std::size_t end) {
        double sum = 0.0;
        for (std::size_t i = begin; i < end; ++i) {
            sum += data[i];
        }
        return sum / static_cast<double>(end - begin);
    };
    for (std::int64_t window = 1; window <= 6; ++window) {
        const std::vector<double> smoothed =
            hicx::chic::smooth_interaction_values(std::span<const double>(data), window);
        const std::int64_t half = window / 2;
        const std::int64_t upstream = window % 2 == 0 ? half - 1 : half;
        const std::int64_t n = static_cast<std::int64_t>(data.size());
        for (std::int64_t i = 0; i < n; ++i) {
            std::int64_t begin = 0;
            std::int64_t end = 0;
            if (i >= n - half) {
                begin = i - half;
                end = n;
            } else {
                begin = std::max<std::int64_t>(i - upstream, 0);
                end = i + half + 1;
            }
            // Sequential and pairwise sums agree for at most seven terms.
            CHECK(smoothed[static_cast<std::size_t>(i)] ==
                  doctest::Approx(mean(static_cast<std::size_t>(begin),
                                       static_cast<std::size_t>(end))).epsilon(1e-15));
        }
    }
    // float32 means are single precision quantities widened afterwards.
    const std::vector<float> single{0.0F, 1.0F, 0.0F, 0.0F, 0.0F};
    const std::vector<double> smoothed =
        hicx::chic::smooth_interaction_values(std::span<const float>(single), 5);
    CHECK(smoothed[2] == static_cast<double>(1.0F / 5.0F));
    CHECK(smoothed[2] != 0.2);
}

TEST_CASE("chic: viewpoints on the real cHi-C matrix") {
    using namespace hicx::chic;
    const ViewpointMatrix matrix = ViewpointMatrix::load(kChic + "FL-E13-5_chr1.cool");
    CHECK(matrix.bin_size() == 1000);
    const ReferencePoint sox17{"chr1", "4487435", "4487435"};
    const ViewpointRange range = calculate_viewpoint_range(matrix, sox17, 500000, 500000);
    CHECK(range.region_start == 3987435);
    CHECK(range.region_end == 4987435);
    const ComputedViewpoint viewpoint =
        compute_viewpoint(matrix, sox17, "chr1", range.region_start, range.region_end);
    REQUIRE(viewpoint.data.size() == 1001);
    CHECK(viewpoint.index_before_viewpoint == 500);
    std::size_t nonzero = 0;
    double total = 0.0;
    for (const double value : viewpoint.data) {
        nonzero += value != 0.0 ? 1 : 0;
        total += value;
    }
    // chicQualityControl's raw filter: 0.14485514485514486
    CHECK(static_cast<double>(nonzero) / 1001.0 == 0.14485514485514486);
    // chicViewpoint's sum_of_interactions for Sox17
    CHECK(total == 810.0);

    // Clipped at the chromosome end: (197195432 - 197145200) + 1000
    const ViewpointRange near_end =
        calculate_viewpoint_range(matrix, ReferencePoint{"chr1", "197145200", "197145200"},
                                  200000, 200000);
    CHECK(near_end.region_end == 197195431);
    CHECK(near_end.downstream == 51232);

    // A position past the chromosome is a TypeError, a start after the end a
    // ValueError (numpy cannot broadcast), an unknown chromosome a ValueError.
    const ReferencePoint beyond{"chr1", "794979306", "794979306"};
    CHECK_THROWS_AS((void)compute_viewpoint(matrix, beyond, "chr1", 3987435, 4987435), TypeError);
    const ReferencePoint inverted{"chr1", "70931031", "70093803"};
    const ViewpointRange inverted_range = calculate_viewpoint_range(matrix, inverted, 500000, 500000);
    CHECK_THROWS_AS((void)compute_viewpoint(matrix, inverted, "chr1", inverted_range.region_start,
                                            inverted_range.region_end),
                    ValueError);
    CHECK_THROWS_AS((void)calculate_viewpoint_range(matrix, ReferencePoint{"chr2", "1", "1"}, 10, 10),
                    ValueError);
}

TEST_CASE("chic: relative values treat a zero denominator as absent") {
    const std::vector<double> data{1.0, 3.0};
    CHECK(hicx::chic::compute_relative_values(data, 8.0)[1] == 0.375);
    CHECK(hicx::chic::compute_relative_values(data, 0.0)[1] == 0.75);
}

TEST_CASE("chic: background model extension and the p-value lookup rule") {
    using namespace hicx::chic;
    const BackgroundModel model =
        read_background_model(kChic + "background.txt", 600000, 700000, 500000, false);
    CHECK(model.min_key() == -600000);
    CHECK(model.max_key() == 700000);
    CHECK(model.at(700000) == model.at(500000));
    CHECK(model.at(-500000)[0] == 58.254577634279);
    const BackgroundModel means =
        read_background_model(kChic + "background.txt", 200000, 200000, 500000, true);
    CHECK(interaction_background_data(means, 200000, 200000).size() == 401);

    const BackgroundModel plain =
        read_background_model(kChic + "background.txt", 200000, 200000, 500000, false);
    const std::vector<double> raw{9.8, 0.0, 26.6};
    const std::vector<double> p = p_values(plain, raw, 1);
    // Index 0 is one bin upstream, which is not a key in base pairs, so the
    // model at -500 kb is used; scipy gives 1 - betainc = 0 exactly here.
    CHECK(p[0] == 1.0 - hicx::scipy::betainc(58.254577634279, 10.8, 0.998604709282));
    CHECK(p[0] == 0.0);
    CHECK(p[1] == 1.0);
    const std::vector<double>& upstream = plain.at(plain.max_key());
    CHECK(p[2] == 1.0 - hicx::scipy::betainc(upstream[0], 27.6, upstream[1]));
}

TEST_CASE("scipy::betainc is Boost's ibeta, exact where cephes is one ulp short") {
    // scipy.special.betainc(58.254577634279, x + 1, 0.998604709282) is 1.0 for
    // all four; cephes incbet returns 1 - MACHEP there.
    for (const double x : {9.8, 26.6, 75.0, 90.0}) {
        CHECK(hicx::scipy::betainc(58.254577634279, x + 1.0, 0.998604709282) == 1.0);
        CHECK(hicx::stats::betainc(58.254577634279, x + 1.0, 0.998604709282) ==
              0.9999999999999999);
    }
    // Values from scipy 1.14.1, compared bit for bit.
    CHECK(hicx::scipy::betainc(2.5, 3.5, 0.5) == 0.6697652726313549);
    CHECK(hicx::scipy::betainc(0.10000000000000001, 0.20000000000000001, 0.5) ==
          0.6705707961028995);
    CHECK(std::isnan(hicx::scipy::betainc(0.0, 1.0, 0.5)));
    CHECK(std::isnan(hicx::scipy::betainc(1.0, 1.0, 1.5)));
    CHECK(hicx::scipy::betainc(2.0, 3.0, 0.0) == 0.0);
    CHECK(hicx::scipy::betainc(2.0, 3.0, 1.0) == 1.0);
}

TEST_CASE("lbfgsb_scipy follows scipy's iterations on smooth problems") {
    using hicx::stats::Bound;
    const double inf = std::numeric_limits<double>::infinity();
    const double eps = std::numeric_limits<double>::epsilon();
    const std::function<double(std::span<const double>)> rosen = [](std::span<const double> x) {
        return (1.0 - x[0]) * (1.0 - x[0]) + 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
    };
    const std::function<double(std::span<const double>)> quad = [](std::span<const double> x) {
        return (x[0] - 3.0) * (x[0] - 3.0) + 10.0 * (x[1] + 0.5) * (x[1] + 0.5) + x[0] * x[1];
    };
    struct Expected {
        const std::function<double(std::span<const double>)>* f;
        std::vector<double> x0;
        std::vector<Bound> bounds;
        double x;
        double y;
        int iterations;
        int evaluations;
    };
    // fmin_l_bfgs_b(f, x0, approx_grad=1, bounds=...): x, nit, funcalls
    const std::vector<Expected> cases{
        {&rosen, {-1.2, 1.0}, {{-inf, inf}, {-inf, inf}}, 0.999996066816395, 0.9999922573958389, 36, 132},
        {&rosen, {-1.2, 1.0}, {{-2.0, 0.8}, {-1.0, 0.5}}, 0.7085594986832366, 0.5, 9, 33},
        {&rosen, {0.3, 2.0}, {{0.5, inf}, {eps, 1.0}}, 0.9999999950019915, 1.0, 14, 54},
        {&quad, {0.0, 0.0}, {{eps, inf}, {eps, 1.0}}, 2.999999858829152, 2.220446049250313e-16, 2, 9},
        {&quad, {5.0, 0.9}, {{-inf, 2.0}, {0.0, 1.0}}, 2.0, 0.0, 1, 6},
    };
    for (const Expected& expected : cases) {
        const auto result = hicx::stats::minimise_lbfgsb_scipy(*expected.f, expected.x0, expected.bounds);
        CHECK(result.status == 0);
        CHECK(result.iterations == expected.iterations);
        CHECK(result.function_evaluations == expected.evaluations);
        // Same iterations; the iterates agree to the level scipy's OpenBLAS
        // kernels and libm leave room for.
        CHECK(result.x[0] == doctest::Approx(expected.x).epsilon(2e-7));
        CHECK(result.x[1] == doctest::Approx(expected.y).epsilon(2e-7));
    }
}

TEST_CASE("fit_nbinom_scipy returns the Python's answer where there is nothing to fit") {
    const std::vector<double> empty;
    const auto none = hicx::stats::fit_nbinom_scipy(empty);
    CHECK(none.size == 10.0);
    CHECK(std::isnan(none.prob));
    // factorial overflows above 170: scipy returns its moment start.
    const std::vector<double> overflow{177.0, 0.0, 2.0, 5.0};
    const auto start = hicx::stats::fit_nbinom_scipy(overflow);
    const double mean = 46.0;
    const double variance = (131.0 * 131.0 + 46.0 * 46.0 + 44.0 * 44.0 + 41.0 * 41.0) / 4.0;
    const double size = (mean * mean) / (variance - mean);
    CHECK(start.size == doctest::Approx(size).epsilon(1e-15));
    CHECK(start.prob == doctest::Approx(size / (size + mean)).epsilon(1e-15));
}

TEST_CASE("chic_hdf5 writes h5py's layout") {
    const std::string path = temporary_path("layout.hdf5");
    {
        hicx::chic::Hdf5Writer writer(path);
        writer.set_attribute("/", "type", std::string("interactions"));
        const std::vector<std::int64_t> range{200000, 100000};
        writer.set_attribute("/", "range", std::span<const std::int64_t>(range));
        writer.create_group("m/chr1/g");
        CHECK_THROWS((void)writer.create_group("m/chr1/g"));
        hicx::chic::InteractionFileData data;
        data.chromosome = "chr1";
        data.gene = "g";
        data.sum_of_interactions = 3.0;
        data.starts = {0, 1000};
        data.ends = {1000, 2000};
        data.relative_positions = {0, 1000};
        data.interaction_data = {0.5, 0.25};
        data.raw = {1.5, 0.75};
        data.pvalues = {1.0, 0.5};
        data.xfold = {2.0, 4.0};
        hicx::chic::write_interaction_datasets(writer, "m/chr1/g", data, 5, 6);
        CHECK(writer.hard_link("m/chr1/g", "m/genes/g"));
        CHECK_FALSE(writer.hard_link("m/chr1/g", "m/genes/g"));
    }
    const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    REQUIRE(file >= 0);
    const hid_t dataset = H5Dopen2(file, "m/genes/g/pvalue", H5P_DEFAULT);
    REQUIRE(dataset >= 0);
    const hid_t plist = H5Dget_create_plist(dataset);
    CHECK(H5Pget_layout(plist) == H5D_CHUNKED);
    CHECK(H5Pget_nfilters(plist) == 1);
    unsigned int flags = 0;
    std::size_t count = 1;
    unsigned int level = 0;
    char name[32];
    CHECK(H5Pget_filter2(plist, 0, &flags, &count, &level, sizeof(name), name, nullptr) ==
          H5Z_FILTER_DEFLATE);
    CHECK(level == 9);
    H5Pclose(plist);
    H5Dclose(dataset);
    const hid_t gene = H5Dopen2(file, "m/chr1/g/gene", H5P_DEFAULT);
    const hid_t type = H5Dget_type(gene);
    CHECK(H5Tis_variable_str(type) > 0);
    CHECK(H5Tget_cset(type) == H5T_CSET_UTF8);
    H5Tclose(type);
    H5Dclose(gene);
    H5Fclose(file);
    std::remove(path.c_str());
}
