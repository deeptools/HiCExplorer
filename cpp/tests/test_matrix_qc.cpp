// Unit tests of the hicMatrixQC statistics on small synthetic matrices whose
// answers are worked out by hand in the comments, and on the shared
// gm12878_chr1.cool.

#include <doctest/doctest.h>

#include <cmath>
#include <string>
#include <tuple>
#include <vector>

#include "hicx/cool_adapter.hpp"
#include "hicx/matrix_qc.hpp"

namespace qc = hicx::matrixqc;

namespace {

using Pixel = std::tuple<std::int64_t, std::int64_t, double>;

qc::Layout make_layout(std::int64_t resolution,
                       const std::vector<std::pair<std::string, std::int64_t>>& chromosomes) {
    qc::Layout layout;
    layout.resolution = resolution;
    layout.first_bin.push_back(0);
    for (const auto& [name, bins] : chromosomes) {
        layout.names.push_back(name);
        layout.lengths.push_back(bins * resolution);
        layout.first_bin.push_back(layout.first_bin.back() + bins);
    }
    return layout;
}

// The pixels in one chunk per worker; `workers` > 1 deals them out in turn so
// that the merge of the per-thread sums is exercised.
qc::PixelSource source_of(const std::vector<Pixel>& pixels, unsigned workers = 1) {
    return [pixels, workers](unsigned threads, const qc::ChunkVisitor& visit) {
        REQUIRE(threads >= workers);
        std::vector<hicx::PixelChunk> chunks(workers);
        for (std::size_t i = 0; i < pixels.size(); ++i) {
            hicx::PixelChunk& chunk = chunks[i % workers];
            chunk.bin1.push_back(std::get<0>(pixels[i]));
            chunk.bin2.push_back(std::get<1>(pixels[i]));
            chunk.count.push_back(std::get<2>(pixels[i]));
        }
        for (unsigned w = 0; w < workers; ++w) {
            visit(w, chunks[w]);
        }
    };
}

qc::Options small_options() {
    qc::Options options;
    options.threads = 2;
    options.rao_min_contacts = 7;
    options.fit_min = 20000;
    options.fit_max = 80000;
    return options;
}

// Two chromosomes, A with 10 bins (0..9) and B with 6 bins (10..15), 10 kb bins.
// pixels: diagonal (0,0)=5 (1,1)=4 (9,9)=1 (10,10)=6; (0,1)=3; (2,5)=2 at offset 3;
// (10,12)=2 at offset 2 on B; trans (3,12)=7.
const std::vector<Pixel> kTwoChromosomes = {{0, 0, 5}, {0, 1, 3}, {1, 1, 4}, {2, 5, 2},
                                            {9, 9, 1}, {3, 12, 7}, {10, 10, 6}, {10, 12, 2}};

}  // namespace

TEST_CASE("totals, cis and trans, short and long range, diagonal") {
    const qc::Layout layout = make_layout(10000, {{"chrA", 10}, {"chrB", 6}});
    for (const unsigned workers : {1U, 2U}) {
        const qc::MatrixStats s =
            qc::compute_matrix_stats(layout, source_of(kTwoChromosomes, workers), std::nullopt, small_options());
        CHECK(s.total_contacts == 30);
        CHECK(s.pixels == 8);
        CHECK(s.cis == 23);
        CHECK(s.trans == 7);
        // short range is below 20 kb: offsets 0 and 1 (19), offsets 2 and 3 are long (2 + 2)
        CHECK(s.intra_short == 19);
        CHECK(s.intra_long == 4);
        CHECK(s.diagonal == 16);
        CHECK(s.first_offdiagonal == 3);
        CHECK(s.checks.diagonal_fraction == doctest::Approx(16.0 / 30.0));
        CHECK(s.checks.first_offdiag_fraction == doctest::Approx(3.0 / 30.0));
    }
}

TEST_CASE("the short range cut-off moves with the option") {
    const qc::Layout layout = make_layout(10000, {{"chrA", 10}, {"chrB", 6}});
    qc::Options options = small_options();
    options.short_range_cutoff = 10000;  // only the diagonal is closer than one bin size
    const qc::MatrixStats s = qc::compute_matrix_stats(layout, source_of(kTwoChromosomes), std::nullopt, options);
    CHECK(s.intra_short == 16);
    CHECK(s.intra_long == 7);
}

TEST_CASE("coverage, empty bins and quantiles") {
    const qc::Layout layout = make_layout(10000, {{"chrA", 10}, {"chrB", 6}});
    const qc::MatrixStats s =
        qc::compute_matrix_stats(layout, source_of(kTwoChromosomes), std::nullopt, small_options());
    // coverage per bin: [8,7,2,7,0,2,0,0,0,1, 8,0,9,0,0,0]: 8 of 16 bins are empty
    CHECK(s.empty_fraction == doctest::Approx(0.5));
    // non-empty, sorted: 1 2 2 7 7 8 8 9
    CHECK(s.coverage.min == 1);
    CHECK(s.coverage.max == 9);
    CHECK(s.coverage.median == doctest::Approx(7.0));
    CHECK(s.coverage.p25 == doctest::Approx(2.0));
    CHECK(s.coverage.p75 == doctest::Approx(8.0));
}

TEST_CASE("the MAD filter of hicCorrectMatrix") {
    // 8 bins, coverage without the diagonal 11 10 13 13 14 15 16 16: median
    // 13.5, MAD 2.0 (absolute deviations .5 .5 .5 1.5 2.5 2.5 2.5 3.5)
    const qc::Layout layout = make_layout(1000, {{"chr1", 8}});
    const std::vector<Pixel> pixels = {{0, 1, 10}, {2, 3, 12}, {4, 5, 14}, {6, 7, 16}, {0, 2, 1}, {3, 5, 1}};
    // coverage: b0 = 11, b1 = 10, b2 = 13, b3 = 13, b4 = 14, b5 = 15, b6 = 16, b7 = 16
    qc::Options options;
    options.threads = 1;
    options.mad_lower = -1.0;
    options.mad_upper = 5.0;
    const qc::MatrixStats s = qc::compute_matrix_stats(layout, source_of(pixels), std::nullopt, options);
    CHECK(s.mad_median == doctest::Approx(13.5));
    CHECK(s.mad_value == doctest::Approx(2.0));
    // z(b1) = 0.6745 * (10 - 13.5) / 2 = -1.18 and z(b0) = 0.6745 * (11 - 13.5) / 2 = -0.84:
    // only bin 1 leaves [-1, 5]
    CHECK(s.mad_removed == 1);
    CHECK(s.mad_removed_fraction == doctest::Approx(1.0 / 8.0));
    CHECK(s.mad_removed_empty == 0);
    CHECK(s.chromosomes.at(0).mad_removed_fraction == doctest::Approx(1.0 / 8.0));
}

TEST_CASE("the contact distance curve, per slot") {
    const qc::Layout layout = make_layout(10000, {{"chrA", 10}, {"chrB", 6}});
    const qc::MatrixStats s =
        qc::compute_matrix_stats(layout, source_of(kTwoChromosomes), std::nullopt, small_options());
    // slot 0: offset 0 (16 contacts, 10 + 6 bins); the slot of a distance is
    // 1 + the QC distance bin: 10 kb in bin 40, 20 kb in 43, 30 kb in 44
    CHECK(s.decay.contacts[0] == 16);
    CHECK(s.decay.pairs[0] == 16);
    CHECK(s.decay.contacts[41] == 3);
    CHECK(s.decay.pairs[41] == 9 + 5);
    CHECK(s.decay.contacts[44] == 2);
    CHECK(s.decay.pairs[44] == 8 + 4);
    CHECK(s.decay.contacts[45] == 2);
    CHECK(s.decay.pairs[45] == 7 + 3);
    CHECK(s.decay.mean_distance[41] == doctest::Approx(10000.0));
    CHECK(s.decay.mean_distance[45] == doctest::Approx(30000.0));
    // chrA alone: its slot 41 holds the (0,1) pixel
    REQUIRE(s.chromosomes.size() == 2);
    CHECK(s.chromosomes[0].name == "chrA");
    CHECK(s.chromosomes[0].decay.contacts[41] == 3);
    CHECK(s.chromosomes[0].decay.pairs[41] == 9);
    CHECK(s.chromosomes[1].decay.contacts[44] == 2);
    CHECK(s.chromosomes[0].cis_contacts == 15);
    CHECK(s.chromosomes[1].cis_contacts == 8);
}

TEST_CASE("the decay exponent of a power law is its exponent") {
    // one chromosome of 2000 bins of 1 kb; the contacts at offset d are
    // chosen so that P(d) = contacts / (n - d) = 1e6 / d
    const std::int64_t n = 2000;
    const qc::Layout layout = make_layout(1000, {{"chr1", n}});
    std::vector<Pixel> pixels;
    for (std::int64_t d = 1; d < n; ++d) {
        pixels.emplace_back(0, d, 1e6 * static_cast<double>(n - d) / static_cast<double>(d));
    }
    qc::Options options;
    options.threads = 1;
    options.fit_min = 100000;
    options.fit_max = 1000000;
    const qc::MatrixStats s = qc::compute_matrix_stats(layout, source_of(pixels), std::nullopt, options);
    CHECK(s.decay.fit_points == 10);
    CHECK(s.decay.exponent == doctest::Approx(-1.0).epsilon(0.02));
    CHECK(s.chromosomes.at(0).decay.exponent == doctest::Approx(-1.0).epsilon(0.02));
    // and a power law with exponent -1.5
    std::vector<Pixel> steeper;
    for (std::int64_t d = 1; d < n; ++d) {
        steeper.emplace_back(0, d, 1e8 * static_cast<double>(n - d) / std::pow(static_cast<double>(d), 1.5));
    }
    const qc::MatrixStats t = qc::compute_matrix_stats(layout, source_of(steeper), std::nullopt, options);
    CHECK(t.decay.exponent == doctest::Approx(-1.5).epsilon(0.02));
}

TEST_CASE("the decay exponent needs at least three groups of distances") {
    const qc::Layout layout = make_layout(500000, {{"chr1", 20}});
    qc::Options options;
    options.threads = 1;  // offsets 1 and 2 only lie in 100 kb to 1 Mb
    const qc::MatrixStats s =
        qc::compute_matrix_stats(layout, source_of({{0, 1, 5}, {0, 2, 3}}), std::nullopt, options);
    CHECK(std::isnan(s.decay.exponent));
}

TEST_CASE("fit_line is ordinary least squares") {
    double slope = 0, intercept = 0;
    CHECK(qc::fit_line({0, 1, 2, 3}, {1, 3, 5, 7}, slope, intercept));
    CHECK(slope == doctest::Approx(2.0));
    CHECK(intercept == doctest::Approx(1.0));
    CHECK_FALSE(qc::fit_line({0, 1}, {1, 3}, slope, intercept));
    CHECK_FALSE(qc::fit_line({1, 1, 1}, {1, 3, 5}, slope, intercept));
}

TEST_CASE("map resolution: bins aggregated to multiples of the bin size") {
    const qc::Layout layout = make_layout(10000, {{"chrA", 10}, {"chrB", 6}});
    qc::Options options = small_options();  // 7 contacts
    options.rao_fraction = 0.5;
    const qc::MatrixStats s =
        qc::compute_matrix_stats(layout, source_of(kTwoChromosomes), std::nullopt, options);
    // coverage [8,7,2,7,0,2,0,0,0,1, 8,0,9,0,0,0]: 5 bins hold 7 or more of 16 (8 non-empty)
    CHECK(s.rao_fraction_all == doctest::Approx(5.0 / 16.0));
    CHECK(s.rao_fraction_nonempty == doctest::Approx(5.0 / 8.0));
    // factor 2 within each chromosome: A 15 9 2 0 1, B 8 9 0: four of 8 groups, of 6 non-empty
    const qc::LadderStep step = qc::aggregate_ladder_step(
        std::vector<double>{8, 7, 2, 7, 0, 2, 0, 0, 0, 1, 8, 0, 9, 0, 0, 0}, layout, 2, options);
    CHECK(step.bin_size == 20000);
    CHECK(step.fraction_all == doctest::Approx(0.5));
    CHECK(step.fraction_nonempty == doctest::Approx(4.0 / 6.0));
    CHECK(step.reached);
    CHECK_FALSE(qc::aggregate_ladder_step(std::vector<double>{8, 7, 2, 7, 0, 2, 0, 0, 0, 1, 8, 0, 9, 0, 0, 0},
                                           layout, 1, options).reached);

    std::vector<qc::LadderStep> steps = {s.ladder.at(0), step};
    const qc::MapResolution map = qc::pick_map_resolution("aggregation", steps);
    CHECK(map.bin_size == 20000);
    CHECK(map.bin_size_nonempty == 10000);  // 5 of the 8 non-empty bins are enough already
    CHECK(qc::pick_map_resolution("aggregation", {s.ladder.at(0)}).bin_size == 0);
}

TEST_CASE("trans contacts: observed over expected per chromosome pair") {
    // three chromosomes of two 1 kb bins: A (0,1), B (2,3), C (4,5)
    const qc::Layout layout = make_layout(1000, {{"A", 2}, {"B", 2}, {"C", 2}});
    // trans contacts: A-B 6, A-C 2, B-C 2; marginals T_A = 8, T_B = 8, T_C = 4, S = 20,
    // total 10; f = 0.4 0.4 0.2, sum of squares 0.36, so E_AB = 10 * 2 * 0.16 / 0.64 = 5,
    // E_AC = E_BC = 10 * 2 * 0.08 / 0.64 = 2.5
    const std::vector<Pixel> pixels = {{0, 2, 4}, {1, 3, 2}, {0, 4, 2}, {2, 5, 2}, {0, 0, 9}};
    qc::Options options;
    options.threads = 1;
    options.trans_threshold = 1.1;
    const qc::MatrixStats s = qc::compute_matrix_stats(layout, source_of(pixels), std::nullopt, options);
    const qc::TransTable& t = s.trans_table;
    REQUIRE(t.names.size() == 3);
    CHECK(t.total == 10);
    CHECK(t.observed[0 * 3 + 1] == 6);
    CHECK(t.observed[1 * 3 + 0] == 6);
    CHECK(t.expected[0 * 3 + 1] == doctest::Approx(5.0));
    CHECK(t.expected[0 * 3 + 2] == doctest::Approx(2.5));
    CHECK(t.expected[1 * 3 + 2] == doctest::Approx(2.5));
    // ratios: A-B 1.2, A-C 0.8, B-C 0.8; only A-B exceeds 1.1
    REQUIRE(t.flagged.size() == 1);
    CHECK(t.flagged[0].a == 0);
    CHECK(t.flagged[0].b == 1);
    CHECK(t.flagged[0].ratio == doctest::Approx(1.2));
    // the expected values add up to the observed total
    CHECK(t.expected[0 * 3 + 1] + t.expected[0 * 3 + 2] + t.expected[1 * 3 + 2] == doctest::Approx(10.0));
}

TEST_CASE("only the longest chromosomes are listed") {
    const qc::Layout layout = make_layout(1000, {{"short", 2}, {"long1", 6}, {"long2", 5}});
    qc::Options options;
    options.threads = 1;
    options.max_chromosomes = 2;
    const qc::MatrixStats s = qc::compute_matrix_stats(layout, source_of({{0, 1, 1}, {2, 3, 4}, {8, 9, 3}}),
                                                       std::nullopt, options);
    REQUIRE(s.chromosomes.size() == 2);
    CHECK(s.chromosomes[0].name == "long1");  // in genome order
    CHECK(s.chromosomes[1].name == "long2");
    CHECK(s.total_contacts == 8);  // the genome-wide statistics still use all
}

TEST_CASE("chromosome classes and the X, Y and M checks") {
    CHECK(qc::classify_chromosome("chrX") == qc::ChromosomeClass::X);
    CHECK(qc::classify_chromosome("Y") == qc::ChromosomeClass::Y);
    CHECK(qc::classify_chromosome("chrM") == qc::ChromosomeClass::M);
    CHECK(qc::classify_chromosome("MT") == qc::ChromosomeClass::M);
    CHECK(qc::classify_chromosome("chr12") == qc::ChromosomeClass::Autosome);
    CHECK(qc::classify_chromosome("2L") == qc::ChromosomeClass::Autosome);
    CHECK(qc::classify_chromosome("chrUn_KI270302v1") == qc::ChromosomeClass::Other);

    const qc::Layout layout = make_layout(1000, {{"chr1", 4}, {"chrX", 4}, {"chrY", 2}, {"chrM", 1}});
    std::vector<Pixel> pixels;
    for (int bin = 0; bin < 4; ++bin) {
        pixels.emplace_back(bin, bin, 10);
    }
    for (int bin = 4; bin < 8; ++bin) {
        pixels.emplace_back(bin, bin, 5);
    }
    pixels.emplace_back(10, 10, 10);
    qc::Options options;
    options.threads = 1;
    const qc::MatrixStats s = qc::compute_matrix_stats(layout, source_of(pixels), std::nullopt, options);
    CHECK(s.checks.has_x);
    CHECK(s.checks.x_autosome_ratio == doctest::Approx(0.5));
    CHECK(s.checks.has_y);
    CHECK(s.checks.y_fraction == doctest::Approx(0.0));
    CHECK(s.checks.has_m);
    CHECK(s.checks.m_fraction == doctest::Approx(10.0 / 70.0));
}

TEST_CASE("the balancing state: NaN weights and the marginals") {
    const qc::Layout layout = make_layout(10000, {{"chrA", 10}, {"chrB", 6}});
    qc::Correction correction;
    correction.factors.assign(16, 1.0);
    correction.factors[4] = std::nan("");
    const qc::MatrixStats s =
        qc::compute_matrix_stats(layout, source_of(kTwoChromosomes), correction, small_options());
    CHECK(s.balance.usable);
    CHECK(s.balance.nan_fraction == doctest::Approx(1.0 / 16.0));
    // weights of 1 leave the counts, so the marginals are the coverage of the 15 finite bins:
    // sum 44 over 15 bins
    CHECK(s.balance.marginal_mean == doctest::Approx(44.0 / 15.0));
    // sum of squares 64 + 49 + 4 + 49 + 4 + 1 + 64 + 81 = 316
    const double mean = 44.0 / 15.0;
    CHECK(s.balance.marginal_cv == doctest::Approx(std::sqrt(316.0 / 15.0 - mean * mean) / mean));

    // balanced counts c * w_i * w_j with w = 2 everywhere: every marginal is 4 times the coverage
    qc::Correction doubled;
    doubled.factors.assign(16, 2.0);
    const qc::MatrixStats d =
        qc::compute_matrix_stats(layout, source_of(kTwoChromosomes), doubled, small_options());
    CHECK(d.balance.marginal_mean == doctest::Approx(4.0 * 44.0 / 16.0));
    // divisive weights divide by w_i * w_j
    doubled.op = '/';
    const qc::MatrixStats v =
        qc::compute_matrix_stats(layout, source_of(kTwoChromosomes), doubled, small_options());
    CHECK(v.balance.marginal_mean == doctest::Approx(0.25 * 44.0 / 16.0));
    // no correction: no marginals
    const qc::MatrixStats none =
        qc::compute_matrix_stats(layout, source_of(kTwoChromosomes), std::nullopt, small_options());
    CHECK_FALSE(none.balance.usable);
    CHECK_FALSE(none.balance.has_marginals);
}

TEST_CASE("the QC statistics of hicBuildMatrix are read back from the metadata") {
    const std::string text =
        "{\"schema_version\": 1, \"tool\": \"hicBuildMatrix\", \"tool_version\": \"4.0.0\", "
        "\"parameters\": {\"min_rest_site_distance\": 300, \"max_library_insert_size\": 1000, "
        "\"restriction_sequences\": [\"GATC\"], \"dangling_sequences\": [\"GATC\"], "
        "\"short_range_cutoff\": 20000, \"keep_self_ligation\": false, \"restriction_cut_file\": true}, "
        "\"counts\": {\"sequenced_reads\": 100, \"pairs_mappable_unique_high_quality\": 70, "
        "\"hic_contacts\": 30, \"one_mate_unmapped\": 10, \"one_mate_not_unique\": 5, "
        "\"low_mapping_quality\": 15, \"dangling_end\": [{\"restriction_sequence\": \"GATC\", "
        "\"dangling_sequence\": \"GATC\", \"count\": 4}], \"self_ligation\": 3, "
        "\"one_mate_not_close_to_rest_site\": 0, \"same_fragment\": 6, \"self_circle\": 2, "
        "\"duplicated_pairs\": 1, \"inter_chromosomal\": 10, \"intra_short_range\": 8, "
        "\"intra_long_range\": 12, \"inward_pairs\": 5, \"outward_pairs\": 6, \"left_pairs\": 7, "
        "\"right_pairs\": 8}, \"distance\": {\"bins_per_decade\": 10, \"lower_edges_bp\": [1, 2, 2], "
        "\"contacts\": [4, 0, 16], \"orientation\": {\"inward\": [1, 0, 2], \"outward\": [1, 0, 3], "
        "\"left\": [1, 0, 5], \"right\": [1, 0, 6]}}}";
    const auto parsed = hicx::json::parse(text);
    REQUIRE(parsed.has_value());
    const auto q = qc::parse_read_qc(*parsed);
    REQUIRE(q.has_value());
    CHECK(q->sequenced_reads == 100);
    CHECK(q->mappable == 70);
    CHECK(q->hic_contacts == 30);
    CHECK(q->dangling_end == std::vector<std::int64_t>{4});
    CHECK(q->restriction_sequences == std::vector<std::string>{"GATC"});
    CHECK(q->min_rest_site_distance.value() == 300);
    CHECK(q->restriction_cut_file);
    CHECK(q->orientation_pairs[3] == 8);
    CHECK(q->distance_contacts == std::vector<std::int64_t>{4, 0, 16});
    CHECK(q->orientation_distance[2] == std::vector<std::int64_t>{1, 0, 5});
    CHECK_FALSE(qc::parse_read_qc(*hicx::json::parse("{\"a\": 1}")).has_value());
}

TEST_CASE("gm12878_chr1.cool: the pixel pass agrees with the file and across threads") {
    const std::string path = std::string(HICX_TEST_DATA_DIR) + "/gm12878_chr1.cool";
    const hicx::CoolFile cool(path);
    const std::optional<qc::Layout> layout = qc::layout_of(cool);
    REQUIRE(layout.has_value());
    qc::Options options;
    options.threads = 1;
    const qc::MatrixStats one = qc::compute_matrix_stats(*layout, qc::cool_pixel_source(cool, *layout), std::nullopt, options);
    options.threads = 4;
    const qc::MatrixStats four = qc::compute_matrix_stats(*layout, qc::cool_pixel_source(cool, *layout), std::nullopt, options);
    CHECK(one.pixels == cool.nnz());
    CHECK(one.total_contacts == cool.info_value("sum")->as_double());
    CHECK(four.total_contacts == one.total_contacts);
    CHECK(four.cis == one.cis);
    CHECK(four.intra_short == one.intra_short);
    CHECK(four.decay.exponent == one.decay.exponent);
    CHECK(four.mad_removed == one.mad_removed);
    CHECK(one.trans == 0);  // a single chromosome
    CHECK(one.decay.exponent < -0.5);
    CHECK(one.decay.exponent > -2.0);
}
