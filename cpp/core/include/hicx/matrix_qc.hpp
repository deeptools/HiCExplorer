// Quality-control statistics of a contact matrix (cool / mcool): the
// computation behind hicMatrixQC.
//
// The statistics are the ones a BAM based QC (hicBuildMatrix' QC.log) cannot
// give once only the matrix is left, plus the QC.log counts themselves when
// hicBuildMatrix stored them in the cool metadata ("hicexplorer_qc"). Every
// number is computed in one streamed pass over the stored pixels (the upper
// triangle, in chunks) with per-thread accumulators; the only objects the
// size of the genome are a few vectors with one double per bin. A matrix is
// never densified.
//
// Definitions that the report repeats:
//   * a contact is a unit of a stored pixel count; the diagonal pixel (i, i)
//     is one pair of bins, so every stored pixel counts once;
//   * the distance of a pixel (i, j), i <= j, is the offset j - i in bins;
//     contacts at offset d are said to be at the distance d * resolution, the
//     distance between the starts of the two bins. A pixel at offset d holds
//     read pairs whose real distance lies between (d - 1) and (d + 1) bin
//     sizes, so every distance threshold is exact only up to one bin size;
//   * the diagonal (offset 0) holds the pairs inside one bin, all closer than
//     one bin size. They are short range whenever the cut-off exceeds 0;
//   * the coverage of a bin is the sum of its row of the symmetric matrix
//     (the diagonal entry counted once).

#ifndef HICX_MATRIX_QC_HPP
#define HICX_MATRIX_QC_HPP

#include <array>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "hicx/cool_adapter.hpp"
#include "hicx/json_lite.hpp"

namespace hicx::matrixqc {

inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

struct Options {
    // Intra-chromosomal contacts closer than this (bp) are short range.
    std::int64_t short_range_cutoff = 20000;
    // The decay exponent is the slope of log10 P(s) against log10 s for the
    // distances s in [fit_min, fit_max] (bp), P(s) taken in ten log-spaced
    // groups of distances.
    std::int64_t fit_min = 100000;
    std::int64_t fit_max = 1000000;
    // Map resolution (Rao et al. 2014): the smallest bin size at which
    // `rao_fraction` of the bins hold at least `rao_min_contacts` contacts.
    double rao_min_contacts = 1000.0;
    double rao_fraction = 0.8;
    // hicCorrectMatrix' modified z-score filter (--filterThreshold).
    double mad_lower = -1.5;
    double mad_upper = 5.0;
    // Chromosome pairs with observed / expected trans contacts above this are
    // flagged.
    double trans_threshold = 2.0;
    // The longest chromosomes that get their own rows (decay, coverage strip,
    // trans enrichment); the genome-wide statistics always use all.
    std::size_t max_chromosomes = 30;
    // Cells of the coverage strip of a chromosome.
    std::size_t strip_cells = 40;
    unsigned threads = 4;
    // The bin column of the balancing weights.
    std::string weight_column = "weight";
};

// Bins of a fixed-size matrix: chromosome c holds the bins
// [first_bin[c], first_bin[c + 1]).
struct Layout {
    std::int64_t resolution = 0;
    std::vector<std::string> names;
    std::vector<std::int64_t> lengths;
    std::vector<std::int64_t> first_bin;  // names.size() + 1 entries

    [[nodiscard]] std::int64_t bins() const { return first_bin.empty() ? 0 : first_bin.back(); }
    [[nodiscard]] std::size_t chromosomes() const { return names.size(); }
};

// The balancing factors a pixel's count is multiplied (op '*') or divided
// (op '/') by to give the balanced count; NaN marks a masked bin.
struct Correction {
    std::vector<double> factors;
    char op = '*';
};

// Delivers every stored pixel exactly once, in chunks, calling
// visit(thread, chunk) with 0 <= thread < threads; calls with one thread
// index never overlap.
using ChunkVisitor = std::function<void(unsigned thread, const PixelChunk& chunk)>;
using PixelSource = std::function<void(unsigned threads, const ChunkVisitor& visit)>;

// --- results ------------------------------------------------------------------

// Contacts per distance on the grid of the QC tables of hicBuildMatrix: slot 0
// is offset 0 (the diagonal), slot k >= 1 is the QC distance bin k - 1 (ten per
// decade, qc_distance_edge). Slots with no bin pair are left out by the
// callers.
inline constexpr std::size_t kDecaySlots = 101;
[[nodiscard]] std::int64_t decay_slot_lower_edge(std::size_t slot);  // bp, slot >= 1

struct DecayCurve {
    std::vector<double> contacts;       // per slot
    std::vector<double> pairs;          // bin pairs per slot
    std::vector<double> mean_distance;  // pair-weighted mean distance in bp
    double exponent = kNaN;             // slope of log10 P(s) against log10 s
    double intercept = kNaN;            // log10 P at s = 1 bp
    std::size_t fit_points = 0;         // groups of distances in the fit
};

struct Quantiles {
    double min = kNaN, p01 = kNaN, p05 = kNaN, p25 = kNaN, median = kNaN, p75 = kNaN,
           p95 = kNaN, p99 = kNaN, max = kNaN;
};

struct ChromosomeSummary {
    std::string name;
    std::int64_t length = 0;
    std::int64_t bins = 0;
    double cis_contacts = 0.0;
    double coverage_sum = 0.0;
    double empty_fraction = kNaN;
    double median_coverage = kNaN;  // of the non-empty bins
    double mad_removed_fraction = kNaN;
    DecayCurve decay;
    // Mean coverage of equal parts of the chromosome relative to the median
    // coverage of the non-empty bins of the genome; -1 for a part with no
    // contacts.
    std::vector<double> strip;
};

struct LadderStep {
    std::int64_t bin_size = 0;
    double fraction_all = kNaN;       // of all bins
    double fraction_nonempty = kNaN;  // of the bins with at least one contact
    bool reached = false;             // fraction_all >= rao_fraction
    bool reached_nonempty = false;    // fraction_nonempty >= rao_fraction
};

struct BalanceInfo {
    bool has_weight_column = false;
    bool usable = false;  // a weight column with at least one finite weight
    double nan_fraction = kNaN;
    // The attributes of the weight dataset (converged, ignore_diags, ...)
    // as text, in name order.
    std::map<std::string, std::string> attributes;
    bool has_marginals = false;
    double marginal_mean = kNaN;
    double marginal_cv = kNaN;
};

struct TransTable {
    std::vector<std::string> names;  // the chromosomes of the table
    // k x k, row major; the diagonal is NaN in expected.
    std::vector<double> observed;
    std::vector<double> expected;
    double total = 0.0;  // trans contacts between the listed chromosomes
    struct Pair {
        std::size_t a = 0, b = 0;
        double observed = 0, expected = 0, ratio = 0;
    };
    std::vector<Pair> flagged;  // a < b, ratio above the threshold, largest first
};

struct CheapChecks {
    double diagonal_fraction = kNaN;       // offset 0 / all contacts
    double first_offdiag_fraction = kNaN;  // offset 1 / all contacts
    bool has_x = false;
    double x_autosome_ratio = kNaN;  // mean coverage per non-empty bin
    bool has_y = false;
    double y_fraction = kNaN;  // of the coverage
    bool has_m = false;
    double m_fraction = kNaN;
};

// The QC.log of hicBuildMatrix, read back from the cool metadata.
struct ReadQc {
    std::string tool;
    std::string tool_version;
    std::optional<std::int64_t> min_rest_site_distance;
    std::int64_t max_library_insert_size = 0;
    std::vector<std::string> restriction_sequences;
    std::vector<std::string> dangling_sequences;
    std::int64_t short_range_cutoff = 20000;
    bool keep_self_ligation = false;
    bool restriction_cut_file = false;

    std::int64_t sequenced_reads = 0, mappable = 0, hic_contacts = 0, one_mate_unmapped = 0,
                 one_mate_not_unique = 0, low_mapping_quality = 0;
    std::vector<std::int64_t> dangling_end;
    std::int64_t self_ligation = 0, not_close_to_rest_site = 0, same_fragment = 0,
                 self_circle = 0, duplicated_pairs = 0, inter_chromosomal = 0,
                 intra_short = 0, intra_long = 0;
    std::array<std::int64_t, 4> orientation_pairs{};  // inward, outward, left, right

    std::vector<std::int64_t> lower_edges;  // bp
    std::vector<std::int64_t> distance_contacts;
    std::array<std::vector<std::int64_t>, 4> orientation_distance;
};

[[nodiscard]] std::optional<ReadQc> parse_read_qc(const json::Value& qc);

struct MatrixStats {
    std::string label;  // the sample name of the report
    std::string uri;
    std::int64_t resolution = 0;  // 0 when the bins are not of a fixed size
    std::int64_t bins = 0;
    std::int64_t pixels = 0;
    std::size_t chromosomes_total = 0;
    // Pixel statistics need fixed-size bins; otherwise only read_qc is set.
    bool pixel_statistics = false;
    std::string note;

    double total_contacts = 0, cis = 0, trans = 0, intra_short = 0, intra_long = 0;
    std::int64_t short_range_cutoff = 0;
    double diagonal = 0, first_offdiagonal = 0;

    DecayCurve decay;
    std::vector<ChromosomeSummary> chromosomes;  // the listed (longest) ones

    Quantiles coverage;  // of the non-empty bins
    double empty_fraction = kNaN;
    double mad_median = kNaN, mad_value = kNaN;
    std::int64_t mad_removed = 0, mad_removed_empty = 0;
    double mad_removed_fraction = kNaN;
    double coverage_median_nonempty = kNaN;

    // Bins that hold >= rao_min_contacts contacts at this resolution.
    double rao_fraction_all = kNaN, rao_fraction_nonempty = kNaN;
    // This matrix' coverage aggregated to multiples of its bin size.
    std::vector<LadderStep> ladder;

    BalanceInfo balance;
    TransTable trans_table;
    CheapChecks checks;
    std::optional<ReadQc> read_qc;
};

struct MapResolution {
    // "resolutions" (the resolutions present in an mcool) or "aggregation"
    // (the bins of a single cool aggregated to multiples of its bin size).
    std::string method;
    std::vector<LadderStep> steps;
    // 0: not reached by any step.
    std::int64_t bin_size = 0;
    std::int64_t bin_size_nonempty = 0;
    // The smallest bin size tried already reaches the fraction: the map
    // resolution is then that bin size or finer, and the file cannot tell.
    bool finest_reached = false;
};

struct FileStats {
    std::string path;
    bool multires = false;
    std::vector<MatrixStats> matrices;  // one per resolution of an mcool
    MapResolution map_resolution;
};

// --- computation --------------------------------------------------------------

// The statistics of one matrix from a stream of its pixels. `correction` (the
// balancing weights) is optional.
[[nodiscard]] MatrixStats compute_matrix_stats(const Layout& layout, const PixelSource& source,
                                               const std::optional<Correction>& correction,
                                               const Options& options);

// Fraction of the bins of the aggregation of `coverage` to multiples of
// `factor` bins that hold >= min_contacts, over all bins and over the
// non-empty bins. Chromosomes are aggregated separately; the last aggregated
// bin of a chromosome may be shorter.
[[nodiscard]] LadderStep aggregate_ladder_step(const std::vector<double>& coverage,
                                               const Layout& layout, std::int64_t factor,
                                               const Options& options);

// The multiples of the bin size tried when a single cool is aggregated.
[[nodiscard]] std::vector<std::int64_t> ladder_factors(const Layout& layout);

// Ordinary least squares of y on x; false with fewer than 3 points.
[[nodiscard]] bool fit_line(const std::vector<double>& x, const std::vector<double>& y,
                            double& slope, double& intercept);

// X, Y, M and autosome names: "chrX", "X", "chrM", "MT", "chr12", "2L".
enum class ChromosomeClass { Autosome, X, Y, M, Other };
[[nodiscard]] ChromosomeClass classify_chromosome(const std::string& name);

// The layout of a fixed-size cool; throws std::runtime_error for variable
// bins (resolution 0 is never returned).
[[nodiscard]] std::optional<Layout> layout_of(const CoolFile& cool);

// The pixels of a cool (raw counts), read in parallel by row ranges.
[[nodiscard]] PixelSource cool_pixel_source(const CoolFile& cool, const Layout& layout);

// The resolutions of an mcool: 0 entries for a plain cool. Throws when the
// file is not a cooler.
[[nodiscard]] std::vector<std::pair<std::int64_t, std::string>> cool_uris(const std::string& path);

// Everything for one file: all its resolutions (or those listed), and the map
// resolution. `label` names the sample.
[[nodiscard]] FileStats analyse_file(const std::string& path, const std::string& label,
                                     const std::vector<std::int64_t>& resolutions,
                                     const Options& options);

// The map resolution from the ladder steps; steps are ordered by bin size.
[[nodiscard]] MapResolution pick_map_resolution(std::string method, std::vector<LadderStep> steps);

}  // namespace hicx::matrixqc

#endif  // HICX_MATRIX_QC_HPP
