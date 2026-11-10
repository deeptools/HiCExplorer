// The computational core of hicDetectStripes, separated from the tool so that
// the unit tests (cpp/tests/test_detect_stripes.cpp) and the oracle
// comparison (cpp/scripts/stripenn_oracle) can reach it directly.
//
// hicDetectStripes implements the Stripenn approach (Yoon et al. 2022,
// Stripenn 1.1.65.22): the maximum-pixel saturation of each chromosome's
// contact image, the sliding frame scan with Canny edge detection and the
// vertical-line and block filters, the pairing of edge columns into stripe
// boxes, the removal of redundant boxes, the empirical background
// distribution and the p-value ranked against it, and the stripiness score.
// Every stage is a port of the corresponding stripenn function in the order
// of its floating point operations, and the acceptance gate is identical
// output to Stripenn on the same input (cpp/scripts/stripenn_oracle compares
// every stage). Where a stage cannot be bit identical the cause is stated at
// the stage.
//
// Stripenn's own structure is kept, including its quirks, because the output
// is defined by them: the null distribution is built from all chromosomes of
// the matrix (not only the requested ones), each chromosome's null sample
// starts from the seeded generator's initial state, the first percentile
// grid is the whole chromosome's positive pixels, and the zero rows and
// columns of a frame are removed before the search.

#ifndef HICX_TOOLS_DETECT_STRIPES_IMPL_HPP
#define HICX_TOOLS_DETECT_STRIPES_IMPL_HPP

#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "stripes_numpy.hpp"

namespace hicx::stripes {

// Stripenn's parameters with its defaults (cli.py).
struct Parameters {
    double canny = 2.0;
    int min_length = 10;  // minL, in bins
    int max_width = 8;    // maxW, in bins
    std::vector<double> maxpixel{0.95, 0.96, 0.97, 0.98, 0.99};
    double pvalue = 0.1;
    std::string mask = "0";
    int bfilter = 3;
    std::uint64_t seed = 123456789;
    // Stripenn run with one core draws every chromosome's null sample from one
    // generator in sequence; with several cores (the default) each chromosome
    // starts from the seeded initial state. See compute_null.
    bool chained_rng = false;
};

struct ChromInfo {
    std::string name;
    std::int64_t size = 0;       // bp
    std::int64_t first_bin = 0;  // index of the chromosome's first bin in the matrix
    std::int64_t nbins = 0;      // ceil(size / resolution)
};

// Python's int(x) of a non-negative double, and ceil division helpers on bp.
[[nodiscard]] inline std::int64_t floor_div(std::int64_t a, std::int64_t b) {
    std::int64_t q = a / b;
    if ((a % b != 0) && ((a < 0) != (b < 0))) {
        --q;
    }
    return q;
}
[[nodiscard]] inline std::int64_t ceil_div(std::int64_t a, std::int64_t b) { return -floor_div(-a, b); }

// Which rows of a chromosome carry a non-zero value in the wide windows that
// Stripenn's null distribution tests (see compute_null); filled while the
// pixels are read.
struct NullGeometry {
    int itera = 0;      // number of units the chromosome is cut into
    int unitsize = 0;   // bins per unit
};
[[nodiscard]] NullGeometry null_geometry(std::int64_t chrom_size, std::int64_t resolution);

// One chromosome's balanced contact values near the diagonal, as cooler's
// dense fetch of the balanced matrix returns them: NaN in every row and
// column of a masked bin, 0 for a pixel that is not stored. Only the upper
// triangle within `depth` bins of the diagonal is kept; asking for a value
// outside it throws.
class BandMatrix {
  public:
    BandMatrix() = default;
    BandMatrix(std::int64_t bins, std::int64_t depth, std::int64_t resolution,
               std::int64_t chrom_size);

    // Adds a stored pixel (bin1 <= bin2, bins relative to the chromosome).
    void add(std::int64_t bin1, std::int64_t bin2, double value);
    void set_bad(const std::vector<std::uint8_t>& bad_bins) { bad_ = bad_bins; }
    // After all pixels: nothing to do, kept so that loaders have a closing call.
    void finish() {}

    [[nodiscard]] std::int64_t bins() const noexcept { return bins_; }
    [[nodiscard]] std::int64_t depth() const noexcept { return depth_; }
    [[nodiscard]] bool bad(std::int64_t bin) const { return !bad_.empty() && bad_[static_cast<std::size_t>(bin)] != 0; }

    // The matrix entry (row, column) with cooler's NaN semantics.
    [[nodiscard]] double get(std::int64_t row, std::int64_t column) const;
    // The same with NaN replaced by 0 (stripenn's nantozero).
    [[nodiscard]] double get_zero(std::int64_t row, std::int64_t column) const {
        const double v = get(row, column);
        return v != v ? 0.0 : v;
    }
    // Cooler's fetch(rows, cols) as a dense row-major array (NaN kept).
    [[nodiscard]] std::vector<double> dense(std::int64_t row0, std::int64_t row1, std::int64_t col0,
                                            std::int64_t col1) const;

    // Row flags for compute_null / available_cols (see below).
    [[nodiscard]] const std::vector<std::uint8_t>& unit_flags() const { return unit_flags_; }
    [[nodiscard]] const std::vector<std::uint8_t>& window_flags() const { return window_flags_; }
    [[nodiscard]] std::int64_t resolution() const noexcept { return resolution_; }
    [[nodiscard]] std::int64_t chrom_size() const noexcept { return chrom_size_; }

  private:
    std::int64_t bins_ = 0;
    std::int64_t depth_ = 0;
    std::int64_t resolution_ = 0;
    std::int64_t chrom_size_ = 0;
    std::vector<double> values_;  // [bin1 * (depth + 1) + (bin2 - bin1)]
    std::vector<std::uint8_t> bad_;
    std::vector<std::uint8_t> unit_flags_;
    std::vector<std::uint8_t> window_flags_;
    NullGeometry geometry_;
    std::vector<std::int64_t> unit_row_end_;  // per unit, the last row (exclusive) of its fetch
};

// Band depth the whole method needs at a resolution (see stripes_impl.cpp).
[[nodiscard]] std::int64_t required_depth(std::int64_t resolution, int max_width);

// ---------------------------------------------------------------------------
// numpy.quantile over the positive entries of a whole symmetric chromosome
// matrix (getQuantile_original), in two passes over the stored upper triangle.
class QuantileSelector {
  public:
    explicit QuantileSelector(std::vector<double> quantiles);
    // Pass 1: every positive stored pixel; off-diagonal ones count twice.
    void count(double value, bool diagonal);
    void prepare();
    // Pass 2: the same pixels again.
    void gather(double value, bool diagonal);
    [[nodiscard]] std::vector<double> result();

  private:
    std::vector<double> quantiles_;
    std::vector<std::uint64_t> histogram_;
    std::uint64_t total_ = 0;
    std::vector<np::QuantileRanks> ranks_;
    std::vector<std::int64_t> wanted_buckets_;
    std::vector<std::vector<double>> bucket_values_;
    std::vector<std::int64_t> bucket_slot_;
    std::vector<std::uint64_t> bucket_start_;
};

// Stripenn.mpmean for one chromosome: the 400 expected values.
[[nodiscard]] std::vector<double> expected_values(const BandMatrix& band);

// ---------------------------------------------------------------------------
// The null distribution (getStripe.nulldist).

struct NullTables {
    // [d * columns + k] for distance d in 0..399 and sample k
    std::int64_t columns = 0;
    std::vector<double> leftup;
    std::vector<double> rightup;
    std::vector<double> leftdown;
    std::vector<double> rightdown;
    // per distance, the non-NaN values sorted (for rank queries)
    std::vector<std::vector<double>> sorted_leftup;
    std::vector<std::vector<double>> sorted_rightup;
    std::vector<std::vector<double>> sorted_leftdown;
    std::vector<std::vector<double>> sorted_rightdown;
    void build_sorted();
};

// nulldist's first stage: the number of rows with data in the chromosome's
// diagonal unit blocks (available_cols). Needs the unit flags of the band.
[[nodiscard]] std::int64_t available_columns(const BandMatrix& band);

// nulldist's allocation of the 1000 samples over the chromosomes, from the
// chromosome sizes and the available column counts of every chromosome of the
// matrix (in matrix order). Returns, per chromosome, its sample count, or -1
// for a chromosome that gets no samples; reproduces nulldist's bookkeeping
// including the way its two filtered lists are indexed.
[[nodiscard]] std::vector<std::int64_t> allocate_samples(
    const std::vector<ChromInfo>& chroms, const std::vector<std::int64_t>& available);

// main_null_calc for one chromosome with `samples` samples. `rng` is advanced.
struct NullColumns {
    std::int64_t columns = 0;
    std::vector<double> leftup;  // [d * columns + k]
    std::vector<double> rightup;
    std::vector<double> leftdown;
    std::vector<double> rightdown;
};
[[nodiscard]] NullColumns compute_null(const BandMatrix& band, std::int64_t samples, np::PyRandom& rng);

// Concatenates the chromosomes' columns in order (np.column_stack).
[[nodiscard]] NullTables join_null(const std::vector<NullColumns>& parts);

// ---------------------------------------------------------------------------
// The frame scan.

// One candidate stripe in Stripenn's result table.
struct Stripe {
    std::int32_t chrom = 0;  // index into the list of requested chromosomes
    std::int64_t pos1 = 0;
    std::int64_t pos2 = 0;
    std::int64_t pos3 = 0;
    std::int64_t pos4 = 0;
    std::int64_t length = 0;
    std::int64_t width = 0;
    double total = 0.0;
    double mean = 0.0;
    double maxpixel = 0.0;   // the percentile as a fraction (0.95)
    int perc_index = 0;
    std::int64_t num = 0;    // frame index
    std::int64_t start = 0;
    std::int64_t end = 0;
    int x = 0;
    int y = 0;
    int h = 0;
    int w = 0;
    double pvalue = std::numeric_limits<double>::quiet_NaN();
    double stripiness = std::numeric_limits<double>::quiet_NaN();
    double oe_mean = std::numeric_limits<double>::quiet_NaN();
};

// Optional capture of a frame's intermediates, for the oracle comparison.
struct FrameTrace {
    std::vector<float> gray;       // [6 * S * S]
    std::vector<std::uint8_t> edges;  // [6 * S * S]
    std::vector<std::uint8_t> vert;   // [6 * S * S]
    int size = 0;
};

// StripeSearch on one compacted frame: `submat` is S x S row major, `M` the
// chromosome's maxpixel value for the percentile. Returns the boxes of
// the six brightness settings in Stripenn's order, before RemoveRedundant.
struct RawStripe {
    int x = 0;
    int y = 0;
    int w = 0;
    int h = 0;
    double total = 0.0;
    double mean = 0.0;
};
[[nodiscard]] std::vector<RawStripe> stripe_search(const std::vector<double>& submat, int S, double M,
                                                   const Parameters& params, FrameTrace* trace);

// search_frame for frame `idx` of a chromosome and one maxpixel value: the
// boxes of stripe_search mapped to genomic positions, before RemoveRedundant.
// `frame_trace` receives the frame's intermediates when not null, and
// `compacted` the frame's start / end arrays (for the dump).
struct FrameResult {
    std::vector<Stripe> stripes;
    std::vector<double> submat;
    std::vector<std::int64_t> start_array;
    std::vector<std::int64_t> end_array;
    bool searched = false;
};
[[nodiscard]] FrameResult search_frame(const BandMatrix& band, const ChromInfo& chrom,
                                       std::int64_t idx, double M, const Parameters& params,
                                       FrameTrace* trace);

// RemoveRedundant(by='size'|'pvalue') over the candidates of one chromosome in
// table order; returns the kept ones in the same order.
enum class RedundantBy { Size, Pvalue };
[[nodiscard]] std::vector<Stripe> remove_redundant(const std::vector<Stripe>& stripes, RedundantBy by);

// getStripe.pvalue for the candidates of one chromosome.
void compute_pvalues(const BandMatrix& band, const ChromInfo& chrom, const NullTables& null,
                     std::vector<Stripe>& stripes, unsigned int threads);

// scoringstripes for the candidates of one chromosome (Stripiness and the
// O/E mean). `expected` is the chromosome's 400 expected values.
void compute_scores(const BandMatrix& band, const ChromInfo& chrom, const std::vector<double>& expected,
                    const Parameters& params, std::vector<Stripe>& stripes, unsigned int threads);

// Python's repr(float) (what pandas writes).
[[nodiscard]] std::string py_repr(double value);

}  // namespace hicx::stripes

#endif  // HICX_TOOLS_DETECT_STRIPES_IMPL_HPP
