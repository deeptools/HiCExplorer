// Contact input of hicDifferentialAnalysis: per chromosome pixel bands of every
// sample, the shared bin mask, the per-sample distance decay, and the two
// calibration devices (replicate split and planted differences).
//
// Counts are the raw `count` column of each cooler, which must hold
// non-negative integers; a balancing weight column is not applied, because the
// model works on counts and takes library size and distance decay as offsets.
//
// Memory. Only the pixels within `band` bins of the diagonal are kept, as a
// CSR matrix of 32-bit column indices and 32-bit counts (8 bytes a pixel),
// built directly from the stored pixel order; the coverage of every bin is
// accumulated from all pixels while they stream past.
//
// Randomness. The replicate split and the planted differences draw binomial
// variates from a counter-based generator keyed by the seed, a stream label,
// and the global bin ids of the pixel, so a pixel's draw depends on nothing
// else: not on the thread count, the chunk size of the reader, or the order in
// which chromosomes are processed.

#ifndef HICX_TOOLS_DIFF_CONTACTS_IMPL_HPP
#define HICX_TOOLS_DIFF_CONTACTS_IMPL_HPP

#include <cstddef>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

#include "hicx/cool_adapter.hpp"

namespace hicx::diffc {

// --------------------------------------------------------------------------
// Randomness

[[nodiscard]] std::uint64_t mix64(std::uint64_t z) noexcept;

// splitmix64 over a state derived from the keys.
class CounterRng {
  public:
    CounterRng(std::initializer_list<std::uint64_t> keys) noexcept;
    [[nodiscard]] std::uint64_t next() noexcept;
    // Uniform on [0, 1), 53 bits.
    [[nodiscard]] double uniform() noexcept;

  private:
    std::uint64_t state_ = 0;
};

// An exact Binomial(n, p) draw: Bernoulli trials for n <= 32, otherwise
// inversion that accumulates the probability mass outwards from the mode.
[[nodiscard]] std::int64_t binomial(std::int64_t n, double p, CounterRng& rng);

inline constexpr std::uint64_t kSplitStream = 0x53504c4954ULL;  // "SPLIT"
inline constexpr std::uint64_t kPlantStream = 0x504c414e54ULL;  // "PLANT"

// --------------------------------------------------------------------------
// Samples

// One sample on one chromosome: the stored pixels with col - row <= band,
// in CSR form over chromosome-local bins, plus the cis coverage of every bin
// (row plus column sums without the diagonal, over all distances).
struct ChromosomeData {
    std::int64_t first_bin = 0;
    std::int64_t bins = 0;
    std::vector<std::int64_t> row_start;
    std::vector<std::int32_t> col;
    std::vector<std::uint32_t> count;
    std::vector<double> coverage;

    // [begin, end) positions of row i's pixels with col in [col_first, col_last).
    [[nodiscard]] std::pair<std::size_t, std::size_t> row_range(std::int64_t i,
                                                                std::int64_t col_first,
                                                                std::int64_t col_last) const;
};

struct LoadRequest {
    // For messages.
    std::string path;
    // Replicate split for calibration: every count is divided into
    // Binomial(count, 1/2) and its complement, keyed by (split_seed, index,
    // bin1, bin2), and both halves are returned.
    bool split = false;
    std::uint64_t split_seed = 0;
    std::size_t index = 0;
};

// Reads bins [first, last) of `file` once: one sample, or with a split its
// two complementary halves. Throws on a count that is not a non-negative
// integer below 2^32, and on pixels that are not sorted by (bin1, bin2), as
// the cooler format requires.
[[nodiscard]] std::vector<ChromosomeData> load_chromosome(const hicx::CoolFile& file,
                                                         std::int64_t first, std::int64_t last,
                                                         std::int64_t band,
                                                         const LoadRequest& request);

// hicCorrectMatrix's MAD outlier rule on one coverage vector: median over the
// positive values, the median absolute deviation over all of them, modified
// z-score 0.6745 (x - median) / MAD. A bin is invalid when its coverage is 0 or
// its z-score lies outside [lower, upper]. 1 = invalid.
[[nodiscard]] std::vector<char> coverage_outliers(const std::vector<double>& coverage,
                                                  double lower, double upper);

// Mean count at each distance 0..band over the pixel pairs whose two bins are
// valid; 0 where no such pair exists.
[[nodiscard]] std::vector<double> distance_decay(const ChromosomeData& data,
                                                 const std::vector<char>& invalid,
                                                 std::int64_t band);

// The number of valid bin pairs (i, i + d) at each distance d in 0..band.
[[nodiscard]] std::vector<double> valid_pairs(const std::vector<char>& invalid, std::int64_t band);

// Planted difference: every stored pixel (i, j) (local bins, i <= j) with
// select(i, j) is replaced by Binomial(count, keep), keyed by (seed,
// kPlantStream, stream, global bin1, global bin2).
void thin_pixels(ChromosomeData& data, const std::function<bool(std::int64_t, std::int64_t)>& select,
                 double keep, std::uint64_t seed, std::uint64_t stream);

// --------------------------------------------------------------------------
// Text input

// Whitespace separated fields of every line that is not empty and does not
// start with '#', "track" or "browser".
[[nodiscard]] std::vector<std::vector<std::string>> read_fields(const std::string& path);

// numpy.median of a copy.
[[nodiscard]] double median(std::vector<double> values);

}  // namespace hicx::diffc

#endif  // HICX_TOOLS_DIFF_CONTACTS_IMPL_HPP
