// Bin and region model.
//
// The Python side keeps the genomic axis of a matrix as a list of
// "cut intervals", tuples (chromosome, start, end, extra), plus two derived
// structures built by hicmatrix.HiCMatrix.intervalListToIntervalTree:
//
//   * one interval tree per chromosome, mapping a genomic position to a bin id
//   * chrBinBoundaries, an ordered mapping chromosome -> [first bin, last bin)
//
// BinTable holds all three. The quirks of the Python implementation are
// reproduced on purpose, including the fact that a chromosome that reappears
// after a different chromosome resets its interval tree and overwrites its
// boundary entry.

#ifndef HICX_BINS_HPP
#define HICX_BINS_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace hicx {

struct CutInterval {
    std::string chrom;
    std::int64_t start = 0;
    std::int64_t end = 0;
    // The fourth element of the Python tuple. It is 1.0 for cool files and the
    // /intervals/extra_list value for h5 files. Some h5 matrices written by
    // hicFindTADs store text there instead of a number, in which case extra is
    // NaN and extra_text holds the value.
    double extra = 1.0;
    std::string extra_text;

    friend bool operator==(const CutInterval& a, const CutInterval& b) {
        return a.chrom == b.chrom && a.start == b.start && a.end == b.end &&
               a.extra == b.extra && a.extra_text == b.extra_text;
    }
};

struct BinRange {
    std::int64_t first = 0;  // inclusive
    std::int64_t last = 0;   // exclusive

    // chrBinBoundaries is compared for equality by hicSumMatrices and
    // hicCompareMatrices, so the ranges need a value comparison.
    friend bool operator==(const BinRange& a, const BinRange& b) {
        return a.first == b.first && a.last == b.last;
    }
};

class BinTable {
  public:
    BinTable() = default;
    explicit BinTable(std::vector<CutInterval> intervals);

    [[nodiscard]] const std::vector<CutInterval>& intervals() const noexcept {
        return intervals_;
    }
    [[nodiscard]] std::size_t size() const noexcept { return intervals_.size(); }
    [[nodiscard]] bool empty() const noexcept { return intervals_.empty(); }
    [[nodiscard]] const CutInterval& bin_pos(std::size_t bin_id) const;

    // chrBinBoundaries, in insertion order.
    [[nodiscard]] const std::vector<std::pair<std::string, BinRange>>&
    chrom_bin_boundaries() const noexcept {
        return boundaries_;
    }
    [[nodiscard]] std::vector<std::string> chrom_names() const;
    [[nodiscard]] std::optional<BinRange> chrom_bin_range(const std::string& chrom) const;

    // get_chromosome_sizes: end coordinate of the last bin of each chromosome.
    [[nodiscard]] std::vector<std::pair<std::string, std::int64_t>>
    chromosome_sizes() const;
    // get_chromosome_sizes_real: end of the last bin minus start of the first
    // bin plus one.
    [[nodiscard]] std::vector<std::pair<std::string, std::int64_t>>
    chromosome_sizes_real() const;

    // getRegionBinRange: first bin overlapping startpos and first bin
    // overlapping endpos. Empty when either position is outside the bins of
    // the chromosome, mirroring the Python IndexError path that returns None.
    [[nodiscard]] std::optional<std::pair<std::int64_t, std::int64_t>>
    region_bin_range(const std::string& chrom, std::int64_t startpos,
                     std::int64_t endpos) const;

    // Bin id covering a single genomic position, or nullopt.
    [[nodiscard]] std::optional<std::int64_t> bin_at(const std::string& chrom,
                                                     std::int64_t position) const;

    // getBinSize: the median of the differences between consecutive bin starts
    // over all chromosomes that hold more than one bin, truncated to an
    // integer. Matches hicmatrix.HiCMatrix.getBinSize including its handling
    // of a single bin.
    [[nodiscard]] std::int64_t bin_size() const;
    // False when more than one percent of the bins deviate from the median.
    [[nodiscard]] bool bin_size_homogeneous() const;

  private:
    void build_index();

    struct BinInterval {
        std::int64_t start = 0;
        std::int64_t end = 0;
        std::int64_t bin_id = 0;
    };

    std::vector<CutInterval> intervals_;
    std::vector<std::pair<std::string, BinRange>> boundaries_;
    std::unordered_map<std::string, std::size_t> boundary_index_;
    std::unordered_map<std::string, std::vector<BinInterval>> trees_;
    mutable std::optional<std::int64_t> bin_size_;
    mutable bool bin_size_homogeneous_ = true;
};

}  // namespace hicx

#endif  // HICX_BINS_HPP
