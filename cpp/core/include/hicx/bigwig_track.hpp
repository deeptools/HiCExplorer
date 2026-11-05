// A bigWig file read as pyBigWig reads it, for the bigwig tracks of
// hicPlotMatrix's plotBigwig (Python API, hicexplorer._hicx.hicPlotMatrix).
// It is libBigWig, the library pyBigWig wraps: chroms() and
// stats(chrom, start, end, nBins=n) with the default type 'mean'.

#ifndef HICX_BIGWIG_TRACK_HPP
#define HICX_BIGWIG_TRACK_HPP

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace hicx::bigwig {

class Track {
  public:
    // pyBigWig.open(path); throws std::runtime_error("Received an error
    // during file opening!") when libBigWig cannot open it.
    explicit Track(const std::string& path);
    ~Track();
    Track(const Track&) = delete;
    Track& operator=(const Track&) = delete;

    // chroms(): name and length, in the file's order
    [[nodiscard]] const std::vector<std::pair<std::string, std::int64_t>>& chroms() const noexcept {
        return chroms_;
    }
    [[nodiscard]] std::optional<std::int64_t> length(const std::string& chrom) const;

    // stats(chrom, start, end, nBins=bins): the mean of every bin, NaN where
    // pyBigWig returns None. A bin count below 1 is 1, as in pyBigWig. Throws
    // std::invalid_argument("Invalid interval bounds!") for an unknown
    // chromosome, an end past it or an empty interval.
    [[nodiscard]] std::vector<double> mean(const std::string& chrom, std::int64_t start, std::int64_t end,
                                           std::int64_t bins) const;

  private:
    void* file_ = nullptr;
    std::vector<std::pair<std::string, std::int64_t>> chroms_;
};

}  // namespace hicx::bigwig

#endif  // HICX_BIGWIG_TRACK_HPP
