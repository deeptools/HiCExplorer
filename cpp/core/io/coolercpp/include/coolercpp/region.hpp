// Genomic regions and cooler URIs (cooler.util.parse_region and friends).

#ifndef COOLERCPP_REGION_HPP
#define COOLERCPP_REGION_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "coolercpp/dtype.hpp"

namespace coolercpp {

// An ordered mapping of sequence names to lengths: the pandas Series
// Cooler.chromsizes returns. Lookup finds the first entry with a name.
class ChromSizes {
  public:
    ChromSizes() = default;
    ChromSizes(std::vector<std::string> names, std::vector<std::int64_t> lengths,
               DType dtype = DType::Int64);

    [[nodiscard]] std::size_t size() const noexcept { return names_.size(); }
    [[nodiscard]] const std::vector<std::string>& names() const noexcept { return names_; }
    [[nodiscard]] const std::vector<std::int64_t>& lengths() const noexcept { return lengths_; }
    // The dtype of the Series values (int32 for a cooler's chroms/length).
    [[nodiscard]] DType dtype() const noexcept { return dtype_; }
    [[nodiscard]] std::optional<std::int64_t> find(std::string_view name) const;
    // chromsizes[name]; throws KeyError.
    [[nodiscard]] std::int64_t operator[](std::string_view name) const;

  private:
    std::vector<std::string> names_;
    std::vector<std::int64_t> lengths_;
    DType dtype_ = DType::Int64;
};

// A region argument: a UCSC style string such as "chr3:10,000,000-12M", or a
// (chrom, start, end) triple whose start and end may be None.
class Region {
  public:
    Region(const char* text) : text_(text) {}               // NOLINT
    Region(std::string text) : text_(std::move(text)) {}     // NOLINT
    Region(std::string_view text) : text_(std::string(text)) {}  // NOLINT
    Region(std::string chrom, std::optional<std::int64_t> start,
           std::optional<std::int64_t> end)
        : chrom_(std::move(chrom)), start_(start), end_(end), is_tuple_(true) {}

    [[nodiscard]] bool is_string() const noexcept { return !is_tuple_; }
    [[nodiscard]] const std::string& text() const noexcept { return text_; }
    [[nodiscard]] const std::string& chrom() const noexcept { return chrom_; }
    [[nodiscard]] std::optional<std::int64_t> start() const noexcept { return start_; }
    [[nodiscard]] std::optional<std::int64_t> end() const noexcept { return end_; }

  private:
    std::string text_;
    std::string chrom_;
    std::optional<std::int64_t> start_;
    std::optional<std::int64_t> end_;
    bool is_tuple_ = false;
};

// (chrom, start or None, end or None)
struct GenomicRange {
    std::string chrom;
    std::optional<std::int64_t> start;
    std::optional<std::int64_t> end;
    friend bool operator==(const GenomicRange&, const GenomicRange&) = default;
};

// A well formed (chrom, start, end) triple.
struct RegionTuple {
    std::string chrom;
    std::int64_t start = 0;
    std::int64_t end = 0;
    friend bool operator==(const RegionTuple&, const RegionTuple&) = default;
};

// cooler.util.parse_cooler_uri: (file path, group path). More than one "::"
// raises ValueError.
[[nodiscard]] std::pair<std::string, std::string> parse_cooler_uri(std::string_view uri);

// cooler.util.parse_humanized: "10,000" -> 10000, "1.5M" -> 1500000.
[[nodiscard]] std::int64_t parse_humanized(std::string_view text);

// cooler.util.parse_region_string.
[[nodiscard]] GenomicRange parse_region_string(std::string_view text);

// cooler.util.parse_region.
[[nodiscard]] RegionTuple parse_region(const Region& region,
                                       const ChromSizes* chromsizes = nullptr);

}  // namespace coolercpp

#endif  // COOLERCPP_REGION_HPP
