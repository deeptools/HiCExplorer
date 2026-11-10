// Reader for the 4DN .pairs format, the input of hicBuildMatrix --pairsFile
// (cpp/PLAN.md tier 9, item 9.2). HiCExplorer 3.7 has no such reader.
//
// The format (4DN pairs specification, v1.0; pairtools writes the same):
//
//   ## pairs format v1.0
//   #sorted: chr1-chr2-pos1-pos2
//   #shape: upper triangle
//   #genome_assembly: mm10
//   #chromsize: chr1 195471971
//   #columns: readID chr1 pos1 chr2 pos2 strand1 strand2 pair_type mapq1 mapq2
//   .	chr1	3000000	chr1	5324596	-	+	UU	60	60
//
// Data lines are tab separated, positions are one-based 5' ends. Files are
// plain, gzip or bgzip compressed; htslib detects which.

#ifndef HICX_PAIRS_FILE_HPP
#define HICX_PAIRS_FILE_HPP

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace hicx {

struct PairsHeader {
    // The text after "## pairs format v" on the first line, which is required.
    std::string format_version;
    // "#sorted:" split at '-', "none" when the file says so, empty when absent.
    std::vector<std::string> sorted_keys;
    std::string shape;  // "#shape:", for example "upper triangle"
    std::vector<std::pair<std::string, std::int64_t>> chrom_sizes;  // "#chromsize:" in order
    std::string genome_assembly;                                    // "#genome_assembly:"
    // "#columns:", or the seven reserved columns of the specification
    // (readID chrom1 pos1 chrom2 pos2 strand1 strand2) when the line is absent.
    // The names chr1 and chr2, which the specification used before it settled
    // on chrom1 and chrom2 and which ENCODE's files still carry, are stored as
    // chrom1 and chrom2; the same holds for the keys of "#sorted:".
    std::vector<std::string> columns;
    bool columns_declared = false;

    // Field indices of the columns a matrix build reads, -1 when absent.
    int chrom1 = -1;
    int pos1 = -1;
    int chrom2 = -1;
    int pos2 = -1;
    int strand1 = -1;
    int strand2 = -1;
    int mapq1 = -1;
    int mapq2 = -1;
    int pair_type = -1;
    // One more than the largest of those indices: a data line needs at least
    // this many fields.
    std::size_t fields_needed = 0;

    // "#sorted: chr1-chr2-pos1-pos2" (or chrom1-chrom2-pos1-pos2).
    [[nodiscard]] bool sorted_by_chroms_then_positions() const;
    // "#shape: upper triangle".
    [[nodiscard]] bool upper_triangle() const;
};

// Parses the header lines, every line up to the first data line. Throws
// std::runtime_error when the first line is not "## pairs format v...", when a
// "#chromsize:" line is malformed or repeats a chromosome, when a column is
// named twice, or when chrom1, pos1, chrom2 or pos2 is missing.
[[nodiscard]] PairsHeader parse_pairs_header(const std::vector<std::string>& lines);

// What a pair's own annotations say about it, before any filter of the tool.
enum class PairKind : std::uint8_t {
    Mapped,           // usable as far as the file says
    Unmapped,         // a side is unmapped: chromosome "!", or a pair_type letter
                      // other than U, R and M (N null, X corrupt, W walk)
    NotUnique,        // a pair_type letter M (multi-mapping), no unmapped side
    MarkedDuplicate,  // pair_type DD, as pairtools dedup --mark-dups writes
};

// One data line, as views into the line.
struct PairsRecord {
    std::string_view chrom1;
    std::string_view chrom2;
    std::int64_t pos1 = 0;
    std::int64_t pos2 = 0;
    std::int32_t mapq1 = -1;  // -1 without the column
    std::int32_t mapq2 = -1;
    std::int8_t strand1 = 0;  // +1 for '+', -1 for '-', 0 otherwise or without the column
    std::int8_t strand2 = 0;
    PairKind kind = PairKind::Mapped;
};

// Parses one data line with the column layout of `header`. `fields` is
// scratch space the caller keeps between calls. Returns the reason when the
// line cannot be parsed (too few fields, a position or mapping quality that is
// not a non-negative integer), std::nullopt otherwise.
[[nodiscard]] std::optional<std::string> parse_pairs_record(
    std::string_view line, const PairsHeader& header, std::vector<std::string_view>& fields,
    PairsRecord& out);

class PairsReader {
  public:
    // Opens the file and reads its header. `threads` above 1 decompresses
    // bgzip blocks on that many threads, which does not change what is read.
    PairsReader(const std::string& path, int threads);
    ~PairsReader();
    PairsReader(const PairsReader&) = delete;
    PairsReader& operator=(const PairsReader&) = delete;

    [[nodiscard]] const PairsHeader& header() const noexcept { return header_; }

    // The next data line, without its line break (and without a trailing
    // carriage return); false at the end of the file. The view is valid until
    // the next call. Empty lines and later lines starting with '#' are
    // skipped, as cooler cload pairs skips comment lines.
    bool next(std::string_view& line);

    // The one-based line number of the line last returned by next().
    [[nodiscard]] std::int64_t line_number() const noexcept { return line_number_; }
    [[nodiscard]] const std::string& path() const noexcept { return path_; }

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    PairsHeader header_;
    std::string path_;
    std::int64_t line_number_ = 0;
};

}  // namespace hicx

#endif  // HICX_PAIRS_FILE_HPP
