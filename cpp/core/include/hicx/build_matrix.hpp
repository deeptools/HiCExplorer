// The machinery behind hicBuildMatrix, hicBuildMatrixMicroC and hicQuickQC:
// hicexplorer/lib/buildMatrixMethods.py.
//
// Everything here reproduces the Python literally, including the places where
// the Python is wrong. Each of those is marked QUIRK with the line it comes
// from, and the list is repeated in the report. The port never silently
// improves on the reference.
//
// Memory (cpp/PLAN.md 2.1). Chromosome names are interned once and every
// interval carries a 32 bit id instead of a std::string. The DpnII cut site
// file for dm3 holds 599,308 intervals; as Python tuples of (str, int, int)
// inside an intervaltree that is hundreds of megabytes, here it is 16 bytes
// each.

#ifndef HICX_BUILD_MATRIX_HPP
#define HICX_BUILD_MATRIX_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace hicx {

// A genomic interval with an interned chromosome name.
struct GenomeInterval {
    std::uint32_t chrom = 0;
    std::int64_t start = 0;
    std::int64_t end = 0;
};

// Interning table for chromosome names. Ids are assigned in order of first
// appearance, which is the order every Python dict over chromosomes has.
class ChromNames {
  public:
    std::uint32_t intern(const std::string& name);
    [[nodiscard]] std::optional<std::uint32_t> lookup(const std::string& name) const;
    [[nodiscard]] const std::string& name(std::uint32_t id) const {
        return names_[id];
    }
    [[nodiscard]] std::size_t size() const noexcept { return names_.size(); }

  private:
    std::vector<std::string> names_;
    std::unordered_map<std::string, std::uint32_t> index_;
};

using ChromSizes = std::vector<std::pair<std::string, std::int64_t>>;

// hicexplorer.utilities.getUserRegion. Throws std::runtime_error with the
// Python's message for an unknown chromosome, which is the NameError that
// hicexplorer/test/trivial_runs/test_hicBuildMatrix_trivial_runs.py hits.
struct UserRegion {
    ChromSizes chrom_sizes;
    std::int64_t start = 0;
    std::int64_t end = 0;
};
[[nodiscard]] UserRegion user_region(const ChromSizes& chrom_sizes,
                                     const std::string& region);

// hicexplorer.utilities.genomicRegion: strip whitespace, drop ",;|!{}()" and
// turn every '-' into ':'.
[[nodiscard]] std::string normalise_region(const std::string& text);

// get_bins: even bins of bin_size, one row per chromosome of chrom_sizes.
//
// QUIRK (buildMatrixMethods.py:137-140). With --region the loop start offset
// taken from getUserRegion is applied to *every* chromosome of the returned
// list, not only to the region's own. getUserRegion returns a single
// chromosome, so this is invisible today, but the offset is a loop-invariant
// that the code treats as per chromosome.
[[nodiscard]] std::vector<GenomeInterval> get_bins(std::int64_t bin_size,
                                                   const ChromSizes& chrom_sizes,
                                                   const std::string& region,
                                                   ChromNames& names);

// bed2interval_list over one BED file, appending to `out`.
//
// QUIRK (buildMatrixMethods.py:180). With --region the filter is
//     chrom == region_chrom and region_start <= start and region_end <= end
// and region_end is the *end of the chromosome* whenever --region names a bare
// chromosome. The second condition therefore demands that the cut site reach
// the end of the chromosome, so the restriction site list comes out empty or
// nearly so. That is why the self circle and same fragment classification
// changes as soon as --region is given.
void bed2interval_list(const std::string& path, const ChromSizes& chrom_sizes,
                       const std::string& region, ChromNames& names,
                       std::vector<GenomeInterval>& out);

// get_rf_bins: one bin per restriction site, sites closer than min_distance
// merged, max_distance added on both sides, bins shorter than min_distance
// dropped.
[[nodiscard]] std::vector<GenomeInterval> get_rf_bins(
    const std::vector<GenomeInterval>& cut_sites, std::int64_t min_distance,
    std::int64_t max_distance);

// enlarge_bins: closes the gaps so consecutive bins touch. Modifies in place,
// as the Python does.
void enlarge_bins(std::vector<GenomeInterval>& bins, const ChromSizes& chrom_sizes,
                  const ChromNames& names);

// The flat array of intervals the Python passes to its workers as a
// multiprocessing RawArray of C_Interval, plus the per chromosome index range.
//
// QUIRK (buildMatrixMethods.py:42-46, 962). C_Interval's fields are c_uint, so
// a negative bin start, which get_rf_bins can produce for a cut site closer to
// the chromosome start than --maxLibraryInsertSize, wraps to a value near
// 2^32 in the *lookup* while the bin table written to the matrix keeps the
// negative number. Reproduced by storing uint32 here and int64 there.
struct SearchInterval {
    std::uint32_t begin = 0;
    std::uint32_t end = 0;
    std::uint32_t data = 0;
};

class BinSearchIndex {
  public:
    BinSearchIndex() = default;
    BinSearchIndex(const std::vector<GenomeInterval>& bins, const ChromNames& names);

    // The bin holding `position` on `chrom`, or nullopt. This is
    // buildMatrixMethods.py:658-688 verbatim, including the fact that the
    // comparison is inclusive at both ends and that the midpoint is computed
    // with int(float) truncation rather than a floor division.
    [[nodiscard]] std::optional<std::uint32_t> bin_at(std::uint32_t chrom,
                                                      std::int64_t position) const;
    // The interval of the bin found by the last successful bin_at, which the
    // coverage code needs.
    [[nodiscard]] const SearchInterval& interval(std::int64_t flat_index) const {
        return flat_[static_cast<std::size_t>(flat_index)];
    }
    [[nodiscard]] std::optional<std::int64_t> flat_index_at(
        std::uint32_t chrom, std::int64_t position) const;
    [[nodiscard]] bool has_chrom(std::uint32_t chrom) const;

  private:
    std::vector<SearchInterval> flat_;
    // chrom id -> inclusive [first, last] index range into flat_.
    std::unordered_map<std::uint32_t, std::pair<std::int64_t, std::int64_t>> range_;
};

// The restriction site positions, queried only for "does anything overlap
// [a, b)". intervaltree's tree[a:b] returns intervals with begin < b and
// end > a, and an empty result for a >= b.
class RestrictionSiteIndex {
  public:
    RestrictionSiteIndex() = default;
    RestrictionSiteIndex(const std::vector<GenomeInterval>& sites,
                         std::size_t chrom_count);

    [[nodiscard]] bool has_chrom(std::uint32_t chrom) const;
    [[nodiscard]] bool overlaps(std::uint32_t chrom, std::int64_t begin,
                                std::int64_t end) const;
    [[nodiscard]] bool empty() const noexcept { return total_ == 0; }

  private:
    struct Chrom {
        std::vector<std::int64_t> begin;
        std::vector<std::int64_t> prefix_max_end;  // size begin.size() + 1
    };
    std::vector<Chrom> per_chrom_;
    std::vector<char> present_;
    std::size_t total_ = 0;
};

// Bio.Seq.Seq(s).reverse_complement() over the IUPAC alphabet.
[[nodiscard]] std::string reverse_complement(const std::string& sequence);

// The counters the QC log is built from. Every one of them is an integer sum
// over read pairs, so combining per thread partitions in index order gives the
// same value for any thread count (cpp/OPTIMIZATION.md 3).
struct QcCounters {
    std::int64_t one_mate_unmapped = 0;
    std::int64_t one_mate_low_quality = 0;
    std::int64_t one_mate_not_unique = 0;
    std::int64_t duplicated_pairs = 0;
    std::vector<std::int64_t> dangling_end;  // one per restriction sequence
    std::int64_t self_circle = 0;
    std::int64_t self_ligation = 0;
    std::int64_t same_fragment = 0;
    std::int64_t mate_not_close_to_rf = 0;
    std::int64_t count_inward = 0;
    std::int64_t count_outward = 0;
    std::int64_t count_left = 0;
    std::int64_t count_right = 0;
    std::int64_t inter_chromosomal = 0;
    std::int64_t short_range = 0;
    std::int64_t long_range = 0;
    std::int64_t pair_added = 0;
    std::int64_t iter_num = 0;

    void add(const QcCounters& other);
};

// The parameters of the QC log header, so that the log can be produced without
// dragging the whole argument set in.
struct QcLogInputs {
    std::string out_file_name;
    std::int64_t min_distance = 0;  // 0 means --minDistance was falsy
    std::int64_t max_library_insert_size = 1000;
    bool keep_self_ligation = false;
    bool has_restriction_cut_file = false;
    // Restriction sequence and its dangling sequence, in the order the Python
    // dict holds them.
    std::vector<std::pair<std::string, std::string>> dangling_sequences;
};

// The exact text buildMatrixMethods.createMatrix writes to QCfolder/QC.log and
// stores in the cool metadata under 'statistics'.
[[nodiscard]] std::string format_qc_log(const QcLogInputs& inputs,
                                        const QcCounters& counters);

// The five tab separated tables hicPrepareQCreport writes next to it. The PNGs
// and hicQC.html of that tool are matplotlib and pandas rendering, drawn by
// the tools through hicx::plot::draw("hicPrepareQCreport") after these tables.
void write_qc_tables(const std::string& folder, const std::string& qc_log_text);

// hicPrepareQCreport.main's tables for one or more QC logs: the rows are
// named by --labels when there are as many labels as logs, and by each log's
// File entry otherwise. Throws std::runtime_error with the reference's error
// where pandas raises (logs of different lengths, a missing or text column,
// duplicate row names). write_qc_tables is this function for one log and no
// labels.
void write_qc_report_tables(const std::string& folder, const std::vector<std::string>& log_texts,
                            const std::optional<std::vector<std::string>>& labels);

}  // namespace hicx

#endif  // HICX_BUILD_MATRIX_HPP
