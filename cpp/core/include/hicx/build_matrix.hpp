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

#include <array>
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
// The loop of get_bins after the region is resolved: range(start, size,
// bin_size) for every chromosome of `chrom_sizes`, with the one start offset.
[[nodiscard]] std::vector<GenomeInterval> get_bins(std::int64_t bin_size,
                                                   const ChromSizes& chrom_sizes,
                                                   std::int64_t start, ChromNames& names);

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
//
// `start_is_python_zero`, when given, receives for every returned bin whether
// its start is the Python int 0 of max(0, start - max_distance) rather than a
// numpy.int64: the Python computes the coordinates on numpy arrays, and the
// public get_rf_bins returns those types.
[[nodiscard]] std::vector<GenomeInterval> get_rf_bins(
    const std::vector<GenomeInterval>& cut_sites, std::int64_t min_distance,
    std::int64_t max_distance, std::vector<bool>* start_is_python_zero = nullptr);

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

// The binary search of buildMatrixMethods.py:658-681 over the flat
// intervals [start, end] (inclusive indices) of one chromosome: the index of
// the interval with begin <= position <= end, or nullopt. `probe(index)`
// returns the interval at an index; the midpoint is int((start + end) / 2).
template <class Probe>
[[nodiscard]] std::optional<std::int64_t> search_bins(std::int64_t start, std::int64_t end,
                                                      std::uint64_t position, Probe&& probe) {
    std::int64_t middle = (start + end) / 2;
    while (!(start > end)) {
        const SearchInterval interval = probe(middle);
        if (static_cast<std::uint64_t>(interval.begin) <= position &&
            position <= static_cast<std::uint64_t>(interval.end)) {
            return middle;
        }
        if (static_cast<std::uint64_t>(interval.begin) > position) {
            end = middle - 1;
        } else {
            start = middle + 1;
        }
        middle = (start + end) / 2;
    }
    return std::nullopt;
}

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

    // Intra-chromosomal contacts per log-spaced distance bin, taken from the
    // read pairs before binning (the distance of the two read start
    // positions, which the short/long range counts of the QC log use as
    // well). qc_distance_bin() defines the bins. Empty until a pair is
    // recorded. `orientation_distance[o]` is the same histogram for the
    // read pair type o (0 inward, 1 outward, 2 left, 3 right).
    std::vector<std::int64_t> distance_contacts;
    std::array<std::vector<std::int64_t>, 4> orientation_distance;

    // Records one intra-chromosomal contact; `orientation` is -1 when the
    // read pair type is unknown (a .pairs file without strands).
    void record_intra(std::int64_t distance, int orientation);

    void add(const QcCounters& other);
};

// The distance bins of the QC tables stored in the cool metadata: 10 bins
// per decade. Bin i holds the integer distances d with
// ceil(10^(i/10)) <= d < ceil(10^((i+1)/10)); bin 0 starts at 1 bp and also
// takes the distance 0. kQcDistanceBins bins reach 10^10 bp.
inline constexpr int kQcDistanceBinsPerDecade = 10;
inline constexpr int kQcDistanceBins = 100;
// The lower edge (bp) of bin i, 0 <= i <= kQcDistanceBins.
[[nodiscard]] std::int64_t qc_distance_edge(int i);
[[nodiscard]] int qc_distance_bin(std::int64_t distance);

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
    // The longest chromosome of the run; the distance tables of the metadata
    // end at the bin holding it. 0 means: at the last non-empty bin.
    std::int64_t max_chromosome_length = 0;
};

// The exact text buildMatrixMethods.createMatrix writes to QCfolder/QC.log and
// stores in the cool metadata under 'statistics'.
[[nodiscard]] std::string format_qc_log(const QcLogInputs& inputs,
                                        const QcCounters& counters);

// The QC statistics as the JSON object stored under the key "hicexplorer_qc"
// of the cool `metadata` attribute: every count of the QC log, the
// parameters of the run, the tool name and version, and the two distance
// tables (contacts per log-spaced distance bin, and the same split by read
// pair type). Schema version 1; see docs/tools/hicMatrixQC.md.
[[nodiscard]] std::string format_qc_metadata_json(const QcLogInputs& inputs,
                                                  const QcCounters& counters,
                                                  const std::string& tool);

// The five tab separated tables hicPrepareQCreport writes next to it. The PNGs
// of that tool are matplotlib and pandas rendering, drawn by the tools through
// hicx::plot::draw("hicPrepareQCreport") after these tables and
// write_qc_report_html.
void write_qc_tables(const std::string& folder, const std::string& qc_log_text);

// hicPrepareQCreport.main's tables for one or more QC logs: the rows are
// named by --labels when there are as many labels as logs, and by each log's
// File entry otherwise. Throws std::runtime_error with the reference's error
// where pandas raises (logs of different lengths, a missing or text column,
// duplicate row names). write_qc_tables is this function for one log and no
// labels.
void write_qc_report_tables(const std::string& folder, const std::vector<std::string>& log_texts,
                            const std::optional<std::vector<std::string>>& labels);

// save_html of hicPrepareQCreport.main for the same logs and labels:
// folder/hicQC.html, the four derived tables and the QC table rendered by
// pandas 2.2's Styler.to_html into hicexplorer/qc_template.html (compiled
// in). Byte identical to the reference except for the five random hexadecimal
// characters of each table id (uuid4().hex[:5] there, random here too).
// Returns the DataFrames the five charts plot (make_figure_*) as the JSON
// object {"pairs_sequenced", "unmappable_and_non_unique", "pairs_discarded",
// "distance", "read_orientation"}, each {"index", "index_name", "columns",
// "columns_name", "values" (row major)} with every label as str() prints it,
// for hicexplorer_plot. Throws as write_qc_report_tables does.
[[nodiscard]] std::string write_qc_report_html(const std::string& folder,
                                               const std::vector<std::string>& log_texts,
                                               const std::optional<std::vector<std::string>>& labels);

// The pass/warn/fail cutoffs hicx::write_qc_pretty_report judges each sample
// against. Two defaults (min_valid_pairs_percent, max_inter_chromosomal_percent)
// come from HiCExplorer's own usage guidance (docs/content/example_usage.rst:
// "usually 25%-40% of the reads are valid"; "good Hi-C libraries have lower
// than 10% inter chromosomal contacts"). The rest have no published HiCExplorer
// guidance behind them; they are conservative, commonly used starting points,
// stated as such in the report itself, and every one is a CLI option so a lab
// can set its own.
struct QcThresholds {
    double min_valid_pairs_percent = 25.0;
    double max_inter_chromosomal_percent = 10.0;
    double max_duplicated_pairs_percent = 20.0;
    double max_dangling_end_percent = 10.0;
    double max_self_ligation_percent = 10.0;
    // A value on the right side of its cutoff, or exactly on it, is PASS. A
    // value on the wrong side by at most this fraction of the cutoff is WARN,
    // and one further out is FAIL. With 0.2, inter-chromosomal contacts are
    // PASS up to 10 %, WARN above 10 % up to 12 % and FAIL above 12 % against
    // the 10 % cutoff; valid pairs are PASS from 25 %, WARN from 20 % up to
    // 25 % and FAIL below 20 % against the 25 % cutoff. The margin is not a
    // command line option.
    double warn_margin_fraction = 0.2;
};

// An HTML report at folder/hicQC_report.html in the design of the
// falcoAmadeus report (light/dark theme, palette picker, a summary sidebar
// with the worst verdict per section, collapsible section cards). The
// overview has one card per sample (row of the --logfiles/--labels input)
// with five metrics, each judged PASS/WARN/FAIL against `thresholds` (see
// QcThresholds::warn_margin_fraction for the WARN band) and shown as a
// coloured badge; the other sections hold write_qc_report_tables's numbers
// as plain HTML tables plus grouped bar charts drawn as HTML tables too (a bar
// is a table cell with a percentage width and a bgcolor). The report needs no
// JavaScript and no external file: Galaxy's HTML sanitizer, which removes
// <script>, <style>, <svg> and inline styles, leaves every chart, table and
// anchor link readable. An inline script only adds the theme toggle, palette
// picker, collapsible cards and sortable tables. C++-only: no Python
// HiCExplorer counterpart, and does not replace hicQC.html, which is
// unchanged and still written alongside it when the caller also draws it.
void write_qc_pretty_report(const std::string& folder, const std::vector<std::string>& log_texts,
                            const std::optional<std::vector<std::string>>& labels,
                            const QcThresholds& thresholds);

}  // namespace hicx

#endif  // HICX_BUILD_MATRIX_HPP
