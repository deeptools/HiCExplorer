// The run of hicBuildMatrix, hicBuildMatrixMicroC and hicQuickQC, that is
// hicexplorer/lib/buildMatrixMethods.createMatrix, plus the read pair
// classification of buildMatrixMethods.process_data.
//
// This used to be tools/build_matrix_impl.hpp. It lives in the core because
// the public Python names buildMatrixMethods.createMatrix and process_data
// call it too (cpp/python/API_BINDINGS.md); the tools keep their argument
// parsers, their error printing and the drawing of the QC report. The design
// notes and the reproduced defects Q1 to Q10, D1 and D2 are in
// tools/hicBuildMatrix.cpp.
//
// Arguments::python_api selects what the Python function does where the
// tools deliberately deviate from a crash of the reference (D1, D2 and the
// short --danglingSequence list): the API raises the exception 3.7.6 raised.

#ifndef HICX_BUILD_MATRIX_RUN_HPP
#define HICX_BUILD_MATRIX_RUN_HPP

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "hicx/bam_file.hpp"
#include "hicx/build_matrix.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/matrix_data.hpp"
#include "hicx/python_error.hpp"
#include "hicx/reduce_matrix.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/version.hpp"

namespace hicx::build_matrix {

struct Arguments {
    std::vector<std::string> sam_files;
    std::string out_file_name;
    std::string qc_folder;
    std::vector<std::string> restriction_cut_files;
    std::vector<std::string> restriction_sequences;
    std::vector<std::string> dangling_sequences;
    std::string out_bam;
    std::vector<std::int64_t> bin_size;
    std::int64_t min_distance = 300;
    std::optional<std::int64_t> max_distance;
    std::int64_t max_library_insert_size = 1000;
    std::string genome_assembly;
    std::string region;
    bool keep_self_ligation = false;
    bool keep_self_circles = false;
    std::int64_t min_mapping_quality = 15;
    std::int64_t threads = 4;
    std::int64_t input_buffer_size = 400000;
    bool do_test_run = false;
    std::int64_t do_test_run_lines = 1000000;
    bool skip_duplication_check = false;
    std::string chromosome_sizes;
    // hicBuildMatrix --pairsFile only (tools/build_matrix_pairs_impl.hpp): the
    // .pairs input, and whether --minMappingQuality was given explicitly,
    // since without mapq columns an explicit value is refused and the default
    // is not applied.
    std::string pairs_file;
    bool min_mapping_quality_given = false;
    // hicBuildMatrix --noPlot: QC.log and the tables without the figures.
    bool no_plot = false;

    // --- buildMatrixMethods.createMatrix called from Python ----------------
    bool python_api = false;
    // pDanglingSequence was given (truthy) and pRestrictionSequence is not
    // None: the dangling sequences dictionary is built and logged.
    bool dangling_given = false;
    // pRestrictionCutFile is not None, which decides the QC log rows.
    bool restriction_cut_file_given = false;
    // pMinDistance is None (min_distance is then 0, which is falsy as None is).
    bool min_distance_none = false;
};

// One restriction site of a BED file, as bed2interval_list returns it.
using RfInterval = std::tuple<std::string, std::int64_t, std::int64_t>;

// What the caller of run provides besides the arguments. Every
// member may be left empty, which gives the tools' behaviour.
struct Hooks {
    // A log record of hicexplorer.lib.buildMatrixMethods ("debug", "info",
    // "error"). Empty: the tools print the one record they always printed,
    // "reading ... to build hic_matrix", to stderr.
    std::function<void(const char* level, const std::string& message)> log;
    // warnings.warn(message) (the tools print the message to stderr).
    std::function<void(const std::string& message)> warn;
    // Called once both BAM files are open, where createMatrix closes the file
    // objects of pSamFiles and pOutBam.
    std::function<void()> bams_opened;
    // createMatrix upper-cases pRestrictionSequence[index] and, with
    // `dangling_too`, pDanglingSequence[index] in place.
    std::function<void(std::size_t index, bool dangling_too)> uppercase;
    // The chromosome sizes of pChromosomeSizes, read by the caller (the
    // tools read Arguments::chromosome_sizes).
    std::function<ChromSizes()> chromosome_sizes;
    // The restriction sites of every pRestrictionCutFile, read by the caller
    // with bed2interval_list (the tools read Arguments::restriction_cut_files).
    std::function<std::vector<RfInterval>(const ChromSizes&)> restriction_sites;
    // getUserRegion(chrom_sizes, pRegion) for get_bins: the chromosome list
    // and the start of the region (the tools resolve Arguments::region).
    std::function<std::pair<ChromSizes, std::int64_t>(const ChromSizes&)> user_region;
    // Called after QC.log is written, where createMatrix calls
    // hicPrepareQCreport.main; the tools write the five tables instead.
    std::function<void(const std::string& qc_log_file)> after_qc_log;
    // Called before the matrix file is written, where createMatrix closes and
    // removes pOutFileName.
    std::function<void()> before_matrix;
};

// An error the tools report with exit status 1 and createMatrix with
// SystemExit: with `logged` the message is a log.error record followed by
// exit(1), otherwise exit(message).
class Exit : public std::runtime_error {
  public:
    Exit(const std::string& message, bool logged) : std::runtime_error(message), logged_(logged) {}
    [[nodiscard]] bool logged() const noexcept { return logged_; }

  private:
    bool logged_;
};

// Runs the whole of createMatrix: reads the BAM pair, classifies, writes
// QC.log and the matrix. Throws Exit, PythonError (python_api only) and
// std::runtime_error.
void run(const Arguments& args, const Hooks& hooks = {});

[[nodiscard]] inline bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

[[nodiscard]] std::string to_upper(std::string text);

// The chromosome sizes file of --chromosomeSizes, as the tools read it.
[[nodiscard]] ChromSizes read_chromosome_sizes(const std::string& path);

// The text of a Python str repr for ASCII text ('...' unless the text holds
// a ' and no "), with \t, \n, \r, \\ and the other control characters escaped
// as Python does; bytes above 0x7f are copied as they are.
[[nodiscard]] std::string python_str_repr(const std::string& text);

// --- the read pair classification -------------------------------------------

// readBamFiles' filters of a primary pair, buildMatrixMethods.py:513-540:
// either mate unmapped (flag 0x4), then a mapping quality below the minimum,
// which counts as "not unique" when the first mate's quality is 0 (Q1: the
// Python's `mate1.mapq == 0 & mate2.mapq == 0` is the chained comparison
// `mate1.mapq == (0 & mate2.mapq) == 0`) and as "low quality" otherwise.
enum class PairVerdict { OneMateUnmapped, OneMateNotUnique, LowQuality, Keep };
[[nodiscard]] inline PairVerdict pair_verdict(std::int64_t flag1, std::int64_t flag2, std::int64_t mapq1,
                                              std::int64_t mapq2, std::int64_t min_mapping_quality) {
    if ((flag1 & 0x4) == 4 || (flag2 & 0x4) == 4) {
        return PairVerdict::OneMateUnmapped;
    }
    if (mapq1 < min_mapping_quality || mapq2 < min_mapping_quality) {
        return mapq1 == 0 ? PairVerdict::OneMateNotUnique : PairVerdict::LowQuality;
    }
    return PairVerdict::Keep;
}

// The parts of process_data's arguments that every pair is classified with.
struct ClassifyOptions {
    // pRestrictionSequence (empty for None)
    std::vector<std::string> restriction_sequences;
    // `if pDanglingSequences`: the dangling sequence dictionary is not empty.
    bool has_dangling = false;
    bool keep_self_circles = false;
    bool keep_self_ligation = false;
    std::int64_t max_insert_size = 1000;
    bool quick_qc_mode = false;
    // pBinsize, the coverage cell width (createMatrix passes 10).
    std::int64_t coverage_bin_size = 10;
};

// process_data's per pair loop, buildMatrixMethods.py:642-861, over the pairs
// [first, last) of a buffer, adding to `counters`. Q2 and Q3 of
// tools/hicBuildMatrix.cpp are reproduced. The policy P reads the pairs and
// receives the results, so that the tools (htslib records) and the Python
// function (pysam objects) run the same decisions:
//
//   std::int64_t tid(i, w), pos(i, w), qlen(i, w), seq_len(i, w)
//   bool is_reverse(i, w)                 mate w (0 or 1) of pair i
//   std::optional<SearchInterval> bin(i, w, read_middle)
//                                         the bin lookup of :652-688, nullopt
//                                         for "not close to a restriction site"
//   bool rf_enabled()                     `pRfPositions and pRestrictionSequence`
//   bool rf_overlaps(i, begin, end)       a restriction site of mate 1's
//                                         chromosome overlaps [begin, end)
//   bool dangling(i, w, s)                check_dangling_end of mate w against
//                                         restriction sequence s
//   std::pair<std::int64_t, std::int64_t> coverage_range(bin_id)
//                                         pCoverageIndex[bin_id].begin, .end
//   void coverage_add(std::int64_t cell)  pCoverage[cell] += 1
//   void set_isize(i)                     :814-815
//   void add_pair(i, bin_of_mate1, bin_of_mate2)   a pair counted
//
// Python's evaluation order is kept, so a policy reading attributes lazily
// fails where process_data fails.
template <class P>
void classify_pairs(const ClassifyOptions& options, std::size_t first, std::size_t last, P& p,
                    QcCounters& counters) {
    if (counters.dangling_end.size() < options.restriction_sequences.size()) {
        counters.dangling_end.resize(options.restriction_sequences.size(), 0);
    }
    const bool has_restriction = !options.restriction_sequences.empty();
    for (std::size_t i = first; i < last; ++i) {
        // Bin lookup, both mates. The interval of the second mate's bin
        // survives the loop and is used for both mates' coverage (Q3).
        std::uint32_t mate_bins[2] = {0, 0};
        SearchInterval mate_bin;
        bool unassigned = false;
        for (int w = 0; w < 2; ++w) {
            const std::int64_t read_middle = p.pos(i, w) + p.qlen(i, w) / 2;
            const std::optional<SearchInterval> found = p.bin(i, w, read_middle);
            if (!found.has_value()) {
                unassigned = true;
                break;
            }
            mate_bin = *found;
            mate_bins[w] = found->data;
        }
        if (unassigned) {
            ++counters.mate_not_close_to_rf;
            continue;
        }

        enum class Orientation { DiffChromosome, Inward, Outward, Left, Right };
        Orientation orientation = Orientation::DiffChromosome;
        if (p.tid(i, 0) != p.tid(i, 1)) {
            orientation = Orientation::DiffChromosome;
        } else {
            const int first_mate = p.pos(i, 0) < p.pos(i, 1) ? 0 : 1;
            const int second_mate = 1 - first_mate;
            const bool first_reverse = p.is_reverse(i, first_mate);
            const bool second_reverse = p.is_reverse(i, second_mate);
            if (!first_reverse && second_reverse) {
                orientation = Orientation::Inward;
            } else if (first_reverse && !second_reverse) {
                orientation = Orientation::Outward;
            } else if (first_reverse && second_reverse) {
                orientation = Orientation::Left;
            } else {
                orientation = Orientation::Right;
            }

            auto separation = [&] {
                const std::int64_t d = p.pos(i, 1) - p.pos(i, 0);
                return d < 0 ? -d : d;
            };
            auto fragment = [&](const std::string& sequence) {
                const auto length = static_cast<std::int64_t>(sequence.size());
                const std::int64_t frag_start = std::min(p.pos(i, 0), p.pos(i, 1)) + length;
                const std::int64_t frag_end =
                    std::max(p.pos(i, 0) + p.qlen(i, 0), p.pos(i, 1) + p.qlen(i, 1)) - length;
                return std::make_pair(frag_start, frag_end);
            };

            if (separation() < 25000 && orientation == Orientation::Outward) {
                bool has_rf = false;
                if (p.rf_enabled()) {
                    for (const std::string& sequence : options.restriction_sequences) {
                        const auto [frag_start, frag_end] = fragment(sequence);
                        if (p.rf_overlaps(i, frag_start, frag_end)) {
                            has_rf = true;
                        }
                        if (!has_rf) {
                            ++counters.self_circle;
                            if (!options.keep_self_circles) {
                                continue;  // Q2: the read pair is NOT skipped
                            }
                        }
                    }
                }
            }
            if (separation() < options.max_insert_size && orientation == Orientation::Inward) {
                if (has_restriction && options.has_dangling) {
                    bool one_match = false;
                    for (std::size_t s = 0; s < options.restriction_sequences.size(); ++s) {
                        if (p.dangling(i, 0, s) || p.dangling(i, 1, s)) {
                            ++counters.dangling_end[s];
                            one_match = true;
                            break;
                        }
                    }
                    if (one_match) {
                        continue;
                    }
                }
                bool has_rf = false;
                if (p.rf_enabled()) {
                    for (const std::string& sequence : options.restriction_sequences) {
                        const auto [frag_start, frag_end] = fragment(sequence);
                        if (p.rf_overlaps(i, frag_start, frag_end)) {
                            has_rf = true;
                        }
                    }
                }
                if (!has_rf) {
                    ++counters.same_fragment;
                    continue;
                }
                ++counters.self_ligation;
                if (!options.keep_self_ligation) {
                    continue;
                }
            }
            p.set_isize(i);
        }

        if (p.tid(i, 0) != p.tid(i, 1)) {
            ++counters.inter_chromosomal;
        } else {
            const std::int64_t d = p.pos(i, 1) - p.pos(i, 0);
            if ((d < 0 ? -d : d) < 20000) {
                ++counters.short_range;
            } else {
                ++counters.long_range;
            }
        }
        switch (orientation) {
            case Orientation::Inward: ++counters.count_inward; break;
            case Orientation::Outward: ++counters.count_outward; break;
            case Orientation::Left: ++counters.count_left; break;
            case Orientation::Right: ++counters.count_right; break;
            case Orientation::DiffChromosome: break;
        }

        // Coverage, with Q3's bin for the first mate.
        const std::uint32_t bin_id = mate_bin.data;
        for (int w = 0; w < 2; ++w) {
            const std::int64_t offset = std::max<std::int64_t>(0, p.pos(i, w) - static_cast<std::int64_t>(mate_bin.begin));
            const std::int64_t bin_size = options.coverage_bin_size;
            if (bin_size == 0) {
                throw PythonError("ZeroDivisionError", "division by zero");
            }
            // int(x / pBinsize) of the Python: true division, then truncation
            const auto truncated_quotient = [bin_size](std::int64_t x) {
                return bin_size > 0 && x >= 0
                           ? x / bin_size
                           : static_cast<std::int64_t>(static_cast<double>(x) / static_cast<double>(bin_size));
            };
            const std::int64_t vec_start = truncated_quotient(offset);
            const auto [cov_begin, cov_end] = p.coverage_range(bin_id);
            const std::int64_t length_coverage = cov_end - cov_begin;
            std::int64_t vec_end = 0;
            if (bin_size > 0) {
                vec_end = std::min(length_coverage, vec_start + p.seq_len(i, w) / bin_size);
            } else {
                vec_end = std::min(length_coverage,
                                   static_cast<std::int64_t>(static_cast<double>(vec_start) +
                                                             static_cast<double>(p.seq_len(i, w)) /
                                                                 static_cast<double>(bin_size)));
            }
            for (std::int64_t cell = cov_begin + vec_start; cell < cov_begin + vec_end; ++cell) {
                p.coverage_add(cell);
            }
        }

        p.add_pair(i, mate_bins[0], mate_bins[1]);
        ++counters.pair_added;
    }
}

// --- the tools' pieces, shared with the pairs route ---------------------------

// The counts of one chunk are sorted and coalesced, then merged into a running
// sorted (key, count) accumulator. Memory is O(distinct pixels) rather than
// O(read pairs), and the merge is order independent because it merges sorted
// keys and adds integers.
class PixelAccumulator {
  public:
    void absorb(std::vector<std::uint64_t>& chunk);
    [[nodiscard]] const std::vector<std::uint64_t>& keys() const { return keys_; }
    [[nodiscard]] const std::vector<std::int64_t>& counts() const { return counts_; }

  private:
    void merge(const std::vector<std::uint64_t>& keys, const std::vector<std::int64_t>& counts);

    std::vector<std::uint64_t> keys_;
    std::vector<std::int64_t> counts_;
};

// The provenance metadata of the matrix file (Q9 and Q10 of
// tools/hicBuildMatrix.cpp).
[[nodiscard]] std::map<std::string, std::string> matrix_metadata(const std::string& qc_log,
                                                                 const std::string& genome_assembly);

// Writes `data` the way createMatrix saves it: an mcool with one group per
// resolution for three or more bin sizes (Q7, Q10), else h5 or cool by the
// name's end, else nothing (hiCMatrix.save's behaviour).
void save_matrix(const std::string& out_file_name, const std::vector<std::int64_t>& bin_size,
                 MatrixData& data, const std::map<std::string, std::string>& metadata);

// The matrix file of a finished run: the upper triangle of the pixel counts
// over the enlarged bins, written as h5, cool or mcool the way
// buildMatrixMethods.createMatrix writes them. Shared by the BAM route and the
// .pairs route (tools/build_matrix_pairs_impl.hpp). There are `pixel_count`
// distinct pixels, which for_each_pixel(sink) hands to sink(key, count) in
// ascending key order, key = min_bin << 32 | max_bin; a source may release its
// storage while it does so. `bin_max` is the coverage column of the h5 bin
// table; `genome_assembly` may be empty.
template <class ForEachPixel>
void write_matrix_outputs(const Arguments& args, const std::vector<GenomeInterval>& bin_intervals,
                          const std::vector<double>& bin_max, const ChromNames& names,
                          std::size_t pixel_count, ForEachPixel&& for_each_pixel,
                          const std::string& qc_log, const std::string& genome_assembly) {
    const std::int64_t matrix_size = static_cast<std::int64_t>(bin_intervals.size());
    MatrixData data;
    data.cut_intervals.reserve(bin_intervals.size());
    for (std::size_t i = 0; i < bin_intervals.size(); ++i) {
        CutInterval interval;
        interval.chrom = names.name(bin_intervals[i].chrom);
        interval.start = bin_intervals[i].start;
        interval.end = bin_intervals[i].end;
        interval.extra = bin_max[i];
        data.cut_intervals.push_back(std::move(interval));
    }

    // The upper triangle of C + C.T - diag(C), built directly: a pair
    // (a, b) is one count at (min, max), and the diagonal is counted once
    // because the Python subtracts the doubled diagonal back off.
    {
        std::vector<std::int64_t> indptr(static_cast<std::size_t>(matrix_size) + 1, 0);
        // Reserved rather than sized, so that pages are touched only as the
        // pixels arrive, while a draining source frees its own blocks.
        std::vector<std::int32_t> indices;
        std::vector<double> values;
        indices.reserve(pixel_count);
        values.reserve(pixel_count);
        for_each_pixel([&](std::uint64_t key, std::uint64_t count) {
            const auto row = static_cast<std::int64_t>(key >> 32);
            indices.push_back(static_cast<std::int32_t>(key & 0xffffffffULL));
            values.push_back(static_cast<double>(count));
            ++indptr[static_cast<std::size_t>(row) + 1];
        });
        for (std::size_t r = 0; r < static_cast<std::size_t>(matrix_size); ++r) {
            indptr[r + 1] += indptr[r];
        }
        data.matrix = CsrMatrix(matrix_size, matrix_size, std::move(indptr), std::move(indices),
                                std::move(values), "int64");
        data.matrix.set_symmetry(Symmetry::UpperTriangle);
    }

    save_matrix(args.out_file_name, args.bin_size, data, matrix_metadata(qc_log, genome_assembly));
}

}  // namespace hicx::build_matrix

#endif  // HICX_BUILD_MATRIX_RUN_HPP
