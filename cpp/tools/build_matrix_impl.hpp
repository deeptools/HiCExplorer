// Shared implementation of hicBuildMatrix and hicBuildMatrixMicroC.
//
// hicexplorer/hicBuildMatrixMicroC.py is hicBuildMatrix.py with the
// restriction enzyme arguments removed and --binSize made required: its main()
// calls the same lib/buildMatrixMethods.createMatrix with pRestrictionCutFile,
// pRestrictionSequence, pDanglingSequence, pMinDistance, pMaxDistance and
// pKeepSelfLigation all None. So the two tools differ in their argument parser
// and in nothing else, and this header is the "nothing else".
//
// The design notes and the list of reproduced defects are in
// hicBuildMatrix.cpp.

#ifndef HICX_TOOLS_BUILD_MATRIX_IMPL_HPP
#define HICX_TOOLS_BUILD_MATRIX_IMPL_HPP

#include <atomic>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#include "hicx/bam_file.hpp"
#include "hicx/build_matrix.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/reduce_matrix.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/version.hpp"

namespace {

// Set by main before anything else, so that the shared code can name the tool
// the user invoked in its usage and error messages.
const char* g_tool = "hicBuildMatrix";

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
    // hicBuildMatrix --pairsFile only (build_matrix_pairs_impl.hpp): the
    // .pairs input, and whether --minMappingQuality was given explicitly,
    // since without mapq columns an explicit value is refused and the default
    // is not applied.
    std::string pairs_file;
    bool min_mapping_quality_given = false;
    // hicBuildMatrix --noPlot: QC.log and the tables without the figures.
    bool no_plot = false;
};

// hicexplorer.utilities.genomicRegion as an argparse type: a value that is
// only whitespace becomes None, one that is only separators is an error.
// hicQuickQC includes this file and has no --region.
[[maybe_unused]] std::optional<std::string> genomic_region_check(const std::string& text) {
    bool blank = true;
    for (const char c : text) {
        if (std::isspace(static_cast<unsigned char>(c)) == 0) {
            blank = false;
            break;
        }
    }
    if (!blank && hicx::normalise_region(text).empty()) {
        return text + " is not a valid region";
    }
    return std::nullopt;
}

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string to_upper(std::string text) {
    for (char& c : text) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return text;
}

// --- the read pair buffer ---------------------------------------------------

struct PairBuffer {
    std::vector<hicx::ReadFields> mate1;
    std::vector<hicx::ReadFields> mate2;
    // Kept only for --outBam.
    std::vector<hicx::BamRecord> record1;
    std::vector<hicx::BamRecord> record2;

    void clear() {
        mate1.clear();
        mate2.clear();
        record1.clear();
        record2.clear();
    }
    [[nodiscard]] std::size_t size() const { return mate1.size(); }
};

// ReadPositionMatrix. The Python key is a string built from the two chromosome
// names ordered by string comparison and the two positions ordered
// numerically. The tuple below is equivalent for every reference name that
// holds no '-', which is every name in the corpus, and costs 16 bytes against
// the string's 80.
struct DuplicateKey {
    std::uint32_t chrom_low = 0;
    std::uint32_t chrom_high = 0;
    std::int32_t pos_low = 0;
    std::int32_t pos_high = 0;

    bool operator==(const DuplicateKey& other) const = default;
};

struct DuplicateKeyHash {
    std::size_t operator()(const DuplicateKey& key) const noexcept {
        std::uint64_t a = (static_cast<std::uint64_t>(key.chrom_low) << 32) |
                          static_cast<std::uint32_t>(key.chrom_high);
        std::uint64_t b = (static_cast<std::uint64_t>(
                               static_cast<std::uint32_t>(key.pos_low))
                           << 32) |
                          static_cast<std::uint32_t>(key.pos_high);
        a ^= b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2);
        return static_cast<std::size_t>(a);
    }
};

// --- classification ---------------------------------------------------------

struct Context {
    const hicx::BinSearchIndex* bins = nullptr;
    const hicx::RestrictionSiteIndex* rf = nullptr;
    bool rf_present = false;  // `if pRfPositions` truthiness
    std::vector<std::uint32_t> tid_to_chrom;
    std::vector<std::string> restriction_sequences;
    std::vector<std::string> pat_forw;
    std::vector<std::string> pat_rev;
    bool has_dangling = false;
    bool keep_self_circles = false;
    bool keep_self_ligation = false;
    std::int64_t max_insert_size = 1000;
    bool quick_qc_mode = false;
    std::vector<std::uint32_t> coverage_begin;
    std::vector<std::uint32_t> coverage_end;
    std::atomic<std::uint32_t>* coverage = nullptr;
    std::size_t coverage_size = 0;
};

bool starts_with_probe(const hicx::ReadFields& read, const std::string& pattern) {
    if (pattern.empty()) {
        return true;
    }
    if (static_cast<int>(pattern.size()) > read.seq_len ||
        static_cast<int>(pattern.size()) > hicx::kDanglingProbe) {
        return false;
    }
    return std::memcmp(read.head, pattern.data(), pattern.size()) == 0;
}

bool ends_with_probe(const hicx::ReadFields& read, const std::string& pattern) {
    if (pattern.empty()) {
        return true;
    }
    const int n = static_cast<int>(pattern.size());
    if (n > read.seq_len || n > hicx::kDanglingProbe) {
        return false;
    }
    const int stored = std::min(read.seq_len, hicx::kDanglingProbe);
    return std::memcmp(read.tail + (stored - n), pattern.data(), pattern.size()) == 0;
}

bool check_dangling_end(const hicx::ReadFields& read, const std::string& forward,
                        const std::string& reverse) {
    if (!read.is_reverse() && starts_with_probe(read, forward)) {
        return true;
    }
    if (read.is_reverse() && ends_with_probe(read, reverse)) {
        return true;
    }
    return false;
}

struct Partition {
    hicx::QcCounters counters;
    std::vector<std::uint64_t> pixels;      // packed (min_bin << 32) | max_bin
    std::vector<std::int64_t> bam_indices;  // indices into the pair buffer
};

void classify_range(const Context& context, const PairBuffer& buffer,
                    std::size_t first, std::size_t last, bool want_bam,
                    Partition& out) {
    out.counters.dangling_end.assign(context.restriction_sequences.size(), 0);
    for (std::size_t index = first; index < last; ++index) {
        const hicx::ReadFields& mate1 = buffer.mate1[index];
        const hicx::ReadFields& mate2 = buffer.mate2[index];

        // Bin lookup, both mates. mate_bin and mate_bin_id survive the loop
        // and are reused by the coverage code below; see Q3.
        std::uint32_t mate_bins[2] = {0, 0};
        int found = 0;
        std::int64_t last_flat = -1;
        bool unassigned = false;
        for (int which = 0; which < 2; ++which) {
            const hicx::ReadFields& mate = which == 0 ? mate1 : mate2;
            const std::uint32_t chrom =
                context.tid_to_chrom[static_cast<std::size_t>(mate.tid)];
            const std::int64_t read_middle = mate.pos + mate.qlen / 2;
            const auto flat = context.bins->flat_index_at(chrom, read_middle);
            if (!flat.has_value()) {
                unassigned = true;
                break;
            }
            last_flat = *flat;
            mate_bins[found++] = context.bins->interval(*flat).data;
        }
        if (unassigned) {
            ++out.counters.mate_not_close_to_rf;
            continue;
        }

        enum class Orientation { DiffChromosome, Inward, Outward, Left, Right };
        Orientation orientation = Orientation::DiffChromosome;
        bool skip_pair = false;
        if (mate1.tid != mate2.tid) {
            orientation = Orientation::DiffChromosome;
        } else {
            const hicx::ReadFields& first_mate = mate1.pos < mate2.pos ? mate1 : mate2;
            const hicx::ReadFields& second_mate = mate1.pos < mate2.pos ? mate2 : mate1;
            if (!first_mate.is_reverse() && second_mate.is_reverse()) {
                orientation = Orientation::Inward;
            } else if (first_mate.is_reverse() && !second_mate.is_reverse()) {
                orientation = Orientation::Outward;
            } else if (first_mate.is_reverse() && second_mate.is_reverse()) {
                orientation = Orientation::Left;
            } else {
                orientation = Orientation::Right;
            }

            const std::int64_t separation =
                std::llabs(static_cast<std::int64_t>(mate2.pos) - mate1.pos);
            const std::uint32_t chrom =
                context.tid_to_chrom[static_cast<std::size_t>(mate1.tid)];

            if (separation < 25000 && orientation == Orientation::Outward) {
                bool has_rf = false;
                if (context.rf_present && !context.restriction_sequences.empty()) {
                    for (const auto& sequence : context.restriction_sequences) {
                        const std::int64_t frag_start =
                            std::min(mate1.pos, mate2.pos) +
                            static_cast<std::int64_t>(sequence.size());
                        const std::int64_t frag_end =
                            std::max(mate1.pos + mate1.qlen, mate2.pos + mate2.qlen) -
                            static_cast<std::int64_t>(sequence.size());
                        if (context.rf->has_chrom(chrom)) {
                            has_rf = has_rf ||
                                     context.rf->overlaps(chrom, frag_start, frag_end);
                        }
                        if (!has_rf) {
                            ++out.counters.self_circle;
                            if (!context.keep_self_circles) {
                                continue;  // Q2: the read pair is NOT skipped
                            }
                        }
                    }
                }
            }
            if (separation < context.max_insert_size &&
                orientation == Orientation::Inward) {
                if (!context.restriction_sequences.empty() && context.has_dangling) {
                    bool one_match = false;
                    for (std::size_t s = 0; s < context.restriction_sequences.size();
                         ++s) {
                        if (check_dangling_end(mate1, context.pat_forw[s],
                                               context.pat_rev[s]) ||
                            check_dangling_end(mate2, context.pat_forw[s],
                                               context.pat_rev[s])) {
                            ++out.counters.dangling_end[s];
                            one_match = true;
                            break;
                        }
                    }
                    if (one_match) {
                        skip_pair = true;
                    }
                }
                if (!skip_pair) {
                    bool has_rf = false;
                    if (context.rf_present && !context.restriction_sequences.empty()) {
                        for (const auto& sequence : context.restriction_sequences) {
                            const std::int64_t frag_start =
                                std::min(mate1.pos, mate2.pos) +
                                static_cast<std::int64_t>(sequence.size());
                            const std::int64_t frag_end =
                                std::max(mate1.pos + mate1.qlen,
                                         mate2.pos + mate2.qlen) -
                                static_cast<std::int64_t>(sequence.size());
                            if (context.rf->has_chrom(chrom)) {
                                has_rf =
                                    has_rf ||
                                    context.rf->overlaps(chrom, frag_start, frag_end);
                            }
                        }
                    }
                    if (!has_rf) {
                        ++out.counters.same_fragment;
                        skip_pair = true;
                    } else {
                        ++out.counters.self_ligation;
                        if (!context.keep_self_ligation) {
                            skip_pair = true;
                        }
                    }
                }
            }
            // Q5: the isize patch of buildMatrixMethods.py:814 happens in the
            // worker's copy and never reaches the written record, so it is
            // deliberately not applied here.
        }
        if (skip_pair) {
            continue;
        }
        if (found != 2) {
            continue;
        }

        if (mate1.tid != mate2.tid) {
            ++out.counters.inter_chromosomal;
        } else if (std::llabs(static_cast<std::int64_t>(mate2.pos) - mate1.pos) <
                   20000) {
            ++out.counters.short_range;
        } else {
            ++out.counters.long_range;
        }
        switch (orientation) {
            case Orientation::Inward: ++out.counters.count_inward; break;
            case Orientation::Outward: ++out.counters.count_outward; break;
            case Orientation::Left: ++out.counters.count_left; break;
            case Orientation::Right: ++out.counters.count_right; break;
            case Orientation::DiffChromosome: break;
        }

        // Coverage, with Q3's wrong bin for the first mate.
        const hicx::SearchInterval& bin = context.bins->interval(last_flat);
        const std::uint32_t bin_id = bin.data;
        const std::uint32_t cov_begin = context.coverage_begin[bin_id];
        const std::uint32_t cov_end = context.coverage_end[bin_id];
        const std::int64_t length_coverage =
            static_cast<std::int64_t>(cov_end) - static_cast<std::int64_t>(cov_begin);
        for (int which = 0; which < 2; ++which) {
            const hicx::ReadFields& mate = which == 0 ? mate1 : mate2;
            const std::int64_t offset =
                std::max<std::int64_t>(0, static_cast<std::int64_t>(mate.pos) -
                                              static_cast<std::int64_t>(bin.begin));
            const std::int64_t vec_start = offset / 10;
            const std::int64_t vec_end =
                std::min<std::int64_t>(length_coverage, vec_start + mate.seq_len / 10);
            const std::int64_t from = static_cast<std::int64_t>(cov_begin) + vec_start;
            const std::int64_t to = static_cast<std::int64_t>(cov_begin) + vec_end;
            for (std::int64_t i = from;
                 i < to && i < static_cast<std::int64_t>(context.coverage_size); ++i) {
                if (i >= 0) {
                    context.coverage[i].fetch_add(1, std::memory_order_relaxed);
                }
            }
        }

        if (!context.quick_qc_mode) {
            const std::uint32_t low = std::min(mate_bins[0], mate_bins[1]);
            const std::uint32_t high = std::max(mate_bins[0], mate_bins[1]);
            out.pixels.push_back((static_cast<std::uint64_t>(low) << 32) | high);
        }
        ++out.counters.pair_added;
        if (want_bam) {
            out.bam_indices.push_back(static_cast<std::int64_t>(index));
        }
    }
}

// --- pixel accumulator ------------------------------------------------------
//
// The counts of one chunk are sorted and coalesced, then merged into a running
// sorted (key, count) accumulator. Memory is O(distinct pixels) rather than
// O(read pairs), which is what keeps the port inside the budget of
// cpp/STATUS.md on a real library, and the merge is order independent because
// it merges sorted keys and adds integers.
class PixelAccumulator {
  public:
    void absorb(std::vector<std::uint64_t>& chunk) {
        if (chunk.empty()) {
            return;
        }
        std::sort(chunk.begin(), chunk.end());
        std::vector<std::uint64_t> keys;
        std::vector<std::int64_t> counts;
        keys.reserve(chunk.size());
        counts.reserve(chunk.size());
        for (const std::uint64_t key : chunk) {
            if (!keys.empty() && keys.back() == key) {
                ++counts.back();
            } else {
                keys.push_back(key);
                counts.push_back(1);
            }
        }
        merge(keys, counts);
        chunk.clear();
        chunk.shrink_to_fit();
    }

    [[nodiscard]] const std::vector<std::uint64_t>& keys() const { return keys_; }
    [[nodiscard]] const std::vector<std::int64_t>& counts() const { return counts_; }

  private:
    void merge(const std::vector<std::uint64_t>& keys,
               const std::vector<std::int64_t>& counts) {
        if (keys_.empty()) {
            keys_ = keys;
            counts_ = counts;
            return;
        }
        std::vector<std::uint64_t> out_keys;
        std::vector<std::int64_t> out_counts;
        out_keys.reserve(keys_.size() + keys.size());
        out_counts.reserve(keys_.size() + keys.size());
        std::size_t a = 0;
        std::size_t b = 0;
        while (a < keys_.size() && b < keys.size()) {
            if (keys_[a] < keys[b]) {
                out_keys.push_back(keys_[a]);
                out_counts.push_back(counts_[a]);
                ++a;
            } else if (keys[b] < keys_[a]) {
                out_keys.push_back(keys[b]);
                out_counts.push_back(counts[b]);
                ++b;
            } else {
                out_keys.push_back(keys_[a]);
                out_counts.push_back(counts_[a] + counts[b]);
                ++a;
                ++b;
            }
        }
        for (; a < keys_.size(); ++a) {
            out_keys.push_back(keys_[a]);
            out_counts.push_back(counts_[a]);
        }
        for (; b < keys.size(); ++b) {
            out_keys.push_back(keys[b]);
            out_counts.push_back(counts[b]);
        }
        keys_.swap(out_keys);
        counts_.swap(out_counts);
    }

    std::vector<std::uint64_t> keys_;
    std::vector<std::int64_t> counts_;
};

hicx::ChromSizes read_chromosome_sizes(const std::string& path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("could not open " + path);
    }
    // The Python builds an OrderedDict, so a repeated name keeps its first
    // position and its last size.
    hicx::ChromSizes sizes;
    std::string line;
    while (std::getline(file, line)) {
        while (!line.empty() &&
               (line.back() == '\r' || line.back() == '\n' || line.back() == ' ' ||
                line.back() == '\t')) {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        const std::size_t tab = line.find('\t');
        if (tab == std::string::npos) {
            continue;
        }
        const std::string name = line.substr(0, tab);
        const std::int64_t size = std::stoll(line.substr(tab + 1));
        const auto found = std::find_if(
            sizes.begin(), sizes.end(),
            [&name](const auto& entry) { return entry.first == name; });
        if (found != sizes.end()) {
            found->second = size;
        } else {
            sizes.emplace_back(name, size);
        }
    }
    return sizes;
}

// The matrix file of a finished run: the upper triangle of the pixel counts
// over the enlarged bins, written as h5, cool or mcool the way
// buildMatrixMethods.createMatrix writes them. Shared by the BAM route and the
// .pairs route (build_matrix_pairs_impl.hpp). There are `pixel_count` distinct
// pixels, which for_each_pixel(sink) hands to sink(key, count) in ascending
// key order, key = min_bin << 32 | max_bin; a source may release its storage
// while it does so. `bin_max` is the coverage column of the h5 bin table;
// `genome_assembly` may be empty.
template <class ForEachPixel>
void write_matrix_outputs(const Arguments& args,
                          const std::vector<hicx::GenomeInterval>& bin_intervals,
                          const std::vector<double>& bin_max, const hicx::ChromNames& names,
                          std::size_t pixel_count, ForEachPixel&& for_each_pixel,
                          const std::string& qc_log, const std::string& genome_assembly) {
    const std::int64_t matrix_size = static_cast<std::int64_t>(bin_intervals.size());
    hicx::MatrixData data;
    data.cut_intervals.reserve(bin_intervals.size());
    for (std::size_t i = 0; i < bin_intervals.size(); ++i) {
        hicx::CutInterval interval;
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
        data.matrix = hicx::CsrMatrix(matrix_size, matrix_size, std::move(indptr),
                                      std::move(indices), std::move(values), "int64");
        data.matrix.set_symmetry(hicx::Symmetry::UpperTriangle);
    }

    // Q9 buildMatrixMethods.py:1371-1376. The three provenance strings are
    // wrapped in np.string_, that is numpy bytes, and hicmatrix's cool
    // writer stores str(value) of them (hicmatrix/lib/cool.py:394, 398,
    // 401). str() of a bytes object is its repr, so what lands in the
    // cool file is the seven extra characters of "b'...'" around the
    // value. Reproduced, not fixed: the attribute is part of the file.
    const auto as_python_bytes_repr = [](const std::string& value) {
        return "b'" + value + "'";
    };
    std::map<std::string, std::string> metadata;
    metadata["statistics"] = qc_log;
    metadata["matrix-generated-by"] =
        as_python_bytes_repr(std::string("HiCExplorer-") + hicx::kVersion);
    metadata["matrix-generated-by-url"] =
        as_python_bytes_repr("https://github.com/deeptools/HiCExplorer");
    if (!genome_assembly.empty()) {
        metadata["genome-assembly"] = as_python_bytes_repr(genome_assembly);
    }

    std::error_code ec;
    std::filesystem::remove(args.out_file_name, ec);

    if (ends_with(args.out_file_name, ".mcool") && args.bin_size.size() > 2) {
        // Q7: only three or more resolutions take this branch.
        hicx::CoolSaveOptions options;
        options.symmetric = true;
        options.apply_correction = false;
        options.hic_metadata = metadata;
        options.has_hic_metadata = true;
        hicx::write_cool(args.out_file_name + "::/resolutions/" +
                             std::to_string(args.bin_size[0]),
                         data, options);
        // Q10 hicmatrix/lib/cool.py:394-402. create_cooler_input deletes
        // 'matrix-generated-by', 'matrix-generated-by-url' and
        // 'genome-assembly' out of the metadata dictionary it was handed,
        // and buildMatrixMethods.py:1381-1404 hands the *same* dictionary
        // to every resolution. So only the first resolution of an mcool
        // carries the provenance and the rest carry none. Reproduced.
        options.hic_metadata.erase("matrix-generated-by");
        options.hic_metadata.erase("matrix-generated-by-url");
        options.hic_metadata.erase("genome-assembly");
        for (std::size_t r = 1; r < args.bin_size.size(); ++r) {
            const std::int64_t factor = args.bin_size[r] / args.bin_size[0];
            hicx::MatrixData merged = hicx::merge_bins(data, factor);
            hicx::CoolSaveOptions append = options;
            append.append = true;
            hicx::write_cool(args.out_file_name + "::/resolutions/" +
                                 std::to_string(args.bin_size[r]),
                             merged, append);
        }
    } else if (ends_with(args.out_file_name, "h5")) {
        hicx::write_hicexplorer_h5(args.out_file_name, data);
    } else if (ends_with(args.out_file_name, "cool")) {
        hicx::CoolSaveOptions options;
        options.symmetric = true;
        options.apply_correction = false;
        options.hic_metadata = metadata;
        options.has_hic_metadata = true;
        hicx::write_cool(args.out_file_name, data, options);
    }
    // A name ending in neither is written nowhere at all, silently, which
    // is hiCMatrix.save's behaviour.
}

int run_build_matrix(const Arguments& args) {
    try {
        if (!ends_with(args.out_file_name, ".h5") &&
            !ends_with(args.out_file_name, ".cool") &&
            args.out_file_name.find(".mcool") == std::string::npos) {
            std::fprintf(stderr,
                         "ERROR:hicexplorer.lib.buildMatrixMethods:Please define the "
                         "file extension. h5 and cool are supported, or the "
                         "specializations of cool, mcool. Given input %s\n",
                         args.out_file_name.c_str());
            return 1;
        }
        std::int64_t max_library_insert_size = args.max_library_insert_size;
        if (args.max_distance.has_value()) {
            max_library_insert_size = *args.max_distance;  // backwards compatibility
        }
        std::error_code ec;
        std::filesystem::create_directories(args.qc_folder, ec);
        if (!std::filesystem::is_directory(args.qc_folder)) {
            std::fprintf(stderr, "Can't open/create QC folder path: %s. Please check\n",
                         args.qc_folder.c_str());
            return 1;
        }
        std::int64_t threads = args.threads;
        if (threads < 2) {
            threads = 2;
            std::fputs("\nAt least two threads need to be defined. Setting "
                       "--threads = 2!s\n",
                       stderr);
        }
        if (args.restriction_sequences.size() != args.dangling_sequences.size()) {
            // buildMatrixMethods.py:968-970 indexes pDanglingSequence with the
            // restriction sequence index, so a short list is an IndexError.
            std::fprintf(stderr,
                         "%s: --danglingSequence has %zu entries but "
                         "--restrictionSequence has %zu; they are matched by "
                         "position\n",
                         g_tool, args.dangling_sequences.size(),
                         args.restriction_sequences.size());
            return 1;
        }

        hicx::BamReader reader1(args.sam_files[0]);
        hicx::BamReader reader2(args.sam_files[1]);
        std::fprintf(stderr,
                     "INFO:hicexplorer.lib.buildMatrixMethods:reading %s and %s to "
                     "build hic_matrix\n\n",
                     args.sam_files[0].c_str(), args.sam_files[1].c_str());

        std::unique_ptr<hicx::BamWriter> out_bam;
        if (!args.do_test_run && !args.out_bam.empty()) {
            out_bam = std::make_unique<hicx::BamWriter>(args.out_bam, reader1);
        }

        const hicx::ChromSizes chrom_sizes =
            args.chromosome_sizes.empty() ? hicx::chrom_sizes_of(reader1)
                                          : read_chromosome_sizes(args.chromosome_sizes);

        hicx::ChromNames names;
        // Intern the reference names first so that tid_to_chrom is a direct
        // lookup and the id order matches the BAM header.
        for (const auto& reference : reader1.references()) {
            names.intern(reference);
        }
        for (const auto& entry : chrom_sizes) {
            names.intern(entry.first);
        }

        std::vector<hicx::GenomeInterval> rf_interval;
        for (const auto& path : args.restriction_cut_files) {
            hicx::bed2interval_list(path, chrom_sizes, args.region, names, rf_interval);
        }
        std::vector<hicx::GenomeInterval> bin_intervals;
        if (!args.bin_size.empty()) {
            bin_intervals =
                hicx::get_bins(args.bin_size[0], chrom_sizes, args.region, names);
        } else {
            bin_intervals =
                hicx::get_rf_bins(rf_interval, args.min_distance, max_library_insert_size);
        }
        const std::int64_t matrix_size = static_cast<std::int64_t>(bin_intervals.size());
        const hicx::BinSearchIndex bin_index(bin_intervals, names);
        const hicx::RestrictionSiteIndex rf_index(rf_interval, names.size());
        const bool rf_present = !rf_interval.empty();
        rf_interval.clear();
        rf_interval.shrink_to_fit();

        Context context;
        context.bins = &bin_index;
        context.rf = &rf_index;
        context.rf_present = rf_present;
        context.tid_to_chrom.reserve(reader1.references().size());
        for (const auto& reference : reader1.references()) {
            context.tid_to_chrom.push_back(*names.lookup(reference));
        }
        for (std::size_t i = 0; i < args.restriction_sequences.size(); ++i) {
            const std::string restriction = to_upper(args.restriction_sequences[i]);
            const std::string dangling = to_upper(args.dangling_sequences[i]);
            context.restriction_sequences.push_back(restriction);
            context.pat_forw.push_back(dangling);
            context.pat_rev.push_back(hicx::reverse_complement(dangling));
        }
        context.has_dangling = !args.dangling_sequences.empty();
        context.keep_self_circles = args.keep_self_circles;
        context.keep_self_ligation = args.keep_self_ligation;
        context.max_insert_size = max_library_insert_size;
        context.quick_qc_mode = args.do_test_run;

        // Coverage vectors, ten base pairs per cell.
        context.coverage_begin.resize(bin_intervals.size());
        context.coverage_end.resize(bin_intervals.size());
        std::int64_t coverage_elements = 0;
        for (std::size_t i = 0; i < bin_intervals.size(); ++i) {
            context.coverage_begin[i] = static_cast<std::uint32_t>(coverage_elements);
            coverage_elements += (bin_intervals[i].end - bin_intervals[i].start) / 10;
            context.coverage_end[i] = static_cast<std::uint32_t>(coverage_elements - 1);
        }
        std::vector<std::atomic<std::uint32_t>> coverage(
            static_cast<std::size_t>(std::max<std::int64_t>(0, coverage_elements)));
        for (auto& cell : coverage) {
            cell.store(0, std::memory_order_relaxed);
        }
        context.coverage = coverage.data();
        context.coverage_size = coverage.size();

        std::int64_t input_buffer_size = args.input_buffer_size;
        if (args.do_test_run) {
            input_buffer_size = args.do_test_run_lines;
        }
        const std::size_t worker_count =
            static_cast<std::size_t>(std::max<std::int64_t>(1, threads - 1));

        hicx::QcCounters total;
        total.dangling_end.assign(context.restriction_sequences.size(), 0);
        PixelAccumulator accumulator;
        std::unordered_set<DuplicateKey, DuplicateKeyHash> seen;

        // The chromosome rank by name, which is the ordering the Python's
        // string key imposes on the duplicate check.
        std::vector<std::uint32_t> name_rank(reader1.references().size());
        {
            std::vector<std::size_t> order(reader1.references().size());
            for (std::size_t i = 0; i < order.size(); ++i) {
                order[i] = i;
            }
            std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
                return reader1.references()[a] < reader1.references()[b];
            });
            for (std::size_t r = 0; r < order.size(); ++r) {
                name_rank[order[r]] = static_cast<std::uint32_t>(r);
            }
        }

        PairBuffer buffer;
        hicx::BamRecord record1;
        hicx::BamRecord record2;
        bool all_data_read = false;
        const bool want_bam = out_bam != nullptr;

        while (!all_data_read) {
            buffer.clear();
            std::int64_t discarded = 0;
            while (static_cast<std::int64_t>(buffer.size()) < input_buffer_size) {
                if (!reader1.read(record1) || !reader2.read(record2)) {
                    all_data_read = true;
                    break;
                }
                ++discarded;  // iter_num in readBamFiles, corrected below
                while ((record1.flag() & 0x100) == 0x100) {
                    if (!reader1.read(record1)) {
                        all_data_read = true;
                        break;
                    }
                }
                while ((record2.flag() & 0x100) == 0x100) {
                    if (!reader2.read(record2)) {
                        all_data_read = true;
                        break;
                    }
                }
                if (record1.qname() != record2.qname()) {
                    throw std::runtime_error(
                        "FATAL ERROR " + record1.qname() + " " + record2.qname() +
                        " Be sure that the sam files have the same read order If "
                        "using Bowtie2 or Hisat2 add the --reorder option");
                }

                // get_supplementary_alignment plus get_correct_map. No record
                // in the test corpus carries an SA tag, so this path has never
                // been exercised against the Python on real data; it is
                // written from the Python and marked in the report.
                auto resolve_supplementary = [](hicx::BamReader& reader,
                                                hicx::BamRecord& primary) {
                    if (!primary.has_tag("SA")) {
                        return;
                    }
                    const std::string tag = primary.tag_string("SA");
                    std::size_t entries = 0;
                    for (const char c : tag) {
                        if (c == ';') {
                            ++entries;
                        }
                    }
                    std::int64_t best = primary.leading_bases_before_match();
                    hicx::BamRecord best_record = primary;
                    hicx::BamRecord candidate;
                    for (std::size_t i = 0; i < entries; ++i) {
                        if (!reader.read(candidate)) {
                            break;
                        }
                        if (candidate.qname() != primary.qname()) {
                            continue;
                        }
                        const std::int64_t value =
                            candidate.leading_bases_before_match();
                        if (value < best) {
                            best = value;
                            best_record = candidate;
                        }
                    }
                    primary = best_record;
                };
                resolve_supplementary(reader1, record1);
                resolve_supplementary(reader2, record2);

                if ((record1.flag() & 0x4) == 4 || (record2.flag() & 0x4) == 4) {
                    ++total.one_mate_unmapped;
                    continue;
                }
                if (record1.mapq() < args.min_mapping_quality ||
                    record2.mapq() < args.min_mapping_quality) {
                    // Q1: the Python's chained comparison looks at mate1 only.
                    if (record1.mapq() == 0) {
                        ++total.one_mate_not_unique;
                        continue;
                    }
                    ++total.one_mate_low_quality;
                    continue;
                }
                if (!args.skip_duplication_check) {
                    const std::uint32_t rank1 =
                        record1.tid() >= 0
                            ? name_rank[static_cast<std::size_t>(record1.tid())]
                            : 0;
                    const std::uint32_t rank2 =
                        record2.tid() >= 0
                            ? name_rank[static_cast<std::size_t>(record2.tid())]
                            : 0;
                    DuplicateKey key;
                    key.chrom_low = std::min(rank1, rank2);
                    key.chrom_high = std::max(rank1, rank2);
                    key.pos_low = std::min(record1.pos(), record2.pos());
                    key.pos_high = std::max(record1.pos(), record2.pos());
                    if (!seen.insert(key).second) {
                        ++total.duplicated_pairs;
                        continue;
                    }
                }
                buffer.mate1.push_back(record1.fields());
                buffer.mate2.push_back(record2.fields());
                if (want_bam) {
                    buffer.record1.push_back(record1);
                    buffer.record2.push_back(record2);
                }
            }
            discarded -= static_cast<std::int64_t>(buffer.size());
            total.iter_num += discarded + static_cast<std::int64_t>(buffer.size());
            if (buffer.size() == 0) {
                break;
            }

            // Fixed contiguous partitions, reduced sequentially inside a
            // partition and combined in index order.
            const std::size_t pairs = buffer.size();
            const std::size_t parts = std::min<std::size_t>(worker_count, pairs);
            std::vector<Partition> partitions(parts);
            std::vector<std::thread> workers;
            workers.reserve(parts);
            for (std::size_t p = 0; p < parts; ++p) {
                const std::size_t first = pairs * p / parts;
                const std::size_t last = pairs * (p + 1) / parts;
                if (parts == 1) {
                    classify_range(context, buffer, first, last, want_bam,
                                   partitions[p]);
                } else {
                    workers.emplace_back([&, p, first, last] {
                        classify_range(context, buffer, first, last, want_bam,
                                       partitions[p]);
                    });
                }
            }
            for (auto& worker : workers) {
                worker.join();
            }

            std::vector<std::uint64_t> chunk_pixels;
            for (auto& partition : partitions) {
                total.add(partition.counters);
                chunk_pixels.insert(chunk_pixels.end(), partition.pixels.begin(),
                                    partition.pixels.end());
                partition.pixels.clear();
                partition.pixels.shrink_to_fit();
            }
            accumulator.absorb(chunk_pixels);

            if (want_bam) {
                for (const auto& partition : partitions) {
                    for (const std::int64_t index : partition.bam_indices) {
                        hicx::BamRecord& mate1 =
                            buffer.record1[static_cast<std::size_t>(index)];
                        hicx::BamRecord& mate2 =
                            buffer.record2[static_cast<std::size_t>(index)];
                        mate1.set_paired_flag_first();
                        mate2.set_paired_flag_second();
                        mate1.set_mate(mate2.tid(), mate2.pos());
                        mate2.set_mate(mate1.tid(), mate1.pos());
                        out_bam->write(mate1);
                        out_bam->write(mate2);
                    }
                }
            }

            if (args.do_test_run && total.iter_num > args.do_test_run_lines) {
                break;
            }
        }
        buffer.clear();
        seen.clear();

        if (out_bam) {
            out_bam->close();
        }

        // --- QC ------------------------------------------------------------
        hicx::QcLogInputs qc_inputs;
        qc_inputs.out_file_name = args.out_file_name;
        qc_inputs.min_distance = args.min_distance;
        qc_inputs.max_library_insert_size = max_library_insert_size;
        qc_inputs.keep_self_ligation = args.keep_self_ligation;
        qc_inputs.has_restriction_cut_file = !args.restriction_cut_files.empty();
        for (std::size_t i = 0; i < context.restriction_sequences.size(); ++i) {
            qc_inputs.dangling_sequences.emplace_back(context.restriction_sequences[i],
                                                      context.pat_forw[i]);
        }
        const std::string qc_log = hicx::format_qc_log(qc_inputs, total);
        {
            std::ofstream log_file(args.qc_folder + "/QC.log", std::ios::binary);
            if (!log_file) {
                throw std::runtime_error("could not write " + args.qc_folder +
                                         "/QC.log");
            }
            log_file << qc_log;
        }
        hicx::write_qc_tables(args.qc_folder, qc_log);

        if (args.do_test_run) {
            hicx::report_resource_usage(g_tool);
            return 0;
        }

        // --- the matrix ------------------------------------------------------
        hicx::enlarge_bins(bin_intervals, chrom_sizes, names);
        std::vector<double> bin_max(bin_intervals.size());
        for (std::size_t i = 0; i < bin_intervals.size(); ++i) {
            std::uint32_t max_element = 0;
            const std::int64_t from = context.coverage_begin[i];
            const std::int64_t to = context.coverage_end[i];  // Q4: exclusive end
            for (std::int64_t k = from;
                 k < to && k < static_cast<std::int64_t>(coverage.size()); ++k) {
                const std::uint32_t value = coverage[static_cast<std::size_t>(k)].load(
                    std::memory_order_relaxed);
                if (value > max_element) {
                    max_element = value;
                }
            }
            bin_max[i] = max_element == 0 ? std::numeric_limits<double>::quiet_NaN()
                                          : static_cast<double>(max_element);
        }
        coverage.clear();

        (void)matrix_size;
        write_matrix_outputs(
            args, bin_intervals, bin_max, names, accumulator.keys().size(),
            [&accumulator](auto&& sink) {
                const auto& keys = accumulator.keys();
                const auto& counts = accumulator.counts();
                for (std::size_t k = 0; k < keys.size(); ++k) {
                    sink(keys[k], static_cast<std::uint64_t>(counts[k]));
                }
            },
            qc_log, args.genome_assembly);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s: %s\n", g_tool, error.what());
        return 1;
    }
    hicx::report_resource_usage(g_tool);
    return 0;
}

// buildMatrixMethods.py:1362, QC.main(["-l", QC.log, "-o", QCfolder]): the
// tables are written by run_build_matrix, and the five PNGs and hicQC.html are
// drawn here with hicPrepareQCreport's default --dpi 200. Called by the tools
// after everything else, because on success it does not return.
int draw_qc_report(const std::string& qc_folder) {
    hicx::plot::JsonObject data;
    data.add("outputFolder", hicx::plot::json_string(qc_folder));
    data.add("dpi", hicx::plot::json_int(200));
    return hicx::plot::draw("hicPrepareQCreport", data.str(), std::nullopt);
}

}  // namespace

#endif  // HICX_TOOLS_BUILD_MATRIX_IMPL_HPP
