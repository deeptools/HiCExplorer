// hicBuildMatrix --pairsFile: a contact matrix from a 4DN or pairtools .pairs
// file instead of two BAM files (cpp/PLAN.md tier 9, item 9.2, class EX).
//
// HiCExplorer 3.7 has no such input, so there is no Python oracle. The route
// is validated E2 against `cooler cload pairs` on the same file and bins, and
// E2 against the BAM route on a .pairs file that carries exactly the read
// pairs the BAM route keeps (cpp/scripts/make_bam_route_pairs.py).
//
// ---------------------------------------------------------------------------
// Format (hicx/pairs_file.hpp)
// ---------------------------------------------------------------------------
//  F1 Plain, gzip or bgzip text, read through htslib. bgzip blocks are
//     decompressed on --threads threads.
//  F2 The first line must be "## pairs format v1...". Of the header lines,
//     #sorted, #shape, #chromsize, #columns and #genome_assembly are read and
//     every other one is ignored. Without #columns the seven reserved columns
//     readID chrom1 pos1 chrom2 pos2 strand1 strand2 are assumed. The column
//     names chr1 and chr2, which ENCODE's files carry, mean chrom1 and chrom2.
//  F3 Columns read: chrom1 pos1 chrom2 pos2 (required), strand1 strand2,
//     mapq1 mapq2, pair_type. frag1, frag2 and every other column are passed
//     over; hicBuildMatrix builds its own restriction fragment bins.
//  F4 A position is the one-based 5' end of its read. The pair is binned at
//     position - 1 into the bin whose half-open interval [start, end) holds
//     it, which is how cooler cload pairs bins a pair. (The BAM route bins the
//     read middle into a closed interval; see make_bam_route_pairs.py for how
//     the two meet at a shared bin end.)
//  F5 Chromosome sizes: --chromosomeSizes, else the #chromsize lines in their
//     order; without either the run is refused. A pair on a chromosome without
//     bins, or at a position outside every bin (in a gap between restriction
//     fragment bins, before --region, beyond the chromosome end), is not
//     counted, as the BAM route drops a mate that finds no bin.
//  F6 Sorted or unsorted files are both accepted; only the duplicate check
//     depends on the order (P4).
//
// ---------------------------------------------------------------------------
// Filters: which of hicBuildMatrix's filters apply to a pair of positions
// ---------------------------------------------------------------------------
//  P1 Unmapped: chromosome "!", or a pair_type letter other than U, R and M
//     (pairtools' N null, X corrupt, W walk), counts as "One mate unmapped".
//  P2 Not unique: a pair_type letter M counts as "One mate not unique".
//  P3 --minMappingQuality applies when the file has mapq1 and mapq2: a pair
//     with either below it counts as "One mate not unique" when either is 0
//     and as "Low mapping quality" otherwise. This deliberately does not copy
//     Q1 of the BAM route, whose chained comparison looks at the first mate
//     only. Without the two columns an explicit --minMappingQuality is refused
//     (exit 1, before any output), and the default is not applied, with a
//     note on stderr.
//  P4 Duplicates. pair_type DD (pairtools dedup --mark-dups) counts as
//     "duplicated pairs" and is removed, also with --skipDuplicationCheck,
//     which skips only the tool's own check. The own check removes a pair
//     whose two ends, each a (chromosome, position), equal those of an earlier
//     pair in either order. (The BAM route's key orders the two chromosomes
//     and the two positions separately, so it also joins pairs whose
//     positions are swapped between chromosomes; the pairs key does not.)
//     A file whose header declares "#sorted: chr1-chr2-pos1-pos2" and
//     "#shape: upper triangle", as pairtools sort writes, is checked in a
//     stream against the previous pair only, in constant memory, and the
//     declared order is verified: a pair out of that order is an error. Any
//     other file is checked against a hash set, whose memory grows with the
//     number of pairs, as the BAM route's does.
//  P5 Refused at the command line (exit 2), because a .pairs file carries
//     neither read sequence nor alignment spans: --restrictionSequence and
//     --danglingSequence (a dangling end is found in the read sequence),
//     --keepSelfLigation and --keepSelfCircles (self ligation, self circle and
//     same fragment are decided on the aligned spans of both reads), --outBam
//     (there are no alignments to write). pairtools select and pairtools
//     restrict are the tools for those decisions on pairs.
//  P6 Bins: --binSize gives fixed bins. Without it, --restrictionCutFile,
//     --minDistance and --maxLibraryInsertSize (or --maxDistance) give the
//     BAM route's restriction fragment bins. With --binSize those four
//     options have no use and are refused (exit 2); with neither, the run is
//     refused.
//  P7 As in the BAM route: --region, --chromosomeSizes, --threads,
//     --inputBufferSize (lines per chunk), --doTestRun (the first
//     --doTestRunLines lines), --skipDuplicationCheck. --genomeAssembly,
//     else the file's #genome_assembly.
//
// ---------------------------------------------------------------------------
// Outputs
// ---------------------------------------------------------------------------
//  O1 h5, cool and mcool through the BAM route's writer
//     (write_matrix_outputs), with its Q7, Q9 and Q10.
//  O2 The h5 bin table's fourth column, the per bin maximum of the BAM
//     route's read coverage, needs read lengths. Every bin gets NaN, the value
//     the BAM route writes for a bin no read covers.
//  O3 QC.log, the five tables and the QC report in the layout of
//     hicBuildMatrixMicroC, which hicPrepareQCreport already reads: no
//     dangling end, self ligation, self circle or rest site lines, and
//     "same fragment" 0. "Sequenced reads" counts pair lines. The distance
//     classes use the pairs' positions; the four orientation counts need
//     strand1 and strand2 and count pairs whose two strands are + or -.
//     "Max library insert size" is printed as in hicBuildMatrixMicroC and "Min
//     rest. site distance" only for restriction fragment bins.
//
// ---------------------------------------------------------------------------
// Structure and determinism
// ---------------------------------------------------------------------------
// A chunk of --inputBufferSize lines is parsed on --threads threads over
// contiguous ranges, filtered and checked for duplicates serially in file
// order, and binned over contiguous ranges whose counters and pixels are
// combined in index order: the result is the same for every thread count.
// Pixels are sorted in runs of at most 2^24 keys, kept as 12 bytes per
// distinct pixel in blocks that merges release as they read them
// (PairsPixelRuns), and drained block by block into the CSR matrix: memory is
// about 24 bytes per distinct pixel at the hand-over, against the 32 bytes of
// an eight byte key and an eight byte count merged twice over.

#ifndef HICX_TOOLS_BUILD_MATRIX_PAIRS_IMPL_HPP
#define HICX_TOOLS_BUILD_MATRIX_PAIRS_IMPL_HPP

#include <deque>
#include <functional>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

#include "build_matrix_impl.hpp"
#include "hicx/pairs_file.hpp"

namespace {

constexpr std::uint32_t kUnknownChrom = std::numeric_limits<std::uint32_t>::max();

// Chromosome name to the id hicx::ChromNames gives it. Read concurrently while
// a chunk is parsed; interned into only between parses.
class PairsChromIds {
  public:
    explicit PairsChromIds(hicx::ChromNames& names) : names_(names) {
        for (std::uint32_t id = 0; id < names.size(); ++id) {
            add(names.name(id), id);
        }
    }

    [[nodiscard]] std::uint32_t find(std::string_view name) const {
        const auto found = ids_.find(name);
        return found == ids_.end() ? kUnknownChrom : found->second;
    }

    std::uint32_t intern(std::string_view name) {
        const std::uint32_t known = find(name);
        if (known != kUnknownChrom) {
            return known;
        }
        const std::uint32_t id = names_.intern(std::string(name));
        add(names_.name(id), id);
        return id;
    }

  private:
    void add(const std::string& name, std::uint32_t id) {
        storage_.push_back(name);
        ids_.emplace(std::string_view(storage_.back()), id);
    }

    hicx::ChromNames& names_;
    std::deque<std::string> storage_;  // stable addresses for the views
    std::unordered_map<std::string_view, std::uint32_t> ids_;
};

// The bins of each chromosome sorted by start, looked up half-open (F4).
class PairsBinIndex {
  public:
    PairsBinIndex(const std::vector<hicx::GenomeInterval>& bins, std::size_t chrom_count) {
        chroms_.resize(chrom_count);
        std::vector<std::vector<std::size_t>> members(chrom_count);
        for (std::size_t i = 0; i < bins.size(); ++i) {
            members[bins[i].chrom].push_back(i);
        }
        for (std::size_t chrom = 0; chrom < chrom_count; ++chrom) {
            std::vector<std::size_t>& group = members[chrom];
            std::sort(group.begin(), group.end(), [&bins](std::size_t a, std::size_t b) {
                if (bins[a].start != bins[b].start) {
                    return bins[a].start < bins[b].start;
                }
                if (bins[a].end != bins[b].end) {
                    return bins[a].end < bins[b].end;
                }
                return a < b;
            });
            Chrom& target = chroms_[chrom];
            for (const std::size_t i : group) {
                target.start.push_back(bins[i].start);
                target.end.push_back(bins[i].end);
                target.id.push_back(static_cast<std::uint32_t>(i));
            }
        }
    }

    // The bin holding the zero-based `position`, or nullopt.
    [[nodiscard]] std::optional<std::uint32_t> bin_of(std::uint32_t chrom,
                                                      std::int64_t position) const {
        if (chrom >= chroms_.size() || position < 0) {
            return std::nullopt;
        }
        const Chrom& target = chroms_[chrom];
        const auto upper = std::upper_bound(target.start.begin(), target.start.end(), position);
        if (upper == target.start.begin()) {
            return std::nullopt;
        }
        const auto k = static_cast<std::size_t>(upper - target.start.begin()) - 1;
        if (position >= target.end[k]) {
            return std::nullopt;
        }
        return target.id[k];
    }

  private:
    struct Chrom {
        std::vector<std::int64_t> start;
        std::vector<std::int64_t> end;
        std::vector<std::uint32_t> id;
    };
    std::vector<Chrom> chroms_;
};

struct ParsedPair {
    std::uint32_t chrom1 = kUnknownChrom;  // kUnknownChrom: not interned yet
    std::uint32_t chrom2 = kUnknownChrom;
    std::int64_t pos1 = 0;  // one-based, as written
    std::int64_t pos2 = 0;
    std::int32_t mapq1 = -1;
    std::int32_t mapq2 = -1;
    std::int8_t strand1 = 0;
    std::int8_t strand2 = 0;
    hicx::PairKind kind = hicx::PairKind::Mapped;
};

// The duplicate check of P4.
class PairsDuplicateCheck {
  public:
    explicit PairsDuplicateCheck(bool streaming) : streaming_(streaming) {}

    [[nodiscard]] bool streaming() const noexcept { return streaming_; }

    // Whether the pair repeats an earlier one. Throws in streaming mode when
    // the pair breaks the order the header declares.
    bool duplicate(const ParsedPair& pair, std::int64_t line_number, const std::string& path) {
        if (!streaming_) {
            Key key{pair.chrom1, pair.chrom2, pair.pos1, pair.pos2};
            if (std::tie(key.chrom2, key.pos2) < std::tie(key.chrom1, key.pos1)) {
                std::swap(key.chrom1, key.chrom2);
                std::swap(key.pos1, key.pos2);
            }
            return !seen_.insert(key).second;
        }
        const auto order_error = [&](const std::string& what) {
            return std::runtime_error(
                path + " line " + std::to_string(line_number) + ": " + what +
                ", although the header declares '#sorted: chr1-chr2-pos1-pos2' and "
                "'#shape: upper triangle'. Sort the file (pairtools sort) or give "
                "--skipDuplicationCheck");
        };
        if (pair.chrom1 == pair.chrom2 && pair.pos1 > pair.pos2) {
            throw order_error("pos1 is larger than pos2 on one chromosome");
        }
        const std::uint64_t block = (static_cast<std::uint64_t>(pair.chrom1) << 32) | pair.chrom2;
        if (!has_previous_ || block != previous_block_) {
            const std::uint64_t flipped =
                (static_cast<std::uint64_t>(pair.chrom2) << 32) | pair.chrom1;
            if (blocks_.count(block) != 0) {
                throw order_error("the chromosome pair of this line appeared before and was "
                                  "followed by another one");
            }
            if (pair.chrom1 != pair.chrom2 && blocks_.count(flipped) != 0) {
                throw order_error("this chromosome pair appeared before in the other order");
            }
            blocks_.insert(block);
            previous_block_ = block;
            previous_pos1_ = pair.pos1;
            previous_pos2_ = pair.pos2;
            has_previous_ = true;
            return false;
        }
        if (std::tie(pair.pos1, pair.pos2) < std::tie(previous_pos1_, previous_pos2_)) {
            throw order_error("the positions are smaller than those of the previous pair");
        }
        const bool repeated = pair.pos1 == previous_pos1_ && pair.pos2 == previous_pos2_;
        previous_pos1_ = pair.pos1;
        previous_pos2_ = pair.pos2;
        return repeated;
    }

  private:
    struct Key {
        std::uint32_t chrom1;
        std::uint32_t chrom2;
        std::int64_t pos1;
        std::int64_t pos2;
        bool operator==(const Key& other) const = default;
    };
    struct KeyHash {
        std::size_t operator()(const Key& key) const noexcept {
            std::uint64_t a = (static_cast<std::uint64_t>(key.chrom1) << 32) | key.chrom2;
            const std::uint64_t b = static_cast<std::uint64_t>(key.pos1) * 0x9e3779b97f4a7c15ULL ^
                                    static_cast<std::uint64_t>(key.pos2);
            a ^= b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2);
            return static_cast<std::size_t>(a);
        }
    };

    bool streaming_;
    std::unordered_set<Key, KeyHash> seen_;
    std::unordered_set<std::uint64_t> blocks_;
    std::uint64_t previous_block_ = 0;
    std::int64_t previous_pos1_ = 0;
    std::int64_t previous_pos2_ = 0;
    bool has_previous_ = false;
};

// Pixel keys (min_bin << 32 | max_bin) collected as sorted, coalesced runs of
// (key, count) pairs, 12 bytes per distinct pixel. A run is stored in blocks
// of at most 2^20 entries that a merge releases as it consumes them, so a
// merge, and the final hand-over to the matrix, needs the distinct pixels once
// plus a few blocks, not twice. Merges follow a binary counter (a run is
// merged into its predecessor while that one is at most twice its size), so
// the number of runs stays logarithmic and every key is merged O(log n) times.
class PairsPixelRuns {
  public:
    void add(const std::vector<std::uint64_t>& keys) {
        if (buffer_.size() + keys.size() > kRunKeys) {
            flush();
        }
        buffer_.insert(buffer_.end(), keys.begin(), keys.end());
    }

    // Merges everything into one run; returns its number of distinct pixels.
    std::size_t finish() {
        flush();
        buffer_ = std::vector<std::uint64_t>();
        while (runs_.size() > 1) {
            merge_last_two();
        }
        return runs_.empty() ? 0 : runs_.back().size;
    }

    // After finish(): sink(key, count) for every pixel in key order, releasing
    // the run as it goes.
    template <class Sink>
    void drain(Sink&& sink) {
        if (runs_.empty()) {
            return;
        }
        Run run = std::move(runs_.back());
        runs_.clear();
        for (Cursor cursor{&run}; !cursor.done(); cursor.advance()) {
            sink(cursor.key(), static_cast<std::uint64_t>(cursor.count()));
        }
    }

  private:
    static constexpr std::size_t kRunKeys = std::size_t{1} << 24;
    static constexpr std::size_t kBlock = std::size_t{1} << 20;

    struct Run {
        std::deque<std::vector<std::uint64_t>> keys;
        std::deque<std::vector<std::uint32_t>> counts;
        std::size_t size = 0;

        void push(std::uint64_t key, std::uint64_t count) {
            if (count > std::numeric_limits<std::uint32_t>::max()) {
                throw std::runtime_error("a pixel holds more than 4,294,967,295 pairs");
            }
            if (keys.empty() || keys.back().size() == kBlock) {
                keys.emplace_back();
                counts.emplace_back();
            }
            keys.back().push_back(key);
            counts.back().push_back(static_cast<std::uint32_t>(count));
            ++size;
        }
    };

    // Reads a run front to back and frees each block once it is read.
    struct Cursor {
        Run* run;
        std::size_t position = 0;

        [[nodiscard]] bool done() const { return run->keys.empty(); }
        [[nodiscard]] std::uint64_t key() const { return run->keys.front()[position]; }
        [[nodiscard]] std::uint32_t count() const { return run->counts.front()[position]; }
        void advance() {
            if (++position == run->keys.front().size()) {
                run->keys.pop_front();
                run->counts.pop_front();
                position = 0;
            }
        }
    };

    void flush() {
        if (buffer_.empty()) {
            return;
        }
        std::sort(buffer_.begin(), buffer_.end());
        Run run;
        std::uint64_t count = 0;
        for (std::size_t i = 0; i < buffer_.size(); ++i) {
            ++count;
            if (i + 1 == buffer_.size() || buffer_[i + 1] != buffer_[i]) {
                run.push(buffer_[i], count);
                count = 0;
            }
        }
        buffer_.clear();
        runs_.push_back(std::move(run));
        while (runs_.size() >= 2 && runs_[runs_.size() - 2].size <= 2 * runs_.back().size) {
            merge_last_two();
        }
    }

    void merge_last_two() {
        Run right = std::move(runs_.back());
        runs_.pop_back();
        Run left = std::move(runs_.back());
        runs_.pop_back();
        Run out;
        Cursor a{&left};
        Cursor b{&right};
        while (!a.done() && !b.done()) {
            if (a.key() < b.key()) {
                out.push(a.key(), a.count());
                a.advance();
            } else if (b.key() < a.key()) {
                out.push(b.key(), b.count());
                b.advance();
            } else {
                out.push(a.key(), static_cast<std::uint64_t>(a.count()) + b.count());
                a.advance();
                b.advance();
            }
        }
        for (; !a.done(); a.advance()) {
            out.push(a.key(), a.count());
        }
        for (; !b.done(); b.advance()) {
            out.push(b.key(), b.count());
        }
        runs_.push_back(std::move(out));
    }

    std::vector<std::uint64_t> buffer_;
    std::vector<Run> runs_;
};

// The lines of one chunk, copied out of the reader.
struct PairsLineChunk {
    std::string text;
    std::vector<std::size_t> offset;
    std::vector<std::uint32_t> length;
    std::vector<std::int64_t> number;

    void clear() {
        text.clear();
        offset.clear();
        length.clear();
        number.clear();
    }
    void add(std::string_view line, std::int64_t line_number) {
        offset.push_back(text.size());
        length.push_back(static_cast<std::uint32_t>(line.size()));
        number.push_back(line_number);
        text.append(line);
    }
    [[nodiscard]] std::size_t size() const noexcept { return offset.size(); }
    [[nodiscard]] std::string_view line(std::size_t i) const {
        return {text.data() + offset[i], length[i]};
    }
};

struct PairsPartition {
    hicx::QcCounters counters;
    std::vector<std::uint64_t> pixels;
};

// Runs body(first, last, part) over `parts` contiguous ranges of [0, count),
// on threads when there is more than one range.
void for_contiguous_ranges(std::size_t workers, std::size_t count,
                           const std::function<void(std::size_t, std::size_t, std::size_t)>& body) {
    const std::size_t parts = std::max<std::size_t>(1, std::min(workers, count));
    if (parts == 1) {
        body(0, count, 0);
        return;
    }
    std::vector<std::thread> threads;
    threads.reserve(parts);
    for (std::size_t p = 0; p < parts; ++p) {
        threads.emplace_back([&body, p, parts, count] {
            body(count * p / parts, count * (p + 1) / parts, p);
        });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

void bin_pairs(const PairsBinIndex& index, const std::vector<ParsedPair>& parsed,
               const std::vector<std::uint32_t>& accepted, std::size_t first, std::size_t last,
               bool quick_qc, PairsPartition& out) {
    for (std::size_t k = first; k < last; ++k) {
        const ParsedPair& pair = parsed[accepted[k]];
        const std::optional<std::uint32_t> bin1 = index.bin_of(pair.chrom1, pair.pos1 - 1);
        const std::optional<std::uint32_t> bin2 =
            bin1.has_value() ? index.bin_of(pair.chrom2, pair.pos2 - 1) : std::nullopt;
        if (!bin1.has_value() || !bin2.has_value()) {
            ++out.counters.mate_not_close_to_rf;
            continue;
        }
        if (pair.chrom1 != pair.chrom2) {
            ++out.counters.inter_chromosomal;
        } else {
            if (std::llabs(pair.pos2 - pair.pos1) < 20000) {
                ++out.counters.short_range;
            } else {
                ++out.counters.long_range;
            }
            // The BAM route's rule: the first read is the one with the smaller
            // position, the second read on a tie.
            const bool one_first = pair.pos1 < pair.pos2;
            const std::int8_t first_strand = one_first ? pair.strand1 : pair.strand2;
            const std::int8_t second_strand = one_first ? pair.strand2 : pair.strand1;
            if (first_strand != 0 && second_strand != 0) {
                if (first_strand > 0 && second_strand < 0) {
                    ++out.counters.count_inward;
                } else if (first_strand < 0 && second_strand > 0) {
                    ++out.counters.count_outward;
                } else if (first_strand < 0) {
                    ++out.counters.count_left;
                } else {
                    ++out.counters.count_right;
                }
            }
        }
        if (!quick_qc) {
            const std::uint32_t low = std::min(*bin1, *bin2);
            const std::uint32_t high = std::max(*bin1, *bin2);
            out.pixels.push_back((static_cast<std::uint64_t>(low) << 32) | high);
        }
        ++out.counters.pair_added;
    }
}

int run_build_matrix_pairs(const Arguments& args) {
    try {
        if (!ends_with(args.out_file_name, ".h5") && !ends_with(args.out_file_name, ".cool") &&
            args.out_file_name.find(".mcool") == std::string::npos) {
            std::fprintf(stderr,
                         "ERROR:hicexplorer.lib.buildMatrixMethods:Please define the "
                         "file extension. h5 and cool are supported, or the "
                         "specializations of cool, mcool. Given input %s\n",
                         args.out_file_name.c_str());
            return 1;
        }
        const std::size_t workers =
            static_cast<std::size_t>(std::max<std::int64_t>(1, args.threads));
        hicx::PairsReader reader(args.pairs_file, static_cast<int>(workers));
        const hicx::PairsHeader& header = reader.header();

        const bool has_mapq = header.mapq1 >= 0 && header.mapq2 >= 0;
        if (!has_mapq && args.min_mapping_quality_given) {
            std::fprintf(stderr,
                         "%s: --minMappingQuality needs the columns mapq1 and mapq2, which %s "
                         "does not have\n",
                         g_tool, args.pairs_file.c_str());
            return 1;
        }
        const hicx::ChromSizes chrom_sizes = args.chromosome_sizes.empty()
                                                 ? header.chrom_sizes
                                                 : read_chromosome_sizes(args.chromosome_sizes);
        if (chrom_sizes.empty()) {
            std::fprintf(stderr,
                         "%s: %s has no #chromsize header lines; give the chromosome sizes "
                         "with --chromosomeSizes\n",
                         g_tool, args.pairs_file.c_str());
            return 1;
        }
        std::int64_t max_library_insert_size = args.max_library_insert_size;
        if (args.max_distance.has_value()) {
            max_library_insert_size = *args.max_distance;
        }

        hicx::ChromNames names;
        for (const auto& entry : chrom_sizes) {
            names.intern(entry.first);
        }
        const bool fragment_bins = args.bin_size.empty();
        std::vector<hicx::GenomeInterval> bin_intervals;
        if (!fragment_bins) {
            bin_intervals = hicx::get_bins(args.bin_size[0], chrom_sizes, args.region, names);
        } else {
            std::vector<hicx::GenomeInterval> rf_interval;
            for (const auto& path : args.restriction_cut_files) {
                hicx::bed2interval_list(path, chrom_sizes, args.region, names, rf_interval);
            }
            bin_intervals =
                hicx::get_rf_bins(rf_interval, args.min_distance, max_library_insert_size);
        }

        std::error_code ec;
        std::filesystem::create_directories(args.qc_folder, ec);
        if (!std::filesystem::is_directory(args.qc_folder)) {
            std::fprintf(stderr, "Can't open/create QC folder path: %s. Please check\n",
                         args.qc_folder.c_str());
            return 1;
        }

        const bool streaming = header.sorted_by_chroms_then_positions() && header.upper_triangle();
        std::fprintf(stderr, "INFO:%s:reading %s to build hic_matrix\n", g_tool,
                     args.pairs_file.c_str());
        if (!has_mapq) {
            std::fprintf(stderr,
                         "INFO:%s:%s has no mapq1 and mapq2 columns; no mapping quality filter "
                         "is applied\n",
                         g_tool, args.pairs_file.c_str());
        }
        if (!args.skip_duplication_check) {
            std::fprintf(stderr, "INFO:%s:duplicate check %s\n", g_tool,
                         streaming ? "in file order (sorted, upper triangle)"
                                   : "with a hash set (the file is not declared sorted and "
                                     "upper triangle)");
        }

        PairsChromIds ids(names);
        const PairsBinIndex bin_index(bin_intervals, names.size());
        PairsDuplicateCheck duplicates(streaming);
        PairsPixelRuns pixel_runs;
        hicx::QcCounters total;
        const std::int64_t chunk_lines = std::max<std::int64_t>(
            1, args.do_test_run ? args.do_test_run_lines : args.input_buffer_size);

        PairsLineChunk chunk;
        std::vector<ParsedPair> parsed;
        std::vector<std::uint32_t> accepted;
        std::string_view line;
        bool all_read = false;
        while (!all_read) {
            chunk.clear();
            while (static_cast<std::int64_t>(chunk.size()) < chunk_lines) {
                if (!reader.next(line)) {
                    all_read = true;
                    break;
                }
                chunk.add(line, reader.line_number());
            }
            const std::size_t count = chunk.size();
            if (count == 0) {
                break;
            }

            // Parse, in contiguous ranges.
            parsed.assign(count, ParsedPair{});
            std::vector<std::optional<std::pair<std::size_t, std::string>>> errors(
                std::min(workers, count));
            for_contiguous_ranges(workers, count, [&](std::size_t first, std::size_t last,
                                                      std::size_t part) {
                std::vector<std::string_view> fields;
                hicx::PairsRecord record;
                for (std::size_t i = first; i < last; ++i) {
                    if (auto error = hicx::parse_pairs_record(chunk.line(i), header, fields,
                                                              record)) {
                        errors[part] = std::make_pair(i, std::move(*error));
                        return;
                    }
                    ParsedPair& out = parsed[i];
                    out.kind = record.kind;
                    out.pos1 = record.pos1;
                    out.pos2 = record.pos2;
                    out.mapq1 = record.mapq1;
                    out.mapq2 = record.mapq2;
                    out.strand1 = record.strand1;
                    out.strand2 = record.strand2;
                    if (record.kind == hicx::PairKind::Mapped) {
                        out.chrom1 = ids.find(record.chrom1);
                        out.chrom2 = ids.find(record.chrom2);
                    }
                }
            });
            for (const auto& error : errors) {
                if (error.has_value()) {
                    throw std::runtime_error(args.pairs_file + " line " +
                                             std::to_string(chunk.number[error->first]) +
                                             ": the line " + error->second);
                }
            }

            // Filters and the duplicate check, serially in file order.
            accepted.clear();
            std::vector<std::string_view> fields;
            hicx::PairsRecord record;
            for (std::size_t i = 0; i < count; ++i) {
                ParsedPair& pair = parsed[i];
                ++total.iter_num;
                if (pair.kind == hicx::PairKind::Unmapped) {
                    ++total.one_mate_unmapped;
                    continue;
                }
                if (pair.kind == hicx::PairKind::NotUnique) {
                    ++total.one_mate_not_unique;
                    continue;
                }
                if (has_mapq && (pair.mapq1 < args.min_mapping_quality ||
                                 pair.mapq2 < args.min_mapping_quality)) {
                    if (pair.mapq1 == 0 || pair.mapq2 == 0) {
                        ++total.one_mate_not_unique;
                    } else {
                        ++total.one_mate_low_quality;
                    }
                    continue;
                }
                if (pair.kind == hicx::PairKind::MarkedDuplicate) {
                    ++total.duplicated_pairs;
                    continue;
                }
                if (pair.chrom1 == kUnknownChrom || pair.chrom2 == kUnknownChrom) {
                    // A chromosome without a size: interned here, serially, so
                    // that the id does not depend on the thread count. It has
                    // no bins, so the pair is not counted (F5).
                    (void)hicx::parse_pairs_record(chunk.line(i), header, fields, record);
                    pair.chrom1 = ids.intern(record.chrom1);
                    pair.chrom2 = ids.intern(record.chrom2);
                }
                if (!args.skip_duplication_check &&
                    duplicates.duplicate(pair, chunk.number[i], args.pairs_file)) {
                    ++total.duplicated_pairs;
                    continue;
                }
                accepted.push_back(static_cast<std::uint32_t>(i));
            }

            // Binning, in contiguous ranges combined in index order.
            std::vector<PairsPartition> partitions(
                std::max<std::size_t>(1, std::min(workers, accepted.size())));
            for_contiguous_ranges(workers, accepted.size(),
                                  [&](std::size_t first, std::size_t last, std::size_t part) {
                                      bin_pairs(bin_index, parsed, accepted, first, last,
                                                args.do_test_run, partitions[part]);
                                  });
            for (auto& partition : partitions) {
                total.add(partition.counters);
                pixel_runs.add(partition.pixels);
                partition.pixels.clear();
                partition.pixels.shrink_to_fit();
            }

            if (args.do_test_run && total.iter_num >= args.do_test_run_lines) {
                break;
            }
        }
        parsed.clear();
        parsed.shrink_to_fit();
        chunk.clear();

        // --- QC --------------------------------------------------------------
        hicx::QcLogInputs qc_inputs;
        qc_inputs.out_file_name = args.out_file_name;
        qc_inputs.min_distance = fragment_bins ? args.min_distance : 0;
        qc_inputs.max_library_insert_size = max_library_insert_size;
        qc_inputs.keep_self_ligation = false;
        qc_inputs.has_restriction_cut_file = false;
        const std::string qc_log = hicx::format_qc_log(qc_inputs, total);
        {
            std::ofstream log_file(args.qc_folder + "/QC.log", std::ios::binary);
            if (!log_file) {
                throw std::runtime_error("could not write " + args.qc_folder + "/QC.log");
            }
            log_file << qc_log;
        }
        hicx::write_qc_tables(args.qc_folder, qc_log);

        if (args.do_test_run) {
            hicx::report_resource_usage(g_tool);
            return 0;
        }

        // --- the matrix --------------------------------------------------------
        const std::size_t pixel_count = pixel_runs.finish();
        hicx::enlarge_bins(bin_intervals, chrom_sizes, names);
        const std::vector<double> bin_max(bin_intervals.size(),
                                          std::numeric_limits<double>::quiet_NaN());
        const std::string& genome_assembly =
            args.genome_assembly.empty() ? header.genome_assembly : args.genome_assembly;
        write_matrix_outputs(
            args, bin_intervals, bin_max, names, pixel_count,
            [&pixel_runs](auto&& sink) { pixel_runs.drain(sink); }, qc_log, genome_assembly);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s: %s\n", g_tool, error.what());
        return 1;
    }
    hicx::report_resource_usage(g_tool);
    return 0;
}

}  // namespace

#endif  // HICX_TOOLS_BUILD_MATRIX_PAIRS_IMPL_HPP
