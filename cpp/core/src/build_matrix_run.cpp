#include "hicx/build_matrix_run.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <deque>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <thread>
#include <unordered_set>

#include "hicx/build_matrix_api.hpp"

namespace hicx::build_matrix {

namespace {

// --- the read pair buffer ---------------------------------------------------

struct PairBuffer {
    std::vector<ReadFields> mate1;
    std::vector<ReadFields> mate2;
    // Kept only for --outBam.
    std::vector<BamRecord> record1;
    std::vector<BamRecord> record2;
    // Pairs read and not accepted, readBamFiles' `iter_num - len(buffer)`.
    std::int64_t discarded = 0;

    void clear() {
        mate1.clear();
        mate2.clear();
        record1.clear();
        record2.clear();
        discarded = 0;
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
        std::uint64_t b = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(key.pos_low)) << 32) |
                          static_cast<std::uint32_t>(key.pos_high);
        a ^= b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2);
        return static_cast<std::size_t>(a);
    }
};

// --- classification of buffered htslib records -------------------------------

struct Context {
    const BinSearchIndex* bins = nullptr;
    const RestrictionSiteIndex* rf = nullptr;
    bool rf_present = false;  // `if pRfPositions` truthiness
    std::vector<std::uint32_t> tid_to_chrom;
    std::vector<std::string> pat_forw;
    std::vector<std::string> pat_rev;
    ClassifyOptions options;
    std::vector<std::uint32_t> coverage_begin;
    std::vector<std::uint32_t> coverage_end;
    std::atomic<std::uint32_t>* coverage = nullptr;
    std::size_t coverage_size = 0;
};

bool starts_with_probe(const ReadFields& read, const std::string& pattern) {
    if (pattern.empty()) {
        return true;
    }
    if (static_cast<int>(pattern.size()) > read.seq_len ||
        static_cast<int>(pattern.size()) > kDanglingProbe) {
        return false;
    }
    return std::memcmp(read.head, pattern.data(), pattern.size()) == 0;
}

bool ends_with_probe(const ReadFields& read, const std::string& pattern) {
    if (pattern.empty()) {
        return true;
    }
    const int n = static_cast<int>(pattern.size());
    if (n > read.seq_len || n > kDanglingProbe) {
        return false;
    }
    const int stored = std::min(read.seq_len, kDanglingProbe);
    return std::memcmp(read.tail + (stored - n), pattern.data(), pattern.size()) == 0;
}

bool check_dangling_end(const ReadFields& read, const std::string& forward, const std::string& reverse) {
    if (!read.is_reverse() && starts_with_probe(read, forward)) {
        return true;
    }
    if (read.is_reverse() && ends_with_probe(read, reverse)) {
        return true;
    }
    return false;
}

struct Partition {
    QcCounters counters;
    std::vector<std::uint64_t> pixels;      // packed (min_bin << 32) | max_bin
    std::vector<std::int64_t> bam_indices;  // indices into the pair buffer
};

// The classify_pairs policy over a PairBuffer of htslib fields.
class RecordPolicy {
  public:
    RecordPolicy(const Context& context, const PairBuffer& buffer, bool want_bam, Partition& out)
        : context_(context), buffer_(buffer), want_bam_(want_bam), out_(out) {}

    [[nodiscard]] const ReadFields& mate(std::size_t i, int w) const {
        return w == 0 ? buffer_.mate1[i] : buffer_.mate2[i];
    }
    [[nodiscard]] std::int64_t tid(std::size_t i, int w) const { return mate(i, w).tid; }
    [[nodiscard]] std::int64_t pos(std::size_t i, int w) const { return mate(i, w).pos; }
    [[nodiscard]] std::int64_t qlen(std::size_t i, int w) const { return mate(i, w).qlen; }
    [[nodiscard]] std::int64_t seq_len(std::size_t i, int w) const { return mate(i, w).seq_len; }
    [[nodiscard]] bool is_reverse(std::size_t i, int w) const { return mate(i, w).is_reverse(); }
    [[nodiscard]] std::optional<SearchInterval> bin(std::size_t i, int w, std::int64_t read_middle) const {
        const std::uint32_t chrom = context_.tid_to_chrom[static_cast<std::size_t>(mate(i, w).tid)];
        const auto flat = context_.bins->flat_index_at(chrom, read_middle);
        if (!flat.has_value()) {
            return std::nullopt;
        }
        return context_.bins->interval(*flat);
    }
    [[nodiscard]] bool rf_enabled() const {
        return context_.rf_present && !context_.options.restriction_sequences.empty();
    }
    [[nodiscard]] bool rf_overlaps(std::size_t i, std::int64_t begin, std::int64_t end) const {
        const std::uint32_t chrom = context_.tid_to_chrom[static_cast<std::size_t>(mate(i, 0).tid)];
        return context_.rf->has_chrom(chrom) && context_.rf->overlaps(chrom, begin, end);
    }
    [[nodiscard]] bool dangling(std::size_t i, int w, std::size_t s) const {
        return check_dangling_end(mate(i, w), context_.pat_forw[s], context_.pat_rev[s]);
    }
    [[nodiscard]] std::pair<std::int64_t, std::int64_t> coverage_range(std::uint32_t bin_id) const {
        return {context_.coverage_begin[bin_id], context_.coverage_end[bin_id]};
    }
    void coverage_add(std::int64_t cell) const {
        if (cell >= 0 && cell < static_cast<std::int64_t>(context_.coverage_size)) {
            context_.coverage[cell].fetch_add(1, std::memory_order_relaxed);
        }
    }
    // Q5: the isize patch of buildMatrixMethods.py:814 happens in the worker's
    // copy and never reaches the written record, so it is not applied.
    void set_isize(std::size_t) const {}
    void add_pair(std::size_t i, std::uint32_t bin1, std::uint32_t bin2) {
        if (!context_.options.quick_qc_mode) {
            const std::uint32_t low = std::min(bin1, bin2);
            const std::uint32_t high = std::max(bin1, bin2);
            out_.pixels.push_back((static_cast<std::uint64_t>(low) << 32) | high);
        }
        if (want_bam_) {
            out_.bam_indices.push_back(static_cast<std::int64_t>(i));
        }
    }

  private:
    const Context& context_;
    const PairBuffer& buffer_;
    bool want_bam_;
    Partition& out_;
};

void classify_range(const Context& context, const PairBuffer& buffer, std::size_t first, std::size_t last,
                    bool want_bam, Partition& out) {
    out.counters.dangling_end.assign(context.options.restriction_sequences.size(), 0);
    RecordPolicy policy(context, buffer, want_bam, out);
    classify_pairs(context.options, first, last, policy, out.counters);
}

// get_supplementary_alignment plus get_correct_map. No record in the test
// corpus carries an SA tag, so this path has never been exercised against the
// Python on real data; it is written from the Python.
void resolve_supplementary(BamReader& reader, BamRecord& primary) {
    if (!primary.has_tag("SA")) {
        return;
    }
    const std::string tag = primary.tag_string("SA");
    const std::size_t entries = static_cast<std::size_t>(std::count(tag.begin(), tag.end(), ';'));
    std::int64_t best = primary.leading_bases_before_match();
    BamRecord best_record = primary;
    BamRecord candidate;
    for (std::size_t i = 0; i < entries; ++i) {
        if (!reader.read(candidate)) {
            break;
        }
        if (candidate.qname() != primary.qname()) {
            continue;
        }
        const std::int64_t value = candidate.leading_bases_before_match();
        if (value < best) {
            best = value;
            best_record = candidate;
        }
    }
    primary = best_record;
}

std::string python_format_fixed(double value, int digits) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.*f", digits, value);
    return buffer;
}

// repr() of a dict_keys of str
std::string dict_keys_repr(const std::vector<std::string>& keys) {
    std::string out = "dict_keys([";
    for (std::size_t i = 0; i < keys.size(); ++i) {
        out += (i > 0 ? ", " : "") + python_str_repr(keys[i]);
    }
    return out + "])";
}

constexpr const char* kUnboundOutBam =
    "cannot access local variable 'out_bam_file' where it is not associated with a value";

}  // namespace

std::string to_upper(std::string text) {
    for (char& c : text) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return text;
}

std::string python_str_repr(const std::string& text) {
    const bool has_single = text.find('\'') != std::string::npos;
    const bool has_double = text.find('"') != std::string::npos;
    const char quote = has_single && !has_double ? '"' : '\'';
    std::string out(1, quote);
    for (const char c : text) {
        const auto u = static_cast<unsigned char>(c);
        if (c == quote || c == '\\') {
            out += '\\';
            out += c;
        } else if (c == '\t') {
            out += "\\t";
        } else if (c == '\n') {
            out += "\\n";
        } else if (c == '\r') {
            out += "\\r";
        } else if (u < 0x20 || u == 0x7f) {
            char buffer[8];
            std::snprintf(buffer, sizeof(buffer), "\\x%02x", u);
            out += buffer;
        } else {
            out += c;
        }
    }
    out += quote;
    return out;
}

ChromSizes read_chromosome_sizes(const std::string& path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("could not open " + path);
    }
    // The Python builds an OrderedDict, so a repeated name keeps its first
    // position and its last size.
    ChromSizes sizes;
    std::string line;
    while (std::getline(file, line)) {
        while (!line.empty() &&
               (line.back() == '\r' || line.back() == '\n' || line.back() == ' ' || line.back() == '\t')) {
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
        const auto found = std::find_if(sizes.begin(), sizes.end(),
                                        [&name](const auto& entry) { return entry.first == name; });
        if (found != sizes.end()) {
            found->second = size;
        } else {
            sizes.emplace_back(name, size);
        }
    }
    return sizes;
}

void PixelAccumulator::absorb(std::vector<std::uint64_t>& chunk) {
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

void PixelAccumulator::merge(const std::vector<std::uint64_t>& keys, const std::vector<std::int64_t>& counts) {
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

std::map<std::string, std::string> matrix_metadata(const std::string& qc_log, const std::string& genome_assembly,
                                                   const std::string& qc_json) {
    // Q9 buildMatrixMethods.py:1371-1376. The three provenance strings are
    // wrapped in np.string_, that is numpy bytes, and hicmatrix's cool
    // writer stores str(value) of them (hicmatrix/lib/cool.py:394, 398,
    // 401). str() of a bytes object is its repr, so what lands in the
    // cool file is the seven extra characters of "b'...'" around the
    // value. Reproduced, not fixed: the attribute is part of the file.
    const auto as_python_bytes_repr = [](const std::string& value) { return "b'" + value + "'"; };
    std::map<std::string, std::string> metadata;
    metadata["statistics"] = qc_log;
    metadata["matrix-generated-by"] = as_python_bytes_repr(std::string("HiCExplorer-") + kVersion);
    metadata["matrix-generated-by-url"] = as_python_bytes_repr("https://github.com/deeptools/HiCExplorer");
    if (!genome_assembly.empty()) {
        metadata["genome-assembly"] = as_python_bytes_repr(genome_assembly);
    }
    if (!qc_json.empty()) {
        // Raw JSON, stored as a nested object under "hicexplorer_qc" in the
        // cool metadata (save_matrix hands it to CoolSaveOptions::metadata_json).
        metadata["hicexplorer_qc"] = qc_json;
    }
    return metadata;
}

namespace {
// The QC statistics travel in the metadata map under "hicexplorer_qc"; the
// cool writer takes them as a nested JSON object of its own field, so that
// they are written for every resolution of an mcool (the provenance keys are
// erased from the shared map after the first resolution, this one is not).
void take_qc_json(CoolSaveOptions& options) {
    const auto found = options.hic_metadata.find("hicexplorer_qc");
    if (found != options.hic_metadata.end()) {
        options.metadata_json.emplace_back("hicexplorer_qc", found->second);
        options.hic_metadata.erase(found);
    }
}
}  // namespace

void save_matrix(const std::string& out_file_name, const std::vector<std::int64_t>& bin_size,
                 MatrixData& data, const std::map<std::string, std::string>& metadata) {
    std::error_code ec;
    std::filesystem::remove(out_file_name, ec);

    if (ends_with(out_file_name, ".mcool") && bin_size.size() > 2) {
        // Q7: only three or more resolutions take this branch.
        CoolSaveOptions options;
        options.symmetric = true;
        options.apply_correction = false;
        options.hic_metadata = metadata;
        options.has_hic_metadata = true;
        take_qc_json(options);
        write_cool(out_file_name + "::/resolutions/" + std::to_string(bin_size[0]), data, options);
        // Q10 hicmatrix/lib/cool.py:394-402. create_cooler_input deletes
        // 'matrix-generated-by', 'matrix-generated-by-url' and
        // 'genome-assembly' out of the metadata dictionary it was handed,
        // and buildMatrixMethods.py:1381-1404 hands the *same* dictionary
        // to every resolution. So only the first resolution of an mcool
        // carries the provenance and the rest carry none. Reproduced.
        options.hic_metadata.erase("matrix-generated-by");
        options.hic_metadata.erase("matrix-generated-by-url");
        options.hic_metadata.erase("genome-assembly");
        for (std::size_t r = 1; r < bin_size.size(); ++r) {
            const std::int64_t factor = bin_size[r] / bin_size[0];
            MatrixData merged = merge_bins(data, factor);
            CoolSaveOptions append = options;
            append.append = true;
            write_cool(out_file_name + "::/resolutions/" + std::to_string(bin_size[r]), merged, append);
        }
    } else if (ends_with(out_file_name, "h5")) {
        write_hicexplorer_h5(out_file_name, data);
    } else if (ends_with(out_file_name, "cool")) {
        CoolSaveOptions options;
        options.symmetric = true;
        options.apply_correction = false;
        options.hic_metadata = metadata;
        options.has_hic_metadata = true;
        take_qc_json(options);
        write_cool(out_file_name, data, options);
    }
    // A name ending in neither is written nowhere at all, silently, which
    // is hiCMatrix.save's behaviour.
}

void run(const Arguments& args, const Hooks& hooks) {
    const bool api = args.python_api;
    const auto log = [&hooks](const char* level, const std::string& message) {
        if (hooks.log) {
            hooks.log(level, message);
        }
    };

    if (!ends_with(args.out_file_name, ".h5") && !ends_with(args.out_file_name, ".cool") &&
        args.out_file_name.find(".mcool") == std::string::npos) {
        throw Exit("Please define the file extension. h5 and cool are supported, or the "
                   "specializations of cool, mcool. Given input " + args.out_file_name,
                   true);
    }
    std::int64_t max_library_insert_size = args.max_library_insert_size;
    if (args.max_distance.has_value()) {
        max_library_insert_size = *args.max_distance;  // backwards compatibility
    }
    {
        std::error_code ec;
        std::filesystem::create_directories(args.qc_folder, ec);
        if (!std::filesystem::is_directory(args.qc_folder)) {
            throw Exit("Can't open/create QC folder path: " + args.qc_folder + ". Please check", false);
        }
    }
    std::int64_t threads = args.threads;
    if (threads < 2) {
        threads = 2;
        const std::string message = "\nAt least two threads need to be defined. Setting --threads = 2!s\n";
        if (hooks.warn) {
            hooks.warn(message);
        } else {
            std::fputs(message.c_str(), stderr);
        }
    }
    if (api && args.dangling_given && args.restriction_sequences.empty()) {
        throw Exit("\nIf --danglingSequence is set, --restrictionSequence needs to be set too.\n", false);
    }
    if (!api && args.restriction_sequences.size() != args.dangling_sequences.size()) {
        // buildMatrixMethods.py:968-970 indexes pDanglingSequence with the
        // restriction sequence index, so a short list is an IndexError; the
        // tools refuse any mismatch.
        throw std::runtime_error("--danglingSequence has " + std::to_string(args.dangling_sequences.size()) +
                                 " entries but --restrictionSequence has " +
                                 std::to_string(args.restriction_sequences.size()) +
                                 "; they are matched by position");
    }

    const std::string reading = "reading " + args.sam_files[0] + " and " + args.sam_files[1] + " to build hic_matrix\n";
    if (hooks.log) {
        hooks.log("info", reading);
    }
    BamReader reader1(args.sam_files[0]);
    BamReader reader2(args.sam_files[1]);
    if (!hooks.log) {
        std::fprintf(stderr, "INFO:hicexplorer.lib.buildMatrixMethods:%s\n", reading.c_str());
    }
    if (hooks.bams_opened) {
        hooks.bams_opened();
    }

    std::unique_ptr<BamWriter> out_bam;
    if (!args.do_test_run && !args.out_bam.empty()) {
        out_bam = std::make_unique<BamWriter>(args.out_bam, reader1);
    }

    ChromSizes chrom_sizes;
    if (hooks.chromosome_sizes) {
        chrom_sizes = hooks.chromosome_sizes();
    } else if (args.chromosome_sizes.empty()) {
        chrom_sizes = chrom_sizes_of(reader1);
    } else {
        chrom_sizes = read_chromosome_sizes(args.chromosome_sizes);
    }

    ChromNames names;
    // Intern the reference names first so that tid_to_chrom is a direct
    // lookup and the id order matches the BAM header.
    for (const auto& reference : reader1.references()) {
        names.intern(reference);
    }
    for (const auto& entry : chrom_sizes) {
        names.intern(entry.first);
    }

    std::vector<GenomeInterval> rf_interval;
    if (hooks.restriction_sites) {
        std::vector<std::string> keys;
        for (const auto& [chrom, start, end] : hooks.restriction_sites(chrom_sizes)) {
            const std::uint32_t id = names.intern(chrom);
            if (std::find(keys.begin(), keys.end(), chrom) == keys.end()) {
                keys.push_back(chrom);
            }
            rf_interval.push_back({id, start, end});
        }
        if (args.restriction_cut_file_given) {
            log("debug", "rf_positions " + dict_keys_repr(keys));
        }
    } else {
        for (const auto& path : args.restriction_cut_files) {
            bed2interval_list(path, chrom_sizes, args.region, names, rf_interval);
        }
    }
    std::vector<GenomeInterval> bin_intervals;
    std::vector<bool> start_is_python_zero;
    if (!args.bin_size.empty()) {
        if (hooks.user_region) {
            const auto [region_sizes, start] = hooks.user_region(chrom_sizes);
            bin_intervals = get_bins(args.bin_size[0], region_sizes, start, names);
        } else {
            bin_intervals = get_bins(args.bin_size[0], chrom_sizes, args.region, names);
        }
    } else {
        if (api) {
            log("info", "Minimum distance considered between restriction sites is " +
                            (args.min_distance_none ? std::string("None") : std::to_string(args.min_distance)) +
                            "\nMax distance: " + std::to_string(max_library_insert_size) + "\n");
            if (rf_interval.empty()) {
                throw PythonError("ValueError", "not enough values to unpack (expected 3, got 0)");
            }
            if (args.min_distance_none) {
                // np.diff(start) - rest_site_len <= None compares element-wise;
                // with a single site the bin filter fails instead
                throw PythonError("TypeError", rf_interval.size() > 1
                                                   ? "'<=' not supported between instances of 'int' and 'NoneType'"
                                                   : "'>=' not supported between instances of 'int' and 'NoneType'");
            }
        }
        bin_intervals = get_rf_bins(rf_interval, args.min_distance, max_library_insert_size);
    }
    const BinSearchIndex bin_index(bin_intervals, names);
    const RestrictionSiteIndex rf_index(rf_interval, names.size());
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

    // The restriction sequences, their dangling sequences and the keys of the
    // Python's dangling_end dictionary (a repeated restriction sequence is one
    // key there, with the dangling sequence given last).
    std::vector<std::string> dangling_keys;
    std::vector<std::size_t> key_of;  // restriction sequence index -> key index
    std::vector<std::string> key_dangling;
    if (api) {
        const bool build = args.dangling_given;
        for (std::size_t i = 0; i < args.restriction_sequences.size(); ++i) {
            if (build && hooks.uppercase) {
                hooks.uppercase(i, i < args.dangling_sequences.size());
            }
            if (build && i >= args.dangling_sequences.size()) {
                throw PythonError("IndexError", "list index out of range");
            }
            const std::string restriction =
                build ? to_upper(args.restriction_sequences[i]) : args.restriction_sequences[i];
            const std::string dangling = build ? to_upper(args.dangling_sequences[i]) : std::string();
            const auto found = std::find(dangling_keys.begin(), dangling_keys.end(), restriction);
            if (found == dangling_keys.end()) {
                key_of.push_back(dangling_keys.size());
                dangling_keys.push_back(restriction);
                key_dangling.push_back(dangling);
            } else {
                key_of.push_back(static_cast<std::size_t>(found - dangling_keys.begin()));
                key_dangling[key_of.back()] = dangling;
            }
            context.options.restriction_sequences.push_back(restriction);
        }
        for (std::size_t i = 0; i < context.options.restriction_sequences.size(); ++i) {
            const std::string& dangling = key_dangling[key_of[i]];
            context.pat_forw.push_back(dangling);
            context.pat_rev.push_back(reverse_complement(dangling));
        }
        if (build) {
            std::string text = "{";
            for (std::size_t k = 0; k < dangling_keys.size(); ++k) {
                text += (k > 0 ? ", " : "") + python_str_repr(dangling_keys[k]) + ": {'pat_forw': " +
                        python_str_repr(key_dangling[k]) +
                        ", 'pat_rev': " + python_str_repr(reverse_complement(key_dangling[k])) + "}";
            }
            log("info", "dangling sequences to check are " + text + "}\n");
        }
        context.options.has_dangling = build && !dangling_keys.empty();
    } else {
        for (std::size_t i = 0; i < args.restriction_sequences.size(); ++i) {
            const std::string restriction = to_upper(args.restriction_sequences[i]);
            const std::string dangling = to_upper(args.dangling_sequences[i]);
            context.options.restriction_sequences.push_back(restriction);
            context.pat_forw.push_back(dangling);
            context.pat_rev.push_back(reverse_complement(dangling));
        }
        context.options.has_dangling = !args.dangling_sequences.empty();
    }
    context.options.keep_self_circles = args.keep_self_circles;
    context.options.keep_self_ligation = args.keep_self_ligation;
    context.options.max_insert_size = max_library_insert_size;
    context.options.quick_qc_mode = args.do_test_run;

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
    const std::size_t worker_count = static_cast<std::size_t>(std::max<std::int64_t>(1, threads - 1));

    QcCounters total;
    total.dangling_end.assign(context.options.restriction_sequences.size(), 0);
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

    BamRecord record1;
    BamRecord record2;
    bool all_data_read = false;
    // --outBam with --doTestRun: the tools ignore it (D2), the Python fails
    // when it writes the first pair.
    const bool want_bam = out_bam != nullptr || (api && args.do_test_run && !args.out_bam.empty());

    // readBamFiles: fills `buffer` with up to input_buffer_size accepted pairs.
    const auto read_buffer = [&](PairBuffer& buffer) {
        buffer.clear();
        std::int64_t read = 0;
        while (static_cast<std::int64_t>(buffer.size()) < input_buffer_size) {
            if (!reader1.read(record1) || !reader2.read(record2)) {
                all_data_read = true;
                break;
            }
            ++read;
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
                const std::string message = "FATAL ERROR " + record1.qname() + " " + record2.qname() +
                                            " Be sure that the sam files have the same read order If "
                                            "using Bowtie2 or Hisat2 add the --reorder option";
                if (api) {
                    throw PythonError("AssertionError", message);
                }
                throw std::runtime_error(message);
            }
            resolve_supplementary(reader1, record1);
            resolve_supplementary(reader2, record2);

            switch (pair_verdict(record1.flag(), record2.flag(), record1.mapq(), record2.mapq(),
                                 args.min_mapping_quality)) {
                case PairVerdict::OneMateUnmapped: ++total.one_mate_unmapped; continue;
                case PairVerdict::OneMateNotUnique: ++total.one_mate_not_unique; continue;
                case PairVerdict::LowQuality: ++total.one_mate_low_quality; continue;
                case PairVerdict::Keep: break;
            }
            if (!args.skip_duplication_check) {
                const std::uint32_t rank1 =
                    record1.tid() >= 0 ? name_rank[static_cast<std::size_t>(record1.tid())] : 0;
                const std::uint32_t rank2 =
                    record2.tid() >= 0 ? name_rank[static_cast<std::size_t>(record2.tid())] : 0;
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
        buffer.discarded = read - static_cast<std::int64_t>(buffer.size());
    };

    // process_data on the threads, then the master's part of the result.
    const auto process_buffer = [&](PairBuffer& buffer) {
        // Fixed contiguous partitions, reduced sequentially inside a
        // partition and combined in index order.
        const std::size_t pairs = buffer.size();
        const std::size_t parts = std::max<std::size_t>(1, std::min<std::size_t>(worker_count, pairs));
        std::vector<Partition> partitions(parts);
        std::vector<std::thread> workers;
        workers.reserve(parts);
        for (std::size_t p = 0; p < parts; ++p) {
            const std::size_t first = pairs * p / parts;
            const std::size_t last = pairs * (p + 1) / parts;
            if (parts == 1) {
                classify_range(context, buffer, first, last, want_bam, partitions[p]);
            } else {
                workers.emplace_back([&, p, first, last] {
                    classify_range(context, buffer, first, last, want_bam, partitions[p]);
                });
            }
        }
        for (auto& worker : workers) {
            worker.join();
        }

        std::vector<std::uint64_t> chunk_pixels;
        for (auto& partition : partitions) {
            total.add(partition.counters);
            chunk_pixels.insert(chunk_pixels.end(), partition.pixels.begin(), partition.pixels.end());
            partition.pixels.clear();
            partition.pixels.shrink_to_fit();
        }
        accumulator.absorb(chunk_pixels);

        if (want_bam) {
            for (const auto& partition : partitions) {
                for (const std::int64_t index : partition.bam_indices) {
                    if (out_bam == nullptr) {
                        throw PythonError("UnboundLocalError", kUnboundOutBam);
                    }
                    BamRecord& mate1 = buffer.record1[static_cast<std::size_t>(index)];
                    BamRecord& mate2 = buffer.record2[static_cast<std::size_t>(index)];
                    mate1.set_paired_flag_first();
                    mate2.set_paired_flag_second();
                    mate1.set_mate(mate2.tid(), mate2.pos());
                    mate2.set_mate(mate1.tid(), mate1.pos());
                    out_bam->write(mate1);
                    out_bam->write(mate2);
                }
            }
        }
    };

    if (!api) {
        PairBuffer buffer;
        while (!all_data_read) {
            read_buffer(buffer);
            total.iter_num += buffer.discarded + static_cast<std::int64_t>(buffer.size());
            if (buffer.size() == 0) {
                break;
            }
            process_buffer(buffer);
            if (args.do_test_run && total.iter_num > args.do_test_run_lines) {
                break;
            }
        }
    } else {
        // createMatrix's master loop: every worker holds a buffer that was
        // read before any result was taken, results are taken in the order
        // the buffers were read, and a buffer is read for a worker as soon as
        // its previous result was taken. readBamFiles' discarded pairs count
        // when a buffer is read, its pairs when its result is taken.
        const auto start_time = std::chrono::system_clock::now();
        std::deque<PairBuffer> in_flight;
        bool stop_reading = false;
        const auto refill = [&] {
            while (!stop_reading && !all_data_read && in_flight.size() < worker_count) {
                PairBuffer buffer;
                read_buffer(buffer);
                total.iter_num += buffer.discarded;
                if (buffer.size() == 0) {
                    // D1: readBamFiles returned (None, None, True, ...)
                    throw PythonError("TypeError", "object of type 'NoneType' has no len()");
                }
                in_flight.push_back(std::move(buffer));
            }
        };
        refill();
        while (!in_flight.empty()) {
            PairBuffer buffer = std::move(in_flight.front());
            in_flight.pop_front();
            process_buffer(buffer);
            total.iter_num += static_cast<std::int64_t>(buffer.size());
            buffer.clear();
            if (std::fmod(static_cast<double>(total.iter_num), 1e6) < 100000) {
                const double elapsed =
                    std::chrono::duration<double>(std::chrono::system_clock::now() - start_time).count();
                log("info", "processing " + std::to_string(total.iter_num) + " lines took " +
                                python_format_fixed(elapsed, 2) + " secs (" +
                                python_format_fixed(static_cast<double>(total.iter_num) / elapsed, 1) +
                                " lines per second)\n");
                log("info", std::to_string(total.pair_added) + " (" +
                                python_format_fixed(static_cast<double>(100 * total.pair_added) /
                                                        static_cast<double>(total.iter_num),
                                                    2) +
                                "%) valid pairs added to matrix\n");
            }
            if (args.do_test_run && total.iter_num > args.do_test_run_lines) {
                log("debug", "\n## *WARNING*. Early exit because of --doTestRun parameter  ##\n\n");
                stop_reading = true;
            }
            refill();
        }
        log("debug", "Parallel stuff done");
    }
    seen.clear();

    if (out_bam) {
        out_bam->close();
    }

    // The Python enlarges the bins before it writes QC.log.
    std::vector<GenomeInterval> enlarged;
    if (api && !args.do_test_run) {
        // dict(chrom_sizes)[chrom], in the order enlarge_bins looks them up
        std::vector<char> known(names.size(), 0);
        for (const auto& entry : chrom_sizes) {
            known[*names.lookup(entry.first)] = 1;
        }
        const auto require = [&](std::uint32_t chrom) {
            if (known[chrom] == 0) {
                throw PythonError("KeyError", python_str_repr(names.name(chrom)));
            }
        };
        for (std::size_t i = 0; i + 1 < bin_intervals.size(); ++i) {
            if (bin_intervals[i].chrom != bin_intervals[i + 1].chrom) {
                require(bin_intervals[i].chrom);
            }
        }
        if (bin_intervals.empty()) {
            throw PythonError("IndexError", "list index out of range");
        }
        require(bin_intervals.back().chrom);
        enlarged = bin_intervals;
        enlarge_bins(enlarged, chrom_sizes, names);
    }

    // --- QC ------------------------------------------------------------
    QcLogInputs qc_inputs;
    qc_inputs.out_file_name = args.out_file_name;
    qc_inputs.min_distance = args.min_distance;
    qc_inputs.max_library_insert_size = max_library_insert_size;
    qc_inputs.keep_self_ligation = args.keep_self_ligation;
    std::string qc_log;
    std::string qc_json;
    for (const auto& entry : chrom_sizes) {
        qc_inputs.max_chromosome_length =
            std::max<std::int64_t>(qc_inputs.max_chromosome_length, static_cast<std::int64_t>(entry.second));
    }
    if (api) {
        qc_inputs.has_restriction_cut_file = args.restriction_cut_file_given;
        QcCounters merged = total;
        merged.dangling_end.assign(dangling_keys.size(), 0);
        for (std::size_t i = 0; i < total.dangling_end.size(); ++i) {
            merged.dangling_end[key_of[i]] += total.dangling_end[i];
        }
        if (!dangling_keys.empty()) {
            std::vector<std::string> known;
            if (args.dangling_given) {
                known = dangling_keys;
            }
            log("debug", "dangling_sequences " + dict_keys_repr(known));
            log("debug", "dangling_end " + dict_keys_repr(dangling_keys));
            if (!args.dangling_given) {
                throw PythonError("KeyError", python_str_repr(dangling_keys.front()));
            }
            for (std::size_t k = 0; k < dangling_keys.size(); ++k) {
                qc_inputs.dangling_sequences.emplace_back(dangling_keys[k], key_dangling[k]);
            }
        }
        qc_log = format_qc_log(qc_inputs, merged);
        qc_json = format_qc_metadata_json(qc_inputs, merged, "hicBuildMatrix");
    } else {
        qc_inputs.has_restriction_cut_file = !args.restriction_cut_files.empty();
        for (std::size_t i = 0; i < context.options.restriction_sequences.size(); ++i) {
            qc_inputs.dangling_sequences.emplace_back(context.options.restriction_sequences[i],
                                                      context.pat_forw[i]);
        }
        qc_log = format_qc_log(qc_inputs, total);
        qc_json = format_qc_metadata_json(qc_inputs, total, args.tool_name.empty() ? "hicBuildMatrix" : args.tool_name);
    }
    const std::string qc_log_file =
        api ? build_matrix_api::join(args.qc_folder, "QC.log") : args.qc_folder + "/QC.log";
    {
        std::ofstream log_file(qc_log_file, std::ios::binary);
        if (!log_file) {
            throw std::runtime_error("could not write " + qc_log_file);
        }
        log_file << qc_log;
    }
    if (api) {
        log("debug", "log_file_name " + qc_log_file);
        log("debug", "pQCfolder " + args.qc_folder);
    }
    if (hooks.after_qc_log) {
        hooks.after_qc_log(qc_log_file);
    } else {
        write_qc_tables(args.qc_folder, qc_log);
    }
    if (hooks.before_matrix) {
        hooks.before_matrix();
    }

    if (args.do_test_run) {
        if (api && ends_with(args.out_file_name, ".mcool") && args.bin_size.size() > 2) {
            throw PythonError("UnboundLocalError",
                              "cannot access local variable 'hic_ma' where it is not associated with a value");
        }
        return;
    }

    // --- the matrix ------------------------------------------------------
    if (!api) {
        enlarged = bin_intervals;
        enlarge_bins(enlarged, chrom_sizes, names);
    }
    std::vector<double> bin_max(enlarged.size());
    for (std::size_t i = 0; i < enlarged.size(); ++i) {
        std::uint32_t max_element = 0;
        const std::int64_t from = context.coverage_begin[i];
        const std::int64_t to = context.coverage_end[i];  // Q4: exclusive end
        for (std::int64_t k = from; k < to && k < static_cast<std::int64_t>(coverage.size()); ++k) {
            const std::uint32_t value = coverage[static_cast<std::size_t>(k)].load(std::memory_order_relaxed);
            if (value > max_element) {
                max_element = value;
            }
        }
        bin_max[i] = max_element == 0 ? std::numeric_limits<double>::quiet_NaN() : static_cast<double>(max_element);
    }
    coverage.clear();

    write_matrix_outputs(
        args, enlarged, bin_max, names, accumulator.keys().size(),
        [&accumulator](auto&& sink) {
            const auto& keys = accumulator.keys();
            const auto& counts = accumulator.counts();
            for (std::size_t k = 0; k < keys.size(); ++k) {
                sink(keys[k], static_cast<std::uint64_t>(counts[k]));
            }
        },
        qc_log, args.genome_assembly, qc_json);
}

}  // namespace hicx::build_matrix
