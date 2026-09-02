#include "hicx/build_matrix.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <limits>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace hicx {
namespace {

// Python's repr(float), which is what pandas.DataFrame.to_csv emits when no
// float_format is given: the shortest decimal string that round trips, with a
// mandatory ".0" on an integral value.
std::string python_float_repr(double value) {
    if (std::isnan(value)) {
        return {};  // pandas na_rep, an empty field
    }
    if (std::isinf(value)) {
        return value > 0 ? "inf" : "-inf";
    }
    char buffer[64];
    const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value);
    std::string text(buffer, result.ptr);
    if (text.find('.') == std::string::npos && text.find('e') == std::string::npos &&
        text.find("inf") == std::string::npos && text.find("nan") == std::string::npos) {
        text += ".0";
    }
    return text;
}

std::string format_two_decimals(double value) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.2f", value);
    return buffer;
}

std::string trim(const std::string& text) {
    const auto begin = text.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) {
        return {};
    }
    const auto end = text.find_last_not_of(" \t\r\n");
    return text.substr(begin, end - begin + 1);
}

std::vector<std::string> split(const std::string& text, char delimiter) {
    std::vector<std::string> parts;
    std::size_t start = 0;
    while (true) {
        const auto position = text.find(delimiter, start);
        if (position == std::string::npos) {
            parts.push_back(text.substr(start));
            break;
        }
        parts.push_back(text.substr(start, position - start));
        start = position + 1;
    }
    return parts;
}

}  // namespace

std::uint32_t ChromNames::intern(const std::string& name) {
    const auto found = index_.find(name);
    if (found != index_.end()) {
        return found->second;
    }
    const auto id = static_cast<std::uint32_t>(names_.size());
    names_.push_back(name);
    index_.emplace(name, id);
    return id;
}

std::optional<std::uint32_t> ChromNames::lookup(const std::string& name) const {
    const auto found = index_.find(name);
    if (found == index_.end()) {
        return std::nullopt;
    }
    return found->second;
}

std::string normalise_region(const std::string& text) {
    std::string region;
    for (const char c : text) {
        if (std::isspace(static_cast<unsigned char>(c)) != 0) {
            continue;
        }
        if (std::string(",;|!{}()").find(c) != std::string::npos) {
            continue;
        }
        region.push_back(c == '-' ? ':' : c);
    }
    return region;
}

UserRegion user_region(const ChromSizes& chrom_sizes, const std::string& region) {
    const std::vector<std::string> parts = split(region, ':');
    const std::string& chrom = parts.front();
    std::int64_t size = -1;
    for (const auto& entry : chrom_sizes) {
        if (entry.first == chrom) {
            size = entry.second;  // dict(chromSizes): the last value wins
        }
    }
    if (size < 0) {
        std::string known;
        for (const auto& entry : chrom_sizes) {
            known += known.empty() ? "'" : ", '";
            known += entry.first + "'";
        }
        throw std::runtime_error("Unknown chromosome: " + chrom +
                                 "\nKnown chromosomes are: [" + known + "] ");
    }
    UserRegion out;
    out.start = 0;
    if (parts.size() > 1 && !parts[1].empty()) {
        out.start = std::stoll(parts[1]);
    }
    out.end = size;
    if (parts.size() > 2 && !parts[2].empty()) {
        const std::int64_t value = std::stoll(parts[2]);
        out.end = value <= size ? value : size;
    }
    if (out.start > out.end || out.start < 0) {
        throw std::runtime_error(region +
                                 " not valid. The format is chrom:start:end. "
                                 "Without comas, dashes or dots. ");
    }
    // The tilesize form of getUserRegion is never reachable from
    // hicBuildMatrix: genomicRegion turns every '-' into ':' but the tool's
    // help documents chr:start-end only, and a fourth field would have to be
    // typed with a colon. It is left out rather than guessed at.
    out.chrom_sizes = {{chrom, out.end}};
    return out;
}

std::vector<GenomeInterval> get_bins(std::int64_t bin_size,
                                     const ChromSizes& chrom_sizes,
                                     const std::string& region, ChromNames& names) {
    std::vector<GenomeInterval> bins;
    std::int64_t start = 0;
    const ChromSizes* sizes = &chrom_sizes;
    ChromSizes region_sizes;
    if (!region.empty()) {
        const UserRegion user = user_region(chrom_sizes, region);
        region_sizes = user.chrom_sizes;
        start = user.start;
        sizes = &region_sizes;
    }
    for (const auto& [chrom, size] : *sizes) {
        const std::uint32_t id = names.intern(chrom);
        for (std::int64_t interval = start; interval < size; interval += bin_size) {
            bins.push_back({id, interval, std::min(size, interval + bin_size)});
        }
    }
    return bins;
}

void bed2interval_list(const std::string& path, const ChromSizes& chrom_sizes,
                       const std::string& region, ChromNames& names,
                       std::vector<GenomeInterval>& out) {
    std::string region_chrom;
    std::int64_t region_start = 0;
    std::int64_t region_end = 0;
    const bool filtered = !region.empty();
    if (filtered) {
        const UserRegion user = user_region(chrom_sizes, region);
        region_chrom = user.chrom_sizes.front().first;
        region_start = user.start;
        region_end = user.end;
    }

    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("could not open " + path);
    }
    std::string line;
    std::size_t count = 0;
    while (std::getline(file, line)) {
        ++count;
        std::istringstream stream(line);
        std::string chrom;
        std::string start_text;
        std::string end_text;
        if (!(stream >> chrom >> start_text >> end_text)) {
            // line.strip().split() with fewer than three fields: the Python
            // logs an error and then falls through with the values of the
            // previous line still bound. An empty trailing line is the common
            // case and is skipped here rather than duplicating the previous
            // interval, which is the only observable difference and is
            // recorded in the report.
            std::fprintf(stderr,
                         "ERROR:hicexplorer.lib.buildMatrixMethods:error reading "
                         "BED file at line %zu\n",
                         count);
            continue;
        }
        const std::int64_t start = std::stoll(start_text);
        const std::int64_t end = std::stoll(end_text);
        if (filtered) {
            if (chrom == region_chrom && region_start <= start && region_end <= end) {
                out.push_back({names.intern(chrom), start, end});
            }
        } else {
            out.push_back({names.intern(chrom), start, end});
        }
    }
}

std::vector<GenomeInterval> get_rf_bins(const std::vector<GenomeInterval>& cut_sites,
                                        std::int64_t min_distance,
                                        std::int64_t max_distance) {
    if (cut_sites.empty()) {
        // zip(*[]) raises "not enough values to unpack"; the port reports it.
        throw std::runtime_error(
            "no restriction cut sites left to build bins from. With --region "
            "this is the empty list bed2interval_list produces; see the QUIRK "
            "note on the region filter");
    }
    const std::size_t n = cut_sites.size();
    const std::int64_t rest_site_len = cut_sites[0].end - cut_sites[0].start;

    // np.flatnonzero(np.diff(start) - rest_site_len <= min_distance) + 1
    std::vector<std::size_t> to_merge;
    for (std::size_t i = 1; i < n; ++i) {
        if (cut_sites[i].start - cut_sites[i - 1].start - rest_site_len <=
            min_distance) {
            to_merge.push_back(i);
        }
    }

    std::vector<std::int64_t> start(n);
    std::vector<std::int64_t> end(n);
    for (std::size_t i = 0; i < n; ++i) {
        start[i] = cut_sites[i].start - max_distance;
        end[i] = cut_sites[i].end + max_distance;
    }

    std::vector<std::int64_t> new_start{std::max<std::int64_t>(0, start[0])};
    std::vector<std::int64_t> new_end;
    std::vector<std::uint32_t> new_chrom{cut_sites[0].chrom};
    std::size_t merge_idx = 0;
    for (std::size_t idx = 1; idx < n; ++idx) {
        if (cut_sites[idx].chrom != cut_sites[idx - 1].chrom) {
            new_start.push_back(std::max<std::int64_t>(0, start[idx]));
            new_end.push_back(end[idx - 1]);
            new_chrom.push_back(cut_sites[idx].chrom);
            ++merge_idx;
            continue;
        }
        if (merge_idx < to_merge.size() && idx == to_merge[merge_idx]) {
            ++merge_idx;
            continue;
        }
        if (end[idx - 1] > start[idx]) {
            const std::int64_t middle = start[idx] + (end[idx - 1] - start[idx]) / 2;
            new_start.push_back(middle);
            new_end.push_back(middle);
        } else {
            new_start.push_back(start[idx]);
            new_end.push_back(end[idx - 1]);
        }
        new_chrom.push_back(cut_sites[idx].chrom);
    }
    new_end.push_back(end[n - 1]);
    if (new_chrom.size() != new_start.size() || new_end.size() != new_start.size()) {
        throw std::runtime_error("error");  // the Python assert
    }

    std::vector<GenomeInterval> bins;
    bins.reserve(new_start.size());
    for (std::size_t i = 0; i < new_start.size(); ++i) {
        if (new_end[i] - new_start[i] >= min_distance) {
            bins.push_back({new_chrom[i], new_start[i], new_end[i]});
        }
    }
    return bins;
}

void enlarge_bins(std::vector<GenomeInterval>& bins, const ChromSizes& chrom_sizes,
                  const ChromNames& names) {
    if (bins.empty()) {
        return;
    }
    std::unordered_map<std::string, std::int64_t> sizes;
    for (const auto& entry : chrom_sizes) {
        sizes[entry.first] = entry.second;
    }
    bool chr_start = true;
    for (std::size_t idx = 0; idx + 1 < bins.size(); ++idx) {
        const std::uint32_t chrom = bins[idx].chrom;
        std::int64_t start = bins[idx].start;
        const std::int64_t end = bins[idx].end;
        const std::uint32_t chrom_next = bins[idx + 1].chrom;
        const std::int64_t start_next = bins[idx + 1].start;
        const std::int64_t end_next = bins[idx + 1].end;
        if (chr_start) {
            start = 0;
            chr_start = false;
        }
        if (chrom == chrom_next && end != start_next) {
            const std::int64_t middle = start_next - (start_next - end) / 2;
            bins[idx] = {chrom, start, middle};
            bins[idx + 1] = {chrom, middle, end_next};
        }
        if (chrom != chrom_next) {
            const auto found = sizes.find(names.name(chrom));
            if (found == sizes.end()) {
                throw std::runtime_error("unknown chromosome in enlarge_bins: " +
                                         names.name(chrom));
            }
            bins[idx] = {chrom, start, found->second};
            bins[idx + 1] = {chrom_next, 0, end_next};
        }
    }
    const auto& last = bins.back();
    const auto found = sizes.find(names.name(last.chrom));
    if (found == sizes.end()) {
        throw std::runtime_error("unknown chromosome in enlarge_bins: " +
                                 names.name(last.chrom));
    }
    bins.back() = {last.chrom, last.start, found->second};
}

BinSearchIndex::BinSearchIndex(const std::vector<GenomeInterval>& bins,
                               const ChromNames& names) {
    (void)names;
    // The Python groups by chromosome in the insertion order of
    // intervalListToIntervalTree's dict, that is by first appearance in the
    // bin list, and sorts (begin, end, data) within a group.
    std::vector<std::uint32_t> order;
    std::unordered_map<std::uint32_t, std::vector<std::size_t>> groups;
    for (std::size_t i = 0; i < bins.size(); ++i) {
        auto& group = groups[bins[i].chrom];
        if (group.empty()) {
            order.push_back(bins[i].chrom);
        }
        group.push_back(i);
    }
    flat_.reserve(bins.size());
    std::int64_t end_index = -1;
    for (const std::uint32_t chrom : order) {
        const std::int64_t start_index = end_index + 1;
        std::vector<std::size_t>& group = groups[chrom];
        std::sort(group.begin(), group.end(),
                  [&bins](std::size_t a, std::size_t b) {
                      if (bins[a].start != bins[b].start) {
                          return bins[a].start < bins[b].start;
                      }
                      if (bins[a].end != bins[b].end) {
                          return bins[a].end < bins[b].end;
                      }
                      return a < b;  // the bin id is the tuple's third element
                  });
        for (const std::size_t i : group) {
            flat_.push_back({static_cast<std::uint32_t>(bins[i].start),
                             static_cast<std::uint32_t>(bins[i].end),
                             static_cast<std::uint32_t>(i)});
        }
        end_index = start_index + static_cast<std::int64_t>(group.size()) - 1;
        range_.emplace(chrom, std::make_pair(start_index, end_index));
    }
}

bool BinSearchIndex::has_chrom(std::uint32_t chrom) const {
    return range_.find(chrom) != range_.end();
}

std::optional<std::int64_t> BinSearchIndex::flat_index_at(
    std::uint32_t chrom, std::int64_t position) const {
    const auto found = range_.find(chrom);
    if (found == range_.end()) {
        return std::nullopt;  // the KeyError branch
    }
    std::int64_t start = found->second.first;
    std::int64_t end = found->second.second;
    // int((start + end) / 2) truncates toward zero, which differs from a floor
    // division as soon as the sum is negative. It cannot be here, but the
    // truncating form is what the Python has.
    std::int64_t middle = (start + end) / 2;
    while (!(start > end)) {
        const SearchInterval& probe = flat_[static_cast<std::size_t>(middle)];
        const std::uint64_t value = static_cast<std::uint64_t>(position);
        if (static_cast<std::uint64_t>(probe.begin) <= value &&
            value <= static_cast<std::uint64_t>(probe.end)) {
            return middle;
        }
        if (static_cast<std::uint64_t>(probe.begin) > value) {
            end = middle - 1;
        } else {
            start = middle + 1;
        }
        middle = (start + end) / 2;
    }
    return std::nullopt;
}

std::optional<std::uint32_t> BinSearchIndex::bin_at(std::uint32_t chrom,
                                                    std::int64_t position) const {
    const auto index = flat_index_at(chrom, position);
    if (!index.has_value()) {
        return std::nullopt;
    }
    return flat_[static_cast<std::size_t>(*index)].data;
}

RestrictionSiteIndex::RestrictionSiteIndex(const std::vector<GenomeInterval>& sites,
                                           std::size_t chrom_count) {
    per_chrom_.resize(chrom_count);
    present_.assign(chrom_count, 0);
    std::vector<std::vector<std::pair<std::int64_t, std::int64_t>>> staged(chrom_count);
    for (const auto& site : sites) {
        staged[site.chrom].emplace_back(site.start, site.end);
        present_[site.chrom] = 1;
    }
    for (std::size_t chrom = 0; chrom < chrom_count; ++chrom) {
        auto& items = staged[chrom];
        std::sort(items.begin(), items.end());
        Chrom& target = per_chrom_[chrom];
        target.begin.reserve(items.size());
        target.prefix_max_end.assign(items.size() + 1,
                                     std::numeric_limits<std::int64_t>::min());
        for (std::size_t i = 0; i < items.size(); ++i) {
            target.begin.push_back(items[i].first);
            target.prefix_max_end[i + 1] =
                std::max(target.prefix_max_end[i], items[i].second);
        }
        total_ += items.size();
    }
}

bool RestrictionSiteIndex::has_chrom(std::uint32_t chrom) const {
    return chrom < present_.size() && present_[chrom] != 0;
}

bool RestrictionSiteIndex::overlaps(std::uint32_t chrom, std::int64_t begin,
                                    std::int64_t end) const {
    if (chrom >= per_chrom_.size() || begin >= end) {
        return false;  // intervaltree returns an empty set for a null slice
    }
    const Chrom& target = per_chrom_[chrom];
    // Intervals with begin < end of the query.
    const auto upper =
        std::lower_bound(target.begin.begin(), target.begin.end(), end);
    const auto count = static_cast<std::size_t>(upper - target.begin.begin());
    return target.prefix_max_end[count] > begin;
}

std::string reverse_complement(const std::string& sequence) {
    static const auto complement = [](char base) -> char {
        switch (base) {
            case 'A': return 'T';
            case 'T': return 'A';
            case 'G': return 'C';
            case 'C': return 'G';
            case 'U': return 'A';
            case 'M': return 'K';
            case 'K': return 'M';
            case 'R': return 'Y';
            case 'Y': return 'R';
            case 'W': return 'W';
            case 'S': return 'S';
            case 'B': return 'V';
            case 'V': return 'B';
            case 'D': return 'H';
            case 'H': return 'D';
            case 'N': return 'N';
            case 'a': return 't';
            case 't': return 'a';
            case 'g': return 'c';
            case 'c': return 'g';
            default: return base;
        }
    };
    std::string out;
    out.reserve(sequence.size());
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        out.push_back(complement(*it));
    }
    return out;
}

void QcCounters::add(const QcCounters& other) {
    one_mate_unmapped += other.one_mate_unmapped;
    one_mate_low_quality += other.one_mate_low_quality;
    one_mate_not_unique += other.one_mate_not_unique;
    duplicated_pairs += other.duplicated_pairs;
    if (dangling_end.size() < other.dangling_end.size()) {
        dangling_end.resize(other.dangling_end.size(), 0);
    }
    for (std::size_t i = 0; i < other.dangling_end.size(); ++i) {
        dangling_end[i] += other.dangling_end[i];
    }
    self_circle += other.self_circle;
    self_ligation += other.self_ligation;
    same_fragment += other.same_fragment;
    mate_not_close_to_rf += other.mate_not_close_to_rf;
    count_inward += other.count_inward;
    count_outward += other.count_outward;
    count_left += other.count_left;
    count_right += other.count_right;
    inter_chromosomal += other.inter_chromosomal;
    short_range += other.short_range;
    long_range += other.long_range;
    pair_added += other.pair_added;
    iter_num += other.iter_num;
}

std::string format_qc_log(const QcLogInputs& inputs, const QcCounters& counters) {
    const double iter_num = static_cast<double>(counters.iter_num);
    const std::int64_t mappable =
        counters.iter_num - (counters.one_mate_unmapped +
                             counters.one_mate_low_quality + counters.one_mate_not_unique);
    const double mappable_f = static_cast<double>(mappable);
    const std::string msg =
        inputs.keep_self_ligation ? " (not removed)" : " (removed)";

    std::string out;
    auto line = [&out](const std::string& text) { out += text; };

    if (inputs.min_distance != 0) {
        line("\nFile\t" + inputs.out_file_name + "\t\t\nSequenced reads\t" +
             std::to_string(counters.iter_num) + "\t\t\nMin rest. site distance\t" +
             std::to_string(inputs.min_distance) + "\t\t\nMax library insert size\t" +
             std::to_string(inputs.max_library_insert_size) + "\t\t\n\n");
    } else {
        line("\nFile\t" + inputs.out_file_name + "\t\t\nSequenced reads\t" +
             std::to_string(counters.iter_num) + "\t\t\nMax library insert size\t" +
             std::to_string(inputs.max_library_insert_size) + "\t\t\n\n");
    }
    line("#\tcount\t(percentage w.r.t. total sequenced reads)\n");
    line("Pairs mappable, unique and high quality\t" + std::to_string(mappable) +
         "\t(" + format_two_decimals(100 * mappable_f / iter_num) + ")\n");
    line("Hi-C contacts\t" + std::to_string(counters.pair_added) + "\t(" +
         format_two_decimals(100 * static_cast<double>(counters.pair_added) / iter_num) +
         ")\n");
    line("One mate unmapped\t" + std::to_string(counters.one_mate_unmapped) + "\t(" +
         format_two_decimals(100 * static_cast<double>(counters.one_mate_unmapped) /
                             iter_num) +
         ")\n");
    line("One mate not unique\t" + std::to_string(counters.one_mate_not_unique) +
         "\t(" +
         format_two_decimals(100 * static_cast<double>(counters.one_mate_not_unique) /
                             iter_num) +
         ")\n");
    line("Low mapping quality\t" + std::to_string(counters.one_mate_low_quality) +
         "\t(" +
         format_two_decimals(100 * static_cast<double>(counters.one_mate_low_quality) /
                             iter_num) +
         ")\n");
    line("\n#\tcount\t(percentage w.r.t. mappable, unique and high quality "
         "pairs)\n");

    for (std::size_t i = 0; i < inputs.dangling_sequences.size(); ++i) {
        const auto& [restriction, dangling] = inputs.dangling_sequences[i];
        const std::int64_t value =
            i < counters.dangling_end.size() ? counters.dangling_end[i] : 0;
        line("dangling end " + dangling + " (restriction sequence " + restriction +
             ")\t" + std::to_string(value) + "\t(" +
             format_two_decimals(100 * static_cast<double>(value) / mappable_f) + ")\n");
    }
    if (inputs.has_restriction_cut_file) {
        line("self ligation" + msg + "\t" + std::to_string(counters.self_ligation) +
             "\t(" +
             format_two_decimals(100 * static_cast<double>(counters.self_ligation) /
                                 mappable_f) +
             ")\n");
        line("One mate not close to rest site\t" +
             std::to_string(counters.mate_not_close_to_rf) + "\t(" +
             format_two_decimals(100 * static_cast<double>(counters.mate_not_close_to_rf) /
                                 mappable_f) +
             ")\n");
    }
    line("same fragment\t" + std::to_string(counters.same_fragment) + "\t(" +
         format_two_decimals(100 * static_cast<double>(counters.same_fragment) /
                             mappable_f) +
         ")\n");
    if (inputs.has_restriction_cut_file) {
        line("self circle\t" + std::to_string(counters.self_circle) + "\t(" +
             format_two_decimals(100 * static_cast<double>(counters.self_circle) /
                                 mappable_f) +
             ")\n");
    }
    line("duplicated pairs\t" + std::to_string(counters.duplicated_pairs) + "\t(" +
         format_two_decimals(100 * static_cast<double>(counters.duplicated_pairs) /
                             mappable_f) +
         ")\n");

    if (counters.pair_added > 0) {
        const double added = static_cast<double>(counters.pair_added);
        line("\n#\tcount\t(percentage w.r.t. total valid pairs used)\n");
        line("inter chromosomal\t" + std::to_string(counters.inter_chromosomal) +
             "\t(" +
             format_two_decimals(100 * static_cast<double>(counters.inter_chromosomal) /
                                 added) +
             ")\n");
        line("Intra short range (< 20kb)\t" + std::to_string(counters.short_range) +
             "\t(" +
             format_two_decimals(100 * static_cast<double>(counters.short_range) / added) +
             ")\n");
        line("Intra long range (>= 20kb)\t" + std::to_string(counters.long_range) +
             "\t(" +
             format_two_decimals(100 * static_cast<double>(counters.long_range) / added) +
             ")\n");
        line("Read pair type: inward pairs\t" + std::to_string(counters.count_inward) +
             "\t(" +
             format_two_decimals(100 * static_cast<double>(counters.count_inward) / added) +
             ")\n");
        line("Read pair type: outward pairs\t" +
             std::to_string(counters.count_outward) + "\t(" +
             format_two_decimals(100 * static_cast<double>(counters.count_outward) / added) +
             ")\n");
        line("Read pair type: left pairs\t" + std::to_string(counters.count_left) +
             "\t(" +
             format_two_decimals(100 * static_cast<double>(counters.count_left) / added) +
             ")\n");
        line("Read pair type: right pairs\t" + std::to_string(counters.count_right) +
             "\t(" +
             format_two_decimals(100 * static_cast<double>(counters.count_right) / added) +
             ")\n");
    }
    return out;
}

namespace {

// hicPrepareQCreport.main's parser: from the first line starting with "File",
// every non empty, non comment line with at least two tab separated fields
// becomes one column.
struct QcTable {
    std::vector<std::string> columns;  // in insertion order, "File" first
    std::unordered_map<std::string, std::string> text;
    std::unordered_map<std::string, double> number;

    [[nodiscard]] bool has(const std::string& column) const {
        return number.find(column) != number.end() ||
               text.find(column) != text.end();
    }
};

QcTable parse_qc_log(const std::string& log_text) {
    QcTable table;
    bool in_log_part = false;
    for (const auto& raw : split(log_text, '\n')) {
        std::string line = trim(raw);
        if (line.rfind("File", 0) == 0) {
            in_log_part = true;
        }
        if (!in_log_part) {
            continue;
        }
        if (line.empty() || line.rfind('#', 0) == 0) {
            continue;
        }
        const std::vector<std::string> fields = split(line, '\t');
        if (fields.size() == 1) {
            continue;
        }
        const std::string& key = fields[0];
        if (!table.has(key)) {
            table.columns.push_back(key);
        }
        const std::string& value = fields[1];
        char* stop = nullptr;
        const long long parsed = std::strtoll(value.c_str(), &stop, 10);
        if (stop != nullptr && *stop == '\0' && !value.empty()) {
            table.number[key] = static_cast<double>(parsed);
        } else {
            table.text[key] = value;
        }
    }
    return table;
}

void write_table(const std::string& path, const std::string& index_value,
                 const std::vector<std::string>& columns,
                 const std::vector<std::string>& values) {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("could not write " + path);
    }
    out << "File";
    for (const auto& column : columns) {
        out << '\t' << column;
    }
    out << '\n' << index_value;
    for (const auto& value : values) {
        out << '\t' << value;
    }
    out << '\n';
}

std::string integer_field(double value) {
    return std::to_string(static_cast<long long>(value));
}

}  // namespace

void write_qc_tables(const std::string& folder, const std::string& qc_log_text) {
    const QcTable table = parse_qc_log(qc_log_text);
    const auto file_entry = table.text.find("File");
    if (file_entry == table.text.end()) {
        throw std::runtime_error("the QC log has no File line");
    }
    const std::string index_value = file_entry->second;

    auto value_of = [&table](const std::string& column) -> double {
        const auto found = table.number.find(column);
        return found == table.number.end() ? 0.0 : found->second;
    };

    // QC_table.txt: every column, in the order of the log.
    {
        std::vector<std::string> columns;
        std::vector<std::string> values;
        for (const auto& column : table.columns) {
            if (column == "File") {
                continue;
            }
            columns.push_back(column);
            values.push_back(integer_field(value_of(column)));
        }
        write_table(folder + "/QC_table.txt", index_value, columns, values);
    }

    const double sequenced = value_of("Sequenced reads");
    const double mappable = value_of("Pairs mappable, unique and high quality");
    const double contacts = value_of("Hi-C contacts");

    // unmapable_table.txt
    {
        static const char* const names[] = {"Hi-C contacts", "Low mapping quality",
                                            "One mate not unique",
                                            "One mate unmapped"};
        std::vector<std::string> columns;
        std::vector<std::string> values;
        for (const char* name : names) {
            columns.emplace_back(name);
            values.push_back(integer_field(value_of(name)));
            columns.emplace_back(std::string(name) + "_%");
            values.push_back(python_float_repr(value_of(name) / sequenced));
        }
        write_table(folder + "/unmapable_table.txt", index_value, columns, values);
    }

    // discarded_table.txt, with the prefix matching of
    // hicPrepareQCreport.make_figure_pairs_discarded: an exact column match
    // emits the count and the percentage, a substring match emits only the
    // percentage.
    {
        static const char* const prefixes[] = {
            "One mate not close to rest site", "dangling end", "duplicated pairs",
            "same fragment", "self circle", "self ligation (removed)"};
        std::vector<std::string> columns;
        std::vector<std::string> values;
        for (const char* prefix : prefixes) {
            const std::string name(prefix);
            if (table.has(name)) {
                columns.push_back(name);
                values.push_back(integer_field(value_of(name)));
                columns.push_back(name + " %");
                values.push_back(python_float_repr(value_of(name) / mappable));
                continue;
            }
            for (const auto& column : table.columns) {
                if (column.find(name) != std::string::npos) {
                    columns.push_back(column + " %");
                    values.push_back(python_float_repr(value_of(column) / mappable));
                }
            }
        }
        write_table(folder + "/discarded_table.txt", index_value, columns, values);
    }

    // distance_table.txt
    {
        static const char* const names[] = {"inter chromosomal",
                                            "Intra short range (< 20kb)",
                                            "Intra long range (>= 20kb)"};
        std::vector<std::string> columns;
        std::vector<std::string> values;
        for (const char* name : names) {
            columns.emplace_back(name);
            values.push_back(integer_field(value_of(name)));
            columns.emplace_back(std::string(name) + " %");
            values.push_back(python_float_repr(value_of(name) / contacts));
        }
        write_table(folder + "/distance_table.txt", index_value, columns, values);
    }

    // read_orientation_table.txt, whose percentages are relative to the sum of
    // the four orientations rather than to the contact count.
    {
        static const char* const names[] = {"Read pair type: inward pairs",
                                            "Read pair type: outward pairs",
                                            "Read pair type: left pairs",
                                            "Read pair type: right pairs"};
        double total = 0.0;
        for (const char* name : names) {
            total += value_of(name);
        }
        std::vector<std::string> columns;
        std::vector<std::string> values;
        for (const char* name : names) {
            columns.emplace_back(name);
            values.push_back(integer_field(value_of(name)));
            columns.emplace_back(std::string(name) + " %");
            values.push_back(python_float_repr(value_of(name) / total));
        }
        write_table(folder + "/read_orientation_table.txt", index_value, columns,
                    values);
    }
}

}  // namespace hicx
