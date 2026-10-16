#include "hicx/build_matrix.hpp"

#include "hicx/argparse.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <limits>
#include <optional>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

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

// A cell of hicPrepareQCreport's DataFrame: int(fields[1]) where Python's
// int() accepts the text, the text otherwise.
struct QcCell {
    bool is_int = false;
    std::int64_t integer = 0;
    std::string text;

    [[nodiscard]] std::string str() const { return is_int ? std::to_string(integer) : text; }
};

// The DataFrame of hicPrepareQCreport.main after set_index: named columns of
// cells in insertion order, and the index.
struct QcFrame {
    std::string index_name;
    std::vector<std::string> index;
    std::vector<std::string> columns;
    std::vector<std::vector<QcCell>> cells;  // per column, per row

    [[nodiscard]] std::ptrdiff_t find(const std::string& name) const {
        const auto it = std::find(columns.begin(), columns.end(), name);
        return it == columns.end() ? -1 : it - columns.begin();
    }
    [[nodiscard]] bool has(const std::string& name) const { return find(name) >= 0; }
    [[nodiscard]] const std::vector<QcCell>& column(const std::string& name) const {
        const std::ptrdiff_t position = find(name);
        if (position < 0) {
            throw std::runtime_error("KeyError: \"['" + name + "'] not in index\"");
        }
        return cells[static_cast<std::size_t>(position)];
    }
    // A value of a column as the number pandas divides with.
    [[nodiscard]] double number(const std::string& name, std::size_t row) const {
        const QcCell& cell = column(name)[row];
        if (!cell.is_int) {
            throw std::runtime_error("TypeError: unsupported operand type(s) for /: 'str' and "
                                     "'int' (the QC column '" + name + "' holds text)");
        }
        return static_cast<double>(cell.integer);
    }
};

// pandas.DataFrame.to_csv quotes a field that holds the separator, a quote or
// a line break.
std::string csv_field(const std::string& text) {
    if (text.find_first_of("\t\"\n\r") == std::string::npos) {
        return text;
    }
    std::string out = "\"";
    for (const char c : text) {
        out += c;
        if (c == '"') {
            out += '"';
        }
    }
    return out + "\"";
}

// hicPrepareQCreport.main from reading the logs to set_index and the two
// derived columns.
QcFrame build_qc_frame(const std::vector<std::string>& log_texts,
                       const std::optional<std::vector<std::string>>& labels) {
    std::vector<std::string> order;
    std::unordered_map<std::string, std::vector<QcCell>> params;
    for (const std::string& log_text : log_texts) {
        bool in_log_part = false;
        std::vector<std::string> lines = split(log_text, '\n');
        // readlines() yields no empty element after a final newline.
        if (!lines.empty() && lines.back().empty()) {
            lines.pop_back();
        }
        for (const std::string& raw : lines) {
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
            auto [entry, inserted] = params.try_emplace(fields[0]);
            if (inserted) {
                order.push_back(fields[0]);
            }
            QcCell cell;
            std::int64_t value = 0;
            if (cli::python_int(fields[1], &value)) {
                cell.is_int = true;
                cell.integer = value;
            } else {
                cell.text = fields[1];
            }
            entry->second.push_back(std::move(cell));
        }
    }

    QcFrame frame;
    std::size_t rows = 0;
    for (std::size_t i = 0; i < order.size(); ++i) {
        const std::size_t length = params[order[i]].size();
        if (i == 0) {
            rows = length;
        } else if (length != rows) {
            throw std::runtime_error("ValueError: All arrays must be of the same length (the QC "
                                     "logs do not hold the same lines)");
        }
    }
    for (const std::string& name : order) {
        frame.columns.push_back(name);
        frame.cells.push_back(std::move(params[name]));
    }

    if (labels.has_value() && labels->size() == log_texts.size()) {
        if (labels->size() != rows) {
            throw std::runtime_error(
                "*ERROR* Some log files may not be valid. Please check that the log files "
                "contain at the end the summary information. The reference exits 0 here "
                "without writing anything; the C++ port exits 1 (cpp/AGENTS_CONTRACT.md rule 7).");
        }
        frame.index_name = "Labels";
        frame.index = *labels;
    } else {
        const std::ptrdiff_t position = frame.find("File");
        if (position < 0) {
            throw std::runtime_error("KeyError: \"None of ['File'] are in the columns\"");
        }
        frame.index_name = "File";
        for (const QcCell& cell : frame.cells[static_cast<std::size_t>(position)]) {
            frame.index.push_back(cell.str());
        }
        frame.columns.erase(frame.columns.begin() + position);
        frame.cells.erase(frame.cells.begin() + position);
    }
    const std::string mappable = "Pairs mappable, unique and high quality";
    if (!frame.has(mappable)) {
        std::vector<QcCell> derived(rows);
        for (std::size_t row = 0; row < rows; ++row) {
            auto integer_of = [&](const char* name) {
                const QcCell& cell = frame.column(name)[row];
                if (!cell.is_int) {
                    throw std::runtime_error(std::string("TypeError: the QC column '") + name +
                                             "' holds text");
                }
                return cell.integer;
            };
            derived[row].is_int = true;
            derived[row].integer =
                integer_of("Sequenced reads") -
                (integer_of("One mate unmapped") + integer_of("One mate not unique") +
                 integer_of("Low mapping quality"));
        }
        frame.columns.push_back(mappable);
        frame.cells.push_back(std::move(derived));
    }
    if (frame.has("same fragment (800 bp)")) {
        std::vector<QcCell> copy = frame.column("same fragment (800 bp)");
        const std::ptrdiff_t position = frame.find("same fragment");
        if (position >= 0) {
            frame.cells[static_cast<std::size_t>(position)] = std::move(copy);
        } else {
            frame.columns.push_back("same fragment");
            frame.cells.push_back(std::move(copy));
        }
    }
    return frame;
}

// One output table: its columns, each already rendered per output row.
struct QcOutputTable {
    std::vector<std::string> columns;
    std::vector<std::vector<std::string>> values;  // per column, per output row

    void add(std::string name, std::vector<std::string> column) {
        columns.push_back(std::move(name));
        values.push_back(std::move(column));
    }
};

// The rows of table[counts].join(prc_table.T, rsuffix=...), a left join on
// the index: for every row, in order, every row with the same name, in order.
// With unique names that is every row once; the same log twice without
// labels gives four rows.
std::vector<std::pair<std::size_t, std::size_t>> joined_rows(const QcFrame& frame) {
    std::vector<std::pair<std::size_t, std::size_t>> rows;
    for (std::size_t left = 0; left < frame.index.size(); ++left) {
        for (std::size_t right = 0; right < frame.index.size(); ++right) {
            if (frame.index[right] == frame.index[left]) {
                rows.emplace_back(left, right);
            }
        }
    }
    return rows;
}

std::vector<std::pair<std::size_t, std::size_t>> plain_rows(const QcFrame& frame) {
    std::vector<std::pair<std::size_t, std::size_t>> rows;
    for (std::size_t row = 0; row < frame.index.size(); ++row) {
        rows.emplace_back(row, row);
    }
    return rows;
}

void write_qc_output(const std::string& path, const QcFrame& frame, const QcOutputTable& table,
                     const std::vector<std::pair<std::size_t, std::size_t>>& rows) {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("could not write " + path);
    }
    out << csv_field(frame.index_name);
    for (const std::string& column : table.columns) {
        out << '\t' << csv_field(column);
    }
    out << '\n';
    for (std::size_t row = 0; row < rows.size(); ++row) {
        out << csv_field(frame.index[rows[row].first]);
        for (const auto& column : table.values) {
            out << '\t' << csv_field(column[row]);
        }
        out << '\n';
    }
    if (!out) {
        throw std::runtime_error("could not write " + path);
    }
}

// The counts of the left rows.
std::vector<std::string> count_column(const QcFrame& frame, const std::string& name,
                                      const std::vector<std::pair<std::size_t, std::size_t>>& rows) {
    std::vector<std::string> out;
    const std::vector<QcCell>& cells = frame.column(name);
    for (const auto& row : rows) {
        out.push_back(cells[row.first].str());
    }
    return out;
}

// The fractions of the right rows, each over its own row's denominator.
std::vector<std::string> ratio_column(const QcFrame& frame, const std::string& name,
                                      const std::vector<double>& denominators,
                                      const std::vector<std::pair<std::size_t, std::size_t>>& rows) {
    std::vector<std::string> out;
    for (const auto& row : rows) {
        out.push_back(python_float_repr(frame.number(name, row.second) / denominators[row.second]));
    }
    return out;
}

std::vector<double> numbers(const QcFrame& frame, const std::string& name) {
    std::vector<double> out;
    for (std::size_t row = 0; row < frame.index.size(); ++row) {
        out.push_back(frame.number(name, row));
    }
    return out;
}

}  // namespace

void write_qc_report_tables(const std::string& folder, const std::vector<std::string>& log_texts,
                            const std::optional<std::vector<std::string>>& labels) {
    const QcFrame frame = build_qc_frame(log_texts, labels);
    const auto plain = plain_rows(frame);
    const auto joined = joined_rows(frame);

    // QC_table.txt: table.to_csv, every column.
    {
        QcOutputTable table;
        for (const std::string& column : frame.columns) {
            table.add(column, count_column(frame, column, plain));
        }
        write_qc_output(folder + "/QC_table.txt", frame, table, plain);
    }

    // make_figure_umappable_non_unique_reads: counts and '_%' over the
    // sequenced reads.
    {
        const std::vector<double> sequenced = numbers(frame, "Sequenced reads");
        QcOutputTable table;
        for (const char* name : {"Hi-C contacts", "Low mapping quality", "One mate not unique",
                                 "One mate unmapped"}) {
            table.add(name, count_column(frame, name, joined));
            table.add(std::string(name) + "_%", ratio_column(frame, name, sequenced, joined));
        }
        write_qc_output(folder + "/unmapable_table.txt", frame, table, joined);
    }

    // make_figure_pairs_discarded: an exact column match gives the count and
    // ' %', every column that only contains the prefix gives ' %' alone.
    {
        const std::vector<double> mappable =
            numbers(frame, "Pairs mappable, unique and high quality");
        QcOutputTable table;
        for (const char* prefix : {"One mate not close to rest site", "dangling end",
                                   "duplicated pairs", "same fragment", "self circle",
                                   "self ligation (removed)"}) {
            const std::string name(prefix);
            if (frame.has(name)) {
                table.add(name, count_column(frame, name, joined));
                table.add(name + " %", ratio_column(frame, name, mappable, joined));
                continue;
            }
            for (const std::string& column : frame.columns) {
                if (column.find(name) != std::string::npos) {
                    table.add(column + " %", ratio_column(frame, column, mappable, joined));
                }
            }
        }
        write_qc_output(folder + "/discarded_table.txt", frame, table, joined);
    }

    // make_figure_distance: over the Hi-C contacts.
    {
        const std::vector<double> contacts = numbers(frame, "Hi-C contacts");
        QcOutputTable table;
        for (const char* name : {"inter chromosomal", "Intra short range (< 20kb)",
                                 "Intra long range (>= 20kb)"}) {
            table.add(name, count_column(frame, name, joined));
            table.add(std::string(name) + " %", ratio_column(frame, name, contacts, joined));
        }
        write_qc_output(folder + "/distance_table.txt", frame, table, joined);
    }

    // make_figure_read_orientation: over the sum of the four orientations.
    {
        static const char* const names[] = {"Read pair type: inward pairs",
                                            "Read pair type: outward pairs",
                                            "Read pair type: left pairs",
                                            "Read pair type: right pairs"};
        std::vector<double> totals(frame.index.size(), 0.0);
        for (const char* name : names) {
            const std::vector<double> values = numbers(frame, name);
            for (std::size_t row = 0; row < totals.size(); ++row) {
                totals[row] += values[row];
            }
        }
        QcOutputTable table;
        for (const char* name : names) {
            table.add(name, count_column(frame, name, joined));
            table.add(std::string(name) + " %", ratio_column(frame, name, totals, joined));
        }
        write_qc_output(folder + "/read_orientation_table.txt", frame, table, joined);
    }
}

void write_qc_tables(const std::string& folder, const std::string& qc_log_text) {
    write_qc_report_tables(folder, {qc_log_text}, std::nullopt);
}

namespace {

// The value of the first column that is exactly `name`, or that contains it
// (the same match write_qc_report_tables's discarded_table block uses for a
// column such as "dangling end" that the log may suffix, e.g. with a cut
// length). Returns nullopt when no such column exists for this row.
std::optional<double> qc_value_by_name_or_prefix(const QcFrame& frame, const std::string& name,
                                                  std::size_t row) {
    if (frame.has(name)) {
        return frame.number(name, row);
    }
    for (const std::string& column : frame.columns) {
        if (column.find(name) != std::string::npos) {
            return frame.number(column, row);
        }
    }
    return std::nullopt;
}

enum class QcStatus { Pass, Warn, Fail, Unknown };

const char* qc_status_label(QcStatus status) {
    switch (status) {
        case QcStatus::Pass:
            return "PASS";
        case QcStatus::Warn:
            return "WARN";
        case QcStatus::Fail:
            return "FAIL";
        default:
            return "N/A";
    }
}

const char* qc_status_class(QcStatus status) {
    switch (status) {
        case QcStatus::Pass:
            return "qc-pass";
        case QcStatus::Warn:
            return "qc-warn";
        case QcStatus::Fail:
            return "qc-fail";
        default:
            return "qc-unknown";
    }
}

// `value` passes when it is at least `minimum`; within `margin` of the
// shortfall below `minimum` (relative to `minimum`) is WARN, further below is
// FAIL. Used for metrics a report wants to see enough of (valid pairs).
QcStatus qc_status_at_least(std::optional<double> value, double minimum, double margin) {
    if (!value.has_value()) {
        return QcStatus::Unknown;
    }
    if (*value >= minimum) {
        return QcStatus::Pass;
    }
    if (*value >= minimum * (1.0 - margin)) {
        return QcStatus::Warn;
    }
    return QcStatus::Fail;
}

// `value` passes when it is at most `maximum`; within `margin` above
// `maximum` (relative to `maximum`) is WARN, further above is FAIL. Used for
// metrics a report wants to see little of (duplicates, dangling ends,
// self-ligation, inter-chromosomal contacts).
QcStatus qc_status_at_most(std::optional<double> value, double maximum, double margin) {
    if (!value.has_value()) {
        return QcStatus::Unknown;
    }
    if (*value <= maximum) {
        return QcStatus::Pass;
    }
    if (*value <= maximum * (1.0 + margin)) {
        return QcStatus::Warn;
    }
    return QcStatus::Fail;
}

std::string qc_html_escape(const std::string& text) {
    std::string out;
    out.reserve(text.size());
    for (const char c : text) {
        switch (c) {
            case '&':
                out += "&amp;";
                break;
            case '<':
                out += "&lt;";
                break;
            case '>':
                out += "&gt;";
                break;
            case '"':
                out += "&quot;";
                break;
            default:
                out += c;
        }
    }
    return out;
}

std::string qc_percent(std::optional<double> value) {
    if (!value.has_value()) {
        return "n/a";
    }
    char buffer[32];
    std::snprintf(buffer, sizeof(buffer), "%.1f%%", *value);
    return buffer;
}

struct QcMetricRow {
    std::string label;
    std::optional<double> value_percent;
    std::string rule;  // shown under the metric, e.g. "cutoff: >= 25%"
    QcStatus status;
};

void qc_write_metric_row(std::ostream& out, const QcMetricRow& row) {
    out << "<div class=\"qc-metric\"><div class=\"qc-metric-top\">"
        << "<span class=\"qc-metric-label\">" << qc_html_escape(row.label) << "</span>"
        << "<span class=\"qc-badge " << qc_status_class(row.status) << "\">"
        << qc_status_label(row.status) << "</span></div>"
        << "<div class=\"qc-metric-value\">" << qc_percent(row.value_percent) << "</div>"
        << "<div class=\"qc-metric-rule\">" << qc_html_escape(row.rule) << "</div></div>\n";
}

// Every column that is exactly `prefix`, or (when none is) every column that
// contains it, each with its own value -- the discarded-pairs columns can be
// per restriction enzyme, e.g. two "dangling end" columns for a two-enzyme
// run, and the original report lists each, not their sum.
std::vector<std::pair<std::string, double>> qc_values_by_prefix_all(const QcFrame& frame,
                                                                     const std::string& prefix,
                                                                     std::size_t row) {
    std::vector<std::pair<std::string, double>> out;
    if (frame.has(prefix)) {
        out.emplace_back(prefix, frame.number(prefix, row));
        return out;
    }
    for (const std::string& column : frame.columns) {
        if (column.find(prefix) != std::string::npos) {
            out.emplace_back(column, frame.number(column, row));
        }
    }
    return out;
}

// 1,234,567-style grouping, matching the original report's "{:,}" format.
std::string qc_format_count(double value) {
    const auto whole = static_cast<long long>(value + (value >= 0 ? 0.5 : -0.5));
    std::string digits = std::to_string(whole < 0 ? -whole : whole);
    std::string grouped;
    int since_sep = 0;
    for (auto it = digits.rbegin(); it != digits.rend(); ++it) {
        if (since_sep == 3) {
            grouped += ',';
            since_sep = 0;
        }
        grouped += *it;
        ++since_sep;
    }
    std::reverse(grouped.begin(), grouped.end());
    return (whole < 0 ? "-" : "") + grouped;
}

// One row of a bar-chart section: a label, its count, its value as a
// fraction of `denominator` (nullopt when there is nothing to divide by),
// and an optional PASS/WARN/FAIL badge (Unknown status renders no badge).
struct QcBarRow {
    std::string label;
    double count = 0.0;
    std::optional<double> denominator;
    QcStatus status = QcStatus::Unknown;
};

void qc_write_bar_section(std::ostream& out, const std::vector<QcBarRow>& rows) {
    double max_fraction = 0.001;
    for (const QcBarRow& row : rows) {
        if (row.denominator.has_value() && *row.denominator > 0) {
            max_fraction = std::max(max_fraction, row.count / *row.denominator);
        }
    }
    out << "<div class=\"qc-bars\">\n";
    for (const QcBarRow& row : rows) {
        const std::optional<double> fraction =
            (row.denominator.has_value() && *row.denominator > 0)
                ? std::optional<double>(row.count / *row.denominator)
                : std::nullopt;
        const double width_percent = fraction.has_value() ? 100.0 * (*fraction) / max_fraction : 0.0;
        out << "<div class=\"qc-bar-row\"><div class=\"qc-bar-label\">"
            << qc_html_escape(row.label);
        if (row.status != QcStatus::Unknown) {
            out << " <span class=\"qc-badge qc-badge-sm " << qc_status_class(row.status) << "\">"
                << qc_status_label(row.status) << "</span>";
        }
        out << "</div><div class=\"qc-bar-track\"><div class=\"qc-bar-fill "
            << (row.status == QcStatus::Fail  ? "qc-bar-fail"
               : row.status == QcStatus::Warn ? "qc-bar-warn"
                                              : "qc-bar-plain")
            << "\" style=\"width:" << width_percent << "%\"></div></div>"
            << "<div class=\"qc-bar-value\">" << qc_format_count(row.count);
        if (fraction.has_value()) {
            out << " (" << qc_percent(100.0 * (*fraction)) << ")";
        }
        out << "</div></div>\n";
    }
    out << "</div>\n";
}

struct QcTableColumn {
    std::string header;
    std::vector<std::string> values;  // one per row, already formatted
};

void qc_write_table(std::ostream& out, const std::string& table_id,
                    const std::vector<std::string>& row_labels,
                    const std::vector<QcTableColumn>& columns) {
    out << "<div class=\"qc-table-wrap\"><table class=\"qc-table\" id=\"" << table_id
        << "\"><thead><tr><th>Sample</th>";
    for (const QcTableColumn& column : columns) {
        out << "<th>" << qc_html_escape(column.header) << "</th>";
    }
    out << "</tr></thead><tbody>\n";
    for (std::size_t row = 0; row < row_labels.size(); ++row) {
        out << "<tr><td class=\"qc-table-rowlabel\">" << qc_html_escape(row_labels[row])
            << "</td>";
        for (const QcTableColumn& column : columns) {
            out << "<td>" << qc_html_escape(column.values[row]) << "</td>";
        }
        out << "</tr>\n";
    }
    out << "</tbody></table></div>\n";
}

}  // namespace

void write_qc_pretty_report(const std::string& folder, const std::vector<std::string>& log_texts,
                            const std::optional<std::vector<std::string>>& labels,
                            const QcThresholds& thresholds) {
    const QcFrame frame = build_qc_frame(log_texts, labels);
    const double margin = thresholds.warn_margin_fraction;
    const std::size_t n = frame.index.size();

    // Per-sample metrics judged against a cutoff, computed once and reused
    // for the overview, the relevant section's bars and its badges.
    struct QcJudged {
        std::optional<double> valid_pct;
        std::optional<double> inter_chrom_pct;
        std::optional<double> duplicated_pct;
        std::optional<double> dangling_pct;
        std::optional<double> self_ligation_pct;
        QcStatus valid_status = QcStatus::Unknown;
        QcStatus inter_chrom_status = QcStatus::Unknown;
        QcStatus duplicated_status = QcStatus::Unknown;
        QcStatus dangling_status = QcStatus::Unknown;
        QcStatus self_ligation_status = QcStatus::Unknown;
    };
    std::vector<QcJudged> judged(n);
    for (std::size_t row = 0; row < n; ++row) {
        const double sequenced = frame.has("Sequenced reads") ? frame.number("Sequenced reads", row) : 0.0;
        const double valid = frame.has("Hi-C contacts") ? frame.number("Hi-C contacts", row) : 0.0;
        const double mappable = frame.has("Pairs mappable, unique and high quality")
                                    ? frame.number("Pairs mappable, unique and high quality", row)
                                    : 0.0;
        QcJudged& j = judged[row];
        j.valid_pct = sequenced > 0 ? std::optional<double>(100.0 * valid / sequenced) : std::nullopt;
        j.inter_chrom_pct = (valid > 0 && frame.has("inter chromosomal"))
                                ? std::optional<double>(100.0 * frame.number("inter chromosomal", row) / valid)
                                : std::nullopt;
        if (mappable > 0) {
            if (const auto value = qc_value_by_name_or_prefix(frame, "duplicated pairs", row)) {
                j.duplicated_pct = 100.0 * (*value) / mappable;
            }
            if (const auto value = qc_value_by_name_or_prefix(frame, "dangling end", row)) {
                j.dangling_pct = 100.0 * (*value) / mappable;
            }
            if (const auto value = qc_value_by_name_or_prefix(frame, "self ligation (removed)", row)) {
                j.self_ligation_pct = 100.0 * (*value) / mappable;
            }
        }
        j.valid_status = qc_status_at_least(j.valid_pct, thresholds.min_valid_pairs_percent, margin);
        j.inter_chrom_status =
            qc_status_at_most(j.inter_chrom_pct, thresholds.max_inter_chromosomal_percent, margin);
        j.duplicated_status =
            qc_status_at_most(j.duplicated_pct, thresholds.max_duplicated_pairs_percent, margin);
        j.dangling_status = qc_status_at_most(j.dangling_pct, thresholds.max_dangling_end_percent, margin);
        j.self_ligation_status =
            qc_status_at_most(j.self_ligation_pct, thresholds.max_self_ligation_percent, margin);
    }

    std::ofstream out(folder + "/hicQC_report.html", std::ios::binary);
    if (!out) {
        throw std::runtime_error("could not write " + folder + "/hicQC_report.html");
    }

    out << R"HTML(<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>HiCExplorer QC report</title>
<style>
:root {
  --qc-pass: #1b7a43; --qc-pass-bg: #e6f4ea;
  --qc-warn: #8a5a00; --qc-warn-bg: #fff3cd;
  --qc-fail: #a3242a; --qc-fail-bg: #fdecea;
  --qc-unknown: #555; --qc-unknown-bg: #eee;
  --qc-ink: #1f2430; --qc-muted: #5b6270; --qc-line: #e3e6eb; --qc-accent: #2a6099;
}
* { box-sizing: border-box; }
html { scroll-behavior: smooth; }
body {
  margin: 0; background: #f6f7f9; color: var(--qc-ink);
  font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
}
a { color: var(--qc-accent); }
#qc-nav {
  position: fixed; top: 0; left: 0; width: 230px; height: 100%; overflow-y: auto;
  background: #fff; border-right: 1px solid var(--qc-line); padding: 24px 18px;
}
#qc-nav h1 { font-size: 1.1rem; margin: 0 0 18px; }
#qc-nav a {
  display: block; padding: 6px 8px; border-radius: 6px; text-decoration: none;
  color: var(--qc-ink); font-size: 0.85rem; margin-bottom: 2px;
}
#qc-nav a:hover { background: #f0f2f5; }
#qc-nav a.qc-active { background: var(--qc-accent); color: #fff; }
#qc-content { margin-left: 230px; padding: 32px 28px 80px; max-width: 1200px; }
section { margin-bottom: 44px; }
section h2 { font-size: 1.25rem; margin: 0 0 4px; scroll-margin-top: 16px; }
.qc-sub { color: var(--qc-muted); margin: 0 0 16px; font-size: 0.88rem; }
.qc-samples { display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 16px; }
.qc-card { background: #fff; border: 1px solid var(--qc-line); border-radius: 10px; padding: 16px 18px; box-shadow: 0 1px 2px rgba(0,0,0,0.04); }
.qc-card h3 { font-size: 1rem; margin: 0 0 2px; }
.qc-card .qc-overall { margin: 0 0 12px; font-size: 0.82rem; color: var(--qc-muted); }
.qc-metric { border-top: 1px solid var(--qc-line); padding: 8px 0; }
.qc-metric:first-of-type { border-top: none; }
.qc-metric-top { display: flex; justify-content: space-between; align-items: center; }
.qc-metric-label { font-size: 0.85rem; }
.qc-metric-value { font-size: 1.15rem; font-weight: 600; margin-top: 2px; }
.qc-metric-rule { font-size: 0.72rem; color: var(--qc-muted); margin-top: 1px; }
.qc-badge { font-size: 0.68rem; font-weight: 700; letter-spacing: 0.03em; padding: 2px 8px; border-radius: 999px; white-space: nowrap; }
.qc-badge-sm { font-size: 0.62rem; padding: 1px 6px; }
.qc-pass { color: var(--qc-pass); background: var(--qc-pass-bg); }
.qc-warn { color: var(--qc-warn); background: var(--qc-warn-bg); }
.qc-fail { color: var(--qc-fail); background: var(--qc-fail-bg); }
.qc-unknown { color: var(--qc-unknown); background: var(--qc-unknown-bg); }
.qc-bars { background: #fff; border: 1px solid var(--qc-line); border-radius: 10px; padding: 14px 18px; margin-bottom: 10px; }
.qc-bar-row { display: grid; grid-template-columns: 220px 1fr 160px; align-items: center; gap: 12px; padding: 7px 0; }
.qc-bar-label { font-size: 0.85rem; }
.qc-bar-track { background: #eef0f3; border-radius: 5px; height: 14px; overflow: hidden; }
.qc-bar-fill { height: 100%; border-radius: 5px; transition: width 0.4s ease; }
.qc-bar-plain { background: var(--qc-accent); }
.qc-bar-warn { background: #d8a200; }
.qc-bar-fail { background: #c0392b; }
.qc-bar-value { font-size: 0.8rem; color: var(--qc-muted); text-align: right; white-space: nowrap; }
.qc-note { font-size: 0.78rem; color: var(--qc-muted); line-height: 1.5; margin-top: 10px; }
dl.qc-defs { font-size: 0.85rem; background: #fff; border: 1px solid var(--qc-line); border-radius: 10px; padding: 14px 18px; margin: 0 0 16px; }
dl.qc-defs dt { font-weight: 700; margin-top: 8px; }
dl.qc-defs dt:first-child { margin-top: 0; }
dl.qc-defs dd { margin: 2px 0 0; color: var(--qc-muted); }
.qc-table-wrap { overflow-x: auto; background: #fff; border: 1px solid var(--qc-line); border-radius: 10px; }
table.qc-table { border-collapse: collapse; width: 100%; font-size: 0.82rem; }
table.qc-table th, table.qc-table td { padding: 8px 12px; text-align: right; white-space: nowrap; border-bottom: 1px solid var(--qc-line); }
table.qc-table th:first-child, table.qc-table td:first-child { text-align: left; position: sticky; left: 0; background: inherit; }
table.qc-table thead th { background: #f0f2f5; cursor: pointer; user-select: none; }
table.qc-table thead th:hover { background: #e4e7eb; }
table.qc-table tbody tr:nth-of-type(even) { background: #fafbfc; }
.qc-download { font-size: 0.78rem; margin-top: 8px; display: inline-block; }
</style>
</head>
<body>
<nav id="qc-nav">
<h1>HiCExplorer QC</h1>
<a href="#overview">Overview</a>
<a href="#sequenced">Pairs sequenced</a>
<a href="#unmappable">Unmappable / non-unique</a>
<a href="#discarded">Pairs discarded</a>
<a href="#distance">Contact distance</a>
<a href="#orientation">Read orientation</a>
<a href="#fulltable">Full table</a>
</nav>
<div id="qc-content">
)HTML";

    // --- Overview --------------------------------------------------------
    out << "<section id=\"overview\"><h2>Overview</h2><p class=\"qc-sub\">" << n << " sample"
        << (n == 1 ? "" : "s") << ", five metrics judged against a cutoff (see the note at the end "
                                   "of this report).</p>\n<div class=\"qc-samples\">\n";
    for (std::size_t row = 0; row < n; ++row) {
        const QcJudged& j = judged[row];
        const std::vector<QcMetricRow> metrics = {
            {"Valid pairs", j.valid_pct,
             "cutoff: >= " + qc_percent(thresholds.min_valid_pairs_percent) + " of sequenced reads",
             j.valid_status},
            {"Inter-chromosomal contacts", j.inter_chrom_pct,
             "cutoff: <= " + qc_percent(thresholds.max_inter_chromosomal_percent) + " of valid pairs",
             j.inter_chrom_status},
            {"Duplicated pairs", j.duplicated_pct,
             "cutoff: <= " + qc_percent(thresholds.max_duplicated_pairs_percent) + " of mappable pairs",
             j.duplicated_status},
            {"Dangling ends", j.dangling_pct,
             "cutoff: <= " + qc_percent(thresholds.max_dangling_end_percent) + " of mappable pairs",
             j.dangling_status},
            {"Self-ligation", j.self_ligation_pct,
             "cutoff: <= " + qc_percent(thresholds.max_self_ligation_percent) + " of mappable pairs",
             j.self_ligation_status},
        };
        int failed = 0;
        int warned = 0;
        for (const QcMetricRow& metric : metrics) {
            if (metric.status == QcStatus::Fail) {
                ++failed;
            } else if (metric.status == QcStatus::Warn) {
                ++warned;
            }
        }
        const char* overall = failed > 0 ? "has metrics that failed their cutoff"
                             : warned > 0 ? "has metrics close to their cutoff"
                                          : "passed every judged metric";
        out << "<div class=\"qc-card\"><h3>" << qc_html_escape(frame.index[row]) << "</h3>"
            << "<p class=\"qc-overall\">" << overall << "</p>\n";
        for (const QcMetricRow& metric : metrics) {
            qc_write_metric_row(out, metric);
        }
        out << "</div>\n";
    }
    out << "</div></section>\n";

    // --- Pairs sequenced ---------------------------------------------------
    out << "<section id=\"sequenced\"><h2>Pairs sequenced</h2>"
           "<p class=\"qc-sub\">Sequenced reads, pairs mappable/unique/high quality, and the valid "
           "Hi-C contacts kept from them.</p>\n";
    {
        std::vector<QcBarRow> rows;
        for (std::size_t row = 0; row < n; ++row) {
            const double sequenced = frame.has("Sequenced reads") ? frame.number("Sequenced reads", row) : 0.0;
            const double mappable = frame.has("Pairs mappable, unique and high quality")
                                        ? frame.number("Pairs mappable, unique and high quality", row)
                                        : 0.0;
            const double valid = frame.has("Hi-C contacts") ? frame.number("Hi-C contacts", row) : 0.0;
            rows.push_back({frame.index[row] + ": sequenced reads", sequenced, sequenced, QcStatus::Unknown});
            rows.push_back(
                {frame.index[row] + ": mappable, unique, high quality", mappable, sequenced, QcStatus::Unknown});
            rows.push_back(
                {frame.index[row] + ": Hi-C contacts (valid pairs)", valid, sequenced, QcStatus::Unknown});
        }
        qc_write_bar_section(out, rows);
    }
    out << "</section>\n";

    // --- Unmappable and non-unique pairs ----------------------------------
    out << "<section id=\"unmappable\"><h2>Unmappable and non-unique pairs</h2>"
           "<p class=\"qc-sub\">Fraction of sequenced reads, per sample, that mapped uniquely and at "
           "high quality (Hi-C contacts) against those that did not.</p>\n";
    {
        std::vector<QcBarRow> rows;
        std::vector<std::string> row_labels;
        std::vector<QcTableColumn> columns = {{"Hi-C contacts", {}}, {"Hi-C contacts %", {}},
                                              {"Low mapping quality", {}}, {"Low mapping quality %", {}},
                                              {"One mate not unique", {}}, {"One mate not unique %", {}},
                                              {"One mate unmapped", {}}, {"One mate unmapped %", {}}};
        for (std::size_t row = 0; row < n; ++row) {
            const double sequenced = frame.has("Sequenced reads") ? frame.number("Sequenced reads", row) : 0.0;
            row_labels.push_back(frame.index[row]);
            for (const char* name :
                {"Hi-C contacts", "Low mapping quality", "One mate not unique", "One mate unmapped"}) {
                const double value = frame.has(name) ? frame.number(name, row) : 0.0;
                rows.push_back({frame.index[row] + ": " + name, value, sequenced, QcStatus::Unknown});
            }
        }
        std::size_t column_index = 0;
        for (const char* name :
            {"Hi-C contacts", "Low mapping quality", "One mate not unique", "One mate unmapped"}) {
            for (std::size_t row = 0; row < n; ++row) {
                const double sequenced =
                    frame.has("Sequenced reads") ? frame.number("Sequenced reads", row) : 0.0;
                const double value = frame.has(name) ? frame.number(name, row) : 0.0;
                columns[column_index].values.push_back(qc_format_count(value));
                columns[column_index + 1].values.push_back(
                    sequenced > 0 ? qc_percent(100.0 * value / sequenced) : "n/a");
            }
            column_index += 2;
        }
        qc_write_bar_section(out, rows);
        qc_write_table(out, "unmappable-table", row_labels, columns);
        out << "<a class=\"qc-download\" href=\"unmapable_table.txt\">Download table (tab separated)"
               "</a>\n";
    }
    out << "</section>\n";

    // --- Pairs discarded ---------------------------------------------------
    out << "<section id=\"discarded\"><h2>Pairs discarded</h2>"
           "<p class=\"qc-sub\">Fraction of mappable, unique, high quality pairs discarded when "
           "building the Hi-C matrix, per sample.</p>\n"
           "<dl class=\"qc-defs\">"
           "<dt>Dangling ends</dt><dd>Reads that start with the restriction site: digested but not "
           "ligated.</dd>"
           "<dt>Same fragment</dt><dd>Inward-facing mates up to 800 bp apart with no restriction "
           "site between them; not a valid Hi-C pair.</dd>"
           "<dt>Self circle</dt><dd>Pairs within 25 kb with outward read orientation.</dd>"
           "<dt>Self ligation</dt><dd>Pairs with a restriction site between them, within 800 bp."
           "</dd></dl>\n";
    {
        std::vector<QcBarRow> rows;
        std::vector<std::string> row_labels;
        static const char* const prefixes[] = {"One mate not close to rest site", "dangling end",
                                                "duplicated pairs", "same fragment", "self circle",
                                                "self ligation (removed)"};
        // Columns are built after the per-sample pass below, since a prefix
        // such as "dangling end" can match more than one real column (one
        // per restriction enzyme) and every sample must get every column
        // that *any* sample has, in a fixed order.
        std::vector<std::string> discovered_columns;
        std::vector<std::unordered_map<std::string, double>> per_row_values(n);
        for (std::size_t row = 0; row < n; ++row) {
            const double mappable = frame.has("Pairs mappable, unique and high quality")
                                        ? frame.number("Pairs mappable, unique and high quality", row)
                                        : 0.0;
            row_labels.push_back(frame.index[row]);
            for (const char* prefix : prefixes) {
                for (const auto& [column, value] : qc_values_by_prefix_all(frame, prefix, row)) {
                    per_row_values[row][column] = value;
                    if (std::find(discovered_columns.begin(), discovered_columns.end(), column) ==
                        discovered_columns.end()) {
                        discovered_columns.push_back(column);
                    }
                    QcStatus status = QcStatus::Unknown;
                    if (column == "duplicated pairs") {
                        status = judged[row].duplicated_status;
                    } else if (column.find("dangling end") != std::string::npos) {
                        status = judged[row].dangling_status;
                    } else if (column == "self ligation (removed)") {
                        status = judged[row].self_ligation_status;
                    }
                    rows.push_back({frame.index[row] + ": " + column, value, mappable, status});
                }
            }
        }
        std::vector<QcTableColumn> columns;
        for (const std::string& column : discovered_columns) {
            columns.push_back({column, {}});
            columns.push_back({column + " %", {}});
        }
        for (std::size_t row = 0; row < n; ++row) {
            const double mappable = frame.has("Pairs mappable, unique and high quality")
                                        ? frame.number("Pairs mappable, unique and high quality", row)
                                        : 0.0;
            for (std::size_t c = 0; c < discovered_columns.size(); ++c) {
                const auto it = per_row_values[row].find(discovered_columns[c]);
                const double value = it == per_row_values[row].end() ? 0.0 : it->second;
                columns[2 * c].values.push_back(qc_format_count(value));
                columns[2 * c + 1].values.push_back(mappable > 0 ? qc_percent(100.0 * value / mappable)
                                                                 : "n/a");
            }
        }
        qc_write_bar_section(out, rows);
        qc_write_table(out, "discarded-table", row_labels, columns);
        out << "<a class=\"qc-download\" href=\"discarded_table.txt\">Download table (tab separated)"
               "</a>\n";
    }
    out << "</section>\n";

    // --- Contact distance ---------------------------------------------------
    out << "<section id=\"distance\"><h2>Contact distance</h2>"
           "<p class=\"qc-sub\">Fraction of valid Hi-C contacts, per sample, by genomic distance. "
           "HiCExplorer's own usage guidance: good Hi-C libraries have less than 10% "
           "inter-chromosomal contacts.</p>\n";
    {
        std::vector<QcBarRow> rows;
        std::vector<std::string> row_labels;
        static const char* const names[] = {"inter chromosomal", "Intra short range (< 20kb)",
                                            "Intra long range (>= 20kb)"};
        std::vector<QcTableColumn> columns;
        for (const char* name : names) {
            columns.push_back({name, {}});
            columns.push_back({std::string(name) + " %", {}});
        }
        for (std::size_t row = 0; row < n; ++row) {
            const double valid = frame.has("Hi-C contacts") ? frame.number("Hi-C contacts", row) : 0.0;
            row_labels.push_back(frame.index[row]);
            std::size_t column_index = 0;
            for (const char* name : names) {
                const double value = frame.has(name) ? frame.number(name, row) : 0.0;
                const QcStatus status =
                    (std::string(name) == "inter chromosomal") ? judged[row].inter_chrom_status
                                                                : QcStatus::Unknown;
                rows.push_back({frame.index[row] + ": " + name, value, valid, status});
                columns[column_index].values.push_back(qc_format_count(value));
                columns[column_index + 1].values.push_back(valid > 0 ? qc_percent(100.0 * value / valid)
                                                                     : "n/a");
                column_index += 2;
            }
        }
        qc_write_bar_section(out, rows);
        qc_write_table(out, "distance-table", row_labels, columns);
        out << "<a class=\"qc-download\" href=\"distance_table.txt\">Download table (tab separated)"
               "</a>\n";
    }
    out << "</section>\n";

    // --- Read orientation ---------------------------------------------------
    out << "<section id=\"orientation\"><h2>Read orientation</h2>"
           "<p class=\"qc-sub\">Fraction of the four read-pair orientations, per sample, informative "
           "for ligation bias.</p>\n";
    {
        std::vector<QcBarRow> rows;
        std::vector<std::string> row_labels;
        static const char* const names[] = {"Read pair type: inward pairs", "Read pair type: outward pairs",
                                            "Read pair type: left pairs", "Read pair type: right pairs"};
        std::vector<QcTableColumn> columns;
        for (const char* name : names) {
            columns.push_back({name, {}});
            columns.push_back({std::string(name) + " %", {}});
        }
        for (std::size_t row = 0; row < n; ++row) {
            double total = 0.0;
            for (const char* name : names) {
                total += frame.has(name) ? frame.number(name, row) : 0.0;
            }
            row_labels.push_back(frame.index[row]);
            std::size_t column_index = 0;
            for (const char* name : names) {
                const double value = frame.has(name) ? frame.number(name, row) : 0.0;
                rows.push_back({frame.index[row] + ": " + name, value, total, QcStatus::Unknown});
                columns[column_index].values.push_back(qc_format_count(value));
                columns[column_index + 1].values.push_back(total > 0 ? qc_percent(100.0 * value / total)
                                                                     : "n/a");
                column_index += 2;
            }
        }
        qc_write_bar_section(out, rows);
        qc_write_table(out, "orientation-table", row_labels, columns);
        out << "<a class=\"qc-download\" href=\"read_orientation_table.txt\">Download table (tab "
               "separated)</a>\n";
    }
    out << "</section>\n";

    // --- Full table ---------------------------------------------------------
    out << "<section id=\"fulltable\"><h2>Full table</h2>"
           "<p class=\"qc-sub\">Every field QC.log reports, one column per sample. Click a header to "
           "sort.</p>\n";
    {
        std::vector<QcTableColumn> columns;
        for (const std::string& name : frame.columns) {
            QcTableColumn column{name, {}};
            for (std::size_t row = 0; row < n; ++row) {
                column.values.push_back(frame.column(name)[row].str());
            }
            columns.push_back(std::move(column));
        }
        qc_write_table(out, "full-table", frame.index, columns);
        out << "<a class=\"qc-download\" href=\"QC_table.txt\">Download table (tab separated)</a>\n";
    }
    out << "</section>\n";

    out << R"HTML(</div>
<p class="qc-note" style="margin-left:258px;max-width:900px;">Valid pairs and inter-chromosomal
cutoffs follow HiCExplorer's own usage guidance. The duplicated pairs, dangling ends and
self-ligation cutoffs have no published HiCExplorer standard behind them; they are conservative
starting points, and every cutoff is a hicPrepareQCreport command line option
(--minValidPairsPercent, --maxInterChromosomalPercent, --maxDuplicatedPairsPercent,
--maxDanglingEndPercent, --maxSelfLigationPercent). A metric within )HTML"
        << static_cast<int>(margin * 100 + 0.5)
        << R"HTML(% of the wrong side of its cutoff is WARN rather than an outright PASS or FAIL.</p>
<script>
(function () {
  var sections = Array.prototype.slice.call(document.querySelectorAll("section[id]"));
  var links = Array.prototype.slice.call(document.querySelectorAll("#qc-nav a"));
  function setActive() {
    var top = window.scrollY + 80;
    var current = sections[0];
    sections.forEach(function (s) { if (s.offsetTop <= top) current = s; });
    links.forEach(function (a) {
      a.classList.toggle("qc-active", a.getAttribute("href") === "#" + current.id);
    });
  }
  window.addEventListener("scroll", setActive, { passive: true });
  setActive();

  document.querySelectorAll("table.qc-table").forEach(function (table) {
    var tbody = table.querySelector("tbody");
    table.querySelectorAll("thead th").forEach(function (th, index) {
      th.addEventListener("click", function () {
        var ascending = th.getAttribute("data-asc") !== "true";
        table.querySelectorAll("thead th").forEach(function (h) { h.removeAttribute("data-asc"); });
        th.setAttribute("data-asc", ascending ? "true" : "false");
        var rows = Array.prototype.slice.call(tbody.querySelectorAll("tr"));
        rows.sort(function (a, b) {
          var av = a.children[index].textContent.replace(/,/g, "");
          var bv = b.children[index].textContent.replace(/,/g, "");
          var an = parseFloat(av), bn = parseFloat(bv);
          var cmp = (!isNaN(an) && !isNaN(bn)) ? an - bn : av.localeCompare(bv);
          return ascending ? cmp : -cmp;
        });
        rows.forEach(function (row) { tbody.appendChild(row); });
      });
    });
  });
})();
</script>
</body>
</html>
)HTML";

    if (!out) {
        throw std::runtime_error("could not write " + folder + "/hicQC_report.html");
    }
}

}  // namespace hicx
