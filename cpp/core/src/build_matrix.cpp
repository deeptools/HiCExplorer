#include "hicx/build_matrix.hpp"

#include "hicx/argparse.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/version.hpp"

#include "qc_template_html.inc"  // kQcTemplateHtml, generated from hicexplorer/qc_template.html

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <limits>
#include <optional>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <random>
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

// A cell of hicPrepareQCreport's DataFrames: int(fields[1]) where Python's
// int() accepts the text, the text otherwise; a fraction of the derived
// tables is a float64.
struct QcCell {
    bool is_int = false;
    bool is_float = false;
    std::int64_t integer = 0;
    double real = 0.0;
    std::string text;

    // As DataFrame.to_csv writes it.
    [[nodiscard]] std::string str() const {
        if (is_int) {
            return std::to_string(integer);
        }
        if (is_float) {
            return python_float_repr(real);
        }
        return text;
    }
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

// One output table: its columns, each a cell per output row.
struct QcOutputTable {
    std::vector<std::string> columns;
    std::vector<std::vector<QcCell>> values;  // per column, per output row

    void add(std::string name, std::vector<QcCell> column) {
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
            out << '\t' << csv_field(column[row].str());
        }
        out << '\n';
    }
    if (!out) {
        throw std::runtime_error("could not write " + path);
    }
}

// The counts of the left rows.
std::vector<QcCell> count_column(const QcFrame& frame, const std::string& name,
                                 const std::vector<std::pair<std::size_t, std::size_t>>& rows) {
    std::vector<QcCell> out;
    const std::vector<QcCell>& cells = frame.column(name);
    for (const auto& row : rows) {
        out.push_back(cells[row.first]);
    }
    return out;
}

// The fractions of the right rows, each over its own row's denominator.
std::vector<QcCell> ratio_column(const QcFrame& frame, const std::string& name,
                                 const std::vector<double>& denominators,
                                 const std::vector<std::pair<std::size_t, std::size_t>>& rows) {
    std::vector<QcCell> out;
    for (const auto& row : rows) {
        QcCell cell;
        cell.is_float = true;
        cell.real = frame.number(name, row.second) / denominators[row.second];
        out.push_back(std::move(cell));
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

// The five tables of hicPrepareQCreport.main: QC_table.txt is the frame
// itself, the four others are what make_figure_* return (the counts joined
// with the fractions on the index). `rows` holds each table's (left, right)
// frame rows.
struct QcTables {
    QcOutputTable all, unmap, discarded, distance, orientation;
    std::vector<std::pair<std::size_t, std::size_t>> plain, joined;
};

QcTables build_qc_tables(const QcFrame& frame) {
    QcTables tables;
    const auto plain = plain_rows(frame);
    const auto joined = joined_rows(frame);
    tables.plain = plain;
    tables.joined = joined;

    // QC_table.txt: table.to_csv, every column.
    {
        QcOutputTable table;
        for (const std::string& column : frame.columns) {
            table.add(column, count_column(frame, column, plain));
        }
        tables.all = std::move(table);
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
        tables.unmap = std::move(table);
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
        tables.discarded = std::move(table);
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
        tables.distance = std::move(table);
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
        tables.orientation = std::move(table);
    }
    return tables;
}

// ---------------------------------------------------------------------------
// hicQC.html (save_html) and the DataFrames of the five charts (make_figure_*).

// f"{x:.<precision>f}": Python prints nan and the infinities as words.
std::string qc_fixed(double value, int precision) {
    if (std::isnan(value)) {
        return "nan";
    }
    if (std::isinf(value)) {
        return value > 0 ? "inf" : "-inf";
    }
    char buffer[400];
    std::snprintf(buffer, sizeof(buffer), "%.*f", precision, value);
    return buffer;
}

// '{:.2%}'.format(x)
std::string qc_format_percent(double value) {
    if (std::isnan(value)) {
        return "nan%";
    }
    return qc_fixed(value * 100.0, 2) + "%";
}

// The thousands separator of '{:,}' on the leading integer digits of `text`.
std::string qc_group_thousands(const std::string& text) {
    const std::size_t start = (!text.empty() && text[0] == '-') ? 1 : 0;
    std::size_t end = start;
    while (end < text.size() && std::isdigit(static_cast<unsigned char>(text[end]))) {
        ++end;
    }
    const std::string digits = text.substr(start, end - start);
    std::string grouped;
    for (std::size_t i = 0; i < digits.size(); ++i) {
        if (i > 0 && (digits.size() - i) % 3 == 0) {
            grouped += ',';
        }
        grouped += digits[i];
    }
    return text.substr(0, start) + grouped + text.substr(end);
}

// The formatter save_html gives the four derived tables:
// lambda x: '{:,}'.format(x) if x > 1 else '{:.2%}'.format(x)
// ('{:,}' of a float is its repr with grouped integer digits, unless the repr
// is in exponent notation or a word).
std::string qc_derived_formatter(const QcCell& cell) {
    if (cell.is_int) {
        return cell.integer > 1 ? qc_group_thousands(std::to_string(cell.integer))
                                : qc_format_percent(static_cast<double>(cell.integer));
    }
    if (cell.is_float) {
        if (!(cell.real > 1)) {
            return qc_format_percent(cell.real);
        }
        const std::string repr = npy::float_repr(cell.real);
        if (repr.find('e') != std::string::npos || repr.find("inf") != std::string::npos) {
            return repr;
        }
        return qc_group_thousands(repr);
    }
    throw std::runtime_error("TypeError: '>' not supported between instances of 'str' and 'int' "
                             "(a QC table column holds text)");
}

// Styler's default formatter: a column of the QC frame holds Python ints and
// strings, which it shows as str() does.
std::string qc_default_formatter(const QcCell& cell) { return cell.str(); }

// Five random hexadecimal characters, as Styler's uuid4().hex[:5].
std::string qc_styler_uuid() {
    static std::mt19937_64 engine{std::random_device{}()};
    static const char* const hex = "0123456789abcdef";
    std::string id;
    for (int i = 0; i < 5; ++i) {
        id += hex[engine() % 16];
    }
    return id;
}

// DataFrame.style.format(formatter).to_html() in pandas 2.2 (html_style.tpl
// and html_table.tpl) for a table without styles, caption or attributes, a
// named single level index and unnamed columns. Styler escapes nothing by
// default.
std::string qc_styler_to_html(const QcFrame& frame, const QcOutputTable& table,
                              const std::vector<std::pair<std::size_t, std::size_t>>& rows,
                              std::string (*format)(const QcCell&)) {
    const std::string uuid = qc_styler_uuid();
    std::string out = "<style type=\"text/css\">\n</style>\n<table id=\"T_" + uuid + "\">\n";
    out += "  <thead>\n    <tr>\n      <th class=\"blank level0\" >&nbsp;</th>\n";
    for (std::size_t c = 0; c < table.columns.size(); ++c) {
        const std::string col = std::to_string(c);
        out += "      <th id=\"T_" + uuid + "_level0_col" + col + "\" class=\"col_heading level0 col" +
               col + "\" >" + table.columns[c] + "</th>\n";
    }
    out += "    </tr>\n";
    out += "    <tr>\n      <th class=\"index_name level0\" >" + frame.index_name + "</th>\n";
    for (std::size_t c = 0; c < table.columns.size(); ++c) {
        out += "      <th class=\"blank col" + std::to_string(c) + "\" >&nbsp;</th>\n";
    }
    out += "    </tr>\n  </thead>\n  <tbody>\n";
    for (std::size_t r = 0; r < rows.size(); ++r) {
        const std::string row = std::to_string(r);
        out += "    <tr>\n      <th id=\"T_" + uuid + "_level0_row" + row +
               "\" class=\"row_heading level0 row" + row + "\" >" + frame.index[rows[r].first] +
               "</th>\n";
        for (std::size_t c = 0; c < table.columns.size(); ++c) {
            const std::string col = std::to_string(c);
            out += "      <td id=\"T_" + uuid + "_row" + row + "_col" + col + "\" class=\"data row" +
                   row + " col" + col + "\" >" + format(table.values[c][r]) + "</td>\n";
        }
        out += "    </tr>\n";
    }
    out += "  </tbody>\n</table>\n";
    return out;
}

// str.replace: every occurrence, left to right.
void qc_replace_all(std::string& text, const std::string& pattern, const std::string& replacement) {
    std::string out;
    std::size_t from = 0;
    for (std::size_t at = text.find(pattern); at != std::string::npos;
         at = text.find(pattern, from)) {
        out.append(text, from, at - from);
        out += replacement;
        from = at + pattern.size();
    }
    out.append(text, from, std::string::npos);
    text = std::move(out);
}

// save_html. The template is hicexplorer/qc_template.html, compiled in.
void qc_save_html(const std::string& filename, const QcFrame& frame, const QcTables& tables) {
    std::string html(reinterpret_cast<const char*>(kQcTemplateHtml), sizeof(kQcTemplateHtml));
    qc_replace_all(html, "%%TABLE_UNMAP%%",
                   qc_styler_to_html(frame, tables.unmap, tables.joined, &qc_derived_formatter));
    qc_replace_all(html, "%%TABLE_DISCARDED%%",
                   qc_styler_to_html(frame, tables.discarded, tables.joined, &qc_derived_formatter));
    qc_replace_all(html, "%%TABLE_DISTANCE%%",
                   qc_styler_to_html(frame, tables.distance, tables.joined, &qc_derived_formatter));
    qc_replace_all(html, "%%TABLE_ORIENTATION%%",
                   qc_styler_to_html(frame, tables.orientation, tables.joined,
                                     &qc_derived_formatter));
    QcOutputTable all = tables.all;
    for (const char* dropped : {"Min rest. site distance", "Max library insert size"}) {
        const auto it = std::find(all.columns.begin(), all.columns.end(), dropped);
        if (it != all.columns.end()) {
            all.values.erase(all.values.begin() + (it - all.columns.begin()));
            all.columns.erase(it);
        }
    }
    qc_replace_all(html, "%%TABLE%%",
                   qc_styler_to_html(frame, all, tables.plain, &qc_default_formatter));
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        throw std::runtime_error("could not write " + filename);
    }
    out << html;
    if (!out) {
        throw std::runtime_error("could not write " + filename);
    }
}

// A chart's DataFrame as JSON: {"index", "index_name", "columns",
// "columns_name", "values" (row major)}. The labels are str() of the index
// values, which are the --labels or the File entries.
std::string qc_chart_frame_json(const std::vector<std::string>& index,
                                const std::optional<std::string>& index_name,
                                const std::vector<std::string>& columns,
                                const std::optional<std::string>& columns_name,
                                const std::vector<double>& values) {
    auto optional_string = [](const std::optional<std::string>& value) {
        return value.has_value() ? plot::json_string(*value) : std::string("null");
    };
    plot::JsonObject object;
    object.add("index", plot::json_strings(index));
    object.add("index_name", optional_string(index_name));
    object.add("columns", plot::json_strings(columns));
    object.add("columns_name", optional_string(columns_name));
    object.add("values", plot::json_numbers(values));
    return object.str();
}

// table[names].T / denominators: a row per name, a column per sample. The
// columns of the transposed frame are the index of the denominator Series,
// so pandas divides position by position, also with repeated names.
std::string qc_transposed_ratio_json(const QcFrame& frame, const std::vector<std::string>& names,
                                     const std::vector<double>& denominators) {
    std::vector<double> values;
    for (const std::string& name : names) {
        const std::vector<double> column = numbers(frame, name);
        for (std::size_t s = 0; s < column.size(); ++s) {
            values.push_back(column[s] / denominators[s]);
        }
    }
    return qc_chart_frame_json(names, std::nullopt, frame.index, frame.index_name, values);
}

// make_figure_*: the DataFrame each chart plots.
std::string qc_chart_json(const QcFrame& frame) {
    const std::size_t samples = frame.index.size();
    plot::JsonObject charts;
    {
        // make_figure_pairs_used: table[[...]] / 1e6.
        const std::vector<std::string> names = {"Hi-C contacts",
                                                "Pairs mappable, unique and high quality",
                                                "Sequenced reads"};
        std::vector<std::vector<double>> columns;
        for (const std::string& name : names) {
            columns.push_back(numbers(frame, name));
        }
        std::vector<double> values;
        for (std::size_t s = 0; s < samples; ++s) {
            for (const std::vector<double>& column : columns) {
                values.push_back(column[s] / 1e6);
            }
        }
        charts.add("pairs_sequenced", qc_chart_frame_json(frame.index, frame.index_name, names,
                                                          std::nullopt, values));
    }
    charts.add("unmappable_and_non_unique",
               qc_transposed_ratio_json(frame,
                                        {"Hi-C contacts", "Low mapping quality",
                                         "One mate not unique", "One mate unmapped"},
                                        numbers(frame, "Sequenced reads")));
    {
        // make_figure_pairs_discarded: an exact column, or every column that
        // contains the prefix.
        std::vector<std::string> names;
        for (const char* prefix : {"One mate not close to rest site", "dangling end",
                                   "duplicated pairs", "same fragment", "self circle",
                                   "self ligation (removed)"}) {
            const std::string name(prefix);
            if (frame.has(name)) {
                names.push_back(name);
                continue;
            }
            for (const std::string& column : frame.columns) {
                if (column.find(name) != std::string::npos) {
                    names.push_back(column);
                }
            }
        }
        charts.add("pairs_discarded",
                   qc_transposed_ratio_json(frame, names,
                                            numbers(frame, "Pairs mappable, unique and high quality")));
    }
    charts.add("distance",
               qc_transposed_ratio_json(frame,
                                        {"inter chromosomal", "Intra short range (< 20kb)",
                                         "Intra long range (>= 20kb)"},
                                        numbers(frame, "Hi-C contacts")));
    {
        // make_figure_read_orientation: _t / _t.sum(axis=0), the integer sum
        // of the four orientations per sample.
        const std::vector<std::string> names = {
            "Read pair type: inward pairs", "Read pair type: outward pairs",
            "Read pair type: left pairs", "Read pair type: right pairs"};
        std::vector<double> totals(samples, 0.0);
        for (const std::string& name : names) {
            const std::vector<double> values = numbers(frame, name);
            for (std::size_t s = 0; s < samples; ++s) {
                totals[s] += values[s];
            }
        }
        charts.add("read_orientation", qc_transposed_ratio_json(frame, names, totals));
    }
    return charts.str();
}

}  // namespace

void write_qc_report_tables(const std::string& folder, const std::vector<std::string>& log_texts,
                            const std::optional<std::vector<std::string>>& labels) {
    const QcFrame frame = build_qc_frame(log_texts, labels);
    const QcTables tables = build_qc_tables(frame);
    write_qc_output(folder + "/QC_table.txt", frame, tables.all, tables.plain);
    write_qc_output(folder + "/unmapable_table.txt", frame, tables.unmap, tables.joined);
    write_qc_output(folder + "/discarded_table.txt", frame, tables.discarded, tables.joined);
    write_qc_output(folder + "/distance_table.txt", frame, tables.distance, tables.joined);
    write_qc_output(folder + "/read_orientation_table.txt", frame, tables.orientation,
                    tables.joined);
}

std::string write_qc_report_html(const std::string& folder,
                                 const std::vector<std::string>& log_texts,
                                 const std::optional<std::vector<std::string>>& labels) {
    const QcFrame frame = build_qc_frame(log_texts, labels);
    const QcTables tables = build_qc_tables(frame);
    std::string charts = qc_chart_json(frame);
    qc_save_html(folder + "/hicQC.html", frame, tables);
    return charts;
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

// PASS when `value` is at least `minimum`; WARN when it is below `minimum`
// but at least minimum * (1 - margin); FAIL below that. Used for metrics a
// report wants to see enough of (valid pairs).
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

// PASS when `value` is at most `maximum`; WARN when it is above `maximum`
// but at most maximum * (1 + margin); FAIL above that. Used for
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

// Section status: the worst verdict of the metrics a section judges. A section
// that judges none stays Unknown and is shown as "not judged".
QcStatus qc_worst(QcStatus a, QcStatus b) {
    auto rank = [](QcStatus s) {
        switch (s) {
            case QcStatus::Fail:
                return 3;
            case QcStatus::Warn:
                return 2;
            case QcStatus::Pass:
                return 1;
            default:
                return 0;
        }
    };
    return rank(a) >= rank(b) ? a : b;
}

// The lower-case word of a verdict, used as the CSS class of the sidebar entry
// and the card title, and as the status key of a chart bar.
const char* qc_status_word(QcStatus status) {
    switch (status) {
        case QcStatus::Pass:
            return "pass";
        case QcStatus::Warn:
            return "warn";
        case QcStatus::Fail:
            return "fail";
        default:
            return "info";
    }
}

// A JSON string literal that is also safe inside an inline <script>: '<' is
// escaped so a sample label can never close the script element.
std::string qc_json_string(const std::string& text) {
    std::string out = "\"";
    for (const char c : text) {
        switch (c) {
            case '"':
                out += "\\\"";
                break;
            case '\\':
                out += "\\\\";
                break;
            case '\n':
                out += "\\n";
                break;
            case '\r':
                out += "\\r";
                break;
            case '\t':
                out += "\\t";
                break;
            case '<':
                out += "\\u003c";
                break;
            case '>':
                out += "\\u003e";
                break;
            case '&':
                out += "\\u0026";
                break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buffer[8];
                    std::snprintf(buffer, sizeof(buffer), "\\u%04x", static_cast<unsigned char>(c));
                    out += buffer;
                } else {
                    out += c;
                }
        }
    }
    out += "\"";
    return out;
}

std::string qc_json_number(std::optional<double> value) {
    if (!value.has_value() || !std::isfinite(*value)) {
        return "null";
    }
    char buffer[40];
    std::snprintf(buffer, sizeof(buffer), "%.10g", *value);
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

// One series of a grouped bar chart: one bar per sample. `y` is what the bar
// height shows (a count or a percentage, nullopt draws no bar), `text` the
// hover text ("1,234 (5.6%)"), and `status` the verdict of a judged metric
// (Unknown for a metric without a cutoff, drawn in a palette colour).
struct QcChartSeries {
    std::string name;
    std::vector<std::optional<double>> y;
    std::vector<std::string> text;
    std::vector<QcStatus> status;
};

struct QcChart {
    std::string id;
    std::string y_title;
    bool percent = true;
    std::vector<QcChartSeries> series;
};

// The hover text of one bar: the count, its share of `denominator` when there
// is one, and the verdict when the metric is judged.
std::string qc_bar_text(double count, std::optional<double> denominator, QcStatus status) {
    std::string text = qc_format_count(count);
    if (denominator.has_value() && *denominator > 0) {
        text += " (" + qc_percent(100.0 * count / *denominator) + ")";
    }
    if (status != QcStatus::Unknown) {
        text += std::string(" ") + qc_status_label(status);
    }
    return text;
}

std::optional<double> qc_share(double count, std::optional<double> denominator) {
    if (denominator.has_value() && *denominator > 0) {
        return 100.0 * count / *denominator;
    }
    return std::nullopt;
}

std::string qc_chart_json(const QcChart& chart, const std::vector<std::string>& samples) {
    std::ostringstream out;
    out << "{\"id\":" << qc_json_string(chart.id) << ",\"ytitle\":" << qc_json_string(chart.y_title)
        << ",\"percent\":" << (chart.percent ? "true" : "false") << ",\"samples\":[";
    for (std::size_t i = 0; i < samples.size(); ++i) {
        out << (i ? "," : "") << qc_json_string(samples[i]);
    }
    out << "],\"series\":[";
    for (std::size_t s = 0; s < chart.series.size(); ++s) {
        const QcChartSeries& series = chart.series[s];
        out << (s ? "," : "") << "{\"name\":" << qc_json_string(series.name) << ",\"y\":[";
        for (std::size_t i = 0; i < series.y.size(); ++i) {
            out << (i ? "," : "") << qc_json_number(series.y[i]);
        }
        out << "],\"text\":[";
        for (std::size_t i = 0; i < series.text.size(); ++i) {
            out << (i ? "," : "") << qc_json_string(series.text[i]);
        }
        out << "],\"status\":[";
        for (std::size_t i = 0; i < series.status.size(); ++i) {
            out << (i ? "," : "")
                << (series.status[i] == QcStatus::Unknown
                        ? std::string("null")
                        : qc_json_string(qc_status_word(series.status[i])));
        }
        out << "]}";
    }
    out << "]}";
    return out.str();
}

// A table column: one formatted value per row and, for a judged metric, the
// verdict of each row, shown as a badge next to the value. `judged_metric`
// names that metric in a visually hidden label before each badge, so a screen
// reader (and the harness validator qc_pretty_report) can tell which metric
// the badge belongs to.
struct QcTableColumn {
    std::string header;
    std::vector<std::string> values;  // one per row, already formatted
    std::vector<QcStatus> statuses = {};  // empty, or one per row
    std::string judged_metric = {};
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
            out << "<td>" << qc_html_escape(column.values[row]);
            if (row < column.statuses.size() && column.statuses[row] != QcStatus::Unknown) {
                out << " <span class=\"qc-sr\">" << qc_html_escape(column.judged_metric) << ":</span>"
                    << "<span class=\"qc-badge qc-badge-sm " << qc_status_class(column.statuses[row])
                    << "\">" << qc_status_label(column.statuses[row]) << "</span>";
            }
            out << "</td>";
        }
        out << "</tr>\n";
    }
    out << "</tbody></table></div>\n";
}

void qc_write_download(std::ostream& out, const char* file) {
    out << "<a class=\"qc-download\" href=\"" << file << "\">Download table (tab separated)</a>\n";
}

// One collapsible card of the report: its body is written later into `body`.
struct QcSection {
    std::string id;
    std::string title;
    QcStatus status = QcStatus::Unknown;
    std::string description;  // HTML, written before the documentation link
    std::string help_page;    // page under the HiCExplorer documentation
    std::ostringstream body;
};

constexpr const char* kQcDocsBase = "https://hicexplorer.readthedocs.io/en/latest/content/tools/";

std::string qc_status_title(QcStatus status) {
    return status == QcStatus::Unknown ? "not judged" : qc_status_word(status);
}

std::string qc_now_iso() {
    const std::time_t now = std::time(nullptr);
    std::tm local{};
    localtime_r(&now, &local);
    char buffer[32];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M", &local);
    return buffer;
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

    auto count_of = [&frame](const std::string& name, std::size_t row) {
        return frame.has(name) ? frame.number(name, row) : 0.0;
    };
    const int margin_percent = static_cast<int>(margin * 100 + 0.5);

    std::vector<QcChart> charts;
    QcSection overview, sequenced_section, unmappable, discarded, distance, orientation, full;

    // --- Overview --------------------------------------------------------
    overview.id = "overview";
    overview.title = "Overview";
    overview.help_page = "hicPrepareQCreport.html";
    overview.description =
        std::to_string(n) + " sample" + (n == 1 ? "" : "s") +
        ", five metrics per sample, each compared with a cutoff. A value on the right side of its "
        "cutoff, or exactly on it, is PASS; a value past the cutoff by at most " +
        std::to_string(margin_percent) + "% of the cutoff is WARN; a value further past it is FAIL.";
    overview.body << "<div class=\"qc-samples\">\n";
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
        QcStatus sample_status = QcStatus::Unknown;
        for (const QcMetricRow& metric : metrics) {
            if (metric.status == QcStatus::Fail) {
                ++failed;
            } else if (metric.status == QcStatus::Warn) {
                ++warned;
            }
            sample_status = qc_worst(sample_status, metric.status);
        }
        overview.status = qc_worst(overview.status, sample_status);
        const char* overall = failed > 0 ? "has metrics past their cutoff by more than the margin"
                             : warned > 0 ? "has metrics past their cutoff, within the margin"
                                          : "passed every judged metric";
        overview.body << "<div class=\"qc-card\"><h3 class=\"" << qc_status_word(sample_status) << "\">"
                      << qc_html_escape(frame.index[row]) << "</h3>"
                      << "<p class=\"qc-overall\">" << overall << "</p>\n";
        for (const QcMetricRow& metric : metrics) {
            qc_write_metric_row(overview.body, metric);
        }
        overview.body << "</div>\n";
    }
    overview.body << "</div>\n<p class=\"qc-note\">Valid pairs and inter-chromosomal cutoffs follow "
                     "HiCExplorer's own usage guidance. The duplicated pairs, dangling ends and "
                     "self-ligation cutoffs have no published HiCExplorer standard behind them; they are "
                     "conservative starting points, and every cutoff is a hicPrepareQCreport command line "
                     "option (--minValidPairsPercent, --maxInterChromosomalPercent, "
                     "--maxDuplicatedPairsPercent, --maxDanglingEndPercent, --maxSelfLigationPercent). A "
                     "metric on the right side of its cutoff, or exactly on it, is PASS. A metric past its "
                     "cutoff by at most "
                  << margin_percent
                  << "% of the cutoff (the margin) is WARN, and one further past it is FAIL.</p>\n";

    // --- Pairs sequenced ---------------------------------------------------
    sequenced_section.id = "sequenced";
    sequenced_section.title = "Pairs sequenced";
    sequenced_section.help_page = "hicQC.html";
    sequenced_section.description =
        "Sequenced reads, pairs mappable, unique and of high quality, and the valid Hi-C contacts kept "
        "from them, per sample. Valid pairs as a percentage of sequenced reads are judged; their bars "
        "are hatched and drawn in the PASS, WARN or FAIL colour.";
    {
        QcChart chart{"plot-sequenced", "Pairs", false, {}};
        QcChartSeries s_seq{"Sequenced reads", {}, {}, {}};
        QcChartSeries s_map{"Mappable, unique, high quality", {}, {}, {}};
        QcChartSeries s_valid{"Hi-C contacts (valid pairs)", {}, {}, {}};
        std::vector<QcTableColumn> columns = {{"Sequenced reads", {}},
                                              {"Mappable, unique, high quality", {}},
                                              {"Mappable, unique, high quality %", {}},
                                              {"Hi-C contacts", {}},
                                              {"Hi-C contacts %", {}, {}}};
        for (std::size_t row = 0; row < n; ++row) {
            const double sequenced = count_of("Sequenced reads", row);
            const double mappable = count_of("Pairs mappable, unique and high quality", row);
            const double valid = count_of("Hi-C contacts", row);
            const QcStatus valid_status = judged[row].valid_status;
            s_seq.y.emplace_back(sequenced);
            s_seq.text.push_back(qc_bar_text(sequenced, sequenced, QcStatus::Unknown));
            s_seq.status.push_back(QcStatus::Unknown);
            s_map.y.emplace_back(mappable);
            s_map.text.push_back(qc_bar_text(mappable, sequenced, QcStatus::Unknown));
            s_map.status.push_back(QcStatus::Unknown);
            s_valid.y.emplace_back(valid);
            s_valid.text.push_back(qc_bar_text(valid, sequenced, valid_status));
            s_valid.status.push_back(valid_status);
            sequenced_section.status = qc_worst(sequenced_section.status, valid_status);
            columns[0].values.push_back(qc_format_count(sequenced));
            columns[1].values.push_back(qc_format_count(mappable));
            columns[2].values.push_back(sequenced > 0 ? qc_percent(100.0 * mappable / sequenced) : "n/a");
            columns[3].values.push_back(qc_format_count(valid));
            columns[4].values.push_back(sequenced > 0 ? qc_percent(100.0 * valid / sequenced) : "n/a");
            columns[4].statuses.push_back(valid_status);
            columns[4].judged_metric = "Valid pairs";
        }
        chart.series = {s_seq, s_map, s_valid};
        sequenced_section.body << "<div class=\"qc-plot\" id=\"" << chart.id << "\"></div>\n";
        charts.push_back(std::move(chart));
        qc_write_table(sequenced_section.body, "sequenced-table", frame.index, columns);
        qc_write_download(sequenced_section.body, "QC_table.txt");
    }

    // --- Unmappable and non-unique pairs ----------------------------------
    unmappable.id = "unmappable";
    unmappable.title = "Unmappable and non-unique pairs";
    unmappable.help_page = "hicQC.html";
    unmappable.description =
        "Fraction of sequenced reads, per sample, that mapped uniquely and at high quality (Hi-C "
        "contacts) against those that did not.";
    {
        static const char* const names[] = {"Hi-C contacts", "Low mapping quality", "One mate not unique",
                                            "One mate unmapped"};
        QcChart chart{"plot-unmappable", "% of sequenced reads", true, {}};
        std::vector<QcTableColumn> columns;
        for (const char* name : names) {
            columns.push_back({name, {}});
            columns.push_back({std::string(name) + " %", {}});
            chart.series.push_back({name, {}, {}, {}});
        }
        for (std::size_t row = 0; row < n; ++row) {
            const double sequenced = count_of("Sequenced reads", row);
            for (std::size_t k = 0; k < 4; ++k) {
                const double value = count_of(names[k], row);
                chart.series[k].y.push_back(qc_share(value, sequenced));
                chart.series[k].text.push_back(qc_bar_text(value, sequenced, QcStatus::Unknown));
                chart.series[k].status.push_back(QcStatus::Unknown);
                columns[2 * k].values.push_back(qc_format_count(value));
                columns[2 * k + 1].values.push_back(sequenced > 0 ? qc_percent(100.0 * value / sequenced)
                                                                  : "n/a");
            }
        }
        unmappable.body << "<div class=\"qc-plot\" id=\"" << chart.id << "\"></div>\n";
        charts.push_back(std::move(chart));
        qc_write_table(unmappable.body, "unmappable-table", frame.index, columns);
        qc_write_download(unmappable.body, "unmapable_table.txt");
    }

    // --- Pairs discarded ---------------------------------------------------
    discarded.id = "discarded";
    discarded.title = "Pairs discarded";
    discarded.help_page = "hicQC.html";
    discarded.description =
        "Fraction of mappable, unique, high quality pairs discarded when building the Hi-C matrix, per "
        "sample. Duplicated pairs, dangling ends and self-ligation are judged; their bars are hatched and "
        "drawn in the PASS, WARN or FAIL colour.";
    {
        static const char* const prefixes[] = {"One mate not close to rest site", "dangling end",
                                                "duplicated pairs", "same fragment", "self circle",
                                                "self ligation (removed)"};
        // Columns are built after the per-sample pass below, since a prefix
        // such as "dangling end" can match more than one real column (one
        // per restriction enzyme) and every sample must get every column
        // that *any* sample has, in a fixed order.
        std::vector<std::string> discovered_columns;
        std::vector<std::unordered_map<std::string, double>> per_row_values(n);
        std::vector<std::unordered_map<std::string, QcStatus>> per_row_status(n);
        for (std::size_t row = 0; row < n; ++row) {
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
                    per_row_status[row][column] = status;
                    discarded.status = qc_worst(discarded.status, status);
                }
            }
        }
        QcChart chart{"plot-discarded", "% of mappable pairs", true, {}};
        std::vector<QcTableColumn> columns;
        for (const std::string& column : discovered_columns) {
            columns.push_back({column, {}});
            columns.push_back({column + " %", {}});
            chart.series.push_back({column, {}, {}, {}});
        }
        for (std::size_t row = 0; row < n; ++row) {
            const double mappable = count_of("Pairs mappable, unique and high quality", row);
            for (std::size_t c = 0; c < discovered_columns.size(); ++c) {
                const auto it = per_row_values[row].find(discovered_columns[c]);
                const double value = it == per_row_values[row].end() ? 0.0 : it->second;
                const auto st = per_row_status[row].find(discovered_columns[c]);
                const QcStatus status = st == per_row_status[row].end() ? QcStatus::Unknown : st->second;
                chart.series[c].y.push_back(qc_share(value, mappable));
                chart.series[c].text.push_back(qc_bar_text(value, mappable, status));
                chart.series[c].status.push_back(status);
                columns[2 * c].values.push_back(qc_format_count(value));
                columns[2 * c + 1].values.push_back(mappable > 0 ? qc_percent(100.0 * value / mappable)
                                                                 : "n/a");
                columns[2 * c + 1].statuses.push_back(status);
                columns[2 * c + 1].judged_metric = discovered_columns[c];
            }
        }
        discarded.body << "<div class=\"qc-plot\" id=\"" << chart.id << "\"></div>\n";
        charts.push_back(std::move(chart));
        discarded.body << "<dl class=\"qc-defs\">"
                          "<dt>Dangling ends</dt><dd>Reads that start with the restriction site: digested "
                          "but not ligated.</dd>"
                          "<dt>Same fragment</dt><dd>Inward-facing mates up to 800 bp apart with no "
                          "restriction site between them; not a valid Hi-C pair.</dd>"
                          "<dt>Self circle</dt><dd>Pairs within 25 kb with outward read orientation.</dd>"
                          "<dt>Self ligation</dt><dd>Pairs with a restriction site between them, within "
                          "800 bp.</dd></dl>\n";
        qc_write_table(discarded.body, "discarded-table", frame.index, columns);
        qc_write_download(discarded.body, "discarded_table.txt");
    }

    // --- Contact distance ---------------------------------------------------
    distance.id = "distance";
    distance.title = "Contact distance";
    distance.help_page = "hicQC.html";
    distance.description =
        "Fraction of valid Hi-C contacts, per sample, by genomic distance. HiCExplorer's own usage "
        "guidance: good Hi-C libraries have less than 10% inter-chromosomal contacts. Inter-chromosomal "
        "contacts are judged; their bars are hatched and drawn in the PASS, WARN or FAIL colour.";
    {
        static const char* const names[] = {"inter chromosomal", "Intra short range (< 20kb)",
                                            "Intra long range (>= 20kb)"};
        QcChart chart{"plot-distance", "% of valid pairs", true, {}};
        std::vector<QcTableColumn> columns;
        for (const char* name : names) {
            columns.push_back({name, {}});
            columns.push_back({std::string(name) + " %", {}});
            chart.series.push_back({name, {}, {}, {}});
        }
        for (std::size_t row = 0; row < n; ++row) {
            const double valid = count_of("Hi-C contacts", row);
            for (std::size_t k = 0; k < 3; ++k) {
                const double value = count_of(names[k], row);
                const QcStatus status = k == 0 ? judged[row].inter_chrom_status : QcStatus::Unknown;
                distance.status = qc_worst(distance.status, status);
                chart.series[k].y.push_back(qc_share(value, valid));
                chart.series[k].text.push_back(qc_bar_text(value, valid, status));
                chart.series[k].status.push_back(status);
                columns[2 * k].values.push_back(qc_format_count(value));
                columns[2 * k + 1].values.push_back(valid > 0 ? qc_percent(100.0 * value / valid) : "n/a");
                columns[2 * k + 1].statuses.push_back(status);
                columns[2 * k + 1].judged_metric = names[k];
            }
        }
        distance.body << "<div class=\"qc-plot\" id=\"" << chart.id << "\"></div>\n";
        charts.push_back(std::move(chart));
        qc_write_table(distance.body, "distance-table", frame.index, columns);
        qc_write_download(distance.body, "distance_table.txt");
    }

    // --- Read orientation ---------------------------------------------------
    orientation.id = "orientation";
    orientation.title = "Read orientation";
    orientation.help_page = "hicQC.html";
    orientation.description =
        "Fraction of the four read-pair orientations, per sample, informative for ligation bias.";
    {
        static const char* const names[] = {"Read pair type: inward pairs", "Read pair type: outward pairs",
                                            "Read pair type: left pairs", "Read pair type: right pairs"};
        QcChart chart{"plot-orientation", "% of read pairs", true, {}};
        std::vector<QcTableColumn> columns;
        for (const char* name : names) {
            columns.push_back({name, {}});
            columns.push_back({std::string(name) + " %", {}});
            chart.series.push_back({std::string(name).substr(16), {}, {}, {}});
        }
        for (std::size_t row = 0; row < n; ++row) {
            double total = 0.0;
            for (const char* name : names) {
                total += count_of(name, row);
            }
            for (std::size_t k = 0; k < 4; ++k) {
                const double value = count_of(names[k], row);
                chart.series[k].y.push_back(qc_share(value, total));
                chart.series[k].text.push_back(qc_bar_text(value, total, QcStatus::Unknown));
                chart.series[k].status.push_back(QcStatus::Unknown);
                columns[2 * k].values.push_back(qc_format_count(value));
                columns[2 * k + 1].values.push_back(total > 0 ? qc_percent(100.0 * value / total) : "n/a");
            }
        }
        orientation.body << "<div class=\"qc-plot\" id=\"" << chart.id << "\"></div>\n";
        charts.push_back(std::move(chart));
        qc_write_table(orientation.body, "orientation-table", frame.index, columns);
        qc_write_download(orientation.body, "read_orientation_table.txt");
    }

    // --- Full table ---------------------------------------------------------
    full.id = "fulltable";
    full.title = "Full table";
    full.help_page = "hicPrepareQCreport.html";
    full.description = "Every field QC.log reports, one row per sample. Click a header to sort.";
    {
        std::vector<QcTableColumn> columns;
        for (const std::string& name : frame.columns) {
            QcTableColumn column{name, {}};
            for (std::size_t row = 0; row < n; ++row) {
                column.values.push_back(frame.column(name)[row].str());
            }
            columns.push_back(std::move(column));
        }
        qc_write_table(full.body, "full-table", frame.index, columns);
        qc_write_download(full.body, "QC_table.txt");
    }

    const std::vector<QcSection*> sections = {&overview,    &sequenced_section, &unmappable, &discarded,
                                              &distance,    &orientation,       &full};

    std::ofstream out(folder + "/hicQC_report.html", std::ios::binary);
    if (!out) {
        throw std::runtime_error("could not write " + folder + "/hicQC_report.html");
    }

    out << R"HTML(<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<script>(function(){var d=document.documentElement,t=null,m=/[?&]theme=(light|dark)/.exec(location.search);if(m){t=m[1];}else{try{t=localStorage.getItem('hicx-theme');}catch(e){}}if(t!=='light'&&t!=='dark')t='dark';d.setAttribute('data-theme',t);})();</script>
<title>HiCExplorer QC report</title>
<style>
:root {
  --bg:#f5f5f7; --surface:#ffffff; --surface-2:#f5f5f7; --border:#e3e3e8;
  --text:#1d1d1f; --muted:#6e6e73; --brand:#0071e3; --brand-2:#42a5ff;
  --pass:#34c759; --warn:#ff9f0a; --fail:#ff3b30;
  --pass-ink:#1b7a43; --pass-bg:#e6f4ea; --warn-ink:#8a5a00; --warn-bg:#fff3cd;
  --fail-ink:#a3242a; --fail-bg:#fdecea; --unknown-ink:#555; --unknown-bg:#eee;
  --header-bg:rgba(255,255,255,.72);
  --shadow:0 8px 30px rgba(0,0,0,.06),0 1px 2px rgba(0,0,0,.05);
  --shadow-lg:0 18px 50px rgba(0,0,0,.10);
  --radius:18px; --sidebar-w:17em;
}
:root[data-theme="dark"] {
  --bg:#000000; --surface:#1c1c1e; --surface-2:#2c2c2e; --border:#38383a;
  --text:#f5f5f7; --muted:#98989d; --brand:#0a84ff; --brand-2:#64d2ff;
  --pass:#30d158; --warn:#ff9f0a; --fail:#ff453a;
  --pass-ink:#30d158; --pass-bg:rgba(48,209,88,.16); --warn-ink:#ffb340; --warn-bg:rgba(255,159,10,.16);
  --fail-ink:#ff6961; --fail-bg:rgba(255,69,58,.16); --unknown-ink:#98989d; --unknown-bg:#2c2c2e;
  --header-bg:rgba(28,28,30,.72);
  --shadow:0 8px 30px rgba(0,0,0,.5),0 1px 2px rgba(0,0,0,.6);
  --shadow-lg:0 18px 50px rgba(0,0,0,.6);
}
* { box-sizing:border-box; }
html { scroll-behavior:smooth; }
body {
  color:var(--text); margin:0; padding:0; background:var(--bg);
  font-family:-apple-system,BlinkMacSystemFont,'SF Pro Text','Segoe UI',Roboto,Helvetica,Arial,sans-serif;
  font-size:15px; line-height:1.55; -webkit-font-smoothing:antialiased; letter-spacing:-.006em;
  -webkit-print-color-adjust:exact; print-color-adjust:exact;
}
a { color:var(--brand); text-decoration:none; }
a:hover { text-decoration:underline; }
div.header {
  position:sticky; top:0; z-index:20; display:flex; align-items:center;
  justify-content:space-between; gap:.6em 1em; flex-wrap:wrap; width:100%; margin:0;
  background:var(--header-bg); border-bottom:1px solid var(--border);
  -webkit-backdrop-filter:saturate(180%) blur(20px); backdrop-filter:saturate(180%) blur(20px);
  padding:.9em 1.6em;
}
div.header-left { display:flex; align-items:center; gap:.5em; flex-wrap:wrap; min-width:0; }
#header_title { font-size:1.6rem; font-weight:700; letter-spacing:-.022em; color:var(--text); margin-right:.3em; }
#header_info {
  text-align:right; font-size:.82rem; color:var(--muted);
  font-family:ui-monospace,SFMono-Regular,Menlo,monospace;
}
div.summary {
  position:fixed; top:4.6em; left:0; width:var(--sidebar-w);
  height:calc(100vh - 6.6em); overflow-y:auto; margin:0; padding:1em .9em;
}
div.summary h2 {
  font-size:.72rem; text-transform:uppercase; letter-spacing:.08em;
  color:var(--muted); font-weight:700; margin:.2em 0 .7em .5em;
}
div.summary ul { padding:0; margin:0; list-style:none; }
div.summary ul li a {
  display:flex; align-items:center; gap:.6em; padding:.5em .7em;
  margin-bottom:2px; border-radius:8px; color:var(--text);
  font-size:.9rem; font-weight:500; transition:background .12s;
}
div.summary ul li a:hover { background:var(--surface-2); text-decoration:none; }
div.summary ul li a::before, div.module h2::before, div.qc-card h3::before {
  content:""; flex:0 0 auto; width:9px; height:9px; border-radius:50%; background:var(--muted);
}
div.summary ul li a.pass::before, div.module h2.pass::before, div.qc-card h3.pass::before { background:var(--pass); box-shadow:0 0 0 3px rgba(52,199,89,.16); }
div.summary ul li a.warn::before, div.module h2.warn::before, div.qc-card h3.warn::before { background:var(--warn); box-shadow:0 0 0 3px rgba(255,159,10,.16); }
div.summary ul li a.fail::before, div.module h2.fail::before, div.qc-card h3.fail::before { background:var(--fail); box-shadow:0 0 0 3px rgba(255,59,48,.16); }
a.toggle-all {
  display:inline-block; margin:.1em 0 .7em .5em; font-size:.72rem;
  color:var(--brand); cursor:pointer; font-weight:600;
}
div.main {
  margin-left:max(var(--sidebar-w), calc((100% - 1180px) / 2)); margin-right:auto;
  padding:1.5em 2em 4.5em; max-width:1180px; min-width:0;
}
div.module {
  background:var(--surface); border:1px solid var(--border); border-radius:var(--radius);
  box-shadow:var(--shadow); padding:1.35em 1.7em 1.6em; margin-bottom:1.5em; min-width:0;
  transition:box-shadow .2s ease;
}
div.module:hover { box-shadow:var(--shadow-lg); }
div.module h2 {
  display:flex; align-items:center; gap:.55em; scroll-margin-top:5.5em;
  font-size:1.12rem; font-weight:650; color:var(--text);
  margin:0 0 1em; padding:0 0 .7em; border-bottom:1px solid var(--border);
}
div.module h2::before { width:11px; height:11px; }
div.module h2 span.title { flex:1 1 auto; min-width:0; }
div.module h2.collapsible { cursor:pointer; user-select:none; }
div.module h2.collapsible::after {
  content:""; margin-left:auto; flex:0 0 auto; width:.5em; height:.5em; margin-top:-.15em;
  border-right:2px solid var(--muted); border-bottom:2px solid var(--muted);
  transform:rotate(45deg); transition:transform .15s;
}
div.module.collapsed h2.collapsible { margin-bottom:0; padding-bottom:0; border-bottom:0; }
div.module.collapsed h2.collapsible::after { transform:rotate(-45deg); margin-top:.15em; }
div.module.collapsed { padding-bottom:1em; }
div.module.collapsed > :not(h2) { display:none; }
p.module-desc { margin:-.3em 0 1.15em; color:var(--muted); font-size:.86rem; line-height:1.5; max-width:78ch; }
p.module-desc a { font-weight:600; white-space:nowrap; }
div.qc-plot { width:100%; min-width:0; margin:0 0 1em; }
.qc-samples { display:grid; grid-template-columns:repeat(auto-fit, minmax(min(100%, 260px), 1fr)); gap:14px; }
.qc-card { background:var(--surface-2); border:1px solid var(--border); border-radius:12px; padding:14px 16px; min-width:0; }
.qc-card h3 { display:flex; align-items:center; gap:.5em; font-size:.98rem; margin:0 0 2px; overflow-wrap:anywhere; }
.qc-card .qc-overall { margin:0 0 10px; font-size:.8rem; color:var(--muted); }
.qc-metric { border-top:1px solid var(--border); padding:7px 0; }
.qc-metric:first-of-type { border-top:none; }
.qc-metric-top { display:flex; justify-content:space-between; align-items:center; gap:.5em; }
.qc-metric-label { font-size:.85rem; }
.qc-metric-value { font-size:1.15rem; font-weight:600; margin-top:2px; font-variant-numeric:tabular-nums; }
.qc-metric-rule { font-size:.72rem; color:var(--muted); margin-top:1px; }
.qc-badge { font-size:.68rem; font-weight:700; letter-spacing:.03em; padding:2px 8px; border-radius:999px; white-space:nowrap; font-family:-apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif; }
.qc-badge-sm { font-size:.62rem; padding:1px 6px; }
.qc-sr { position:absolute; width:1px; height:1px; overflow:hidden; clip:rect(0 0 0 0); white-space:nowrap; }
.qc-pass { color:var(--pass-ink); background:var(--pass-bg); }
.qc-warn { color:var(--warn-ink); background:var(--warn-bg); }
.qc-fail { color:var(--fail-ink); background:var(--fail-bg); }
.qc-unknown { color:var(--unknown-ink); background:var(--unknown-bg); }
.qc-note { font-size:.78rem; color:var(--muted); line-height:1.5; margin:14px 0 0; max-width:90ch; }
dl.qc-defs { font-size:.85rem; background:var(--surface-2); border:1px solid var(--border); border-radius:12px; padding:12px 16px; margin:0 0 14px; }
dl.qc-defs dt { font-weight:700; margin-top:8px; }
dl.qc-defs dt:first-child { margin-top:0; }
dl.qc-defs dd { margin:2px 0 0; color:var(--muted); }
.qc-table-wrap { position:relative; overflow-x:auto; border:1px solid var(--border); border-radius:10px; max-width:100%; }
table.qc-table { border-collapse:collapse; width:100%; font-size:.84rem; }
table.qc-table th {
  background:var(--surface-2); color:var(--muted); text-transform:uppercase; font-size:.69rem;
  letter-spacing:.04em; font-weight:700; padding:.6em .9em; border-bottom:1px solid var(--border);
  text-align:right; white-space:nowrap; cursor:pointer; user-select:none;
}
table.qc-table th:hover { color:var(--text); }
table.qc-table td {
  font-family:ui-monospace,SFMono-Regular,Menlo,monospace; background:var(--surface); color:var(--text);
  padding:.5em .9em; border-bottom:1px solid var(--border); text-align:right; white-space:nowrap;
}
table.qc-table tr:last-child td { border-bottom:0; }
table.qc-table tr:hover td { background:var(--surface-2); }
table.qc-table th:first-child, table.qc-table td:first-child { text-align:left; position:sticky; left:0; z-index:1; }
table.qc-table td:first-child { font-family:inherit; font-weight:600; }
.qc-download { font-size:.8rem; margin-top:8px; display:inline-block; font-weight:600; }
div.footer {
  position:fixed; bottom:0; left:0; width:100%; z-index:20; margin:0;
  background:var(--header-bg); border-top:1px solid var(--border);
  -webkit-backdrop-filter:saturate(180%) blur(20px); backdrop-filter:saturate(180%) blur(20px);
  color:var(--muted); font-size:.8rem; padding:.55em 1.6em;
}
button.hx-btn, select.palette-select {
  border:1px solid var(--border); background:var(--surface-2); color:var(--text);
  border-radius:999px; height:2.1em; cursor:pointer; flex:0 0 auto; font-family:inherit;
}
button.hx-btn {
  width:2.1em; font-size:1rem; line-height:1; display:inline-flex; align-items:center;
  justify-content:center; padding:0; transition:background .15s, transform .15s;
}
button.hx-btn svg { width:1em; height:1em; display:block; }
button.hx-btn:hover { background:var(--brand); color:#fff; transform:scale(1.08); }
select.palette-select {
  font-size:.82rem; font-weight:600; line-height:1; padding:0 1.9em 0 .9em;
  -webkit-appearance:none; appearance:none;
  background-image:url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='10' height='6'><path d='M1 1l4 4 4-4' fill='none' stroke='%236e6e73' stroke-width='1.6'/></svg>");
  background-repeat:no-repeat; background-position:right .8em center;
}
select.palette-select:hover { border-color:var(--brand); }
[hidden] { display:none !important; }
body.nav-collapsed div.summary { display:none; }
body.nav-collapsed div.main { margin-left:max(0px, calc((100% - 1180px) / 2)); }
.js-plotly-plot .modebar-btn { display:inline-flex !important; align-items:center; justify-content:center; }
.js-plotly-plot .modebar-btn svg { vertical-align:middle; display:block; }
@media (max-width:820px) {
  div.header { padding:.7em 1em; }
  #header_title { font-size:1.3rem; }
  #header_info { text-align:left; }
  div.summary { position:static; width:auto; height:auto; margin:.8em 16px;
    background:var(--surface); border:1px solid var(--border); border-radius:var(--radius); }
  div.summary ul { display:flex; flex-wrap:wrap; gap:.25em; }
  div.main { margin-left:0; padding:.5em 16px 4.5em; }
  body.nav-collapsed div.main { margin-left:0; }
  div.module { padding:1.1em 1em 1.2em; }
}
@media print {
  div.summary, div.footer, div.header { position:static; width:auto; height:auto; }
  div.main { margin-left:0; } div.module { break-inside:avoid; box-shadow:none; }
  .modebar, button.hx-btn, select.palette-select { display:none !important; }
}
</style>
</head>
<body>
<div class="header">
<div class="header-left">
<button class="hx-btn" id="nav-toggle" type="button" title="Show or hide the summary sidebar" aria-label="Toggle sidebar" hidden>&#9776;</button>
<div id="header_title">HiCExplorer QC</div>
<button class="hx-btn" id="theme-toggle" type="button" title="Toggle light or dark theme" aria-label="Toggle theme" hidden></button>
<select class="palette-select" id="palette-select" title="Plot colour scheme" aria-label="Plot colour scheme" hidden>
<option value="ocean">Ocean</option><option value="sunset">Sunset</option><option value="ink">Ink</option>
</select>
</div>
)HTML";
    out << "<div id=\"header_info\">" << qc_now_iso() << "<br>" << n << " sample" << (n == 1 ? "" : "s")
        << "</div>\n</div>\n";

    // --- Sidebar: every section with its worst verdict ---------------------
    out << "<div class=\"summary\"><h2>Summary</h2>\n<ul>\n";
    for (const QcSection* section : sections) {
        out << "<li><a class=\"" << qc_status_word(section->status) << "\" href=\"#" << section->id
            << "\" title=\"" << qc_status_title(section->status) << "\">" << qc_html_escape(section->title)
            << "</a></li>\n";
    }
    out << "</ul>\n</div>\n<div class=\"main\">\n";

    for (const QcSection* section : sections) {
        out << "<div class=\"module\" id=\"module-" << section->id << "\">\n<h2 class=\""
            << qc_status_word(section->status) << "\" id=\"" << section->id << "\">"
            << "<span class=\"title\">" << qc_html_escape(section->title) << ": "
            << qc_status_title(section->status) << "</span></h2>\n"
            << "<p class=\"module-desc\">" << section->description << " <a href=\"" << kQcDocsBase
            << section->help_page << "\" target=\"_blank\" rel=\"noopener\">HiCExplorer documentation</a></p>\n"
            << section->body.str() << "</div>\n";
    }
    out << "</div>\n<div class=\"footer\">HiCExplorer " << kVersion << "</div>\n";

    out << "<script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\" charset=\"utf-8\"></script>\n"
        << "<script>\nvar HICX_CHARTS = [\n";
    for (std::size_t c = 0; c < charts.size(); ++c) {
        out << (c ? ",\n" : "") << qc_chart_json(charts[c], frame.index);
    }
    out << "\n];\n";
    out << R"HTML((function () {
  var root = document.documentElement, body = document.body;
  function store(key, value) { try { localStorage.setItem(key, value); } catch (e) {} }
  function load(key) { try { return localStorage.getItem(key); } catch (e) { return null; } }

  // Colourblind-safe categorical palettes (as in falcoAmadeus), each with a
  // variant per background; none mixes red with green.
  var PALETTES = {
    ocean:  {light: ['#0071e3', '#00c7be', '#5856d6', '#ff9500', '#1d1d1f', '#8e8e93'],
             dark:  ['#4aa3ff', '#40e0d6', '#a29bff', '#ffb340', '#f2f2f7', '#c7c7cc']},
    sunset: {light: ['#ff9500', '#ff2d55', '#af52de', '#a2845e', '#1d1d1f', '#0071e3'],
             dark:  ['#ffb340', '#ff6482', '#d68bff', '#d6b48c', '#f2f2f7', '#64d2ff']},
    ink:    {light: ['#1d1d1f', '#0071e3', '#ff9500', '#af52de', '#8a6d1f', '#6e6e73'],
             dark:  ['#f2f2f7', '#4aa3ff', '#ffb340', '#d68bff', '#d4b06a', '#aeaeb2']}
  };
  var STATUS = {light: {pass: '#34c759', warn: '#ff9f0a', fail: '#ff3b30'},
                dark:  {pass: '#30d158', warn: '#ff9f0a', fail: '#ff453a'}};
  var PATTERNS = ['/', '\\', 'x', '.', '+', '-'];
  var SUN = '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><circle cx="12" cy="12" r="4"/><path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4"/></svg>';
  var MOON = '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"><path d="M21 12.8A9 9 0 1 1 11.2 3a7 7 0 0 0 9.8 9.8z"/></svg>';

  var palette = 'ocean';
  var qp = /[?&]palette=(ocean|sunset|ink)/.exec(location.search);
  if (qp) { palette = qp[1]; } else { var sp = load('hicx-palette'); if (sp && PALETTES[sp]) palette = sp; }

  function legendRows(c, perRow) { return Math.ceil(c.series.length / perRow); }
  function dark() { return root.getAttribute('data-theme') === 'dark'; }

  function drawCharts() {
    if (!window.Plotly) return;
    var d = dark(), mode = d ? 'dark' : 'light';
    var colours = PALETTES[palette][mode], st = STATUS[mode];
    var fg = d ? '#e6ebf2' : '#2a3441', grid = d ? '#3a4451' : '#d6dbe3',
        line = d ? '#4a5462' : '#b7bec9', tick = d ? '#9aa6b2' : '#697586';
    HICX_CHARTS.forEach(function (c) {
      var el = document.getElementById(c.id);
      if (!el) return;
      var k = 0, j = 0, traces = [];
      c.series.forEach(function (s) {
        var judged = s.status.some(function (x) { return x; });
        var base = colours[k % colours.length];
        if (!judged) k++;
        var marker = {color: judged ? s.status.map(function (x) { return x ? st[x] : base; }) : base};
        var group = 'g' + traces.length;
        if (judged) {
          var pattern = {shape: PATTERNS[j++ % PATTERNS.length], size: 7, solidity: 0.3,
                         fgcolor: d ? 'rgba(0,0,0,0.55)' : 'rgba(255,255,255,0.75)'};
          marker.pattern = pattern;
          // The bars carry one verdict colour each, so the legend gets a
          // neutral hatched swatch of its own.
          traces.push({type: 'bar', name: s.name, x: [c.samples[0]], y: [null], legendgroup: group, xaxis: 'x2',
                       hoverinfo: 'skip',
                       marker: {color: d ? '#636366' : '#d1d1d6',
                                pattern: {shape: pattern.shape, size: 5, solidity: 0.45, fgcolor: fg}}});
        }
        traces.push({type: 'bar', name: s.name, x: c.samples, y: s.y, text: s.text, textposition: 'none',
                     legendgroup: group, showlegend: !judged,
                     hovertemplate: '%{fullData.name}: %{text}<extra></extra>', marker: marker});
      });
      var axis = {gridcolor: grid, zerolinecolor: grid, linecolor: line, tickcolor: line, ticks: 'outside',
                  tickfont: {size: 11, color: tick}, automargin: true, showline: true};
      var layout = {
        barmode: 'group', bargap: 0.25, bargroupgap: 0.06,
        // Room for the legend rows below the axis: one entry per row on a narrow screen.
        height: el.clientWidth < 600 ? 340 + 24 * legendRows(c, 1) : 340 + 22 * legendRows(c, 4),
        paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
        font: {family: '-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif', size: 12, color: fg},
        margin: {t: 34, r: 16, b: 16, l: 60},
        xaxis: Object.assign({}, axis, {type: 'category', showgrid: false}),
        xaxis2: {type: 'category', overlaying: 'x', visible: false},
        yaxis: Object.assign({}, axis, {title: {text: c.ytitle, font: {color: fg}}, rangemode: 'tozero',
                                        ticksuffix: c.percent ? '%' : '', separatethousands: true}),
        legend: {orientation: 'h', x: 0, xanchor: 'left', yref: 'container', y: 0, yanchor: 'bottom',
                 font: {size: 11, color: fg}},
        hovermode: 'x unified',
        hoverlabel: {bgcolor: d ? '#2c2c2e' : '#ffffff', bordercolor: d ? '#48484a' : '#d6dbe3',
                     font: {color: d ? '#f5f5f7' : '#1d1d1f', size: 12}}
      };
      Plotly.react(el, traces, layout, {displaylogo: false, responsive: true, displayModeBar: 'hover',
        modeBarButtonsToRemove: ['select2d', 'lasso2d', 'autoScale2d', 'toggleSpikelines',
                                 'hoverClosestCartesian', 'hoverCompareCartesian']});
    });
  }
  function resizeCharts() {
    if (!window.Plotly) return;
    document.querySelectorAll('.js-plotly-plot').forEach(function (el) {
      if (el.offsetParent) Plotly.Plots.resize(el);
    });
  }

  // Theme toggle (default dark) and palette picker.
  var themeBtn = document.getElementById('theme-toggle');
  function setTheme(t) {
    root.setAttribute('data-theme', t);
    store('hicx-theme', t);
    themeBtn.innerHTML = t === 'dark' ? SUN : MOON;
    drawCharts();
  }
  themeBtn.onclick = function () { setTheme(dark() ? 'light' : 'dark'); };
  themeBtn.innerHTML = dark() ? SUN : MOON;
  themeBtn.hidden = false;
  var sel = document.getElementById('palette-select');
  sel.value = palette;
  sel.onchange = function () { palette = sel.value; store('hicx-palette', palette); drawCharts(); };
  sel.hidden = false;

  // Collapsible summary sidebar; the state is kept in localStorage.
  var navBtn = document.getElementById('nav-toggle');
  function setNav(collapsed) {
    body.classList.toggle('nav-collapsed', collapsed);
    store('hicx-nav', collapsed ? '1' : '0');
    resizeCharts();
  }
  navBtn.onclick = function () { setNav(!body.classList.contains('nav-collapsed')); };
  navBtn.hidden = false;
  body.classList.toggle('nav-collapsed', load('hicx-nav') === '1');

  // Collapsible cards and "collapse all".
  var modules = document.querySelectorAll('div.module');
  modules.forEach(function (m) {
    var h = m.querySelector('h2');
    h.classList.add('collapsible');
    h.addEventListener('click', function () {
      m.classList.toggle('collapsed');
      if (!m.classList.contains('collapsed')) resizeCharts();
    });
  });
  var sumTitle = document.querySelector('div.summary h2');
  var all = document.createElement('a');
  all.className = 'toggle-all';
  all.textContent = 'collapse all';
  all.addEventListener('click', function () {
    var collapse = all.textContent === 'collapse all';
    modules.forEach(function (m) { m.classList.toggle('collapsed', collapse); });
    all.textContent = collapse ? 'expand all' : 'collapse all';
    if (!collapse) resizeCharts();
  });
  sumTitle.parentNode.insertBefore(all, sumTitle.nextSibling);
  // A sidebar link opens its card when it is collapsed.
  document.querySelectorAll('div.summary a[href^="#"]').forEach(function (a) {
    a.addEventListener('click', function () {
      var t = document.getElementById(a.getAttribute('href').slice(1));
      if (t && t.parentNode.classList.contains('collapsed')) {
        t.parentNode.classList.remove('collapsed');
        resizeCharts();
      }
    });
  });

  // Sortable tables: click a header to sort by that column.
  document.querySelectorAll('table.qc-table').forEach(function (table) {
    var tbody = table.querySelector('tbody');
    table.querySelectorAll('thead th').forEach(function (th, index) {
      th.addEventListener('click', function () {
        var ascending = th.getAttribute('data-asc') !== 'true';
        table.querySelectorAll('thead th').forEach(function (h) { h.removeAttribute('data-asc'); });
        th.setAttribute('data-asc', ascending ? 'true' : 'false');
        var rows = Array.prototype.slice.call(tbody.querySelectorAll('tr'));
        rows.sort(function (a, b) {
          var av = a.children[index].textContent.replace(/,/g, '');
          var bv = b.children[index].textContent.replace(/,/g, '');
          var an = parseFloat(av), bn = parseFloat(bv);
          var cmp = (!isNaN(an) && !isNaN(bn)) ? an - bn : av.localeCompare(bv);
          return ascending ? cmp : -cmp;
        });
        rows.forEach(function (row) { tbody.appendChild(row); });
      });
    });
  });

  drawCharts();
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
