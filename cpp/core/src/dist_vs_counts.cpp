#include "hicx/dist_vs_counts.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <optional>
#include <string>
#include <unordered_map>

#include "hicx/python_error.hpp"

namespace hicx::dist_counts {

namespace {

bool starts_with(const std::string& text, const char* prefix) {
    return text.rfind(prefix, 0) == 0;
}

std::int64_t python_mod(std::int64_t a, std::int64_t m) {
    const std::int64_t r = a % m;
    return (r != 0 && ((r < 0) != (m < 0))) ? r + m : r;
}

}  // namespace

std::vector<std::int64_t> fitted_starts(const std::vector<CutInterval>& intervals, bool* snapped,
                                        std::int64_t* median_out) {
    if (snapped != nullptr) {
        *snapped = false;
    }
    std::vector<std::int64_t> starts;
    starts.reserve(intervals.size());
    for (const CutInterval& interval : intervals) {
        starts.push_back(interval.start);
    }
    if (intervals.size() <= 1) {
        return starts;
    }
    std::vector<std::string> order;
    std::unordered_map<std::string, std::vector<std::int64_t>> by_chrom;
    for (const CutInterval& interval : intervals) {
        auto [it, inserted] = by_chrom.try_emplace(interval.chrom);
        if (inserted) {
            order.push_back(interval.chrom);
        }
        it->second.push_back(interval.start);
    }
    std::vector<std::int64_t> diffs;
    for (const std::string& chrom : order) {
        const std::vector<std::int64_t>& values = by_chrom[chrom];
        for (std::size_t i = 1; i < values.size(); ++i) {
            diffs.push_back(values[i] - values[i - 1]);
        }
    }
    if (diffs.empty()) {
        throw PythonError("ValueError", "need at least one array to concatenate");
    }
    std::sort(diffs.begin(), diffs.end());
    const std::size_t n = diffs.size();
    const double median_value =
        n % 2 == 1 ? static_cast<double>(diffs[n / 2])
                   : (static_cast<double>(diffs[n / 2 - 1]) + static_cast<double>(diffs[n / 2])) / 2.0;
    const auto median = static_cast<std::int64_t>(median_value);
    if (median_out != nullptr) {
        *median_out = median;
    }
    std::size_t deviating = 0;
    for (const CutInterval& interval : intervals) {
        if (interval.end - interval.start != median) {
            ++deviating;
        }
    }
    if (static_cast<double>(deviating) > static_cast<double>(intervals.size()) * 0.01) {
        if (median == 0) {
            throw PythonError("ZeroDivisionError", "integer modulo by zero");
        }
        for (std::int64_t& start : starts) {
            const std::int64_t down = -python_mod(start, median);
            const std::int64_t up = python_mod(-start, median);
            start += std::llabs(down) <= std::llabs(up) ? down : up;
        }
        if (snapped != nullptr) {
            *snapped = true;
        }
    }
    return starts;
}

std::vector<std::int64_t> unit_sizes(const std::vector<CutInterval>& intervals) {
    std::vector<std::string> order;
    std::unordered_map<std::string, std::pair<std::int64_t, std::int64_t>> bounds;
    if (intervals.empty()) {
        return {};
    }
    std::int64_t chr_start_id = 0;
    std::optional<std::string> previous;
    auto assign = [&](const std::string& name, std::int64_t first, std::int64_t last) {
        if (bounds.find(name) == bounds.end()) {
            order.push_back(name);
        }
        bounds[name] = {first, last};
    };
    std::int64_t id = 0;
    for (const CutInterval& interval : intervals) {
        if (!previous.has_value() || *previous != interval.chrom) {
            if (!previous.has_value()) {
                previous = interval.chrom;
            }
            assign(*previous, chr_start_id, id);
            chr_start_id = id;
            previous = interval.chrom;
        }
        ++id;
    }
    assign(intervals.back().chrom, chr_start_id, id);
    std::vector<std::int64_t> sizes;
    for (const std::string& name : order) {
        if (!starts_with(name, "_ignore_")) {
            sizes.push_back(bounds[name].second - bounds[name].first);
        }
    }
    return sizes;
}

std::vector<Distance> scan(const std::vector<Entry>& entries, std::size_t begin, std::size_t end,
                           std::int32_t offset, std::int64_t unit_bins, const std::vector<CutInterval>& intervals,
                           std::int64_t binsize, bool maxdepth_set) {
    const std::vector<std::int64_t> starts = fitted_starts(intervals);
    std::unordered_map<std::string, std::int64_t> ids;
    std::vector<std::int64_t> chrom_id(intervals.size());
    std::vector<char> ignored(intervals.size());
    for (std::size_t i = 0; i < intervals.size(); ++i) {
        chrom_id[i] = ids.try_emplace(intervals[i].chrom, static_cast<std::int64_t>(ids.size())).first->second;
        ignored[i] = starts_with(intervals[i].chrom, "_ignore_") ? 1 : 0;
    }

    std::vector<double> sums;
    std::vector<std::int64_t> counts;
    for (std::size_t k = begin; k < end; ++k) {
        const auto row = static_cast<std::size_t>(entries[k].row - offset);
        const auto col = static_cast<std::size_t>(entries[k].col - offset);
        std::int64_t dist = 0;
        if (chrom_id[row] != chrom_id[col]) {
            dist = -1;  // inter unit; its chromosome name is '' and it is kept
        } else {
            if (ignored[row]) {
                continue;
            }
            dist = starts[col] - starts[row];
        }
        if (dist == -1) {
            dist = -binsize;
        }
        const std::int64_t bin = static_cast<std::int64_t>(static_cast<double>(dist) / static_cast<double>(binsize)) + 1;
        if (bin < 0) {
            throw PythonError("ValueError", "'list' argument must have no negative elements");
        }
        const auto index = static_cast<std::size_t>(bin);
        if (index >= sums.size()) {
            sums.resize(index + 1, 0.0);
            counts.resize(index + 1, 0);
        }
        sums[index] += entries[k].value;
        ++counts[index];
    }

    const std::vector<std::int64_t> sizes = unit_sizes(intervals);
    std::vector<Distance> out;
    std::vector<std::size_t> zero_value_bins;
    std::size_t adjacent = 0;
    for (std::size_t k = 0; k < sums.size(); ++k) {
        Distance d;
        d.k = k;
        d.sum = sums[k];
        if (maxdepth_set && k == 0) {
            d.skipped = true;
            out.push_back(d);
            continue;
        }
        if (k == 0) {
            std::int64_t total_intra = unit_bins * unit_bins;
            for (const std::int64_t size : sizes) {
                total_intra -= size * size;
            }
            d.diagonal_length = static_cast<double>(total_intra) / 2.0;
            d.diagonal_is_float = true;
        } else {
            const auto offset_k = static_cast<std::int64_t>(k - 1);
            std::int64_t length = 0;
            for (const std::int64_t size : sizes) {
                if (size > offset_k) {
                    length += size - offset_k;
                }
            }
            d.diagonal_length = static_cast<double>(length);
        }
        if (static_cast<double>(counts[k]) > d.diagonal_length) {
            d.diagonal_length = static_cast<double>(counts[k]);
            d.diagonal_is_float = false;
        }
        if (d.diagonal_length != 0.0) {
            d.defined = true;
            d.mean = d.sum / d.diagonal_length;
            if (d.sum == 0.0) {
                d.zero_value = true;
                if (!zero_value_bins.empty() && k - zero_value_bins.back() == 1) {
                    ++adjacent;
                }
                zero_value_bins.push_back(k);
            }
            if (zero_value_bins.size() > 10 && adjacent > 10) {
                d.stop = true;
                out.push_back(d);
                break;
            }
        } else {
            d.mean = std::nan("");
        }
        d.nan_sum = std::isnan(d.sum);
        out.push_back(d);
    }
    return out;
}

}  // namespace hicx::dist_counts
