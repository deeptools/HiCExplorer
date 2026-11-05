#include "hicx/bedtools_ops.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <unordered_map>

namespace hicx::bedtools {

namespace {

struct StartAndRow {
    std::int64_t start = 0;
    std::size_t row = 0;
};

}  // namespace

std::vector<std::size_t> sort_order(const std::vector<std::string>& chrom,
                                    const std::vector<std::int64_t>& start) {
    // bedtools' BedFile::loadBedFileIntoMapNoBin, which appends every record to
    // the vector of its chromosome in file order, and whose map iterates in
    // lexicographic key order.
    std::map<std::string, std::vector<StartAndRow>> by_chrom;
    for (std::size_t row = 0; row < chrom.size(); ++row) {
        by_chrom[chrom[row]].push_back(StartAndRow{start[row], row});
    }

    std::vector<std::size_t> order;
    order.reserve(chrom.size());
    for (auto& entry : by_chrom) {
        // The comparator bedtools uses, start only. Not a total order, and
        // std::sort is not stable, which is exactly the point: see the header.
        std::sort(entry.second.begin(), entry.second.end(),
                  [](const StartAndRow& a, const StartAndRow& b) {
                      return a.start < b.start;
                  });
        for (const StartAndRow& record : entry.second) {
            order.push_back(record.row);
        }
    }
    return order;
}

Intervals merge(const Intervals& input) {
    Intervals output;
    std::size_t row = 0;
    while (row < input.size()) {
        const std::string& chrom = input.chrom[row];
        std::int64_t open_start = input.start[row];
        std::int64_t open_end = input.end[row];
        std::size_t next = row + 1;
        while (next < input.size() && input.chrom[next] == chrom) {
            if (input.start[next] < input.start[next - 1]) {
                throw std::runtime_error(
                    "ERROR: input file is not sorted by chromosome then by start "
                    "position; bedtools merge requires a sorted input");
            }
            if (input.start[next] <= open_end) {
                // -d 0 merges overlapping and book-ended intervals.
                open_end = std::max(open_end, input.end[next]);
            } else {
                output.push_back(chrom, open_start, open_end);
                open_start = input.start[next];
                open_end = input.end[next];
            }
            ++next;
        }
        output.push_back(chrom, open_start, open_end);
        row = next;
    }
    return output;
}

std::vector<std::int64_t> intersect_count(const Intervals& a, const Intervals& b) {
    struct ChromIndex {
        std::vector<std::int64_t> start;
        std::vector<std::int64_t> end;
        bool disjoint_sorted = true;
    };
    std::unordered_map<std::string, ChromIndex> index;
    for (std::size_t row = 0; row < b.size(); ++row) {
        ChromIndex& entry = index[b.chrom[row]];
        if (!entry.start.empty() &&
            (b.start[row] < entry.start.back() || b.start[row] < entry.end.back())) {
            entry.disjoint_sorted = false;
        }
        entry.start.push_back(b.start[row]);
        entry.end.push_back(b.end[row]);
    }

    std::vector<std::int64_t> counts(a.size(), 0);
    for (std::size_t row = 0; row < a.size(); ++row) {
        const auto found = index.find(a.chrom[row]);
        if (found == index.end()) {
            continue;
        }
        const ChromIndex& entry = found->second;
        const std::int64_t query_start = a.start[row];
        const std::int64_t query_end = a.end[row];
        if (entry.disjoint_sorted) {
            // Both start[] and end[] are ascending, so the intervals that
            // overlap [query_start, query_end) are a contiguous run.
            const std::size_t high = static_cast<std::size_t>(
                std::lower_bound(entry.start.begin(), entry.start.end(), query_end) -
                entry.start.begin());
            const std::size_t low = static_cast<std::size_t>(
                std::upper_bound(entry.end.begin(), entry.end.begin() +
                                                        static_cast<std::ptrdiff_t>(high),
                                 query_start) -
                entry.end.begin());
            counts[row] = high > low ? static_cast<std::int64_t>(high - low) : 0;
            continue;
        }
        std::int64_t count = 0;
        for (std::size_t k = 0; k < entry.start.size(); ++k) {
            if (entry.start[k] < query_end && query_start < entry.end[k]) {
                ++count;
            }
        }
        counts[row] = count;
    }
    return counts;
}

}  // namespace hicx::bedtools
