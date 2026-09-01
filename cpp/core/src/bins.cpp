#include "hicx/bins.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace hicx {

BinTable::BinTable(std::vector<CutInterval> intervals)
    : intervals_(std::move(intervals)) {
    build_index();
}

void BinTable::build_index() {
    boundaries_.clear();
    boundary_index_.clear();
    trees_.clear();
    if (intervals_.empty()) {
        return;
    }

    // Literal transcription of intervalListToIntervalTree.
    std::int64_t intval_id = 0;
    std::int64_t chr_start_id = 0;
    const std::string* previous_chrom = nullptr;
    std::string current_chrom;

    auto set_boundary = [this](const std::string& chrom, BinRange range) {
        auto it = boundary_index_.find(chrom);
        if (it == boundary_index_.end()) {
            boundary_index_.emplace(chrom, boundaries_.size());
            boundaries_.emplace_back(chrom, range);
        } else {
            boundaries_[it->second].second = range;
        }
    };

    for (const CutInterval& interval : intervals_) {
        if (previous_chrom == nullptr || *previous_chrom != interval.chrom) {
            if (previous_chrom == nullptr) {
                current_chrom = interval.chrom;
                previous_chrom = &current_chrom;
            }
            set_boundary(*previous_chrom, BinRange{chr_start_id, intval_id});
            chr_start_id = intval_id;
            // A repeated chromosome resets its tree, exactly as the Python
            // implementation does with cut_int_tree[chrom] = IntervalTree().
            trees_[interval.chrom].clear();
            current_chrom = interval.chrom;
            previous_chrom = &current_chrom;
        }
        trees_[interval.chrom].push_back(
            BinInterval{interval.start, interval.end, intval_id});
        ++intval_id;
    }
    set_boundary(*previous_chrom, BinRange{chr_start_id, intval_id});

    for (auto& [chrom, tree] : trees_) {
        (void)chrom;
        std::sort(tree.begin(), tree.end(),
                  [](const BinInterval& a, const BinInterval& b) {
                      if (a.start != b.start) {
                          return a.start < b.start;
                      }
                      if (a.end != b.end) {
                          return a.end < b.end;
                      }
                      return a.bin_id < b.bin_id;
                  });
    }
}

const CutInterval& BinTable::bin_pos(std::size_t bin_id) const {
    if (bin_id >= intervals_.size()) {
        throw std::out_of_range("binIndex: " + std::to_string(bin_id) + " not found");
    }
    return intervals_[bin_id];
}

std::vector<std::string> BinTable::chrom_names() const {
    std::vector<std::string> names;
    names.reserve(boundaries_.size());
    for (const auto& [chrom, range] : boundaries_) {
        (void)range;
        names.push_back(chrom);
    }
    return names;
}

std::optional<BinRange> BinTable::chrom_bin_range(const std::string& chrom) const {
    const auto it = boundary_index_.find(chrom);
    if (it == boundary_index_.end()) {
        return std::nullopt;
    }
    return boundaries_[it->second].second;
}

std::vector<std::pair<std::string, std::int64_t>> BinTable::chromosome_sizes() const {
    std::vector<std::pair<std::string, std::int64_t>> sizes;
    sizes.reserve(boundaries_.size());
    for (const auto& [chrom, range] : boundaries_) {
        if (range.last <= 0 || static_cast<std::size_t>(range.last) > intervals_.size()) {
            continue;
        }
        const CutInterval& last = intervals_[static_cast<std::size_t>(range.last) - 1];
        sizes.emplace_back(last.chrom, last.end);
    }
    return sizes;
}

std::vector<std::pair<std::string, std::int64_t>> BinTable::chromosome_sizes_real() const {
    std::vector<std::pair<std::string, std::int64_t>> sizes;
    sizes.reserve(boundaries_.size());
    for (const auto& [chrom, range] : boundaries_) {
        if (range.last <= 0 || static_cast<std::size_t>(range.last) > intervals_.size()) {
            continue;
        }
        const CutInterval& first = intervals_[static_cast<std::size_t>(range.first)];
        const CutInterval& last = intervals_[static_cast<std::size_t>(range.last) - 1];
        sizes.emplace_back(last.chrom, last.end - first.start + 1);
    }
    return sizes;
}

std::optional<std::int64_t> BinTable::bin_at(const std::string& chrom,
                                             std::int64_t position) const {
    const auto it = trees_.find(chrom);
    if (it == trees_.end()) {
        return std::nullopt;
    }
    const std::vector<BinInterval>& tree = it->second;
    // Bins are non overlapping, so the first interval with start <= position
    // and end > position is the answer. sorted(...)[0] in Python picks the
    // smallest start among the overlapping intervals, which is the same.
    auto upper = std::upper_bound(tree.begin(), tree.end(), position,
                                  [](std::int64_t value, const BinInterval& item) {
                                      return value < item.start;
                                  });
    for (auto candidate = tree.begin(); candidate != upper; ++candidate) {
        if (candidate->start <= position && position < candidate->end) {
            return candidate->bin_id;
        }
    }
    return std::nullopt;
}

std::optional<std::pair<std::int64_t, std::int64_t>> BinTable::region_bin_range(
    const std::string& chrom, std::int64_t startpos, std::int64_t endpos) const {
    const std::optional<std::int64_t> start_bin = bin_at(chrom, startpos);
    const std::optional<std::int64_t> end_bin = bin_at(chrom, endpos);
    if (!start_bin.has_value() || !end_bin.has_value()) {
        return std::nullopt;
    }
    return std::make_pair(*start_bin, *end_bin);
}

std::int64_t BinTable::bin_size() const {
    if (bin_size_.has_value()) {
        return *bin_size_;
    }
    if (intervals_.empty()) {
        throw std::runtime_error("cannot compute the bin size of an empty bin table");
    }
    if (intervals_.size() == 1) {
        bin_size_ = intervals_[0].end - intervals_[0].start;
        bin_size_homogeneous_ = true;
        return *bin_size_;
    }

    // Starts grouped by chromosome, in first appearance order, over the whole
    // list (not only over contiguous blocks), as the Python comprehension does.
    std::vector<std::string> order;
    std::unordered_map<std::string, std::vector<std::int64_t>> starts;
    for (const CutInterval& interval : intervals_) {
        auto [it, inserted] = starts.try_emplace(interval.chrom);
        if (inserted) {
            order.push_back(interval.chrom);
        }
        it->second.push_back(interval.start);
    }

    std::vector<std::int64_t> diffs;
    for (const std::string& chrom : order) {
        const std::vector<std::int64_t>& values = starts[chrom];
        if (values.size() < 2) {
            continue;
        }
        for (std::size_t i = 1; i < values.size(); ++i) {
            diffs.push_back(values[i] - values[i - 1]);
        }
    }
    if (diffs.empty()) {
        throw std::runtime_error(
            "cannot compute the bin size: every chromosome holds a single bin");
    }

    std::sort(diffs.begin(), diffs.end());
    const std::size_t n = diffs.size();
    double median = 0.0;
    if (n % 2 == 1) {
        median = static_cast<double>(diffs[n / 2]);
    } else {
        // np.median averages the two central values in float64.
        median = (static_cast<double>(diffs[n / 2 - 1]) +
                  static_cast<double>(diffs[n / 2])) /
                 2.0;
    }
    const std::int64_t median_int = static_cast<std::int64_t>(median);  // int() truncates

    std::size_t deviating = 0;
    for (const CutInterval& interval : intervals_) {
        if ((interval.end - interval.start) != median_int) {
            ++deviating;
        }
    }
    bin_size_homogeneous_ =
        !(static_cast<double>(deviating) >
          static_cast<double>(intervals_.size()) * 0.01);
    bin_size_ = median_int;
    return *bin_size_;
}

bool BinTable::bin_size_homogeneous() const {
    (void)bin_size();
    return bin_size_homogeneous_;
}

}  // namespace hicx
