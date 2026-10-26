// CPython's set and intervaltree 3.1.0's IntervalTree, reproduced so that
// the iteration order of IntervalTree.overlap's result is CPython's.
//
// hicMergeLoops decides ties by that order (see cpp/tools/hicMergeLoops.cpp,
// property 3, for the measurements); the tool and the public Python API's
// hicexplorer.hicMergeLoops.mergeLoops (hicx/merge_loops.hpp) both use
// these classes. Moved here unchanged from cpp/tools/hicMergeLoops.cpp.

#ifndef HICX_PY_INTERVALTREE_HPP
#define HICX_PY_INTERVALTREE_HPP

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace hicx::pytree {

// CPython's hash of an int that fits in 64 bits: the value modulo the
// Mersenne prime 2^61 - 1, with the sign kept, and -1 replaced by -2.
[[nodiscard]] inline std::int64_t python_hash_of_int(std::int64_t value) {
    constexpr std::uint64_t kModulus = (1ULL << 61) - 1;
    const bool negative = value < 0;
    const std::uint64_t magnitude =
        negative ? static_cast<std::uint64_t>(0) - static_cast<std::uint64_t>(value)
                 : static_cast<std::uint64_t>(value);
    std::int64_t hash = static_cast<std::int64_t>(magnitude % kModulus);
    if (negative) {
        hash = -hash;
    }
    return hash == -1 ? -2 : hash;
}

// CPython's tuple hash, the xxHash based one introduced in 3.8, for a pair of
// ints. intervaltree's Interval.__hash__ is hash((begin, end)); the data field
// is left out of it.
[[nodiscard]] inline std::int64_t python_hash_of_pair(std::int64_t first, std::int64_t second) {
    constexpr std::uint64_t kPrime1 = 11400714785074694791ULL;
    constexpr std::uint64_t kPrime2 = 14029467366897019727ULL;
    constexpr std::uint64_t kPrime5 = 2870177450012600261ULL;
    std::uint64_t acc = kPrime5;
    const std::int64_t items[2] = {first, second};
    for (const std::int64_t item : items) {
        acc += static_cast<std::uint64_t>(python_hash_of_int(item)) * kPrime2;
        acc = (acc << 31) | (acc >> 33);
        acc *= kPrime1;
    }
    acc += 2ULL ^ (kPrime5 ^ 3527539ULL);
    if (acc == static_cast<std::uint64_t>(-1)) {
        return 1546275796;
    }
    return static_cast<std::int64_t>(acc);
}

// A CPython 3.12 set (Objects/setobject.c) of Interval keys, reproduced slot
// for slot so that it iterates in the same order: open addressing with nine
// linear probes then perturbation, dummies on deletion that a later insertion
// reuses (the last one passed), resizing at fill*5 >= mask*3, and the fast
// paths of set_merge. A key is an interval's index in its Tree; two keys are
// equal only when they are the same index, as Interval equality includes the
// data field and every id is distinct.
class PySet {
public:
    static constexpr std::int32_t kEmpty = -1;
    static constexpr std::int32_t kDummy = -2;

    struct Entry {
        std::int32_t key = kEmpty;
        std::int64_t hash = 0;
    };

    PySet() : table_(8) {}

    [[nodiscard]] std::size_t size() const { return used_; }
    [[nodiscard]] bool empty() const { return used_ == 0; }

    // Keys in iteration order.
    [[nodiscard]] std::vector<std::int32_t> items() const {
        std::vector<std::int32_t> out;
        out.reserve(used_);
        for (const Entry& entry : table_) {
            if (entry.key >= 0) {
                out.push_back(entry.key);
            }
        }
        return out;
    }

    // set_add_entry
    void add(std::int32_t key, std::int64_t hash) {
        const std::size_t mask = mask_;
        std::size_t perturb = static_cast<std::size_t>(hash);
        std::size_t i = static_cast<std::size_t>(hash) & mask;
        Entry* freeslot = nullptr;
        while (true) {
            std::size_t probes = (i + kLinearProbes <= mask) ? kLinearProbes : 0;
            std::size_t j = i;
            while (true) {
                Entry& entry = table_[j];
                if (entry.key == kEmpty) {
                    if (freeslot != nullptr) {
                        ++used_;
                        freeslot->key = key;
                        freeslot->hash = hash;
                        return;
                    }
                    ++fill_;
                    ++used_;
                    entry.key = key;
                    entry.hash = hash;
                    if (fill_ * 5 < mask * 3) {
                        return;
                    }
                    resize(used_ > 50000 ? used_ * 2 : used_ * 4);
                    return;
                }
                if (entry.key == key) {
                    return;
                }
                if (entry.key == kDummy) {
                    freeslot = &entry;
                }
                ++j;
                if (probes == 0) {
                    break;
                }
                --probes;
            }
            perturb >>= kPerturbShift;
            i = (i * 5 + 1 + perturb) & mask;
        }
    }

    [[nodiscard]] bool contains(std::int32_t key, std::int64_t hash) const {
        return find(key, hash) != kNotFound;
    }

    // set_discard_entry; returns whether the key was present.
    bool discard(std::int32_t key, std::int64_t hash) {
        const std::size_t slot = find(key, hash);
        if (slot == kNotFound) {
            return false;
        }
        table_[slot].key = kDummy;
        table_[slot].hash = -1;
        --used_;
        return true;
    }

    // set_merge (s.update(other) or set(other) for a set argument)
    void merge(const PySet& other) {
        if (&other == this || other.used_ == 0) {
            return;
        }
        if ((fill_ + other.used_) * 5 >= mask_ * 3) {
            resize((used_ + other.used_) * 2);
        }
        if (fill_ == 0 && mask_ == other.mask_ && other.fill_ == other.used_) {
            table_ = other.table_;
            fill_ = other.fill_;
            used_ = other.used_;
            return;
        }
        if (fill_ == 0) {
            fill_ = other.used_;
            used_ = other.used_;
            for (const Entry& entry : other.table_) {
                if (entry.key >= 0) {
                    insert_clean(table_, mask_, entry.key, entry.hash);
                }
            }
            return;
        }
        for (const Entry& entry : other.table_) {
            if (entry.key >= 0) {
                add(entry.key, entry.hash);
            }
        }
    }

    // set_difference_update_internal for a set argument that is a subset of
    // this one (the only use, `s_center -= child.s_center`), so the
    // intersection shortcut for a much larger argument never applies.
    void difference_update(const PySet& other) {
        for (const Entry& entry : other.table_) {
            if (entry.key >= 0) {
                discard(entry.key, entry.hash);
            }
        }
        if (fill_ - used_ <= mask_ / 4) {
            return;
        }
        resize(used_ > 50000 ? used_ * 2 : used_ * 4);
    }

private:
    static constexpr std::size_t kLinearProbes = 9;
    static constexpr int kPerturbShift = 5;
    static constexpr std::size_t kNotFound = static_cast<std::size_t>(-1);

    [[nodiscard]] std::size_t find(std::int32_t key, std::int64_t hash) const {
        const std::size_t mask = mask_;
        std::size_t perturb = static_cast<std::size_t>(hash);
        std::size_t i = static_cast<std::size_t>(hash) & mask;
        while (true) {
            std::size_t probes = (i + kLinearProbes <= mask) ? kLinearProbes : 0;
            std::size_t j = i;
            while (true) {
                const Entry& entry = table_[j];
                if (entry.key == kEmpty) {
                    return kNotFound;
                }
                if (entry.key == key) {
                    return j;
                }
                ++j;
                if (probes == 0) {
                    break;
                }
                --probes;
            }
            perturb >>= kPerturbShift;
            i = (i * 5 + 1 + perturb) & mask;
        }
    }

    static void insert_clean(std::vector<Entry>& table, std::size_t mask, std::int32_t key,
                             std::int64_t hash) {
        std::size_t perturb = static_cast<std::size_t>(hash);
        std::size_t i = static_cast<std::size_t>(hash) & mask;
        while (true) {
            if (table[i].key == kEmpty) {
                table[i].key = key;
                table[i].hash = hash;
                return;
            }
            if (i + kLinearProbes <= mask) {
                for (std::size_t j = 1; j <= kLinearProbes; ++j) {
                    if (table[i + j].key == kEmpty) {
                        table[i + j].key = key;
                        table[i + j].hash = hash;
                        return;
                    }
                }
            }
            perturb >>= kPerturbShift;
            i = (i * 5 + 1 + perturb) & mask;
        }
    }

    // set_table_resize
    void resize(std::size_t minused) {
        std::size_t newsize = 8;
        while (newsize <= minused) {
            newsize <<= 1;
        }
        if (newsize == 8 && mask_ == 7 && fill_ == used_) {
            return;  // the small table is kept as it is
        }
        std::vector<Entry> old = std::move(table_);
        table_.assign(newsize, Entry{});
        mask_ = newsize - 1;
        fill_ = used_;
        for (const Entry& entry : old) {
            if (entry.key >= 0) {
                insert_clean(table_, mask_, entry.key, entry.hash);
            }
        }
    }

    std::vector<Entry> table_;
    std::size_t mask_ = 7;
    std::size_t fill_ = 0;
    std::size_t used_ = 0;
};

// One chromosome's intervaltree 3.1.0 IntervalTree, reproduced node for node:
// the centered AVL tree of intervaltree/node.py (insertion, rotations with
// promotion of intervals to the new top, removal with pruning and
// pop_greatest_child) with every s_center a PySet, and the sorted boundary
// table. The tree's shape and the sets' slot layouts depend on the order of
// earlier insertions and removals, and IntervalTree.overlap returns a set
// whose iteration order decides ties in hicMergeLoops, so all of it is kept.
class Tree {
public:
    std::vector<std::int64_t> begin;
    std::vector<std::int64_t> end;
    std::vector<std::size_t> id;
    std::vector<std::int64_t> hash;
    std::vector<bool> present;

    // IntervalTree.add for a new interval (every id is distinct, so the
    // `interval in self` early return never fires).
    void add(std::int64_t b, std::int64_t e, std::size_t i) {
        const auto key = static_cast<std::int32_t>(begin.size());
        begin.push_back(b);
        end.push_back(e);
        id.push_back(i);
        hash.push_back(python_hash_of_pair(b, e));
        present.push_back(true);
        if (top_ == kNone) {
            top_ = from_interval(key);
        } else {
            top_ = node_add(top_, key);
        }
        ++boundaries_[b];
        ++boundaries_[e];
    }

    // IntervalTree.remove; false when the interval is not in the tree, where
    // the reference raises ValueError.
    bool remove(std::int32_t key) {
        if (!present[static_cast<std::size_t>(key)]) {
            return false;
        }
        bool done = false;
        top_ = remove_helper(top_, key, done);
        present[static_cast<std::size_t>(key)] = false;
        for (const std::int64_t boundary : {begin[static_cast<std::size_t>(key)],
                                            end[static_cast<std::size_t>(key)]}) {
            auto found = boundaries_.find(boundary);
            if (--found->second == 0) {
                boundaries_.erase(found);
            }
        }
        return true;
    }

    // IntervalTree.overlap(begin, end): the keys of the returned set in its
    // iteration order.
    [[nodiscard]] std::vector<std::int32_t> overlap(std::int64_t b, std::int64_t e) const {
        if (top_ == kNone || b >= e) {
            return {};
        }
        PySet result;
        search_point(top_, b, result);
        PySet more;
        for (auto it = boundaries_.lower_bound(b); it != boundaries_.end() && it->first < e;
             ++it) {
            search_point(top_, it->first, more);
        }
        result.merge(more);
        return result.items();
    }

private:
    static constexpr std::int32_t kNone = -1;

    struct Node {
        std::int64_t x_center = 0;
        PySet s_center;
        std::int32_t child[2] = {kNone, kNone};
        int depth = 0;
        int balance = 0;
    };

    std::vector<Node> nodes_;
    std::int32_t top_ = kNone;
    std::map<std::int64_t, int> boundaries_;

    Node& node(std::int32_t index) { return nodes_[static_cast<std::size_t>(index)]; }
    [[nodiscard]] const Node& node(std::int32_t index) const {
        return nodes_[static_cast<std::size_t>(index)];
    }
    [[nodiscard]] bool contains_point(std::int32_t key, std::int64_t point) const {
        const auto k = static_cast<std::size_t>(key);
        return begin[k] <= point && point < end[k];
    }
    void set_add(std::int32_t n, std::int32_t key) {
        node(n).s_center.add(key, hash[static_cast<std::size_t>(key)]);
    }

    // Node(x_center, s_center) for a node without children: the constructor's
    // rotate() only refreshes depth and balance.
    std::int32_t new_node(std::int64_t x_center) {
        nodes_.emplace_back();
        nodes_.back().x_center = x_center;
        const auto index = static_cast<std::int32_t>(nodes_.size() - 1);
        refresh_balance(index);
        return index;
    }

    std::int32_t from_interval(std::int32_t key) {
        const std::int32_t n = new_node(begin[static_cast<std::size_t>(key)]);
        set_add(n, key);
        return n;
    }

    void refresh_balance(std::int32_t n) {
        Node& self = node(n);
        const int left = self.child[0] != kNone ? node(self.child[0]).depth : 0;
        const int right = self.child[1] != kNone ? node(self.child[1]).depth : 0;
        self.depth = 1 + std::max(left, right);
        self.balance = right - left;
    }

    std::int32_t rotate(std::int32_t n) {
        refresh_balance(n);
        if (std::abs(node(n).balance) < 2) {
            return n;
        }
        const int my_heavy = node(n).balance > 0 ? 1 : 0;
        const std::int32_t heavy_child = node(n).child[my_heavy];
        const int child_heavy = node(heavy_child).balance > 0 ? 1 : 0;
        if (my_heavy == child_heavy || node(heavy_child).balance == 0) {
            return srotate(n);
        }
        return drotate(n);
    }

    std::int32_t srotate(std::int32_t n) {
        const int heavy = node(n).balance > 0 ? 1 : 0;
        const int light = 1 - heavy;
        const std::int32_t save = node(n).child[heavy];
        node(n).child[heavy] = node(save).child[light];
        const std::int32_t rotated = rotate(n);
        node(save).child[light] = rotated;

        // Promote the intervals of the new light child's s_center that overlap
        // save.x_center, in that set's iteration order.
        std::vector<std::int32_t> promotees;
        const std::int64_t center = node(save).x_center;
        for (const std::int32_t key : node(node(save).child[light]).s_center.items()) {
            if (contains_point(key, center)) {
                promotees.push_back(key);
            }
        }
        if (!promotees.empty()) {
            for (const std::int32_t key : promotees) {
                bool done = false;
                node(save).child[light] = remove_helper(node(save).child[light], key, done);
            }
            for (const std::int32_t key : promotees) {
                set_add(save, key);
            }
        }
        refresh_balance(save);
        return save;
    }

    std::int32_t drotate(std::int32_t n) {
        const int my_heavy = node(n).balance > 0 ? 1 : 0;
        node(n).child[my_heavy] = srotate(node(n).child[my_heavy]);
        refresh_balance(n);
        return srotate(n);
    }

    std::int32_t node_add(std::int32_t n, std::int32_t key) {
        if (contains_point(key, node(n).x_center)) {
            set_add(n, key);
            return n;
        }
        const int direction = begin[static_cast<std::size_t>(key)] > node(n).x_center ? 1 : 0;
        if (node(n).child[direction] == kNone) {
            const std::int32_t leaf = from_interval(key);
            node(n).child[direction] = leaf;
            refresh_balance(n);
            return n;
        }
        const std::int32_t sub = node_add(node(n).child[direction], key);
        node(n).child[direction] = sub;
        return rotate(n);
    }

    std::int32_t remove_helper(std::int32_t n, std::int32_t key, bool& done) {
        if (contains_point(key, node(n).x_center)) {
            if (!node(n).s_center.discard(key, hash[static_cast<std::size_t>(key)])) {
                throw std::runtime_error("intervaltree: interval missing from its node");
            }
            if (!node(n).s_center.empty()) {
                done = true;
                return n;
            }
            return prune(n);
        }
        const int direction = begin[static_cast<std::size_t>(key)] > node(n).x_center ? 1 : 0;
        if (node(n).child[direction] == kNone) {
            throw std::runtime_error("intervaltree: interval missing from the tree");
        }
        const std::int32_t sub = remove_helper(node(n).child[direction], key, done);
        node(n).child[direction] = sub;
        if (!done) {
            return rotate(n);
        }
        return n;
    }

    std::int32_t prune(std::int32_t n) {
        if (node(n).child[0] == kNone || node(n).child[1] == kNone) {
            const int direction = node(n).child[0] == kNone ? 1 : 0;
            return node(n).child[direction];
        }
        const auto [heir, left] = pop_greatest_child(node(n).child[0]);
        node(n).child[0] = left;
        node(heir).child[0] = node(n).child[0];
        node(heir).child[1] = node(n).child[1];
        refresh_balance(heir);
        return rotate(heir);
    }

    std::pair<std::int32_t, std::int32_t> pop_greatest_child(std::int32_t n) {
        if (node(n).child[1] == kNone) {
            const std::vector<std::int32_t> keys = node(n).s_center.items();
            std::int64_t max_end = std::numeric_limits<std::int64_t>::min();
            for (const std::int32_t key : keys) {
                max_end = std::max(max_end, end[static_cast<std::size_t>(key)]);
            }
            // sorted by (end, begin), pop the last, then every other interval
            // whose end differs from the popped one's raises new_x_center.
            std::int64_t new_x_center = node(n).x_center;
            for (const std::int32_t key : keys) {
                const std::int64_t e = end[static_cast<std::size_t>(key)];
                if (e == max_end) {
                    continue;
                }
                new_x_center = std::max(new_x_center, e);
            }
            const std::int32_t child = new_node(new_x_center);
            for (const std::int32_t key : node(n).s_center.items()) {
                if (contains_point(key, new_x_center)) {
                    set_add(child, key);
                }
            }
            // self.s_center -= child.s_center
            const PySet taken = node(child).s_center;
            node(n).s_center.difference_update(taken);
            if (!node(n).s_center.empty()) {
                return {child, n};
            }
            return {child, node(n).child[0]};
        }
        const auto [greatest, right] = pop_greatest_child(node(n).child[1]);
        node(n).child[1] = right;
        // for iv in set(self.s_center): move the ones overlapping the new
        // centre into greatest_child.
        PySet copy;
        copy.merge(node(n).s_center);
        const std::int64_t center = node(greatest).x_center;
        for (const std::int32_t key : copy.items()) {
            if (contains_point(key, center)) {
                node(n).s_center.discard(key, hash[static_cast<std::size_t>(key)]);
                set_add(greatest, key);  // center_hit holds, so Node.add stores it here
            }
        }
        if (!node(n).s_center.empty()) {
            refresh_balance(n);
            return {greatest, rotate(n)};
        }
        return {greatest, prune(n)};
    }

    void search_point(std::int32_t n, std::int64_t point, PySet& result) const {
        while (n != kNone) {
            const Node& self = node(n);
            for (const std::int32_t key : self.s_center.items()) {
                if (contains_point(key, point)) {
                    result.add(key, hash[static_cast<std::size_t>(key)]);
                }
            }
            if (point < self.x_center && self.child[0] != kNone) {
                n = self.child[0];
            } else if (point > self.x_center && self.child[1] != kNone) {
                n = self.child[1];
            } else {
                return;
            }
        }
    }
};

}  // namespace hicx::pytree

#endif  // HICX_PY_INTERVALTREE_HPP
