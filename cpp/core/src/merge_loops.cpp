// See hicx/merge_loops.hpp.

#include "hicx/merge_loops.hpp"

#include <algorithm>
#include <cstdlib>
#include <optional>
#include <tuple>
#include <unordered_set>

#include "hicx/py_intervaltree.hpp"

namespace hicx::merge_loops {

namespace {

using hicx::pytree::Tree;

std::vector<Tree> build(const Intervals& intervals) {
    std::vector<Tree> trees(static_cast<std::size_t>(intervals.tree_count));
    for (std::size_t k = 0; k < intervals.tree.size(); ++k) {
        trees[static_cast<std::size_t>(intervals.tree[k])].add(intervals.begin[k], intervals.end[k], k);
    }
    return trees;
}

// The intervals of an overlap set, in its iteration order, with the tree
// they came from.
struct Hits {
    std::int64_t tree = -1;
    std::vector<std::int32_t> keys;
};

}  // namespace

Result merge(const Intervals& x, const Intervals& y, const std::vector<Loop>& loops) {
    std::vector<Tree> trees_x = build(x);
    std::vector<Tree> trees_y = build(y);
    std::optional<Hits> x_interval;
    std::optional<Hits> y_interval;
    Result result;
    for (std::size_t r = 0; r < loops.size(); ++r) {
        const Loop& loop = loops[r];
        if (loop.x_tree >= 0) {
            x_interval = Hits{loop.x_tree, trees_x[static_cast<std::size_t>(loop.x_tree)].overlap(loop.x_begin,
                                                                                                    loop.x_end)};
        }
        if (loop.y_guard) {
            if (loop.y_tree < 0) {
                result.error = Error::MissingYTree;
                result.loop = r;
                return result;
            }
            y_interval = Hits{loop.y_tree, trees_y[static_cast<std::size_t>(loop.y_tree)].overlap(loop.y_begin,
                                                                                                    loop.y_end)};
        }
        if (!x_interval) {
            result.error = Error::UnboundX;
            result.loop = r;
            return result;
        }
        if (x_interval->keys.size() <= 1) {
            continue;
        }
        if (!y_interval) {
            result.error = Error::UnboundY;
            result.loop = r;
            return result;
        }
        if (y_interval->keys.size() <= 1) {
            continue;
        }
        const Tree& tx = trees_x[static_cast<std::size_t>(x_interval->tree)];
        const Tree& ty = trees_y[static_cast<std::size_t>(y_interval->tree)];
        // snapshots: the removals below change the trees
        struct Item {
            std::int32_t key;
            std::size_t id;
        };
        std::vector<Item> xs;
        std::vector<Item> ys;
        std::unordered_set<std::size_t> interest_x;
        for (const std::int32_t key : x_interval->keys) {
            xs.push_back({key, tx.id[static_cast<std::size_t>(key)]});
            interest_x.insert(tx.id[static_cast<std::size_t>(key)]);
        }
        std::size_t max_index = 0;
        std::int64_t max_distance = 0;
        std::unordered_set<std::size_t> all_ids;
        for (const std::int32_t key : y_interval->keys) {
            const std::size_t slot = static_cast<std::size_t>(key);
            ys.push_back({key, ty.id[slot]});
            if (interest_x.find(ty.id[slot]) == interest_x.end()) {
                continue;
            }
            const std::int64_t width = std::llabs(ty.begin[slot] - ty.end[slot]);
            if (width > max_distance) {
                max_distance = width;
                max_index = ty.id[slot];
            }
            all_ids.insert(ty.id[slot]);
        }
        for (const Item& item : xs) {
            if (item.id == max_index || all_ids.find(item.id) == all_ids.end()) {
                continue;
            }
            if (loop.x_tree < 0) {
                result.error = Error::MissingXRemove;
                result.loop = r;
                return result;
            }
            if (loop.x_tree != x_interval->tree ||
                !trees_x[static_cast<std::size_t>(loop.x_tree)].remove(item.key)) {
                result.error = Error::NotInXTree;
                result.loop = r;
                return result;
            }
        }
        for (const Item& item : ys) {
            if (item.id == max_index || all_ids.find(item.id) == all_ids.end()) {
                continue;
            }
            if (loop.y_tree < 0) {
                result.error = Error::MissingYRemove;
                result.loop = r;
                return result;
            }
            if (loop.y_tree != y_interval->tree ||
                !trees_y[static_cast<std::size_t>(loop.y_tree)].remove(item.key)) {
                result.error = Error::NotInYTree;
                result.loop = r;
                return result;
            }
        }
    }
    // zip(X, Y) over the dictionaries in key order; the ids of the sorted x
    // intervals that are also in the y tree
    const std::size_t pairs = std::min(x.dictionary.size(), y.dictionary.size());
    for (std::size_t k = 0; k < pairs; ++k) {
        const Tree& tx = trees_x[static_cast<std::size_t>(x.dictionary[k])];
        const Tree& ty = trees_y[static_cast<std::size_t>(y.dictionary[k])];
        std::unordered_set<std::size_t> alive_y;
        for (std::size_t slot = 0; slot < ty.id.size(); ++slot) {
            if (ty.present[slot]) {
                alive_y.insert(ty.id[slot]);
            }
        }
        std::vector<std::tuple<std::int64_t, std::int64_t, std::size_t>> sorted_x;
        for (std::size_t slot = 0; slot < tx.id.size(); ++slot) {
            if (tx.present[slot]) {
                sorted_x.emplace_back(tx.begin[slot], tx.end[slot], tx.id[slot]);
            }
        }
        std::sort(sorted_x.begin(), sorted_x.end());
        for (const auto& [begin, end, id] : sorted_x) {
            if (alive_y.find(id) != alive_y.end()) {
                result.ids.push_back(static_cast<std::int64_t>(id));
            }
        }
    }
    return result;
}

}  // namespace hicx::merge_loops
