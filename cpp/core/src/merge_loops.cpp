// See hicx/merge_loops.hpp.

#include "hicx/merge_loops.hpp"

#include <algorithm>
#include <cstdlib>
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

}  // namespace

Result merge(const Intervals& x, const Intervals& y, const std::vector<Loop>& loops) {
    std::vector<Tree> trees_x = build(x);
    std::vector<Tree> trees_y = build(y);
    Result result;
    for (std::size_t r = 0; r < loops.size(); ++r) {
        const Loop& loop = loops[r];
        // HiCExplorer 4: the x anchor is looked up in the x tree of loop[0],
        // the y anchor in the y tree of loop[3], and a loop whose chromosome
        // has no tree is skipped. 3.7.6 looked the y tree up under loop[0] and
        // evaluated such a loop against the intervals of the previous one.
        if (loop.x_tree < 0 || loop.y_tree < 0) {
            continue;
        }
        const std::size_t x_index = static_cast<std::size_t>(loop.x_tree);
        const std::size_t y_index = static_cast<std::size_t>(loop.y_tree);
        const std::vector<std::int32_t> x_keys = trees_x[x_index].overlap(loop.x_begin, loop.x_end);
        const std::vector<std::int32_t> y_keys = trees_y[y_index].overlap(loop.y_begin, loop.y_end);
        if (x_keys.size() <= 1 || y_keys.size() <= 1) {
            continue;
        }
        const Tree& tx = trees_x[x_index];
        const Tree& ty = trees_y[y_index];
        // snapshots: the removals below change the trees
        struct Item {
            std::int32_t key;
            std::size_t id;
        };
        std::vector<Item> xs;
        std::vector<Item> ys;
        std::unordered_set<std::size_t> interest_x;
        for (const std::int32_t key : x_keys) {
            xs.push_back({key, tx.id[static_cast<std::size_t>(key)]});
            interest_x.insert(tx.id[static_cast<std::size_t>(key)]);
        }
        std::size_t max_index = 0;
        std::int64_t max_distance = 0;
        std::unordered_set<std::size_t> all_ids;
        for (const std::int32_t key : y_keys) {
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
            trees_x[x_index].remove(item.key);
        }
        for (const Item& item : ys) {
            if (item.id == max_index || all_ids.find(item.id) == all_ids.end()) {
                continue;
            }
            trees_y[y_index].remove(item.key);
        }
    }
    // The ids of the x intervals, tree by tree in the dictionary's key order
    // and sorted within a tree, that are also alive in a y tree. 3.7.6 paired
    // the x and y trees positionally (zip), which matches the tree of a loop's
    // x chromosome with another chromosome's y tree as soon as a loop is
    // inter-chromosomal; HiCExplorer 4 asks whether the loop is alive in both.
    std::unordered_set<std::size_t> alive_y;
    for (const std::int64_t index : y.dictionary) {
        const Tree& ty = trees_y[static_cast<std::size_t>(index)];
        for (std::size_t slot = 0; slot < ty.id.size(); ++slot) {
            if (ty.present[slot]) {
                alive_y.insert(ty.id[slot]);
            }
        }
    }
    for (const std::int64_t index : x.dictionary) {
        const Tree& tx = trees_x[static_cast<std::size_t>(index)];
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
