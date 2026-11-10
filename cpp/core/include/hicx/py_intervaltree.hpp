// A centered interval tree whose overlap queries return their hits in the
// order of the Python package intervaltree 3.1.0.
//
// hicMergeLoops picks among equally wide anchors by the iteration order of the
// set that IntervalTree.overlap returns (cpp/tools/hicMergeLoops.cpp, property
// 3), and that order depends on the exact shape of the tree and on the slot
// layout of the sets held in its nodes. The tree is therefore built the way
// the package documents it: a centered interval tree (Edelsbrunner 1980;
// Berg et al., Computational Geometry, section 10.1) in which every node keeps
// the intervals that contain its centre in a Python set, kept height balanced
// as an AVL tree (Adelson-Velsky and Landis 1962), with intervals promoted to
// the new top node when a rotation moves their centre. Used by
// core/src/merge_loops.cpp, tools/hicMergeLoops.cpp and the Python API.
//
// Identity with the package was verified on 56 random insert / remove /
// overlap sequences (about 5000 queries, ties and repeated coordinates
// included) against intervaltree 3.1.0 on Python 3.12.7:
// cpp/tests/test_cpython_oracles.cpp replays
// cpp/tests/fixtures/cpython/intervaltree.txt, written by
// cpp/scripts/gen_intervaltree_fixture.py.

#ifndef HICX_PY_INTERVALTREE_HPP
#define HICX_PY_INTERVALTREE_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include "hicx/py_set_order.hpp"

namespace hicx::pytree {

// The intervals of one chromosome. An interval is identified by its slot, the
// order in which it was added; its element identity in the Python set is
// (begin, end, slot), which hashes as the tuple (begin, end).
class Tree {
public:
    std::vector<std::int64_t> begin;
    std::vector<std::int64_t> end;
    std::vector<std::size_t> id;
    std::vector<bool> present;

    // Adds [b, e) with the caller's identifier i.
    void add(std::int64_t b, std::int64_t e, std::size_t i) {
        const auto slot = static_cast<std::int32_t>(begin.size());
        begin.push_back(b);
        end.push_back(e);
        id.push_back(i);
        hash_.push_back(pyorder::hash_of_int_pair(b, e));
        present.push_back(true);
        root_ = (root_ == kNil) ? make_node(slot) : insert(root_, slot);
        ++endpoint_count_[b];
        ++endpoint_count_[e];
    }

    // Removes the interval in slot; false when it is not in the tree.
    bool remove(std::int32_t slot) {
        if (!present[idx(slot)]) {
            return false;
        }
        bool stable = false;
        root_ = erase(root_, slot, stable);
        present[idx(slot)] = false;
        for (const std::int64_t endpoint : {begin[idx(slot)], end[idx(slot)]}) {
            auto found = endpoint_count_.find(endpoint);
            if (--found->second == 0) {
                endpoint_count_.erase(found);
            }
        }
        return true;
    }

    // The slots of the intervals overlapping [b, e), in the iteration order of
    // the set the package returns: the hits at b first, then the hits at every
    // interval endpoint in [b, e) in ascending order, merged into that set.
    [[nodiscard]] std::vector<std::int32_t> overlap(std::int64_t b, std::int64_t e) const {
        if (root_ == kNil || b >= e) {
            return {};
        }
        pyorder::OrderedSet hits;
        collect_at(b, hits);
        pyorder::OrderedSet inner;
        for (auto it = endpoint_count_.lower_bound(b); it != endpoint_count_.end() && it->first < e;
             ++it) {
            collect_at(it->first, inner);
        }
        hits.update(inner);
        std::vector<std::int32_t> out;
        for (const std::int64_t key : hits.items()) {
            out.push_back(static_cast<std::int32_t>(key));
        }
        return out;
    }

private:
    static constexpr std::int32_t kNil = -1;
    enum Side { Left = 0, Right = 1 };

    struct Node {
        std::int64_t center = 0;
        pyorder::OrderedSet members;  // the intervals that contain center
        std::array<std::int32_t, 2> link = {kNil, kNil};
        int height = 0;
        int skew = 0;  // height(right) - height(left)
    };

    std::vector<Node> nodes_;
    std::int32_t root_ = kNil;
    std::vector<std::int64_t> hash_;
    std::map<std::int64_t, int> endpoint_count_;  // multiset of all endpoints

    static std::size_t idx(std::int32_t i) { return static_cast<std::size_t>(i); }
    Node& at(std::int32_t n) { return nodes_[idx(n)]; }
    [[nodiscard]] const Node& at(std::int32_t n) const { return nodes_[idx(n)]; }
    [[nodiscard]] bool holds_point(std::int32_t slot, std::int64_t point) const {
        return begin[idx(slot)] <= point && point < end[idx(slot)];
    }
    [[nodiscard]] Side side_of(std::int32_t slot, std::int32_t n) const {
        return begin[idx(slot)] > at(n).center ? Right : Left;
    }
    void put(std::int32_t n, std::int32_t slot) {
        at(n).members.add(slot, hash_[idx(slot)]);
    }
    void take(std::int32_t n, std::int32_t slot) {
        at(n).members.discard(slot, hash_[idx(slot)]);
    }

    std::int32_t new_node(std::int64_t center) {
        nodes_.emplace_back();
        nodes_.back().center = center;
        const auto n = static_cast<std::int32_t>(nodes_.size() - 1);
        settle(n);
        return n;
    }

    // A leaf centred on the begin of the interval, holding it.
    std::int32_t make_node(std::int32_t slot) {
        const std::int32_t n = new_node(begin[idx(slot)]);
        put(n, slot);
        return n;
    }

    // Recomputes height and skew from the children.
    void settle(std::int32_t n) {
        Node& node = at(n);
        const int left = node.link[Left] != kNil ? at(node.link[Left]).height : 0;
        const int right = node.link[Right] != kNil ? at(node.link[Right]).height : 0;
        node.height = 1 + std::max(left, right);
        node.skew = right - left;
    }

    // Restores the AVL condition at n; returns the new top of the subtree.
    std::int32_t rebalance(std::int32_t n) {
        settle(n);
        if (std::abs(at(n).skew) < 2) {
            return n;
        }
        const Side heavy = at(n).skew > 0 ? Right : Left;
        const std::int32_t child = at(n).link[heavy];
        const Side child_heavy = at(child).skew > 0 ? Right : Left;
        if (heavy == child_heavy || at(child).skew == 0) {
            return rotate_once(n);
        }
        return rotate_twice(n);
    }

    // Single rotation towards the light side. Intervals of the demoted
    // subtree that contain the new top's centre move up into the new top, in
    // the iteration order of the set they leave.
    std::int32_t rotate_once(std::int32_t n) {
        const Side heavy = at(n).skew > 0 ? Right : Left;
        const Side light = heavy == Right ? Left : Right;
        const std::int32_t top = at(n).link[heavy];
        at(n).link[heavy] = at(top).link[light];
        const std::int32_t demoted = rebalance(n);
        at(top).link[light] = demoted;

        const std::int64_t centre = at(top).center;
        std::vector<std::int32_t> promoted;
        for (const std::int64_t key : at(at(top).link[light]).members.items()) {
            if (holds_point(static_cast<std::int32_t>(key), centre)) {
                promoted.push_back(static_cast<std::int32_t>(key));
            }
        }
        for (const std::int32_t slot : promoted) {
            bool stable = false;
            at(top).link[light] = erase(at(top).link[light], slot, stable);
        }
        for (const std::int32_t slot : promoted) {
            put(top, slot);
        }
        settle(top);
        return top;
    }

    std::int32_t rotate_twice(std::int32_t n) {
        const Side heavy = at(n).skew > 0 ? Right : Left;
        at(n).link[heavy] = rotate_once(at(n).link[heavy]);
        settle(n);
        return rotate_once(n);
    }

    std::int32_t insert(std::int32_t n, std::int32_t slot) {
        if (holds_point(slot, at(n).center)) {
            put(n, slot);
            return n;
        }
        const Side side = side_of(slot, n);
        if (at(n).link[side] == kNil) {
            at(n).link[side] = make_node(slot);
            settle(n);
            return n;
        }
        at(n).link[side] = insert(at(n).link[side], slot);
        return rebalance(n);
    }

    // Removes slot from the subtree at n; `stable` is set once the removal
    // left a node non-empty, from where on no node above is rebalanced.
    std::int32_t erase(std::int32_t n, std::int32_t slot, bool& stable) {
        if (holds_point(slot, at(n).center)) {
            if (!at(n).members.discard(slot, hash_[idx(slot)])) {
                throw std::runtime_error("intervaltree: interval missing from its node");
            }
            if (!at(n).members.empty()) {
                stable = true;
                return n;
            }
            return collapse(n);
        }
        const Side side = side_of(slot, n);
        if (at(n).link[side] == kNil) {
            throw std::runtime_error("intervaltree: interval missing from the tree");
        }
        at(n).link[side] = erase(at(n).link[side], slot, stable);
        return stable ? n : rebalance(n);
    }

    // The subtree that replaces node n after its last interval left.
    std::int32_t collapse(std::int32_t n) {
        if (at(n).link[Left] == kNil || at(n).link[Right] == kNil) {
            return at(n).link[at(n).link[Left] == kNil ? Right : Left];
        }
        // Two children: the rightmost node of the left subtree takes n's place.
        const auto [heir, left_rest] = detach_rightmost(at(n).link[Left]);
        at(n).link[Left] = left_rest;
        at(heir).link = at(n).link;
        settle(heir);
        return rebalance(heir);
    }

    // Detaches the rightmost node of the subtree at n and returns it together
    // with the root of what remains. The detached node is a new node: its
    // centre is the largest interval end below the largest one of that node,
    // and it takes over the intervals that contain this centre.
    std::pair<std::int32_t, std::int32_t> detach_rightmost(std::int32_t n) {
        if (at(n).link[Right] == kNil) {
            const std::vector<std::int64_t> keys = at(n).members.items();
            std::int64_t top_end = std::numeric_limits<std::int64_t>::min();
            for (const std::int64_t key : keys) {
                top_end = std::max(top_end, end[idx(static_cast<std::int32_t>(key))]);
            }
            std::int64_t centre = at(n).center;
            for (const std::int64_t key : keys) {
                const std::int64_t e = end[idx(static_cast<std::int32_t>(key))];
                if (e != top_end) {
                    centre = std::max(centre, e);
                }
            }
            const std::int32_t detached = new_node(centre);
            for (const std::int64_t key : at(n).members.items()) {
                if (holds_point(static_cast<std::int32_t>(key), centre)) {
                    put(detached, static_cast<std::int32_t>(key));
                }
            }
            const pyorder::OrderedSet moved = at(detached).members;
            at(n).members.difference_update(moved);
            if (!at(n).members.empty()) {
                return {detached, n};
            }
            return {detached, at(n).link[Left]};
        }
        const auto [detached, right_rest] = detach_rightmost(at(n).link[Right]);
        at(n).link[Right] = right_rest;
        // Intervals of n that contain the detached node's centre move there;
        // they are visited in the order of a copy of n's set.
        pyorder::OrderedSet snapshot;
        snapshot.update(at(n).members);
        const std::int64_t centre = at(detached).center;
        for (const std::int64_t key : snapshot.items()) {
            const auto slot = static_cast<std::int32_t>(key);
            if (holds_point(slot, centre)) {
                take(n, slot);
                put(detached, slot);
            }
        }
        if (!at(n).members.empty()) {
            settle(n);
            return {detached, rebalance(n)};
        }
        return {detached, collapse(n)};
    }

    // Adds to out the intervals that contain point, walking from the root
    // towards the node whose centre is point.
    void collect_at(std::int64_t point, pyorder::OrderedSet& out) const {
        std::int32_t n = root_;
        while (n != kNil) {
            const Node& node = at(n);
            for (const std::int64_t key : node.members.items()) {
                if (holds_point(static_cast<std::int32_t>(key), point)) {
                    out.add(key, hash_[idx(static_cast<std::int32_t>(key))]);
                }
            }
            if (point < node.center && node.link[Left] != kNil) {
                n = node.link[Left];
            } else if (point > node.center && node.link[Right] != kNil) {
                n = node.link[Right];
            } else {
                return;
            }
        }
    }
};

}  // namespace hicx::pytree

#endif  // HICX_PY_INTERVALTREE_HPP
