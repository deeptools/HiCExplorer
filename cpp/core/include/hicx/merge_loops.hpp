// hicexplorer.hicMergeLoops.mergeLoops for the public Python API.
//
// The Python builds two dictionaries of intervaltree IntervalTrees with
// hiCMatrix.intervalListToIntervalTree (one over pTupleX, one over pTupleY,
// each interval carrying its index in the list), walks the loops of the data
// frame removing the losers from both, and returns the ids that survive in
// both trees, chromosome pair by chromosome pair. The trees are
// hicx::pytree::Tree (hicx/py_intervaltree.hpp), so the order of every
// overlap set, which breaks ties between equally wide intervals, is
// CPython's. The dictionary lookups are the caller's (they follow Python's
// key equality): it hands over, per tuple, the tree it goes into, and per
// loop, which trees loop[0] and loop[3] select.
//
// Fixed in HiCExplorer 4, as in the tool (cpp/tools/hicMergeLoops.cpp): the y
// anchor is looked up in the y tree of loop[3] (3.7.6 used loop[0]), a loop
// whose chromosome has no tree is skipped (3.7.6 reused the intervals of the
// previous loop, or raised for the first), and the result lists the loops alive
// in an x and in a y tree (3.7.6 paired the trees positionally).

#ifndef HICX_MERGE_LOOPS_HPP
#define HICX_MERGE_LOOPS_HPP

#include <cstddef>
#include <cstdint>
#include <vector>

namespace hicx::merge_loops {

struct Intervals {
    // per interval (tuple): the tree it is added to, begin and end
    std::vector<std::int64_t> tree;
    std::vector<std::int64_t> begin;
    std::vector<std::int64_t> end;
    std::int64_t tree_count = 0;
    // the trees of the dictionary in its key order
    std::vector<std::int64_t> dictionary;
};

struct Loop {
    std::int64_t x_tree = -1;  // X[loop[0]], -1 when loop[0] is not a key of X
    std::int64_t y_tree = -1;  // Y[loop[3]], -1 when loop[3] is not a key of Y
    std::int64_t x_begin = 0;  // the overlap queries
    std::int64_t x_end = 0;
    std::int64_t y_begin = 0;
    std::int64_t y_end = 0;
};

struct Result {
    std::vector<std::int64_t> ids;
};

[[nodiscard]] Result merge(const Intervals& x, const Intervals& y, const std::vector<Loop>& loops);

}  // namespace hicx::merge_loops

#endif  // HICX_MERGE_LOOPS_HPP
