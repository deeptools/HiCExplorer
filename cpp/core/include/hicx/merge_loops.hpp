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
// Reproduced from the Python, as the tool does (cpp/tools/hicMergeLoops.cpp):
// the y tree is looked up under loop[0]; x_interval and y_interval survive
// across loops; the losers are removed from the trees of loop[0], which for a
// stale interval set is not the tree it came from.

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
    bool y_guard = false;      // loop[3] in Y
    std::int64_t y_tree = -1;  // Y[loop[0]], -1 when loop[0] is not a key of Y
    std::int64_t x_begin = 0;  // the overlap queries
    std::int64_t x_end = 0;
    std::int64_t y_begin = 0;
    std::int64_t y_end = 0;
};

enum class Error {
    None,
    UnboundX,         // len(x_interval) before any assignment
    UnboundY,         // len(y_interval) before any assignment
    MissingYTree,     // Y[loop[0]] for the overlap: KeyError(loop[0])
    MissingXRemove,   // X[loop[0]].remove: KeyError(loop[0])
    MissingYRemove,   // Y[loop[0]].remove: KeyError(loop[0])
    NotInXTree,       // IntervalTree.remove of an interval it does not hold: ValueError
    NotInYTree,
};

struct Result {
    Error error = Error::None;
    std::size_t loop = 0;  // the loop at which the error happened
    std::vector<std::int64_t> ids;
};

[[nodiscard]] Result merge(const Intervals& x, const Intervals& y, const std::vector<Loop>& loops);

}  // namespace hicx::merge_loops

#endif  // HICX_MERGE_LOOPS_HPP
