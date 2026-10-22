// Port of hicexplorer/hicMergeLoops.py.
//
// Merges loop calls made at several resolutions: a loop is dropped when it
// overlaps, in both its x and its y anchor, a loop of a lower resolution, and
// the lower-resolution call is the one kept.
//
// The reference builds two intervaltree dictionaries, one over the x anchors
// and one over the y anchors, both keyed by the *x* chromosome, and then walks
// the loops removing the losers from both trees. Four properties of that code
// are behaviour and are reproduced here rather than repaired. All four are
// pinned by the characterization tests added to
// hicexplorer/test/general/test_hicMergeLoops.py.
//
//  1. **The y tree is indexed by the x chromosome** (hicMergeLoops.py:81:
//     `target_regions_intervaltree_y[loop[0]]`, where loop[3] is the y
//     chromosome). Every loop in the corpus is intra-chromosomal so the two
//     agree, but on an inter-chromosomal call the lookup either finds the
//     wrong chromosome's tree or raises KeyError, since the guard on line 80
//     tests loop[3] while line 81 uses loop[0].
//  2. **x_interval and y_interval survive across iterations.** They are only
//     assigned inside `if loop[0] in ...` (":77-82"), so a loop on a
//     chromosome that is not in the tree is evaluated against the previous
//     loop's intervals, and a first loop in that state raises TypeError on
//     len(None).
//  3. **The winner among equally wide overlapping anchors is decided by
//     CPython's set iteration order**, and the port reproduces that order
//     exactly. `IntervalTree.overlap` returns a set, the scan at ":98-102"
//     keeps the first interval whose width is *strictly* greater than the
//     best so far, so equally wide intervals (every overlap when all inputs
//     share one resolution) break the tie by the order that set iterates in.
//     That order is deterministic: Interval hashes as `hash((begin, end))`,
//     a tuple of ints, and neither int nor tuple-of-int hashes are salted by
//     PYTHONHASHSEED (3.7.6 gives byte identical output under seeds 0, 1, 7,
//     12345 and random on ext.hicMergeLoops.gse234292_wt_ko.r25000). The
//     chromosome names are str, but they are only dict keys, and dicts
//     iterate in insertion order. The order does depend on the whole history
//     of the trees: the slot layout of a set depends on the order of its
//     insertions and deletions, and the order of insertion into the overlap
//     set depends on the shape of the tree and on the layout of every node's
//     s_center set, which earlier removals (rotations, promotions, pruning)
//     changed. So `Tree` below is intervaltree 3.1.0's node.py reproduced
//     node for node, and `PySet` is CPython 3.12's Objects/setobject.c
//     reproduced slot for slot (linear probing, perturbation, dummy reuse,
//     resize thresholds, the set_merge fast paths). Verified byte identical
//     against 3.7.6 on every case and on 400 random multi-resolution inputs
//     with many ties and repeated coordinates; `PySet` alone was checked
//     against real CPython sets over random add, discard, update, copy and
//     difference sequences. An earlier version ordered the hits by
//     `(hash & mask, id)` and diverged whenever two tied intervals collided.
//  4. **drop_duplicates(keep=False)** removes every copy of a repeated row, not
//     all but one (":165"), so a loop called identically at two resolutions
//     disappears from the merge entirely.
//
// Threading and SIMD: none. The tool is 21,524 loops against a few thousand
// intervals per chromosome and the C++ run costs 0.03 s of CPU against the
// Python's 15 s. cpp/OPTIMIZATION.md 6 requires a measurement first, and the
// measurement leaves nothing to gain.

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/bedtools_ops.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/text_table.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicMergeLoops --inputFiles INPUTFILES [INPUTFILES ...] --outFileName\n"
    "                     OUTFILENAME --lowestResolution LOWESTRESOLUTION [--help]\n"
    "                     [--version]\n";

const char* const kHelp =
    "\n"
    "This script merges the locations of loops detected at several resolutions.\n"
    "\n"
    "Loops in the inputFiles need to have the following format:\n"
    "\n"
    "chr start end chr start end\n"
    "\n"
    "Loops are merged if the x and y position of a loop overlap with the x and y "
    "position of another loop; all loops are considered as an overlap within +/- "
    "the bin size of the lowest resolution.\n"
    "I.e. for a loop with coordinates x and y, the overlap with all other loops is "
    "checked for (x - lowest resolution) and (y + lowest resolution).\n"
    "If two or more locations are to be merged, the loop at the lowest resolution "
    "is taken as the merged loop.\n"
    "\n"
    "Required arguments:\n"
    "  --inputFiles INPUTFILES [INPUTFILES ...], -i INPUTFILES [INPUTFILES ...]\n"
    "                        The loop files from hicDetectLoops. To use files from\n"
    "                        other sources, please follow 'chr start end chr start\n"
    "                        end' format and remove any header.\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        The name of the merged loop file.\n"
    "  --lowestResolution LOWESTRESOLUTION, -r LOWESTRESOLUTION\n"
    "                        The lowest resolution of all loop files, i.e. 5kb,\n"
    "                        10kb and 25kb, please use 25000.\n"
    "\n"
    "Optional arguments:\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::vector<std::string> input_files;
    std::string out_file_name;
    std::int64_t lowest_resolution = 0;
};

// hicMergeLoops.py parse_arguments.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicMergeLoops",
                       "This script merges the locations of loops detected at several resolutions.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--inputFiles", "-i"})
        .nargs("+")
        .required()
        .input({"bedgraph", "txt"})
        .help("The loop files from hicDetectLoops.");
    required.add({"--outFileName", "-o"})
        .required()
        .output({"bedgraph"})
        .help("The name of the merged loop file.");
    required.add({"--lowestResolution", "-r"})
        .type("int")
        .required()
        .help("The lowest resolution of all loop files, i.e. 5kb, 10kb and 25kb, please use 25000.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.input_files = ns.strs("inputFiles");
    args.out_file_name = ns.str("outFileName");
    args.lowest_resolution = ns.integer("lowestResolution");
    return args;
}

// CPython's hash of an int that fits in 64 bits: the value modulo the
// Mersenne prime 2^61 - 1, with the sign kept, and -1 replaced by -2.
[[nodiscard]] std::int64_t python_hash_of_int(std::int64_t value) {
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
[[nodiscard]] std::int64_t python_hash_of_pair(std::int64_t first, std::int64_t second) {
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

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    try {
        hicx::TextTable table;
        bool first = true;
        for (const std::string& path : args.input_files) {
            hicx::TextTable part = hicx::TextTable::read_tsv(path);
            if (first) {
                table = std::move(part);
                first = false;
            } else {
                table.append(part);
            }
        }
        if (table.cols() < 6) {
            std::fprintf(stderr,
                         "hicMergeLoops: the loop files need six columns "
                         "'chr start end chr start end', found %zu\n",
                         table.cols());
            return 1;
        }

        // BedTool.from_dataframe(dataframe).sort().to_dataframe(...)
        {
            std::vector<std::string> chrom(table.rows());
            std::vector<std::int64_t> start(table.rows());
            for (std::size_t r = 0; r < table.rows(); ++r) {
                chrom[r] = table.column(0).text(r, ".");
                start[r] = table.column(1).as_int(r);
            }
            table = table.select_rows(hicx::bedtools::sort_order(chrom, start)).reparse();
        }
        // drop_duplicates(keep=False): every copy of a repeated row goes.
        table.drop_duplicates_drop_all();

        // intervalListToIntervalTree over the x anchors and over the y anchors.
        // The dictionary is keyed by chromosome and reset whenever the
        // chromosome of the input changes, so a chromosome that reappears after
        // another one loses everything it had; the input is sorted, so that
        // does not happen here, and the behaviour is reproduced anyway.
        std::vector<std::string> tree_order_x;
        std::vector<std::string> tree_order_y;
        std::map<std::string, Tree> tree_x;
        std::map<std::string, Tree> tree_y;
        {
            std::string previous;
            for (std::size_t r = 0; r < table.rows(); ++r) {
                const std::string chrom = table.column(0).text(r, ".");
                if (r == 0 || chrom != previous) {
                    if (tree_x.find(chrom) == tree_x.end()) {
                        tree_order_x.push_back(chrom);
                    }
                    tree_x[chrom] = Tree();
                    previous = chrom;
                }
                tree_x[chrom].add(table.column(1).as_int(r), table.column(2).as_int(r), r);
            }
            previous.clear();
            for (std::size_t r = 0; r < table.rows(); ++r) {
                const std::string chrom = table.column(3).text(r, ".");
                if (r == 0 || chrom != previous) {
                    if (tree_y.find(chrom) == tree_y.end()) {
                        tree_order_y.push_back(chrom);
                    }
                    tree_y[chrom] = Tree();
                    previous = chrom;
                }
                tree_y[chrom].add(table.column(4).as_int(r), table.column(5).as_int(r), r);
            }
        }

        std::vector<std::int32_t> x_interval;
        std::vector<std::int32_t> y_interval;
        bool x_valid = false;
        bool y_valid = false;
        std::string x_chrom_of_hits;
        std::string y_chrom_of_hits;

        for (std::size_t r = 0; r < table.rows(); ++r) {
            const std::string chrom_x = table.column(0).text(r, ".");
            const std::string chrom_y = table.column(3).text(r, ".");
            const std::int64_t x_start = table.column(1).as_int(r);
            const std::int64_t x_end = table.column(2).as_int(r);
            const std::int64_t y_start = table.column(4).as_int(r);
            const std::int64_t y_end = table.column(5).as_int(r);

            const std::int64_t factor_x =
                args.lowest_resolution - std::abs(x_end - x_start);
            const std::int64_t factor_y =
                args.lowest_resolution - std::abs(y_end - y_start);

            if (tree_x.find(chrom_x) != tree_x.end()) {
                x_interval = tree_x[chrom_x].overlap(x_start - factor_x - 1,
                                                     x_end + factor_x + 1);
                x_chrom_of_hits = chrom_x;
                x_valid = true;
            }
            if (tree_y.find(chrom_y) != tree_y.end()) {
                // hicMergeLoops.py:81 queries the y tree of the *x* chromosome.
                const auto found = tree_y.find(chrom_x);
                if (found == tree_y.end()) {
                    std::fprintf(stderr,
                                 "hicMergeLoops: loop %zu has y chromosome '%s' while "
                                 "hicMergeLoops.py:81 looks the y tree up under the x "
                                 "chromosome '%s', which does not exist; the reference "
                                 "raises KeyError here\n",
                                 r, chrom_y.c_str(), chrom_x.c_str());
                    return 1;
                }
                y_interval = found->second.overlap(y_start - factor_y - 1,
                                                   y_end + factor_y + 1);
                y_chrom_of_hits = chrom_x;
                y_valid = true;
            }
            if (!x_valid || !y_valid) {
                std::fprintf(stderr,
                             "hicMergeLoops: the first loop is on a chromosome that is "
                             "not in the interval trees; hicMergeLoops.py:84 then calls "
                             "len() on None and raises TypeError\n");
                return 1;
            }

            if (x_interval.size() <= 1 || y_interval.size() <= 1) {
                continue;
            }

            const Tree& tx = tree_x[x_chrom_of_hits];
            const Tree& ty = tree_y[y_chrom_of_hits];

            std::unordered_set<std::size_t> interest_x;
            interest_x.reserve(x_interval.size() * 2);
            for (const std::int32_t slot : x_interval) {
                interest_x.insert(tx.id[static_cast<std::size_t>(slot)]);
            }

            std::size_t max_index = 0;
            std::int64_t max_distance = 0;
            std::unordered_set<std::size_t> all_ids;
            for (const std::int32_t key : y_interval) {
                const auto slot = static_cast<std::size_t>(key);
                if (interest_x.find(ty.id[slot]) == interest_x.end()) {
                    continue;
                }
                const std::int64_t width = std::abs(ty.begin[slot] - ty.end[slot]);
                if (width > max_distance) {
                    max_distance = width;
                    max_index = ty.id[slot];
                }
                all_ids.insert(ty.id[slot]);
            }

            // The losers leave both trees in the iteration order of the
            // overlap sets, which shapes the trees for every later query.
            for (const std::int32_t key : x_interval) {
                const std::size_t identifier = tx.id[static_cast<std::size_t>(key)];
                if (identifier == max_index || all_ids.find(identifier) == all_ids.end()) {
                    continue;
                }
                if (!tree_x[x_chrom_of_hits].remove(key)) {
                    std::fprintf(stderr, "hicMergeLoops: interval removed twice; the reference "
                                         "raises ValueError here\n");
                    return 1;
                }
            }
            for (const std::int32_t key : y_interval) {
                const std::size_t identifier = ty.id[static_cast<std::size_t>(key)];
                if (identifier == max_index || all_ids.find(identifier) == all_ids.end()) {
                    continue;
                }
                if (!tree_y[y_chrom_of_hits].remove(key)) {
                    std::fprintf(stderr, "hicMergeLoops: interval removed twice; the reference "
                                         "raises ValueError here\n");
                    return 1;
                }
            }
        }

        // zip(tree_x, tree_y) pairs the chromosome keys positionally, in
        // dictionary insertion order, and takes the ids that survive in both.
        std::vector<std::size_t> result;
        const std::size_t pairs = std::min(tree_order_x.size(), tree_order_y.size());
        for (std::size_t k = 0; k < pairs; ++k) {
            const Tree& tx = tree_x[tree_order_x[k]];
            const Tree& ty = tree_y[tree_order_y[k]];
            std::unordered_set<std::size_t> alive_y;
            for (std::size_t slot = 0; slot < ty.id.size(); ++slot) {
                if (ty.present[slot]) {
                    alive_y.insert(ty.id[slot]);
                }
            }
            for (std::size_t slot = 0; slot < tx.id.size(); ++slot) {
                if (!tx.present[slot]) {
                    continue;
                }
                if (alive_y.find(tx.id[slot]) != alive_y.end()) {
                    result.push_back(tx.id[slot]);
                }
            }
        }
        std::sort(result.begin(), result.end());

        table.select_rows(result).write_tsv(args.out_file_name);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicMergeLoops: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicMergeLoops");
    return 0;
}
