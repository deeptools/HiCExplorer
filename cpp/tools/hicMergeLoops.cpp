// Port of hicexplorer/hicMergeLoops.py.
//
// Merges loop calls made at several resolutions: a loop is dropped when it
// overlaps, in both its x and its y anchor, a loop of a lower resolution, and
// the lower-resolution call is the one kept.
//
// The reference builds two intervaltree dictionaries, one over the x anchors
// and one over the y anchors, and walks the loops removing the losers from
// both trees. Three defects of that code are fixed in HiCExplorer 4 (the tests
// in hicexplorer/test/general/test_hicMergeLoops.py and the API cases pin the
// corrected behaviour), and one property is kept:
//
//  1. **The y tree was indexed by the x chromosome** (hicMergeLoops.py:81:
//     `target_regions_intervaltree_y[loop[0]]`, where loop[3] is the y
//     chromosome). On an inter-chromosomal loop the lookup found the wrong
//     chromosome's tree or raised KeyError. HiCExplorer 4 queries the y tree
//     of the y chromosome and removes the losers from it. The final zip over
//     the two dictionaries paired the x and y trees by position, which is
//     wrong as soon as the chromosomes of the x and y anchors differ; a loop
//     now survives when it is still present in an x tree and in a y tree.
//  2. **x_interval and y_interval survived across iterations.** They were only
//     assigned inside `if loop[0] in ...` (":77-82"), so a loop on a
//     chromosome that was not in the tree was evaluated against the previous
//     loop's intervals, and a first such loop raised TypeError on len(None).
//     HiCExplorer 4 evaluates every loop against its own intervals and leaves
//     a loop whose chromosome has no tree alone.
//  3. **The winner among equally wide overlapping anchors is decided by
//     CPython's set iteration order**, which the port reproduces. This is
//     kept: `IntervalTree.overlap` returns a set, the scan at ":98-102" keeps
//     the first interval whose width is *strictly* greater than the best so
//     far, so equally wide intervals (every overlap when all inputs share one
//     resolution) break the tie by the order that set iterates in. That order
//     is deterministic (Interval hashes as `hash((begin, end))`, a tuple of
//     ints, which PYTHONHASHSEED does not salt), the help text does not say
//     which of equally wide loops is kept, and choosing another rule would be
//     a design decision that moves rows of the merge. `Tree`
//     (hicx/py_intervaltree.hpp) is our centered AVL interval tree with the
//     shape and the member sets of intervaltree 3.1.0, and `OrderedSet`
//     (hicx/py_set_order.hpp) is the slot order of a CPython set. Both are
//     checked against the real packages on random sequences
//     (tests/test_cpython_oracles.cpp).
//  4. **drop_duplicates(keep=False)** removed every copy of a repeated row
//     (":165"), so a loop called identically at two resolutions, or in the
//     same file given twice, disappeared from the merge entirely. HiCExplorer 4
//     keeps the first copy.
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
#include "hicx/py_intervaltree.hpp"
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

}  // namespace

// Tree, the reproduction of intervaltree's query order that property 3 needs,
// lives in hicx/py_intervaltree.hpp, which the public Python API shares.
using hicx::pytree::Tree;

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
        // drop_duplicates(): the first copy of a repeated row stays (3.7.6 used
        // keep=False and dropped every copy).
        table.drop_duplicates();

        // intervalListToIntervalTree over the x anchors and over the y anchors,
        // one tree per chromosome. hiCMatrix.intervalListToIntervalTree reset
        // the dictionary entry whenever the chromosome of the list changed, so
        // a chromosome that reappeared after another one lost everything it
        // had. The x anchors are sorted, but the y chromosomes of
        // inter-chromosomal loops are not, so HiCExplorer 4 keeps adding to the
        // tree of the chromosome.
        std::vector<std::string> tree_order_x;
        std::vector<std::string> tree_order_y;
        std::map<std::string, Tree> tree_x;
        std::map<std::string, Tree> tree_y;
        for (std::size_t r = 0; r < table.rows(); ++r) {
            const std::string chrom = table.column(0).text(r, ".");
            if (tree_x.find(chrom) == tree_x.end()) {
                tree_order_x.push_back(chrom);
            }
            tree_x[chrom].add(table.column(1).as_int(r), table.column(2).as_int(r), r);
        }
        for (std::size_t r = 0; r < table.rows(); ++r) {
            const std::string chrom = table.column(3).text(r, ".");
            if (tree_y.find(chrom) == tree_y.end()) {
                tree_order_y.push_back(chrom);
            }
            tree_y[chrom].add(table.column(4).as_int(r), table.column(5).as_int(r), r);
        }

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

            // A loop whose chromosome has no tree has nothing to merge with
            // (3.7.6 reused the previous loop's intervals here).
            const auto found_x = tree_x.find(chrom_x);
            const auto found_y = tree_y.find(chrom_y);
            if (found_x == tree_x.end() || found_y == tree_y.end()) {
                continue;
            }
            const std::vector<std::int32_t> x_interval =
                found_x->second.overlap(x_start - factor_x - 1, x_end + factor_x + 1);
            const std::vector<std::int32_t> y_interval =
                found_y->second.overlap(y_start - factor_y - 1, y_end + factor_y + 1);

            if (x_interval.size() <= 1 || y_interval.size() <= 1) {
                continue;
            }

            Tree& tx = found_x->second;
            Tree& ty = found_y->second;

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
                tx.remove(key);
            }
            for (const std::int32_t key : y_interval) {
                const std::size_t identifier = ty.id[static_cast<std::size_t>(key)];
                if (identifier == max_index || all_ids.find(identifier) == all_ids.end()) {
                    continue;
                }
                ty.remove(key);
            }
        }

        // A loop survives when it is present in an x tree and in a y tree
        // (3.7.6 paired the trees of the two dictionaries by position).
        std::vector<std::size_t> result;
        std::unordered_set<std::size_t> alive_y;
        for (const std::string& name : tree_order_y) {
            const Tree& ty = tree_y[name];
            for (std::size_t slot = 0; slot < ty.id.size(); ++slot) {
                if (ty.present[slot]) {
                    alive_y.insert(ty.id[slot]);
                }
            }
        }
        for (const std::string& name : tree_order_x) {
            const Tree& tx = tree_x[name];
            for (std::size_t slot = 0; slot < tx.id.size(); ++slot) {
                if (tx.present[slot] && alive_y.find(tx.id[slot]) != alive_y.end()) {
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
