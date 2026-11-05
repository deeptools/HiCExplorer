// Port of hicexplorer/hicMergeTADbins.py.
//
// Merges the bins of a matrix per TAD, using a BED file of domains to decide
// where one merged bin ends and the next begins. It is the second caller of
// reduceMatrix.reduce_matrix with pDiagonal=True, so it shares
// hicx::reduce_matrix with hicMergeMatrixBins rather than reimplementing it.
//
// Three behaviours are reproduced rather than fixed, all pinned by
// hicexplorer/test/general/test_hicMergeTADbins.py:
//
//  * hicMergeTADbins.py:87 sets correction_factors to None before saving, by
//    design: the factors belong to the unmerged bins and would no longer match
//    the size of the matrix. Note that for an h5 input the factors are not in
//    correction_factors to begin with, because the hicmatrix loader unpacks
//    correction_factors and distance_counts the wrong way round
//    (cpp/PLAN.md 2.7 quirk 1); both end up absent from the output.
//  * the merged matrix does not conserve the total count. reduce_matrix sums
//    the upper triangle block by block and then rebuilds the symmetric matrix
//    as R + R.T - diag(R). With pDiagonal=True the subtracted diagonal is the
//    whole within-TAD block sum rather than the original main diagonal, so on
//    Li_et_al_2015.h5 the symmetric total falls from 30,482,969.637651745 to
//    23,695,524.863065504, a loss of 22 percent. The upper triangle *is*
//    conserved. See test_..._conserves_the_upper_triangle_but_not_the_full_sum.
//  * hicMergeTADbins.py:129 unpacks getRegionBinRange without checking for
//    None, so a domain that lies outside the matrix ends the run with a bare
//    TypeError and no diagnostic. This port prints a message naming the domain
//    and exits 1, which is the status the Python traceback produces, because a
//    silent crash carries no information a caller can act on.

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/adjust_ops.hpp"
#include "hicx/bins.hpp"
#include "hicx/reduce_matrix.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicMergeTADbins [-h] --matrix MATRIX --domains DOMAINS --outFile\n"
    "                       OUTFILE [--version]\n";

const char* const kHelp =
    "\n"
    "Uses a BED file of domains or TAD boundaries to merge\n"
    "the bin counts of a Hi-C matrix per TAD.\n"
    "\n"
    "The output matrix contains the total counts per TAD and\n"
    "the total contacts with all other TADs.\n"
    "\n"
    "options:\n"
    "  -h, --help            show this help message and exit\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        Path to Hi-C matrix to use.\n"
    "  --domains DOMAINS     Path to a bed file containing the domains.\n"
    "  --outFile OUTFILE, -o OUTFILE\n"
    "                        Name for the resulting matrix file.\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::string matrix;
    std::string domains;
    std::string out_file;
};

// hicMergeTADbins.py parse_arguments: no argument groups, argparse's own
// -h/--help first.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicMergeTADbins",
                       "Uses a BED file of domains or TAD boundaries to merge the bin counts of a "
                       "Hi-C matrix per TAD.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& options = parser.group("options");
    options.add({"-h", "--help"}).action(cli::Action::Help).help("show this help message and exit");
    options.add({"--matrix", "-m"})
        .required()
        .input({"h5", "cool", "mcool"})
        .help("Path to Hi-C matrix to use.");
    options.add({"--domains"})
        .file_type("r")
        .required()
        .input({"bed"})
        .help("Path to a bed file containing the domains.");
    options.add({"--outFile", "-o"})
        .required()
        .output({"h5", "cool"})
        .help("Name for the resulting matrix file.");
    options.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.matrix = ns.str("matrix");
    args.domains = ns.str("domains");
    args.out_file = ns.str("outFile");
    return args;
}

std::string trim_right(const std::string& text) {
    std::size_t end = text.size();
    while (end > 0 && std::isspace(static_cast<unsigned char>(text[end - 1])) != 0) {
        --end;
    }
    return text.substr(0, end);
}

std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (true) {
        const std::size_t marker = line.find('\t', start);
        if (marker == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, marker - start));
        start = marker + 1;
    }
    return fields;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    try {
        hicx::ToolMatrix hic = hicx::ToolMatrix::load(args.matrix);
        // hicMergeTADbins.py:140 and :57 both call restoreMaskedBins. Nothing
        // is masked straight after a load, so both are no operations here.

        const hicx::BinTable bins(hic.data().cut_intervals);

        std::ifstream bed(args.domains);
        if (!bed) {
            std::fprintf(stderr,
                         "hicMergeTADbins: error: argument --domains: can't open '%s'\n",
                         args.domains.c_str());
            return 2;
        }
        // get_boundary_bin_id: the first and the last bin of every domain, as
        // a sorted set.
        std::set<std::int64_t> boundary_set;
        std::string line;
        int line_number = 0;
        while (std::getline(bed, line)) {
            ++line_number;
            if (line.rfind("browser", 0) == 0 || line.rfind("track", 0) == 0 ||
                line.rfind("#", 0) == 0) {
                continue;
            }
            const std::vector<std::string> fields = split_tabs(trim_right(line));
            if (fields.size() < 3) {
                std::fprintf(stderr, "hicMergeTADbins: could not read line\n%s\n",
                             line.c_str());
                return 1;
            }
            std::int64_t start = 0;
            std::int64_t end = 0;
            try {
                start = std::stoll(fields[1]);
                end = std::stoll(fields[2]);
            } catch (const std::exception&) {
                std::fprintf(stderr,
                             "hicMergeTADbins: error reading line: %d. One of the "
                             "fields is not an integer.\n",
                             line_number);
                return 1;
            }
            if (start > end) {
                std::fprintf(stderr,
                             "hicMergeTADbins: error in line #%d, end1 larger than "
                             "start1 in %s\n",
                             line_number, line.c_str());
                return 1;
            }
            const std::optional<std::pair<std::int64_t, std::int64_t>> range =
                bins.region_bin_range(fields[0], start, end);
            if (!range.has_value()) {
                // hicMergeTADbins.py:129 unpacks None here and dies with a
                // TypeError. Same exit status, with the domain named.
                std::fprintf(stderr,
                             "hicMergeTADbins: domain %s:%lld-%lld of line %d is not "
                             "inside the matrix; the Python reference raises a bare "
                             "TypeError here (hicMergeTADbins.py:129)\n",
                             fields[0].c_str(), static_cast<long long>(start),
                             static_cast<long long>(end), line_number);
                return 1;
            }
            boundary_set.insert(range->first);
            boundary_set.insert(range->second);
        }
        const std::vector<std::int64_t> boundaries(boundary_set.begin(),
                                                   boundary_set.end());

        hicx::BinMergePlan plan =
            hicx::plan_tad_merge(hic.data().cut_intervals, boundaries);
        if (plan.bins_to_merge.empty()) {
            std::fputs("INFO:hicexplorer.hicMergeTADbins:Nothing to merge.\n", stderr);
            hicx::report_resource_usage("hicMergeTADbins");
            return 0;
        }

        // hicMergeTADbins.py:87: the factors no longer match the merged bins.
        hic.data().correction_factors.reset();
        hic.data().correction_factors_are_column = false;

        // update_matrix(reduce_matrix(matrix, bins_to_merge, diagonal=True),
        //               new_bins)
        hic.matrix() = hicx::reduce_matrix(hic.matrix(), plan.bins_to_merge, true, true);
        hic.data().cut_intervals = std::move(plan.intervals);
        hic.refresh_boundaries();
        hic.data().nan_bins = hicx::empty_column_bins(hic.matrix());

        hic.save(args.out_file);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicMergeTADbins: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicMergeTADbins");
    return 0;
}
