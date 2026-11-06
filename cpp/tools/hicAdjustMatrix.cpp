// Port of hicexplorer/hicAdjustMatrix.py.
//
// Keeps, removes or masks a list of chromosomes or of BED regions, and can
// additionally blank the inter- or intra-chromosomal contacts.
//
// The tool is a thin shell over four hiCMatrix mutators, all of which live in
// the core: reorderChromosomes and reorderBins (hicx::reorder_bins), maskBins
// followed by the restoreMaskedBins that save performs
// (hicx::mask_and_restore_bins), and maskBins with the restore deliberately
// disabled (hicx::delete_bins).
//
// Five behaviours of the 3.7.6 Python are fixed in HiCExplorer 4:
//
//  * --maskBadRegions never opened the file it was given. hicAdjustMatrix.py:
//    161-165 loaded the matrix and did nothing else (an h5 matrix was written
//    back unchanged, a cool input died with a TypeError from
//    len(pArgs.chromosomes), which is always None). HiCExplorer 4 masks the
//    bins of the regions in the BED file, as --regions does with --action mask.
//  * with none of --chromosomes, --regions and --maskBadRegions the Python
//    logged one line, wrote nothing and exited 0. HiCExplorer 4 exits 1 with
//    an error that names the three options.
//  * a chromosome name that is not in the matrix was only warned about, so a
//    misspelt name silently changed the selection (with --action remove the
//    chromosome stayed in the output). HiCExplorer 4 exits 1 and names every
//    chromosome that is not in the matrix. Chromosomes in a --regions or
//    --maskBadRegions BED file that are not in the matrix are still only
//    warned about, because a BED file may cover more than one assembly part.
//  * --regions --action remove cleared orig_bin_ids, orig_cut_intervals and
//    nan_bins after masking (:157-159), so the bins of the regions went, and
//    so did every empty (NaN) bin of the matrix, which nobody had asked to
//    remove. HiCExplorer 4 removes exactly the bins of the regions; the empty
//    bins that remain stay in the matrix, as they do with --chromosomes
//    --action remove.
//  * the "entry may generate discounted regions" warning of --action remove
//    named the chromosome of the last BED line (a leaked loop variable) and
//    looked up its range, instead of those of the region it reports on.
//    HiCExplorer 4 uses the chromosome of the region.
//
// Reproduced because it is the documented output and not a defect:
//  * --interIntraHandling inter zeroes only the blocks above the diagonal. The
//    mirror entries stay in memory and are dropped by the triu(k=0) both
//    writers apply, so the file is the one the help text promises.
//
// This tool takes the both-triangles exemption of cpp/PLAN.md 4.4 rule 2,
// because reorderBins permutes rows and columns independently. It takes it
// only when it has to: hicx::select_bins keeps the halved storage whenever the
// index list is strictly increasing, which covers --action remove, the
// single-chromosome cooler fast path and every --regions BED whose intervals
// are in file order.

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "hicx/adjust_ops.hpp"
#include "hicx/bins.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/matrix_ops.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/argparse.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicAdjustMatrix --matrix MATRIX --outFileName OUTFILENAME\n"
    "                       [--chromosomes CHROMOSOMES [CHROMOSOMES ...] |\n"
    "                       --regions REGIONS | --maskBadRegions MASKBADREGIONS]\n"
    "                       [--action {keep,remove,mask}]\n"
    "                       [--interIntraHandling {inter,intra}] [--help]\n"
    "                       [--version]\n";

const char* const kHelp =
    "\n"
    "                    This tool adjusts hic matrices by keeping, removing or\n"
    "                    masking a given list of regions or chromosmes.\n"
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        The Hi-C matrix to adjust. HiCExplorer supports the\n"
    "                        following file formats: h5 (native HiCExplorer format)\n"
    "                        and cool.\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        File name to save the adjusted matrix.\n"
    "\n"
    "Optional arguments:\n"
    "  --action {keep,remove,mask}, -a {keep,remove,mask}\n"
    "                        Keep, remove or mask the list of specified\n"
    "                        chromosomes/regions. keep/remove: These options\n"
    "                        keep/remove bins of matrix by deleting them. This may\n"
    "                        cause issue plotting the matrix if several parts of a\n"
    "                        single chromosome are going to be deleted. In that\n"
    "                        case, one may consider using the mask option (Default:\n"
    "                        keep).\n"
    "  --interIntraHandling {inter,intra}, -iih {inter,intra}\n"
    "                        Remove the inter- or intra-chromosomal contacts of the\n"
    "                        given chromosomes. (Default: None).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n"
    "\n"
    "  --chromosomes CHROMOSOMES [CHROMOSOMES ...], -c CHROMOSOMES [CHROMOSOMES ...]\n"
    "                        List of chromosomes to keep/remove.\n"
    "  --regions REGIONS, -r REGIONS\n"
    "                        BED file which stores a list of regions to keep/remove.\n"
    "  --maskBadRegions MASKBADREGIONS, -mbr MASKBADREGIONS\n"
    "                        Bad regions are identified and masked.\n";

enum class Action { Keep, Remove, Mask };

struct Arguments {
    std::string matrix;
    std::string out_file_name;
    std::vector<std::string> chromosomes;
    bool has_chromosomes = false;
    std::optional<std::string> regions;
    std::optional<std::string> mask_bad_regions;
    Action action = Action::Keep;
    std::optional<hicx::InterIntra> inter_intra;
};

// hicAdjustMatrix.py parse_arguments. The mutually exclusive group hangs off
// the parser itself, so argparse lists its options under the parser's own
// "options" group.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicAdjustMatrix",
                       "This tool adjusts hic matrices by keeping, removing or masking a given list "
                       "of regions or chromosmes.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrix", "-m"})
        .required()
        .input({"h5", "cool", "mcool"})
        .help("The Hi-C matrix to adjust.");
    required.add({"--outFileName", "-o"})
        .required()
        .output({"h5", "cool"})
        .help("File name to save the adjusted matrix.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    cli::ArgumentGroup& options = parser.group("options");
    cli::MutuallyExclusiveGroup& selection = parser.mutually_exclusive(options);
    selection.add({"--chromosomes", "-c"})
        .nargs("+")
        .help("List of chromosomes to keep/remove.");
    selection.add({"--regions", "-r"})
        .input({"bed"})
        .help("BED file which stores a list of regions to keep/remove.");
    selection.add({"--maskBadRegions", "-mbr"})
        .input({"bed"})
        .help("Bad regions are identified and masked.");
    optional.add({"--action", "-a"})
        .default_value("keep")
        .choices({"keep", "remove", "mask"})
        .help("Keep, remove or mask the list of specified chromosomes/regions.");
    optional.add({"--interIntraHandling", "-iih"})
        .choices({hicx::json::Value::null(), hicx::json::Value::string("inter"),
                  hicx::json::Value::string("intra")})
        .help("Remove the inter- or intra-chromosomal contacts of the given chromosomes.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.matrix = ns.str("matrix");
    args.out_file_name = ns.str("outFileName");
    args.has_chromosomes = ns.given("chromosomes");
    args.chromosomes = ns.strs("chromosomes");
    args.regions = ns.opt_str("regions");
    args.mask_bad_regions = ns.opt_str("maskBadRegions");
    const std::string action = ns.str("action");
    args.action = action == "remove" ? Action::Remove
                  : action == "mask" ? Action::Mask
                                     : Action::Keep;
    if (const std::optional<std::string> handling = ns.opt_str("interIntraHandling")) {
        args.inter_intra = *handling == "inter" ? hicx::InterIntra::Inter : hicx::InterIntra::Intra;
    }
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

std::string chromosome_list_repr(const std::vector<std::string>& names) {
    std::string text = "[";
    for (std::size_t i = 0; i < names.size(); ++i) {
        if (i != 0) {
            text += ", ";
        }
        text += "'" + names[i] + "'";
    }
    text += "]";
    return text;
}

std::vector<std::string> boundary_names(const hicx::ToolMatrix& hic) {
    std::vector<std::string> names;
    names.reserve(hic.boundaries().size());
    for (const auto& entry : hic.boundaries()) {
        names.push_back(entry.first);
    }
    return names;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    try {
        std::optional<hicx::ToolMatrix> hic;

        if (args.has_chromosomes) {
            // hicAdjustMatrix.py:73-79. One chromosome, --action keep and a
            // cool input load only that chromosome's block, which gives a
            // different NaN bin list from a whole file load followed by a
            // selection, so the branch is observable and is reproduced.
            if (hicx::check_cooler(args.matrix) && args.chromosomes.size() == 1 &&
                args.action == Action::Keep) {
                const hicx::CoolFile cool(args.matrix);
                const std::vector<std::string>& available = cool.chrom_names();
                if (std::find(available.begin(), available.end(), args.chromosomes[0]) ==
                    available.end()) {
                    std::fprintf(stderr,
                                 "ERROR:hicexplorer.hicAdjustMatrix:Chromosome not "
                                 "available in matrix: %s %s\n",
                                 args.matrix.c_str(), args.chromosomes[0].c_str());
                    return 1;
                }
                hic = hicx::ToolMatrix::load(args.matrix, args.chromosomes[0]);
            } else {
                hic = hicx::ToolMatrix::load(args.matrix);
            }

            std::vector<std::string> chromosomes_list = boundary_names(*hic);
            std::vector<std::string> to_operate_on;
            std::vector<std::string> unknown;
            for (const std::string& chromosome : args.chromosomes) {
                if (std::find(chromosomes_list.begin(), chromosomes_list.end(),
                              chromosome) != chromosomes_list.end()) {
                    to_operate_on.push_back(chromosome);
                } else {
                    unknown.push_back(chromosome);
                }
            }
            if (!unknown.empty()) {
                // Fixed in HiCExplorer 4: 3.7.6 only warned and went on with
                // the names it knew.
                std::fprintf(stderr,
                             "ERROR:hicexplorer.hicAdjustMatrix:Chromosome not "
                             "available in matrix: %s %s. Available: %s\n",
                             args.matrix.c_str(), chromosome_list_repr(unknown).c_str(),
                             chromosome_list_repr(chromosomes_list).c_str());
                return 1;
            }

            const auto bins_of = [&](const std::vector<std::string>& names) {
                std::vector<std::int64_t> order;
                for (const std::string& chromosome : names) {
                    for (const auto& entry : hic->boundaries()) {
                        if (entry.first != chromosome) {
                            continue;
                        }
                        for (std::int64_t bin = entry.second.first;
                             bin < entry.second.last; ++bin) {
                            order.push_back(bin);
                        }
                        break;
                    }
                }
                return order;
            };

            if (args.action == Action::Keep) {
                hicx::reorder_bins(hic->data(), bins_of(to_operate_on));
                hic->refresh_boundaries();
            } else if (args.action == Action::Remove) {
                // The Python removes the named chromosomes from the list of
                // all chromosomes and reorders to what is left, so the
                // surviving chromosomes keep their original order.
                std::vector<std::string> remaining;
                for (const std::string& chromosome : chromosomes_list) {
                    if (std::find(to_operate_on.begin(), to_operate_on.end(),
                                  chromosome) == to_operate_on.end()) {
                        remaining.push_back(chromosome);
                    }
                }
                hicx::reorder_bins(hic->data(), bins_of(remaining));
                hic->refresh_boundaries();
            } else {
                // maskChromosomes, then the restoreMaskedBins that save does.
                hicx::mask_and_restore_bins(hic->data(), bins_of(to_operate_on));
                hic->refresh_boundaries();
            }
        } else if (args.regions.has_value() || args.mask_bad_regions.has_value()) {
            // --maskBadRegions masks the bins of its BED file (fixed in
            // HiCExplorer 4, see the top of this file), whatever --action says.
            const bool bad_regions = args.mask_bad_regions.has_value();
            const std::string& regions_file = bad_regions ? *args.mask_bad_regions : *args.regions;
            const Action action = bad_regions ? Action::Mask : args.action;
            hic = hicx::ToolMatrix::load(args.matrix);
            const std::vector<std::string> chromosomes_list = boundary_names(*hic);
            const hicx::BinTable bins(hic->data().cut_intervals);

            std::ifstream bed(regions_file);
            if (!bed) {
                std::fprintf(stderr,
                             "hicAdjustMatrix: error: argument %s: can't "
                             "open '%s'\n",
                             bad_regions ? "--maskBadRegions/-mbr" : "--regions/-r",
                             regions_file.c_str());
                return 2;
            }
            struct Region {
                std::string chrom;
                std::int64_t start = 0;
                std::int64_t end = 0;
            };
            std::vector<Region> genomic_regions;
            std::string line;
            while (std::getline(bed, line)) {
                const std::string stripped = trim_right(line);
                // hicAdjustMatrix.py:120 tests len(line), the raw line, not
                // len(_line), the fields. readlines keeps the newline, so the
                // threshold is two characters of content.
                if (line.size() + 1 < 3) {
                    std::fputs("WARNING:hicexplorer.hicAdjustMatrix:An entry shorter "
                               "than 3 columns has been found!\n",
                               stderr);
                    continue;
                }
                const std::vector<std::string> fields = split_tabs(stripped);
                if (fields.size() < 3) {
                    continue;
                }
                Region region;
                region.chrom = fields[0];
                try {
                    region.start = std::stoll(fields[1]);
                    region.end = std::stoll(fields[2]);
                } catch (const std::exception&) {
                    std::fprintf(stderr,
                                 "hicAdjustMatrix: %s is not an integer in line: %s\n",
                                 fields[1].c_str(), stripped.c_str());
                    return 1;
                }
                if (std::find(chromosomes_list.begin(), chromosomes_list.end(),
                              region.chrom) != chromosomes_list.end()) {
                    genomic_regions.push_back(region);
                } else {
                    std::fprintf(stderr,
                                 "WARNING:hicexplorer.hicAdjustMatrix:Chromosome not "
                                 "available in matrix, ignoring regions: %s %s\n",
                                 args.matrix.c_str(), region.chrom.c_str());
                }
            }
            if (genomic_regions.empty()) {
                std::fprintf(stderr,
                             "ERROR:hicexplorer.hicAdjustMatrix:No valid chromosome "
                             "given. Available: %s\n",
                             chromosome_list_repr(chromosomes_list).c_str());
                return 1;
            }

            std::vector<std::int64_t> indices;
            for (const Region& region : genomic_regions) {
                const std::optional<std::pair<std::int64_t, std::int64_t>> range =
                    bins.region_bin_range(region.chrom, region.start, region.end);
                if (!range.has_value()) {
                    continue;  // getRegionBinRange returned None and is skipped
                }
                for (std::int64_t bin = range->first; bin <= range->second; ++bin) {
                    indices.push_back(bin);  // the end bin is inclusive
                }
                if (action == Action::Remove) {
                    // Fixed in HiCExplorer 4: 3.7.6 looked up the chromosome of
                    // the last BED line here (a leaked loop variable).
                    const std::optional<hicx::BinRange> chrom_range =
                        bins.chrom_bin_range(region.chrom);
                    if (chrom_range.has_value() && range->first > chrom_range->first &&
                        range->second < chrom_range->last - 1) {
                        std::fprintf(stderr,
                                     "WARNING:hicexplorer.hicAdjustMatrix:%s:%lld-%lld "
                                     "entry may generate discounted regions on a "
                                     "chromosome.Please consider using `mask` action "
                                     "to deal with that.\n",
                                     region.chrom.c_str(),
                                     static_cast<long long>(range->first),
                                     static_cast<long long>(range->second));
                    }
                }
            }

            if (action == Action::Keep) {
                hicx::reorder_bins(hic->data(), indices);
            } else if (action == Action::Mask) {
                hicx::mask_and_restore_bins(hic->data(), indices);
            } else {
                // Fixed in HiCExplorer 4: exactly the bins of the regions are
                // removed. 3.7.6 also dropped every NaN bin of the matrix.
                const std::int64_t rows = hic->data().matrix.rows();
                std::vector<char> drop(static_cast<std::size_t>(rows), 0);
                for (const std::int64_t bin : indices) {
                    drop[static_cast<std::size_t>(bin)] = 1;
                }
                std::vector<std::int64_t> survivors;
                for (std::int64_t bin = 0; bin < rows; ++bin) {
                    if (drop[static_cast<std::size_t>(bin)] == 0) {
                        survivors.push_back(bin);
                    }
                }
                hicx::reorder_bins(hic->data(), survivors);
            }
            hic->refresh_boundaries();
        } else {
            std::fputs("ERROR:hicexplorer.hicAdjustMatrix:No data to adjust given. "
                       "Please specify one of --chromosomes, --regions or "
                       "--maskBadRegions.\n",
                       stderr);
            return 1;
        }

        if (args.inter_intra.has_value()) {
            if (!hic.has_value()) {
                // hic_matrix is None and the Python dies on
                // None.chrBinBoundaries with an AttributeError.
                std::fputs("hicAdjustMatrix: AttributeError: 'NoneType' object has no "
                           "attribute 'chrBinBoundaries'\n",
                           stderr);
                return 1;
            }
            hicx::zero_inter_or_intra(hic->matrix(), hic->boundaries(),
                                      *args.inter_intra);
        }

        if (hic.has_value()) {
            hic->save(args.out_file_name);
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicAdjustMatrix: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicAdjustMatrix");
    return 0;
}
