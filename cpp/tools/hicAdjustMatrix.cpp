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
// Five behaviours are reproduced rather than fixed, all pinned by
// hicexplorer/test/general/test_hicAdjustMatrix.py:
//
//  * --maskBadRegions never opens the file it is given. hicAdjustMatrix.py:
//    161-165 loads the matrix and does nothing else, because :162 evaluates
//    len(pArgs.chromosomes) while --chromosomes is mutually exclusive with the
//    option and therefore always None. For an h5 input check_cooler is false
//    and short circuits the expression, so the matrix is written back
//    unchanged and the tool exits 0; for a cool input the len(None) raises a
//    TypeError. Both are reproduced, the second as an error message and exit
//    status 1.
//  * with none of --chromosomes, --regions and --maskBadRegions the tool logs
//    one line, writes nothing and exits 0. A caller cannot tell that from
//    success.
//  * a chromosome name that is not in the matrix is warned about and ignored,
//    and only an empty result exits 1. That is what makes the stray BED path
//    of test_remove_inter harmless: --chromosomes takes nargs='+', so the path
//    is swallowed as a fourth chromosome name.
//  * --regions --action remove clears orig_bin_ids, orig_cut_intervals and
//    nan_bins after masking (:157-159), so restoreMaskedBins has nothing to
//    put back and the bins, including the NaN bins that were folded into the
//    mask, really disappear.
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

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicAdjustMatrix: error: %s\n", message.c_str());
    std::exit(2);
}

Arguments parse_arguments(int argc, char** argv) {
    Arguments args;
    bool matrix_seen = false;
    bool out_seen = false;
    std::string* pending = nullptr;
    bool collecting_chromosomes = false;
    std::string exclusive_seen;

    const auto claim_exclusive = [&](const std::string& option) {
        if (!exclusive_seen.empty() && exclusive_seen != option) {
            fail("argument " + option + ": not allowed with argument " + exclusive_seen);
        }
        exclusive_seen = option;
    };

    for (int i = 1; i < argc; ++i) {
        const std::string token(argv[i]);
        if (pending != nullptr) {
            *pending = token;
            pending = nullptr;
            continue;
        }
        const bool is_option = token.size() > 1 && token[0] == '-' &&
                               std::isdigit(static_cast<unsigned char>(token[1])) == 0;
        if (!is_option) {
            if (collecting_chromosomes) {
                args.chromosomes.push_back(token);
                continue;
            }
            fail("unrecognized arguments: " + token);
        }
        collecting_chromosomes = false;

        std::string name = token;
        std::optional<std::string> inline_value;
        const std::size_t equals = token.find('=');
        if (equals != std::string::npos && token.rfind("--", 0) == 0) {
            name = token.substr(0, equals);
            inline_value = token.substr(equals + 1);
        }
        if (name == "-h" || name == "--help") {
            std::fputs(kUsage, stdout);
            std::fputs(kHelp, stdout);
            std::exit(0);
        }
        if (name == "--version") {
            std::printf("hicAdjustMatrix %s\n", hicx::kVersion);
            std::exit(0);
        }
        if (name == "-c" || name == "--chromosomes") {
            claim_exclusive("--chromosomes/-c");
            args.has_chromosomes = true;
            if (inline_value.has_value()) {
                args.chromosomes.push_back(*inline_value);
            } else {
                collecting_chromosomes = true;
            }
            continue;
        }
        std::string value_text;
        std::string* target = nullptr;
        if (name == "-m" || name == "--matrix") {
            target = &args.matrix;
            matrix_seen = true;
        } else if (name == "-o" || name == "--outFileName") {
            target = &args.out_file_name;
            out_seen = true;
        } else if (name == "-r" || name == "--regions") {
            claim_exclusive("--regions/-r");
            target = &value_text;
        } else if (name == "-mbr" || name == "--maskBadRegions") {
            claim_exclusive("--maskBadRegions/-mbr");
            target = &value_text;
        } else if (name == "-a" || name == "--action") {
            target = &value_text;
        } else if (name == "-iih" || name == "--interIntraHandling") {
            target = &value_text;
        } else {
            fail("unrecognized arguments: " + token);
        }

        if (inline_value.has_value()) {
            *target = *inline_value;
        } else {
            if (i + 1 >= argc) {
                fail("argument " + name + ": expected one argument");
            }
            *target = std::string(argv[++i]);
        }
        if (name == "-r" || name == "--regions") {
            args.regions = value_text;
        } else if (name == "-mbr" || name == "--maskBadRegions") {
            args.mask_bad_regions = value_text;
        } else if (name == "-a" || name == "--action") {
            if (value_text == "keep") {
                args.action = Action::Keep;
            } else if (value_text == "remove") {
                args.action = Action::Remove;
            } else if (value_text == "mask") {
                args.action = Action::Mask;
            } else {
                fail("argument --action/-a: invalid choice: '" + value_text +
                     "' (choose from 'keep', 'remove', 'mask')");
            }
        } else if (name == "-iih" || name == "--interIntraHandling") {
            if (value_text == "inter") {
                args.inter_intra = hicx::InterIntra::Inter;
            } else if (value_text == "intra") {
                args.inter_intra = hicx::InterIntra::Intra;
            } else {
                fail("argument --interIntraHandling/-iih: invalid choice: '" +
                     value_text + "' (choose from None, 'inter', 'intra')");
            }
        }
    }
    if (pending != nullptr) {
        fail("expected one argument");
    }
    std::string missing;
    if (!matrix_seen) {
        missing += "--matrix/-m";
    }
    if (!out_seen) {
        missing += missing.empty() ? "--outFileName/-o" : ", --outFileName/-o";
    }
    if (!missing.empty()) {
        fail("the following arguments are required: " + missing);
    }
    if (args.has_chromosomes && args.chromosomes.empty()) {
        fail("argument --chromosomes/-c: expected at least one argument");
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
            for (const std::string& chromosome : args.chromosomes) {
                if (std::find(chromosomes_list.begin(), chromosomes_list.end(),
                              chromosome) != chromosomes_list.end()) {
                    to_operate_on.push_back(chromosome);
                } else {
                    std::fprintf(stderr,
                                 "WARNING:hicexplorer.hicAdjustMatrix:Chromosome not "
                                 "available in matrix: %s %s\n",
                                 args.matrix.c_str(), chromosome.c_str());
                }
            }
            if (to_operate_on.empty()) {
                std::fprintf(stderr,
                             "ERROR:hicexplorer.hicAdjustMatrix:No valid chromosome "
                             "given: %s. Available: %s\n",
                             chromosome_list_repr(args.chromosomes).c_str(),
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
        } else if (args.regions.has_value()) {
            hic = hicx::ToolMatrix::load(args.matrix);
            const std::vector<std::string> chromosomes_list = boundary_names(*hic);
            const hicx::BinTable bins(hic->data().cut_intervals);

            std::ifstream bed(*args.regions);
            if (!bed) {
                std::fprintf(stderr,
                             "hicAdjustMatrix: error: argument --regions/-r: can't "
                             "open '%s'\n",
                             args.regions->c_str());
                return 2;
            }
            struct Region {
                std::string chrom;
                std::int64_t start = 0;
                std::int64_t end = 0;
            };
            std::vector<Region> genomic_regions;
            // hicAdjustMatrix.py:124 leaks the loop variable `chrom`, and the
            // --action remove branch at :144 then looks up whichever
            // chromosome the *last* readable BED line named, not the one of
            // the region it is reporting on. It only feeds a warning, so it is
            // reproduced here as one.
            std::string last_chrom;
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
                last_chrom = fields[0];
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
                if (args.action == Action::Remove) {
                    const std::optional<hicx::BinRange> chrom_range =
                        bins.chrom_bin_range(last_chrom);
                    if (chrom_range.has_value() && range->first > chrom_range->first &&
                        range->second < chrom_range->last - 1) {
                        std::fprintf(stderr,
                                     "WARNING:hicexplorer.hicAdjustMatrix:%s:%lld-%lld "
                                     "entry may generate discounted regions on a "
                                     "chromosome.Please consider using `mask` action "
                                     "to deal with that.\n",
                                     last_chrom.c_str(),
                                     static_cast<long long>(range->first),
                                     static_cast<long long>(range->second));
                    }
                }
            }

            if (args.action == Action::Keep) {
                hicx::reorder_bins(hic->data(), indices);
            } else if (args.action == Action::Mask) {
                hicx::mask_and_restore_bins(hic->data(), indices);
            } else {
                // maskBins, then orig_bin_ids, orig_cut_intervals and nan_bins
                // are cleared, so the bins are gone for good.
                hicx::delete_bins(hic->data(), indices);
            }
            hic->refresh_boundaries();
        } else if (args.mask_bad_regions.has_value()) {
            if (hicx::check_cooler(args.matrix)) {
                // hicAdjustMatrix.py:162 evaluates len(pArgs.chromosomes) with
                // chromosomes always None, which check_cooler short circuits
                // only for h5. Same exit status as the Python traceback.
                std::fputs("hicAdjustMatrix: TypeError: object of type 'NoneType' has "
                           "no len(). hicAdjustMatrix.py:162 takes len() of "
                           "--chromosomes, which is mutually exclusive with "
                           "--maskBadRegions and therefore always None\n",
                           stderr);
                return 1;
            }
            // The BED file is never opened; the matrix is written back as it
            // was loaded.
            hic = hicx::ToolMatrix::load(args.matrix);
        } else {
            std::fputs("INFO:hicexplorer.hicAdjustMatrix:No data to adjust given. "
                       "Please specify either --chromosomes or --region parameter.\n",
                       stderr);
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
