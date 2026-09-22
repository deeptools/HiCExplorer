// chicChicagoBackgroundModel: estimates CHiCAGO's own background model
// (Cairns et al. 2016, Genome Biology 17:127) from CHiCAGO's own input
// files, as a standalone C++-only tool alongside chicViewpointBackgroundModel
// (cpp/PLAN.md 9.15; the CHiCAGO pipeline is exposed as its own
// chicChicago<Stage> tools rather than as a mode flag on the existing chic*
// tools, per the project owner's 2026-09-22 correction to the original 9.15
// design).
//
// Mirrors chicViewpointBackgroundModel's role (fit once, reuse for every
// scoring run) but on CHiCAGO's own design: a .rmap/.baitmap pair, the
// genome-wide design tables (.npb/.nbpb/.poe, produced by CHiCAGO's own
// makeDesignFiles.py) and one .chinput file of observed interaction counts.
//
// Statistics implemented in hicx::chicago (see chicago.hpp/.cpp): readSample
// filtering, addTLB pooling, normaliseBaits/normaliseOtherEnds (non-shrunken,
// Chicago's own default), estimateTechnicalNoise, estimateDistFun (cubic) and
// estimateBrownianComponent's dispersion (MASS::glm.nb's theta, via
// MASS::theta.ml's own Newton iteration).
//
// Scope, honestly stated (see the final report handed to the orchestrating
// session for the measured validation): this tool ingests ONE already-merged
// .chinput file, not several replicate files merged the way R's
// readAndMerge/mergeSamples does; a caller with several replicates sums their
// N columns per (baitID, otherEndID) pair before running this tool.
// estimateBrownianComponent's dispersion is always fit on the full .poe
// design (R's own behaviour when brownianNoise.subset, 1000 baits by
// default, is not smaller than the actual bait count); when the input has
// more baits than that, R itself subsamples baits and averages several
// stochastic glm.nb fits (no fixed seed by default), so its own reported
// dispersion is not reproducible run to run. This tool always reports the
// full-dataset fit instead and flags in its output whether R's subsampling
// would have triggered.

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/chicago.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: chicChicagoBackgroundModel --rmap RMAP --baitmap BAITMAP\n"
    "                                  (--chinput CHINPUT | --matrices MATRICES [MATRICES ...])\n"
    "                                  --nperbin NPERBIN --nbaitsperbin NBAITSPERBIN\n"
    "                                  --proxOE PROXOE [--outFileName OUTFILENAME]\n"
    "                                  [--minFragLen MINFRAGLEN] [--maxFragLen MAXFRAGLEN]\n"
    "                                  [--minNPerBait MINNPERBAIT] [--maxLBrownEst MAXLBROWNEST]\n"
    "                                  [--binsize BINSIZE] [--noRemoveAdjacent]\n"
    "                                  [--tlbFilterTopPercent TLBFILTERTOPPERCENT]\n"
    "                                  [--tlbMinProxOEPerBin TLBMINPROXOEPERBIN]\n"
    "                                  [--tlbMinProxB2BPerBin TLBMINPROXB2BPERBIN]\n"
    "                                  [--techNoiseMinBaitsPerBin TECHNOISEMINBAITSPERBIN]\n"
    "                                  [--brownianNoiseSubset BROWNIANNOISESUBSET] [--help]\n"
    "                                  [--version]\n";

const char* const kHelp =
    "\n"
    "chicChicagoBackgroundModel estimates CHiCAGO's own background model (Cairns et al. 2016) "
    "from CHiCAGO's own .rmap/.baitmap/.chinput input and genome-wide\n"
    "design tables (.npb/.nbpb/.poe): the cubic distance function, the bait and other-end "
    "scaling factors, the negative-binomial dispersion and the Poisson\n"
    "technical-noise table per (bait, other-end) pool. The output feeds chicChicagoScores.\n"
    "\n"
    "This tool has no Python HiCExplorer counterpart; it is C++ only (cpp/PLAN.md 9.15).\n"
    "\n"
    "An example usage is:\n"
    "\n"
    "$ chicChicagoBackgroundModel --rmap design.rmap --baitmap design.baitmap --chinput "
    "sample.chinput --nperbin design.npb --nbaitsperbin design.nbpb --proxOE design.poe "
    "--outFileName background_model.chicago.txt\n"
    "\n"
    "Required arguments:\n"
    "  --rmap RMAP           CHiCAGO .rmap restriction fragment file.\n"
    "  --baitmap BAITMAP     CHiCAGO .baitmap baited fragment file.\n"
    "  --chinput CHINPUT     One CHiCAGO .chinput interaction count file. Several\n"
    "                        replicates must be pre-merged (N summed per bait/\n"
    "                        other-end pair) before this tool; it reads one file.\n"
    "                        Mutually exclusive with --matrices; exactly one of\n"
    "                        the two is required.\n"
    "  --matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]\n"
    "                        One or more Hi-C matrices (cool, h5 or .hic) to derive\n"
    "                        the per (bait, other-end) fragment N counts from\n"
    "                        directly, instead of a pre-built .chinput file. Several\n"
    "                        matrices are summed per fragment pair, the matrix\n"
    "                        equivalent of pre-merging several replicate .chinput\n"
    "                        files. Only cis (bait, other-end) fragment pairs within\n"
    "                        --maxLBrownEst of each other are derived; trans and\n"
    "                        farther-cis counts, used elsewhere for technical-noise\n"
    "                        estimation, are not available this way and still need\n"
    "                        --chinput. Mutually exclusive with --chinput; exactly\n"
    "                        one of the two is required.\n"
    "  --nperbin NPERBIN     CHiCAGO .npb NPerBin design table.\n"
    "  --nbaitsperbin NBAITSPERBIN\n"
    "                        CHiCAGO .nbpb NBaitsPerBin design table.\n"
    "  --proxOE PROXOE       CHiCAGO .poe ProxOE design table.\n"
    "\n"
    "Optional arguments:\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        The name of the background model file (Default:\n"
    "                        chicago_background_model.txt).\n"
    "  --minFragLen MINFRAGLEN\n"
    "                        Minimum other-end fragment length (Default: 150).\n"
    "  --maxFragLen MAXFRAGLEN\n"
    "                        Maximum other-end fragment length (Default: 40000).\n"
    "  --minNPerBait MINNPERBAIT\n"
    "                        Minimum total read count per bait (Default: 250).\n"
    "  --maxLBrownEst MAXLBROWNEST\n"
    "                        Maximum distance the Brownian component is estimated\n"
    "                        over (Default: 1500000).\n"
    "  --binsize BINSIZE     Distance bin width (Default: 20000).\n"
    "  --noRemoveAdjacent    Keep interactions with fragments immediately adjacent\n"
    "                        to their bait (Default: they are removed).\n"
    "  --tlbFilterTopPercent TLBFILTERTOPPERCENT\n"
    "                        Percent of other ends with the highest trans-counts\n"
    "                        dropped before pooling (Default: 0.01).\n"
    "  --tlbMinProxOEPerBin TLBMINPROXOEPERBIN\n"
    "                        Minimum non-bait2bait other ends per trans-count pool\n"
    "                        (Default: 50000).\n"
    "  --tlbMinProxB2BPerBin TLBMINPROXB2BPERBIN\n"
    "                        Minimum bait2bait other ends per trans-count pool\n"
    "                        (Default: 2500).\n"
    "  --techNoiseMinBaitsPerBin TECHNOISEMINBAITSPERBIN\n"
    "                        Minimum baits per technical-noise pool (Default: 1000).\n"
    "  --brownianNoiseSubset BROWNIANNOISESUBSET\n"
    "                        Reported only: the bait count above which R's own\n"
    "                        estimateBrownianComponent subsamples (Default: 1000).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

}  // namespace

int main(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("chicChicagoBackgroundModel",
                       "Estimates CHiCAGO's own background model from CHiCAGO's own input files.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--rmap"}).required().input({"rmap", "txt"}).help("CHiCAGO .rmap file.");
    required.add({"--baitmap"}).required().input({"baitmap", "txt"}).help("CHiCAGO .baitmap file.");
    cli::MutuallyExclusiveGroup& input_source = parser.mutually_exclusive(required, /*required=*/true);
    input_source.add({"--chinput"}).input({"chinput", "txt"}).help("CHiCAGO .chinput file.");
    input_source.add({"--matrices", "-m"})
        .nargs("+")
        .input({"cool", "h5", "hic"})
        .help("Hi-C matrices to derive per-fragment-pair counts from directly.");
    required.add({"--nperbin"}).required().input({"npb", "txt"}).help("CHiCAGO .npb design table.");
    required.add({"--nbaitsperbin"})
        .required()
        .input({"nbpb", "txt"})
        .help("CHiCAGO .nbpb design table.");
    required.add({"--proxOE"}).required().input({"poe", "txt"}).help("CHiCAGO .poe design table.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--outFileName", "-o"})
        .default_value("chicago_background_model.txt")
        .output({"txt"})
        .help("The name of the background model file.");
    optional.add({"--minFragLen"}).type("int").default_value(150).help("Minimum other-end length.");
    optional.add({"--maxFragLen"}).type("int").default_value(40000).help("Maximum other-end length.");
    optional.add({"--minNPerBait"}).type("int").default_value(250).help("Minimum reads per bait.");
    optional.add({"--maxLBrownEst"})
        .type("int")
        .default_value(1500000)
        .help("Brownian estimation distance cutoff.");
    optional.add({"--binsize"}).type("int").default_value(20000).help("Distance bin width.");
    optional.add({"--noRemoveAdjacent"})
        .action(cli::Action::StoreTrue)
        .help("Keep bait-adjacent interactions.");
    optional.add({"--tlbFilterTopPercent"})
        .type("float")
        .default_value(0.01)
        .help("Top trans-count percent dropped before pooling.");
    optional.add({"--tlbMinProxOEPerBin"})
        .type("int")
        .default_value(50000)
        .help("Minimum non-bait2bait other ends per pool.");
    optional.add({"--tlbMinProxB2BPerBin"})
        .type("int")
        .default_value(2500)
        .help("Minimum bait2bait other ends per pool.");
    optional.add({"--techNoiseMinBaitsPerBin"})
        .type("int")
        .default_value(1000)
        .help("Minimum baits per technical-noise pool.");
    optional.add({"--brownianNoiseSubset"})
        .type("int")
        .default_value(1000)
        .help("Bait count above which R subsamples (reported only).");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace args = parser.parse(argc, argv);

    try {
        using namespace hicx::chicago;
        FilterSettings fs;
        fs.min_frag_len = args.integer("minFragLen");
        fs.max_frag_len = args.integer("maxFragLen");
        fs.min_n_per_bait = args.integer("minNPerBait");
        fs.max_l_brown_est = args.integer("maxLBrownEst");
        fs.binsize = args.integer("binsize");
        fs.remove_adjacent = !args.flag("noRemoveAdjacent");

        auto rmap = read_rmap(args.str("rmap"));
        auto baitmap = read_baitmap(args.str("baitmap"));
        const std::vector<ChinputRecord> raw = args.given("chinput")
            ? read_chinput(args.str("chinput"))
            : chinput_from_matrices(args.strs("matrices"), rmap, baitmap, fs);

        const BackgroundModel model = fit_chicago_background(
            raw, rmap, baitmap, args.str("nperbin"), args.str("nbaitsperbin"), args.str("proxOE"), fs,
            args.integer("techNoiseMinBaitsPerBin"), args.real("tlbFilterTopPercent"),
            args.integer("tlbMinProxOEPerBin"), args.integer("tlbMinProxB2BPerBin"),
            args.integer("brownianNoiseSubset"));

        std::ofstream out(args.str("outFileName"), std::ios::binary);
        if (!out) {
            throw std::runtime_error("[Errno 2] No such file or directory: '" + args.str("outFileName") +
                                     "'");
        }
        out.precision(15);
        out << "# chicChicagoBackgroundModel output (hicexplorer v4, CHiCAGO background model, "
               "PLAN.md 9.15)\n";
        out << "[Settings]\n";
        out << "minFragLen\t" << fs.min_frag_len << "\n";
        out << "maxFragLen\t" << fs.max_frag_len << "\n";
        out << "minNPerBait\t" << fs.min_n_per_bait << "\n";
        out << "maxLBrownEst\t" << fs.max_l_brown_est << "\n";
        out << "binsize\t" << fs.binsize << "\n";
        out << "removeAdjacent\t" << (fs.remove_adjacent ? 1 : 0) << "\n";
        out << "nNonB2BPools\t" << model.tlb.n_non_b2b_pools << "\n";
        out << "dispersionNPairs\t" << model.dispersion_n_pairs << "\n";
        out << "subsetWouldTriggerInR\t" << (model.subset_would_trigger_in_r ? 1 : 0) << "\n";

        out << "[DistFun]\n";
        out << "cubic0\t" << model.dist_fun.cubic[0] << "\n";
        out << "cubic1\t" << model.dist_fun.cubic[1] << "\n";
        out << "cubic2\t" << model.dist_fun.cubic[2] << "\n";
        out << "cubic3\t" << model.dist_fun.cubic[3] << "\n";
        out << "obsMinLog\t" << model.dist_fun.obs_min_log << "\n";
        out << "obsMaxLog\t" << model.dist_fun.obs_max_log << "\n";
        out << "headA\t" << model.dist_fun.head_coef[0] << "\n";
        out << "headB\t" << model.dist_fun.head_coef[1] << "\n";
        out << "tailA\t" << model.dist_fun.tail_coef[0] << "\n";
        out << "tailB\t" << model.dist_fun.tail_coef[1] << "\n";

        out << "[Dispersion]\n";
        out << "alpha\t" << model.dispersion << "\n";

        out << "[BaitFactors]\n";
        out << "baitID\ts_j\n";
        for (const auto& [bait_id, s_j] : model.bait_factors.s_j) {
            out << bait_id << "\t" << s_j << "\n";
        }

        out << "[OtherEndPools]\n";
        out << "otherEndID\tpool\n";
        for (const auto& [oe_id, pool] : model.tlb.pool_of_other_end) {
            out << oe_id << "\t" << pool << "\n";
        }

        out << "[PoolFactors]\n";
        out << "pool\ts_i\n";
        for (const auto& [pool, s_i] : model.s_i_by_tlb_pool) {
            out << pool << "\t" << s_i << "\n";
        }

        out << "[BaitPools]\n";
        out << "baitID\ttblb\n";
        for (const auto& [bait_id, tblb] : model.tech_noise.tblb_of_bait) {
            out << bait_id << "\t" << tblb << "\n";
        }

        out << "[TechnicalNoise]\n";
        out << "tlb\ttblb\tTmean\n";
        for (const auto& [key, tmean] : model.tech_noise.tmean_by_pool) {
            out << key.first << "\t" << key.second << "\t" << tmean << "\n";
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chicChicagoBackgroundModel: %s\n", error.what());
        return 1;
    }
    return 0;
}
