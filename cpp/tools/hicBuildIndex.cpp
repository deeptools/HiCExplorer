// hicBuildIndex: PLAN.md tier 13, "custom index generation". No Python
// counterpart, so its own name follows this project's tool naming rather
// than porting a Python module.
//
// A thin wrapper around `minibwa index`, the conda-packaged aligner at
// ~/miniconda3/envs/__minibwa@0.5/bin/minibwa (HICX_MINIBWA_BIN overrides
// the path; unset, "minibwa" is looked up on PATH). It exists so a user can
// build a custom minibwa index from any FASTA -- not only the baked-in
// genomes a later Docker-image task adds (PLAN.md tier 13 explicitly scopes
// that out of this tool) -- from the same C++ tool directory and the same
// GUI catalog as every other tool, with the index prefix stored per project
// or in the shared cache directory the GUI's Settings dialog exposes
// (gui/hicexplorer_gui/settings.py, index_cache_dir).
//
// The wrapper adds nothing of its own: every option below is minibwa
// index's own option, passed straight through, and --outPrefix is
// minibwa's own optional positional out.prefix argument, made an explicit,
// named, required argument here because this project's file-role
// convention (hicx/argparse.hpp) needs a single argument to hang the
// "output, kind=prefix" tag on for the GUI form and for
// cpp/scripts/equiv.py's file discovery.
//
// Validation: cpp/PLAN.md 5.1 class E0 (byte-identical), against running
// `minibwa index` directly on the same FASTA with the same options --
// there is no algorithm here to diverge from, so anything other than
// identical bytes would be a bug in how this wrapper builds minibwa's own
// command line.

#include <cstdio>

#include "hicx/argparse.hpp"
#include "hicx/minibwa_bridge.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicBuildIndex --inFile FASTA --outPrefix PREFIX\n"
    "                     [--seed SEED] [--saSampleRate SASAMPLERATE]\n"
    "                     [--lowMemory] [--blockSize BLOCKSIZE]\n"
    "                     [--threads THREADS] [--methylation] [--help]\n"
    "                     [--version]\n";

const char* const kHelp =
    "\n"
    "New in HiCExplorer v4 (PLAN.md tier 13). Builds a minibwa index from a\n"
    "FASTA reference, for hicAlignReads and for hicBuildMatrix /\n"
    "hicBuildMatrixMicroC once the resulting BAM is written. A thin wrapper\n"
    "around `minibwa index`; every option below is minibwa's own option.\n"
    "\n"
    "Required arguments:\n"
    "  --inFile FASTA, -i FASTA\n"
    "                        The reference FASTA to index.\n"
    "  --outPrefix PREFIX, -o PREFIX\n"
    "                        Prefix of the index files minibwa writes\n"
    "                        (minibwa's own out.prefix argument).\n"
    "\n"
    "Optional arguments:\n"
    "  --seed SEED           Random seed for ambiguous bases (minibwa's -s).\n"
    "                        (Default: 11).\n"
    "  --saSampleRate SASAMPLERATE\n"
    "                        SA sample rate at 1/(1<<INT) (minibwa's -u).\n"
    "                        (Default: 4).\n"
    "  --lowMemory           Low-memory GPL'd BWT construction algorithm\n"
    "                        (minibwa's -l). (Default: False).\n"
    "  --blockSize BLOCKSIZE\n"
    "                        Block size, effective with --lowMemory\n"
    "                        (minibwa's -b). (Default: 10m, minibwa's own\n"
    "                        syntax, e.g. \"10m\").\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads, effective without\n"
    "                        --lowMemory (minibwa's -t). (Default: 4).\n"
    "  --methylation         Build an FM-index for BS-seq mapping\n"
    "                        (minibwa's --meth). (Default: False).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::string in_file;
    std::string out_prefix;
    std::int64_t seed = 11;
    std::int64_t sa_sample_rate = 4;
    bool low_memory = false;
    std::string block_size = "10m";
    bool block_size_given = false;
    std::int64_t threads = 4;
    bool methylation = false;
};

Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicBuildIndex", "Builds a minibwa index from a FASTA reference.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--inFile", "-i"})
        .required()
        .input({"fasta"})
        .help("The reference FASTA to index.");
    required.add({"--outPrefix", "-o"})
        .required()
        .output({"bwaidx"}, "prefix")
        .help("Prefix of the index files minibwa writes.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--seed"}).type("int").default_value(11).help(
        "Random seed for ambiguous bases.");
    optional.add({"--saSampleRate"}).type("int").default_value(4).help(
        "SA sample rate at 1/(1<<INT).");
    optional.add({"--lowMemory"}).action(cli::Action::StoreTrue).help(
        "Low-memory GPL'd BWT construction algorithm.");
    optional.add({"--blockSize"}).default_value("10m").help(
        "Block size, effective with --lowMemory.");
    optional.add({"--threads", "-t"}).type("int").default_value(4).help(
        "Number of threads, effective without --lowMemory.");
    optional.add({"--methylation"}).action(cli::Action::StoreTrue).help(
        "Build an FM-index for BS-seq mapping.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.in_file = ns.str("inFile");
    args.out_prefix = ns.str("outPrefix");
    args.seed = ns.integer("seed");
    args.sa_sample_rate = ns.integer("saSampleRate");
    args.low_memory = ns.flag("lowMemory");
    args.block_size = ns.str("blockSize");
    args.block_size_given = ns.given("blockSize");
    args.threads = ns.integer("threads");
    args.methylation = ns.flag("methylation");
    return args;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);

    std::vector<std::string> minibwa_args = {"index"};
    minibwa_args.push_back("-s");
    minibwa_args.push_back(std::to_string(args.seed));
    minibwa_args.push_back("-u");
    minibwa_args.push_back(std::to_string(args.sa_sample_rate));
    if (args.low_memory) {
        minibwa_args.push_back("-l");
    }
    if (args.block_size_given) {
        minibwa_args.push_back("-b");
        minibwa_args.push_back(args.block_size);
    }
    minibwa_args.push_back("-t");
    minibwa_args.push_back(std::to_string(args.threads));
    if (args.methylation) {
        minibwa_args.push_back("--meth");
    }
    minibwa_args.push_back(args.in_file);
    minibwa_args.push_back(args.out_prefix);

    const int status = hicx::minibwa::run("hicBuildIndex", minibwa_args);
    if (status != 0) {
        return status;
    }
    hicx::report_resource_usage("hicBuildIndex");
    return 0;
}
