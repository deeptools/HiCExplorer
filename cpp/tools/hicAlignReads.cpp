// hicAlignReads: PLAN.md tier 13, the FASTQ-to-BAM alignment step that runs
// ahead of hicBuildMatrix / hicBuildMatrixMicroC. No Python counterpart, so
// its own name follows this project's tool naming rather than porting a
// Python module.
//
// A thin wrapper around `minibwa map`, the conda-packaged aligner at
// ~/miniconda3/envs/__minibwa@0.5/bin/minibwa (HICX_MINIBWA_BIN overrides
// the path; unset, "minibwa" is looked up on PATH). minibwa map aligns one
// FASTQ file against one index and writes SAM (its `map [options] <in.idx>
// <in.fastq>`, no paired-input syntax); this wrapper adds a BAM conversion
// step (hicx::minibwa::run_map_to_bam, core/src/minibwa_bridge.cpp) so the
// output is what hicBuildMatrix's --samFiles actually reads, and nothing
// else. It is run once per mate: hicBuildMatrix already takes two
// independently aligned BAMs, one per mate (its own --samFiles takes
// "two sam files"), which is standard Hi-C practice because a ligation
// junction read is not a simple concordant pair; minibwa's own --hic flag
// ("map Hi-C reads; equivalent to option -5P": take the alignment with the
// smallest query position as primary, and skip pairing and mate rescue) is
// built for exactly this per-mate, non-paired alignment mode, so
// hicAlignReads passes it through rather than inventing any pairing
// behaviour of its own.
//
// Every option below is minibwa map's own option, passed straight through.
// Not exposed, all deliberately, because they do not fit this wrapper's
// fixed SAM-to-BAM output path or add nothing for the Hi-C use case:
//   -f (PAF output instead of SAM) -- the wrapper always converts to BAM;
//   -o (output file) -- fixed to a private temporary file internally;
//   --mmap's "=lite" argument -- exposed as a plain on/off flag here
//     (--memoryMap), since the "lite" variant is a memory/speed trade-off
//     minibwa's own --help does not explain further and no test in this
//     project exercises;
//   --help/--version of minibwa itself -- this tool has its own.
//
// Validation: cpp/PLAN.md 5.1 class E0 (byte-identical), against running
// `minibwa map ... -o x.sam` directly on the same FASTQ and index and then
// converting with `samtools view -b --no-PG x.sam -o x.bam`. Both routes go
// through the same htslib (the conda env's htslib 1.21, which samtools
// 1.21 also links, cpp/AGENTS_CONTRACT.md section 2). Measured on a 500-read
// synthetic FASTQ against a 20 kb synthetic reference at --threads 1 and 4:
// `samtools view` reports byte-for-byte identical alignment records at both
// thread counts, and the SAM/BAM headers agree except for one thing that
// cannot agree by construction -- the @PG CL line minibwa itself writes
// quotes the SAM path it was told to write to (its own -o argument), and
// this wrapper always passes a freshly generated private temporary path
// there, never the same string twice. That is the file-format analogue of
// PLAN.md 5.1 class E1's own normalisation of "creation-date, generated-by
// ... tool-url": a provenance field that legitimately differs by
// construction, not a divergence in what was computed. --no-PG on the
// reference's samtools step is needed for the comparison itself: without
// it, samtools stamps its own extra @PG record on conversion, which this
// wrapper's direct htslib conversion (hicx::minibwa::run_map_to_bam) never
// adds, so an unqualified `cmp` against plain `samtools view -b` is not the
// right reference command; `--no-PG` makes the two routes actually
// comparable. Full commands and results are in this task's report.

#include <cstdio>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/minibwa_bridge.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicAlignReads --inFile FASTQ --index PREFIX --outFileName BAM\n"
    "                     [--hic] [--methylation] [--threads THREADS]\n"
    "                     [--shortReadThreshold SHORTREADTHRESHOLD]\n"
    "                     [--readGroup READGROUP] [--baseTag {cs,ds,MD}]\n"
    "                     [--minSeedLength MINSEEDLENGTH]\n"
    "                     [--maxSeedOccurrences MAXSEEDOCCURRENCES]\n"
    "                     [--maxGapSize MAXGAPSIZE] [--bandwidth BANDWIDTH]\n"
    "                     [--longBandwidth LONGBANDWIDTH]\n"
    "                     [--minChainingScore MINCHAININGSCORE]\n"
    "                     [--minSecondaryRatio MINSECONDARYRATIO]\n"
    "                     [--maxSecondary MAXSECONDARY] [--chainOnly]\n"
    "                     [--preset {sr,lr,adap}] [--matchScore MATCHSCORE]\n"
    "                     [--mismatchPenalty MISMATCHPENALTY]\n"
    "                     [--gapOpenPenalty GAPOPENPENALTY]\n"
    "                     [--gapExtendPenalty GAPEXTENDPENALTY]\n"
    "                     [--minDpScore MINDPSCORE] [--skipPairing]\n"
    "                     [--mateRescue MATERESCUE]\n"
    "                     [--isizeDistribution ISIZEDISTRIBUTION]\n"
    "                     [--noUnmapped] [--outSecondaryCount OUTSECONDARYCOUNT]\n"
    "                     [--outSecondaryScore OUTSECONDARYSCORE]\n"
    "                     [--xaThreshold XATHRESHOLD] [--copyComments]\n"
    "                     [--softClipSupplementary] [--headerInsert HEADERINSERT]\n"
    "                     [--smallestPrimary] [--batchSize BATCHSIZE]\n"
    "                     [--memoryMap] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "New in HiCExplorer v4 (PLAN.md tier 13). Aligns one FASTQ file against a\n"
    "minibwa index and writes a BAM ready for hicBuildMatrix's --samFiles (run\n"
    "once per mate: hicBuildMatrix reads two independently aligned BAMs, one\n"
    "per mate). A thin wrapper around `minibwa map`; every option below is\n"
    "minibwa's own option, except --index and --outFileName, which name this\n"
    "wrapper's fixed BAM output path.\n"
    "\n"
    "Required arguments:\n"
    "  --inFile FASTQ, -i FASTQ\n"
    "                        The FASTQ file to align (one mate).\n"
    "  --index PREFIX        Prefix of a minibwa index (hicBuildIndex's\n"
    "                        --outPrefix).\n"
    "  --outFileName BAM, -o BAM\n"
    "                        Output BAM file.\n"
    "\n"
    "Optional arguments:\n"
    "  --hic                 Map Hi-C reads: take the alignment with the\n"
    "                        smallest query position as primary and skip\n"
    "                        pairing and mate rescue (minibwa's --hic, equal\n"
    "                        to -5P). (Default: False).\n"
    "  --methylation         Map directional bisulfite sequencing reads\n"
    "                        (minibwa's --meth). (Default: False).\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of worker threads. (Default: 1).\n"
    "  --shortReadThreshold SHORTREADTHRESHOLD\n"
    "                        Treat reads shorter than this as short reads in\n"
    "                        the default adaptive mode (minibwa's -l).\n"
    "                        (Default: 325).\n"
    "  --readGroup READGROUP\n"
    "                        SAM read group line, e.g. '@RG\\tID:foo\\tSM:bar'\n"
    "                        (minibwa's -R).\n"
    "  --baseTag {cs,ds,MD}  Output a base alignment tag (minibwa's -b).\n"
    "  --minSeedLength MINSEEDLENGTH\n"
    "                        Minimum seed length (minibwa's -k).\n"
    "                        (Default: 19).\n"
    "  --maxSeedOccurrences MAXSEEDOCCURRENCES\n"
    "                        Maximum seed occurrences (minibwa's -c).\n"
    "                        (Default: 250).\n"
    "  --maxGapSize MAXGAPSIZE\n"
    "                        Maximum gap size, controlling extension and\n"
    "                        chain breaking (minibwa's -g). (Default: 100).\n"
    "  --bandwidth BANDWIDTH\n"
    "                        Bandwidth (minibwa's -w). (Default: 100).\n"
    "  --longBandwidth LONGBANDWIDTH\n"
    "                        Long bandwidth, for long reads or the adaptive\n"
    "                        mode (minibwa's -W). (Default: 20000).\n"
    "  --minChainingScore MINCHAININGSCORE\n"
    "                        Minimum chaining score (minibwa's -m).\n"
    "                        (Default: 25).\n"
    "  --minSecondaryRatio MINSECONDARYRATIO\n"
    "                        Minimum secondary-to-primary score ratio\n"
    "                        (minibwa's -p). (Default: 0.5).\n"
    "  --maxSecondary MAXSECONDARY\n"
    "                        Retain at most this many secondary alignments\n"
    "                        (minibwa's -N). (Default: 50).\n"
    "  --chainOnly           Perform chaining only, without base alignment\n"
    "                        (minibwa's --chain-only). (Default: False).\n"
    "  --preset {sr,lr,adap}\n"
    "                        Preset: short reads, long reads, or mixed\n"
    "                        short/long (minibwa's -x). (Default: adap).\n"
    "  --matchScore MATCHSCORE\n"
    "                        Matching score (minibwa's -A). (Default: 2).\n"
    "  --mismatchPenalty MISMATCHPENALTY\n"
    "                        Mismatch penalty (minibwa's -B). (Default: 8).\n"
    "  --gapOpenPenalty GAPOPENPENALTY\n"
    "                        Gap open penalty, INT1[,INT2] (minibwa's -O).\n"
    "                        (Default: 12,23).\n"
    "  --gapExtendPenalty GAPEXTENDPENALTY\n"
    "                        Gap extension penalty, INT1[,INT2] (minibwa's\n"
    "                        -E). (Default: 2,1).\n"
    "  --minDpScore MINDPSCORE\n"
    "                        Suppress alignments with a DP score lower than\n"
    "                        this times the matching score (minibwa's -s).\n"
    "                        (Default: 30).\n"
    "  --skipPairing         Skip pairing and mate rescue (minibwa's -P).\n"
    "                        (Default: False).\n"
    "  --mateRescue MATERESCUE\n"
    "                        Mate rescue for up to this many candidates, 0 to\n"
    "                        skip (minibwa's --rescue). (Default: 10).\n"
    "  --isizeDistribution ISIZEDISTRIBUTION\n"
    "                        Mean, stddev, max and min of the insert size\n"
    "                        distribution, INT[,INT[,INT[,INT]]] (minibwa's\n"
    "                        -I). (Default: inferred).\n"
    "  --noUnmapped          Do not output unmapped reads (minibwa's -u).\n"
    "                        (Default: False).\n"
    "  --outSecondaryCount OUTSECONDARYCOUNT\n"
    "                        Output up to this many secondary alignments\n"
    "                        (minibwa's --outn). (Default: 0).\n"
    "  --outSecondaryScore OUTSECONDARYSCORE\n"
    "                        Output a secondary hit if its score is at least\n"
    "                        this fraction of the best score (minibwa's\n"
    "                        --outs). (Default: 0.8).\n"
    "  --xaThreshold XATHRESHOLD\n"
    "                        If at most this many hits score above 80%% of\n"
    "                        the best hit, output them to the XA tag\n"
    "                        (minibwa's --xa). (Default: 5).\n"
    "  --copyComments        Copy FASTA/Q comments to the output (minibwa's\n"
    "                        -y). (Default: False).\n"
    "  --softClipSupplementary\n"
    "                        Use soft clipping for supplementary alignments\n"
    "                        (minibwa's -Y). (Default: False).\n"
    "  --headerInsert HEADERINSERT\n"
    "                        A string starting with @ to insert into the\n"
    "                        header, or a file of lines to insert (minibwa's\n"
    "                        -H).\n"
    "  --smallestPrimary     Take the alignment with the smallest query\n"
    "                        position as primary (minibwa's -5), implied by\n"
    "                        --hic. (Default: False).\n"
    "  --batchSize BATCHSIZE\n"
    "                        Process NUM1-NUM2 bp of query sequences per\n"
    "                        batch, NUM1[,NUM2] (minibwa's -K). (Default:\n"
    "                        100m,1g).\n"
    "  --memoryMap           Load the index via memory mapped files, slower\n"
    "                        mapping (minibwa's --mmap). (Default: False).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::string in_file;
    std::string index;
    std::string out_file_name;
    bool hic = false;
    bool methylation = false;
    std::int64_t threads = 1;
    std::optional<std::int64_t> short_read_threshold;
    std::optional<std::string> read_group;
    std::optional<std::string> base_tag;
    std::optional<std::int64_t> min_seed_length;
    std::optional<std::int64_t> max_seed_occurrences;
    std::optional<std::int64_t> max_gap_size;
    std::optional<std::int64_t> bandwidth;
    std::optional<std::int64_t> long_bandwidth;
    std::optional<std::int64_t> min_chaining_score;
    std::optional<double> min_secondary_ratio;
    std::optional<std::int64_t> max_secondary;
    bool chain_only = false;
    std::optional<std::string> preset;
    std::optional<std::int64_t> match_score;
    std::optional<std::int64_t> mismatch_penalty;
    std::optional<std::string> gap_open_penalty;
    std::optional<std::string> gap_extend_penalty;
    std::optional<std::int64_t> min_dp_score;
    bool skip_pairing = false;
    std::optional<std::int64_t> mate_rescue;
    std::optional<std::string> isize_distribution;
    bool no_unmapped = false;
    std::optional<std::int64_t> out_secondary_count;
    std::optional<double> out_secondary_score;
    std::optional<std::int64_t> xa_threshold;
    bool copy_comments = false;
    bool soft_clip_supplementary = false;
    std::optional<std::string> header_insert;
    bool smallest_primary = false;
    std::optional<std::string> batch_size;
    bool memory_map = false;
};

Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicAlignReads",
                       "Aligns one FASTQ file against a minibwa index and writes a BAM.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--inFile", "-i"})
        .required()
        .input({"fastq", "fastq.gz"})
        .help("The FASTQ file to align (one mate).");
    required.add({"--index"})
        .required()
        .input({"bwaidx"}, "prefix")
        .help("Prefix of a minibwa index (hicBuildIndex's --outPrefix).");
    required.add({"--outFileName", "-o"})
        .required()
        .output({"bam"})
        .help("Output BAM file.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--hic"}).action(cli::Action::StoreTrue).help(
        "Map Hi-C reads (equal to minibwa's -5P).");
    optional.add({"--methylation"}).action(cli::Action::StoreTrue).help(
        "Map directional bisulfite sequencing reads.");
    optional.add({"--threads", "-t"}).type("int").default_value(1).help("Number of worker threads.");
    optional.add({"--shortReadThreshold"}).type("int").help(
        "Treat reads shorter than this as short reads in the default adaptive mode.");
    optional.add({"--readGroup"}).help("SAM read group line.");
    optional.add({"--baseTag"}).choices({"cs", "ds", "MD"}).help("Output a base alignment tag.");
    optional.add({"--minSeedLength"}).type("int").help("Minimum seed length.");
    optional.add({"--maxSeedOccurrences"}).type("int").help("Maximum seed occurrences.");
    optional.add({"--maxGapSize"}).type("int").help(
        "Maximum gap size, controlling extension and chain breaking.");
    optional.add({"--bandwidth"}).type("int").help("Bandwidth.");
    optional.add({"--longBandwidth"}).type("int").help(
        "Long bandwidth, for long reads or the adaptive mode.");
    optional.add({"--minChainingScore"}).type("int").help("Minimum chaining score.");
    optional.add({"--minSecondaryRatio"}).type("float").help(
        "Minimum secondary-to-primary score ratio.");
    optional.add({"--maxSecondary"}).type("int").help("Retain at most this many secondary alignments.");
    optional.add({"--chainOnly"}).action(cli::Action::StoreTrue).help(
        "Perform chaining only, without base alignment.");
    optional.add({"--preset"}).choices({"sr", "lr", "adap"}).help("Alignment preset.");
    optional.add({"--matchScore"}).type("int").help("Matching score.");
    optional.add({"--mismatchPenalty"}).type("int").help("Mismatch penalty.");
    optional.add({"--gapOpenPenalty"}).help("Gap open penalty, INT1[,INT2].");
    optional.add({"--gapExtendPenalty"}).help("Gap extension penalty, INT1[,INT2].");
    optional.add({"--minDpScore"}).type("int").help(
        "Suppress alignments with a DP score lower than this times the matching score.");
    optional.add({"--skipPairing"}).action(cli::Action::StoreTrue).help("Skip pairing and mate rescue.");
    optional.add({"--mateRescue"}).type("int").help(
        "Mate rescue for up to this many candidates, 0 to skip.");
    optional.add({"--isizeDistribution"}).help(
        "Mean, stddev, max and min of the insert size distribution, INT[,INT[,INT[,INT]]].");
    optional.add({"--noUnmapped"}).action(cli::Action::StoreTrue).help("Do not output unmapped reads.");
    optional.add({"--outSecondaryCount"}).type("int").help("Output up to this many secondary alignments.");
    optional.add({"--outSecondaryScore"}).type("float").help(
        "Output a secondary hit if its score is at least this fraction of the best score.");
    optional.add({"--xaThreshold"}).type("int").help(
        "If at most this many hits score above 80%% of the best hit, output them to the XA tag.");
    optional.add({"--copyComments"}).action(cli::Action::StoreTrue).help(
        "Copy FASTA/Q comments to the output.");
    optional.add({"--softClipSupplementary"}).action(cli::Action::StoreTrue).help(
        "Use soft clipping for supplementary alignments.");
    optional.add({"--headerInsert"}).help(
        "A string starting with @ to insert into the header, or a file of lines to insert.");
    optional.add({"--smallestPrimary"}).action(cli::Action::StoreTrue).help(
        "Take the alignment with the smallest query position as primary, implied by --hic.");
    optional.add({"--batchSize"}).help("Process NUM1-NUM2 bp of query sequences per batch, NUM1[,NUM2].");
    optional.add({"--memoryMap"}).action(cli::Action::StoreTrue).help(
        "Load the index via memory mapped files, slower mapping.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.in_file = ns.str("inFile");
    args.index = ns.str("index");
    args.out_file_name = ns.str("outFileName");
    args.hic = ns.flag("hic");
    args.methylation = ns.flag("methylation");
    args.threads = ns.integer("threads");
    if (ns.given("shortReadThreshold")) args.short_read_threshold = ns.integer("shortReadThreshold");
    if (ns.given("readGroup")) args.read_group = ns.str("readGroup");
    if (ns.given("baseTag")) args.base_tag = ns.str("baseTag");
    if (ns.given("minSeedLength")) args.min_seed_length = ns.integer("minSeedLength");
    if (ns.given("maxSeedOccurrences")) args.max_seed_occurrences = ns.integer("maxSeedOccurrences");
    if (ns.given("maxGapSize")) args.max_gap_size = ns.integer("maxGapSize");
    if (ns.given("bandwidth")) args.bandwidth = ns.integer("bandwidth");
    if (ns.given("longBandwidth")) args.long_bandwidth = ns.integer("longBandwidth");
    if (ns.given("minChainingScore")) args.min_chaining_score = ns.integer("minChainingScore");
    if (ns.given("minSecondaryRatio")) args.min_secondary_ratio = ns.real("minSecondaryRatio");
    if (ns.given("maxSecondary")) args.max_secondary = ns.integer("maxSecondary");
    args.chain_only = ns.flag("chainOnly");
    if (ns.given("preset")) args.preset = ns.str("preset");
    if (ns.given("matchScore")) args.match_score = ns.integer("matchScore");
    if (ns.given("mismatchPenalty")) args.mismatch_penalty = ns.integer("mismatchPenalty");
    if (ns.given("gapOpenPenalty")) args.gap_open_penalty = ns.str("gapOpenPenalty");
    if (ns.given("gapExtendPenalty")) args.gap_extend_penalty = ns.str("gapExtendPenalty");
    if (ns.given("minDpScore")) args.min_dp_score = ns.integer("minDpScore");
    args.skip_pairing = ns.flag("skipPairing");
    if (ns.given("mateRescue")) args.mate_rescue = ns.integer("mateRescue");
    if (ns.given("isizeDistribution")) args.isize_distribution = ns.str("isizeDistribution");
    args.no_unmapped = ns.flag("noUnmapped");
    if (ns.given("outSecondaryCount")) args.out_secondary_count = ns.integer("outSecondaryCount");
    if (ns.given("outSecondaryScore")) args.out_secondary_score = ns.real("outSecondaryScore");
    if (ns.given("xaThreshold")) args.xa_threshold = ns.integer("xaThreshold");
    args.copy_comments = ns.flag("copyComments");
    args.soft_clip_supplementary = ns.flag("softClipSupplementary");
    if (ns.given("headerInsert")) args.header_insert = ns.str("headerInsert");
    args.smallest_primary = ns.flag("smallestPrimary");
    if (ns.given("batchSize")) args.batch_size = ns.str("batchSize");
    args.memory_map = ns.flag("memoryMap");
    return args;
}

// A two-token option: "-k" "19".
template <typename T>
void push_opt(std::vector<std::string>& out, const char* flag, const std::optional<T>& value) {
    if (value.has_value()) {
        out.push_back(flag);
        if constexpr (std::is_same_v<T, std::string>) {
            out.push_back(*value);
        } else {
            out.push_back(std::to_string(*value));
        }
    }
}

// A one-token, '='-attached option: "--rescue=10". minibwa's own --help
// shows --rescue, --outn, --outs and --xa exactly this way, unlike its
// other long and short options.
template <typename T>
void push_opt_eq(std::vector<std::string>& out, const char* flag, const std::optional<T>& value) {
    if (value.has_value()) {
        std::string token = flag;
        if constexpr (std::is_same_v<T, std::string>) {
            token += *value;
        } else {
            token += std::to_string(*value);
        }
        out.push_back(token);
    }
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);

    std::vector<std::string> map_args;
    map_args.push_back("-t");
    map_args.push_back(std::to_string(args.threads));
    if (args.hic) map_args.push_back("--hic");
    if (args.methylation) map_args.push_back("--meth");
    push_opt(map_args, "-l", args.short_read_threshold);
    push_opt(map_args, "-R", args.read_group);
    push_opt(map_args, "-b", args.base_tag);
    push_opt(map_args, "-k", args.min_seed_length);
    push_opt(map_args, "-c", args.max_seed_occurrences);
    push_opt(map_args, "-g", args.max_gap_size);
    push_opt(map_args, "-w", args.bandwidth);
    push_opt(map_args, "-W", args.long_bandwidth);
    push_opt(map_args, "-m", args.min_chaining_score);
    push_opt(map_args, "-p", args.min_secondary_ratio);
    push_opt(map_args, "-N", args.max_secondary);
    if (args.chain_only) map_args.push_back("--chain-only");
    push_opt(map_args, "-x", args.preset);
    push_opt(map_args, "-A", args.match_score);
    push_opt(map_args, "-B", args.mismatch_penalty);
    push_opt(map_args, "-O", args.gap_open_penalty);
    push_opt(map_args, "-E", args.gap_extend_penalty);
    push_opt(map_args, "-s", args.min_dp_score);
    if (args.skip_pairing) map_args.push_back("-P");
    push_opt_eq(map_args, "--rescue=", args.mate_rescue);
    push_opt(map_args, "-I", args.isize_distribution);
    if (args.no_unmapped) map_args.push_back("-u");
    push_opt_eq(map_args, "--outn=", args.out_secondary_count);
    push_opt_eq(map_args, "--outs=", args.out_secondary_score);
    push_opt_eq(map_args, "--xa=", args.xa_threshold);
    if (args.copy_comments) map_args.push_back("-y");
    if (args.soft_clip_supplementary) map_args.push_back("-Y");
    push_opt(map_args, "-H", args.header_insert);
    if (args.smallest_primary) map_args.push_back("-5");
    push_opt(map_args, "-K", args.batch_size);
    if (args.memory_map) map_args.push_back("--mmap");

    const int status = hicx::minibwa::run_map_to_bam("hicAlignReads", map_args, args.index, args.in_file,
                                                      args.out_file_name);
    if (status != 0) {
        return status;
    }
    hicx::report_resource_usage("hicAlignReads");
    return 0;
}
