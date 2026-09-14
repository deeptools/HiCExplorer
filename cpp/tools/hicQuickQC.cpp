// Port of hicexplorer/hicQuickQC.py.
//
// hicQuickQC is a thin front end: it rewrites its arguments into a
// hicBuildMatrix command line with --doTestRun and calls hicBuildMatrix.main.
// The port does the same against the shared BAM subsystem in
// build_matrix_impl.hpp, so the quality control counters, QC.log and the five
// *_table.txt files come from exactly the code hicBuildMatrix runs.
//
// What hicQuickQC.py:100-127 hands to hicBuildMatrix, and therefore what is
// fixed here rather than exposed:
//
//   --outFileName <NamedTemporaryFile .h5>  --QCfolder <QCfolder>
//   --doTestRun --doTestRunLines <lines>  --threads 1  --binSize 10000
//   --restrictionSequence ...  --danglingSequence ...  --restrictionCutFile ...
//
// Everything else takes hicBuildMatrix's default: minimum distance 300,
// maximum library insert size 1000, minimum mapping quality 15, duplicate
// check on. --threads 1 is raised to 2 by buildMatrixMethods.py:889 with a
// warning, which the shared code reproduces, so one worker classifies the
// pairs.
//
// ---------------------------------------------------------------------------
// Reproduced
// ---------------------------------------------------------------------------
//  R1 The temporary matrix name is part of the output. It is the first line of
//     QC.log and the first column of every table, so the port creates the same
//     kind of name, tmp plus eight characters plus .h5, in the directory
//     Python's tempfile would choose. The harness compares these files byte
//     for byte after the one named normalisation `quickqc_temporary_matrix`
//     (cpp/scripts/comparators/text.py), which accepts that name and nothing
//     else.
//  R2 `--restrictionCutFile` is appended only `if args.danglingSequence`
//     (hicQuickQC.py:120), a copy and paste of the block above it. Both options
//     are required, so the condition always holds; kept for the record.
//
// ---------------------------------------------------------------------------
// Deliberate deviations
// ---------------------------------------------------------------------------
//  D1 hicQuickQC.py:101-105 builds the hicBuildMatrix command line with
//     str.format and .split(), so a BAM file, QC folder or temporary directory
//     whose path contains whitespace is split into several arguments and the
//     run fails in hicBuildMatrix's argument parser. The port passes the paths
//     through unchanged.
//  D2 The temporary file is removed whether or not the build succeeds. The
//     Python removes it only after a successful QC report (:1366), so a failed
//     run leaves an empty tmp*.h5 behind.
//
// ---------------------------------------------------------------------------
// Not produced
// ---------------------------------------------------------------------------
// As for hicBuildMatrix: the five PNG figures and hicQC.html that
// hicPrepareQCreport renders into the QC folder are matplotlib and a pandas
// Styler page, tier 7 of cpp/PLAN.md. QC.log and the five tables are written.

#include <cstdio>

#include "build_matrix_impl.hpp"
#include "hicx/argparse.hpp"
#include "quick_qc_impl.hpp"

namespace {

const char* const kUsage =
    "usage: hicQuickQC --samFiles two sam files two sam files --QCfolder FOLDER\n"
    "                  --restrictionCutFile BED file [BED file ...]\n"
    "                  --restrictionSequence RESTRICTIONSEQUENCE\n"
    "                  [RESTRICTIONSEQUENCE ...] --danglingSequence DANGLINGSEQUENCE\n"
    "                  [DANGLINGSEQUENCE ...] [--lines LINES] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "The tool hicQuickQC considers the first n lines of two bam/sam files to get a\n"
    "first estimate of the quality of the data. It is highly recommended to set the\n"
    "restriction enzyme and dangling end parameter to get a good quality report.\n"
    "\n"
    "Required arguments:\n"
    "  --samFiles two sam files two sam files, -s two sam files two sam files\n"
    "                        The two PE alignment sam files to process.\n"
    "  --QCfolder FOLDER     Path of folder to save the quality control data of the\n"
    "                        matrix.\n"
    "  --restrictionCutFile BED file [BED file ...], -rs BED file [BED file ...]\n"
    "                        BED file(s) with all restriction cut places.\n"
    "  --restrictionSequence RESTRICTIONSEQUENCE [RESTRICTIONSEQUENCE ...], -seq ...\n"
    "                        Sequence of the restriction site.\n"
    "  --danglingSequence DANGLINGSEQUENCE [DANGLINGSEQUENCE ...]\n"
    "                        Sequence left by the restriction enzyme after cutting.\n"
    "\n"
    "Optional arguments:\n"
    "  --lines LINES         Number of lines to consider for the QC test run (Default:\n"
    "                        1000000).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

// hicQuickQC.py parse_arguments.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicQuickQC",
                       "The tool hicQuickQC considers the first n lines of two bam/sam files to get "
                       "a first estimate of the quality of the data. It is highly recommended to set "
                       "the restriction enzyme and dangling end parameter to get a good quality "
                       "report.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--samFiles", "-s"})
        .metavar("two sam files")
        .nargs(2)
        .required()
        .input({"bam", "sam"})
        .help("The two PE alignment sam files to process.");
    required.add({"--QCfolder"})
        .metavar("FOLDER")
        .required()
        .output({}, "directory")
        .help("Path of folder to save the quality control data of the matrix.");
    // argparse.FileType('r') opens the file while parsing, so a missing cut
    // site file is an argument error with exit status 2.
    required.add({"--restrictionCutFile", "-rs"})
        .file_type("r")
        .metavar("BED file")
        .nargs("+")
        .required()
        .input({"bed"})
        .help("BED file(s) with all restriction cut places.");
    required.add({"--restrictionSequence", "-seq"})
        .type("str")
        .nargs("+")
        .required()
        .help("Sequence of the restriction site.");
    required.add({"--danglingSequence"})
        .type("str")
        .nargs("+")
        .required()
        .help("Sequence left by the restriction enzyme after cutting.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--lines"})
        .default_value(1000000)
        .type("int")
        .help("Number of lines to consider for the QC test run.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.sam_files = ns.strs("samFiles");
    args.qc_folder = ns.str("QCfolder");
    args.restriction_cut_files = ns.strs("restrictionCutFile");
    args.restriction_sequences = ns.strs("restrictionSequence");
    args.dangling_sequences = ns.strs("danglingSequence");
    args.do_test_run_lines = ns.integer("lines");

    // hicQuickQC.py:101-123, the fixed part of the hicBuildMatrix command line.
    args.do_test_run = true;
    args.threads = 1;
    args.bin_size = {10000};
    return args;
}

}  // namespace

int main(int argc, char** argv) {
    g_tool = "hicQuickQC";
    Arguments args = parse_arguments(argc, argv);

    // hicQuickQC.py:94-99 creates the QC folder before anything else.
    std::error_code ec;
    std::filesystem::create_directories(args.qc_folder, ec);

    std::string temporary;
    try {
        temporary = hicx::quick_qc::create_named_temporary_file(
            hicx::quick_qc::python_gettempdir(), ".h5");
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicQuickQC: %s\n", error.what());
        return 1;
    }
    args.out_file_name = temporary;
    const int status = run_build_matrix(args);
    std::filesystem::remove(temporary, ec);  // D2
    return status;
}
