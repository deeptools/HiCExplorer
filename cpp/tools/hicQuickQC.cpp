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

#include <cerrno>
#include <cstdio>
#include <cstring>

#include "build_matrix_impl.hpp"
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

// argparse.FileType('r') opens the file while parsing, so a missing cut site
// file is an argument error with exit status 2.
void require_readable(const std::string& option, const std::string& path) {
    std::FILE* handle = std::fopen(path.c_str(), "r");
    if (handle == nullptr) {
        fail("argument " + option + ": can't open '" + path + "': [Errno " +
             std::to_string(errno) + "] " + std::strerror(errno) + ": '" + path + "'");
    }
    std::fclose(handle);
}

Arguments parse_arguments(int argc, char** argv) {
    Arguments args;
    bool sam_seen = false;
    bool qc_seen = false;
    bool cut_seen = false;
    bool seq_seen = false;
    bool dangling_seen = false;

    for (int i = 1; i < argc; ++i) {
        std::string name(argv[i]);
        std::optional<std::string> inline_value;
        const std::size_t equals = name.find('=');
        if (equals != std::string::npos && name.rfind("--", 0) == 0) {
            inline_value = name.substr(equals + 1);
            name = name.substr(0, equals);
        }
        auto take_one = [&](const std::string& option) -> std::string {
            if (inline_value.has_value()) {
                return *inline_value;
            }
            if (i + 1 >= argc) {
                fail("argument " + option + ": expected one argument");
            }
            return std::string(argv[++i]);
        };
        auto take_many = [&](const std::string& option) -> std::vector<std::string> {
            std::vector<std::string> values;
            if (inline_value.has_value()) {
                values.push_back(*inline_value);
            }
            while (i + 1 < argc && !is_option(argv[i + 1])) {
                values.emplace_back(argv[++i]);
            }
            if (values.empty()) {
                fail("argument " + option + ": expected at least one argument");
            }
            return values;
        };

        if (name == "-h" || name == "--help") {
            std::fputs(kUsage, stdout);
            std::fputs(kHelp, stdout);
            std::exit(0);
        }
        if (name == "--version") {
            std::printf("hicQuickQC %s\n", hicx::kVersion);
            std::exit(0);
        }
        if (name == "-s" || name == "--samFiles") {
            args.sam_files = take_many("--samFiles/-s");
            if (args.sam_files.size() != 2) {
                fail("argument --samFiles/-s: expected 2 arguments");
            }
            sam_seen = true;
        } else if (name == "--QCfolder") {
            args.qc_folder = take_one("--QCfolder");
            qc_seen = true;
        } else if (name == "-rs" || name == "--restrictionCutFile") {
            args.restriction_cut_files = take_many("--restrictionCutFile/-rs");
            for (const auto& path : args.restriction_cut_files) {
                require_readable("--restrictionCutFile/-rs", path);
            }
            cut_seen = true;
        } else if (name == "-seq" || name == "--restrictionSequence") {
            args.restriction_sequences = take_many("--restrictionSequence/-seq");
            seq_seen = true;
        } else if (name == "--danglingSequence") {
            args.dangling_sequences = take_many("--danglingSequence");
            dangling_seen = true;
        } else if (name == "--lines") {
            args.do_test_run_lines = parse_int(name, take_one(name));
        } else {
            fail("unrecognized arguments: " + name);
        }
    }

    std::string missing;
    auto require = [&missing](bool seen, const char* text) {
        if (!seen) {
            missing += missing.empty() ? text : (std::string(", ") + text);
        }
    };
    require(sam_seen, "--samFiles/-s");
    require(qc_seen, "--QCfolder");
    require(cut_seen, "--restrictionCutFile/-rs");
    require(seq_seen, "--restrictionSequence/-seq");
    require(dangling_seen, "--danglingSequence");
    if (!missing.empty()) {
        fail("the following arguments are required: " + missing);
    }

    // hicQuickQC.py:101-123, the fixed part of the hicBuildMatrix command line.
    args.do_test_run = true;
    args.threads = 1;
    args.bin_size = {10000};
    return args;
}

}  // namespace

int main(int argc, char** argv) {
    g_tool = "hicQuickQC";
    g_usage = kUsage;
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
