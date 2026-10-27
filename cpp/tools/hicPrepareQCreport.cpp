// Port of hicexplorer/hicPrepareQCreport.py (cpp/PLAN.md tier 7, option (a)),
// built twice: as hicPrepareQCreport and as hicQC, the name bin/hicQC gives
// the same main().
//
// The table aggregation runs here (hicx::write_qc_report_tables, shared with
// hicBuildMatrix and hicQuickQC): the QC logs become pandas' DataFrame, named
// by --labels or by each log's File entry, and the five tab separated tables
// are written as DataFrame.to_csv writes them. hicQC.html (Styler.to_html
// into qc_template.html) and the DataFrames of the five bar charts come from
// the same logs (hicx::write_qc_report_html); the charts are drawn by
// plot/hicexplorer_plot/hicPrepareQCreport.py with the reference's pandas
// plotting calls.
//
// Deliberate deviation: when --labels has as many entries as there are logs
// but the logs yield a different number of rows, the reference logs an error
// and exits 0 without writing anything; the port exits 1
// (cpp/AGENTS_CONTRACT.md rule 7).

#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/build_matrix.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/version.hpp"

#ifndef HICX_QC_PROG
#define HICX_QC_PROG "hicPrepareQCreport"
#endif

namespace {

namespace cli = hicx::cli;
namespace plot = hicx::plot;

const char* const kProg = HICX_QC_PROG;

const char* const kHelp =
    "\n"
    "Tabulates and plots QC measures from hicBuildMatrix log files within an HTML\n"
    "output\n"
    "\n"
    "Required arguments:\n"
    "  --logfiles LOGFILES [LOGFILES ...], -l LOGFILES [LOGFILES ...]\n"
    "                        Path to the log files to be processed\n"
    "  --labels LABELS [LABELS ...]\n"
    "                        Label to assign to each log file. Each label should be\n"
    "                        separated by a space. Quote labels that contain\n"
    "                        spaces: E.g. --labels label1 \"labels 2\"\n"
    "  --outputFolder OUTPUTFOLDER, -o OUTPUTFOLDER\n"
    "                        Several files with be saved under this folder: A table\n"
    "                        containing the results and a html file with several\n"
    "                        images.\n"
    "\n"
    "Optional arguments:\n"
    "  --dpi DPI             Image resolution. By default high resolution png\n"
    "                        images with a 200 dpi are created.\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n"
    "\n"
    "C++ port: the tables and hicQC.html are written in C++, and the figures are\n"
    "drawn by the hicexplorer_plot drawing layer with the calls of the Python tool\n"
    "(HICX_PLOT_PYTHON names the interpreter). The C++-only option --plotData FILE\n"
    "writes the tables, hicQC.html and the data of the figures to FILE instead of\n"
    "drawing them.\n"
    "\n"
    "  --prettyReport        C++-only. Also write OUTPUTFOLDER/hicQC_report.html, a\n"
    "                        self-contained, card-based report, independent of\n"
    "                        hicQC.html, that judges five metrics per sample as\n"
    "                        PASS, WARN or FAIL. Two cutoffs (valid pairs, inter-\n"
    "                        chromosomal contacts) follow HiCExplorer's own usage\n"
    "                        guidance; the rest (duplicated pairs, dangling ends,\n"
    "                        self-ligation) are conservative defaults with no\n"
    "                        published HiCExplorer standard behind them. Every\n"
    "                        cutoff below is a command line option. A metric on\n"
    "                        the right side of its cutoff, or exactly on it, is\n"
    "                        PASS; one past the cutoff by at most 20 % of the\n"
    "                        cutoff is WARN; one further past it is FAIL.\n"
    "  --minValidPairsPercent PCT\n"
    "                        C++-only. Valid (Hi-C contact) pairs, as a percentage\n"
    "                        of sequenced reads. PASS at or above PCT, WARN below\n"
    "                        PCT down to 0.8 * PCT, FAIL below 0.8 * PCT.\n"
    "                        (Default: 25.0).\n"
    "  --maxInterChromosomalPercent PCT\n"
    "                        C++-only. Inter-chromosomal contacts, as a percentage\n"
    "                        of valid pairs. PASS at or below PCT, WARN above PCT\n"
    "                        up to 1.2 * PCT, FAIL above 1.2 * PCT. (Default:\n"
    "                        10.0).\n"
    "  --maxDuplicatedPairsPercent PCT\n"
    "                        C++-only. Duplicated pairs, as a percentage of\n"
    "                        mappable pairs. PASS at or below PCT, WARN above PCT\n"
    "                        up to 1.2 * PCT, FAIL above 1.2 * PCT. (Default:\n"
    "                        20.0).\n"
    "  --maxDanglingEndPercent PCT\n"
    "                        C++-only. Dangling ends, as a percentage of mappable\n"
    "                        pairs. PASS at or below PCT, WARN above PCT up to\n"
    "                        1.2 * PCT, FAIL above 1.2 * PCT. (Default: 10.0).\n"
    "  --maxSelfLigationPercent PCT\n"
    "                        C++-only. Self-ligated pairs, as a percentage of\n"
    "                        mappable pairs. PASS at or below PCT, WARN above PCT\n"
    "                        up to 1.2 * PCT, FAIL above 1.2 * PCT. (Default:\n"
    "                        10.0).\n";

std::string read_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("could not read " + path);
    }
    std::ostringstream text;
    text << in.rdbuf();
    return text.str();
}

}  // namespace

int main(int argc, char** argv) {
    cli::Parser parser(kProg,
                       "Tabulates and plots QC measures from  hicBuildMatrix log files within an "
                       "HTML output");
    const std::string usage =
        std::string("usage: ") + kProg +
        " --logfiles matrix1_QCfolder/QC.log matrix2_QCfolder/QC.log --labels \"sample 1\" "
        "\"sample 2\" --outputFolder QC_all_samples\n";
    parser.set_usage(usage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--logfiles", "-l"})
        .file_type("r")
        .nargs("+")
        .required()
        .input({"log", "txt"})
        .help("Path to the log files to be processed");
    required.add({"--labels"})
        .nargs("+")
        .help("Label to assign to each log file. Each label should be separated by a space. "
              "Quote labels that contain spaces: E.g. --labels label1 \"labels 2\"");
    required.add({"--outputFolder", "-o"})
        .required()
        .output({}, "directory")
        .help("Several files with be saved under this folder: A table containing the results "
              "and a html file with several images.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--dpi"})
        .type("int")
        .default_value(std::int64_t{200})
        .help("Image resolution. By default high resolution png images with a 200 dpi are "
              "created.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    optional.add({"--plotData"})
        .metavar("FILE")
        .output({"json"})
        .cpp_only("The data of the figures as JSON, without drawing them (cpp/PLAN.md tier 7); "
                  "the tables and hicQC.html are written as always.")
        .help("Write the data the figures are drawn from as JSON to this file and do not draw "
              "them.");
    optional.add({"--prettyReport"})
        .action(cli::Action::StoreTrue)
        .cpp_only("A second, card-based HTML report with PASS/WARN/FAIL badges per metric, "
                  "written alongside the unchanged hicQC.html.")
        .help("Also write OUTPUTFOLDER/hicQC_report.html, judging five metrics per sample as "
              "PASS, WARN or FAIL.");
    optional.add({"--minValidPairsPercent"})
        .type("float")
        .default_value(25.0)
        .cpp_only("Only read with --prettyReport.")
        .help("Valid pairs as a percentage of sequenced reads: PASS at or above this value, "
              "WARN below it down to 0.8 times it, FAIL below that.");
    optional.add({"--maxInterChromosomalPercent"})
        .type("float")
        .default_value(10.0)
        .cpp_only("Only read with --prettyReport.")
        .help("Inter-chromosomal contacts as a percentage of valid pairs: PASS at or below this value, "
              "WARN above it up to 1.2 times it, FAIL above that.");
    optional.add({"--maxDuplicatedPairsPercent"})
        .type("float")
        .default_value(20.0)
        .cpp_only("Only read with --prettyReport.")
        .help("Duplicated pairs as a percentage of mappable pairs: PASS at or below this value, "
              "WARN above it up to 1.2 times it, FAIL above that.");
    optional.add({"--maxDanglingEndPercent"})
        .type("float")
        .default_value(10.0)
        .cpp_only("Only read with --prettyReport.")
        .help("Dangling ends as a percentage of mappable pairs: PASS at or below this value, "
              "WARN above it up to 1.2 times it, FAIL above that.");
    optional.add({"--maxSelfLigationPercent"})
        .type("float")
        .default_value(10.0)
        .cpp_only("Only read with --prettyReport.")
        .help("Self-ligated pairs as a percentage of mappable pairs: PASS at or below this value, "
              "WARN above it up to 1.2 times it, FAIL above that.");

    const cli::Namespace ns = parser.parse(argc, argv);
    if (const int refused = hicx::plot::preflight("hicPrepareQCreport", !ns.given("plotData")); refused != 0) {
        return refused;
    }
    const std::string folder = ns.str("outputFolder");
    std::string charts;
    try {
        std::error_code error;
        std::filesystem::create_directories(folder, error);
        if (error && !std::filesystem::is_directory(folder)) {
            throw std::runtime_error("OSError: cannot create " + folder + ": " + error.message());
        }
        std::vector<std::string> logs;
        for (const std::string& path : ns.strs("logfiles")) {
            logs.push_back(read_file(path));
        }
        std::optional<std::vector<std::string>> labels;
        if (ns.given("labels")) {
            labels = ns.strs("labels");
        }
        hicx::write_qc_report_tables(folder, logs, labels);
        charts = hicx::write_qc_report_html(folder, logs, labels);
        if (ns.flag("prettyReport")) {
            hicx::QcThresholds thresholds;
            thresholds.min_valid_pairs_percent = ns.real("minValidPairsPercent");
            thresholds.max_inter_chromosomal_percent = ns.real("maxInterChromosomalPercent");
            thresholds.max_duplicated_pairs_percent = ns.real("maxDuplicatedPairsPercent");
            thresholds.max_dangling_end_percent = ns.real("maxDanglingEndPercent");
            thresholds.max_self_ligation_percent = ns.real("maxSelfLigationPercent");
            hicx::write_qc_pretty_report(folder, logs, labels, thresholds);
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s: %s\n", kProg, error.what());
        return 1;
    }

    plot::JsonObject data;
    data.add("outputFolder", plot::json_string(folder));
    data.add("dpi", plot::json_int(ns.integer("dpi")));
    data.add("charts", charts);
    hicx::report_resource_usage(kProg);
    return plot::draw("hicPrepareQCreport", data.str(), ns.opt_str("plotData"));
}
