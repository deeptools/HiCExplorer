// hicMatrixQC: a quality-control report for contact matrices (cool, mcool).
//
// No Python HiCExplorer counterpart. The computation is hicx::matrixqc
// (cpp/core/include/hicx/matrix_qc.hpp): one streamed pass over the stored
// pixels of each matrix, read through the cool layer of the matrix library;
// the report (matrix_qc_report.hpp) is an HTML file in the style of hicQC.html
// that Galaxy's HTML sanitizer leaves readable, and a table of all statistics
// with one column per sample.

#include <cstdio>
#include <exception>
#include <fstream>
#include <string>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/matrix_qc.hpp"
#include "hicx/matrix_qc_report.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/version.hpp"

namespace {

namespace cli = hicx::cli;

const char* const kUsage =
    "usage: hicMatrixQC --matrices MATRICES [MATRICES ...] --outFileName OUTFILENAME\n"
    "                   [--tableFileName TABLEFILENAME] [--labels LABELS [LABELS ...]]\n"
    "                   [--resolutions RESOLUTIONS [RESOLUTIONS ...]]\n"
    "                   [--shortRangeCutoff BP] [--decayRange LOW HIGH]\n"
    "                   [--mapResolutionContacts N] [--mapResolutionFraction F]\n"
    "                   [--madThreshold LOWER UPPER] [--transThreshold X]\n"
    "                   [--maxChromosomes N] [--weightColumn NAME] [--threads INT]\n"
    "                   [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Quality control of contact matrices (cool or mcool). The report is the\n"
    "counterpart of the BAM-based QC of hicBuildMatrix for the finished matrix:\n"
    "\n"
    "  - the counts of QC.log and the distance plots, when hicBuildMatrix stored\n"
    "    them in the metadata of the matrix (key hicexplorer_qc),\n"
    "  - contacts, cis and trans, short and long range, diagonal fractions,\n"
    "  - the contact distance curve P(s) and its decay exponent,\n"
    "  - the map resolution after Rao et al. 2014,\n"
    "  - the coverage per bin and the bins hicCorrectMatrix' MAD filter removes,\n"
    "  - the state of the balancing weights,\n"
    "  - observed over expected trans contacts per chromosome pair,\n"
    "  - X over autosome coverage, chromosome Y and chromosome M.\n"
    "\n"
    "Every matrix is read once, in chunks, so that memory use depends on the\n"
    "number of bins and not on the number of contacts. The report is an HTML\n"
    "file without scripts; the statistics are also written as a tab separated\n"
    "table with one column per sample. An mcool file gives one sample per\n"
    "resolution. The statistics that need the bin size are not computed for a\n"
    "matrix with bins of different sizes.\n"
    "\n"
    "An example usage is:\n"
    "$ hicMatrixQC --matrices matrix.cool --outFileName qc.html\n"
    "\n"
    "Required arguments:\n"
    "  --matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]\n"
    "                        The cool or mcool matrices.\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        File name of the HTML report.\n"
    "\n"
    "Optional arguments:\n"
    "  --tableFileName TABLEFILENAME, -t TABLEFILENAME\n"
    "                        File name of the table of all statistics (tab\n"
    "                        separated). Default: the report's name with the\n"
    "                        extension .tsv.\n"
    "  --labels LABELS [LABELS ...]\n"
    "                        One label per matrix. Default: the file names.\n"
    "  --resolutions RESOLUTIONS [RESOLUTIONS ...], -r RESOLUTIONS [RESOLUTIONS ...]\n"
    "                        The resolutions of an mcool file to analyse.\n"
    "                        Default: all resolutions present. Ignored for a cool\n"
    "                        file.\n"
    "  --shortRangeCutoff BP\n"
    "                        Cis contacts closer than this distance (bp) are short\n"
    "                        range. The distance is the distance between the starts\n"
    "                        of the two bins, so the cut-off is exact up to one bin\n"
    "                        size; contacts inside one bin (the main diagonal) are\n"
    "                        short range. (Default: 20000).\n"
    "  --decayRange LOW HIGH\n"
    "                        The decay exponent is the slope of log10 P(s) against\n"
    "                        log10 s for distances s from LOW to HIGH (bp), fitted\n"
    "                        to ten groups of distances equally spaced on the log\n"
    "                        scale. (Default: 100000 1000000).\n"
    "  --mapResolutionContacts N\n"
    "                        Minimum number of contacts of a bin for the map\n"
    "                        resolution. (Default: 1000).\n"
    "  --mapResolutionFraction F\n"
    "                        Fraction of the bins that must reach that number.\n"
    "                        (Default: 0.8).\n"
    "  --madThreshold LOWER UPPER\n"
    "                        Modified z-scores of the coverage outside this range\n"
    "                        are removed by the MAD filter of hicCorrectMatrix\n"
    "                        (--filterThreshold there). (Default: -1.5 5.0).\n"
    "  --transThreshold X    Chromosome pairs with an observed over expected ratio of\n"
    "                        trans contacts above X are flagged. (Default: 2.0).\n"
    "  --maxChromosomes N    The N longest chromosomes get rows of their own (decay\n"
    "                        exponent, coverage, trans enrichment); the genome-wide\n"
    "                        statistics use all chromosomes. (Default: 30).\n"
    "  --weightColumn NAME   The bin column of the balancing weights. (Default:\n"
    "                        weight).\n"
    "  --threads INT, -p INT\n"
    "                        Number of threads. Each thread holds two vectors with\n"
    "                        one value per bin. (Default: 4).\n"
    "  --help, -h            Show this help message and exit.\n"
    "  --version             Show program's version number and exit.\n";

struct Arguments {
    std::vector<std::string> matrices;
    std::string out_file_name;
    std::string table_file_name;
    std::vector<std::string> labels;
    std::vector<std::int64_t> resolutions;
    hicx::matrixqc::Options options;
};

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicMatrixQC: error: %s\n", message.c_str());
    std::exit(2);
}

Arguments parse_arguments(int argc, char** argv) {
    cli::Parser parser("hicMatrixQC",
                       "Quality control of contact matrices: contacts, cis and trans, the contact "
                       "distance curve, the map resolution, the coverage, the balancing state, "
                       "and the statistics stored by hicBuildMatrix.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrices", "-m"})
        .nargs("+")
        .required()
        .input({"cool", "mcool"})
        .help("The cool or mcool matrices.");
    required.add({"--outFileName", "-o"})
        .required()
        .output({"html"})
        .help("File name of the HTML report.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--tableFileName", "-t"})
        .output({"tabular"})
        .help("File name of the table of all statistics. Default: the report's name with the "
              "extension .tsv.");
    optional.add({"--labels"}).nargs("+").help("One label per matrix. Default: the file names.");
    optional.add({"--resolutions", "-r"})
        .nargs("+")
        .type("int")
        .help("The resolutions of an mcool file to analyse. Default: all resolutions present.");
    optional.add({"--shortRangeCutoff"})
        .type("int")
        .default_value(std::int64_t{20000})
        .help("Cis contacts closer than this distance (bp) are short range.");
    optional.add({"--decayRange"})
        .type("int")
        .nargs(2)
        .help("Distances (bp) LOW HIGH of the decay exponent fit. Default: 100000 1000000.");
    optional.add({"--mapResolutionContacts"})
        .type("float")
        .default_value(1000.0)
        .help("Minimum number of contacts of a bin for the map resolution.");
    optional.add({"--mapResolutionFraction"})
        .type("float")
        .default_value(0.8)
        .help("Fraction of the bins that must reach that number.");
    optional.add({"--madThreshold"})
        .type("float")
        .nargs(2)
        .help("Lower and upper modified z-score of the MAD filter. Default: -1.5 5.0.");
    optional.add({"--transThreshold"})
        .type("float")
        .default_value(2.0)
        .help("Chromosome pairs with an observed over expected ratio of trans contacts above this "
              "value are flagged.");
    optional.add({"--maxChromosomes"})
        .type("int")
        .default_value(std::int64_t{30})
        .help("The longest chromosomes that get rows of their own.");
    optional.add({"--weightColumn"})
        .default_value(hicx::json::Value::string("weight"))
        .help("The bin column of the balancing weights.");
    optional.add({"--threads", "-p"})
        .type("int")
        .metavar("INT")
        .default_value(std::int64_t{4})
        .help("Number of threads.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("Show this help message and exit.");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.matrices = ns.strs("matrices");
    args.out_file_name = ns.str("outFileName");
    if (ns.given("tableFileName")) {
        args.table_file_name = ns.str("tableFileName");
    } else {
        std::string stem = args.out_file_name;
        const std::size_t dot = stem.find_last_of('.');
        const std::size_t slash = stem.find_last_of('/');
        if (dot != std::string::npos && (slash == std::string::npos || dot > slash)) {
            stem.resize(dot);
        }
        args.table_file_name = stem + ".tsv";
    }
    if (ns.given("labels")) {
        args.labels = ns.strs("labels");
        if (args.labels.size() != args.matrices.size()) {
            fail("argument --labels: expected " + std::to_string(args.matrices.size()) +
                 " labels, one per matrix");
        }
    } else {
        args.labels = args.matrices;
    }
    if (ns.given("resolutions")) {
        args.resolutions = ns.integers("resolutions");
    }
    hicx::matrixqc::Options& o = args.options;
    o.short_range_cutoff = ns.integer("shortRangeCutoff");
    if (ns.given("decayRange")) {
        const std::vector<std::int64_t> range = ns.integers("decayRange");
        o.fit_min = range[0];
        o.fit_max = range[1];
    }
    if (!(o.fit_min > 0 && o.fit_min < o.fit_max)) {
        fail("argument --decayRange: LOW and HIGH need 0 < LOW < HIGH");
    }
    o.rao_min_contacts = ns.real("mapResolutionContacts");
    o.rao_fraction = ns.real("mapResolutionFraction");
    if (!(o.rao_fraction > 0 && o.rao_fraction <= 1)) {
        fail("argument --mapResolutionFraction: expected a value in (0, 1]");
    }
    if (ns.given("madThreshold")) {
        const std::vector<double> threshold = ns.reals("madThreshold");
        o.mad_lower = threshold[0];
        o.mad_upper = threshold[1];
    }
    o.trans_threshold = ns.real("transThreshold");
    const std::int64_t chromosomes = ns.integer("maxChromosomes");
    if (chromosomes < 2) {
        fail("argument --maxChromosomes: must be at least 2");
    }
    o.max_chromosomes = static_cast<std::size_t>(chromosomes);
    o.weight_column = ns.str("weightColumn");
    const std::int64_t threads = ns.integer("threads");
    if (threads < 1) {
        fail("argument --threads: must be at least 1");
    }
    o.threads = static_cast<unsigned>(std::min<std::int64_t>(threads, 256));
    return args;
}

void write_file(const std::string& path, const std::string& text) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("cannot write " + path);
    }
    out << text;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    try {
        std::vector<hicx::matrixqc::FileStats> files;
        for (std::size_t i = 0; i < args.matrices.size(); ++i) {
            files.push_back(hicx::matrixqc::analyse_file(args.matrices[i], args.labels[i],
                                                         args.resolutions, args.options));
        }
        write_file(args.out_file_name, hicx::matrixqc::render_html(files, args.options, hicx::kVersion));
        write_file(args.table_file_name, hicx::matrixqc::render_table(files));
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicMatrixQC: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicMatrixQC");
    return 0;
}
