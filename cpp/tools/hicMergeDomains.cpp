// Port of hicexplorer/hicMergeDomains.py.
//
// Merges TAD domain files called at several resolutions into one list, filters
// boundaries against protein peaks, derives parent and child relations between
// overlapping TADs and draws one relation tree per chromosome. The algorithm is
// in cpp/core/src/merge_domains.cpp; this file is the command line and the graphviz
// rendering.
//
// Two facts about the reference that decide what "equivalent" means here:
//
//  1. **No clustering.** hicMergeDomains.py:2 imports scipy's linkage and
//     dendrogram and never calls either.
//  2. **The relation tree is rendered by graphviz's `dot` binary, not drawn
//     by the tool.** create_tree builds a graphviz.Digraph and calls
//     render(cleanup=True), which writes the DOT source to <prefix>_<chrom>,
//     runs `dot -Kdot -T<format> -O <name>` in that file's directory, and
//     deletes the source again, leaving <prefix>_<chrom>.<format>. The port
//     writes a DOT source byte identical to Digraph.source and runs the same
//     command, so the rendered file is equivalent by construction: the same
//     input to the same program. Decided by the orchestrating session
//     2026-09-13.
//
// Deliberate deviations, all on failure paths:
//
//  * With more than one domain file the tree is a required output, so an
//    unknown --outputTreePlotFormat and a missing `dot` are refused before any
//    file is written. The reference writes the merged list and the relation
//    list first and only then fails with ValueError or ExecutableNotFound.
//  * A domain file whose first 20 start coordinates include one made only of
//    zeros gives a bin size of 0, and merge_protein then loops forever
//    (hicMergeDomains.py:359). The port exits 1 with a message instead.
//
// Threading: none. The Python takes 1.0 to 1.9 s of CPU on the designated
// inputs, almost all of it interpreter start up and list membership tests;
// the port's cost is reading the files.

#include <cstdio>
#include <fstream>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/version.hpp"
#include "hicx/merge_domains.hpp"

namespace {

namespace md = hicx::merge_domains;

const char* const kUsage =
    "usage: hicMergeDomains --domainFiles DOMAINFILES [DOMAINFILES ...]\n"
    "                       [--proteinFile PROTEINFILE]\n"
    "                       [--minimumNumberOfPeaks MINIMUMNUMBEROFPEAKS]\n"
    "                       [--value VALUE] [--percent PERCENT]\n"
    "                       [--outputMergedList OUTPUTMERGEDLIST]\n"
    "                       [--outputRelationList OUTPUTRELATIONLIST]\n"
    "                       [--outputTreePlotPrefix OUTPUTTREEPLOTPREFIX]\n"
    "                       [--outputTreePlotFormat OUTPUTTREEPLOTFORMAT] [--help]\n"
    "                       [--version]\n";

const char* const kHelp =
    "\n"
    "hicMergeDomains takes as input multiple TAD domain files from hicFindTads. It "
    "merges TADs from different resolutions to one TAD domains file,\n"
    "considers protein peaks from known TAD binding sites and computes a dependency "
    "graph of the TADs.\n"
    "\n"
    "Two TADs are considered as one if they don't overlap at x bins given by "
    "`--value`; TAD borders need to match the protein peaks given by "
    "`--proteinFile`;\n"
    "a relation between two TADs is given by their overlap of area in percent, "
    "parameter `--percent`. The protein peaks are only considered if in one bin at "
    "least `--minPeak`.\n"
    "\n"
    "An example usage is:\n"
    "\n"
    "`$ hicMergeDomains --domainFiles 10kbtad_domains.bed 50kbtad_domains.bed "
    "--proteinFile ctcf_sorted.bed --outputMergedList two_files_ctcf "
    "--outputRelationList two_files_relation_ctcf --outputTreePlotPrefix "
    "two_files_plot_ctcf --outputTreePlotFormat pdf`\n"
    "\n"
    "Required arguments:\n"
    "  --domainFiles DOMAINFILES [DOMAINFILES ...], -d DOMAINFILES [DOMAINFILES ...]\n"
    "                        The domain files of the different resolutions is\n"
    "                        required\n"
    "\n"
    "Optional arguments:\n"
    "  --proteinFile PROTEINFILE, -p PROTEINFILE\n"
    "                        In order to be able to better assess the relationship\n"
    "                        between TADs, the associated protein file (e.g. CTCF\n"
    "                        for mammals) can be included. The protein file is\n"
    "                        required in broadpeak format\n"
    "  --minimumNumberOfPeaks MINIMUMNUMBEROFPEAKS, -m MINIMUMNUMBEROFPEAKS\n"
    "                        Optional parameter to adjust the number of protein\n"
    "                        peaks when adapting the resolution to the domain\n"
    "                        files. At least minimumNumberOfPeaks of unique peaks\n"
    "                        must be in a bin to considered. Otherwise the bin is\n"
    "                        treated like it has no peaks (Default: 1).\n"
    "  --value VALUE, -v VALUE\n"
    "                        Determine a value by how much the boundaries of two\n"
    "                        TADs must at least differ to consider them as two\n"
    "                        separate TADs (Default: 5000).\n"
    "  --percent PERCENT, -pe PERCENT\n"
    "                        For the relationship determination, a percentage is\n"
    "                        required from which area coverage the TADs are related\n"
    "                        to each other.For example, a relationship should be\n"
    "                        entered from 5 percent area coverage -p 0.05 (Default:\n"
    "                        0.5).\n"
    "  --outputMergedList OUTPUTMERGEDLIST, -om OUTPUTMERGEDLIST\n"
    "                        File name for the merged domains list (Default:\n"
    "                        mergedDomains.bed).\n"
    "  --outputRelationList OUTPUTRELATIONLIST, -or OUTPUTRELATIONLIST\n"
    "                        File name for the relationship list of the TADs\n"
    "                        (Default: relationList.txt).\n"
    "  --outputTreePlotPrefix OUTPUTTREEPLOTPREFIX, -ot OUTPUTTREEPLOTPREFIX\n"
    "                        File name prefix for the relationship tree of the TADs\n"
    "                        (Default: relationship_tree_).\n"
    "  --outputTreePlotFormat OUTPUTTREEPLOTFORMAT, -of OUTPUTTREEPLOTFORMAT\n"
    "                        File format of the relationship tree. Supported\n"
    "                        formats are listed on:\n"
    "                        https://www.graphviz.org/doc/info/output.html\n"
    "                        (Default: pdf).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n"
    "\n"
    "C++ port: the relation tree is rendered by running graphviz's `dot`, which\n"
    "must be on PATH whenever more than one domain file is given.\n";

struct Arguments {
    std::vector<std::string> domain_files;
    std::optional<std::string> protein_file;
    std::int64_t minimum_number_of_peaks = 1;
    std::int64_t value = 5000;
    double percent = 0.5;
    std::string output_merged_list = "mergedDomains.bed";
    std::string output_relation_list = "relationList.txt";
    std::string output_tree_plot_prefix = "relationship_tree_";
    std::string output_tree_plot_format = "pdf";
};

// hicMergeDomains.py parse_arguments. The Python parser uses
// conflict_handler='resolve' and redefines --help/-h, which removes argparse's
// own help action, so there is exactly one.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicMergeDomains",
                       "hicMergeDomains takes as input multiple TAD domain files from hicFindTads. "
                       "It merges TADs from different resolutions to one TAD domains file, considers "
                       "protein peaks from known TAD binding sites and computes a dependency graph "
                       "of the TADs.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--domainFiles", "-d"})
        .nargs("+")
        .required()
        .input({"bed"})
        .help("The domain files of the different resolutions is required");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--proteinFile", "-p"})
        .input({"bed"})
        .help("The associated protein file (e.g. CTCF for mammals), in broadpeak format.");
    optional.add({"--minimumNumberOfPeaks", "-m"})
        .type("int")
        .default_value(1)
        .help("At least minimumNumberOfPeaks unique peaks must be in a bin for it to be "
              "considered.");
    optional.add({"--value", "-v"})
        .type("int")
        .default_value(5000)
        .help("How much the boundaries of two TADs must at least differ to consider them as two "
              "separate TADs.");
    optional.add({"--percent", "-pe"})
        .type("float")
        .default_value(0.5)
        .help("The area coverage from which on two TADs are related to each other.");
    optional.add({"--outputMergedList", "-om"})
        .default_value("mergedDomains.bed")
        .output({"bed"})
        .help("File name for the merged domains list");
    optional.add({"--outputRelationList", "-or"})
        .default_value("relationList.txt")
        .output({"txt"})
        .help("File name for the relationship list of the TADs");
    optional.add({"--outputTreePlotPrefix", "-ot"})
        .default_value("relationship_tree_")
        .output({"pdf", "png", "svg"}, "prefix")
        .note("The C++ port renders the tree with graphviz's dot, which must be on PATH when "
              "more than one domain file is given.")
        .help("File name prefix for the relationship tree of the TADs");
    optional.add({"--outputTreePlotFormat", "-of"})
        .default_value("pdf")
        .help("File format of the relationship tree, one of the graphviz output formats.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.domain_files = ns.strs("domainFiles");
    args.protein_file = ns.opt_str("proteinFile");
    args.minimum_number_of_peaks = ns.integer("minimumNumberOfPeaks");
    args.value = ns.integer("value");
    args.percent = ns.real("percent");
    args.output_merged_list = ns.str("outputMergedList");
    args.output_relation_list = ns.str("outputRelationList");
    args.output_tree_plot_prefix = ns.str("outputTreePlotPrefix");
    args.output_tree_plot_format = ns.str("outputTreePlotFormat");
    return args;
}

void write_file(const std::string& path, const std::function<void(std::ostream&)>& body) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw md::PythonError("FileNotFoundError: [Errno 2] No such file or directory: '" +
                              path + "'");
    }
    body(out);
    if (!out.flush()) {
        throw md::PythonError("cannot write " + path);
    }
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);

    if (args.domain_files.size() == 1 && !args.protein_file.has_value()) {
        std::fputs("ERROR:hicexplorer.hicMergeDomains:Please use multiple or domain files "
                   "or at least one domain file and one protein file.\n",
                   stderr);
        return 1;
    }
    const bool with_tree = args.domain_files.size() > 1;
    if (with_tree) {
        if (!md::is_graphviz_format(args.output_tree_plot_format)) {
            std::fprintf(stderr,
                         "hicMergeDomains: unknown --outputTreePlotFormat '%s'; graphviz "
                         "accepts the formats listed on "
                         "https://www.graphviz.org/doc/info/output.html. Nothing was "
                         "written.\n",
                         args.output_tree_plot_format.c_str());
            return 1;
        }
        if (!md::find_on_path("dot").has_value()) {
            std::fputs("hicMergeDomains: graphviz's dot executable was not found on PATH. "
                       "It renders the relation tree, which is written whenever more than "
                       "one domain file is given. Nothing was written.\n",
                       stderr);
            return 1;
        }
    }

    try {
        std::optional<md::ProteinList> proteins;
        if (args.protein_file.has_value()) {
            proteins = md::read_protein(*args.protein_file);
        }
        md::MergedDomains domains = md::merge_domain_files(
            args.domain_files, proteins.has_value() ? &*proteins : nullptr,
            args.minimum_number_of_peaks, args.value);
        md::RowPool& pool = domains.pool;
        const std::vector<std::size_t>& merged = domains.merged;
        write_file(args.output_merged_list,
                   [&](std::ostream& out) { md::write_domain_list(out, pool, merged); });
        if (with_tree) {
            const std::vector<md::Relation> relations =
                md::create_relationship_list(pool, merged, args.percent);
            write_file(args.output_relation_list,
                       [&](std::ostream& out) { md::write_relation_list(out, relations); });
            md::create_tree(relations, pool, merged, args.output_tree_plot_prefix,
                            [&](const md::TreeGraph& graph) {
                                std::fputs(("Saved relation tree of " + graph.chrom + "\n").c_str(),
                                           stdout);
                                std::fflush(stdout);
                                md::render_graph(graph, args.output_tree_plot_format);
                            });
        }
    } catch (const md::ReferenceNeverTerminates& error) {
        std::fprintf(stderr, "hicMergeDomains: %s\n", error.what());
        return 1;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicMergeDomains: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicMergeDomains");
    return 0;
}
