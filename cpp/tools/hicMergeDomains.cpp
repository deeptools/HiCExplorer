// Port of hicexplorer/hicMergeDomains.py.
//
// Merges TAD domain files called at several resolutions into one list, filters
// boundaries against protein peaks, derives parent and child relations between
// overlapping TADs and draws one relation tree per chromosome. The algorithm is
// in merge_domains_impl.cpp; this file is the command line and the graphviz
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

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

#include "hicx/resource_usage.hpp"
#include "hicx/version.hpp"
#include "merge_domains_impl.hpp"

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

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicMergeDomains: error: %s\n", message.c_str());
    std::exit(2);
}

// argparse treats a token as an option when it starts with '-' and is not a
// negative number (the parser defines no option that looks like one).
bool looks_like_option(const std::string& token) {
    if (token.size() < 2 || token[0] != '-') {
        return false;
    }
    std::size_t i = 1;
    bool digits = false;
    while (i < token.size() && std::isdigit(static_cast<unsigned char>(token[i])) != 0) {
        ++i;
        digits = true;
    }
    if (i < token.size() && token[i] == '.') {
        ++i;
        digits = false;
        while (i < token.size() && std::isdigit(static_cast<unsigned char>(token[i])) != 0) {
            ++i;
            digits = true;
        }
    }
    const bool negative_number = digits && i == token.size();
    return !negative_number && token.find(' ') == std::string::npos;
}

Arguments parse_arguments(int argc, char** argv) {
    Arguments args;
    bool domains_seen = false;
    std::vector<std::string> tokens(argv + 1, argv + argc);

    for (std::size_t i = 0; i < tokens.size(); ++i) {
        std::string name = tokens[i];
        std::optional<std::string> inline_value;
        const std::size_t equals = name.find('=');
        if (name.rfind("--", 0) == 0 && equals != std::string::npos) {
            inline_value = name.substr(equals + 1);
            name = name.substr(0, equals);
        }
        if (name == "-h" || name == "--help") {
            std::fputs(kUsage, stdout);
            std::fputs(kHelp, stdout);
            std::exit(0);
        }
        if (name == "--version") {
            std::printf("hicMergeDomains %s\n", hicx::kVersion);
            std::exit(0);
        }
        const auto single = [&](const std::string& display) -> std::string {
            if (inline_value.has_value()) {
                return *inline_value;
            }
            if (i + 1 >= tokens.size() || looks_like_option(tokens[i + 1])) {
                fail("argument " + display + ": expected one argument");
            }
            return tokens[++i];
        };
        const auto as_int = [&](const std::string& display, const std::string& text) {
            try {
                return md::py_int(text);
            } catch (const md::PythonError&) {
                fail("argument " + display + ": invalid int value: '" + text + "'");
            }
        };
        if (name == "-d" || name == "--domainFiles") {
            args.domain_files.clear();
            if (inline_value.has_value()) {
                args.domain_files.push_back(*inline_value);
            }
            while (i + 1 < tokens.size() && !looks_like_option(tokens[i + 1])) {
                args.domain_files.push_back(tokens[++i]);
            }
            if (args.domain_files.empty()) {
                fail("argument --domainFiles/-d: expected at least one argument");
            }
            domains_seen = true;
        } else if (name == "-p" || name == "--proteinFile") {
            args.protein_file = single("--proteinFile/-p");
        } else if (name == "-m" || name == "--minimumNumberOfPeaks") {
            const std::string display = "--minimumNumberOfPeaks/-m";
            args.minimum_number_of_peaks = as_int(display, single(display));
        } else if (name == "-v" || name == "--value") {
            args.value = as_int("--value/-v", single("--value/-v"));
        } else if (name == "-pe" || name == "--percent") {
            const std::string text = single("--percent/-pe");
            try {
                args.percent = md::py_float(text);
            } catch (const md::PythonError&) {
                fail("argument --percent/-pe: invalid float value: '" + text + "'");
            }
        } else if (name == "-om" || name == "--outputMergedList") {
            args.output_merged_list = single("--outputMergedList/-om");
        } else if (name == "-or" || name == "--outputRelationList") {
            args.output_relation_list = single("--outputRelationList/-or");
        } else if (name == "-ot" || name == "--outputTreePlotPrefix") {
            args.output_tree_plot_prefix = single("--outputTreePlotPrefix/-ot");
        } else if (name == "-of" || name == "--outputTreePlotFormat") {
            args.output_tree_plot_format = single("--outputTreePlotFormat/-of");
        } else {
            fail("unrecognized arguments: " + tokens[i]);
        }
    }
    if (!domains_seen) {
        fail("the following arguments are required: --domainFiles/-d");
    }
    return args;
}

// shutil.which('dot') as subprocess resolves it: the first executable regular
// file named dot in a PATH entry.
std::optional<std::string> find_on_path(const std::string& program) {
    const char* path = std::getenv("PATH");
    if (path == nullptr) {
        return std::nullopt;
    }
    std::string entries(path);
    std::size_t start = 0;
    while (start <= entries.size()) {
        std::size_t colon = entries.find(':', start);
        if (colon == std::string::npos) {
            colon = entries.size();
        }
        std::string directory = entries.substr(start, colon - start);
        if (directory.empty()) {
            directory = ".";
        }
        const std::string candidate = directory + "/" + program;
        std::error_code error;
        if (std::filesystem::is_regular_file(candidate, error) &&
            ::access(candidate.c_str(), X_OK) == 0) {
            return candidate;
        }
        start = colon + 1;
    }
    return std::nullopt;
}

std::string lower(const std::string& text) {
    std::string out;
    for (const char c : text) {
        out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    return out;
}

// graphviz.Digraph.render(filename, cleanup=True), graphviz 0.20.3:
// saving.save (creating missing parent directories), backend.rendering.render
// with outfile=None, which runs [dot, -Kdot, -T<format>, -O, <file name>] with
// the working directory set to the file's parent when it has one, then
// os.remove of the source.
void render_graph(const md::TreeGraph& graph, const std::string& format) {
    std::fputs(("Saved relation tree of " + graph.chrom + "\n").c_str(), stdout);
    std::fflush(stdout);

    const std::filesystem::path source(graph.filename);
    if (source.has_parent_path()) {
        std::filesystem::create_directories(source.parent_path());
    }
    {
        std::ofstream out(graph.filename, std::ios::binary | std::ios::trunc);
        if (!out) {
            throw md::PythonError("cannot write the DOT source " + graph.filename);
        }
        out << graph.source;
        if (!out.flush()) {
            throw md::PythonError("cannot write the DOT source " + graph.filename);
        }
    }

    const std::string type_flag = "-T" + lower(format);
    const std::string file_name = source.filename().string();
    const std::string directory = source.has_parent_path() ? source.parent_path().string()
                                                           : std::string();
    std::fflush(stderr);
    const pid_t child = ::fork();
    if (child < 0) {
        throw md::PythonError(std::string("fork failed: ") + std::strerror(errno));
    }
    if (child == 0) {
        if (!directory.empty() && ::chdir(directory.c_str()) != 0) {
            std::fprintf(stderr, "hicMergeDomains: cannot change into %s: %s\n",
                         directory.c_str(), std::strerror(errno));
            ::_exit(127);
        }
        const char* const dot_argv[] = {"dot", "-Kdot", type_flag.c_str(), "-O",
                                        file_name.c_str(), nullptr};
        ::execvp("dot", const_cast<char* const*>(dot_argv));
        std::fprintf(stderr, "hicMergeDomains: cannot execute dot: %s\n",
                     std::strerror(errno));
        ::_exit(127);
    }
    int status = 0;
    while (::waitpid(child, &status, 0) < 0) {
        if (errno != EINTR) {
            throw md::PythonError(std::string("waitpid failed: ") + std::strerror(errno));
        }
    }
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        // graphviz raises CalledProcessError before the cleanup, so the source
        // stays on disk.
        throw md::PythonError(
            "graphviz.backend.execute.CalledProcessError: dot -Kdot " + type_flag + " -O " +
            file_name + " returned non-zero exit status " +
            std::to_string(WIFEXITED(status) ? WEXITSTATUS(status) : -1) + " for " +
            graph.filename);
    }
    std::error_code error;
    std::filesystem::remove(source, error);
    if (error) {
        throw md::PythonError("cannot remove the DOT source " + graph.filename + ": " +
                              error.message());
    }
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
        if (!find_on_path("dot").has_value()) {
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
                                render_graph(graph, args.output_tree_plot_format);
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
