// The `hicexplorer` command: prints the version, the citation and the list of
// the HiCExplorer tools, as hicexplorer/list_tools.py of HiCExplorer 3.7.6
// does (argparse with RawDescriptionHelpFormatter, help printed when the
// command is run without arguments).
//
// The tool set is the one this build installs: every executable target of
// cpp/tools and the Python drivers of pyproject.toml's console scripts. The
// list and each tool's description are generated at build time by
// cpp/cmake/HicxToolList.cmake into hicexplorer_tool_list.inc, from each
// executable's own --help-json and from each driver's argparse description,
// so they cannot drift from the tools. The line shown per tool is the first
// sentence of the first line of that description, as the Python v4
// list_tools.py printed it.

#include <algorithm>
#include <cstdio>
#include <string>
#include <string_view>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/version.hpp"

namespace {

struct ToolEntry {
    const char* name;
    const char* description;
};

#include "hicexplorer_tool_list.inc"

const char* const kCitation =
    "HiCExplorer addresses the common tasks of Hi-C analysis from processing to visualization.\n"
    "Each tool should be called by its own name as in the following example:\n"
    "\n"
    " $ hicPlotMatrix -m hic_matrix.h5 -o plot.pdf\n"
    "\n"
    "If you find HiCExplorer useful for your research please cite as:\n"
    "\n"
    "Fidel Ramirez, Vivek Bhardwaj, Jose Villaveces, Laura Arrigoni, Bjoern A Gruening, Kin Chung Lam,\n"
    "Bianca Habermann, Asifa Akhtar, Thomas Manke.\n"
    "\"High-resolution TADs reveal DNA sequences underlying genome organization in flies\".\n"
    "Nature Communications, Volume 9, Article number: 189 (2018), doi: "
    "https://doi.org/10.1038/s41467-017-02525-w\n"
    "\n"
    "Joachim Wolff, Vivek Bhardwaj, Stephan Nothjunge, Gautier Richard, Gina Renschler, Ralf Gilsbach,\n"
    "Thomas Manke, Rolf Backofen, Fidel Ram\xc3\xadrez, Bj\xc3\xb6rn A Gr\xc3\xbcning.\n"
    "\"Galaxy HiCExplorer: a web server for reproducible Hi-C data analysis, quality control and "
    "visualization\",\n"
    "Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W11\xe2\x80\x93W16, doi: "
    "https://doi.org/10.1093/nar/gky504\n";

const char* const kUsage = "usage: hicexplorer [-h] [--version]\n";

// Python's str.strip() for the ASCII whitespace characters.
bool is_space(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\x0b' || c == '\x0c' ||
           (c >= '\x1c' && c <= '\x1f');
}

std::string_view strip(std::string_view text) {
    while (!text.empty() && is_space(text.front())) {
        text.remove_prefix(1);
    }
    while (!text.empty() && is_space(text.back())) {
        text.remove_suffix(1);
    }
    return text;
}

// description.strip().split('\n')[0].split('. ')[0].rstrip('.')
std::string summary(std::string_view description) {
    std::string_view first = strip(description);
    first = first.substr(0, first.find('\n'));
    first = first.substr(0, first.find(". "));
    while (!first.empty() && first.back() == '.') {
        first.remove_suffix(1);
    }
    return std::string(first);
}

// '{:<{}}{}\n'.format(name, width, summary) per tool, names sorted, width
// the longest name plus one.
std::string tool_list() {
    std::vector<ToolEntry> tools(std::begin(kTools), std::end(kTools));
    std::sort(tools.begin(), tools.end(), [](const ToolEntry& a, const ToolEntry& b) {
        return std::string_view(a.name) < std::string_view(b.name);
    });
    std::size_t width = 0;
    for (const ToolEntry& tool : tools) {
        width = std::max(width, std::string_view(tool.name).size() + 1);
    }
    std::string out;
    for (const ToolEntry& tool : tools) {
        std::string name(tool.name);
        name.resize(std::max(width, name.size()), ' ');
        out += name + summary(tool.description) + "\n";
    }
    return out;
}

// What argparse prints after the usage line: the raw description, the
// options section; argparse collapses the description's leading newline and
// separates the sections by one empty line.
std::string help_text() {
    return std::string("\n") + kCitation + "\nThe following is the list of tools:\n\n" + tool_list() +
           "\nFor more information visit: http://hicexplorer.readthedocs.org\n"
           "\n"
           "options:\n"
           "  -h, --help  show this help message and exit\n"
           "  --version   show program's version number and exit\n";
}

}  // namespace

int main(int argc, char** argv) {
    namespace cli = hicx::cli;
    const std::string description = std::string("\n") + kCitation;
    cli::Parser parser("hicexplorer", description);
    parser.set_usage(kUsage).set_help(help_text()).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& options = parser.group("options");
    options.add({"-h", "--help"}).action(cli::Action::Help).help("show this help message and exit");
    options.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    if (argc == 1) {
        // list_tools.main: without arguments the help is printed.
        char help[] = "--help";
        char* help_argv[] = {argv[0], help, nullptr};
        (void)parser.parse(2, help_argv);
        return 0;
    }
    (void)parser.parse(argc, argv);
    return 0;
}
