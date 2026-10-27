// The tool list of the `hicexplorer` command and of hicexplorer.list_tools:
// which installed commands are HiCExplorer tools, the one line summary of a
// tool's description, and the aligned listing.
//
// Shared by the hicexplorer executable (cpp/tools/hicexplorer.cpp) and by the
// extension (hicexplorer._hicx.list_tools).

#ifndef HICX_TOOL_LIST_HPP
#define HICX_TOOL_LIST_HPP

#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace hicx::toollist {

// A file the installed distribution lists (importlib.metadata PackagePath).
struct DistributionFile {
    std::vector<std::string> parts;  // path components
    std::string name;                // last component
    std::string stem;                // name without the last suffix
    std::string suffix;              // ".exe" or ""
};

// The commands of the distribution, sorted: its console scripts and the files
// of its bin/ or Scripts/ directory (stem for .exe files), without
// `hicexplorer` itself, only names starting with hic or chic.
[[nodiscard]] std::vector<std::string> command_names(
    const std::vector<std::string>& console_scripts, const std::vector<DistributionFile>& files);

// description.strip().split('\n')[0].split('. ')[0].rstrip('.')
[[nodiscard]] std::string summary(std::string_view description);

// '{:<{}}{}\n'.format(name, width, summary) per entry in the given order,
// width the longest name plus one (in code points).
[[nodiscard]] std::string format_tool_list(
    const std::vector<std::pair<std::string, std::string>>& entries);

}  // namespace hicx::toollist

#endif  // HICX_TOOL_LIST_HPP
