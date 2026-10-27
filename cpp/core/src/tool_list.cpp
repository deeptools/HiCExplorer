#include "hicx/tool_list.hpp"

#include <algorithm>
#include <set>

namespace hicx::toollist {

namespace {

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

std::size_t code_points(const std::string& text) {
    std::size_t n = 0;
    for (const char c : text) {
        n += (static_cast<unsigned char>(c) & 0xC0) != 0x80 ? 1 : 0;
    }
    return n;
}

bool starts_with(const std::string& text, std::string_view prefix) {
    return text.size() >= prefix.size() && text.compare(0, prefix.size(), prefix) == 0;
}

}  // namespace

std::vector<std::string> command_names(const std::vector<std::string>& console_scripts,
                                       const std::vector<DistributionFile>& files) {
    std::set<std::string> names(console_scripts.begin(), console_scripts.end());
    for (const DistributionFile& file : files) {
        if (file.parts.size() > 1) {
            const std::string& directory = file.parts[file.parts.size() - 2];
            if (directory == "bin" || directory == "Scripts") {
                names.insert(file.suffix == ".exe" ? file.stem : file.name);
            }
        }
    }
    names.erase("hicexplorer");
    std::vector<std::string> out;
    for (const std::string& name : names) {
        if (starts_with(name, "hic") || starts_with(name, "chic")) {
            out.push_back(name);
        }
    }
    // std::set orders the UTF-8 bytes, which is the code point order of sorted()
    return out;
}

std::string summary(std::string_view description) {
    std::string_view first = strip(description);
    first = first.substr(0, first.find('\n'));
    first = first.substr(0, first.find(". "));
    while (!first.empty() && first.back() == '.') {
        first.remove_suffix(1);
    }
    return std::string(first);
}

std::string format_tool_list(const std::vector<std::pair<std::string, std::string>>& entries) {
    std::size_t width = 0;
    for (const auto& entry : entries) {
        width = std::max(width, code_points(entry.first) + 1);
    }
    std::string out;
    for (const auto& entry : entries) {
        out += entry.first;
        out.append(width - code_points(entry.first), ' ');
        out += entry.second + "\n";
    }
    return out;
}

}  // namespace hicx::toollist
