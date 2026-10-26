#include "hicx/build_matrix_api.hpp"

#include <string>
#include <vector>

namespace hicx::build_matrix_api {

std::string normpath(const std::string& path) {
    if (path.empty()) {
        return ".";
    }
    int initial_slashes = path[0] == '/' ? 1 : 0;
    if (initial_slashes == 1 && path.size() > 1 && path[1] == '/' && (path.size() < 3 || path[2] != '/')) {
        initial_slashes = 2;
    }
    std::vector<std::string> parts;
    std::size_t begin = 0;
    while (begin <= path.size()) {
        std::size_t end = path.find('/', begin);
        if (end == std::string::npos) {
            end = path.size();
        }
        const std::string part = path.substr(begin, end - begin);
        if (!part.empty() && part != ".") {
            if (part != ".." || (initial_slashes == 0 && parts.empty()) ||
                (!parts.empty() && parts.back() == "..")) {
                parts.push_back(part);
            } else if (!parts.empty()) {
                parts.pop_back();
            }
        }
        begin = end + 1;
    }
    std::string out(static_cast<std::size_t>(initial_slashes), '/');
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (i > 0) {
            out += '/';
        }
        out += parts[i];
    }
    return out.empty() ? "." : out;
}

std::string abspath(const std::string& path, const std::string& cwd) {
    return normpath(path.empty() || path[0] != '/' ? join(cwd, path) : path);
}

std::string dirname(const std::string& path) {
    const std::size_t slash = path.rfind('/');
    std::string head = slash == std::string::npos ? std::string() : path.substr(0, slash + 1);
    if (!head.empty() && head.find_first_not_of('/') != std::string::npos) {
        while (!head.empty() && head.back() == '/') {
            head.pop_back();
        }
    }
    return head;
}

std::string join(const std::string& a, const std::string& b) {
    if (!b.empty() && b[0] == '/') {
        return b;
    }
    if (a.empty() || a.back() == '/') {
        return a + b;
    }
    return a + "/" + b;
}

TesterPaths tester_paths(const std::optional<std::string>& test_data_dir, const std::string& module_file,
                         const std::string& cwd, int levels) {
    TesterPaths paths;
    if (test_data_dir.has_value() && !test_data_dir->empty()) {
        paths.root = *test_data_dir;
    } else {
        std::string directory = abspath(module_file, cwd);
        for (int i = 0; i < levels; ++i) {
            directory = dirname(directory);
        }
        paths.root = directory + "/test/test_data/";
    }
    paths.bam_file_1 = join(paths.root, "hic.bam");
    return paths;
}

}  // namespace hicx::build_matrix_api
