#include "hicx/build_matrix_api.hpp"

#include <algorithm>
#include <cctype>
#include <string>
#include <unordered_map>
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

std::string read_position_key(const std::string& chrom1, const std::string& chrom2, bool chrom1_first,
                              const std::string& start1, const std::string& start2, bool start1_first) {
    std::string key = chrom1_first ? chrom1 + "-" + chrom2 : chrom2 + "-" + chrom1;
    key += start1_first ? "-" + start1 + "-" + start2 : "-" + start2 + "-" + start1;
    return key;
}

std::optional<std::vector<std::int64_t>> python_range(std::int64_t start, std::int64_t stop, std::int64_t step) {
    if (step == 0) {
        return std::nullopt;
    }
    std::vector<std::int64_t> values;
    if (step > 0) {
        for (std::int64_t v = start; v < stop; v += step) {
            values.push_back(v);
            if (v > stop - step) {
                break;
            }
        }
    } else {
        for (std::int64_t v = start; v > stop; v += step) {
            values.push_back(v);
            if (v < stop - step) {
                break;
            }
        }
    }
    return values;
}

std::vector<std::pair<std::size_t, std::size_t>> ordered_dict_items(const std::vector<std::string>& names) {
    std::vector<std::pair<std::size_t, std::size_t>> items;
    std::unordered_map<std::string, std::size_t> position;
    for (std::size_t i = 0; i < names.size(); ++i) {
        const auto [it, inserted] = position.emplace(names[i], items.size());
        if (inserted) {
            items.emplace_back(i, i);
        } else {
            items[it->second].second = i;
        }
    }
    return items;
}

std::string ascii_upper(std::string text) {
    for (char& c : text) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return text;
}

bool starts_with(const std::string& text, const std::string& prefix) {
    return text.size() >= prefix.size() && text.compare(0, prefix.size(), prefix) == 0;
}

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() && text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::size_t supplementary_count(const std::string& sa_tag) {
    return static_cast<std::size_t>(std::count(sa_tag.begin(), sa_tag.end(), ';'));
}

std::int64_t bases_before_match(const std::vector<CigarOperation>& cigar) {
    std::int64_t sum = 0;
    for (const auto& [op, count] : cigar) {
        if (op == 0) {  // pysam.CMATCH
            break;
        }
        sum += count;
    }
    return sum;
}

std::size_t first_index_of_minimum(const std::vector<std::int64_t>& values) {
    return static_cast<std::size_t>(std::min_element(values.begin(), values.end()) - values.begin());
}

std::vector<EnlargeAssignment> enlarge_bins(const std::vector<std::int64_t>& starts,
                                            const std::vector<std::int64_t>& ends,
                                            const std::vector<bool>& same_chrom) {
    using Kind = EnlargeField::Kind;
    const std::size_t n = starts.size();
    struct Current {
        std::size_t chrom_from;
        EnlargeField start;
        EnlargeField end;
    };
    std::vector<Current> current(n);
    for (std::size_t i = 0; i < n; ++i) {
        current[i] = {i, {Kind::OriginalStart, starts[i], i}, {Kind::OriginalEnd, ends[i], i}};
    }
    std::vector<EnlargeAssignment> out;
    const auto assign = [&](std::size_t index, std::size_t chrom_from, const EnlargeField& start,
                            const EnlargeField& end) {
        current[index] = {chrom_from, start, end};
        out.push_back({index, chrom_from, start, end});
    };
    bool chr_start = true;
    for (std::size_t idx = 0; idx + 1 < n; ++idx) {
        const Current here = current[idx];
        const Current next = current[idx + 1];
        EnlargeField start = here.start;
        if (chr_start) {
            start = {Kind::Zero, 0, 0};
            chr_start = false;
        }
        if (same_chrom[idx] && here.end.value != next.start.value) {
            // start_next - int((start_next - end) / 2), truncating toward zero
            const std::int64_t middle = next.start.value - (next.start.value - here.end.value) / 2;
            const EnlargeField field{Kind::Middle, middle, idx + 1};
            assign(idx, here.chrom_from, start, field);
            assign(idx + 1, here.chrom_from, field, next.end);
        }
        if (!same_chrom[idx]) {
            assign(idx, here.chrom_from, start, {Kind::Size, 0, here.chrom_from});
            assign(idx + 1, next.chrom_from, {Kind::Zero, 0, 0}, next.end);
        }
    }
    if (n > 0) {
        const Current last = current[n - 1];
        assign(n - 1, last.chrom_from, last.start, {Kind::Size, 0, last.chrom_from});
    }
    return out;
}

}  // namespace hicx::build_matrix_api
