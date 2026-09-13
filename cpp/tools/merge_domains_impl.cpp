// See merge_domains_impl.hpp. Line references are to
// hicexplorer/hicMergeDomains.py.

#include "merge_domains_impl.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <limits>
#include <set>
#include <sstream>
#include <unordered_set>
#include <utility>

namespace hicx::merge_domains {

namespace {

const char* const kIndexError = "IndexError: list index out of range";

template <class T>
const T& at(const std::vector<T>& items, std::size_t index) {
    if (index >= items.size()) {
        throw PythonError(kIndexError);
    }
    return items[index];
}

template <class T>
T& at(std::vector<T>& items, std::size_t index) {
    if (index >= items.size()) {
        throw PythonError(kIndexError);
    }
    return items[index];
}

bool is_py_space(char c) {
    const auto u = static_cast<unsigned char>(c);
    return u == ' ' || (u >= 0x09 && u <= 0x0d) || (u >= 0x1c && u <= 0x1f);
}

std::string py_strip(const std::string& text) {
    std::size_t begin = 0;
    std::size_t end = text.size();
    while (begin < end && is_py_space(text[begin])) {
        ++begin;
    }
    while (end > begin && is_py_space(text[end - 1])) {
        --end;
    }
    return text.substr(begin, end - begin);
}

// repr() of a str, for the ValueError messages.
std::string py_str_repr(const std::string& text) {
    const char quote = text.find('\'') != std::string::npos &&
                               text.find('"') == std::string::npos
                           ? '"'
                           : '\'';
    std::string out(1, quote);
    for (const char c : text) {
        if (c == quote || c == '\\') {
            out.push_back('\\');
        }
        out.push_back(c);
    }
    out.push_back(quote);
    return out;
}

// Removes single underscores between digits, as int() and float() accept
// them. Returns false for any other placement.
bool strip_digit_underscores(const std::string& text, std::string* out) {
    out->clear();
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] == '_') {
            const bool digit_before =
                i > 0 && std::isdigit(static_cast<unsigned char>(text[i - 1])) != 0;
            const bool digit_after =
                i + 1 < text.size() &&
                std::isdigit(static_cast<unsigned char>(text[i + 1])) != 0;
            if (!digit_before || !digit_after) {
                return false;
            }
            continue;
        }
        out->push_back(text[i]);
    }
    return true;
}

// Python's s[start:] for a non negative start.
std::string slice_from(const std::string& text, std::size_t start) {
    return start >= text.size() ? std::string() : text.substr(start);
}

// Python's s[len(s) - 1:], the last character or the empty string.
std::string last_char(const std::string& text) {
    return text.empty() ? std::string() : text.substr(text.size() - 1);
}

// int(a) > float(b) and friends compare exactly in Python, without rounding
// the integer to a double. long double holds every int64 exactly on x86-64.
bool int_greater_than_float(std::int64_t a, double b) {
    if (std::isnan(b)) {
        return false;
    }
    return static_cast<long double>(a) > static_cast<long double>(b);
}

std::string join_row(const Row& row) {
    std::string joined;
    for (std::size_t i = 0; i < row.size(); ++i) {
        if (i > 0) {
            joined.push_back('\t');
        }
        joined.append(row[i]);
    }
    return joined;
}

}  // namespace

// ---------------------------------------------------------------------------

std::string py_rstrip(const std::string& text) {
    std::size_t end = text.size();
    while (end > 0 && is_py_space(text[end - 1])) {
        --end;
    }
    return text.substr(0, end);
}

std::vector<std::string> split_tab(const std::string& text) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (true) {
        const std::size_t tab = text.find('\t', start);
        if (tab == std::string::npos) {
            fields.push_back(text.substr(start));
            return fields;
        }
        fields.push_back(text.substr(start, tab - start));
        start = tab + 1;
    }
}

std::vector<std::string> read_lines(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw PythonError("FileNotFoundError: [Errno 2] No such file or directory: " +
                          py_str_repr(path));
    }
    const std::string content((std::istreambuf_iterator<char>(in)),
                              std::istreambuf_iterator<char>());
    std::vector<std::string> lines;
    std::string current;
    bool pending = false;
    for (std::size_t i = 0; i < content.size(); ++i) {
        const char c = content[i];
        if (c == '\n' || c == '\r') {
            lines.push_back(py_rstrip(current));
            current.clear();
            pending = false;
            if (c == '\r' && i + 1 < content.size() && content[i + 1] == '\n') {
                ++i;
            }
            continue;
        }
        current.push_back(c);
        pending = true;
    }
    if (pending) {
        lines.push_back(py_rstrip(current));
    }
    return lines;
}

std::int64_t py_int(const std::string& text) {
    const std::string stripped = py_strip(text);
    std::string body;
    std::size_t start = 0;
    bool negative = false;
    if (!stripped.empty() && (stripped[0] == '+' || stripped[0] == '-')) {
        negative = stripped[0] == '-';
        start = 1;
    }
    const std::string digits_part = stripped.substr(start);
    bool valid = !digits_part.empty() && strip_digit_underscores(digits_part, &body) &&
                 !body.empty() &&
                 std::all_of(body.begin(), body.end(), [](char c) {
                     return std::isdigit(static_cast<unsigned char>(c)) != 0;
                 });
    if (!valid) {
        throw PythonError("ValueError: invalid literal for int() with base 10: " +
                          py_str_repr(text));
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long long magnitude = std::strtoull(body.c_str(), &end, 10);
    const unsigned long long limit =
        negative ? static_cast<unsigned long long>(std::numeric_limits<std::int64_t>::max()) + 1ULL
                 : static_cast<unsigned long long>(std::numeric_limits<std::int64_t>::max());
    if (errno == ERANGE || magnitude > limit) {
        throw PythonError("an integer field does not fit into 64 bits: " +
                          py_str_repr(text));
    }
    if (negative) {
        return magnitude == limit ? std::numeric_limits<std::int64_t>::min()
                                  : -static_cast<std::int64_t>(magnitude);
    }
    return static_cast<std::int64_t>(magnitude);
}

double py_float(const std::string& text) {
    const std::string stripped = py_strip(text);
    const auto fail = [&text]() -> double {
        throw PythonError("ValueError: could not convert string to float: " +
                          py_str_repr(text));
    };
    std::string lowered;
    for (const char c : stripped) {
        lowered.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    std::string unsigned_part = lowered;
    bool negative = false;
    if (!unsigned_part.empty() && (unsigned_part[0] == '+' || unsigned_part[0] == '-')) {
        negative = unsigned_part[0] == '-';
        unsigned_part.erase(0, 1);
    }
    if (unsigned_part == "inf" || unsigned_part == "infinity") {
        return negative ? -std::numeric_limits<double>::infinity()
                        : std::numeric_limits<double>::infinity();
    }
    if (unsigned_part == "nan") {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::string body;
    if (!strip_digit_underscores(stripped, &body) || body.empty()) {
        return fail();
    }
    // [sign] (digits [. [digits]] | . digits) [e [sign] digits], nothing else:
    // strtod would also take hexadecimal floats, which float() rejects.
    std::size_t i = 0;
    if (body[i] == '+' || body[i] == '-') {
        ++i;
    }
    std::size_t mantissa_digits = 0;
    while (i < body.size() && std::isdigit(static_cast<unsigned char>(body[i])) != 0) {
        ++i;
        ++mantissa_digits;
    }
    if (i < body.size() && body[i] == '.') {
        ++i;
        while (i < body.size() && std::isdigit(static_cast<unsigned char>(body[i])) != 0) {
            ++i;
            ++mantissa_digits;
        }
    }
    if (mantissa_digits == 0) {
        return fail();
    }
    if (i < body.size() && (body[i] == 'e' || body[i] == 'E')) {
        ++i;
        if (i < body.size() && (body[i] == '+' || body[i] == '-')) {
            ++i;
        }
        std::size_t exponent_digits = 0;
        while (i < body.size() && std::isdigit(static_cast<unsigned char>(body[i])) != 0) {
            ++i;
            ++exponent_digits;
        }
        if (exponent_digits == 0) {
            return fail();
        }
    }
    if (i != body.size()) {
        return fail();
    }
    return std::strtod(body.c_str(), nullptr);
}

// ---------------------------------------------------------------------------

DomainList create_list_of_file(RowPool& pool, const std::string& path) {
    DomainList result;
    const std::vector<std::string> lines = read_lines(path);
    result.rows.reserve(lines.size());
    std::int64_t bin_size = 10000000;
    int position = 0;
    for (const std::string& line : lines) {
        pool.rows.push_back(split_tab(line));
        result.rows.push_back(pool.rows.size() - 1);
        const Row& x = pool.rows.back();
        if (position < 20) {
            ++position;
            int zeros = 0;
            std::string num = at(x, 1);
            while (last_char(num) == "0") {
                ++zeros;
                num = num.substr(0, num.size() - 1);
            }
            num = last_char(num);
            num.append(static_cast<std::size_t>(zeros), '0');
            const std::int64_t value = py_int(num);
            if (value < bin_size) {
                bin_size = value;
            }
        }
    }
    result.bin_size = bin_size;
    return result;
}

ProteinList read_protein(const std::string& path) {
    const std::vector<std::string> lines = read_lines(path);
    ProteinList grouped(1);
    std::string actual_chr = split_tab(at(lines, 0))[0];
    for (const std::string& line : lines) {
        const Row x = split_tab(line);
        Row first_three(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(
                                                   std::min<std::size_t>(3, x.size())));
        if (x[0] == actual_chr) {
            grouped.back().push_back(std::move(first_three));
        } else {
            actual_chr = x[0];
            grouped.push_back({std::move(first_three)});
        }
    }
    return grouped;
}

MergedProtein merge_protein(const ProteinList& proteins, std::int64_t bin_size,
                            std::int64_t min_peak) {
    MergedProtein merged;
    for (const auto& chromosome : proteins) {
        merged.emplace_back();
        std::int64_t left = 0;
        std::int64_t right = bin_size;
        std::int64_t count = 0;
        for (const Row& peak : chromosome) {
            const std::int64_t position = py_int(at(peak, 1));
            if (position <= right) {
                ++count;
                continue;
            }
            if (count >= min_peak) {
                merged.back().push_back(ProteinBin{at(peak, 0), left, right, count});
            }
            left = right;
            right = left + bin_size;
            count = 0;
            if (position < right) {
                count += 1;
            } else {
                // :359-360 advances by bin_size until it passes the peak, so a
                // bin size of 0 never gets there.
                if (bin_size <= 0 && position > left) {
                    throw ReferenceNeverTerminates(
                        "the protein peaks are binned with a bin size of " +
                        std::to_string(bin_size) +
                        ", derived from a start coordinate made only of zeros in the "
                        "first 20 lines of a domain file; hicMergeDomains.py:359 loops "
                        "forever on this input");
                }
                while (position > left) {
                    left += bin_size;
                }
                right = left + bin_size;
                count = 1;
            }
        }
    }
    return merged;
}

void compare_boundaries_protein(const RowPool& pool, std::vector<std::size_t>& b_list,
                                const MergedProtein& c_list, double para_score) {
    std::size_t pos_tad = 0;
    std::size_t pos_peak = 0;
    std::size_t chrom_position = 0;
    const auto row = [&](std::size_t index) -> const Row& {
        return pool.rows[at(b_list, index)];
    };
    while (pos_tad < b_list.size()) {
        if (chrom_position < c_list.size() &&
            at(row(pos_tad), 0) != at(c_list[chrom_position], 0).chrom) {
            chrom_position = 0;
            while (pos_tad < b_list.size() && chrom_position < c_list.size() &&
                   at(row(pos_tad), 0) != at(c_list[chrom_position], 0).chrom) {
                ++chrom_position;
            }
            pos_peak = 0;
        } else if (chrom_position < c_list.size()) {
            const std::vector<ProteinBin>& peaks = c_list[chrom_position];
            while (pos_peak + 1 < peaks.size() &&
                   py_int(at(row(pos_tad), 1)) > at(peaks, pos_peak).right) {
                ++pos_peak;
            }
            if (py_int(at(row(pos_tad), 1)) < at(peaks, pos_peak).left) {
                if (pos_tad != 0) {
                    const double previous = py_float(at(row(pos_tad - 1), 4));
                    const double current = py_float(at(row(pos_tad), 4));
                    if (std::abs(previous - current) < para_score &&
                        pos_tad != b_list.size() - 1) {
                        const double next = py_float(at(row(pos_tad + 1), 4));
                        if (std::abs(next - py_float(at(row(pos_tad), 4))) < para_score) {
                            // bList.remove(bList[posTad]) removes the first row
                            // with equal content, which need not be this one.
                            const Row& target = row(pos_tad);
                            for (std::size_t j = 0; j < b_list.size(); ++j) {
                                if (pool.rows[b_list[j]] == target) {
                                    b_list.erase(b_list.begin() +
                                                 static_cast<std::ptrdiff_t>(j));
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        ++pos_tad;
    }
}

std::vector<std::size_t> merge_list(const RowPool& pool, const std::vector<std::size_t>& d1,
                                    const std::vector<std::size_t>& d2,
                                    std::int64_t p_value) {
    const auto r1 = [&](std::size_t index) -> const Row& { return pool.rows[at(d1, index)]; };
    const auto r2 = [&](std::size_t index) -> const Row& { return pool.rows[at(d2, index)]; };
    std::size_t pos1 = 0;
    std::size_t pos2 = 0;
    std::vector<std::size_t> merged;
    // Content of every merged row, for the `not in merged_list` test at :179.
    std::unordered_set<std::string> merged_content;
    const auto append = [&](std::size_t row_index) {
        merged.push_back(row_index);
        merged_content.insert(join_row(pool.rows[row_index]));
    };
    const auto far_apart = [&](const Row& a, const Row& b) {
        const std::int64_t left_a = py_int(at(a, 1));
        const std::int64_t left_b = py_int(at(b, 1));
        if (std::llabs(left_a - left_b) > p_value) {
            return true;
        }
        return std::llabs(py_int(at(a, 2)) - py_int(at(b, 2))) > p_value;
    };

    while (true) {
        if (pos1 == d1.size()) {
            break;
        }
        while (at(r1(pos1), 0) == at(r2(pos2), 0)) {
            if (py_int(at(r1(pos1), 1)) <= py_int(at(r2(pos2), 1))) {
                if (far_apart(r1(pos1), r2(pos2))) {
                    append(d1[pos1]);
                }
                if (pos1 + 1 != d1.size()) {
                    ++pos1;
                } else {
                    while (pos2 < d2.size() && at(r1(pos1), 0) == at(r2(pos2), 0)) {
                        append(d2[pos2]);
                        ++pos2;
                    }
                    break;
                }
            } else if (py_int(at(r1(pos1), 1)) > py_int(at(r2(pos2), 1))) {
                if (far_apart(r1(pos1), r2(pos2))) {
                    append(d2[pos2]);
                }
                if (pos2 + 1 != d2.size()) {
                    ++pos2;
                } else {
                    while (pos1 < d1.size() && at(r1(pos1), 0) == at(r2(pos2), 0)) {
                        append(d1[pos1]);
                        ++pos1;
                    }
                    break;
                }
            }
        }
        if (pos1 == d1.size()) {
            break;
        }
        const std::size_t old_pos2 = pos2;
        while (pos2 < d2.size()) {
            if (at(r1(pos1), 0) != at(r2(pos2), 0)) {
                ++pos2;
            } else {
                break;
            }
        }
        if (pos2 == d2.size()) {
            pos2 = old_pos2;
            const std::string chrom = at(r1(pos1), 0);
            while (at(r1(pos1), 0) == chrom) {
                // The same row object can be appended here a second time
                // (:171-173 after :139-142 consumed the last row of d1).
                append(d1[pos1]);
                ++pos1;
                if (pos1 == d1.size()) {
                    break;
                }
            }
        }
    }
    for (std::size_t index = 0; index < d2.size(); ++index) {
        std::string content = join_row(pool.rows[d2[index]]);
        if (merged_content.find(content) == merged_content.end()) {
            merged.push_back(d2[index]);
            merged_content.insert(std::move(content));
        }
    }
    return merged;
}

void add_id(RowPool& pool, const std::vector<std::size_t>& list) {
    std::int64_t id_number = 1;
    for (const std::size_t index : list) {
        at(pool.rows[index], 3) = "ID_" + std::to_string(id_number);
        ++id_number;
    }
}

namespace {

void add_relation_to_list(std::vector<Relation>& relations, const std::string& chrom,
                          const std::string& parent, const std::string& child) {
    if (relations.empty()) {
        relations.push_back(Relation{chrom, parent, {child}});
        return;
    }
    // pos runs below zero and then indexes from the end, as a negative Python
    // index does, until it passes -len and raises IndexError.
    const auto size = static_cast<std::int64_t>(relations.size());
    std::int64_t pos = size - 1;
    while (true) {
        const std::int64_t index = pos >= 0 ? pos : size + pos;
        if (index < 0) {
            throw PythonError(kIndexError);
        }
        Relation& entry = relations[static_cast<std::size_t>(index)];
        if (entry.parent == parent) {
            entry.children.push_back(child);
            return;
        }
        if (py_int(slice_from(entry.parent, 3)) < py_int(slice_from(parent, 3))) {
            relations.push_back(Relation{chrom, parent, {child}});
            return;
        }
        --pos;
    }
}

}  // namespace

std::vector<Relation> create_relationship_list(const RowPool& pool,
                                               const std::vector<std::size_t>& list,
                                               double percent) {
    std::vector<Relation> relations;
    const auto row = [&](std::size_t index) -> const Row& { return pool.rows[list[index]]; };
    std::size_t tad1 = 0;
    std::size_t tad2 = 1;
    while (tad1 < list.size()) {
        while (tad2 < list.size()) {
            const Row& t1 = row(tad1);
            const Row& t2 = row(tad2);
            if (py_int(at(t1, 2)) < py_int(at(t2, 1)) || at(t1, 0) != at(t2, 0)) {
                break;
            }
            const double min_area = (py_float(at(t1, 2)) - py_float(at(t1, 1))) * percent;
            if ((py_float(at(t1, 2)) - min_area) > py_float(at(t2, 1)) &&
                (py_float(at(t1, 2)) + min_area) <= py_float(at(t2, 2))) {
                // float(t2[2]) - int(t2[1]): the int is converted to a double.
                if ((py_float(at(t1, 2)) - py_float(at(t1, 1))) >
                    (py_float(at(t2, 2)) - static_cast<double>(py_int(at(t2, 1))))) {
                    add_relation_to_list(relations, at(t1, 0), at(t1, 3), at(t2, 3));
                } else {
                    add_relation_to_list(relations, at(t2, 0), at(t2, 3), at(t1, 3));
                }
            } else if (!(py_float(at(t2, 1)) < py_float(at(t1, 1))) &&
                       !int_greater_than_float(py_int(at(t2, 2)), py_float(at(t1, 2)))) {
                add_relation_to_list(relations, at(t1, 0), at(t1, 3), at(t2, 3));
            }
            ++tad2;
        }
        tad2 = tad1 + 2;
        ++tad1;
    }
    return relations;
}

void write_domain_list(std::ostream& out, const RowPool& pool,
                       const std::vector<std::size_t>& list) {
    for (const std::size_t index : list) {
        out << join_row(pool.rows[index]) << '\n';
    }
}

void write_relation_list(std::ostream& out, const std::vector<Relation>& relations) {
    for (const Relation& relation : relations) {
        for (const std::string& child : relation.children) {
            out << relation.chrom << '\t' << relation.parent << '\t' << child << '\n';
        }
    }
}

namespace {

// graphviz.quoting.quote for the identifiers create_tree can produce. Every
// node name is an ID written by add_id, "ID_" followed by a decimal number,
// which graphviz leaves unquoted. Anything else would need graphviz's quoting
// rules, which this port does not implement, so it is refused rather than
// written differently.
const std::string& dot_id(const std::string& name) {
    static const std::set<std::string> keywords = {"node",    "edge",     "graph",
                                                   "digraph", "subgraph", "strict"};
    bool plain = !name.empty() &&
                 (std::isalpha(static_cast<unsigned char>(name[0])) != 0 || name[0] == '_');
    for (const char c : name) {
        if (std::isalnum(static_cast<unsigned char>(c)) == 0 && c != '_') {
            plain = false;
        }
    }
    if (!plain || keywords.count(name) != 0 || name.find(':') != std::string::npos) {
        throw std::logic_error("unexpected graph node name '" + name + "'");
    }
    return name;
}

class Digraph {
  public:
    Digraph(std::string filename, bool strict)
        : filename_(std::move(filename)), strict_(strict) {}

    void node(const std::string& name) { body_ += "\t" + dot_id(name) + "\n"; }
    void edge(const std::string& tail, const std::string& head) {
        body_ += "\t" + dot_id(tail) + " -> " + dot_id(head) + "\n";
    }
    [[nodiscard]] TreeGraph graph(const std::string& chrom) const {
        return TreeGraph{chrom, filename_,
                         std::string(strict_ ? "strict digraph {\n" : "digraph {\n") + body_ +
                             "}\n"};
    }

  private:
    std::string filename_;
    bool strict_;
    std::string body_;
};

void remove_first(std::vector<std::string>& items, const std::string& value) {
    const auto found = std::find(items.begin(), items.end(), value);
    if (found == items.end()) {
        throw PythonError("ValueError: list.remove(x): x not in list");
    }
    items.erase(found);
}

bool contains(const std::vector<std::string>& items, const std::string& value) {
    return std::find(items.begin(), items.end(), value) != items.end();
}

}  // namespace

void create_tree(const std::vector<Relation>& relations, const RowPool& pool,
                 const std::vector<std::size_t>& list, const std::string& prefix,
                 const std::function<void(const TreeGraph&)>& render) {
    std::string name = prefix + "_" + at(relations, 0).chrom;
    // :271 is the only Digraph built with strict=True; the ones for every
    // further chromosome at :297 are not strict. Pinned, not repaired.
    Digraph graph(name, true);

    // create_small_list (:306-319): IDs grouped by runs of one chromosome.
    std::vector<std::pair<std::string, std::vector<std::string>>> small;
    {
        const Row& first = pool.rows[at(list, 0)];
        std::string chrom = at(first, 0);
        small.push_back({at(first, 0), {}});
        for (const std::size_t index : list) {
            const Row& row = pool.rows[index];
            if (at(row, 0) == chrom) {
                small.back().second.push_back(at(row, 3));
            } else {
                chrom = at(row, 0);
                small.push_back({at(row, 0), {at(row, 3)}});
            }
        }
    }

    std::string chrom = relations[0].chrom;
    std::size_t pos_s = 0;
    while (at(small, pos_s).first != chrom) {
        ++pos_s;
    }
    for (std::size_t pos_r = 0; pos_r < relations.size(); ++pos_r) {
        const Relation& relation = relations[pos_r];
        if (chrom == relation.chrom) {
            std::vector<std::string>& ids = small[pos_s].second;
            while (py_int(slice_from(at(ids, 0), 3)) < py_int(slice_from(relation.parent, 3))) {
                graph.node(ids[0]);
                ids.erase(ids.begin());
            }
            if (contains(ids, relation.parent)) {
                remove_first(ids, relation.parent);
            }
            for (const std::string& child : relation.children) {
                graph.edge(relation.parent, child);
                if (contains(ids, child)) {
                    remove_first(ids, child);
                }
            }
        } else {
            // The IDs left over on the finished chromosome become lone nodes.
            // The relation that opened the new chromosome is not added to any
            // graph, because the loop moves past it (:301); and the IDs left
            // over on the last chromosome are never flushed (:302-303). Both
            // pinned.
            std::vector<std::string>& ids = small[pos_s].second;
            while (!ids.empty()) {
                graph.node(ids[0]);
                ids.erase(ids.begin());
            }
            render(graph.graph(chrom));
            chrom = relation.chrom;
            name = prefix + "_" + chrom;
            graph = Digraph(name, false);
            pos_s = 0;
            while (at(small, pos_s).first != chrom) {
                ++pos_s;
            }
        }
    }
    render(graph.graph(chrom));
}

MergedDomains merge_domain_files(const std::vector<std::string>& domain_files,
                                 const ProteinList* proteins,
                                 std::int64_t minimum_number_of_peaks, std::int64_t value) {
    MergedDomains result;
    // create_list_with_protein (:394-400). merge_protein rebinds only its
    // local name, so every domain file bins the original protein list with
    // its own bin size.
    const auto load = [&](const std::string& path) {
        DomainList list = create_list_of_file(result.pool, path);
        if (proteins != nullptr) {
            const MergedProtein binned =
                merge_protein(*proteins, list.bin_size, minimum_number_of_peaks);
            compare_boundaries_protein(result.pool, list.rows, binned);
        }
        return list.rows;
    };
    result.merged = load(at(domain_files, 0));
    if (domain_files.size() > 1) {
        std::vector<std::vector<std::size_t>> others;
        for (std::size_t i = 1; i < domain_files.size(); ++i) {
            others.push_back(load(domain_files[i]));
        }
        for (const auto& other : others) {
            result.merged = merge_list(result.pool, result.merged, other, value);
        }
    }
    add_id(result.pool, result.merged);
    return result;
}

bool is_graphviz_format(const std::string& format) {
    static const std::set<std::string> formats = {
        "bmp",      "canon",     "dot",      "gv",        "xdot",     "xdot1.2",
        "xdot1.4",  "cgimage",   "cmap",     "eps",       "exr",      "fig",
        "gd",       "gd2",       "gif",      "gtk",       "ico",      "imap",
        "cmapx",    "imap_np",   "cmapx_np", "ismap",     "jp2",      "jpg",
        "jpeg",     "jpe",       "json",     "json0",     "dot_json", "xdot_json",
        "pct",      "pict",      "pdf",      "pic",       "plain",    "plain-ext",
        "png",      "pov",       "ps",       "ps2",       "psd",      "sgi",
        "svg",      "svgz",      "tga",      "tif",       "tiff",     "tk",
        "vml",      "vmlz",      "vrml",     "wbmp",      "webp",     "xlib",
        "x11"};
    std::string lowered;
    for (const char c : format) {
        lowered.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    return formats.count(lowered) != 0;
}

}  // namespace hicx::merge_domains
