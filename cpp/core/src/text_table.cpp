#include "hicx/text_table.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <limits>
#include <unordered_set>

#include "hicx/numpy_compat.hpp"

namespace hicx {

namespace {

// The power of ten table precise_xstrtod scales by. pandas writes it as 309
// decimal literals; building it with strtod at static initialisation gives the
// same correctly rounded doubles the compiler would produce for the literals,
// without 309 lines of them.
const std::array<double, 309>& power_of_ten_table() {
    static const std::array<double, 309> table = [] {
        std::array<double, 309> values{};
        for (int i = 0; i <= 308; ++i) {
            const std::string literal = "1e" + std::to_string(i);
            values[static_cast<std::size_t>(i)] = std::strtod(literal.c_str(), nullptr);
        }
        return values;
    }();
    return table;
}

inline bool is_digit(char c) { return c >= '0' && c <= '9'; }

// pandas.io.parsers.STR_NA_VALUES.
const std::unordered_set<std::string>& na_values() {
    static const std::unordered_set<std::string> values = {
        "",     "#N/A", "#N/A N/A", "#NA",  "-1.#IND", "-1.#QNAN", "-NaN",
        "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A",     "NA",       "NULL",
        "NaN",  "None", "n/a",      "nan",  "null"};
    return values;
}

bool parse_int64(const std::string& text, std::int64_t* out) {
    if (text.empty()) {
        return false;
    }
    std::size_t i = 0;
    bool negative = false;
    if (text[0] == '+' || text[0] == '-') {
        negative = text[0] == '-';
        i = 1;
    }
    if (i >= text.size()) {
        return false;
    }
    unsigned long long magnitude = 0;
    for (; i < text.size(); ++i) {
        if (!is_digit(text[i])) {
            return false;
        }
        const unsigned long long digit = static_cast<unsigned long long>(text[i] - '0');
        if (magnitude > (0xFFFFFFFFFFFFFFFFULL - digit) / 10ULL) {
            return false;  // does not fit, pandas falls through to float64
        }
        magnitude = magnitude * 10ULL + digit;
    }
    if (negative) {
        if (magnitude > 9223372036854775808ULL) {
            return false;
        }
        *out = magnitude == 9223372036854775808ULL
                   ? std::numeric_limits<std::int64_t>::min()
                   : -static_cast<std::int64_t>(magnitude);
        return true;
    }
    if (magnitude > 9223372036854775807ULL) {
        return false;
    }
    *out = static_cast<std::int64_t>(magnitude);
    return true;
}

std::string int_to_string(std::int64_t value) { return std::to_string(value); }

}  // namespace

double pandas_strtod(const std::string& text, bool* ok) {
    // pandas/_libs/src/parser/tokenizer.c, precise_xstrtod, with tsep = '\0',
    // decimal = '.', sci = 'E' and skip_trailing = 1.
    const std::array<double, 309>& e = power_of_ten_table();
    const char* p = text.c_str();
    if (ok != nullptr) {
        *ok = false;
    }

    while (std::isspace(static_cast<unsigned char>(*p)) != 0) {
        ++p;
    }
    bool negative = false;
    if (*p == '-') {
        negative = true;
        ++p;
    } else if (*p == '+') {
        ++p;
    }

    double number = 0.0;
    int exponent = 0;
    int num_digits = 0;
    int num_decimals = 0;
    const int max_digits = 17;

    while (is_digit(*p)) {
        if (num_digits < max_digits) {
            number = number * 10.0 + static_cast<double>(*p - '0');
            ++num_digits;
        } else {
            ++exponent;
        }
        ++p;
    }

    if (*p == '.') {
        ++p;
        while (num_digits < max_digits && is_digit(*p)) {
            number = number * 10.0 + static_cast<double>(*p - '0');
            ++p;
            ++num_digits;
            ++num_decimals;
        }
        if (num_digits >= max_digits) {
            while (is_digit(*p)) {
                ++p;
            }
        }
        exponent -= num_decimals;
    }

    if (num_digits == 0) {
        return 0.0;
    }
    if (negative) {
        number = -number;
    }

    if (*p == 'e' || *p == 'E') {
        ++p;
        bool exponent_negative = false;
        if (*p == '-') {
            exponent_negative = true;
            ++p;
        } else if (*p == '+') {
            ++p;
        }
        int exponent_digits = 0;
        int n = 0;
        while (is_digit(*p)) {
            n = n * 10 + (*p - '0');
            ++exponent_digits;
            ++p;
        }
        if (exponent_digits == 0) {
            return 0.0;
        }
        exponent += exponent_negative ? -n : n;
    }

    if (exponent > 308) {
        return HUGE_VAL;
    }
    if (exponent > 0) {
        number *= e[static_cast<std::size_t>(exponent)];
    } else if (exponent < -308) {
        if (exponent < -616) {
            number = 0.0;
        } else {
            number /= e[static_cast<std::size_t>(-308 - exponent)];
            number /= e[308];
        }
    } else {
        number /= e[static_cast<std::size_t>(-exponent)];
    }

    while (std::isspace(static_cast<unsigned char>(*p)) != 0) {
        ++p;
    }
    if (ok != nullptr) {
        *ok = *p == '\0';
    }
    return number;
}

bool is_pandas_na(const std::string& text) {
    return na_values().find(text) != na_values().end();
}

std::size_t TableColumn::size() const {
    switch (type) {
        case ColumnType::Int64:
            return ints.size();
        case ColumnType::Float64:
            return floats.size();
        default:
            return strings.size();
    }
}

std::string TableColumn::text(std::size_t row, const std::string& na_rep) const {
    switch (type) {
        case ColumnType::Int64:
            return int_to_string(ints[row]);
        case ColumnType::Float64: {
            const double value = floats[row];
            if (std::isnan(value)) {
                return na_rep;
            }
            return npy::float_repr(value);
        }
        default:
            return strings[row];
    }
}

std::string TableColumn::as_str(std::size_t row) const {
    switch (type) {
        case ColumnType::Int64:
            return int_to_string(ints[row]);
        case ColumnType::Float64:
            return npy::float_repr(floats[row]);
        default:
            return strings[row];
    }
}

std::int64_t TableColumn::as_int(std::size_t row) const {
    switch (type) {
        case ColumnType::Int64:
            return ints[row];
        case ColumnType::Float64:
            return static_cast<std::int64_t>(floats[row]);
        default:
            throw std::runtime_error("column holds text where a number is required: '" +
                                     strings[row] + "'");
    }
}

TextTable TextTable::with_columns(std::vector<TableColumn> columns) {
    TextTable table;
    table.rows_ = columns.empty() ? 0 : columns.front().size();
    table.columns_ = std::move(columns);
    return table;
}

namespace {

// One column per field position, with pandas' dtype inference: int64 if every
// field parses as an integer, otherwise float64 if every field parses as a
// float or is one of pandas' NA spellings, otherwise object. That is the order
// the C parser tries them in.
std::vector<TableColumn> infer_columns(
    const std::vector<std::vector<std::string>>& fields, std::size_t ncols) {
    std::vector<TableColumn> columns(ncols);
    for (std::size_t c = 0; c < ncols; ++c) {
        TableColumn& column = columns[c];
        bool all_int = true;
        for (const std::vector<std::string>& row : fields) {
            std::int64_t value = 0;
            if (!parse_int64(row[c], &value)) {
                all_int = false;
                break;
            }
        }
        if (all_int) {
            column.type = ColumnType::Int64;
            column.ints.reserve(fields.size());
            for (const std::vector<std::string>& row : fields) {
                std::int64_t value = 0;
                parse_int64(row[c], &value);
                column.ints.push_back(value);
            }
            continue;
        }
        bool all_float = true;
        for (const std::vector<std::string>& row : fields) {
            if (is_pandas_na(row[c])) {
                continue;
            }
            bool ok = false;
            const double parsed = pandas_strtod(row[c], &ok);
            (void)parsed;
            if (!ok) {
                all_float = false;
                break;
            }
        }
        if (all_float) {
            column.type = ColumnType::Float64;
            column.floats.reserve(fields.size());
            for (const std::vector<std::string>& row : fields) {
                if (is_pandas_na(row[c])) {
                    column.floats.push_back(std::nan(""));
                    continue;
                }
                bool ok = false;
                column.floats.push_back(pandas_strtod(row[c], &ok));
            }
            continue;
        }
        column.type = ColumnType::String;
        column.strings.reserve(fields.size());
        for (const std::vector<std::string>& row : fields) {
            column.strings.push_back(row[c]);
        }
    }
    return columns;
}

std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> row;
    std::size_t start = 0;
    while (true) {
        const std::size_t marker = line.find('\t', start);
        if (marker == std::string::npos) {
            row.push_back(line.substr(start));
            break;
        }
        row.push_back(line.substr(start, marker - start));
        start = marker + 1;
    }
    return row;
}

}  // namespace

TextTable TextTable::read_tsv(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open '" + path + "'");
    }

    std::vector<std::vector<std::string>> fields;
    std::string line;
    std::size_t ncols = 0;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;  // skip_blank_lines, pandas' default
        }
        std::vector<std::string> row = split_tabs(line);
        if (ncols == 0) {
            ncols = row.size();
        } else if (row.size() != ncols) {
            throw std::runtime_error(
                "Error tokenizing data. C error: Expected " + std::to_string(ncols) +
                " fields in line " + std::to_string(fields.size() + 1) + ", saw " +
                std::to_string(row.size()));
        }
        fields.push_back(std::move(row));
    }

    TextTable table;
    table.rows_ = fields.size();
    table.columns_ = infer_columns(fields, ncols);
    return table;
}

TextTable TextTable::reparse(const std::string& na_rep) const {
    if (rows_ == 0) {
        return *this;
    }
    const std::string text = to_tsv(na_rep);
    std::vector<std::vector<std::string>> fields;
    fields.reserve(rows_);
    std::size_t begin = 0;
    while (begin < text.size()) {
        const std::size_t line_end = text.find('\n', begin);
        const std::string line = text.substr(begin, line_end - begin);
        begin = line_end == std::string::npos ? text.size() : line_end + 1;
        if (line.empty()) {
            continue;
        }
        fields.push_back(split_tabs(line));
    }
    TextTable table;
    table.rows_ = fields.size();
    table.columns_ = infer_columns(fields, fields.empty() ? 0 : fields.front().size());
    return table;
}

void TextTable::append(const TextTable& other) {
    if (other.cols() != cols()) {
        throw std::runtime_error("cannot concatenate tables with different widths");
    }
    for (std::size_t c = 0; c < columns_.size(); ++c) {
        TableColumn& target = columns_[c];
        const TableColumn& source = other.columns_[c];
        if (target.type == source.type) {
            switch (target.type) {
                case ColumnType::Int64:
                    target.ints.insert(target.ints.end(), source.ints.begin(),
                                       source.ints.end());
                    break;
                case ColumnType::Float64:
                    target.floats.insert(target.floats.end(), source.floats.begin(),
                                         source.floats.end());
                    break;
                default:
                    target.strings.insert(target.strings.end(), source.strings.begin(),
                                          source.strings.end());
                    break;
            }
            continue;
        }
        // Mixed dtypes: int64 with float64 promotes to float64, anything with
        // a string column becomes object.
        if (target.type != ColumnType::String && source.type != ColumnType::String) {
            if (target.type == ColumnType::Int64) {
                std::vector<double> promoted;
                promoted.reserve(target.ints.size());
                for (std::int64_t value : target.ints) {
                    promoted.push_back(static_cast<double>(value));
                }
                target.type = ColumnType::Float64;
                target.floats = std::move(promoted);
                target.ints.clear();
            }
            for (std::size_t r = 0; r < other.rows_; ++r) {
                target.floats.push_back(static_cast<double>(source.ints[r]));
            }
            continue;
        }
        if (target.type != ColumnType::String) {
            std::vector<std::string> promoted;
            promoted.reserve(target.size());
            for (std::size_t r = 0; r < rows_; ++r) {
                promoted.push_back(target.as_str(r));
            }
            target.type = ColumnType::String;
            target.strings = std::move(promoted);
            target.ints.clear();
            target.floats.clear();
        }
        for (std::size_t r = 0; r < other.rows_; ++r) {
            target.strings.push_back(source.as_str(r));
        }
    }
    rows_ += other.rows_;
}

TextTable TextTable::select_rows(const std::vector<std::size_t>& rows) const {
    TextTable result;
    result.rows_ = rows.size();
    result.columns_.resize(columns_.size());
    for (std::size_t c = 0; c < columns_.size(); ++c) {
        const TableColumn& source = columns_[c];
        TableColumn& target = result.columns_[c];
        target.type = source.type;
        switch (source.type) {
            case ColumnType::Int64:
                target.ints.reserve(rows.size());
                for (std::size_t r : rows) {
                    target.ints.push_back(source.ints[r]);
                }
                break;
            case ColumnType::Float64:
                target.floats.reserve(rows.size());
                for (std::size_t r : rows) {
                    target.floats.push_back(source.floats[r]);
                }
                break;
            default:
                target.strings.reserve(rows.size());
                for (std::size_t r : rows) {
                    target.strings.push_back(source.strings[r]);
                }
                break;
        }
    }
    return result;
}

TextTable TextTable::select_columns(const std::vector<std::size_t>& columns) const {
    TextTable result;
    result.rows_ = rows_;
    result.columns_.reserve(columns.size());
    for (std::size_t c : columns) {
        result.columns_.push_back(columns_[c]);
    }
    return result;
}

std::string TextTable::to_tsv(const std::string& na_rep) const {
    std::string out;
    for (std::size_t r = 0; r < rows_; ++r) {
        for (std::size_t c = 0; c < columns_.size(); ++c) {
            if (c != 0) {
                out.push_back('\t');
            }
            out += columns_[c].text(r, na_rep);
        }
        out.push_back('\n');
    }
    return out;
}

void TextTable::write_tsv(const std::string& path, const std::string& na_rep) const {
    std::ofstream output(path, std::ios::binary);
    if (!output) {
        throw std::runtime_error("cannot write '" + path + "'");
    }
    output << to_tsv(na_rep);
}

void TextTable::add_chr_prefix(std::size_t index) {
    TableColumn& column = columns_[index];
    std::vector<std::string> prefixed;
    prefixed.reserve(rows_);
    for (std::size_t r = 0; r < rows_; ++r) {
        prefixed.push_back("chr" + column.as_str(r));
    }
    column.type = ColumnType::String;
    column.strings = std::move(prefixed);
    column.ints.clear();
    column.floats.clear();
}

void TextTable::remove_chr_prefix(std::size_t index) {
    TableColumn& column = columns_[index];
    if (column.type != ColumnType::String) {
        // Series.str on a non-object column raises
        // AttributeError: Can only use .str accessor with string values.
        throw std::runtime_error(
            "Can only use .str accessor with string values; the chromosome column "
            "of this file was inferred as a numeric dtype, so --chrPrefixProtein "
            "remove cannot be applied to it");
    }
    for (std::string& value : column.strings) {
        // str.lstrip('chr') strips any leading character out of the set
        // {'c','h','r'}, not the prefix "chr".
        std::size_t first = 0;
        while (first < value.size() &&
               (value[first] == 'c' || value[first] == 'h' || value[first] == 'r')) {
            ++first;
        }
        value.erase(0, first);
    }
}

namespace {

std::string row_key(const TextTable& table, std::size_t row) {
    std::string key;
    for (std::size_t c = 0; c < table.cols(); ++c) {
        key += table.column(c).text(row, "<NA>");
        key.push_back('\x01');
    }
    return key;
}

}  // namespace

void TextTable::drop_duplicates() {
    std::unordered_set<std::string> seen;
    std::vector<std::size_t> keep;
    keep.reserve(rows_);
    for (std::size_t r = 0; r < rows_; ++r) {
        if (seen.insert(row_key(*this, r)).second) {
            keep.push_back(r);
        }
    }
    *this = select_rows(keep);
}

void TextTable::drop_duplicates_drop_all() {
    std::unordered_map<std::string, std::size_t> counts;
    for (std::size_t r = 0; r < rows_; ++r) {
        ++counts[row_key(*this, r)];
    }
    std::vector<std::size_t> keep;
    keep.reserve(rows_);
    for (std::size_t r = 0; r < rows_; ++r) {
        if (counts[row_key(*this, r)] == 1) {
            keep.push_back(r);
        }
    }
    *this = select_rows(keep);
}

}  // namespace hicx
