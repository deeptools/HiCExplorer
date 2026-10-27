#include "hicx/qc_html.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <random>

#include "hicx/numpy_compat.hpp"
#include "hicx/python_error.hpp"

#include "qc_template_html.inc"  // kQcTemplateHtml, generated from hicexplorer/qc_template.html

namespace hicx {

namespace {

// f"{x:.<precision>f}": Python prints nan and the infinities as words.
std::string qc_fixed(double value, int precision) {
    if (std::isnan(value)) {
        return "nan";
    }
    if (std::isinf(value)) {
        return value > 0 ? "inf" : "-inf";
    }
    char buffer[400];
    std::snprintf(buffer, sizeof(buffer), "%.*f", precision, value);
    return buffer;
}

// '{:.2%}'.format(x)
std::string qc_format_percent(double value) {
    if (std::isnan(value)) {
        return "nan%";
    }
    return qc_fixed(value * 100.0, 2) + "%";
}

// The thousands separator of '{:,}' on the leading integer digits of `text`.
std::string qc_group_thousands(const std::string& text) {
    const std::size_t start = (!text.empty() && text[0] == '-') ? 1 : 0;
    std::size_t end = start;
    while (end < text.size() && std::isdigit(static_cast<unsigned char>(text[end]))) {
        ++end;
    }
    const std::string digits = text.substr(start, end - start);
    std::string grouped;
    for (std::size_t i = 0; i < digits.size(); ++i) {
        if (i > 0 && (digits.size() - i) % 3 == 0) {
            grouped += ',';
        }
        grouped += digits[i];
    }
    return text.substr(0, start) + grouped + text.substr(end);
}

// The formatter save_html gives the four derived tables:
// lambda x: '{:,}'.format(x) if x > 1 else '{:.2%}'.format(x)
// ('{:,}' of a float is its repr with grouped integer digits, unless the repr
// is in exponent notation or a word).
std::string derived_formatter(const QcStylerCell& cell) {
    if (cell.is_int) {
        return cell.integer > 1 ? qc_group_thousands(std::to_string(cell.integer))
                                : qc_format_percent(static_cast<double>(cell.integer));
    }
    if (cell.is_float) {
        if (!(cell.real > 1)) {
            return qc_format_percent(cell.real);
        }
        const std::string repr = npy::float_repr(cell.real);
        if (repr.find('e') != std::string::npos || repr.find("inf") != std::string::npos) {
            return repr;
        }
        return qc_group_thousands(repr);
    }
    throw PythonError("TypeError", "'>' not supported between instances of 'str' and 'int'");
}

// Styler's default formatter (pandas.io.formats.style_render._default_formatter
// with precision 6): floats with six decimals, integers and text as they are.
std::string plain_formatter(const QcStylerCell& cell) {
    if (cell.is_int) {
        return std::to_string(cell.integer);
    }
    if (cell.is_float) {
        return qc_fixed(cell.real, 6);
    }
    return cell.text;
}

// Five random hexadecimal characters, as Styler's uuid4().hex[:5].
std::string qc_styler_uuid() {
    static std::mt19937_64 engine{std::random_device{}()};
    static const char* const hex = "0123456789abcdef";
    std::string id;
    for (int i = 0; i < 5; ++i) {
        id += hex[engine() % 16];
    }
    return id;
}

void append_header_cell(std::string& out, const std::string& uuid, const std::string& id,
                        const std::string& css, const std::string& text) {
    out += "      <th id=\"T_" + uuid + "_" + id + "\" class=\"" + css + "\" >" + text + "</th>\n";
}

// str.replace: every occurrence, left to right.
void qc_replace_all(std::string& text, const std::string& pattern, const std::string& replacement) {
    std::string out;
    std::size_t from = 0;
    for (std::size_t at = text.find(pattern); at != std::string::npos;
         at = text.find(pattern, from)) {
        out.append(text, from, at - from);
        out += replacement;
        from = at + pattern.size();
    }
    out.append(text, from, std::string::npos);
    text = std::move(out);
}

}  // namespace

std::string qc_styler_html(const QcDisplayTable& table, QcFormatter formatter) {
    const std::string uuid = qc_styler_uuid();
    std::string out = "<style type=\"text/css\">\n</style>\n<table id=\"T_" + uuid + "\">\n";
    out += "  <thead>\n";
    if (!table.columns.empty()) {
        // a table without columns has no column header row
        out += "    <tr>\n";
        if (table.columns_name.has_value()) {
            out += "      <th class=\"index_name level0\" >" + *table.columns_name + "</th>\n";
        } else {
            out += "      <th class=\"blank level0\" >&nbsp;</th>\n";
        }
        for (std::size_t c = 0; c < table.columns.size(); ++c) {
            const std::string col = std::to_string(c);
            append_header_cell(out, uuid, "level0_col" + col, "col_heading level0 col" + col,
                               table.columns[c]);
        }
        out += "    </tr>\n";
    }
    if (table.index_name.has_value()) {
        out += "    <tr>\n      <th class=\"index_name level0\" >" + *table.index_name + "</th>\n";
        for (std::size_t c = 0; c < table.columns.size(); ++c) {
            out += "      <th class=\"blank col" + std::to_string(c) + "\" >&nbsp;</th>\n";
        }
        out += "    </tr>\n";
    }
    out += "  </thead>\n  <tbody>\n";
    for (std::size_t r = 0; r < table.index.size(); ++r) {
        const std::string row = std::to_string(r);
        out += "    <tr>\n";
        append_header_cell(out, uuid, "level0_row" + row, "row_heading level0 row" + row,
                           table.index[r]);
        for (std::size_t c = 0; c < table.columns.size(); ++c) {
            const std::string col = std::to_string(c);
            const QcStylerCell& cell = table.cells[c][r];
            const std::string shown =
                formatter == QcFormatter::Derived ? derived_formatter(cell) : plain_formatter(cell);
            out += "      <td id=\"T_" + uuid + "_row" + row + "_col" + col + "\" class=\"data row" +
                   row + " col" + col + "\" >" + shown + "</td>\n";
        }
        out += "    </tr>\n";
    }
    out += "  </tbody>\n</table>\n";
    return out;
}

std::string qc_html_document(const QcDisplayTable& unmapped, const QcDisplayTable& discarded,
                             const QcDisplayTable& distance, const QcDisplayTable& orientation,
                             const QcDisplayTable& all) {
    std::string html(reinterpret_cast<const char*>(kQcTemplateHtml), sizeof(kQcTemplateHtml));
    qc_replace_all(html, "%%TABLE_UNMAP%%", qc_styler_html(unmapped, QcFormatter::Derived));
    qc_replace_all(html, "%%TABLE_DISCARDED%%", qc_styler_html(discarded, QcFormatter::Derived));
    qc_replace_all(html, "%%TABLE_DISTANCE%%", qc_styler_html(distance, QcFormatter::Derived));
    qc_replace_all(html, "%%TABLE_ORIENTATION%%",
                   qc_styler_html(orientation, QcFormatter::Derived));
    QcDisplayTable kept = all;
    for (const char* dropped : {"Min rest. site distance", "Max library insert size"}) {
        const auto it = std::find(kept.columns.begin(), kept.columns.end(), dropped);
        if (it != kept.columns.end()) {
            kept.cells.erase(kept.cells.begin() + (it - kept.columns.begin()));
            kept.columns.erase(it);
        }
    }
    qc_replace_all(html, "%%TABLE%%", qc_styler_html(kept, QcFormatter::Plain));
    return html;
}

}  // namespace hicx
