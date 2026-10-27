// hicQC.html of hicPrepareQCreport: pandas' Styler.to_html of the QC tables
// inside hicexplorer/qc_template.html (compiled in).
//
// Shared by the hicPrepareQCreport executable (through build_matrix.cpp) and
// by hicexplorer.hicPrepareQCreport.save_html, whose DataFrames the extension
// reads into QcDisplayTable.

#ifndef HICX_QC_HTML_HPP
#define HICX_QC_HTML_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace hicx {

// One cell of a DataFrame: a Python/numpy integer, a float, or text.
struct QcStylerCell {
    bool is_int = false;
    bool is_float = false;
    std::int64_t integer = 0;
    double real = 0.0;
    std::string text;
};

// A DataFrame with a single level index and single level columns. The labels
// are the strings pandas' default formatter makes of them.
struct QcDisplayTable {
    std::optional<std::string> index_name;
    std::optional<std::string> columns_name;
    std::vector<std::string> index;
    std::vector<std::string> columns;
    std::vector<std::vector<QcStylerCell>> cells;  // per column, per row
};

enum class QcFormatter {
    // save_html's lambda x: '{:,}'.format(x) if x > 1 else '{:.2%}'.format(x)
    Derived,
    // Styler's default: integers and text as they are, floats with 6 decimals
    Plain
};

// DataFrame.style.format(formatter).to_html() in pandas 2.2 (html_style.tpl
// and html_table.tpl) for a table without styles, caption or attributes.
// Styler escapes nothing by default. The table id carries five random
// hexadecimal characters, as uuid4().hex[:5] does there. Throws PythonError
// (TypeError) where the Python formatter would raise for a text cell of the
// Derived formatter.
[[nodiscard]] std::string qc_styler_html(const QcDisplayTable& table, QcFormatter formatter);

// hicexplorer/qc_template.html with the four derived tables and the QC table
// put in (save_html; the QC table loses 'Min rest. site distance' and 'Max
// library insert size').
[[nodiscard]] std::string qc_html_document(const QcDisplayTable& unmapped,
                                           const QcDisplayTable& discarded,
                                           const QcDisplayTable& distance,
                                           const QcDisplayTable& orientation,
                                           const QcDisplayTable& all);

}  // namespace hicx

#endif  // HICX_QC_HTML_HPP
