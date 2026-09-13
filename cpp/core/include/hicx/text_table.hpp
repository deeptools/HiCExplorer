// A pandas.read_csv(sep='\t', header=None) table, and the writer that matches
// DataFrame.to_csv(sep='\t', header=False, index=False).
//
// Three tier 2 tools (hicValidateLocations, hicMergeLoops and, through
// pybedtools, anything that goes near BedTool.from_dataframe) move their data
// through pandas rather than through the matrix layer, and pandas is
// observable in the output in two places that a naive TSV reader gets wrong:
//
//  1. **Column dtype inference.** read_csv assigns one dtype per column: int64
//     if every field parses as an integer, otherwise float64 if every field
//     parses as a float or is one of pandas' NA spellings, otherwise object.
//     The dtype decides how to_csv prints the column again: an integer column
//     prints without a decimal point, a float column prints repr(float), an
//     object column prints the text unchanged. A loop file whose score column
//     happens to hold only whole numbers therefore comes back as integers.
//
//  2. **The float parser is not strtod.** The C parser's default converter is
//     `precise_xstrtod` (pandas/_libs/src/parser/tokenizer.c), which
//     accumulates at most 17 significant decimal digits into a double and then
//     scales by an exact power of ten from a table. That is not the correctly
//     rounded conversion strtod performs, and the difference is visible:
//     "0.0019430592210407135" in hicValidateLocations' loop file comes back as
//     the double whose repr is "0.0019430592210407", and that shorter string
//     is what the tool writes out. Using strtod here would reproduce the
//     input's digits instead, which is arguably more correct and is certainly
//     not what the reference does.
//
// So `pandas_strtod` below is precise_xstrtod, not a wrapper over std::strtod,
// and the difference is the reason hicValidateLocations can be compared byte
// for byte at E0 rather than only at a float tolerance.

#ifndef HICX_TEXT_TABLE_HPP
#define HICX_TEXT_TABLE_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace hicx {

// pandas' precise_xstrtod. `ok` reports whether the whole string was consumed
// as a number; a partially consumed string is not a number to pandas either.
[[nodiscard]] double pandas_strtod(const std::string& text, bool* ok);

// True for the strings pandas treats as missing in a numeric column
// (pandas.io.parsers STR_NA_VALUES).
[[nodiscard]] bool is_pandas_na(const std::string& text);

enum class ColumnType { Int64, Float64, String };

struct TableColumn {
    ColumnType type = ColumnType::String;
    std::vector<std::int64_t> ints;
    std::vector<double> floats;
    std::vector<std::string> strings;

    [[nodiscard]] std::size_t size() const;
    // The text to_csv would write for one cell. `na_rep` is written for a NaN
    // in a float column, which is how pandas spells a missing value; the
    // default is the empty string that DataFrame.to_csv uses, while
    // BedTool.from_dataframe passes ".".
    [[nodiscard]] std::string text(std::size_t row, const std::string& na_rep) const;
    // str(value), which is what 'chr' + column.astype(str) concatenates.
    [[nodiscard]] std::string as_str(std::size_t row) const;
    // The value as an integer. Throws for a string column, and for a float
    // column truncates towards zero the way astype(int) does.
    [[nodiscard]] std::int64_t as_int(std::size_t row) const;
};

class TextTable {
  public:
    // pandas.read_csv(path, sep='\t', header=None). Blank lines are skipped,
    // every line must hold the same number of fields, and there is no comment
    // character, all three as pandas has them with these arguments.
    [[nodiscard]] static TextTable read_tsv(const std::string& path);

    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t cols() const noexcept { return columns_.size(); }
    [[nodiscard]] TableColumn& column(std::size_t index) { return columns_[index]; }
    [[nodiscard]] const TableColumn& column(std::size_t index) const {
        return columns_[index];
    }
    [[nodiscard]] const std::vector<TableColumn>& columns() const noexcept {
        return columns_;
    }

    // Appends the rows of `other`, which must have the same number of columns.
    // pd.concat re-infers nothing: a column that is int64 in one frame and
    // float64 in the other becomes float64, which is numpy's promotion and is
    // what pd.concat does for these dtypes.
    void append(const TextTable& other);

    // The table restricted to `rows`, in that order.
    [[nodiscard]] TextTable select_rows(const std::vector<std::size_t>& rows) const;
    // The table restricted to `columns`, in that order.
    [[nodiscard]] TextTable select_columns(const std::vector<std::size_t>& columns) const;

    // The round trip BedTool.from_dataframe(df).<op>().to_dataframe() performs:
    // the frame is written as text and read back, so every column's dtype is
    // inferred again from that text. It is not the identity. A float column
    // holding NaN is written as `na_rep` and read back as an object column,
    // because read_csv does not recognise pybedtools' '.' as a missing value.
    [[nodiscard]] TextTable reparse(const std::string& na_rep = ".") const;

    // to_csv(sep='\t', header=False, index=False).
    void write_tsv(const std::string& path, const std::string& na_rep = "") const;
    // The same, rendered into a string.
    [[nodiscard]] std::string to_tsv(const std::string& na_rep = "") const;

    // 'chr' + column.astype(str), which turns the column into an object column.
    void add_chr_prefix(std::size_t index);
    // column.str.lstrip('chr'), a character-set strip, not a prefix strip, and
    // an AttributeError on a non-object column. The exception is reproduced as
    // a std::runtime_error.
    void remove_chr_prefix(std::size_t index);

    // drop_duplicates(): keeps the first row of every group of identical rows.
    void drop_duplicates();
    // drop_duplicates(keep=False): drops every row that occurs more than once.
    void drop_duplicates_drop_all();

    static TextTable with_columns(std::vector<TableColumn> columns);

  private:
    std::vector<TableColumn> columns_;
    std::size_t rows_ = 0;
};

}  // namespace hicx

#endif  // HICX_TEXT_TABLE_HPP
