#include "hicx/text_formats.hpp"

#include <zlib.h>

#include <algorithm>
#include <array>
#include <charconv>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "hicx/bins.hpp"
#include "hicx/numpy_compat.hpp"

namespace hicx {

namespace {

// --------------------------------------------------------------------------
// small text helpers

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

// A buffered line reader that transparently handles a gzipped and a plain
// file, which is what hicmatrix.utilities.opener does by sniffing the two
// magic bytes. zlib's gz layer already reads an uncompressed stream verbatim,
// so one code path covers both.
class LineReader {
  public:
    explicit LineReader(const std::string& path) : path_(path) {
        file_ = gzopen(path.c_str(), "rb");
        if (file_ == nullptr) {
            fail("cannot open " + path);
        }
        gzbuffer(file_, 1u << 20);
    }
    LineReader(const LineReader&) = delete;
    LineReader& operator=(const LineReader&) = delete;
    ~LineReader() {
        if (file_ != nullptr) {
            gzclose(file_);
        }
    }

    // Returns false at end of file. The trailing newline is stripped; a final
    // line without one is still returned.
    bool next(std::string& line) {
        line.clear();
        while (true) {
            if (position_ < filled_) {
                const char* start = buffer_.data() + position_;
                const char* found = static_cast<const char*>(
                    std::memchr(start, '\n', filled_ - position_));
                if (found != nullptr) {
                    line.append(start, static_cast<std::size_t>(found - start));
                    position_ += static_cast<std::size_t>(found - start) + 1;
                    return true;
                }
                line.append(start, filled_ - position_);
                position_ = filled_;
            }
            const int read = gzread(file_, buffer_.data(),
                                    static_cast<unsigned>(buffer_.size()));
            if (read < 0) {
                fail("cannot read " + path_);
            }
            if (read == 0) {
                return !line.empty();
            }
            filled_ = static_cast<std::size_t>(read);
            position_ = 0;
        }
    }

  private:
    std::string path_;
    gzFile file_ = nullptr;
    std::vector<char> buffer_ = std::vector<char>(1u << 20);
    std::size_t filled_ = 0;
    std::size_t position_ = 0;
};

// A writer that either compresses with gzip or does not, so that the two text
// writers share one buffering path.
class TextWriter {
  public:
    TextWriter(const std::string& path, bool compressed) : path_(path) {
        if (compressed) {
            // gzip.open(..., 'wt') uses compresslevel 9.
            gz_ = gzopen(path.c_str(), "wb9");
            if (gz_ == nullptr) {
                fail("cannot write " + path);
            }
        } else {
            plain_ = std::fopen(path.c_str(), "wb");
            if (plain_ == nullptr) {
                fail("cannot write " + path);
            }
        }
        buffer_.reserve(kFlushAt + 4096);
    }
    TextWriter(const TextWriter&) = delete;
    TextWriter& operator=(const TextWriter&) = delete;
    ~TextWriter() {
        try {
            flush();
        } catch (const std::exception&) {  // NOLINT: a destructor cannot throw
        }
        if (gz_ != nullptr) {
            gzclose(gz_);
        }
        if (plain_ != nullptr) {
            std::fclose(plain_);
        }
    }

    void write(const std::string& text) {
        buffer_ += text;
        if (buffer_.size() >= kFlushAt) {
            flush();
        }
    }
    void write(char character) {
        buffer_.push_back(character);
        if (buffer_.size() >= kFlushAt) {
            flush();
        }
    }
    void flush() {
        if (buffer_.empty()) {
            return;
        }
        if (gz_ != nullptr) {
            if (gzwrite(gz_, buffer_.data(), static_cast<unsigned>(buffer_.size())) <= 0) {
                fail("cannot write " + path_);
            }
        } else if (std::fwrite(buffer_.data(), 1, buffer_.size(), plain_) !=
                   buffer_.size()) {
            fail("cannot write " + path_);
        }
        buffer_.clear();
    }

  private:
    static constexpr std::size_t kFlushAt = 1u << 20;
    std::string path_;
    gzFile gz_ = nullptr;
    std::FILE* plain_ = nullptr;
    std::string buffer_;
};

std::vector<std::string> split(const std::string& text, char delimiter) {
    std::vector<std::string> parts;
    std::size_t start = 0;
    while (true) {
        const std::size_t found = text.find(delimiter, start);
        if (found == std::string::npos) {
            parts.push_back(text.substr(start));
            return parts;
        }
        parts.push_back(text.substr(start, found - start));
        start = found + 1;
    }
}

// str.strip(): Python strips ASCII whitespace, which includes the tab that
// ends the homer header line.
std::string strip(const std::string& text) {
    const auto is_space = [](unsigned char c) {
        return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f';
    };
    std::size_t first = 0;
    while (first < text.size() && is_space(static_cast<unsigned char>(text[first]))) {
        ++first;
    }
    std::size_t last = text.size();
    while (last > first && is_space(static_cast<unsigned char>(text[last - 1]))) {
        --last;
    }
    return text.substr(first, last - first);
}

std::int64_t parse_int(const std::string& text, const std::string& context) {
    const std::string trimmed = strip(text);
    std::size_t consumed = 0;
    std::int64_t value = 0;
    try {
        value = std::stoll(trimmed, &consumed);
    } catch (const std::exception&) {
        fail(context + ": '" + trimmed + "' is not an integer");
    }
    if (consumed != trimmed.size()) {
        fail(context + ": '" + trimmed + "' is not an integer");
    }
    return value;
}

double parse_double(const std::string& text, const std::string& context) {
    const std::string trimmed = strip(text);
    if (trimmed.empty()) {
        fail(context + ": empty value");
    }
    errno = 0;
    char* end = nullptr;
    // strtod accepts exactly what Python's float() accepts for the forms that
    // occur here, including 'nan' and 'inf'.
    const double value = std::strtod(trimmed.c_str(), &end);
    if (end != trimmed.c_str() + trimmed.size()) {
        fail(context + ": '" + trimmed + "' is not a number");
    }
    return value;
}

// str() of a numpy.float32 scalar: the shortest decimal that round trips as a
// float, laid out with CPython's repr rules. No matrix in the corpus reaches
// a text writer as float32, so this is the untested branch of value_repr.
std::string float32_repr(double value) {
    const float single = static_cast<float>(value);
    if (std::isnan(single)) {
        return "nan";
    }
    if (std::isinf(single)) {
        return single < 0 ? "-inf" : "inf";
    }
    std::array<char, 64> buffer{};
    const auto [end, ec] = std::to_chars(buffer.data(), buffer.data() + buffer.size(),
                                         single);
    if (ec != std::errc()) {
        return npy::float_repr(value);
    }
    std::string text(buffer.data(), end);
    // to_chars omits the '.0' that Python keeps on a whole number.
    if (text.find('.') == std::string::npos && text.find('e') == std::string::npos &&
        text.find("inf") == std::string::npos && text.find("nan") == std::string::npos) {
        text += ".0";
    }
    return text;
}

}  // namespace

std::string value_repr(double value, const std::string& dtype) {
    switch (dtype_from_name(dtype)) {
        case DType::Integer:
            if (!std::isfinite(value)) {
                return npy::float_repr(value);
            }
            return std::to_string(static_cast<std::int64_t>(value));
        case DType::Float32:
            return float32_repr(value);
        case DType::Float64:
        default:
            return npy::float_repr(value);
    }
}

CsrMatrix maximum_with_transpose(const CsrMatrix& matrix) {
    const std::int64_t rows = matrix.rows();
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& values = matrix.data();

    std::vector<std::int32_t> out_rows;
    std::vector<std::int32_t> out_cols;
    std::vector<double> out_values;
    for (std::int64_t row = 0; row < rows; ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t col = indices[k];
            if (col < row) {
                continue;  // triu(m, k=0)
            }
            const double value = values[k];
            if (col == row) {
                // max(v, v) is v, sign and all.
                out_rows.push_back(static_cast<std::int32_t>(row));
                out_cols.push_back(static_cast<std::int32_t>(col));
                out_values.push_back(value);
                continue;
            }
            // The mirror cell of the triangle is a structural zero, so the
            // elementwise maximum of the two is max(value, 0).
            const double mirrored = std::max(value, 0.0);
            if (mirrored == 0.0) {
                continue;  // scipy prunes the result of maximum
            }
            out_rows.push_back(static_cast<std::int32_t>(row));
            out_cols.push_back(static_cast<std::int32_t>(col));
            out_values.push_back(mirrored);
            out_rows.push_back(static_cast<std::int32_t>(col));
            out_cols.push_back(static_cast<std::int32_t>(row));
            out_values.push_back(mirrored);
        }
    }
    CsrMatrix result = CsrMatrix::from_coo(rows, matrix.cols(), out_rows, out_cols,
                                           std::move(out_values), matrix.dtype());
    return result;
}

CsrMatrix upper_triangle_after_maximum(const CsrMatrix& matrix) {
    const std::int64_t rows = matrix.rows();
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& values = matrix.data();

    std::vector<std::int64_t> out_indptr(static_cast<std::size_t>(rows) + 1, 0);
    std::vector<std::int32_t> out_indices;
    std::vector<double> out_values;
    for (std::int64_t row = 0; row < rows; ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t col = indices[k];
            if (col < row) {
                continue;
            }
            const double value = col == row ? values[k] : std::max(values[k], 0.0);
            if (col != row && value == 0.0) {
                continue;
            }
            out_indices.push_back(indices[k]);
            out_values.push_back(value);
        }
        out_indptr[static_cast<std::size_t>(row) + 1] =
            static_cast<std::int64_t>(out_indices.size());
    }
    // The input rows are sorted by column, so the selection is too.
    return CsrMatrix(rows, matrix.cols(), std::move(out_indptr), std::move(out_indices),
                     std::move(out_values), matrix.dtype());
}

// --------------------------------------------------------------------------
// homer

MatrixData read_homer(const std::string& path) {
    LineReader reader(path);
    std::string line;
    if (!reader.next(line)) {
        fail(path + ": the homer file is empty");
    }
    const std::vector<std::string> header = split(strip(line), '\t');
    if (header.size() < 4) {
        // values[2] and values[3] are read unconditionally to derive the bin
        // size, so fewer than two bins is not a readable file.
        fail(path + ": the homer header names fewer than two bins");
    }

    const auto split_name = [&](const std::string& name) {
        const std::vector<std::string> parts = split(strip(name), '-');
        if (parts.size() != 2) {
            fail(path + ": the bin name '" + strip(name) +
                 "' does not split into a chromosome and a start position");
        }
        return std::pair<std::string, std::int64_t>{
            parts[0], parse_int(parts[1], path + ": bin name '" + strip(name) + "'")};
    };

    const std::int64_t bin_size =
        split_name(header[3]).second - split_name(header[2]).second;

    MatrixData data;
    data.cut_intervals.reserve(header.size() - 2);
    for (std::size_t i = 2; i < header.size(); ++i) {
        const auto [chrom, start] = split_name(header[i]);
        // The Python stores the coverage as the integer 1 here, not 1.0.
        data.cut_intervals.push_back(CutInterval{chrom, start, start + bin_size, 1.0, ""});
    }

    const std::int64_t bins = static_cast<std::int64_t>(data.cut_intervals.size());
    std::vector<std::int32_t> rows;
    std::vector<std::int32_t> cols;
    std::vector<double> values;
    std::int64_t row = 0;
    while (reader.next(line)) {
        const std::vector<std::string> fields = split(line, '\t');
        if (static_cast<std::int64_t>(fields.size()) != bins + 2) {
            fail(path + ": row " + std::to_string(row + 1) + " has " +
                 std::to_string(fields.size() >= 2 ? fields.size() - 2 : 0) +
                 " values, the header names " + std::to_string(bins) + " bins");
        }
        if (row >= bins) {
            fail(path + ": more data rows than the " + std::to_string(bins) +
                 " bins the header names");
        }
        for (std::int64_t column = 0; column < bins; ++column) {
            const double value = parse_double(fields[static_cast<std::size_t>(column) + 2],
                                              path + ": row " + std::to_string(row + 1));
            if (value == 0.0) {
                continue;  // csr_matrix(dense) stores only the nonzeros
            }
            rows.push_back(static_cast<std::int32_t>(row));
            cols.push_back(static_cast<std::int32_t>(column));
            values.push_back(value);
        }
        ++row;
    }
    if (row != bins) {
        fail(path + ": " + std::to_string(row) + " data rows for " +
             std::to_string(bins) + " bins");
    }

    data.matrix = CsrMatrix::from_coo(bins, bins, rows, cols, std::move(values),
                                      "float64");
    return data;
}

void write_homer(const std::string& path, const MatrixData& data) {
    const CsrMatrix& matrix = data.matrix;
    const std::int64_t rows = matrix.rows();
    if (static_cast<std::int64_t>(data.cut_intervals.size()) != rows) {
        fail("the bin table has " + std::to_string(data.cut_intervals.size()) +
             " entries but the matrix has " + std::to_string(rows) + " rows");
    }

    std::vector<std::string> names;
    names.reserve(data.cut_intervals.size());
    for (const CutInterval& bin : data.cut_intervals) {
        names.push_back(bin.chrom + "-" + std::to_string(bin.start));
    }

    TextWriter out(path, true);
    out.write(std::string("HiCMatrix (directory=.)\tRegions\t"));
    for (const std::string& name : names) {
        out.write(name);
        out.write('\t');
    }
    out.write('\n');

    const std::string zero = value_repr(0.0, matrix.dtype());
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& values = matrix.data();
    std::vector<double> dense(static_cast<std::size_t>(matrix.cols()), 0.0);
    for (std::int64_t row = 0; row < rows; ++row) {
        std::fill(dense.begin(), dense.end(), 0.0);
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            dense[static_cast<std::size_t>(indices[k])] = values[k];
        }
        out.write(names[static_cast<std::size_t>(row)]);
        out.write('\t');
        out.write(names[static_cast<std::size_t>(row)]);
        out.write('\t');
        for (std::size_t column = 0; column < dense.size(); ++column) {
            if (column > 0) {
                out.write('\t');
            }
            out.write(dense[column] == 0.0 ? zero
                                           : value_repr(dense[column], matrix.dtype()));
        }
        // Homer.save writes the newline before the next row, so the file has
        // no trailing newline.
        if (row < rows - 1) {
            out.write('\n');
        }
    }
    out.flush();
}

// --------------------------------------------------------------------------
// ginteractions

void write_ginteractions(const std::string& path, const MatrixData& data) {
    const CsrMatrix& matrix = data.matrix;
    if (static_cast<std::int64_t>(data.cut_intervals.size()) != matrix.rows()) {
        fail("the bin table has " + std::to_string(data.cut_intervals.size()) +
             " entries but the matrix has " + std::to_string(matrix.rows()) + " rows");
    }
    TextWriter out(path + ".tsv", false);
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& values = matrix.data();
    for (std::int64_t row = 0; row < matrix.rows(); ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        const CutInterval& first = data.cut_intervals[static_cast<std::size_t>(row)];
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t col = indices[k];
            if (col < row) {
                continue;  // triu(matrix, k=0)
            }
            const CutInterval& second = data.cut_intervals[static_cast<std::size_t>(col)];
            out.write(first.chrom);
            out.write('\t');
            out.write(std::to_string(first.start));
            out.write('\t');
            out.write(std::to_string(first.end));
            out.write('\t');
            out.write(second.chrom);
            out.write('\t');
            out.write(std::to_string(second.start));
            out.write('\t');
            out.write(std::to_string(second.end));
            out.write('\t');
            out.write(value_repr(values[k], matrix.dtype()));
            out.write('\n');
        }
    }
    out.flush();
}

// --------------------------------------------------------------------------
// hicpro

MatrixData read_hicpro(const std::string& matrix_path, const std::string& bed_path) {
    MatrixData data;
    {
        LineReader bed(bed_path);
        std::string line;
        while (bed.next(line)) {
            const std::string trimmed = strip(line);
            const std::vector<std::string> fields = split(trimmed, '\t');
            if (fields.size() != 4) {
                fail(bed_path + ": expected four columns, found " +
                     std::to_string(fields.size()));
            }
            const double coverage = static_cast<double>(
                parse_int(fields[3], bed_path + ": fourth column"));
            data.cut_intervals.push_back(
                CutInterval{fields[0], parse_int(fields[1], bed_path + ": start"),
                            parse_int(fields[2], bed_path + ": end"), coverage, ""});
        }
    }
    const std::int64_t bins = static_cast<std::int64_t>(data.cut_intervals.size());
    if (bins == 0) {
        fail(bed_path + ": the bed file is empty");
    }

    std::vector<std::int32_t> rows;
    std::vector<std::int32_t> cols;
    std::vector<double> values;
    {
        LineReader matrix(matrix_path);
        std::string line;
        while (matrix.next(line)) {
            const std::string trimmed = strip(line);
            const std::vector<std::string> fields = split(trimmed, '\t');
            if (fields.size() != 3) {
                fail(matrix_path + ": expected three columns, found " +
                     std::to_string(fields.size()));
            }
            const std::int64_t row = parse_int(fields[0], matrix_path + ": first column") - 1;
            const std::int64_t col = parse_int(fields[1], matrix_path + ": second column") - 1;
            if (row < 0 || row >= bins || col < 0 || col >= bins) {
                fail(matrix_path + ": bin pair (" + fields[0] + ", " + fields[1] +
                     ") is outside the " + std::to_string(bins) + " bins of " + bed_path);
            }
            rows.push_back(static_cast<std::int32_t>(row));
            cols.push_back(static_cast<std::int32_t>(col));
            values.push_back(parse_double(fields[2], matrix_path + ": third column"));
        }
    }
    // csr_matrix((data, (i, j))) sums duplicate coordinates.
    data.matrix = CsrMatrix::from_coo(bins, bins, rows, cols, std::move(values),
                                      "float64");
    return data;
}

void write_hicpro(const std::string& matrix_path, const std::string& bed_path,
                  const MatrixData& data) {
    CsrMatrix matrix = data.matrix;
    matrix.eliminate_zeros();
    if (static_cast<std::int64_t>(data.cut_intervals.size()) != matrix.rows()) {
        fail("the bin table has " + std::to_string(data.cut_intervals.size()) +
             " entries but the matrix has " + std::to_string(matrix.rows()) + " rows");
    }
    {
        TextWriter out(matrix_path, false);
        const std::vector<std::int64_t>& indptr = matrix.indptr();
        const std::vector<std::int32_t>& indices = matrix.indices();
        const std::vector<double>& values = matrix.data();
        for (std::int64_t row = 0; row < matrix.rows(); ++row) {
            const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
            const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
            for (std::size_t k = begin; k < end; ++k) {
                // Hicpro.save writes every stored entry, both triangles
                // included: it ignores pSymmetric.
                out.write(std::to_string(row + 1));
                out.write('\t');
                out.write(std::to_string(indices[k] + 1));
                out.write('\t');
                out.write(value_repr(values[k], matrix.dtype()));
                out.write('\n');
            }
        }
        out.flush();
    }
    {
        TextWriter out(bed_path, false);
        for (std::size_t i = 0; i < data.cut_intervals.size(); ++i) {
            const CutInterval& bin = data.cut_intervals[i];
            out.write(bin.chrom);
            out.write('\t');
            out.write(std::to_string(bin.start));
            out.write('\t');
            out.write(std::to_string(bin.end));
            out.write('\t');
            out.write(std::to_string(i + 1));
            out.write('\n');
        }
        out.flush();
    }
}

// --------------------------------------------------------------------------
// 2D text

std::vector<std::pair<std::string, std::int64_t>> read_chromosome_sizes(
    const std::string& path) {
    LineReader reader(path);
    std::string line;
    std::vector<std::pair<std::string, std::int64_t>> sizes;
    std::unordered_map<std::string, std::size_t> seen;
    while (reader.next(line)) {
        const std::string trimmed = strip(line);
        if (trimmed.empty()) {
            // The Python loop stops at the first empty line, which for a well
            // formed file is the end of it.
            break;
        }
        const std::vector<std::string> fields = split(trimmed, '\t');
        if (fields.size() < 2) {
            fail(path + ": expected 'name<TAB>size', found '" + trimmed + "'");
        }
        const std::int64_t size = parse_int(fields[1], path + ": size of " + fields[0]);
        const auto [it, inserted] = seen.emplace(fields[0], sizes.size());
        if (inserted) {
            sizes.emplace_back(fields[0], size);
        } else {
            sizes[it->second].second = size;
        }
    }
    return sizes;
}

MatrixData read_two_dimensional_text(
    const std::string& path,
    const std::vector<std::pair<std::string, std::int64_t>>& chromosome_sizes,
    std::int64_t resolution) {
    if (resolution < 1) {
        fail("the resolution must be positive");
    }
    MatrixData data;
    for (const auto& [chrom, size] : chromosome_sizes) {
        for (std::int64_t start = 0; start < size; start += resolution) {
            data.cut_intervals.push_back(
                CutInterval{chrom, start, std::min(size, start + resolution), 1.0, ""});
        }
    }
    const std::int64_t bins = static_cast<std::int64_t>(data.cut_intervals.size());
    if (bins == 0) {
        fail(path + ": the chromosome sizes yield no bins");
    }

    // The Python builds a hiCMatrix over these bins and queries its interval
    // trees. The bins are generated here at a fixed resolution, so the tree
    // answer is arithmetic: bin = first_bin_of(chrom) + position / resolution
    // when 0 <= position < size, and no bin otherwise, which is exactly the
    // IndexError that makes getRegionBinRange return None. Doing the
    // arithmetic instead of building a BinTable avoids a second copy of the
    // bin table and the per bin interval records; on the 313,762 bin hg19
    // input that is 35 MB of the peak.
    std::unordered_map<std::string, std::pair<std::int64_t, std::int64_t>> offsets;
    {
        std::int64_t first = 0;
        for (const auto& [chrom, size] : chromosome_sizes) {
            offsets.emplace(chrom, std::pair<std::int64_t, std::int64_t>{first, size});
            first += (size + resolution - 1) / resolution;
        }
    }
    const auto bin_at = [&](const std::string& chrom,
                            std::int64_t position) -> std::optional<std::int64_t> {
        const auto found = offsets.find(chrom);
        if (found == offsets.end()) {
            // hicConvertFormat.py:210 calls getRegionBinRange outside its own
            // try block, so an unknown chromosome raises there. Reporting it
            // is the closest honest equivalent.
            fail(path + ": chromosome '" + chrom +
                 "' is not in the chromosome sizes file");
        }
        if (position < 0 || position >= found->second.second) {
            return std::nullopt;
        }
        return found->second.first + position / resolution;
    };

    // The Python fills a lil_matrix by assignment, so a cell keeps the last
    // value written to it and a value of exactly zero removes it again. The
    // cells are collected in file order here and resolved afterwards, which
    // needs one entry per assignment instead of a dense row index.
    std::vector<std::pair<std::int64_t, double>> assignments;
    LineReader reader(path);
    std::string line;
    std::int64_t line_number = 0;
    while (reader.next(line)) {
        ++line_number;
        if (line.empty()) {
            continue;
        }
        const std::vector<std::string> fields = split(line, '\t');
        if (fields.size() < 7) {
            fail(path + ": line " + std::to_string(line_number) + " has " +
                 std::to_string(fields.size()) + " columns, seven are needed");
        }
        const std::string where = path + ": line " + std::to_string(line_number);
        const std::int64_t start1 = parse_int(fields[1], where);
        const std::int64_t end1 = parse_int(fields[2], where);
        const std::int64_t start2 = parse_int(fields[4], where);
        const std::int64_t end2 = parse_int(fields[5], where);
        const double value = parse_double(fields[6], where);

        // getRegionBinRange returns the pair (bin of startpos, bin of endpos).
        const std::optional<std::int64_t> row_start = bin_at(fields[0], start1);
        const std::optional<std::int64_t> row_end = bin_at(fields[0], end1);
        const std::optional<std::int64_t> col_start = bin_at(fields[3], start2);
        const std::optional<std::int64_t> col_end = bin_at(fields[3], end2);
        if (!row_start.has_value() || !row_end.has_value() ||
            !col_start.has_value() || !col_end.has_value()) {
            // getRegionBinRange returned None, the assignment raised and the
            // Python swallowed the exception (hicConvertFormat.py:212-215).
            continue;
        }
        // matrix[(a, b), (c, d)] = value assigns to (a, c) and (b, d): two
        // cells whenever the interval straddles a bin boundary.
        assignments.emplace_back((*row_start << 32) | *col_start, value);
        assignments.emplace_back((*row_end << 32) | *col_end, value);
    }

    std::vector<std::size_t> order(assignments.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return assignments[a].first < assignments[b].first;
    });
    std::vector<std::int32_t> rows;
    std::vector<std::int32_t> cols;
    std::vector<double> values;
    for (std::size_t position = 0; position < order.size();) {
        const std::int64_t key = assignments[order[position]].first;
        std::size_t last = position;
        while (position < order.size() && assignments[order[position]].first == key) {
            last = position;
            ++position;
        }
        const double value = assignments[order[last]].second;
        if (value == 0.0) {
            continue;
        }
        rows.push_back(static_cast<std::int32_t>(key >> 32));
        cols.push_back(static_cast<std::int32_t>(
            static_cast<std::uint32_t>(key & 0xFFFFFFFF)));
        values.push_back(value);
    }

    data.matrix = CsrMatrix::from_coo(bins, bins, rows, cols, std::move(values),
                                      "float64");
    return data;
}

}  // namespace hicx
