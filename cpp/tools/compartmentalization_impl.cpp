#include "compartmentalization_impl.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <unordered_set>

#include "hicx/numpy_compat.hpp"

namespace hicx::compartments {

namespace {

std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true) {
        const std::size_t tab = line.find('\t', begin);
        if (tab == std::string::npos) {
            fields.push_back(line.substr(begin));
            break;
        }
        fields.push_back(line.substr(begin, tab - begin));
        begin = tab + 1;
    }
    return fields;
}

// pandas' default na_values (pandas/_libs/parsers.pyx STR_NA_VALUES).
bool is_pandas_na(const std::string& text) {
    static const std::unordered_set<std::string> kNa = {
        "",     "#N/A", "#N/A N/A", "#NA",  "-1.#IND", "-1.#QNAN", "-NaN", "-nan",
        "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA",      "NULL",     "NaN",  "None",
        "n/a",  "nan",  "null"};
    return kNa.count(text) != 0;
}

std::int64_t parse_int64(const std::string& text, const std::string& path, int line) {
    errno = 0;
    char* end = nullptr;
    const long long value = std::strtoll(text.c_str(), &end, 10);
    if (text.empty() || end != text.c_str() + text.size() || errno != 0) {
        throw std::runtime_error(path + ":" + std::to_string(line) + ": '" + text +
                                 "' is not an integer coordinate");
    }
    return static_cast<std::int64_t>(value);
}

// numpy's ordering for floating point sort and search: NaN after everything.
bool npy_less(double a, double b) { return a < b || (b != b && a == a); }

}  // namespace

std::vector<PcaRow> read_pca_bedgraph(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("[Errno 2] No such file or directory: '" + path + "'");
    }
    std::vector<PcaRow> rows;
    std::string line;
    int line_number = 0;
    std::size_t expected_fields = 0;
    while (std::getline(file, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;  // skip_blank_lines
        }
        const std::vector<std::string> fields = split_tabs(line);
        if (expected_fields == 0) {
            expected_fields = fields.size();
        } else if (fields.size() > expected_fields) {
            throw std::runtime_error("Error tokenizing data. C error: Expected " +
                                     std::to_string(expected_fields) + " fields in line " +
                                     std::to_string(line_number) + ", saw " +
                                     std::to_string(fields.size()));
        }
        if (fields.size() < 3) {
            throw std::runtime_error(path + ":" + std::to_string(line_number) +
                                     ": fewer than three columns");
        }
        if (expected_fields < 4) {
            // pc1.rename leaves no column 3, and the Python fails on pc1['pc1'].
            throw std::runtime_error("KeyError: 'pc1' (the bedgraph has fewer than four "
                                     "columns)");
        }
        PcaRow row;
        row.chrom = fields[0];
        row.start = parse_int64(fields[1], path, line_number);
        row.end = parse_int64(fields[2], path, line_number);
        const std::string value_text = fields.size() > 3 ? fields[3] : std::string();
        if (is_pandas_na(value_text)) {
            row.pc1 = std::numeric_limits<float>::quiet_NaN();
        } else {
            errno = 0;
            char* end = nullptr;
            const double value = std::strtod(value_text.c_str(), &end);
            if (end != value_text.c_str() + value_text.size()) {
                throw std::runtime_error(path + ":" + std::to_string(line_number) +
                                         ": could not convert string to float: '" +
                                         value_text + "'");
            }
            row.pc1 = static_cast<float>(value);
        }
        rows.push_back(std::move(row));
    }
    return rows;
}

std::vector<double> nanquantile_linear(std::vector<double> values,
                                       const std::vector<double>& quantiles) {
    for (const double q : quantiles) {
        if (!(q >= 0.0 && q <= 1.0)) {
            throw std::invalid_argument("Quantiles must be in the range [0, 1]");
        }
    }
    values.erase(std::remove_if(values.begin(), values.end(),
                                [](double v) { return std::isnan(v); }),
                 values.end());
    std::vector<double> result(quantiles.size(), std::numeric_limits<double>::quiet_NaN());
    const auto n = static_cast<std::int64_t>(values.size());
    if (n == 0) {
        return result;  // numpy warns "All-NaN slice encountered" and returns NaN
    }
    std::sort(values.begin(), values.end());
    for (std::size_t k = 0; k < quantiles.size(); ++k) {
        // _QuantileMethods['linear'].get_virtual_index: (n - 1) * quantiles
        const double virtual_index = static_cast<double>(n - 1) * quantiles[k];
        // _get_indexes
        double previous = std::floor(virtual_index);
        double next = previous + 1.0;
        if (virtual_index >= static_cast<double>(n - 1)) {
            previous = -1.0;
            next = -1.0;
        }
        if (virtual_index < 0.0) {
            previous = 0.0;
            next = 0.0;
        }
        if (std::isnan(virtual_index)) {
            previous = -1.0;
            next = -1.0;
        }
        const auto previous_index = static_cast<std::int64_t>(previous);
        const auto next_index = static_cast<std::int64_t>(next);
        // _get_gamma: virtual_indexes - previous_indexes, the latter as intp.
        const double gamma = virtual_index - static_cast<double>(previous_index);
        const double a =
            values[static_cast<std::size_t>(previous_index < 0 ? n + previous_index
                                                               : previous_index)];
        const double b =
            values[static_cast<std::size_t>(next_index < 0 ? n + next_index : next_index)];
        // _lerp
        const double diff_b_a = b - a;
        double lerp = a + diff_b_a * gamma;
        if (gamma >= 0.5) {
            lerp = b - diff_b_a * (1.0 - gamma);
        }
        result[k] = lerp;
    }
    return result;
}

std::vector<double> linspace(double start, double stop, std::int64_t num) {
    if (num < 0) {
        throw std::invalid_argument("Number of samples, " + std::to_string(num) +
                                    ", must be non-negative.");
    }
    const std::int64_t div = num - 1;
    const double delta = stop - start;
    std::vector<double> y(static_cast<std::size_t>(num));
    for (std::int64_t i = 0; i < num; ++i) {
        y[static_cast<std::size_t>(i)] = static_cast<double>(i);
    }
    if (div > 0) {
        const double step = delta / static_cast<double>(div);
        if (step == 0.0) {
            for (double& value : y) {
                value /= static_cast<double>(div);
                value *= delta;
            }
        } else {
            for (double& value : y) {
                value *= step;
            }
        }
    } else {
        for (double& value : y) {
            value = value * delta;
        }
    }
    for (double& value : y) {
        value += start;
    }
    if (num > 1) {
        y.back() = stop;
    }
    return y;
}

std::vector<std::int64_t> searchsorted_right(const std::vector<double>& sorted,
                                             const std::vector<double>& keys) {
    // numpy/_core/src/npysort/binsearch.cpp, side right: cmp(a, b) is
    // !less(b, a).
    const auto cmp = [](double a, double b) { return !npy_less(b, a); };
    std::vector<std::int64_t> result(keys.size());
    if (keys.empty()) {
        return result;
    }
    const auto arr_len = static_cast<std::int64_t>(sorted.size());
    std::int64_t min_idx = 0;
    std::int64_t max_idx = arr_len;
    double last_key = keys[0];
    for (std::size_t k = 0; k < keys.size(); ++k) {
        const double key = keys[k];
        if (cmp(last_key, key)) {
            max_idx = arr_len;
        } else {
            min_idx = 0;
            max_idx = (max_idx < arr_len) ? (max_idx + 1) : arr_len;
        }
        last_key = key;
        while (min_idx < max_idx) {
            const std::int64_t mid = min_idx + ((max_idx - min_idx) >> 1);
            if (cmp(sorted[static_cast<std::size_t>(mid)], key)) {
                min_idx = mid + 1;
            } else {
                max_idx = mid;
            }
        }
        result[k] = min_idx;
    }
    return result;
}

std::vector<double> quantile_boundaries(const std::vector<PcaRow>& rows,
                                        std::int64_t quantiles, double outliers) {
    std::vector<double> values;
    values.reserve(rows.size());
    for (const auto& row : rows) {
        values.push_back(static_cast<double>(row.pc1));
    }
    if (outliers != 0) {
        const std::vector<double> quantile = {outliers / 100, (100 - outliers) / 100};
        const std::vector<double> boundaries = nanquantile_linear(values, quantile);
        return linspace(boundaries[0], boundaries[1], quantiles);
    }
    if (quantiles == 1) {
        // [j / (args.quantile - 1) for j in range(0, args.quantile)]
        throw std::domain_error("ZeroDivisionError: division by zero");
    }
    std::vector<double> quantile;
    for (std::int64_t j = 0; j < quantiles; ++j) {
        quantile.push_back(static_cast<double>(j) / static_cast<double>(quantiles - 1));
    }
    return nanquantile_linear(values, quantile);
}

SymmetricRows::SymmetricRows(const CsrMatrix& matrix)
    : matrix_(&matrix),
      n_(matrix.rows()),
      upper_(matrix.symmetry() == Symmetry::UpperTriangle) {
    if (matrix.rows() != matrix.cols()) {
        throw std::runtime_error("the matrix is not square");
    }
    if (!upper_) {
        return;
    }
    const auto& indptr = matrix.indptr();
    const auto& indices = matrix.indices();
    column_indptr_.assign(static_cast<std::size_t>(n_) + 1, 0);
    for (std::int64_t k = 0; k < n_; ++k) {
        for (auto p = indptr[static_cast<std::size_t>(k)];
             p < indptr[static_cast<std::size_t>(k) + 1]; ++p) {
            const std::int32_t c = indices[static_cast<std::size_t>(p)];
            if (c > k) {
                ++column_indptr_[static_cast<std::size_t>(c) + 1];
            }
        }
    }
    for (std::size_t c = 0; c < static_cast<std::size_t>(n_); ++c) {
        column_indptr_[c + 1] += column_indptr_[c];
    }
    column_rows_.resize(static_cast<std::size_t>(column_indptr_.back()));
    std::vector<std::int64_t> fill(column_indptr_.begin(), column_indptr_.end() - 1);
    for (std::int64_t k = 0; k < n_; ++k) {
        for (auto p = indptr[static_cast<std::size_t>(k)];
             p < indptr[static_cast<std::size_t>(k) + 1]; ++p) {
            const std::int32_t c = indices[static_cast<std::size_t>(p)];
            if (c > k) {
                column_rows_[static_cast<std::size_t>(fill[static_cast<std::size_t>(c)]++)] =
                    static_cast<std::int32_t>(k);
            }
        }
    }
}

std::size_t SymmetricRows::position_in_row(std::int32_t row, std::int64_t column) const {
    const auto& indptr = matrix_->indptr();
    const auto& indices = matrix_->indices();
    const auto first = indices.begin() + indptr[static_cast<std::size_t>(row)];
    const auto last = indices.begin() + indptr[static_cast<std::size_t>(row) + 1];
    const auto found = std::lower_bound(first, last, static_cast<std::int32_t>(column));
    return static_cast<std::size_t>(found - indices.begin());
}

namespace {

// np.sum over a contiguous float64 array, fed one value at a time. numpy
// reduces through its 8192 element ufunc buffer and adds the pairwise sums of
// consecutive buffers sequentially, starting from zero, so filling a buffer of
// that size and handing each full one to npy::pairwise_sum is the same
// arithmetic without materialising the array.
class PairwiseStream {
  public:
    static constexpr std::size_t kBuffer = 8192;

    void push(double value) {
        buffer_[length_++] = value;
        ++count_;
        if (length_ == kBuffer) {
            accumulated_ += npy::pairwise_sum(buffer_.data(), length_);
            length_ = 0;
        }
    }
    [[nodiscard]] double finish() {
        if (length_ != 0) {
            accumulated_ += npy::pairwise_sum(buffer_.data(), length_);
            length_ = 0;
        }
        return accumulated_;
    }
    [[nodiscard]] std::int64_t count() const noexcept { return count_; }

  private:
    std::array<double, kBuffer> buffer_{};
    std::size_t length_ = 0;
    double accumulated_ = 0.0;
    std::int64_t count_ = 0;
};

}  // namespace

std::vector<double> normalised_sum_per_quantile(
    const CsrMatrix& matrix, const std::vector<std::vector<std::int64_t>>& bin_ids,
    const std::vector<std::int64_t>& quantile_of_row, std::int64_t quantiles,
    const std::vector<std::int64_t>& offsets, std::int64_t chromosome_count) {
    const auto Q = static_cast<std::size_t>(quantiles);
    const SymmetricRows rows(matrix);
    const std::int64_t n = rows.size();

    // The bins of every quantile, pca rows in file order, each row's bins in
    // ascending order: the row_indices and col_indices of :112-120. A
    // quantile with no pca row is skipped by `if row_indices.empty: continue`,
    // which is not the same as a quantile whose rows cover no bin.
    std::vector<std::vector<std::int64_t>> members(Q);
    std::vector<bool> present(Q, false);
    for (std::size_t k = 0; k < quantile_of_row.size(); ++k) {
        const std::int64_t q = quantile_of_row[k];
        if (q < 0 || q >= quantiles) {
            continue;
        }
        present[static_cast<std::size_t>(q)] = true;
        auto& list = members[static_cast<std::size_t>(q)];
        list.insert(list.end(), bin_ids[k].begin(), bin_ids[k].end());
    }

    // The sum and the finite count of every block obs_exp[np.ix_(rows, cols)],
    // in the row major order of the dense block, which is the order np.sum
    // sees after the NaN and infinity filters flatten it.
    std::vector<double> block_sum(Q * Q, 0.0);
    std::vector<std::int64_t> block_count(Q * Q, 0);
    std::vector<double> dense(static_cast<std::size_t>(n), 0.0);
    std::vector<std::int64_t> touched;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (std::size_t qi = 0; qi < Q; ++qi) {
        if (!present[qi]) {
            continue;
        }
        std::vector<PairwiseStream> streams(Q);
        for (const std::int64_t r : members[qi]) {
            rows.for_each_in_row(r, [&](std::int64_t column, double value) {
                dense[static_cast<std::size_t>(column)] = value;
                touched.push_back(column);
            });
            // --offset writes NaN into both triangles of the named diagonals
            // (:100-105) before any block is read.
            for (const std::int64_t dist : offsets) {
                if (dist >= n) {
                    continue;
                }
                if (r + dist < n) {
                    dense[static_cast<std::size_t>(r + dist)] = nan;
                    touched.push_back(r + dist);
                }
                if (r - dist >= 0) {
                    dense[static_cast<std::size_t>(r - dist)] = nan;
                    touched.push_back(r - dist);
                }
            }
            for (std::size_t qj = 0; qj < Q; ++qj) {
                if (!present[qj]) {
                    continue;
                }
                PairwiseStream& stream = streams[qj];
                for (const std::int64_t c : members[qj]) {
                    const double value = dense[static_cast<std::size_t>(c)];
                    if (std::isfinite(value)) {
                        stream.push(value);
                    }
                }
            }
            for (const std::int64_t column : touched) {
                dense[static_cast<std::size_t>(column)] = 0.0;
            }
            touched.clear();
        }
        for (std::size_t qj = 0; qj < Q; ++qj) {
            if (present[qj]) {
                block_sum[qi * Q + qj] = streams[qj].finish();
                block_count[qi * Q + qj] = streams[qj].count();
            }
        }
    }

    // The accumulation of :125-128, repeated once per chromosome of the pca
    // file because the loop at :106 does not use its variable. Every
    // repetition adds the same block values in the same order, so replaying
    // the stored sums is the same floating point sequence without recomputing
    // the blocks.
    std::vector<double> interaction_sum(Q * Q, 0.0);
    std::vector<double> number_of_bins(Q * Q, 0.0);
    for (std::int64_t chrom = 0; chrom < chromosome_count; ++chrom) {
        for (std::size_t qi = 0; qi < Q; ++qi) {
            if (!present[qi]) {
                continue;
            }
            for (std::size_t qj = 0; qj < Q; ++qj) {
                if (!present[qj]) {
                    continue;
                }
                const double sum = block_sum[qi * Q + qj];
                const auto count = static_cast<double>(block_count[qi * Q + qj]);
                interaction_sum[qi * Q + qj] += sum;
                interaction_sum[qj * Q + qi] += sum;
                number_of_bins[qi * Q + qj] += count;
                number_of_bins[qj * Q + qi] += count;
            }
        }
    }
    std::vector<double> result(Q * Q);
    for (std::size_t k = 0; k < Q * Q; ++k) {
        result[k] = interaction_sum[k] / number_of_bins[k];
    }
    nan_to_num_in_place(result);
    return result;
}

void nan_to_num_in_place(std::vector<double>& values) {
    for (double& value : values) {
        if (std::isnan(value)) {
            value = 0.0;
        } else if (std::isinf(value)) {
            value = value > 0 ? std::numeric_limits<double>::max()
                              : std::numeric_limits<double>::lowest();
        }
    }
}

namespace {

double slice_sum(const std::vector<double>& matrix, std::size_t quantiles,
                 std::size_t row_begin, std::size_t row_end, std::size_t col_begin,
                 std::size_t col_end) {
    // numpy reduces a non contiguous slice through a contiguous buffered copy,
    // measured to agree bit for bit with the pairwise sum of that copy.
    std::vector<double> copy;
    copy.reserve((row_end - row_begin) * (col_end - col_begin));
    for (std::size_t i = row_begin; i < row_end; ++i) {
        for (std::size_t j = col_begin; j < col_end; ++j) {
            copy.push_back(matrix[i * quantiles + j]);
        }
    }
    return npy::pairwise_sum(copy);
}

}  // namespace

std::vector<double> within_vs_between(const std::vector<double>& normalised,
                                      std::int64_t quantiles) {
    const auto Q = static_cast<std::size_t>(quantiles);
    std::vector<double> ratios;
    for (std::size_t q = 1; q < Q; ++q) {
        const double within = slice_sum(normalised, Q, 0, q, 0, q) +
                              slice_sum(normalised, Q, Q - q, Q, Q - q, Q);
        const double between = slice_sum(normalised, Q, 0, q, Q - q, Q) +
                               slice_sum(normalised, Q, Q - q, Q, 0, q);
        ratios.push_back(within / between);
    }
    return ratios;
}

std::string savetxt_line(const std::vector<double>& row) {
    std::string line;
    for (std::size_t k = 0; k < row.size(); ++k) {
        if (k != 0) {
            line += ' ';
        }
        const double value = row[k];
        if (std::isnan(value)) {
            line += "nan";
        } else if (std::isinf(value)) {
            line += value > 0 ? "inf" : "-inf";
        } else {
            std::array<char, 64> buffer{};
            const int written = std::snprintf(buffer.data(), buffer.size(), "%.18e", value);
            line.append(buffer.data(), static_cast<std::size_t>(written));
        }
    }
    line += '\n';
    return line;
}

}  // namespace hicx::compartments
