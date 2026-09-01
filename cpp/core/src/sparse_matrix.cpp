#include "hicx/sparse_matrix.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "hicx/numpy_compat.hpp"

namespace hicx {

CsrMatrix::CsrMatrix(std::int64_t rows, std::int64_t cols,
                     std::vector<std::int64_t> indptr,
                     std::vector<std::int32_t> indices, std::vector<double> data,
                     std::string dtype)
    : rows_(rows),
      cols_(cols),
      indptr_(std::move(indptr)),
      indices_(std::move(indices)),
      data_(std::move(data)),
      dtype_(std::move(dtype)) {
    if (static_cast<std::int64_t>(indptr_.size()) != rows_ + 1) {
        throw std::invalid_argument("indptr length does not match the row count");
    }
    if (indices_.size() != data_.size()) {
        throw std::invalid_argument("indices and data have different lengths");
    }
}

DType dtype_from_name(const std::string& name) {
    if (name == "float32") {
        return DType::Float32;
    }
    if (name.rfind("int", 0) == 0 || name.rfind("uint", 0) == 0 || name == "bool") {
        return DType::Integer;
    }
    return DType::Float64;
}

CsrMatrix CsrMatrix::from_coo(std::int64_t rows, std::int64_t cols,
                              const std::vector<std::int32_t>& row,
                              const std::vector<std::int32_t>& col,
                              std::vector<double> data, std::string dtype) {
    if (row.size() != col.size() || row.size() != data.size()) {
        throw std::invalid_argument("coordinate arrays have different lengths");
    }
    const std::size_t n = row.size();

    std::vector<std::int64_t> indptr(static_cast<std::size_t>(rows) + 1, 0);
    for (std::size_t k = 0; k < n; ++k) {
        ++indptr[static_cast<std::size_t>(row[k]) + 1];
    }
    for (std::size_t i = 1; i < indptr.size(); ++i) {
        indptr[i] += indptr[i - 1];
    }

    std::vector<std::int32_t> indices(n);
    std::vector<double> values(n);
    std::vector<std::int64_t> cursor(indptr.begin(), indptr.end() - 1);
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t position = static_cast<std::size_t>(cursor[static_cast<std::size_t>(row[k])]++);
        indices[position] = col[k];
        values[position] = data[k];
    }

    // Sort every row by column index and sum duplicates, like coo_tocsr
    // followed by sum_duplicates.
    std::vector<std::size_t> order;
    std::vector<std::int32_t> sorted_indices;
    std::vector<double> sorted_values;
    sorted_indices.reserve(n);
    sorted_values.reserve(n);
    std::vector<std::int64_t> new_indptr(indptr.size(), 0);

    for (std::int64_t i = 0; i < rows; ++i) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(i)]);
        const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(i) + 1]);
        order.resize(end - begin);
        for (std::size_t k = 0; k < order.size(); ++k) {
            order[k] = begin + k;
        }
        std::stable_sort(order.begin(), order.end(),
                         [&](std::size_t a, std::size_t b) {
                             return indices[a] < indices[b];
                         });
        for (std::size_t k = 0; k < order.size();) {
            const std::int32_t column = indices[order[k]];
            double accumulated = values[order[k]];
            ++k;
            while (k < order.size() && indices[order[k]] == column) {
                accumulated += values[order[k]];
                ++k;
            }
            sorted_indices.push_back(column);
            sorted_values.push_back(accumulated);
        }
        new_indptr[static_cast<std::size_t>(i) + 1] =
            static_cast<std::int64_t>(sorted_indices.size());
    }

    return CsrMatrix(rows, cols, std::move(new_indptr), std::move(sorted_indices),
                     std::move(sorted_values), std::move(dtype));
}

void CsrMatrix::eliminate_zeros() {
    std::vector<std::int64_t> new_indptr(indptr_.size(), 0);
    std::size_t write = 0;
    for (std::int64_t i = 0; i < rows_; ++i) {
        const std::size_t begin = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i)]);
        const std::size_t end = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            if (data_[k] != 0.0) {
                indices_[write] = indices_[k];
                data_[write] = data_[k];
                ++write;
            }
        }
        new_indptr[static_cast<std::size_t>(i) + 1] = static_cast<std::int64_t>(write);
    }
    indices_.resize(write);
    data_.resize(write);
    indptr_ = std::move(new_indptr);
}

std::size_t CsrMatrix::nnz() const {
    if (symmetry_ == Symmetry::Full) {
        return data_.size();
    }
    std::size_t count = 0;
    for (std::int64_t i = 0; i < rows_; ++i) {
        const std::size_t begin = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i)]);
        const std::size_t end = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            if (data_[k] == 0.0) {
                continue;  // the CSR addition drops entries that are zero
            }
            count += (indices_[k] == i) ? 1 : 2;
        }
    }
    return count;
}

bool CsrMatrix::lower_triangle_is_zero() const {
    double total = 0.0;
    for (std::int64_t i = 0; i < rows_; ++i) {
        const std::size_t begin = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i)]);
        const std::size_t end = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            if (indices_[k] < i) {
                total += data_[k];
            }
        }
    }
    return total == 0.0;
}

bool CsrMatrix::has_lower_entries() const {
    for (std::int64_t i = 0; i < rows_; ++i) {
        const std::size_t begin = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i)]);
        const std::size_t end = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            if (indices_[k] < i) {
                return true;
            }
        }
    }
    return false;
}

void CsrMatrix::symmetrize_in_place() {
    if (symmetry_ == Symmetry::UpperTriangle) {
        return;
    }
    if (!has_lower_entries()) {
        symmetry_ = Symmetry::UpperTriangle;
        return;
    }
    if (lower_triangle_is_zero()) {
        // Entries below the diagonal exist but cancel out. The Python code
        // still adds the transpose, so the cheap relabelling is not valid and
        // the full matrix has to be built.
        materialize_full();
    }
}

void CsrMatrix::materialize_full() {
    if (symmetry_ == Symmetry::Full && !lower_triangle_is_zero()) {
        return;
    }
    // result = self + triu(self, 1).T
    // Every stored entry (r, c) with c > r contributes to row c, column r.
    std::vector<std::int64_t> lower_count(static_cast<std::size_t>(rows_) + 1, 0);
    for (std::int64_t i = 0; i < rows_; ++i) {
        const std::size_t begin = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i)]);
        const std::size_t end = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            if (indices_[k] > i) {
                ++lower_count[static_cast<std::size_t>(indices_[k]) + 1];
            }
        }
    }
    std::vector<std::int64_t> lower_offset(lower_count);
    for (std::size_t i = 1; i < lower_offset.size(); ++i) {
        lower_offset[i] += lower_offset[i - 1];
    }
    const std::size_t lower_total = static_cast<std::size_t>(lower_offset.back());
    std::vector<std::int32_t> lower_indices(lower_total);
    std::vector<double> lower_data(lower_total);
    {
        std::vector<std::int64_t> cursor(lower_offset.begin(), lower_offset.end() - 1);
        // Rows are visited in ascending order, so the transposed entries land
        // in each target row already sorted by column.
        for (std::int64_t i = 0; i < rows_; ++i) {
            const std::size_t begin = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i)]);
            const std::size_t end = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i) + 1]);
            for (std::size_t k = begin; k < end; ++k) {
                if (indices_[k] > i) {
                    const std::size_t target = static_cast<std::size_t>(indices_[k]);
                    const std::size_t position = static_cast<std::size_t>(cursor[target]++);
                    lower_indices[position] = static_cast<std::int32_t>(i);
                    lower_data[position] = data_[k];
                }
            }
        }
    }

    std::vector<std::int64_t> new_indptr(indptr_.size(), 0);
    std::vector<std::int32_t> new_indices;
    std::vector<double> new_data;
    new_indices.reserve(indices_.size() + lower_total);
    new_data.reserve(data_.size() + lower_total);
    const bool single_precision = dtype_kind() == DType::Float32;

    for (std::int64_t i = 0; i < rows_; ++i) {
        std::size_t a = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i)]);
        const std::size_t a_end = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i) + 1]);
        std::size_t b = static_cast<std::size_t>(lower_offset[static_cast<std::size_t>(i)]);
        const std::size_t b_end = static_cast<std::size_t>(lower_offset[static_cast<std::size_t>(i) + 1]);
        while (a < a_end || b < b_end) {
            std::int32_t column = 0;
            double value = 0.0;
            const bool take_a = (b >= b_end) || (a < a_end && indices_[a] <= lower_indices[b]);
            const bool take_b = (a >= a_end) || (b < b_end && lower_indices[b] <= indices_[a]);
            if (take_a && take_b) {
                column = indices_[a];
                value = single_precision
                            ? static_cast<double>(static_cast<float>(data_[a]) +
                                                  static_cast<float>(lower_data[b]))
                            : data_[a] + lower_data[b];
                ++a;
                ++b;
            } else if (take_a) {
                column = indices_[a];
                value = data_[a];
                ++a;
            } else {
                column = lower_indices[b];
                value = lower_data[b];
                ++b;
            }
            // scipy's csr_binop_csr only stores non zero results.
            if (value != 0.0) {
                new_indices.push_back(column);
                new_data.push_back(value);
            }
        }
        new_indptr[static_cast<std::size_t>(i) + 1] =
            static_cast<std::int64_t>(new_indices.size());
    }

    indptr_ = std::move(new_indptr);
    indices_ = std::move(new_indices);
    data_ = std::move(new_data);
    symmetry_ = Symmetry::Full;
}

namespace {

// Row sums of the symmetric matrix represented by an upper triangle, in
// scipy's accumulation order.
//
// Row i of the full matrix holds, in ascending column order, first the mirrored
// entries taken from column i of the upper triangle (their row index ascends
// with the scan) and then the entries stored in row i itself. Scanning the
// upper triangle row by row therefore delivers every accumulator its values in
// exactly the order the sparse matrix vector product would.
template <typename T>
std::vector<T> symmetric_row_sums(std::int64_t rows,
                                  const std::vector<std::int64_t>& indptr,
                                  const std::vector<std::int32_t>& indices,
                                  const std::vector<double>& data) {
    std::vector<T> row_sums(static_cast<std::size_t>(rows), static_cast<T>(0));
    for (std::int64_t i = 0; i < rows; ++i) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(i)]);
        const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(i) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const double raw = data[k];
            if (raw == 0.0) {
                continue;  // dropped by the CSR addition
            }
            const T value = static_cast<T>(raw);
            const std::int32_t column = indices[k];
            if (column != i) {
                row_sums[static_cast<std::size_t>(column)] += value;
            }
            row_sums[static_cast<std::size_t>(i)] += value;
        }
    }
    return row_sums;
}

}  // namespace

Scalar CsrMatrix::sum() const {
    if (symmetry_ == Symmetry::UpperTriangle) {
        switch (dtype_kind()) {
            case DType::Integer: {
                const std::vector<std::int64_t> row_sums =
                    symmetric_row_sums<std::int64_t>(rows_, indptr_, indices_, data_);
                std::int64_t total = 0;
                for (const std::int64_t value : row_sums) {
                    total += value;
                }
                return Scalar::from_int(total);
            }
            case DType::Float32: {
                const std::vector<float> row_sums =
                    symmetric_row_sums<float>(rows_, indptr_, indices_, data_);
                return Scalar::from_float(npy::pairwise_sum(row_sums));
            }
            case DType::Float64:
            default: {
                const std::vector<double> row_sums =
                    symmetric_row_sums<double>(rows_, indptr_, indices_, data_);
                return Scalar::from_double(npy::pairwise_sum(row_sums));
            }
        }
    }
    // scipy multiplies by a dense vector of ones, so every row is accumulated
    // sequentially in the dtype of the matrix, and numpy then reduces the
    // dense result pairwise.
    switch (dtype_kind()) {
        case DType::Integer: {
            std::int64_t total = 0;
            for (std::int64_t i = 0; i < rows_; ++i) {
                const std::size_t begin = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i)]);
                const std::size_t end = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i) + 1]);
                std::int64_t row_sum = 0;
                for (std::size_t k = begin; k < end; ++k) {
                    row_sum += static_cast<std::int64_t>(data_[k]);
                }
                total += row_sum;
            }
            return Scalar::from_int(total);
        }
        case DType::Float32: {
            std::vector<float> row_sums(static_cast<std::size_t>(rows_), 0.0F);
            for (std::int64_t i = 0; i < rows_; ++i) {
                const std::size_t begin = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i)]);
                const std::size_t end = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i) + 1]);
                float row_sum = 0.0F;
                for (std::size_t k = begin; k < end; ++k) {
                    row_sum += static_cast<float>(data_[k]);
                }
                row_sums[static_cast<std::size_t>(i)] = row_sum;
            }
            return Scalar::from_float(npy::pairwise_sum(row_sums));
        }
        case DType::Float64:
        default: {
            std::vector<double> row_sums(static_cast<std::size_t>(rows_), 0.0);
            for (std::int64_t i = 0; i < rows_; ++i) {
                const std::size_t begin = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i)]);
                const std::size_t end = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i) + 1]);
                double row_sum = 0.0;
                for (std::size_t k = begin; k < end; ++k) {
                    row_sum += data_[k];
                }
                row_sums[static_cast<std::size_t>(i)] = row_sum;
            }
            return Scalar::from_double(npy::pairwise_sum(row_sums));
        }
    }
}

std::vector<double> CsrMatrix::diagonal() const {
    const std::int64_t length = std::min(rows_, cols_);
    std::vector<double> result(static_cast<std::size_t>(std::max<std::int64_t>(length, 0)), 0.0);
    const bool single = dtype_kind() == DType::Float32;
    for (std::int64_t i = 0; i < length; ++i) {
        const std::size_t begin = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i)]);
        const std::size_t end = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(i) + 1]);
        double value = 0.0;
        for (std::size_t k = begin; k < end; ++k) {
            if (indices_[k] == i) {
                // csr_diagonal sums duplicate entries, in the matrix dtype.
                value = single ? static_cast<double>(static_cast<float>(value) +
                                                     static_cast<float>(data_[k]))
                               : value + data_[k];
            }
        }
        result[static_cast<std::size_t>(i)] = value;
    }
    return result;
}

Scalar CsrMatrix::diagonal_sum() const {
    const std::vector<double> diag = diagonal();
    switch (dtype_kind()) {
        case DType::Integer: {
            std::int64_t total = 0;
            for (const double value : diag) {
                total += static_cast<std::int64_t>(value);
            }
            return Scalar::from_int(total);
        }
        case DType::Float32: {
            std::vector<float> single(diag.size());
            for (std::size_t i = 0; i < diag.size(); ++i) {
                single[i] = static_cast<float>(diag[i]);
            }
            return Scalar::from_float(npy::pairwise_sum(single));
        }
        case DType::Float64:
        default: return Scalar::from_double(npy::pairwise_sum(diag));
    }
}

namespace {

Scalar as_scalar(double value, DType kind) {
    switch (kind) {
        case DType::Integer: return Scalar::from_int(static_cast<std::int64_t>(value));
        case DType::Float32: return Scalar::from_float(static_cast<float>(value));
        case DType::Float64:
        default: return Scalar::from_double(value);
    }
}

}  // namespace

namespace {

// Extremum over the data array of the represented matrix. For an upper
// triangle the mirrored copies hold the same values, and the exact zeros that
// the CSR addition drops are skipped.
Scalar data_extremum(const std::vector<double>& data, DType kind, bool skip_zeros,
                     bool want_max) {
    bool seen = false;
    double best = 0.0;
    for (const double value : data) {
        if (skip_zeros && value == 0.0) {
            continue;
        }
        // numpy's minimum.reduce and maximum.reduce propagate NaN.
        if (std::isnan(value)) {
            return as_scalar(value, kind);
        }
        if (!seen || (want_max ? value > best : value < best)) {
            best = value;
            seen = true;
        }
    }
    if (!seen) {
        throw std::runtime_error("zero-size array reduction: the matrix is empty");
    }
    return as_scalar(best, kind);
}

}  // namespace

Scalar CsrMatrix::data_min() const {
    return data_extremum(data_, dtype_kind(), symmetry_ == Symmetry::UpperTriangle,
                         false);
}

Scalar CsrMatrix::data_max() const {
    return data_extremum(data_, dtype_kind(), symmetry_ == Symmetry::UpperTriangle,
                         true);
}

double CsrMatrix::at(std::int64_t row, std::int64_t col) const {
    if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
        throw std::out_of_range("matrix index out of range");
    }
    if (symmetry_ == Symmetry::UpperTriangle && row > col) {
        std::swap(row, col);
    }
    const std::size_t begin = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(row)]);
    const std::size_t end = static_cast<std::size_t>(indptr_[static_cast<std::size_t>(row) + 1]);
    double value = 0.0;
    for (std::size_t k = begin; k < end; ++k) {
        if (indices_[k] == col) {
            value += data_[k];
        }
    }
    return value;
}

}  // namespace hicx
