#include "hicx/matrix_ops.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <unordered_map>

#include "hicx/numpy_compat.hpp"

namespace hicx {

namespace {

// numpy's ufunc reduction buffer, the same constant numpy_compat uses. A
// pairwise sum over a whole array is the sequential accumulation of the
// pairwise sums of its 8192 element blocks, so the reduction can be fed one
// value at a time without holding the array.
constexpr std::size_t kReduceBufferSize = 8192;

template <typename T>
class StreamingPairwiseSum {
  public:
    StreamingPairwiseSum() { buffer_.reserve(kReduceBufferSize); }

    void push(T value) {
        buffer_.push_back(value);
        if (buffer_.size() == kReduceBufferSize) {
            flush();
        }
    }

    [[nodiscard]] T finish() {
        flush();
        return total_;
    }

  private:
    void flush() {
        if (buffer_.empty()) {
            return;
        }
        total_ += npy::pairwise_sum(buffer_.data(), buffer_.size());
        buffer_.clear();
    }

    std::vector<T> buffer_;
    T total_ = static_cast<T>(0);
};

// How the arithmetic of a binary operation is carried out, which follows the
// promoted dtype and not the storage type.
enum class Arithmetic { Integer, Single, Double };

Arithmetic arithmetic_of(const std::string& dtype) {
    switch (dtype_from_name(dtype)) {
        case DType::Integer: return Arithmetic::Integer;
        case DType::Float32: return Arithmetic::Single;
        case DType::Float64:
        default: return Arithmetic::Double;
    }
}

double apply(Arithmetic kind, char op, double left, double right) {
    switch (kind) {
        case Arithmetic::Integer: {
            const std::int64_t a = static_cast<std::int64_t>(left);
            const std::int64_t b = static_cast<std::int64_t>(right);
            switch (op) {
                case '+': return static_cast<double>(a + b);
                case '-': return static_cast<double>(a - b);
                default: return static_cast<double>(a * b);
            }
        }
        case Arithmetic::Single: {
            const float a = static_cast<float>(left);
            const float b = static_cast<float>(right);
            switch (op) {
                case '+': return static_cast<double>(a + b);
                case '-': return static_cast<double>(a - b);
                default: return static_cast<double>(a * b);
            }
        }
        case Arithmetic::Double:
        default:
            switch (op) {
                case '+': return left + right;
                case '-': return left - right;
                default: return left * right;
            }
    }
}

// scipy's csr_binop_csr: walk the union of the two patterns row by row, apply
// the operation with a zero for a missing operand, and store the result only
// when it is not exactly zero. Entries that are stored as an exact zero are
// skipped on the way in, because hicmatrix's fillLowerTriangle has already
// dropped them on the Python side.
//
// The merge is written once and run twice, first counting and then filling, so
// that the destination arrays are allocated at exactly the size they need. A
// reserve() guess would be wrong in both directions here: the union of the two
// patterns is larger than either operand, so the vectors would reallocate and
// briefly hold two copies, while a difference that cancels can be far smaller
// than either. Both matter at the scale these tools are budgeted for.
template <class OnEntry>
void merge_rows(const CsrMatrix& a, const CsrMatrix& b, Arithmetic kind, char op,
                OnEntry&& on_entry) {
    const std::vector<std::int64_t>& a_indptr = a.indptr();
    const std::vector<std::int32_t>& a_indices = a.indices();
    const std::vector<double>& a_data = a.data();
    const std::vector<std::int64_t>& b_indptr = b.indptr();
    const std::vector<std::int32_t>& b_indices = b.indices();
    const std::vector<double>& b_data = b.data();

    for (std::int64_t row = 0; row < a.rows(); ++row) {
        std::size_t i = static_cast<std::size_t>(a_indptr[static_cast<std::size_t>(row)]);
        const std::size_t i_end =
            static_cast<std::size_t>(a_indptr[static_cast<std::size_t>(row) + 1]);
        std::size_t j = static_cast<std::size_t>(b_indptr[static_cast<std::size_t>(row)]);
        const std::size_t j_end =
            static_cast<std::size_t>(b_indptr[static_cast<std::size_t>(row) + 1]);
        while (i < i_end && a_data[i] == 0.0) {
            ++i;
        }
        while (j < j_end && b_data[j] == 0.0) {
            ++j;
        }
        while (i < i_end || j < j_end) {
            std::int32_t column = 0;
            double left = 0.0;
            double right = 0.0;
            if (j >= j_end || (i < i_end && a_indices[i] < b_indices[j])) {
                column = a_indices[i];
                left = a_data[i];
                ++i;
            } else if (i >= i_end || b_indices[j] < a_indices[i]) {
                column = b_indices[j];
                right = b_data[j];
                ++j;
            } else {
                column = a_indices[i];
                left = a_data[i];
                right = b_data[j];
                ++i;
                ++j;
            }
            while (i < i_end && a_data[i] == 0.0) {
                ++i;
            }
            while (j < j_end && b_data[j] == 0.0) {
                ++j;
            }
            const double value = apply(kind, op, left, right);
            if (value != 0.0) {
                on_entry(row, column, value);
            }
        }
    }
}

CsrMatrix binary_op(const CsrMatrix& a, const CsrMatrix& b, char op) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        // Unreachable from hicSumMatrices and hicCompareMatrices: both compare
        // chrBinBoundaries first, and the boundaries encode the shape.
        throw std::invalid_argument("inconsistent shapes: (" + std::to_string(a.rows()) +
                                    ", " + std::to_string(a.cols()) + ") and (" +
                                    std::to_string(b.rows()) + ", " +
                                    std::to_string(b.cols()) + ")");
    }
    if (a.symmetry() != b.symmetry()) {
        throw std::invalid_argument(
            "both matrices must be stored with the same symmetry");
    }

    const std::string dtype = promote_dtype(a.dtype(), b.dtype());
    const Arithmetic kind = arithmetic_of(dtype);
    const std::int64_t rows = a.rows();

    std::vector<std::int64_t> indptr(static_cast<std::size_t>(rows) + 1, 0);
    merge_rows(a, b, kind, op, [&](std::int64_t row, std::int32_t, double) {
        ++indptr[static_cast<std::size_t>(row) + 1];
    });
    for (std::size_t row = 1; row < indptr.size(); ++row) {
        indptr[row] += indptr[row - 1];
    }

    const std::size_t nnz = static_cast<std::size_t>(indptr.back());
    std::vector<std::int32_t> indices(nnz);
    std::vector<double> data(nnz);
    std::size_t write = 0;
    merge_rows(a, b, kind, op, [&](std::int64_t, std::int32_t column, double value) {
        indices[write] = column;
        data[write] = value;
        ++write;
    });

    CsrMatrix result(rows, a.cols(), std::move(indptr), std::move(indices),
                     std::move(data), dtype);
    if (a.symmetry() == Symmetry::UpperTriangle) {
        // Nothing below the diagonal can have appeared, so this only restores
        // the flag; it never allocates.
        result.symmetrize_in_place();
    }
    return result;
}

// Visits every value of the data array of the represented matrix, in scipy's
// storage order, without materialising it.
//
// For a matrix held as an upper triangle, row i of the symmetric matrix holds
// first the mirrored entries taken from column i of the triangle, in ascending
// row order, and then the entries stored in row i itself. The mirrored entries
// are found with one cursor per row plus a per column list of the rows whose
// cursor currently sits on that column, so the traversal costs O(nbins) memory
// and no copy of the matrix.
template <class Emit>
void for_each_symmetric_value(const CsrMatrix& matrix, Emit&& emit) {
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& data = matrix.data();
    const std::int64_t rows = matrix.rows();
    const std::size_t n = static_cast<std::size_t>(rows);

    std::vector<std::int64_t> cursor(n, 0);
    std::vector<std::int64_t> column_head(n, -1);
    std::vector<std::int64_t> row_next(n, -1);

    // Moves row `r`'s cursor to its next stored entry strictly above the
    // diagonal and files the row under that column.
    const auto advance = [&](std::int64_t r) {
        const std::size_t index = static_cast<std::size_t>(r);
        std::int64_t k = cursor[index];
        const std::int64_t end = indptr[index + 1];
        while (k < end && (static_cast<std::int64_t>(indices[static_cast<std::size_t>(k)]) <= r ||
                           data[static_cast<std::size_t>(k)] == 0.0)) {
            ++k;
        }
        cursor[index] = k;
        if (k < end) {
            const std::size_t column =
                static_cast<std::size_t>(indices[static_cast<std::size_t>(k)]);
            row_next[index] = column_head[column];
            column_head[column] = r;
        }
    };

    for (std::int64_t r = 0; r < rows; ++r) {
        cursor[static_cast<std::size_t>(r)] = indptr[static_cast<std::size_t>(r)];
        advance(r);
    }

    std::vector<std::pair<std::int64_t, double>> mirrored;
    for (std::int64_t row = 0; row < rows; ++row) {
        mirrored.clear();
        for (std::int64_t r = column_head[static_cast<std::size_t>(row)]; r != -1;) {
            const std::int64_t next = row_next[static_cast<std::size_t>(r)];
            mirrored.emplace_back(r, data[static_cast<std::size_t>(cursor[static_cast<std::size_t>(r)])]);
            r = next;
        }
        column_head[static_cast<std::size_t>(row)] = -1;
        // The list is built in an arbitrary order; the data array needs the
        // mirrored entries in ascending row order, which is ascending column
        // order in the symmetric matrix.
        std::sort(mirrored.begin(), mirrored.end(),
                  [](const std::pair<std::int64_t, double>& x,
                     const std::pair<std::int64_t, double>& y) { return x.first < y.first; });
        for (const auto& entry : mirrored) {
            emit(entry.second);
        }
        for (const auto& entry : mirrored) {
            ++cursor[static_cast<std::size_t>(entry.first)];
            advance(entry.first);
        }
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            if (indices[k] < row || data[k] == 0.0) {
                continue;
            }
            emit(data[k]);
        }
    }
}

template <class Emit>
void for_each_data_value(const CsrMatrix& matrix, Emit&& emit) {
    if (matrix.symmetry() == Symmetry::UpperTriangle) {
        for_each_symmetric_value(matrix, std::forward<Emit>(emit));
        return;
    }
    for (const double value : matrix.data()) {
        emit(value);
    }
}

}  // namespace

std::string promote_dtype(const std::string& a, const std::string& b) {
    if (a == b) {
        return a;
    }
    const DType left = dtype_from_name(a);
    const DType right = dtype_from_name(b);
    if (left == DType::Integer && right == DType::Integer) {
        // The only integer widths the corpus carries are int32 and int64.
        return (a == "int64" || b == "int64") ? "int64" : "int32";
    }
    if (left == DType::Float32 && right == DType::Float32) {
        return "float32";
    }
    // numpy.promote_types(int32, float32) and promote_types(int64, float32)
    // are both float64, so every remaining mixture ends up there.
    return "float64";
}

CsrMatrix add(const CsrMatrix& a, const CsrMatrix& b) { return binary_op(a, b, '+'); }

CsrMatrix subtract(const CsrMatrix& a, const CsrMatrix& b) { return binary_op(a, b, '-'); }

CsrMatrix multiply_elementwise(const CsrMatrix& a, const CsrMatrix& b) {
    return binary_op(a, b, '*');
}

Scalar data_sum(const CsrMatrix& matrix) {
    switch (matrix.dtype_kind()) {
        case DType::Integer: {
            // np.add.reduce over an integer array accumulates in the default
            // integer type, so the result is exact and order independent.
            std::int64_t total = 0;
            for_each_data_value(matrix, [&](double value) {
                total += static_cast<std::int64_t>(value);
            });
            return Scalar::from_int(total);
        }
        case DType::Float32: {
            StreamingPairwiseSum<float> accumulator;
            for_each_data_value(matrix,
                                [&](double value) { accumulator.push(static_cast<float>(value)); });
            return Scalar::from_float(accumulator.finish());
        }
        case DType::Float64:
        default: {
            StreamingPairwiseSum<double> accumulator;
            for_each_data_value(matrix, [&](double value) { accumulator.push(value); });
            return Scalar::from_double(accumulator.finish());
        }
    }
}

void divide_data_in_place(CsrMatrix& matrix, const Scalar& divisor) {
    const double value = divisor.as_double();
    for (double& entry : matrix.mutable_data()) {
        entry /= value;
    }
    matrix.set_dtype("float64");
}

void reciprocal_data_in_place(CsrMatrix& matrix) {
    if (matrix.dtype_kind() == DType::Float32) {
        for (double& entry : matrix.mutable_data()) {
            entry = static_cast<double>(1.0F / static_cast<float>(entry));
        }
        return;
    }
    for (double& entry : matrix.mutable_data()) {
        entry = 1.0 / entry;
    }
    matrix.set_dtype("float64");
}

void log2_data_in_place(CsrMatrix& matrix) {
    if (matrix.dtype_kind() == DType::Float32) {
        for (double& entry : matrix.mutable_data()) {
            entry = static_cast<double>(std::log2(static_cast<float>(entry)));
        }
        return;
    }
    for (double& entry : matrix.mutable_data()) {
        entry = std::log2(entry);
    }
    matrix.set_dtype("float64");
}

std::vector<std::pair<std::string, BinRange>> chrom_bin_boundaries(
    const std::vector<CutInterval>& cut_intervals) {
    std::vector<std::pair<std::string, BinRange>> boundaries;
    if (cut_intervals.empty()) {
        return boundaries;
    }
    std::unordered_map<std::string, std::size_t> position;
    const auto assign = [&](const std::string& chrom, BinRange range) {
        const auto found = position.find(chrom);
        if (found != position.end()) {
            // An OrderedDict keeps the position of an existing key on
            // reassignment, which is what a chromosome appearing twice does.
            boundaries[found->second].second = range;
            return;
        }
        position.emplace(chrom, boundaries.size());
        boundaries.emplace_back(chrom, range);
    };

    std::int64_t interval_id = 0;
    std::int64_t chrom_start_id = 0;
    const std::string* previous = nullptr;
    for (const CutInterval& interval : cut_intervals) {
        if (previous == nullptr) {
            previous = &interval.chrom;
        }
        if (*previous != interval.chrom) {
            assign(*previous, BinRange{chrom_start_id, interval_id});
            chrom_start_id = interval_id;
            previous = &interval.chrom;
        }
        ++interval_id;
    }
    assign(cut_intervals.back().chrom, BinRange{chrom_start_id, interval_id});
    return boundaries;
}

void mask_and_restore_bins(MatrixData& data, const std::vector<std::int64_t>& bin_ids) {
    if (bin_ids.empty()) {
        // maskBins returns before doing anything, so restoreMaskedBins has
        // nothing to restore and the dtype is not touched either.
        return;
    }
    const std::int64_t rows = data.matrix.rows();
    std::vector<char> masked(static_cast<std::size_t>(rows), 0);
    const auto mark = [&](std::int64_t bin) {
        if (bin >= 0 && bin < rows) {
            masked[static_cast<std::size_t>(bin)] = 1;
        }
    };
    for (const std::int64_t bin : bin_ids) {
        mark(bin);
    }
    // maskBins folds the NaN bins the matrix already carries into the mask.
    for (const std::int64_t bin : data.nan_bins) {
        mark(bin);
    }

    const std::vector<std::int64_t>& indptr = data.matrix.indptr();
    const std::vector<std::int32_t>& indices = data.matrix.indices();
    const std::vector<double>& values = data.matrix.data();
    // Counted first and allocated exactly, so the surviving entries never cost
    // a second copy of the whole matrix on the way out.
    std::vector<std::int64_t> new_indptr(static_cast<std::size_t>(rows) + 1, 0);
    const auto visit_kept = [&](auto&& handle) {
        for (std::int64_t row = 0; row < rows; ++row) {
            if (masked[static_cast<std::size_t>(row)] != 0) {
                continue;
            }
            const std::size_t begin =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
            const std::size_t end =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
            for (std::size_t k = begin; k < end; ++k) {
                if (masked[static_cast<std::size_t>(indices[k])] != 0) {
                    continue;
                }
                handle(row, k);
            }
        }
    };
    visit_kept([&](std::int64_t row, std::size_t) {
        ++new_indptr[static_cast<std::size_t>(row) + 1];
    });
    for (std::size_t row = 1; row < new_indptr.size(); ++row) {
        new_indptr[row] += new_indptr[row - 1];
    }
    std::vector<std::int32_t> new_indices(static_cast<std::size_t>(new_indptr.back()));
    std::vector<double> new_data(static_cast<std::size_t>(new_indptr.back()));
    std::size_t write = 0;
    visit_kept([&](std::int64_t, std::size_t k) {
        new_indices[write] = indices[k];
        new_data[write] = values[k];
        ++write;
    });

    const Symmetry symmetry = data.matrix.symmetry();
    // restoreMaskedBins vstacks an empty csr_matrix, which numpy creates as
    // float64, so an integer matrix comes out of the round trip as float64.
    CsrMatrix result(rows, data.matrix.cols(), std::move(new_indptr),
                     std::move(new_indices), std::move(new_data), "float64");
    if (symmetry == Symmetry::UpperTriangle) {
        result.symmetrize_in_place();
    }
    data.matrix = std::move(result);

    std::vector<std::int64_t> new_nan_bins;
    for (std::int64_t bin = 0; bin < rows; ++bin) {
        if (masked[static_cast<std::size_t>(bin)] != 0) {
            new_nan_bins.push_back(bin);
        }
    }
    data.nan_bins = std::move(new_nan_bins);

    if (data.correction_factors.has_value()) {
        std::vector<double>& factors = *data.correction_factors;
        for (std::size_t bin = 0; bin < factors.size() && bin < masked.size(); ++bin) {
            if (masked[bin] != 0) {
                factors[bin] = std::nan("");
            }
        }
    }
}

}  // namespace hicx
