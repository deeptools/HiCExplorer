#include "hicx/correct_ops.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "hicx/math/sparse_kernels.hpp"

namespace hicx::correct {

namespace {

std::vector<double> block_diagonal(const kernels::DiagonalBlock& block) {
    std::vector<double> diagonal(static_cast<std::size_t>(block.size()), 0.0);
    const std::int32_t* indices = block.indices();
    const double* data = block.data();
    const std::int64_t first = block.first();
    for (std::int64_t local = 0; local < block.size(); ++local) {
        const std::int64_t begin = block.row_begin(local);
        if (begin < block.row_end(local) && indices[begin] - first == local) {
            diagonal[static_cast<std::size_t>(local)] = data[begin];
        }
    }
    return diagonal;
}

}  // namespace

double numpy_median(std::vector<double> values) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const std::size_t half = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(half),
                     values.end());
    const double upper = values[half];
    if (values.size() % 2 == 1) {
        return upper;
    }
    const double lower =
        *std::max_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(half));
    return (lower + upper) / 2.0;
}

Mad::Mad(const std::vector<double>& points) {
    std::vector<double> positive;
    positive.reserve(points.size());
    for (const double value : points) {
        if (value > 0.0) {
            positive.push_back(value);
        }
    }
    median_ = numpy_median(std::move(positive));

    std::vector<double> difference(points.size());
    std::vector<double> magnitude(points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        difference[i] = points[i] - median_;
        magnitude[i] = std::fabs(difference[i]);
    }
    mad_ = numpy_median(std::move(magnitude));

    z_.resize(points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        // np.multiply(0.6745, np.divide(diff, mad)): the division happens
        // first, so the rounding is that of a quotient scaled afterwards.
        z_[i] = 0.6745 * (difference[i] / mad_);
    }
}

std::vector<double> row_sums(CsrMatrix& matrix, int threads) {
    kernels::DiagonalBlock block(matrix);
    std::vector<double> sums(static_cast<std::size_t>(matrix.rows()), 0.0);
    kernels::symmetric_marginals(block, sums.data(), threads);
    return sums;
}

namespace {

std::vector<double> block_coverage_without_diagonal(const kernels::DiagonalBlock& block,
                                                    int threads) {
    const std::size_t size = static_cast<std::size_t>(block.size());
    std::vector<double> coverage(size, 0.0);
    kernels::symmetric_marginals(block, coverage.data(), threads);
    const std::vector<double> diagonal = block_diagonal(block);
    for (std::size_t i = 0; i < size; ++i) {
        coverage[i] -= diagonal[i];
    }
    return coverage;
}

}  // namespace

std::vector<double> coverage_without_diagonal(CsrMatrix& matrix, std::int64_t first,
                                              std::int64_t last, int threads) {
    const kernels::DiagonalBlock block(matrix, first, last);
    return block_coverage_without_diagonal(block, threads);
}

std::vector<std::int64_t> zero_coverage_bins(CsrMatrix& matrix, int threads) {
    const std::vector<double> sums = row_sums(matrix, threads);
    std::vector<std::int64_t> bins;
    for (std::size_t i = 0; i < sums.size(); ++i) {
        if (sums[i] == 0.0) {
            bins.push_back(static_cast<std::int64_t>(i));
        }
    }
    return bins;
}

std::vector<std::int64_t> filter_by_zscore(
    CsrMatrix& matrix, const std::vector<std::pair<std::string, BinRange>>& boundaries,
    double lower_threshold, double upper_threshold, bool perchr, int threads,
    std::vector<std::string>* empty_chromosome_warnings) {
    std::vector<std::int64_t> to_remove;

    const auto collect = [&](const kernels::DiagonalBlock& block, std::int64_t offset,
                             const std::string* chromosome) {
        const std::size_t size = static_cast<std::size_t>(block.size());
        const std::vector<double> coverage = block_coverage_without_diagonal(block, threads);
        const Mad mad(coverage);
        const std::vector<double>& z = mad.modified_z_scores();
        std::size_t found = 0;
        for (std::size_t i = 0; i < size; ++i) {
            if (z[i] < lower_threshold || z[i] > upper_threshold) {
                to_remove.push_back(offset + static_cast<std::int64_t>(i));
                ++found;
            }
        }
        if (found == 0 && chromosome != nullptr && empty_chromosome_warnings != nullptr) {
            empty_chromosome_warnings->push_back(*chromosome);
        }
    };

    if (perchr) {
        for (const auto& entry : boundaries) {
            kernels::DiagonalBlock block(matrix, entry.second.first, entry.second.last);
            collect(block, entry.second.first, &entry.first);
        }
        std::sort(to_remove.begin(), to_remove.end());
    } else {
        kernels::DiagonalBlock block(matrix);
        collect(block, 0, nullptr);
    }
    return to_remove;
}

void remove_diagonal(CsrMatrix& matrix) {
    CsrMatrix::Arrays arrays = matrix.release();
    std::size_t write = 0;
    std::vector<std::int64_t> indptr(arrays.indptr.size(), 0);
    for (std::int64_t row = 0; row < arrays.rows; ++row) {
        const std::size_t begin = static_cast<std::size_t>(arrays.indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(arrays.indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            if (arrays.indices[k] == static_cast<std::int32_t>(row)) {
                continue;
            }
            arrays.indices[write] = arrays.indices[k];
            arrays.data[write] = arrays.data[k];
            ++write;
        }
        indptr[static_cast<std::size_t>(row) + 1] = static_cast<std::int64_t>(write);
    }
    arrays.indices.resize(write);
    arrays.data.resize(write);
    arrays.indptr = std::move(indptr);
    matrix = CsrMatrix::adopt(std::move(arrays));
}

MaskState mask_bins_in_place(MatrixData& data, const std::vector<char>& masked) {
    MaskState state;
    state.original_intervals = data.cut_intervals;

    const std::int64_t rows = data.matrix.rows();
    std::vector<std::int32_t> renumbered(static_cast<std::size_t>(rows), -1);
    for (std::int64_t bin = 0; bin < rows; ++bin) {
        if (masked[static_cast<std::size_t>(bin)] != 0) {
            state.removed.push_back(bin);
        } else {
            renumbered[static_cast<std::size_t>(bin)] =
                static_cast<std::int32_t>(state.kept.size());
            state.kept.push_back(bin);
        }
    }

    CsrMatrix::Arrays arrays = data.matrix.release();
    const std::int64_t kept_rows = static_cast<std::int64_t>(state.kept.size());
    std::vector<std::int64_t> indptr(static_cast<std::size_t>(kept_rows) + 1, 0);
    std::size_t write = 0;
    for (std::int64_t new_row = 0; new_row < kept_rows; ++new_row) {
        const std::int64_t row = state.kept[static_cast<std::size_t>(new_row)];
        const std::size_t begin =
            static_cast<std::size_t>(arrays.indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(arrays.indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int32_t column = renumbered[static_cast<std::size_t>(arrays.indices[k])];
            if (column < 0) {
                continue;
            }
            arrays.indices[write] = column;
            arrays.data[write] = arrays.data[k];
            ++write;
        }
        indptr[static_cast<std::size_t>(new_row) + 1] = static_cast<std::int64_t>(write);
    }
    arrays.indices.resize(write);
    arrays.data.resize(write);
    arrays.indptr = std::move(indptr);
    arrays.rows = kept_rows;
    arrays.cols = kept_rows;
    data.matrix = CsrMatrix::adopt(std::move(arrays));

    std::vector<CutInterval> intervals;
    intervals.reserve(state.kept.size());
    for (const std::int64_t bin : state.kept) {
        intervals.push_back(state.original_intervals[static_cast<std::size_t>(bin)]);
    }
    data.cut_intervals = std::move(intervals);
    data.nan_bins.clear();
    if (data.correction_factors.has_value()) {
        std::vector<double> factors;
        factors.reserve(state.kept.size());
        for (const std::int64_t bin : state.kept) {
            factors.push_back((*data.correction_factors)[static_cast<std::size_t>(bin)]);
        }
        data.correction_factors = std::move(factors);
    }
    return state;
}

void restore_masked_bins(MatrixData& data, const MaskState& state) {
    if (state.removed.empty()) {
        return;
    }
    const std::int64_t rows =
        static_cast<std::int64_t>(state.kept.size() + state.removed.size());
    CsrMatrix::Arrays arrays = data.matrix.release();

    // The surviving entries are already in ascending original row order, so
    // only the column numbers and the row offsets change.
    for (std::size_t k = 0; k < arrays.indices.size(); ++k) {
        arrays.indices[k] =
            static_cast<std::int32_t>(state.kept[static_cast<std::size_t>(arrays.indices[k])]);
    }
    std::vector<std::int64_t> indptr(static_cast<std::size_t>(rows) + 1, 0);
    std::size_t kept_index = 0;
    for (std::int64_t row = 0; row < rows; ++row) {
        std::int64_t end = indptr[static_cast<std::size_t>(row)];
        if (kept_index < state.kept.size() &&
            state.kept[kept_index] == row) {
            end = arrays.indptr[kept_index + 1];
            ++kept_index;
        }
        indptr[static_cast<std::size_t>(row) + 1] = end;
    }
    arrays.indptr = std::move(indptr);
    arrays.rows = rows;
    arrays.cols = rows;
    // restoreMaskedBins stacks an empty float64 block onto the matrix, so an
    // integer matrix comes out of the round trip as float64.
    arrays.dtype = "float64";
    data.matrix = CsrMatrix::adopt(std::move(arrays));

    data.cut_intervals = state.original_intervals;
    data.nan_bins = state.removed;

    if (data.correction_factors.has_value()) {
        std::vector<double> factors(static_cast<std::size_t>(rows),
                                    std::numeric_limits<double>::quiet_NaN());
        const std::vector<double>& compact = *data.correction_factors;
        for (std::size_t i = 0; i < state.kept.size() && i < compact.size(); ++i) {
            factors[static_cast<std::size_t>(state.kept[i])] = compact[i];
        }
        data.correction_factors = std::move(factors);
    }
}

std::size_t add_missing_diagonal(CsrMatrix& matrix, double value) {
    const std::int64_t rows = matrix.rows();
    std::size_t missing = 0;
    {
        const std::vector<std::int64_t>& indptr = matrix.indptr();
        const std::vector<std::int32_t>& indices = matrix.indices();
        for (std::int64_t row = 0; row < rows; ++row) {
            const std::size_t begin =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
            const std::size_t end =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
            if (begin == end || indices[begin] != static_cast<std::int32_t>(row)) {
                ++missing;
            }
        }
    }
    if (missing == 0) {
        return 0;
    }

    CsrMatrix::Arrays arrays = matrix.release();
    const std::size_t total = arrays.data.size() + missing;
    std::vector<std::int32_t> indices(total);
    std::vector<std::int64_t> indptr(static_cast<std::size_t>(rows) + 1, 0);
    std::size_t write = 0;
    for (std::int64_t row = 0; row < rows; ++row) {
        const std::size_t begin =
            static_cast<std::size_t>(arrays.indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(arrays.indptr[static_cast<std::size_t>(row) + 1]);
        if (begin == end || arrays.indices[begin] != static_cast<std::int32_t>(row)) {
            indices[write++] = static_cast<std::int32_t>(row);
        }
        for (std::size_t k = begin; k < end; ++k) {
            indices[write++] = arrays.indices[k];
        }
        indptr[static_cast<std::size_t>(row) + 1] = static_cast<std::int64_t>(write);
    }
    // The column array is rebuilt first and the original released before the
    // value array is rebuilt, so only one of the two oversized buffers exists
    // at a time.
    arrays.indices = std::move(indices);

    // Walk the original rows again and copy the values, inserting `value`
    // wherever the row gained a column. A row gained one exactly when its new
    // extent is one longer than its old one.
    std::vector<double> data(total);
    std::size_t out = 0;
    for (std::int64_t row = 0; row < rows; ++row) {
        const std::size_t begin =
            static_cast<std::size_t>(arrays.indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(arrays.indptr[static_cast<std::size_t>(row) + 1]);
        const std::size_t out_begin =
            static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t out_end =
            static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        if (out_end - out_begin != end - begin) {
            data[out++] = value;
        }
        for (std::size_t k = begin; k < end; ++k) {
            data[out++] = arrays.data[k];
        }
    }
    arrays.data = std::move(data);
    arrays.indptr = std::move(indptr);
    matrix = CsrMatrix::adopt(std::move(arrays));
    return missing;
}

namespace {

// CPython 3.12 Objects/setobject.c, for a set of distinct non-negative Python
// ints. -1 marks an empty slot; a set built by set(iterable) never holds a
// dummy, because nothing is ever removed from it.
constexpr std::size_t kLinearProbes = 9;
constexpr std::size_t kPerturbShift = 5;
constexpr std::size_t kMinSize = 8;

struct PySet {
    std::vector<std::int64_t> table = std::vector<std::int64_t>(kMinSize, -1);
    std::size_t mask = kMinSize - 1;
    std::size_t fill = 0;
    std::size_t used = 0;
};

// set_insert_clean: no equality checks, because a resize only moves keys that
// are already known to be distinct.
void set_insert_clean(std::vector<std::int64_t>& table, std::size_t mask,
                      std::int64_t key) {
    std::size_t perturb = static_cast<std::size_t>(key);
    std::size_t i = static_cast<std::size_t>(key) & mask;
    while (true) {
        if (table[i] < 0) {
            table[i] = key;
            return;
        }
        if (i + kLinearProbes <= mask) {
            for (std::size_t j = 1; j <= kLinearProbes; ++j) {
                if (table[i + j] < 0) {
                    table[i + j] = key;
                    return;
                }
            }
        }
        perturb >>= kPerturbShift;
        i = (i * 5 + 1 + perturb) & mask;
    }
}

// set_table_resize: the smallest power of two strictly greater than minused,
// with the old table rehashed in slot order.
void set_table_resize(PySet& set, std::size_t minused) {
    std::size_t newsize = kMinSize;
    while (newsize <= minused) {
        newsize <<= 1;
    }
    std::vector<std::int64_t> table(newsize, -1);
    const std::size_t newmask = newsize - 1;
    for (const std::int64_t key : set.table) {
        if (key >= 0) {
            set_insert_clean(table, newmask, key);
        }
    }
    set.table = std::move(table);
    set.mask = newmask;
    set.fill = set.used;
}

void set_add_entry(PySet& set, std::int64_t key) {
    std::size_t perturb = static_cast<std::size_t>(key);
    std::size_t i = static_cast<std::size_t>(key) & set.mask;
    while (true) {
        // do { ... } while (probes--) runs ten times when the whole linear run
        // fits below the mask and once when it does not.
        std::size_t probes = (i + kLinearProbes <= set.mask) ? kLinearProbes : 0;
        std::size_t slot = i;
        while (true) {
            if (set.table[slot] < 0) {
                set.table[slot] = key;
                ++set.fill;
                ++set.used;
                if (set.fill * 5 >= set.mask * 3) {
                    set_table_resize(set, set.used > 50000 ? set.used * 2 : set.used * 4);
                }
                return;
            }
            if (set.table[slot] == key) {
                return;
            }
            if (probes-- == 0) {
                break;
            }
            ++slot;
        }
        perturb >>= kPerturbShift;
        i = (i * 5 + 1 + perturb) & set.mask;
    }
}

}  // namespace

std::vector<std::int64_t> cpython_set_order(const std::vector<std::int64_t>& keys) {
    PySet set;
    for (const std::int64_t key : keys) {
        set_add_entry(set, key);
    }
    std::vector<std::int64_t> order;
    order.reserve(set.used);
    for (const std::int64_t key : set.table) {
        if (key >= 0) {
            order.push_back(key);
        }
    }
    return order;
}

void keep_only_diagonal_blocks(
    CsrMatrix& matrix, const std::vector<std::pair<std::string, BinRange>>& boundaries) {
    const std::int64_t rows = matrix.rows();
    std::vector<std::int32_t> block_end(static_cast<std::size_t>(rows), 0);
    for (const auto& entry : boundaries) {
        for (std::int64_t bin = entry.second.first; bin < entry.second.last && bin < rows;
             ++bin) {
            block_end[static_cast<std::size_t>(bin)] = static_cast<std::int32_t>(entry.second.last);
        }
    }

    CsrMatrix::Arrays arrays = matrix.release();
    std::vector<std::int64_t> indptr(arrays.indptr.size(), 0);
    std::size_t write = 0;
    for (std::int64_t row = 0; row < rows; ++row) {
        const std::size_t begin =
            static_cast<std::size_t>(arrays.indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(arrays.indptr[static_cast<std::size_t>(row) + 1]);
        const std::int32_t limit = block_end[static_cast<std::size_t>(row)];
        for (std::size_t k = begin; k < end; ++k) {
            if (arrays.indices[k] >= limit) {
                continue;
            }
            arrays.indices[write] = arrays.indices[k];
            arrays.data[write] = arrays.data[k];
            ++write;
        }
        indptr[static_cast<std::size_t>(row) + 1] = static_cast<std::int64_t>(write);
    }
    arrays.indices.resize(write);
    arrays.data.resize(write);
    arrays.indptr = std::move(indptr);
    matrix = CsrMatrix::adopt(std::move(arrays));
}

}  // namespace hicx::correct
