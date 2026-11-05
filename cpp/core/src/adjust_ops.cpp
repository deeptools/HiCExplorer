#include "hicx/adjust_ops.hpp"
#include "hicx/matrix_ops.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace hicx {

namespace {

// Where each input index ends up in the output, as a CSR style pair of arrays
// rather than a vector of vectors: an index may be selected more than once,
// and a per-bin vector would cost a heap block per bin.
struct PositionIndex {
    std::vector<std::int64_t> start;    // rows + 1
    std::vector<std::int32_t> position;  // order.size()
};

PositionIndex build_positions(std::int64_t rows, const std::vector<std::int64_t>& order) {
    PositionIndex index;
    index.start.assign(static_cast<std::size_t>(rows) + 1, 0);
    for (const std::int64_t bin : order) {
        ++index.start[static_cast<std::size_t>(bin) + 1];
    }
    for (std::size_t i = 1; i < index.start.size(); ++i) {
        index.start[i] += index.start[i - 1];
    }
    index.position.resize(order.size());
    std::vector<std::int64_t> cursor(index.start.begin(), index.start.end() - 1);
    for (std::size_t i = 0; i < order.size(); ++i) {
        const std::size_t slot =
            static_cast<std::size_t>(cursor[static_cast<std::size_t>(order[i])]++);
        index.position[slot] = static_cast<std::int32_t>(i);
    }
    return index;
}

// The general gather, for a matrix whose stored entries are the whole matrix:
// one output row is one input row with its columns mapped.
CsrMatrix gather_select(const CsrMatrix& matrix, const std::vector<std::int64_t>& order,
                        const PositionIndex& index) {
    const std::size_t out_size = order.size();
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& values = matrix.data();

    std::vector<std::int64_t> new_indptr(out_size + 1, 0);
    std::vector<std::int32_t> new_indices;
    std::vector<double> new_values;
    std::vector<std::pair<std::int32_t, double>> row_buffer;
    for (std::size_t i = 0; i < out_size; ++i) {
        const std::size_t source_row = static_cast<std::size_t>(order[i]);
        const std::size_t begin = static_cast<std::size_t>(indptr[source_row]);
        const std::size_t end = static_cast<std::size_t>(indptr[source_row + 1]);
        row_buffer.clear();
        for (std::size_t k = begin; k < end; ++k) {
            const std::size_t column = static_cast<std::size_t>(indices[k]);
            const std::size_t first = static_cast<std::size_t>(index.start[column]);
            const std::size_t last = static_cast<std::size_t>(index.start[column + 1]);
            for (std::size_t p = first; p < last; ++p) {
                row_buffer.emplace_back(index.position[p], values[k]);
            }
        }
        // scipy's fancy indexing returns sorted column indices. There can be no
        // duplicate column here: an output cell (i, p) has exactly one source.
        std::sort(row_buffer.begin(), row_buffer.end(),
                  [](const std::pair<std::int32_t, double>& a,
                     const std::pair<std::int32_t, double>& b) {
                      return a.first < b.first;
                  });
        for (const auto& entry : row_buffer) {
            new_indices.push_back(entry.first);
            new_values.push_back(entry.second);
        }
        new_indptr[i + 1] = static_cast<std::int64_t>(new_indices.size());
    }
    CsrMatrix result(static_cast<std::int64_t>(out_size),
                     static_cast<std::int64_t>(out_size), std::move(new_indptr),
                     std::move(new_indices), std::move(new_values), matrix.dtype());
    result.set_symmetry(matrix.symmetry());
    return result;
}

// The selection of a matrix held as an upper triangle, producing the upper
// triangle of the result and nothing else.
//
// P A P^T is symmetric whenever A is, and the row selection and the column
// selection are the same list, so the mirror never has to exist. A stored
// entry (r, c) with r <= c stands for both A[r][c] and A[c][r], so it feeds
// every output cell whose unordered index pair maps to {r, c}: for a in the
// positions of r and b in the positions of c, the cell (min(a,b), max(a,b)).
// A cell has exactly one source, so nothing is summed and nothing is
// duplicated.
//
// This is what withdraws hicAdjustMatrix's claim on the both-triangles
// exemption of cpp/PLAN.md 4.4 rule 2: a reordering no longer needs the
// mirror, only a scatter.
CsrMatrix symmetric_select(const CsrMatrix& matrix, const std::vector<std::int64_t>& order,
                           const PositionIndex& index) {
    const std::int64_t rows = matrix.rows();
    const std::size_t out_size = order.size();
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& values = matrix.data();

    // The output row of every contribution is not produced in ascending order,
    // so the entries are counted first and scattered second. Nothing beyond
    // the result itself is allocated.
    std::vector<std::int64_t> new_indptr(out_size + 1, 0);
    const auto for_each_target = [&](auto&& emit) {
        for (std::int64_t row = 0; row < rows; ++row) {
            const std::size_t begin =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
            const std::size_t end =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
            const std::size_t row_first = static_cast<std::size_t>(index.start[static_cast<std::size_t>(row)]);
            const std::size_t row_last = static_cast<std::size_t>(index.start[static_cast<std::size_t>(row) + 1]);
            if (row_first == row_last) {
                continue;
            }
            for (std::size_t k = begin; k < end; ++k) {
                const std::size_t column = static_cast<std::size_t>(indices[k]);
                const std::size_t column_first = static_cast<std::size_t>(index.start[column]);
                const std::size_t column_last = static_cast<std::size_t>(index.start[column + 1]);
                // On the diagonal the two position lists are the same one, so
                // the unordered pair {a, b} would be produced twice. The
                // positions of a bin are stored in ascending output order, so
                // starting b at a keeps each pair exactly once.
                const bool on_diagonal = column == static_cast<std::size_t>(row);
                for (std::size_t a = row_first; a < row_last; ++a) {
                    for (std::size_t b = on_diagonal ? a : column_first;
                         b < column_last; ++b) {
                        const std::int32_t first = index.position[a];
                        const std::int32_t second = index.position[b];
                        const std::int32_t low = first < second ? first : second;
                        const std::int32_t high = first < second ? second : first;
                        emit(low, high, values[k]);
                    }
                }
            }
        }
    };

    for_each_target([&](std::int32_t low, std::int32_t, double) {
        ++new_indptr[static_cast<std::size_t>(low) + 1];
    });
    for (std::size_t i = 1; i < new_indptr.size(); ++i) {
        new_indptr[i] += new_indptr[i - 1];
    }
    const std::size_t total = static_cast<std::size_t>(new_indptr.back());
    std::vector<std::int32_t> new_indices(total);
    std::vector<double> new_values(total);
    {
        std::vector<std::int64_t> cursor(new_indptr.begin(), new_indptr.end() - 1);
        for_each_target([&](std::int32_t low, std::int32_t high, double value) {
            const std::size_t slot =
                static_cast<std::size_t>(cursor[static_cast<std::size_t>(low)]++);
            new_indices[slot] = high;
            new_values[slot] = value;
        });
    }
    // Each row still has to be put into column order, which scipy's fancy
    // indexing also guarantees. The rows are sorted one at a time, so the sort
    // buffer is the length of the longest row rather than of the matrix.
    std::vector<std::pair<std::int32_t, double>> row_buffer;
    for (std::size_t i = 0; i < out_size; ++i) {
        const std::size_t begin = static_cast<std::size_t>(new_indptr[i]);
        const std::size_t end = static_cast<std::size_t>(new_indptr[i + 1]);
        bool sorted = true;
        for (std::size_t k = begin + 1; k < end; ++k) {
            if (new_indices[k] < new_indices[k - 1]) {
                sorted = false;
                break;
            }
        }
        if (sorted) {
            continue;
        }
        row_buffer.clear();
        for (std::size_t k = begin; k < end; ++k) {
            row_buffer.emplace_back(new_indices[k], new_values[k]);
        }
        std::sort(row_buffer.begin(), row_buffer.end(),
                  [](const std::pair<std::int32_t, double>& a,
                     const std::pair<std::int32_t, double>& b) {
                      return a.first < b.first;
                  });
        for (std::size_t k = begin; k < end; ++k) {
            new_indices[k] = row_buffer[k - begin].first;
            new_values[k] = row_buffer[k - begin].second;
        }
    }

    CsrMatrix result(static_cast<std::int64_t>(out_size),
                     static_cast<std::int64_t>(out_size), std::move(new_indptr),
                     std::move(new_indices), std::move(new_values), matrix.dtype());
    result.set_symmetry(Symmetry::UpperTriangle);
    return result;
}

}  // namespace

CsrMatrix select_bins(const CsrMatrix& matrix, const std::vector<std::int64_t>& order) {
    const std::int64_t rows = matrix.rows();
    for (const std::int64_t index : order) {
        if (index < 0 || index >= rows) {
            throw std::out_of_range("bin index " + std::to_string(index) +
                                    " is outside a matrix with " +
                                    std::to_string(rows) + " rows");
        }
    }
    const PositionIndex index = build_positions(rows, order);
    if (matrix.symmetry() == Symmetry::UpperTriangle) {
        return symmetric_select(matrix, order, index);
    }
    return gather_select(matrix, order, index);
}

void reorder_bins(MatrixData& data, const std::vector<std::int64_t>& order) {
    const std::int64_t original_rows = data.matrix.rows();
    data.matrix = select_bins(data.matrix, order);

    std::vector<CutInterval> intervals;
    intervals.reserve(order.size());
    for (const std::int64_t index : order) {
        intervals.push_back(data.cut_intervals[static_cast<std::size_t>(index)]);
    }
    data.cut_intervals = std::move(intervals);

    // The NaN bins follow the permutation, and a NaN bin that is not selected
    // disappears with its row. An empty list stays empty rather than being
    // recomputed, which is what HiCMatrix.py:735-742 does.
    std::vector<std::int64_t> nan_bins;
    if (!data.nan_bins.empty()) {
        std::vector<char> was_nan(static_cast<std::size_t>(original_rows), 0);
        for (const std::int64_t bin : data.nan_bins) {
            if (bin >= 0 && bin < original_rows) {
                was_nan[static_cast<std::size_t>(bin)] = 1;
            }
        }
        for (std::size_t i = 0; i < order.size(); ++i) {
            if (was_nan[static_cast<std::size_t>(order[i])] != 0) {
                nan_bins.push_back(static_cast<std::int64_t>(i));
            }
        }
    }
    data.nan_bins = std::move(nan_bins);

    // reorderBins does not touch the correction factors. cpp/PLAN.md 2.7
    // quirk 3: reproduced, not fixed.
}

void delete_bins(MatrixData& data, const std::vector<std::int64_t>& bin_ids) {
    if (bin_ids.empty()) {
        return;  // maskBins returns before touching anything
    }
    const std::int64_t rows = data.matrix.rows();
    std::vector<char> drop(static_cast<std::size_t>(rows), 0);
    for (const std::int64_t bin : bin_ids) {
        if (bin < 0 || bin >= rows) {
            throw std::out_of_range("bin index " + std::to_string(bin) +
                                    " is outside a matrix with " +
                                    std::to_string(rows) + " rows");
        }
        drop[static_cast<std::size_t>(bin)] = 1;
    }
    // The join with the existing NaN bins at HiCMatrix.py:794-799.
    for (const std::int64_t bin : data.nan_bins) {
        if (bin >= 0 && bin < rows) {
            drop[static_cast<std::size_t>(bin)] = 1;
        }
    }

    std::vector<std::int64_t> keep;
    keep.reserve(static_cast<std::size_t>(rows));
    for (std::int64_t bin = 0; bin < rows; ++bin) {
        if (drop[static_cast<std::size_t>(bin)] == 0) {
            keep.push_back(bin);
        }
    }

    data.matrix = select_bins(data.matrix, keep);

    std::vector<CutInterval> intervals;
    intervals.reserve(keep.size());
    for (const std::int64_t index : keep) {
        intervals.push_back(data.cut_intervals[static_cast<std::size_t>(index)]);
    }
    data.cut_intervals = std::move(intervals);

    if (data.correction_factors.has_value()) {
        std::vector<double> factors;
        factors.reserve(keep.size());
        for (const std::int64_t index : keep) {
            factors.push_back((*data.correction_factors)[static_cast<std::size_t>(index)]);
        }
        data.correction_factors = std::move(factors);
    }

    // maskBins empties nan_bins when it folds them into the mask, and both
    // callers that skip restoreMaskedBins clear it explicitly anyway.
    data.nan_bins.clear();
}

std::vector<std::int64_t> empty_column_bins(const CsrMatrix& matrix) {
    const std::int64_t rows = matrix.rows();
    std::vector<double> column_sum(static_cast<std::size_t>(matrix.cols()), 0.0);
    const bool upper_only = matrix.symmetry() == Symmetry::UpperTriangle;
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& values = matrix.data();
    for (std::int64_t row = 0; row < rows; ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t column = indices[k];
            column_sum[static_cast<std::size_t>(column)] += values[k];
            if (upper_only && column != row) {
                column_sum[static_cast<std::size_t>(row)] += values[k];
            }
        }
    }
    std::vector<std::int64_t> empty;
    for (std::size_t column = 0; column < column_sum.size(); ++column) {
        if (column_sum[column] == 0.0) {
            empty.push_back(static_cast<std::int64_t>(column));
        }
    }
    return empty;
}

void zero_inter_or_intra(CsrMatrix& matrix,
                         const std::vector<std::pair<std::string, BinRange>>& boundaries,
                         InterIntra mode) {
    const std::int64_t rows = matrix.rows();
    // The block a bin belongs to. Written in boundary order, so a chromosome
    // whose range overlaps an earlier one wins, exactly as the Python's
    // successive slice assignments do.
    std::vector<std::int64_t> block_start(static_cast<std::size_t>(rows), -1);
    std::vector<std::int64_t> block_end(static_cast<std::size_t>(rows), -1);
    for (const auto& entry : boundaries) {
        const std::int64_t first = entry.second.first;
        const std::int64_t last = entry.second.last;
        for (std::int64_t bin = first; bin < last && bin < rows; ++bin) {
            if (bin < 0) {
                continue;
            }
            block_start[static_cast<std::size_t>(bin)] = first;
            block_end[static_cast<std::size_t>(bin)] = last;
        }
    }

    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    std::vector<double>& values = matrix.mutable_data();
    for (std::int64_t row = 0; row < rows; ++row) {
        const std::int64_t start = block_start[static_cast<std::size_t>(row)];
        const std::int64_t end = block_end[static_cast<std::size_t>(row)];
        if (end < 0) {
            continue;
        }
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t stop = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < stop; ++k) {
            const std::int64_t column = indices[k];
            if (column < row) {
                // Only the entries at or above the diagonal are addressed by
                // the Python's slice assignments, and the mirror entries are
                // dropped by the triu(k=0) the writers apply.
                continue;
            }
            const bool hit = (mode == InterIntra::Inter) ? (column >= end)
                                                         : (column >= start && column < end);
            if (hit) {
                values[k] = 0.0;
            }
        }
    }
    matrix.eliminate_zeros();
}

// Moved here from transform_ops.cpp, where it pulled LAPACK into every tool
// that restricts chromosomes.
void keep_only_chromosomes(MatrixData& data, const std::vector<std::string>& chromosomes) {
    const std::vector<std::pair<std::string, BinRange>> boundaries =
        chrom_bin_boundaries(data.cut_intervals);
    for (const std::string& name : chromosomes) {
        bool found = false;
        for (const std::pair<std::string, BinRange>& entry : boundaries) {
            if (entry.first == name) {
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::runtime_error("Chromosome name not in matrix. '" + name + "'");
        }
    }

    // The Python builds a boolean mask and takes np.flatnonzero of it, so the
    // selection is in ascending bin order whatever order the names were given
    // in. That is also what keeps the upper triangle representation valid.
    std::vector<bool> selected(data.cut_intervals.size(), false);
    for (const std::pair<std::string, BinRange>& entry : boundaries) {
        if (std::find(chromosomes.begin(), chromosomes.end(), entry.first) ==
            chromosomes.end()) {
            continue;
        }
        for (std::int64_t bin = entry.second.first; bin < entry.second.last; ++bin) {
            selected[static_cast<std::size_t>(bin)] = true;
        }
    }
    std::vector<std::int64_t> order;
    order.reserve(selected.size());
    for (std::size_t bin = 0; bin < selected.size(); ++bin) {
        if (selected[bin]) {
            order.push_back(static_cast<std::int64_t>(bin));
        }
    }

    data.matrix = select_bins(data.matrix, order);

    std::vector<CutInterval> kept;
    kept.reserve(order.size());
    for (std::int64_t bin : order) {
        kept.push_back(data.cut_intervals[static_cast<std::size_t>(bin)]);
    }
    data.cut_intervals = std::move(kept);

    if (data.correction_factors.has_value()) {
        std::vector<double> factors;
        factors.reserve(order.size());
        for (std::int64_t bin : order) {
            factors.push_back((*data.correction_factors)[static_cast<std::size_t>(bin)]);
        }
        data.correction_factors = std::move(factors);
    }

    if (!data.nan_bins.empty()) {
        std::vector<std::int64_t> mapped;
        std::vector<std::int64_t> position(selected.size(), -1);
        for (std::size_t k = 0; k < order.size(); ++k) {
            position[static_cast<std::size_t>(order[k])] = static_cast<std::int64_t>(k);
        }
        for (std::int64_t bin : data.nan_bins) {
            if (bin >= 0 && static_cast<std::size_t>(bin) < position.size() &&
                position[static_cast<std::size_t>(bin)] >= 0) {
                mapped.push_back(position[static_cast<std::size_t>(bin)]);
            }
        }
        std::sort(mapped.begin(), mapped.end());
        data.nan_bins = std::move(mapped);
    }

    // HiCMatrix.py:677 clears distance_counts unconditionally.
    data.distance_counts.reset();
}

}  // namespace hicx
