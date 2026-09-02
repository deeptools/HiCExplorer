#include "hicx/obsexp_ops.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "hicx/matrix_ops.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/simd_reduce.hpp"

namespace hicx {

namespace {

// np.median over a vector of integers, then int(), which truncates.
std::int64_t median_int(std::vector<std::int64_t> values) {
    if (values.empty()) {
        throw std::runtime_error(
            "cannot determine a bin size: no chromosome has more than one bin");
    }
    std::sort(values.begin(), values.end());
    const std::size_t n = values.size();
    double median = 0.0;
    if (n % 2 == 1) {
        median = static_cast<double>(values[n / 2]);
    } else {
        median = 0.5 * (static_cast<double>(values[n / 2 - 1]) +
                        static_cast<double>(values[n / 2]));
    }
    return static_cast<std::int64_t>(median);
}

// snap_nearest_multiple in HiCMatrix.fit_cut_intervals. Python's % is never
// negative for a positive modulus, and np.argmin returns the first minimum, so
// a tie snaps downwards.
std::int64_t snap_nearest_multiple(std::int64_t value, std::int64_t modulus) {
    if (modulus == 0) {
        return value;
    }
    const std::int64_t remainder = ((value % modulus) + modulus) % modulus;
    const std::int64_t down = -remainder;
    const std::int64_t up = (modulus - remainder) % modulus;
    return value + (std::llabs(down) <= std::llabs(up) ? down : up);
}

}  // namespace

std::vector<CutInterval> fit_cut_intervals(const std::vector<CutInterval>& intervals) {
    if (intervals.size() <= 1) {
        return intervals;
    }
    // Counter(chrom).items() keeps first-seen order, and only a chromosome
    // with more than one bin contributes differences.
    std::vector<std::string> order;
    std::vector<std::size_t> counts;
    for (const CutInterval& interval : intervals) {
        const auto it = std::find(order.begin(), order.end(), interval.chrom);
        if (it == order.end()) {
            order.push_back(interval.chrom);
            counts.push_back(1);
        } else {
            counts[static_cast<std::size_t>(it - order.begin())] += 1;
        }
    }
    std::vector<std::int64_t> differences;
    for (std::size_t c = 0; c < order.size(); ++c) {
        if (counts[c] <= 1) {
            continue;
        }
        // np.diff over the starts of that chromosome, in file order.
        bool have_previous = false;
        std::int64_t previous = 0;
        for (const CutInterval& interval : intervals) {
            if (interval.chrom != order[c]) {
                continue;
            }
            if (have_previous) {
                differences.push_back(interval.start - previous);
            }
            previous = interval.start;
            have_previous = true;
        }
    }
    const std::int64_t median = median_int(std::move(differences));

    std::size_t deviating = 0;
    for (const CutInterval& interval : intervals) {
        if (interval.end - interval.start != median) {
            ++deviating;
        }
    }
    if (static_cast<double>(deviating) <=
        static_cast<double>(intervals.size()) * 0.01) {
        return intervals;
    }
    std::vector<CutInterval> fitted = intervals;
    for (CutInterval& interval : fitted) {
        interval.start = snap_nearest_multiple(interval.start, median);
        interval.end = snap_nearest_multiple(interval.end, median);
    }
    return fitted;
}

void enlarge_bins(std::vector<CutInterval>& intervals) {
    if (intervals.empty()) {
        return;
    }
    bool chrom_start = true;
    for (std::size_t idx = 0; idx + 1 < intervals.size(); ++idx) {
        CutInterval& current = intervals[idx];
        const CutInterval next = intervals[idx + 1];
        if (chrom_start) {
            current.start = 0;
            chrom_start = false;
        }
        if (current.chrom == next.chrom && current.end != next.start) {
            const std::int64_t middle =
                next.start - (next.start - current.end) / 2;
            current.end = middle;
            intervals[idx + 1].start = middle;
        }
        if (current.chrom != next.chrom) {
            chrom_start = true;
        }
    }
    // The Python rewrites the last interval with its own values, which is a no
    // operation, so a single chromosome whose last bin is also its first keeps
    // its start. Reproduced by doing nothing here.
}

namespace {

// Walks the union of the stored entries of one row and the band positions
// [row, min(row + width, n)), in increasing column order, calling
// visit(column, value). This is what
//
//     (matrix + diags(ones, range(depth))).tocoo() ... data -= 1
//
// produces, without allocating the band. A position whose value plus one is
// exactly zero is skipped, because scipy's sparse addition stores only results
// that are not exactly zero.
template <class Visit>
void for_each_band_position(const CsrMatrix& matrix, bool integral,
                            std::int64_t first_bin, std::int64_t local_row,
                            std::int64_t block_size, std::int64_t width,
                            Visit&& visit) {
    // The band of ones is added to the matrix and subtracted from the COO
    // data, so a float64 value is round tripped through v + 1 - 1 and loses
    // its lowest bits. On Li_et_al_2015.h5, whose values reach 1914, that
    // shifts the per diagonal mean by one unit in the last place and the
    // z-scores by up to 5e-11 relative. An integer matrix is unaffected,
    // because the sum stays integral, so the round trip is applied only where
    // the Python applies it.
    const auto round_trip = [integral](double value) {
        return integral ? value : (value + 1.0) - 1.0;
    };
    const std::int64_t global_row = first_bin + local_row;
    const std::size_t begin =
        static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(global_row)]);
    const std::size_t end = static_cast<std::size_t>(
        matrix.indptr()[static_cast<std::size_t>(global_row) + 1]);
    const std::int64_t band_end = std::min(local_row + width, block_size);

    constexpr std::int64_t kNone = std::numeric_limits<std::int64_t>::max();
    std::size_t k = begin;
    // Skip stored entries left of the block, which the block slice removes.
    while (k < end &&
           static_cast<std::int64_t>(matrix.indices()[k]) - first_bin < local_row) {
        ++k;
    }
    std::int64_t column = local_row;
    while (true) {
        const std::int64_t band_column = column < band_end ? column : kNone;
        std::int64_t stored_column = kNone;
        if (k < end) {
            const std::int64_t candidate =
                static_cast<std::int64_t>(matrix.indices()[k]) - first_bin;
            if (candidate < block_size) {
                stored_column = candidate;
            } else {
                k = end;  // sorted, so everything after is outside the block
            }
        }
        if (band_column == kNone && stored_column == kNone) {
            return;
        }
        if (stored_column < band_column) {
            const double value = matrix.data()[k];
            if (value + 1.0 != 0.0) {
                visit(stored_column, round_trip(value));
            }
            ++k;
        } else if (stored_column == band_column) {
            const double value = matrix.data()[k];
            if (value + 1.0 != 0.0) {
                visit(band_column, round_trip(value));
            }
            ++k;
            ++column;
        } else {
            visit(band_column, 0.0);
            ++column;
        }
    }
}

}  // namespace

void convert_to_obs_exp_matrix(MatrixData& data, std::int64_t bin_size,
                               const ObsExpOptions& options) {
    if (options.max_depth_bp <= 0.0) {
        throw std::runtime_error(
            "convert_to_obs_exp_matrix requires a maxdepth; the unbounded "
            "branch of the Python densifies the whole matrix and no ported "
            "tool uses it");
    }
    if (options.max_depth_bp < static_cast<double>(bin_size)) {
        throw std::runtime_error("Please specify a maxDepth larger than bin size (" +
                                 std::to_string(bin_size) + ")");
    }

    CsrMatrix& matrix = data.matrix;
    const bool integral = matrix.integral_dtype();
    const std::int64_t n = matrix.rows();
    const std::int64_t width = static_cast<std::int64_t>(
        static_cast<double>(options.max_depth_bp) * 1.5 / static_cast<double>(bin_size));

    // Step 1: triu(m, 0) - triu(m, width). The stored entries are already the
    // upper triangle, so this keeps column - row < width, plus any NaN further
    // out, which the subtraction cannot cancel.
    {
        std::vector<std::int64_t> indptr(static_cast<std::size_t>(n) + 1, 0);
        std::vector<std::int32_t> indices;
        std::vector<double> values;
        indices.reserve(matrix.stored_nnz());
        values.reserve(matrix.stored_nnz());
        for (std::int64_t row = 0; row < n; ++row) {
            const std::size_t begin =
                static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row)]);
            const std::size_t end = static_cast<std::size_t>(
                matrix.indptr()[static_cast<std::size_t>(row) + 1]);
            for (std::size_t k = begin; k < end; ++k) {
                const std::int64_t column = static_cast<std::int64_t>(matrix.indices()[k]);
                if (column < row) {
                    continue;  // triu(k=0)
                }
                const double value = matrix.data()[k];
                if (column - row >= width && !std::isnan(value)) {
                    continue;  // cancels against triu(m, width)
                }
                if (value == 0.0) {
                    continue;  // eliminate_zeros
                }
                indices.push_back(static_cast<std::int32_t>(column));
                values.push_back(value);
            }
            indptr[static_cast<std::size_t>(row) + 1] =
                static_cast<std::int64_t>(values.size());
        }
        matrix = CsrMatrix(n, n, std::move(indptr), std::move(indices),
                           std::move(values), matrix.dtype());
        matrix.set_symmetry(Symmetry::Full);
    }

    // The chromosome blocks, and the bin counts the diagonal lengths use.
    const std::vector<std::pair<std::string, BinRange>> boundaries =
        chrom_bin_boundaries(data.cut_intervals);

    struct Block {
        std::int64_t first = 0;
        std::int64_t last = 0;
        std::vector<CutInterval> intervals;
        std::vector<std::int64_t> chrom_sizes;
    };
    std::vector<Block> blocks;
    if (options.perchr) {
        for (const auto& [chrom, range] : boundaries) {
            Block block;
            block.first = range.first;
            block.last = range.last;
            block.intervals.assign(
                data.cut_intervals.begin() + static_cast<std::ptrdiff_t>(range.first),
                data.cut_intervals.begin() + static_cast<std::ptrdiff_t>(range.last));
            block.chrom_sizes.push_back(range.last - range.first);
            blocks.push_back(std::move(block));
        }
    } else {
        Block block;
        block.first = 0;
        block.last = n;
        block.intervals = data.cut_intervals;
        for (const auto& [chrom, range] : boundaries) {
            block.chrom_sizes.push_back(range.last - range.first);
        }
        blocks.push_back(std::move(block));
    }

    // The transformed matrix is streamed into one pair of arrays in genome
    // order rather than assembled per chromosome and concatenated.
    //
    // The per chromosome blocks are independent, so processing them on several
    // threads is easy and was tried: every block writes into its own slot and
    // the slots are concatenated in block order, which keeps the result
    // identical for every thread count. It was reverted, because it was
    // measured rather than assumed. On small_test_matrix.h5 the wall clock at
    // 32 threads moved from 0.69 s to 0.67 s, three percent, while the peak RSS
    // rose from 112 MB to 187 MB, because every block's output is live at once
    // instead of one growing array. 187 MB is over this case's 184 MB budget.
    // cpp/OPTIMIZATION.md 6 says to revert a change that does not measurably
    // help, and the memory budget is a hard gate, so the sequential form
    // stands. The threading that does pay is one level up, over the bins of
    // the TAD-separation score, where it is a 3.5x wall clock gain at no
    // memory cost.
    std::vector<std::int64_t> out_indptr(static_cast<std::size_t>(n) + 1, 0);
    std::vector<std::int32_t> out_indices;
    std::vector<double> out_values;
    out_indices.reserve(matrix.stored_nnz());
    out_values.reserve(matrix.stored_nnz());
    std::int64_t written_rows = 0;

    for (const Block& block : blocks) {
        if (block.first < written_rows) {
            throw std::runtime_error(
                "chromosome blocks are not in genome order; the port cannot "
                "stream the transformed matrix");
        }
        for (std::int64_t row = written_rows; row < block.first; ++row) {
            out_indptr[static_cast<std::size_t>(row) + 1] =
                static_cast<std::int64_t>(out_values.size());
        }
        written_rows = block.first;

        const std::int64_t block_size = block.last - block.first;
        if (block_size <= 0) {
            continue;
        }
        const std::vector<CutInterval> fitted = fit_cut_intervals(block.intervals);
        std::vector<std::int64_t> starts(static_cast<std::size_t>(block_size));
        for (std::int64_t i = 0; i < block_size; ++i) {
            starts[static_cast<std::size_t>(i)] = fitted[static_cast<std::size_t>(i)].start;
        }
        // For perchr=False the chromosome id of every bin decides whether a
        // pair is inter-chromosomal, which np.unique assigns in sorted order;
        // only the equality matters, so the block-local chromosome name is
        // enough.
        const std::vector<CutInterval>& names = block.intervals;

        const auto distance_index = [&](std::int64_t i, std::int64_t j) -> std::int64_t {
            std::int64_t distance = 0;
            if (names[static_cast<std::size_t>(i)].chrom !=
                names[static_cast<std::size_t>(j)].chrom) {
                distance = -bin_size;  // dist_list[dist_list == -1] = -binsize
            } else {
                distance = starts[static_cast<std::size_t>(j)] -
                           starts[static_cast<std::size_t>(i)];
            }
            const double scaled = static_cast<double>(distance) /
                                  static_cast<double>(bin_size);
            return static_cast<std::int64_t>(scaled) + 1;  // astype(int) truncates
        };

        // Pass A: np.bincount over the band, in COO order.
        std::vector<double> sum_counts;
        std::vector<std::int64_t> distance_len;
        const auto grow = [&](std::int64_t index) {
            if (index < 0) {
                throw std::runtime_error("negative bin distance");
            }
            if (static_cast<std::size_t>(index) >= sum_counts.size()) {
                sum_counts.resize(static_cast<std::size_t>(index) + 1, 0.0);
                distance_len.resize(static_cast<std::size_t>(index) + 1, 0);
            }
        };
        for (std::int64_t i = 0; i < block_size; ++i) {
            for_each_band_position(
                matrix, integral, block.first, i, block_size, width,
                [&](std::int64_t j, double value) {
                    const std::int64_t d = distance_index(i, j);
                    grow(d);
                    sum_counts[static_cast<std::size_t>(d)] += value;
                    distance_len[static_cast<std::size_t>(d)] += 1;
                });
        }
        const std::size_t distances = sum_counts.size();

        // The diagonal lengths and the per distance means.
        std::vector<double> mu(distances, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> sigma(distances, std::numeric_limits<double>::quiet_NaN());
        std::vector<std::int64_t> lengths(distances, 0);
        for (std::size_t d = 0; d < distances; ++d) {
            if (d == 0) {
                // maxdepth is always set here, so the inter-chromosomal bucket
                // is NaN by construction (HiCMatrix.py:481-486).
                continue;
            }
            std::int64_t length = 0;
            for (const std::int64_t size : block.chrom_sizes) {
                if (size > static_cast<std::int64_t>(d) - 1) {
                    length += size - (static_cast<std::int64_t>(d) - 1);
                }
            }
            length = std::max(length, distance_len[d]);
            lengths[d] = length;
            if (length != 0) {
                mu[d] = sum_counts[d] / static_cast<double>(length);
            }
        }

        if (options.zscore) {
            // Pass B: gather (v - mu)^2 per distance, in COO order, so that
            // each group can be reduced with numpy's pairwise scheme.
            std::vector<std::size_t> offset(distances + 1, 0);
            for (std::size_t d = 0; d < distances; ++d) {
                offset[d + 1] = offset[d] + static_cast<std::size_t>(distance_len[d]);
            }
            std::vector<std::size_t> cursor(offset.begin(), offset.end() - 1);
            std::vector<double> deviations(offset[distances], 0.0);
            for (std::int64_t i = 0; i < block_size; ++i) {
                for_each_band_position(
                    matrix, integral, block.first, i, block_size, width,
                    [&](std::int64_t j, double value) {
                        const std::size_t d =
                            static_cast<std::size_t>(distance_index(i, j));
                        const double difference = value - mu[d];
                        deviations[cursor[d]++] = std::abs(difference * difference);
                    });
            }
            for (std::size_t d = 0; d < distances; ++d) {
                if (lengths[d] == 0) {
                    continue;
                }
                // Bit-identical to npy::pairwise_sum, vectorised; see
                // simd_reduce.hpp.
                const double stored =
                    simd::pairwise_sum(deviations.data() + offset[d],
                                       static_cast<std::size_t>(distance_len[d]));
                // (diagonal_length - len(values)) * mu ** 2: the square is
                // formed first, as in the Python, because ((n * mu) * mu) and
                // (n * (mu * mu)) round differently.
                const double missing =
                    static_cast<double>(lengths[d] - distance_len[d]) * (mu[d] * mu[d]);
                sigma[d] = std::sqrt((stored + missing) / static_cast<double>(lengths[d]));
            }
        }

        // Pass C: transform and emit. A LIL assignment drops exact zeros, so a
        // value that comes out zero is simply not stored; NaN is.
        const std::int64_t depth_limit = width;
        for (std::int64_t i = 0; i < block_size; ++i) {
            for_each_band_position(
                matrix, integral, block.first, i, block_size, width,
                [&](std::int64_t j, double value) {
                    const std::size_t d = static_cast<std::size_t>(distance_index(i, j));
                    double transformed = 0.0;
                    if (static_cast<std::int64_t>(d) <= depth_limit + 1) {
                        if (options.zscore) {
                            transformed = sigma[d] == 0.0
                                              ? std::numeric_limits<double>::quiet_NaN()
                                              : (value - mu[d]) / sigma[d];
                        } else {
                            transformed = value / mu[d];
                        }
                    }
                    if (transformed == 0.0) {
                        return;
                    }
                    out_indices.push_back(
                        static_cast<std::int32_t>(block.first + j));
                    out_values.push_back(transformed);
                });
            out_indptr[static_cast<std::size_t>(block.first + i) + 1] =
                static_cast<std::int64_t>(out_values.size());
        }
        written_rows = block.last;
    }
    for (std::int64_t row = written_rows; row < n; ++row) {
        out_indptr[static_cast<std::size_t>(row) + 1] =
            static_cast<std::int64_t>(out_values.size());
    }

    matrix = CsrMatrix(n, n, std::move(out_indptr), std::move(out_indices),
                       std::move(out_values), "float64");
    matrix.set_symmetry(Symmetry::Full);
}

}  // namespace hicx
