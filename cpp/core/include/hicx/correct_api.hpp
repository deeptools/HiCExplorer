// The computations of hicexplorer.hicCorrectMatrix's public functions (MAD,
// filter_by_zscore, fill_gaps, the data of plot_total_contact_dist) in the
// dtype rules of numpy and scipy, for the public Python API
// (cpp/python/src/api/hicCorrectMatrix.cpp). Header only, templated on the
// element type of the matrix.
//
//   * row sums: scipy sums the rows of a CSR with np.add.reduceat over each
//     row's stored values, the first value plus numpy's pairwise sum of the
//     others, in the data dtype; the diagonal is csr.diagonal()
//   * np.median: the middle element, or the mean of the two middle ones, and
//     NaN as soon as one value is NaN (hicx::correct::numpy_median)
//   * MAD: median of the positive points, diff the row sums of points -
//     median, med_abs_deviation the median of |diff|, z = 0.6745 * (diff / mad)
//   * fill_gaps: the lil assignments of the Python on a copy of the matrix,
//     every value cast to the dtype, an assigned zero deleting the entry, the
//     reads with scipy's index checks

#ifndef HICX_CORRECT_API_HPP
#define HICX_CORRECT_API_HPP

#include <cmath>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <map>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

#include "hicx/correct_ops.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/python_error.hpp"
#include "hicx/scipy_coo.hpp"

namespace hicx::correct_api {

inline constexpr double kMadB = 0.6745;

// np.median with numpy's NaN rule.
inline double median(std::vector<double> values) {
    for (const double v : values) {
        if (std::isnan(v)) {
            return v;
        }
    }
    return hicx::correct::numpy_median(std::move(values));
}

struct MadStats {
    double median = 0.0;
    double med_abs_deviation = 0.0;
    std::vector<double> modified_z_score;
};

// MAD.__init__ of a (rows, cols) array of points in C order.
inline MadStats mad(const std::vector<double>& points, std::int64_t rows, std::int64_t cols) {
    MadStats out;
    std::vector<double> positive;
    for (const double v : points) {
        if (v > 0) {
            positive.push_back(v);
        }
    }
    out.median = median(positive);
    std::vector<double> diff(static_cast<std::size_t>(rows));
    std::vector<double> row(static_cast<std::size_t>(cols));
    for (std::int64_t i = 0; i < rows; ++i) {
        for (std::int64_t j = 0; j < cols; ++j) {
            row[static_cast<std::size_t>(j)] = points[static_cast<std::size_t>(i * cols + j)] - out.median;
        }
        diff[static_cast<std::size_t>(i)] = cols == 1 ? row[0] : hicx::npy::pairwise_sum(row);
    }
    std::vector<double> magnitude(diff.size());
    for (std::size_t i = 0; i < diff.size(); ++i) {
        magnitude[i] = std::fabs(diff[i]);
    }
    out.med_abs_deviation = median(magnitude);
    out.modified_z_score.resize(diff.size());
    for (std::size_t i = 0; i < diff.size(); ++i) {
        out.modified_z_score[i] = kMadB * (diff[i] / out.med_abs_deviation);
    }
    return out;
}

// matrix[first:last, first:last].sum(axis=1) - .diagonal() of a CSR given as
// indptr, column indices and data; NaN values are set to zero first with
// `zero_nan`. Integer matrices give int64, float64 ones double.
template <class T>
std::vector<std::conditional_t<std::is_integral_v<T>, std::int64_t, double>> coverage(
    const std::vector<std::int64_t>& indptr, const std::vector<std::int64_t>& columns, const std::vector<T>& data,
    std::int64_t rows, std::int64_t first, std::int64_t last, bool zero_nan) {
    using Acc = std::conditional_t<std::is_integral_v<T>, std::int64_t, double>;
    const std::int64_t n = std::max<std::int64_t>(0, std::min(last, rows) - first);
    std::vector<Acc> sums(static_cast<std::size_t>(n), Acc(0));
    std::vector<T> segment;
    for (std::int64_t i = 0; i < n; ++i) {
        const std::int64_t r = first + i;
        segment.clear();
        T diagonal = T(0);
        for (std::int64_t k = indptr[static_cast<std::size_t>(r)]; k < indptr[static_cast<std::size_t>(r) + 1]; ++k) {
            const std::int64_t c = columns[static_cast<std::size_t>(k)];
            if (c < first || c >= first + n) {
                continue;
            }
            T value = data[static_cast<std::size_t>(k)];
            if constexpr (std::is_floating_point_v<T>) {
                if (zero_nan && std::isnan(value)) {
                    value = T(0);
                }
            }
            segment.push_back(value);
            if (c == r) {
                diagonal = hicx::scipy_coo::add(diagonal, value);
            }
        }
        const T sum = segment.empty() ? T(0) : hicx::scipy_coo::reduceat_run(segment);
        sums[static_cast<std::size_t>(i)] = static_cast<Acc>(sum) - static_cast<Acc>(diagonal);
    }
    return sums;
}

// numpy's float64 -> T cast.
template <class T>
T cast_from_double(double value) {
    if constexpr (std::is_floating_point_v<T>) {
        return static_cast<T>(value);
    } else {
        if (std::isnan(value) || value < -9.2233720368547758e18 || value >= 9.2233720368547758e18) {
            return std::numeric_limits<T>::min();
        }
        return static_cast<T>(static_cast<std::int64_t>(value));
    }
}

// A matrix as rows of sorted cells: what a lil_matrix holds.
template <class T>
class Grid {
  public:
    Grid(std::int64_t rows, std::int64_t cols, const std::vector<std::int64_t>& row,
         const std::vector<std::int64_t>& col, const std::vector<T>& data)
        : rows_(rows), cols_(cols), cells_(static_cast<std::size_t>(rows)) {
        for (std::size_t k = 0; k < data.size(); ++k) {
            auto& cells = cells_[static_cast<std::size_t>(row[k])];
            auto [it, inserted] = cells.emplace(col[k], data[k]);
            if (!inserted) {
                it->second = hicx::scipy_coo::add(it->second, data[k]);
            }
        }
    }

    // matrix[i, j] with scipy's bounds check
    [[nodiscard]] T at(std::int64_t i, std::int64_t j) const {
        if (i < -rows_ || i >= rows_) {
            throw hicx::PythonError("IndexError", "row index (" + std::to_string(i) + ") out of range");
        }
        if (j < -cols_ || j >= cols_) {
            throw hicx::PythonError("IndexError", "column index (" + std::to_string(j) + ") out of range");
        }
        if (i < 0) {
            i += rows_;
        }
        if (j < 0) {
            j += cols_;
        }
        const auto& cells = cells_[static_cast<std::size_t>(i)];
        const auto it = cells.find(j);
        return it == cells.end() ? T(0) : it->second;
    }

    // lil assignment of one value: a zero deletes the entry
    void set(std::int64_t i, std::int64_t j, T value) {
        auto& cells = cells_[static_cast<std::size_t>(i)];
        if (hicx::scipy_coo::nonzero(value)) {
            cells[j] = value;
        } else {
            cells.erase(j);
        }
    }

    [[nodiscard]] hicx::scipy_coo::Csr<T> csr() const {
        hicx::scipy_coo::Csr<T> out;
        out.rows = rows_;
        out.cols = cols_;
        out.indptr.push_back(0);
        for (const auto& cells : cells_) {
            for (const auto& [j, v] : cells) {
                out.indices.push_back(j);
                out.data.push_back(v);
            }
            out.indptr.push_back(static_cast<std::int64_t>(out.data.size()));
        }
        return out;
    }

  private:
    std::int64_t rows_;
    std::int64_t cols_;
    std::vector<std::map<std::int64_t, T>> cells_;
};

// np.mean of the four neighbours, as the dtype the lil assignment casts it
// to: float32 values are averaged in float32, others in float64.
template <class T>
T fill_value(T a, T b, T c, T d) {
    if constexpr (std::is_same_v<T, float>) {
        const std::vector<float> v{a, b, c, d};
        return hicx::npy::pairwise_sum(v) / 4.0F;
    } else {
        const std::vector<double> v{static_cast<double>(a), static_cast<double>(b), static_cast<double>(c),
                                    static_cast<double>(d)};
        const double mean = hicx::npy::pairwise_sum(v) / 4.0;
        if constexpr (std::is_floating_point_v<T>) {
            return mean;
        } else {
            if (std::isnan(mean)) {
                throw hicx::PythonError("ValueError", "cannot convert float NaN to integer");
            }
            return cast_from_double<T>(mean);
        }
    }
}

// fill_gaps on the matrix given as a COO in storage order (duplicates
// summed), with the failed bins that are not part of a run.
template <class T>
hicx::scipy_coo::Csr<T> fill_gaps(std::int64_t rows, std::int64_t cols, const std::vector<std::int64_t>& row,
                                  const std::vector<std::int64_t>& col, const std::vector<T>& data,
                                  const std::vector<std::int64_t>& discontinuous) {
    const Grid<T> original(rows, cols, row, col, data);
    Grid<T> fill(rows, cols, row, col, data);
    const std::int64_t mat_size = rows;
    for (const std::int64_t missing_bin : discontinuous) {
        if (0 < missing_bin && missing_bin < mat_size - 1) {
            // fill_ma[missing_bin, 1:n - 1] = (M[missing_bin - 1, :n - 2] + M[missing_bin + 1, 2:]) / 2
            std::vector<T> values;
            for (std::int64_t j = 1; j < mat_size - 1; ++j) {
                const T sum = hicx::scipy_coo::add(original.at(missing_bin - 1, j - 1), original.at(missing_bin + 1, j + 1));
                values.push_back(cast_from_double<T>(static_cast<double>(sum) * 0.5));
            }
            for (std::int64_t j = 1; j < mat_size - 1; ++j) {
                fill.set(missing_bin, j, values[static_cast<std::size_t>(j - 1)]);
            }
            // the same for the column
            values.clear();
            for (std::int64_t i = 1; i < mat_size - 1; ++i) {
                const T sum = hicx::scipy_coo::add(original.at(i - 1, missing_bin - 1), original.at(i + 1, missing_bin + 1));
                values.push_back(cast_from_double<T>(static_cast<double>(sum) * 0.5));
            }
            for (std::int64_t i = 1; i < mat_size - 1; ++i) {
                fill.set(i, missing_bin, values[static_cast<std::size_t>(i - 1)]);
            }
        }
    }
    for (const std::int64_t a : discontinuous) {
        for (const std::int64_t b : discontinuous) {
            if (0 < a && a < mat_size && 0 < b && b < mat_size) {
                const T v1 = original.at(a - 1, b - 1);
                const T v2 = original.at(a - 1, b + 1);
                const T v3 = original.at(a + 1, b - 1);
                const T v4 = original.at(a + 1, b + 1);
                const T value = fill_value(v1, v2, v3, v4);
                fill.set(a - 1, b, value);
                fill.set(a + 1, b, value);
                fill.set(a, b - 1, value);
                fill.set(a, b + 1, value);
            }
        }
    }
    return fill.csr();
}

// The left edge of the first local minimum of the 100 bin histogram
// plot_total_contact_dist draws, ax1.hist(values, 100): x with
// 1 <= x < 99 and dist[x - 1] > dist[x] < dist[x + 1], or nothing. The
// histogram is np.histogram(values, 100, range=(min, max)) as matplotlib
// calls it: np.linspace edges and the uniform-bin index computation with
// numpy's one-ulp corrections. Shared with the hicCorrectMatrix tool's
// diagnostic_plot.
inline std::optional<double> histogram_threshold(const std::vector<double>& values) {
    constexpr int bins = 100;
    double first = 0.0;
    double last = 1.0;
    if (!values.empty()) {
        first = *std::min_element(values.begin(), values.end());
        last = *std::max_element(values.begin(), values.end());
    }
    if (first == last) {
        first = first - 0.5;
        last = last + 0.5;
    }
    // np.linspace(first, last, bins + 1)
    std::vector<double> edges(bins + 1);
    const double delta = last - first;
    const double step = delta / bins;
    for (int i = 0; i <= bins; ++i) {
        double y = static_cast<double>(i);
        if (step == 0.0) {
            y = y / bins;
            y = y * delta;
        } else {
            y = y * step;
        }
        edges[static_cast<std::size_t>(i)] = y + first;
    }
    edges[bins] = last;
    std::vector<std::int64_t> counts(bins, 0);
    const double norm_denom = last - first;
    for (const double a : values) {
        if (!(a >= first && a <= last)) {
            continue;
        }
        auto index = static_cast<std::int64_t>(((a - first) / norm_denom) * bins);
        if (index == bins) {
            index -= 1;
        }
        if (a < edges[static_cast<std::size_t>(index)]) {
            index -= 1;
        }
        if (a >= edges[static_cast<std::size_t>(index + 1)] && index != bins - 1) {
            index += 1;
        }
        counts[static_cast<std::size_t>(index)] += 1;
    }
    for (int x = 1; x < bins - 1; ++x) {
        const std::size_t i = static_cast<std::size_t>(x);
        if (counts[i - 1] > counts[i] && counts[i] < counts[i + 1]) {
            return edges[i];
        }
    }
    return std::nullopt;
}

}  // namespace hicx::correct_api

#endif  // HICX_CORRECT_API_HPP
