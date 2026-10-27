// See hicx/iterative_correction.hpp.

#include "hicx/iterative_correction.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <utility>

#include "hicx/numpy_compat.hpp"

namespace hicx::iterative_correction {

namespace {

double numpy_mean(const std::vector<double>& values) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return hicx::npy::pairwise_sum(values) / static_cast<double>(values.size());
}

// np.max with NaN propagation
double numpy_max(const std::vector<double>& values) {
    double best = -std::numeric_limits<double>::infinity();
    for (const double v : values) {
        if (std::isnan(v)) {
            return v;
        }
        best = std::max(best, v);
    }
    return best;
}

// The mean of the dense matrix whose entries are `cells`, as scipy's mean
// computes it: every entry divided by rows * cols, the row sums accumulated
// in column order and the rows reduced pairwise.
double sparse_mean(const std::map<std::pair<std::int64_t, std::int64_t>, double>& cells, std::int64_t rows,
                   std::int64_t cols) {
    const double size = static_cast<double>(rows * cols);
    std::vector<double> row_sums(static_cast<std::size_t>(rows), 0.0);
    for (const auto& [position, value] : cells) {
        row_sums[static_cast<std::size_t>(position.first)] += value / size;
    }
    return hicx::npy::pairwise_sum(row_sums);
}

}  // namespace

double asymmetry(const Coo& matrix) {
    std::map<std::pair<std::int64_t, std::int64_t>, double> cells;
    std::map<std::pair<std::int64_t, std::int64_t>, double> difference;
    for (std::size_t k = 0; k < matrix.data.size(); ++k) {
        cells[{matrix.row[k], matrix.col[k]}] += matrix.data[k];
    }
    for (const auto& [position, value] : cells) {
        difference[position] += value;
        difference[{position.second, position.first}] -= value;
    }
    for (auto& [position, value] : difference) {
        value = std::fabs(value);
    }
    const double numerator = sparse_mean(difference, matrix.cols, matrix.rows);
    const double mean = sparse_mean(cells, matrix.rows, matrix.cols);
    return numerator / (1.0 * std::fabs(mean));
}

Result correct(Coo& matrix, std::int64_t max_iterations, double tolerance, const Progress& progress) {
    const std::size_t n = static_cast<std::size_t>(matrix.rows);
    Result result;
    result.total_bias.assign(n, 1.0);
    std::vector<double>& data = matrix.data;
    std::vector<double> s(n);
    for (std::int64_t iternum = 1; iternum <= max_iterations; ++iternum) {
        std::fill(s.begin(), s.end(), 0.0);
        for (std::size_t k = 0; k < data.size(); ++k) {
            s[static_cast<std::size_t>(matrix.row[k])] += data[k] * 1.0;
        }
        std::vector<double> nonzero;
        nonzero.reserve(n);
        for (const double v : s) {
            if (!(v == 0.0)) {
                nonzero.push_back(v);
            }
        }
        const double mean = numpy_mean(nonzero);
        std::vector<double> deviations(n);
        for (std::size_t i = 0; i < n; ++i) {
            s[i] = s[i] / mean;
            result.total_bias[i] *= s[i];
            deviations[i] = std::fabs(s[i] - 1.0);
        }
        const double deviation = n == 0 ? std::numeric_limits<double>::quiet_NaN() : numpy_max(deviations);
        for (std::size_t i = 0; i < n; ++i) {
            s[i] = 1.0 / s[i];
        }
        bool large = false;
        for (std::size_t k = 0; k < data.size(); ++k) {
            data[k] *= s[static_cast<std::size_t>(matrix.row[k])];
        }
        for (std::size_t k = 0; k < data.size(); ++k) {
            data[k] *= s[static_cast<std::size_t>(matrix.col[k])];
            large = large || data[k] > 1e100;
        }
        if (large) {
            result.outcome = Outcome::LargeValuesInLoop;
            return result;
        }
        if (progress) {
            progress(iternum, deviation);
        }
        if (deviation < tolerance) {
            result.converged_at = iternum;
            break;
        }
    }
    std::vector<double> nonzero;
    for (const double v : result.total_bias) {
        if (v != 0.0) {
            nonzero.push_back(v);
        }
    }
    const double corr = numpy_mean(nonzero);
    for (double& v : result.total_bias) {
        v = v / corr;
    }
    bool large = false;
    for (double& v : data) {
        v = v * corr * corr;
        large = large || v > 1e10;
    }
    if (large) {
        result.outcome = Outcome::LargeValuesAtEnd;
    }
    return result;
}

}  // namespace hicx::iterative_correction
