#include "hicx/math/ice.hpp"

#include <cmath>
#include <vector>

#include "hicx/numpy_compat.hpp"

namespace hicx::ice {

namespace {

// np.mean over the entries of `values` selected by `keep`. numpy reduces with
// its pairwise scheme over the *compacted* array that fancy indexing produces,
// so the compaction is part of the semantics and not an implementation detail.
double mean_of_selected(const std::vector<double>& values, const std::vector<char>& keep,
                        std::vector<double>& scratch) {
    scratch.clear();
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (keep[i] != 0) {
            scratch.push_back(values[i]);
        }
    }
    if (scratch.empty()) {
        return std::nan("");  // numpy's mean of an empty slice
    }
    return npy::pairwise_sum(scratch) / static_cast<double>(scratch.size());
}

}  // namespace

Result correct(const kernels::DiagonalBlock& block, const Options& options) {
    Result result;
    const std::int64_t rows = block.size();
    const std::size_t n = static_cast<std::size_t>(rows);

    result.correction_factors.assign(n, 1.0);
    std::vector<double> s(n, 0.0);
    std::vector<double> s_next(n, 0.0);
    std::vector<double> inverse(n, 0.0);
    std::vector<char> keep(n, 0);
    std::vector<double> scratch;
    scratch.reserve(n);

    kernels::symmetric_marginals(block, s.data(), options.threads);

    for (std::int64_t iteration = 1; iteration <= options.max_iterations; ++iteration) {
        result.iterations = iteration;

        // mask = (s == 0); s = s / mean(s[~mask])
        for (std::size_t i = 0; i < n; ++i) {
            keep[i] = s[i] == 0.0 ? 0 : 1;
        }
        const double mean = mean_of_selected(s, keep, scratch);
        double deviation = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            s[i] /= mean;
            result.correction_factors[i] *= s[i];
            const double distance = std::fabs(s[i] - 1.0);
            if (distance > deviation) {
                deviation = distance;
            }
            inverse[i] = 1.0 / s[i];
        }

        const double largest =
            kernels::scale_and_marginals(block, inverse.data(), s_next.data(), options.threads);
        if (largest > 1e100) {
            result.failed = true;
            result.message =
                "matrix correction is producing extremely large values. This is often "
                "caused by bins of low counts. Use a more stringent filtering of bins.";
            return result;
        }
        if (deviation < options.tolerance) {
            result.converged = true;
            break;
        }
        s.swap(s_next);
    }

    // corr = total_bias[total_bias != 0].mean(); total_bias /= corr
    for (std::size_t i = 0; i < n; ++i) {
        keep[i] = result.correction_factors[i] == 0.0 ? 0 : 1;
    }
    const double corr = mean_of_selected(result.correction_factors, keep, scratch);
    for (std::size_t i = 0; i < n; ++i) {
        result.correction_factors[i] /= corr;
    }

    // W.data = W.data * corr * corr, two multiplications in that order.
    double largest = 0.0;
    double* data = block.data();
    for (std::int64_t local = 0; local < rows; ++local) {
        for (std::int64_t k = block.row_begin(local); k < block.row_end(local); ++k) {
            const double value = data[k] * corr * corr;
            data[k] = value;
            const double magnitude = std::fabs(value);
            if (magnitude > largest) {
                largest = magnitude;
            }
        }
    }
    if (largest > 1e10) {
        result.failed = true;
        result.message =
            "matrix correction produced extremely large values. This is often caused by "
            "bins of low counts. Use a more stringent filtering of bins.";
    }
    return result;
}

}  // namespace hicx::ice
