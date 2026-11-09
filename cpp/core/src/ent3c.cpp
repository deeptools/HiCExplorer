#include "hicx/ent3c.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "hicx/lapack_shim.hpp"
#include "hicx/transform_ops.hpp"

namespace hicx::ent3c {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Python's floor division for a positive divisor.
std::int64_t floor_div(std::int64_t a, std::int64_t b) {
    std::int64_t q = a / b;
    if ((a % b != 0) && (a < 0)) {
        --q;
    }
    return q;
}

// The eigenvalues of a symmetric k x k matrix whose upper triangle is stored
// column major in `matrix` (overwritten).
std::vector<double> symmetric_eigenvalues(std::vector<double>& matrix, int k) {
    const char jobz = 'N';
    const char range = 'A';
    const char uplo = 'U';
    const double vl_unused = 0.0;
    const double vu_unused = 0.0;
    const int il_unused = 0;
    const int iu_unused = 0;
    const char safe_minimum = 'S';
    const double abstol = 2.0 * dlamch_(&safe_minimum);
    int found = 0;
    std::vector<double> values(static_cast<std::size_t>(k));
    double z_unused = 0.0;
    const int ldz = 1;
    std::vector<int> isuppz(static_cast<std::size_t>(2 * k));
    int info = 0;
    int lwork = -1;
    int liwork = -1;
    double work_query = 0.0;
    int iwork_query = 0;
    dsyevr_(&jobz, &range, &uplo, &k, matrix.data(), &k, &vl_unused, &vu_unused, &il_unused,
            &iu_unused, &abstol, &found, values.data(), &z_unused, &ldz, isuppz.data(), &work_query,
            &lwork, &iwork_query, &liwork, &info);
    if (info != 0) {
        throw std::runtime_error("ENT3C: the eigenvalue workspace query failed (LAPACK info " +
                                 std::to_string(info) + ")");
    }
    lwork = static_cast<int>(work_query);
    liwork = iwork_query;
    std::vector<double> work(static_cast<std::size_t>(std::max(lwork, 1)));
    std::vector<int> iwork(static_cast<std::size_t>(std::max(liwork, 1)));
    dsyevr_(&jobz, &range, &uplo, &k, matrix.data(), &k, &vl_unused, &vu_unused, &il_unused,
            &iu_unused, &abstol, &found, values.data(), &z_unused, &ldz, isuppz.data(), work.data(),
            &lwork, iwork.data(), &liwork, &info);
    if (info != 0) {
        throw std::runtime_error("ENT3C: the eigenvalue computation did not converge (LAPACK info " +
                                 std::to_string(info) + ")");
    }
    values.resize(static_cast<std::size_t>(found));
    return values;
}

}  // namespace

Windows windows(std::int64_t bins, const Options& options) {
    if (options.submatrix_size < 0) {
        throw std::invalid_argument("the ENT3C submatrix size must be 0 (derived) or positive");
    }
    if (options.submatrix_size == 0 && options.chromosome_split < 1) {
        throw std::invalid_argument("the ENT3C chromosome split must be at least 1");
    }
    if (options.window_shift < 1) {
        throw std::invalid_argument("the ENT3C window shift must be at least 1");
    }
    if (options.max_windows < 0) {
        throw std::invalid_argument("the ENT3C maximum number of windows must be 0 (no limit) or positive");
    }
    Windows w;
    w.size = options.submatrix_size > 0
                 ? options.submatrix_size
                 : static_cast<std::int64_t>(std::floor(static_cast<double>(bins) /
                                                            static_cast<double>(options.chromosome_split) +
                                                        0.5));
    w.shift = options.window_shift;
    if (w.size < 1 || w.size > bins) {
        w.count = 0;
        return w;
    }
    w.count = 1 + floor_div(bins - w.size, w.shift);
    while (options.max_windows > 0 && w.count > options.max_windows) {
        ++w.shift;
        w.count = 1 + floor_div(bins - w.size, w.shift);
    }
    return w;
}

std::vector<bool> empty_bins(const ChromosomeEntries& entries) {
    std::vector<double> sums(static_cast<std::size_t>(entries.nbins), 0.0);
    for (std::size_t i = 0; i < entries.values.size(); ++i) {
        const double v = entries.values[i];
        if (std::isnan(v)) {
            continue;  // numpy.nansum
        }
        sums[static_cast<std::size_t>(entries.rows[i])] += v;
        if (entries.rows[i] != entries.cols[i]) {
            sums[static_cast<std::size_t>(entries.cols[i])] += v;
        }
    }
    std::vector<bool> empty(sums.size());
    for (std::size_t i = 0; i < sums.size(); ++i) {
        empty[i] = sums[i] == 0.0;
    }
    return empty;
}

double window_entropy(const std::vector<double>& counts, std::int64_t n) {
    const auto size = static_cast<std::size_t>(n);
    // Logarithms; a zero (and anything else without a finite logarithm) is
    // missing.
    std::vector<double> logs(size * size);
    for (std::size_t i = 0; i < logs.size(); ++i) {
        const double v = counts[i] > 0.0 ? std::log(counts[i]) : kNaN;
        logs[i] = std::isfinite(v) ? v : kNaN;
    }
    // A bin whose column of logarithms sums to 0, missing values left out,
    // leaves the window as row and column: a column of zeros and ones.
    std::vector<std::size_t> kept;
    for (std::size_t col = 0; col < size; ++col) {
        double sum = 0.0;
        for (std::size_t row = 0; row < size; ++row) {
            const double v = logs[row * size + col];
            if (!std::isnan(v)) {
                sum += v;
            }
        }
        if (sum != 0.0) {
            kept.push_back(col);
        }
    }
    const int k = static_cast<int>(kept.size());
    const auto kk = static_cast<std::size_t>(k);
    // The remaining missing values become the smallest remaining value.
    double smallest = std::numeric_limits<double>::infinity();
    for (const std::size_t r : kept) {
        for (const std::size_t c : kept) {
            const double v = logs[r * size + c];
            if (!std::isnan(v)) {
                smallest = std::min(smallest, v);
            }
        }
    }
    if (k < 2 || !std::isfinite(smallest)) {
        return kNaN;
    }
    std::vector<double> window(kk * kk);
    for (std::size_t r = 0; r < kk; ++r) {
        for (std::size_t c = 0; c < kk; ++c) {
            const double v = logs[kept[r] * size + kept[c]];
            window[r * kk + c] = std::isnan(v) ? smallest : v;
        }
    }

    // numpy.corrcoef(x) with the k remaining rows as variables: centre each
    // row, form the Gram matrix, scale by the standard deviations. The window
    // is symmetric, so x's rows are its columns. `centred` is x row major,
    // which is x^T column major, so dsyrk with trans = 'T' forms x x^T.
    std::vector<double> centred(kk * kk);
    for (std::size_t r = 0; r < kk; ++r) {
        double mean = 0.0;
        for (std::size_t c = 0; c < kk; ++c) {
            mean += window[r * kk + c];
        }
        mean /= static_cast<double>(k);
        for (std::size_t c = 0; c < kk; ++c) {
            centred[r * kk + c] = window[r * kk + c] - mean;
        }
    }
    std::vector<double> gram(kk * kk, 0.0);
    const char uplo = 'U';
    const char trans = 'T';
    const double one = 1.0;
    const double zero = 0.0;
    dsyrk_(&uplo, &trans, &k, &k, &one, centred.data(), &k, &zero, gram.data(), &k);
    std::vector<double> deviation(kk);
    for (std::size_t i = 0; i < kk; ++i) {
        deviation[i] = std::sqrt(gram[i + i * kk]);
    }
    // rho = P / k, upper triangle, column major: element (i, j), i <= j, at
    // i + j * k. numpy clips the correlations to [-1, 1]. A column that is
    // still constant has no correlation (0/0); it correlates 0 with the
    // others and 1 with itself, as every diagonal entry is set to 1.
    const double scale = 1.0 / static_cast<double>(k);
    for (std::size_t j = 0; j < kk; ++j) {
        for (std::size_t i = 0; i < j; ++i) {
            double p = gram[i + j * kk] / deviation[i] / deviation[j];
            p = std::isnan(p) ? 0.0 : std::clamp(p, -1.0, 1.0);
            gram[i + j * kk] = p * scale;
        }
        gram[j + j * kk] = scale;
    }
    const std::vector<double> lambdas = symmetric_eigenvalues(gram, k);
    double entropy = 0.0;
    for (const double lambda : lambdas) {
        if (lambda > 0.0) {
            entropy -= lambda * std::log(lambda);
        }
    }
    return entropy;
}

std::vector<double> entropy_signal(const ChromosomeEntries& entries,
                                   const std::vector<std::int64_t>& keep, const Windows& windows,
                                   int threads) {
    std::vector<double> signal(static_cast<std::size_t>(std::max<std::int64_t>(windows.count, 0)), kNaN);
    if (signal.empty()) {
        return signal;
    }
    // Both triangles of the matrix over the kept bins, by row.
    std::vector<std::int64_t> local(static_cast<std::size_t>(entries.nbins), -1);
    for (std::size_t i = 0; i < keep.size(); ++i) {
        local[static_cast<std::size_t>(keep[i])] = static_cast<std::int64_t>(i);
    }
    struct Contact {
        std::int64_t col;
        double value;
    };
    std::vector<std::vector<Contact>> rows(keep.size());
    for (std::size_t i = 0; i < entries.values.size(); ++i) {
        const std::int64_t r = local[static_cast<std::size_t>(entries.rows[i])];
        const std::int64_t c = local[static_cast<std::size_t>(entries.cols[i])];
        if (r < 0 || c < 0) {
            continue;
        }
        rows[static_cast<std::size_t>(r)].push_back({c, entries.values[i]});
        if (r != c) {
            rows[static_cast<std::size_t>(c)].push_back({r, entries.values[i]});
        }
    }

    const std::int64_t n = windows.size;
    const auto size = static_cast<std::size_t>(n);
    const auto compute = [&](std::int64_t w) {
        const std::int64_t first = w * windows.shift;
        std::vector<double> counts(size * size, 0.0);
        for (std::int64_t r = first; r < first + n; ++r) {
            for (const Contact& contact : rows[static_cast<std::size_t>(r)]) {
                if (contact.col >= first && contact.col < first + n) {
                    counts[static_cast<std::size_t>(r - first) * size +
                           static_cast<std::size_t>(contact.col - first)] = contact.value;
                }
            }
        }
        signal[static_cast<std::size_t>(w)] = window_entropy(counts, n);
    };

    hicx::pin_blas_to_one_thread();
    const int workers =
        static_cast<int>(std::min<std::int64_t>(std::max(threads, 1), windows.count));
    if (workers <= 1) {
        for (std::int64_t w = 0; w < windows.count; ++w) {
            compute(w);
        }
        return signal;
    }
    std::atomic<std::int64_t> next{0};
    std::vector<std::thread> pool;
    std::exception_ptr failure;
    std::atomic<bool> failed{false};
    for (int t = 0; t < workers; ++t) {
        pool.emplace_back([&] {
            try {
                for (std::int64_t w = next++; w < windows.count && !failed; w = next++) {
                    compute(w);
                }
            } catch (...) {
                if (!failed.exchange(true)) {
                    failure = std::current_exception();
                }
            }
        });
    }
    for (std::thread& thread : pool) {
        thread.join();
    }
    if (failure) {
        std::rethrow_exception(failure);
    }
    return signal;
}

double similarity(const std::vector<double>& s1, const std::vector<double>& s2) {
    if (s1.size() != s2.size()) {
        return kNaN;
    }
    // Only the windows with an entropy in both signals.
    std::vector<double> x;
    std::vector<double> y;
    for (std::size_t i = 0; i < s1.size(); ++i) {
        if (!std::isnan(s1[i]) && !std::isnan(s2[i])) {
            x.push_back(s1[i]);
            y.push_back(s2[i]);
        }
    }
    const std::size_t n = x.size();
    if (n < 2) {
        return kNaN;
    }
    double mean1 = 0.0;
    double mean2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        mean1 += x[i];
        mean2 += y[i];
    }
    mean1 /= static_cast<double>(n);
    mean2 /= static_cast<double>(n);
    double sxy = 0.0;
    double sxx = 0.0;
    double syy = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double dx = x[i] - mean1;
        const double dy = y[i] - mean2;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    const double r = sxy / std::sqrt(sxx) / std::sqrt(syy);
    return std::isnan(r) ? kNaN : std::clamp(r, -1.0, 1.0);
}

}  // namespace hicx::ent3c
