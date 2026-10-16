#include "hicx/genomedisco.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

namespace hicx::genomedisco {

namespace {

// A plain, general (not upper-triangle-only) CSR matrix for the random-walk
// arithmetic: hicx::CsrMatrix is built around the upper-triangle Hi-C
// convention, and GenomeDISCO's transition matrix is not symmetric once
// to_transition divides unevenly sized rows, so this stays separate rather
// than bending that type's invariants.
struct FullCsr {
    std::int64_t n = 0;
    std::vector<std::int64_t> indptr;
    std::vector<std::int32_t> indices;
    std::vector<double> data;
};

FullCsr from_triplets(std::int64_t n, std::vector<std::int32_t> rows,
                      std::vector<std::int32_t> cols, std::vector<double> vals) {
    std::vector<std::size_t> order(rows.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return rows[a] != rows[b] ? rows[a] < rows[b] : cols[a] < cols[b];
    });
    FullCsr out;
    out.n = n;
    out.indptr.assign(static_cast<std::size_t>(n) + 1, 0);
    out.indices.reserve(order.size());
    out.data.reserve(order.size());
    std::int32_t prev_row = -1;
    std::int32_t prev_col = -1;
    for (const std::size_t idx : order) {
        const std::int32_t r = rows[idx];
        const std::int32_t c = cols[idx];
        if (r == prev_row && c == prev_col) {
            out.data.back() += vals[idx];
            continue;
        }
        while (prev_row < r) {
            ++prev_row;
            out.indptr[static_cast<std::size_t>(prev_row)] = static_cast<std::int64_t>(out.indices.size());
        }
        out.indices.push_back(c);
        out.data.push_back(vals[idx]);
        prev_col = c;
    }
    while (prev_row < n) {
        ++prev_row;
        out.indptr[static_cast<std::size_t>(prev_row)] = static_cast<std::int64_t>(out.indices.size());
    }
    return out;
}

// m1up + m1down with m1down.setdiag(0): mirrors the (diagonal-free) upper
// triangle entries into a full symmetric matrix.
FullCsr symmetrize(const ChromosomeEntries& e, const std::vector<double>& values) {
    std::vector<std::int32_t> rows;
    std::vector<std::int32_t> cols;
    std::vector<double> vals;
    rows.reserve(e.rows.size() * 2);
    cols.reserve(e.rows.size() * 2);
    vals.reserve(e.rows.size() * 2);
    for (std::size_t i = 0; i < e.rows.size(); ++i) {
        rows.push_back(e.rows[i]);
        cols.push_back(e.cols[i]);
        vals.push_back(values[i]);
        rows.push_back(e.cols[i]);
        cols.push_back(e.rows[i]);
        vals.push_back(values[i]);
    }
    return from_triplets(e.nbins, std::move(rows), std::move(cols), std::move(vals));
}

// data_operations.sqrtvc: each stored (upper-triangle, diagonal-free) value
// divided by sqrt(row coverage) * sqrt(column coverage), where coverage is
// the row sum of the symmetrized matrix (every entry touching that row,
// counted once per entry it appears in).
std::vector<double> sqrtvc_values(const ChromosomeEntries& e) {
    std::vector<double> coverage(static_cast<std::size_t>(e.nbins), 0.0);
    for (std::size_t i = 0; i < e.rows.size(); ++i) {
        coverage[static_cast<std::size_t>(e.rows[i])] += e.values[i];
        coverage[static_cast<std::size_t>(e.cols[i])] += e.values[i];
    }
    std::vector<double> sqrt_coverage(coverage.size());
    for (std::size_t i = 0; i < coverage.size(); ++i) {
        sqrt_coverage[i] = coverage[i] == 0.0 ? 1.0 : std::sqrt(coverage[i]);
    }
    std::vector<double> out(e.values.size());
    for (std::size_t i = 0; i < e.values.size(); ++i) {
        out[i] = e.values[i] / (sqrt_coverage[static_cast<std::size_t>(e.rows[i])] *
                                sqrt_coverage[static_cast<std::size_t>(e.cols[i])]);
    }
    return out;
}

// to_transition: row-normalize so each row sums to 1; a zero row is divided
// by 1 instead (left all-zero), matching data_operations.to_transition's
// `sums[sums == 0.0] = 1.0`.
FullCsr to_transition(FullCsr m) {
    for (std::int64_t row = 0; row < m.n; ++row) {
        const auto begin = static_cast<std::size_t>(m.indptr[static_cast<std::size_t>(row)]);
        const auto end = static_cast<std::size_t>(m.indptr[static_cast<std::size_t>(row) + 1]);
        double sum = 0.0;
        for (std::size_t k = begin; k < end; ++k) sum += m.data[k];
        const double divisor = sum == 0.0 ? 1.0 : sum;
        for (std::size_t k = begin; k < end; ++k) m.data[k] /= divisor;
    }
    return m;
}

std::int64_t count_nonzero_rows(const FullCsr& m) {
    std::int64_t count = 0;
    for (std::int64_t row = 0; row < m.n; ++row) {
        const auto begin = static_cast<std::size_t>(m.indptr[static_cast<std::size_t>(row)]);
        const auto end = static_cast<std::size_t>(m.indptr[static_cast<std::size_t>(row) + 1]);
        double sum = 0.0;
        for (std::size_t k = begin; k < end; ++k) sum += m.data[k];
        if (sum > 0.0) ++count;
    }
    return count;
}

// C = A * B, the standard sparse matrix product (Gustavson's algorithm): for
// each row of A, accumulate into a dense row buffer via every nonzero's
// target row in B, then compact. This is what random_walk's m_input ** t
// computes (scipy's CSR __pow__ is repeated matrix multiplication, not
// elementwise), not hicx::matrix_ops's elementwise multiply().
FullCsr multiply(const FullCsr& a, const FullCsr& b) {
    const std::int64_t n = a.n;
    std::vector<double> acc(static_cast<std::size_t>(n), 0.0);
    std::vector<std::int32_t> touched;
    std::vector<std::int32_t> rows;
    std::vector<std::int32_t> cols;
    std::vector<double> vals;
    for (std::int64_t row = 0; row < n; ++row) {
        touched.clear();
        const auto a_begin = static_cast<std::size_t>(a.indptr[static_cast<std::size_t>(row)]);
        const auto a_end = static_cast<std::size_t>(a.indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t ak = a_begin; ak < a_end; ++ak) {
            const std::int32_t k = a.indices[ak];
            const double a_val = a.data[ak];
            const auto b_begin = static_cast<std::size_t>(b.indptr[static_cast<std::size_t>(k)]);
            const auto b_end = static_cast<std::size_t>(b.indptr[static_cast<std::size_t>(k) + 1]);
            for (std::size_t bk = b_begin; bk < b_end; ++bk) {
                const std::int32_t j = b.indices[bk];
                if (acc[static_cast<std::size_t>(j)] == 0.0) {
                    touched.push_back(j);
                }
                acc[static_cast<std::size_t>(j)] += a_val * b.data[bk];
            }
        }
        std::sort(touched.begin(), touched.end());
        for (const std::int32_t j : touched) {
            const double v = acc[static_cast<std::size_t>(j)];
            if (v != 0.0) {
                rows.push_back(static_cast<std::int32_t>(row));
                cols.push_back(j);
                vals.push_back(v);
            }
            acc[static_cast<std::size_t>(j)] = 0.0;
        }
    }
    return from_triplets(n, std::move(rows), std::move(cols), std::move(vals));
}

// abs(rw1 - rw2).sum(): the L1 distance over every entry of the union of the
// two sparsity patterns (a position stored on only one side counts via the
// other's implicit zero), matching scipy's sparse subtraction + abs + sum.
double l1_distance(const FullCsr& a, const FullCsr& b) {
    double total = 0.0;
    for (std::int64_t row = 0; row < a.n; ++row) {
        auto ak = static_cast<std::size_t>(a.indptr[static_cast<std::size_t>(row)]);
        const auto a_end = static_cast<std::size_t>(a.indptr[static_cast<std::size_t>(row) + 1]);
        auto bk = static_cast<std::size_t>(b.indptr[static_cast<std::size_t>(row)]);
        const auto b_end = static_cast<std::size_t>(b.indptr[static_cast<std::size_t>(row) + 1]);
        while (ak < a_end && bk < b_end) {
            if (a.indices[ak] == b.indices[bk]) {
                total += std::fabs(a.data[ak] - b.data[bk]);
                ++ak;
                ++bk;
            } else if (a.indices[ak] < b.indices[bk]) {
                total += std::fabs(a.data[ak]);
                ++ak;
            } else {
                total += std::fabs(b.data[bk]);
                ++bk;
            }
        }
        while (ak < a_end) {
            total += std::fabs(a.data[ak]);
            ++ak;
        }
        while (bk < b_end) {
            total += std::fabs(b.data[bk]);
            ++bk;
        }
    }
    return total;
}

// data_operations.subsample_to_depth_csr_upperTri: each stored count is
// resampled with numpy.random.binomial(count, target / depth), applied
// per-entry to the upper-triangle values before any normalization. Not
// numpy-equivalent (see genomedisco.hpp): std::binomial_distribution with a
// std::mt19937_64 reproduces the per-entry binomial procedure, not numpy's
// exact draws.
std::vector<double> subsample_values(const ChromosomeEntries& e, double probability,
                                     std::mt19937_64& rng) {
    std::vector<double> out(e.values.size());
    for (std::size_t i = 0; i < e.values.size(); ++i) {
        const auto count = static_cast<std::int64_t>(std::llround(e.values[i]));
        if (count <= 0 || probability >= 1.0) {
            out[i] = e.values[i];
            continue;
        }
        std::binomial_distribution<std::int64_t> dist(count, probability);
        out[i] = static_cast<double>(dist(rng));
    }
    return out;
}

double total(const ChromosomeEntries& e) {
    double s = 0.0;
    for (const double v : e.values) s += v;
    return s;
}

}  // namespace

double chromosome_score(const ChromosomeEntries& e1, const ChromosomeEntries& e2,
                        const Options& options) {
    if (e1.nbins != e2.nbins || e1.nbins <= 0 || e1.rows.empty() || e2.rows.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    std::vector<double> raw1 = e1.values;
    std::vector<double> raw2 = e2.values;
    if (options.subsample) {
        const double depth1 = total(e1);
        const double depth2 = total(e2);
        const double target = std::min(depth1, depth2);
        std::mt19937_64 rng(options.subsample_seed);
        if (depth1 > target) raw1 = subsample_values(e1, target / depth1, rng);
        if (depth2 > target) raw2 = subsample_values(e2, target / depth2, rng);
    }

    ChromosomeEntries s1 = e1;
    s1.values = raw1;
    ChromosomeEntries s2 = e2;
    s2.values = raw2;

    const std::vector<double> sq1 = sqrtvc_values(s1);
    const std::vector<double> sq2 = sqrtvc_values(s2);
    const FullCsr m1f = to_transition(symmetrize(s1, sq1));
    const FullCsr m2f = to_transition(symmetrize(s2, sq2));
    const double nonzero_total =
        0.5 * static_cast<double>(count_nonzero_rows(m1f) + count_nonzero_rows(m2f));
    if (nonzero_total == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    std::vector<double> scores;
    FullCsr rw1 = m1f;
    FullCsr rw2 = m2f;
    for (int t = 1; t <= options.tmax; ++t) {
        if (t > 1) {
            rw1 = multiply(rw1, m1f);
            rw2 = multiply(rw2, m2f);
        }
        if (t >= options.tmin) {
            scores.push_back(l1_distance(rw1, rw2) / nonzero_total);
        }
    }

    double auc;
    if (options.tmin == options.tmax || scores.size() == 1) {
        auc = scores.front();
    } else {
        // sklearn.metrics.auc(range(len(scores)), scores) / (len(scores) - 1):
        // the trapezoidal rule with unit spacing, averaged over the steps.
        double trapezoid = 0.0;
        for (std::size_t i = 1; i < scores.size(); ++i) {
            trapezoid += 0.5 * (scores[i] + scores[i - 1]);
        }
        auc = trapezoid / static_cast<double>(scores.size() - 1);
    }
    return 1.0 - auc;
}

}  // namespace hicx::genomedisco
