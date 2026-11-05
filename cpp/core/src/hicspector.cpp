#include "hicx/hicspector.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "hicx/lapack_shim.hpp"

namespace hicx::hicspector {

namespace {

// A dense, row-major (== column-major for a symmetric matrix) n x n matrix.
struct Dense {
    int n = 0;
    std::vector<double> v;

    explicit Dense(int size) : n(size), v(static_cast<std::size_t>(size) * static_cast<std::size_t>(size), 0.0) {}
    [[nodiscard]] double& at(int r, int c) { return v[static_cast<std::size_t>(r) * static_cast<std::size_t>(n) + c]; }
    [[nodiscard]] double at(int r, int c) const {
        return v[static_cast<std::size_t>(r) * static_cast<std::size_t>(n) + c];
    }
};

Dense build_symmetric(const ChromosomeEntries& e) {
    Dense m(static_cast<int>(e.nbins));
    for (std::size_t i = 0; i < e.rows.size(); ++i) {
        const int r = e.rows[i];
        const int c = e.cols[i];
        m.at(r, c) = e.values[i];
        m.at(c, r) = e.values[i];
    }
    return m;
}

std::vector<double> row_sums(const Dense& m) {
    std::vector<double> sums(static_cast<std::size_t>(m.n), 0.0);
    for (int i = 0; i < m.n; ++i) {
        double s = 0.0;
        for (int j = 0; j < m.n; ++j) s += m.at(i, j);
        sums[static_cast<std::size_t>(i)] = s;
    }
    return sums;
}

// numpy.sign(M).sum(1): the count of nonzero entries per row (Hi-C counts
// are never negative, so sign is 0 or 1).
std::vector<int> nonzero_counts(const Dense& m) {
    std::vector<int> counts(static_cast<std::size_t>(m.n), 0);
    for (int i = 0; i < m.n; ++i) {
        int c = 0;
        for (int j = 0; j < m.n; ++j) {
            if (m.at(i, j) != 0.0) ++c;
        }
        counts[static_cast<std::size_t>(i)] = c;
    }
    return counts;
}

// Restrict a dense matrix to the given index subset, both axes.
Dense restrict_to(const Dense& m, const std::vector<int>& idx) {
    Dense out(static_cast<int>(idx.size()));
    for (std::size_t i = 0; i < idx.size(); ++i) {
        for (std::size_t j = 0; j < idx.size(); ++j) {
            out.at(static_cast<int>(i), static_cast<int>(j)) =
                m.at(idx[i], idx[j]);
        }
    }
    return out;
}

// get_Laplacian: restrict to nonzero-row-sum indices, normalize by
// 1/sqrt(rowsum) on both sides, L = I - D^-1/2 A D^-1/2, symmetrized.
// `kept` receives, in order, the positions (within `m`'s own index space)
// that survived the nonzero-row-sum filter -- the index map the caller
// needs to zero-pad eigenvectors back out.
Dense laplacian(const Dense& m, std::vector<int>& kept) {
    const std::vector<double> sums = row_sums(m);
    kept.clear();
    for (int i = 0; i < m.n; ++i) {
        if (sums[static_cast<std::size_t>(i)] > 0.0) kept.push_back(i);
    }
    Dense restricted = restrict_to(m, kept);
    const int n = restricted.n;
    std::vector<double> inv_sqrt(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        inv_sqrt[static_cast<std::size_t>(i)] = 1.0 / std::sqrt(sums[static_cast<std::size_t>(kept[static_cast<std::size_t>(i)])]);
    }
    Dense lap(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double normalized = inv_sqrt[static_cast<std::size_t>(i)] * restricted.at(i, j) *
                                      inv_sqrt[static_cast<std::size_t>(j)];
            lap.at(i, j) = (i == j ? 1.0 : 0.0) - normalized;
        }
    }
    // Symmetrize: (L + L^T) / 2.
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double avg = (lap.at(i, j) + lap.at(j, i)) / 2.0;
            lap.at(i, j) = avg;
            lap.at(j, i) = avg;
        }
    }
    return lap;
}

// The `k` smallest eigenpairs of a dense symmetric matrix, via LAPACK's
// dsyevr with an ascending index range (RANGE='I', IL=1, IU=k): an exact
// dense computation, not the real tool's iterative Lanczos approximation
// (scipy.sparse.linalg.eigsh), but for the already-filtered, modest sizes
// this runs on, strictly more precise rather than a shortcut (see
// hicspector.hpp). `eigenvectors` is filled column-major-equivalent: the
// i-th returned eigenvector is `n` consecutive entries at i*n.
struct EigenResult {
    std::vector<double> values;
    std::vector<double> vectors;  // vectors[i * n + j] = eigenvector i, component j
    int n = 0;
    int found = 0;
};

EigenResult smallest_eigenpairs(Dense matrix, int k) {
    const int n = matrix.n;
    EigenResult result;
    result.n = n;
    if (n == 0 || k <= 0) {
        return result;
    }
    const int wanted = std::min(k, n);

    const char jobz = 'V';
    const char range = 'I';
    const char uplo = 'U';
    const double vl_unused = 0.0;
    const double vu_unused = 0.0;
    const int il = 1;
    const int iu = wanted;
    const char safe_minimum = 'S';
    const double abstol = 2.0 * dlamch_(&safe_minimum);
    int found = 0;
    std::vector<double> values(static_cast<std::size_t>(n));
    std::vector<double> vectors_colmajor(static_cast<std::size_t>(n) * static_cast<std::size_t>(wanted));
    std::vector<int> isuppz(static_cast<std::size_t>(2 * std::max(wanted, 1)));

    int info = 0;
    int lwork = -1;
    int liwork = -1;
    double work_query = 0.0;
    int iwork_query = 0;
    dsyevr_(&jobz, &range, &uplo, &n, matrix.v.data(), &n, &vl_unused, &vu_unused, &il, &iu, &abstol,
           &found, values.data(), vectors_colmajor.data(), &n, isuppz.data(), &work_query, &lwork,
           &iwork_query, &liwork, &info);
    if (info != 0) {
        return result;
    }
    lwork = static_cast<int>(work_query);
    liwork = iwork_query;
    std::vector<double> work(static_cast<std::size_t>(std::max(lwork, 1)));
    std::vector<int> iwork(static_cast<std::size_t>(std::max(liwork, 1)));
    dsyevr_(&jobz, &range, &uplo, &n, matrix.v.data(), &n, &vl_unused, &vu_unused, &il, &iu, &abstol,
           &found, values.data(), vectors_colmajor.data(), &n, isuppz.data(), work.data(), &lwork,
           iwork.data(), &liwork, &info);
    if (info != 0) {
        return result;
    }

    result.found = found;
    result.values.assign(values.begin(), values.begin() + found);
    result.vectors.assign(static_cast<std::size_t>(found) * static_cast<std::size_t>(n), 0.0);
    // vectors_colmajor is LAPACK column-major: eigenvector i's j-th
    // component is at vectors_colmajor[i * n + j] (lda == n), which is
    // already the layout `result.vectors` uses.
    for (int i = 0; i < found; ++i) {
        for (int j = 0; j < n; ++j) {
            result.vectors[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) + static_cast<std::size_t>(j)] =
                vectors_colmajor[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) + static_cast<std::size_t>(j)];
        }
    }
    return result;
}

double ipr(const std::vector<double>& v) {
    double sum4 = 0.0;
    for (const double x : v) sum4 += x * x * x * x;
    return sum4 == 0.0 ? 0.0 : 1.0 / sum4;
}

double evec_distance(const std::vector<double>& a, const std::vector<double>& b) {
    double d1 = 0.0, d2 = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double diff = a[i] - b[i];
        const double sum = a[i] + b[i];
        d1 += diff * diff;
        d2 += sum * sum;
    }
    return std::sqrt(std::min(d1, d2));
}

}  // namespace

double chromosome_score(const ChromosomeEntries& e1, const ChromosomeEntries& e2, const Options& options) {
    if (e1.nbins != e2.nbins || e1.nbins <= 0 || e1.rows.empty() || e2.rows.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const Dense m1 = build_symmetric(e1);
    const Dense m2 = build_symmetric(e2);
    const int n = m1.n;

    const std::vector<int> k1 = nonzero_counts(m1);
    const std::vector<int> k2 = nonzero_counts(m2);

    std::vector<int> iz;
    for (int i = 0; i < n; ++i) {
        const bool kd1 = !(k1[static_cast<std::size_t>(i)] == 1 && m1.at(i, i) > 0.0);
        const bool kd2 = !(k2[static_cast<std::size_t>(i)] == 1 && m2.at(i, i) > 0.0);
        const bool any_signal = k1[static_cast<std::size_t>(i)] + k2[static_cast<std::size_t>(i)] > 0;
        if (any_signal && kd1 && kd2) iz.push_back(i);
    }
    if (iz.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const Dense m1b = restrict_to(m1, iz);
    const Dense m2b = restrict_to(m2, iz);

    std::vector<int> kept1, kept2;  // unused beyond sizing; laplacian() recomputes kept internally
    const Dense l1 = laplacian(m1b, kept1);
    const Dense l2 = laplacian(m2b, kept2);

    const EigenResult eig1 = smallest_eigenpairs(l1, options.num_evec);
    const EigenResult eig2 = smallest_eigenpairs(l2, options.num_evec);
    if (eig1.found == 0 || eig2.found == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const int m1b_n = m1b.n;
    const int m2b_n = m2b.n;
    std::vector<std::vector<double>> b1_extend(static_cast<std::size_t>(eig1.found),
                                               std::vector<double>(static_cast<std::size_t>(m1b_n), 0.0));
    for (int i = 0; i < eig1.found; ++i) {
        for (std::size_t p = 0; p < kept1.size(); ++p) {
            b1_extend[static_cast<std::size_t>(i)][static_cast<std::size_t>(kept1[p])] =
                eig1.vectors[static_cast<std::size_t>(i) * static_cast<std::size_t>(l1.n) + p];
        }
    }
    std::vector<std::vector<double>> b2_extend(static_cast<std::size_t>(eig2.found),
                                               std::vector<double>(static_cast<std::size_t>(m2b_n), 0.0));
    for (int i = 0; i < eig2.found; ++i) {
        for (std::size_t p = 0; p < kept2.size(); ++p) {
            b2_extend[static_cast<std::size_t>(i)][static_cast<std::size_t>(kept2[p])] =
                eig2.vectors[static_cast<std::size_t>(i) * static_cast<std::size_t>(l2.n) + p];
        }
    }

    std::vector<std::vector<double>> b1_eff;
    for (const auto& v : b1_extend) {
        if (ipr(v) > options.ipr_cut) b1_eff.push_back(v);
    }
    std::vector<std::vector<double>> b2_eff;
    for (const auto& v : b2_extend) {
        if (ipr(v) > options.ipr_cut) b2_eff.push_back(v);
    }
    const int num_evec_eff = static_cast<int>(std::min(b1_eff.size(), b2_eff.size()));
    if (num_evec_eff <= 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double sum_distance = 0.0;
    for (int i = 0; i < num_evec_eff; ++i) {
        sum_distance += evec_distance(b1_eff[static_cast<std::size_t>(i)], b2_eff[static_cast<std::size_t>(i)]);
    }
    const double mean_distance = sum_distance / static_cast<double>(num_evec_eff);
    const double sqrt2 = std::sqrt(2.0);
    return std::fabs(sqrt2 - mean_distance) / sqrt2;
}

}  // namespace hicx::hicspector
