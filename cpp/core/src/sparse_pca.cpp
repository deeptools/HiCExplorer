// See hicx/sparse_pca.hpp.

#include "hicx/sparse_pca.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>

#include "hicx/numpy_compat.hpp"

namespace hicx {

std::int64_t lanczos_subspace_size(std::int64_t nev, std::int64_t n) {
    return std::min(n, std::max<std::int64_t>(2 * nev + 1, 32));
}

CovarianceOperator::CovarianceOperator(const CsrMatrix& obs_exp, int threads)
    : matrix_(&obs_exp), n_(obs_exp.rows()) {
    if (obs_exp.rows() != obs_exp.cols()) {
        throw std::runtime_error("CovarianceOperator: the matrix is not square");
    }
    if (obs_exp.symmetry() != Symmetry::Full) {
        throw std::runtime_error("CovarianceOperator: the matrix must carry both triangles");
    }
    const std::vector<std::int64_t>& indptr = obs_exp.indptr();
    const std::vector<double>& data = obs_exp.data();
    const std::size_t n = static_cast<std::size_t>(n_);
    good_.assign(n, 1);
    means_.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t begin = static_cast<std::size_t>(indptr[i]);
        const std::size_t end = static_cast<std::size_t>(indptr[i + 1]);
        bool finite = true;
        for (std::size_t k = begin; k < end; ++k) {
            finite = finite && std::isfinite(data[k]);
        }
        if (!finite) {
            good_[i] = 0;
            ++masked_;
            continue;
        }
        // The row mean over all n bins, zeros included, as np.cov's
        // X.mean(axis=1) is; the stored values are the only nonzero terms.
        means_[i] = npy::pairwise_sum(data.data() + begin, end - begin) /
                    static_cast<double>(n_);
    }

    // Split the rows into ranges of about equal stored entries. Below a few
    // hundred thousand entries a product is cheaper than waking threads.
    constexpr std::int64_t kMinimumEntriesPerPiece = 200000;
    const std::int64_t stored = indptr.empty() ? 0 : indptr.back();
    const std::int64_t pieces = std::max<std::int64_t>(
        1, std::min<std::int64_t>(threads, stored / kMinimumEntriesPerPiece));
    row_ranges_.push_back(0);
    for (std::int64_t p = 1; p < pieces; ++p) {
        const std::int64_t target = stored * p / pieces;
        const auto position = std::lower_bound(indptr.begin(), indptr.end(), target);
        const std::int64_t row =
            std::clamp<std::int64_t>(position - indptr.begin(), row_ranges_.back(), n_);
        row_ranges_.push_back(row);
    }
    row_ranges_.push_back(n_);
    if (pieces > 1) {
        pool_ = std::make_unique<WorkerPool>(static_cast<int>(pieces));
    }
    u_.assign(n, 0.0);
    w_.assign(n, 0.0);
}

void CovarianceOperator::sparse_product(const double* x, double* y, bool masked_columns,
                                        bool skip_masked_rows) const {
    const std::vector<std::int64_t>& indptr = matrix_->indptr();
    const std::vector<std::int32_t>& indices = matrix_->indices();
    const std::vector<double>& data = matrix_->data();
    const auto rows = [&](std::int64_t first, std::int64_t last) {
        for (std::int64_t i = first; i < last; ++i) {
            const std::size_t row = static_cast<std::size_t>(i);
            if (skip_masked_rows && good_[row] == 0) {
                y[row] = 0.0;
                continue;
            }
            // Sequential within the row, whichever thread runs it.
            double sum = 0.0;
            const std::size_t end = static_cast<std::size_t>(indptr[row + 1]);
            for (std::size_t k = static_cast<std::size_t>(indptr[row]); k < end; ++k) {
                const std::size_t column = static_cast<std::size_t>(indices[k]);
                if (masked_columns && good_[column] == 0) {
                    continue;
                }
                sum += data[k] * x[column];
            }
            y[row] = sum;
        }
    };
    const int pieces = static_cast<int>(row_ranges_.size()) - 1;
    if (pool_ == nullptr || pieces <= 1) {
        rows(0, n_);
        return;
    }
    pool_->run(pieces, [&](int piece) {
        rows(row_ranges_[static_cast<std::size_t>(piece)],
             row_ranges_[static_cast<std::size_t>(piece) + 1]);
    });
}

void CovarianceOperator::apply(const double* x, double* y) const {
    const std::size_t n = static_cast<std::size_t>(n_);
    // u = D x, and the centering scalar m . u, summed in index order.
    double mean_dot_u = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        u_[i] = good_[i] != 0 ? x[i] : 0.0;
        mean_dot_u += means_[i] * u_[i];
    }
    // w = X^T u = A u - (m . u) 1, with masked columns left out of A u.
    sparse_product(u_.data(), w_.data(), true, false);
    double sum_w = 0.0;
    for (std::size_t l = 0; l < n; ++l) {
        w_[l] -= mean_dot_u;
        sum_w += w_[l];
    }
    // y = D X w / (n - 1) = D (A w - m (1 . w)) / (n - 1).
    sparse_product(w_.data(), y, false, true);
    const double factor = 1.0 / static_cast<double>(n_ - 1);
    for (std::size_t i = 0; i < n; ++i) {
        y[i] = good_[i] != 0 ? (y[i] - means_[i] * sum_w) * factor : 0.0;
    }
}

namespace {

// The adapter Spectra's solvers call.
class SpectraCovariance {
  public:
    using Scalar = double;
    explicit SpectraCovariance(const CovarianceOperator& op) : op_(op) {}
    [[nodiscard]] Eigen::Index rows() const { return static_cast<Eigen::Index>(op_.size()); }
    [[nodiscard]] Eigen::Index cols() const { return static_cast<Eigen::Index>(op_.size()); }
    void perform_op(const double* x_in, double* y_out) const { op_.apply(x_in, y_out); }

  private:
    const CovarianceOperator& op_;
};

// A fixed start vector: splitmix64 from a constant seed, mapped to (-0.5,
// 0.5), masked rows zero. It depends on nothing but n and the mask.
std::vector<double> start_vector(const std::vector<unsigned char>& good) {
    std::vector<double> start(good.size(), 0.0);
    std::uint64_t state = 0x9E3779B97F4A7C15ULL;
    for (std::size_t i = 0; i < good.size(); ++i) {
        state += 0x9E3779B97F4A7C15ULL;
        std::uint64_t z = state;
        z = (z ^ (z >> 30U)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27U)) * 0x94D049BB133111EBULL;
        z = z ^ (z >> 31U);
        const double uniform = static_cast<double>(z >> 11U) * 0x1.0p-53;
        start[i] = good[i] != 0 ? uniform - 0.5 : 0.0;
    }
    return start;
}

// The sign convention of --compatMode v4: the largest magnitude component is
// positive, the first such component deciding a tie.
void fix_sign(std::vector<double>& vector) {
    std::size_t extreme = 0;
    double best = -1.0;
    for (std::size_t j = 0; j < vector.size(); ++j) {
        const double magnitude = std::fabs(vector[j]);
        if (magnitude > best) {
            best = magnitude;
            extreme = j;
        }
    }
    if (!vector.empty() && vector[extreme] < 0.0) {
        for (double& value : vector) {
            value = -value;
        }
    }
}

EigenResult dense_fallback(const CovarianceOperator& op, const std::vector<int>& which) {
    const std::int64_t n = op.size();
    DenseSymmetric covariance(n);
    std::vector<double> unit(static_cast<std::size_t>(n), 0.0);
    std::vector<double> column(static_cast<std::size_t>(n), 0.0);
    for (std::int64_t j = 0; j < n; ++j) {
        unit[static_cast<std::size_t>(j)] = 1.0;
        op.apply(unit.data(), column.data());
        unit[static_cast<std::size_t>(j)] = 0.0;
        for (std::int64_t i = 0; i < n; ++i) {
            covariance.at(i, j) = column[static_cast<std::size_t>(i)];
        }
    }
    // The operator is symmetric only to rounding; dsyevr reads one triangle.
    return leading_eigenvectors(covariance, which, EigenSolver::Dsyevr);
}

}  // namespace

EigenResult covariance_eigenvectors_lanczos(const CsrMatrix& obs_exp, const std::vector<int>& which,
                                            int threads, LanczosDiagnostics* diagnostics) {
    const std::int64_t n = obs_exp.rows();
    EigenResult result;
    result.values.assign(which.size(), 0.0);
    result.vectors.assign(which.size(), {});
    LanczosDiagnostics local;
    LanczosDiagnostics& report = diagnostics != nullptr ? *diagnostics : local;
    report = LanczosDiagnostics();

    std::int64_t nev = 0;
    for (const int index : which) {
        if (index >= 1 && index <= n) {
            nev = std::max<std::int64_t>(nev, index);
        }
    }
    if (nev == 0) {
        return result;
    }
    if (n == 1) {
        // np.cov of one observation is 0 * inf, NaN, which hicPCA zeroes; the
        // eigenvector of the 1 by 1 zero matrix is the unit vector.
        for (std::size_t entry = 0; entry < which.size(); ++entry) {
            if (which[entry] == 1) {
                result.vectors[entry] = {1.0};
            }
        }
        report.dense_fallback = true;
        return result;
    }

    const CovarianceOperator op(obs_exp, threads);
    report.masked_rows = op.masked_rows();
    if (n <= kLanczosDenseLimit || nev >= n || op.masked_rows() == n) {
        report.dense_fallback = true;
        return dense_fallback(op, which);
    }

    SpectraCovariance adapter(op);
    const std::vector<double> start = start_vector(op.mask());
    std::int64_t ncv = lanczos_subspace_size(nev, n);
    constexpr Eigen::Index kMaxRestarts = 1000;
    for (int attempt = 0; attempt < 2; ++attempt) {
        Spectra::SymEigsSolver<SpectraCovariance> solver(adapter, static_cast<Eigen::Index>(nev),
                                                        static_cast<Eigen::Index>(ncv));
        solver.init(start.data());
        solver.compute(Spectra::SortRule::LargestAlge, kMaxRestarts, kLanczosTolerance,
                       Spectra::SortRule::LargestAlge);
        report.ncv = ncv;
        report.iterations += solver.num_iterations();
        report.products += solver.num_operations();
        if (solver.info() == Spectra::CompInfo::Successful) {
            const Eigen::VectorXd values = solver.eigenvalues();
            const Eigen::MatrixXd vectors = solver.eigenvectors();
            for (std::size_t entry = 0; entry < which.size(); ++entry) {
                const int index = which[entry] - 1;
                if (index < 0 || index >= n) {
                    continue;
                }
                result.values[entry] = values(index);
                std::vector<double> vector(static_cast<std::size_t>(n));
                for (std::int64_t i = 0; i < n; ++i) {
                    vector[static_cast<std::size_t>(i)] = vectors(i, index);
                }
                fix_sign(vector);
                result.vectors[entry] = std::move(vector);
            }
            return result;
        }
        // A larger Krylov space for the second and last attempt.
        ncv = std::min(n, 2 * ncv);
        if (ncv <= nev) {
            break;
        }
    }
    throw std::runtime_error(
        "the Lanczos eigensolver did not converge on a chromosome of " + std::to_string(n) +
        " bins; --eigenSolver dense computes the same eigenvectors with a dense matrix");
}

}  // namespace hicx
