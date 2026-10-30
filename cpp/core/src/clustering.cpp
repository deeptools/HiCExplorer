#include "hicx/clustering.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>
#include <string>
#include <utility>

#include "hicx/numpy_compat.hpp"
#include "hicx/numpy_random.hpp"
#include "hicx/numpy_sort.hpp"

// The OpenBLAS build numpy, scipy and scikit-learn are linked against, loaded
// on first use rather than linked. OpenBLAS starts its worker threads when the
// library is loaded and they spin before main() could lower their number:
// linked, it cost hicAggregateContacts 1.2 s of CPU on 32 cores for a run that
// needs 0.05 s of wall time, and every run paid it, clustering or not.
// OPENBLAS_NUM_THREADS is read at load, so it is set first (unless the user
// set it) and the thread count is pinned to one afterwards either way, which
// is also what keeps the results independent of the machine.
// HICX_OPENBLAS_LIBRARY is the library CMake found; the bare file name is
// tried first so that the executable's runpath is honoured.

#include <cstdlib>
#include <dlfcn.h>

namespace {

struct Blas {
    void (*dgemv)(int, int, int, int, double, const double*, int, const double*, int, double,
                  double*, int) = nullptr;
    void (*dgemm_rowmajor)(int, int, int, int, int, int, double, const double*, int,
                           const double*, int, double, double*, int) = nullptr;
    double (*ddot)(int, const double*, int, const double*, int) = nullptr;
    void (*dgemm)(const char*, const char*, const int*, const int*, const int*, const double*,
                  const double*, const int*, const double*, const int*, const double*, double*,
                  const int*) = nullptr;
};

const Blas& blas() {
    static const Blas loaded = [] {
        setenv("OPENBLAS_NUM_THREADS", "1", 0);
        const std::string full = HICX_OPENBLAS_LIBRARY;
        const std::string base = full.substr(full.find_last_of('/') + 1);
        void* handle = dlopen(base.c_str(), RTLD_NOW | RTLD_LOCAL);
        if (handle == nullptr) {
            handle = dlopen(full.c_str(), RTLD_NOW | RTLD_LOCAL);
        }
        if (handle == nullptr) {
            throw std::runtime_error("cannot load the BLAS library " + full + ": " +
                                     std::string(dlerror()));
        }
        const auto symbol = [handle](const char* name) {
            void* address = dlsym(handle, name);
            if (address == nullptr) {
                throw std::runtime_error(std::string("the BLAS library lacks ") + name);
            }
            return address;
        };
        Blas table;
        table.dgemv = reinterpret_cast<decltype(table.dgemv)>(symbol("cblas_dgemv"));
        table.dgemm_rowmajor = reinterpret_cast<decltype(table.dgemm_rowmajor)>(symbol("cblas_dgemm"));
        table.ddot = reinterpret_cast<decltype(table.ddot)>(symbol("cblas_ddot"));
        table.dgemm = reinterpret_cast<decltype(table.dgemm)>(symbol("dgemm_"));
        if (void* set_threads = dlsym(handle, "openblas_set_num_threads")) {
            reinterpret_cast<void (*)(int)>(set_threads)(1);
        }
        return table;
    }();
    return loaded;
}

void cblas_dgemv(int order, int trans, int m, int n, double alpha, const double* a, int lda,
                 const double* x, int incx, double beta, double* y, int incy) {
    blas().dgemv(order, trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}
void cblas_dgemm(int order, int transa, int transb, int m, int n, int k, double alpha,
                 const double* a, int lda, const double* b, int ldb, double beta, double* c,
                 int ldc) {
    blas().dgemm_rowmajor(order, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}
double cblas_ddot(int n, const double* x, int incx, const double* y, int incy) {
    return blas().ddot(n, x, incx, y, incy);
}
void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
            const double* alpha, const double* a, const int* lda, const double* b,
            const int* ldb, const double* beta, double* c, const int* ldc) {
    blas().dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

}  // namespace

namespace hicx::cluster {

namespace {

constexpr int kCblasRowMajor = 101;
constexpr int kCblasColMajor = 102;
constexpr int kCblasNoTrans = 111;
constexpr int kCblasTrans = 112;

// Samples per block of the assignment step. The centre sums of a block are
// formed on their own and the blocks are added in order, so the block size
// fixes the rounding of the sums; 256 is the value that reproduces
// scikit-learn 1.3.2's bits (test_sklearn_oracle.cpp, datasets of 255 to 700
// samples).
constexpr std::int64_t kBlockSamples = 256;
// scikit-learn 1.3 defaults of KMeans: max_iter=300, n_init=10, tol=1e-4.
constexpr int kMaxIterations = 300;
constexpr int kRestarts = 10;
constexpr double kRelativeTolerance = 1e-4;

int as_int(std::int64_t value) { return static_cast<int>(value); }

// Squared Euclidean distance, accumulated four coordinates at a time (the
// partial sum of a group of four is formed first and then added to the
// running total); the tail is added one coordinate at a time.
double squared_distance(const double* a, const double* b, std::int64_t n_features) {
    const std::int64_t groups = n_features / 4;
    const std::int64_t tail = n_features % 4;
    double total = 0.0;
    for (std::int64_t g = 0; g < groups; ++g) {
        total += ((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                  (a[2] - b[2]) * (a[2] - b[2]) + (a[3] - b[3]) * (a[3] - b[3]));
        a += 4;
        b += 4;
    }
    for (std::int64_t i = 0; i < tail; ++i) {
        total += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return total;
}

std::vector<double> row_norms(const double* values, std::int64_t rows, std::int64_t features) {
    std::vector<double> norms(static_cast<std::size_t>(rows));
    for (std::int64_t i = 0; i < rows; ++i) {
        norms[static_cast<std::size_t>(i)] =
            detail::einsum_row_norm_squared(values + i * features, features);
    }
    return norms;
}

// Sum of one feature over the samples. numpy reduces a (n, d) array along
// axis 0 by adding row after row (sequentially) when d > 1; with a single
// feature the data are one contiguous run and it uses the pairwise sum.
double sum_over_samples(const std::vector<double>& column, std::int64_t n_features) {
    if (n_features == 1) {
        return npy::pairwise_sum(column.data(), column.size());
    }
    double total = 0.0;
    for (const double value : column) {
        total += value;
    }
    return total;
}

std::vector<double> feature_column(const Samples& X, std::int64_t j) {
    std::vector<double> column(static_cast<std::size_t>(X.samples));
    for (std::int64_t i = 0; i < X.samples; ++i) {
        column[static_cast<std::size_t>(i)] = X.row(i)[j];
    }
    return column;
}

// The convergence threshold of KMeans: tol times the mean over the features of
// the per-feature variance (documented for the `tol` parameter). The variance
// of a column is summed sequentially over the samples; the mean over the
// features is numpy's pairwise sum divided by the count.
double convergence_threshold(const Samples& X) {
    const std::int64_t n = X.samples;
    const std::int64_t d = X.features;
    std::vector<double> variances(static_cast<std::size_t>(d));
    for (std::int64_t j = 0; j < d; ++j) {
        std::vector<double> column = feature_column(X, j);
        const double mean = sum_over_samples(column, d) / static_cast<double>(n);
        for (double& value : column) {
            const double deviation = value - mean;
            value = deviation * deviation;
        }
        variances[static_cast<std::size_t>(j)] = sum_over_samples(column, d) / static_cast<double>(n);
    }
    const double mean_variance =
        npy::pairwise_sum(variances.data(), variances.size()) / static_cast<double>(d);
    return mean_variance * kRelativeTolerance;
}

// Squared distances from every sample to one point through the expansion
// |x|^2 - 2 x.c + |c|^2, clipped at zero, as scikit-learn evaluates them in
// the seeding. Rounding of the expansion decides ties between candidates, so
// the order of the three terms is kept.
inline double expanded_distance(double dot, double point_norm, double sample_norm) {
    double value = -2.0 * dot;
    value += point_norm;
    value += sample_norm;
    return std::max(value, 0.0);
}

// k-means++ seeding (Arthur and Vassilvitskii 2007, SODA, with the greedy
// variant of several local trials per centre that scikit-learn documents: the
// number of trials is 2 + int(ln k)). The first centre is drawn uniformly. Each
// further centre is chosen among the trial candidates, drawn with probability
// proportional to the squared distance to the nearest centre so far, as the
// one that gives the smallest total squared distance.
std::vector<double> kmeans_plusplus_centers(const Samples& X,
                                            const std::vector<double>& x_squared_norms,
                                            std::int64_t n_clusters,
                                            npy::RandomState& random_state) {
    const std::int64_t n = X.samples;
    const std::int64_t d = X.features;
    std::vector<double> centers(static_cast<std::size_t>(n_clusters * d));
    const std::int64_t n_trials =
        2 + static_cast<std::int64_t>(std::log(static_cast<double>(n_clusters)));

    // Uniform first centre: choice(n, p=1/n).
    std::vector<double> uniform(static_cast<std::size_t>(n), 1.0 / static_cast<double>(n));
    const std::int64_t first = random_state.choice(uniform);
    std::copy(X.row(first), X.row(first) + d, centers.begin());

    // Dot products of every sample with the first centre: numpy evaluates a
    // row vector times X^T with a matrix-vector product on X in column-major
    // order.
    const std::vector<double> ones(static_cast<std::size_t>(n), 1.0);
    std::vector<double> nearest(static_cast<std::size_t>(n));
    cblas_dgemv(kCblasColMajor, kCblasTrans, as_int(d), as_int(n), 1.0, X.values.data(),
                as_int(d), centers.data(), 1, 0.0, nearest.data(), 1);
    {
        const double centre_norm = detail::einsum_row_norm_squared(centers.data(), d);
        for (std::int64_t j = 0; j < n; ++j) {
            nearest[static_cast<std::size_t>(j)] = expanded_distance(
                nearest[static_cast<std::size_t>(j)], centre_norm,
                x_squared_norms[static_cast<std::size_t>(j)]);
        }
    }
    // Total squared distance: a (1, n) times (n,) product is numpy's dot of two
    // vectors, a BLAS dot added to 0.0.
    double potential = 0.0 + cblas_ddot(as_int(n), nearest.data(), 1, ones.data(), 1);

    std::vector<double> candidates(static_cast<std::size_t>(n_trials * d));
    std::vector<double> trial_distances(static_cast<std::size_t>(n_trials * n));
    std::vector<double> trial_potentials(static_cast<std::size_t>(n_trials));
    std::vector<double> cumulative(static_cast<std::size_t>(n));
    std::vector<std::int64_t> candidate_ids(static_cast<std::size_t>(n_trials));

    for (std::int64_t c = 1; c < n_clusters; ++c) {
        // Inverse-transform sampling on the running sum of the squared
        // distances: a uniform draw scaled by the potential, located by a
        // left-side binary search.
        std::vector<double> targets = random_state.uniform(static_cast<std::size_t>(n_trials));
        for (double& value : targets) {
            value *= potential;
        }
        double running = 0.0;
        for (std::int64_t j = 0; j < n; ++j) {
            running += 1.0 * nearest[static_cast<std::size_t>(j)];
            cumulative[static_cast<std::size_t>(j)] = running;
        }
        for (std::int64_t t = 0; t < n_trials; ++t) {
            auto position = static_cast<std::int64_t>(
                std::lower_bound(cumulative.begin(), cumulative.end(),
                                 targets[static_cast<std::size_t>(t)]) -
                cumulative.begin());
            position = std::min(position, n - 1);  // rounding can leave the last sum short
            candidate_ids[static_cast<std::size_t>(t)] = position;
            std::copy(X.row(position), X.row(position) + d, candidates.begin() + t * d);
        }
        // Candidates times X^T: a matrix-matrix product, row major, second
        // operand transposed.
        cblas_dgemm(kCblasRowMajor, kCblasNoTrans, kCblasTrans, as_int(n_trials), as_int(n),
                    as_int(d), 1.0, candidates.data(), as_int(d), X.values.data(), as_int(d),
                    0.0, trial_distances.data(), as_int(n));
        for (std::int64_t t = 0; t < n_trials; ++t) {
            const double candidate_norm =
                detail::einsum_row_norm_squared(candidates.data() + t * d, d);
            for (std::int64_t j = 0; j < n; ++j) {
                double& value = trial_distances[static_cast<std::size_t>(t * n + j)];
                value = expanded_distance(value, candidate_norm,
                                          x_squared_norms[static_cast<std::size_t>(j)]);
                // distance to the nearer of the old centres and the candidate
                const double current = nearest[static_cast<std::size_t>(j)];
                value = current <= value ? current : value;
            }
        }
        // Potential of every trial: matrix times a column of ones, a
        // column-major matrix-vector product.
        cblas_dgemv(kCblasColMajor, kCblasTrans, as_int(n), as_int(n_trials), 1.0,
                    trial_distances.data(), as_int(n), ones.data(), 1, 0.0,
                    trial_potentials.data(), 1);
        std::int64_t best = 0;
        for (std::int64_t t = 1; t < n_trials; ++t) {
            if (trial_potentials[static_cast<std::size_t>(t)] <
                trial_potentials[static_cast<std::size_t>(best)]) {
                best = t;  // the first of equal potentials wins
            }
        }
        potential = trial_potentials[static_cast<std::size_t>(best)];
        std::copy(trial_distances.begin() + best * n, trial_distances.begin() + (best + 1) * n,
                  nearest.begin());
        const std::int64_t chosen = candidate_ids[static_cast<std::size_t>(best)];
        std::copy(X.row(chosen), X.row(chosen) + d, centers.begin() + c * d);
    }
    return centers;
}

// A cluster that lost all its samples in an update is given the sample that is
// farthest from its (old) centre, taking the farthest samples in decreasing
// order of distance for several empty clusters; the sample leaves its previous
// cluster. The partial selection of the farthest samples is numpy's
// argpartition (numpy_sort.hpp), whose order among equal distances is part of
// the result.
void reseat_empty_clusters(const Samples& X, const std::vector<double>& centers_old,
                           std::vector<double>& sums, std::vector<double>& counts,
                           const std::vector<std::int32_t>& labels, std::int64_t n_clusters) {
    std::vector<std::int64_t> empty;
    for (std::int64_t j = 0; j < n_clusters; ++j) {
        if (counts[static_cast<std::size_t>(j)] == 0.0) {
            empty.push_back(j);
        }
    }
    if (empty.empty()) {
        return;
    }
    const std::int64_t n = X.samples;
    const std::int64_t d = X.features;
    const auto n_empty = static_cast<std::int64_t>(empty.size());
    // Squared distance of every sample to its centre, summed over the
    // coordinates with numpy's pairwise sum.
    std::vector<double> distances(static_cast<std::size_t>(n));
    std::vector<double> squares(static_cast<std::size_t>(d));
    for (std::int64_t i = 0; i < n; ++i) {
        const double* center = centers_old.data() + labels[static_cast<std::size_t>(i)] * d;
        for (std::int64_t k = 0; k < d; ++k) {
            const double difference = X.row(i)[k] - center[k];
            squares[static_cast<std::size_t>(k)] = difference * difference;
        }
        distances[static_cast<std::size_t>(i)] = npy::pairwise_sum(squares.data(), squares.size());
    }
    const std::vector<std::int64_t> partition = npy::argpartition(distances, -n_empty);
    for (std::int64_t idx = 0; idx < n_empty; ++idx) {
        const std::int64_t new_cluster = empty[static_cast<std::size_t>(idx)];
        const std::int64_t far = partition[static_cast<std::size_t>(n - 1 - idx)];
        const std::int64_t old_cluster = labels[static_cast<std::size_t>(far)];
        for (std::int64_t k = 0; k < d; ++k) {
            sums[static_cast<std::size_t>(old_cluster * d + k)] -= X.row(far)[k];
            sums[static_cast<std::size_t>(new_cluster * d + k)] = X.row(far)[k];
        }
        counts[static_cast<std::size_t>(new_cluster)] = 1.0;
        counts[static_cast<std::size_t>(old_cluster)] -= 1.0;
    }
}

// One pass of Lloyd's algorithm (Lloyd 1982, IEEE Trans. Inf. Theory 28(2)):
// assign every sample to its nearest centre, then, when `update_centers` is
// set, move each centre to the mean of its samples and report how far each
// moved. The assignment works on blocks of kBlockSamples samples: the squared
// distances of a block to all centres come from |c|^2 - 2 x.c (|x|^2 does not
// change the nearest centre) with one matrix product, the nearest centre is the
// first of equal minima, and the block's per-centre sums and counts are added to
// the totals in block order.
void lloyd_pass(const Samples& X, const std::vector<double>& centers_old,
                std::vector<double>& centers_new, std::vector<double>& counts,
                std::vector<std::int32_t>& labels, std::vector<double>& center_shift,
                std::int64_t n_clusters, bool update_centers) {
    const std::int64_t n = X.samples;
    const std::int64_t d = X.features;
    const std::int64_t block = n > kBlockSamples ? kBlockSamples : n;
    const std::int64_t last_block_size = n % block;
    const std::int64_t blocks = n / block + (last_block_size != 0 ? 1 : 0);
    const std::vector<double> centre_norms = row_norms(centers_old.data(), n_clusters, d);

    if (update_centers) {
        std::fill(centers_new.begin(), centers_new.end(), 0.0);
        std::fill(counts.begin(), counts.end(), 0.0);
    }
    std::vector<double> block_sums(static_cast<std::size_t>(n_clusters * d));
    std::vector<double> block_counts(static_cast<std::size_t>(n_clusters));
    std::vector<double> scores(static_cast<std::size_t>(block * n_clusters));

    for (std::int64_t b = 0; b < blocks; ++b) {
        const std::int64_t start = b * block;
        const std::int64_t m = (b == blocks - 1 && last_block_size > 0) ? last_block_size : block;
        std::fill(block_sums.begin(), block_sums.end(), 0.0);
        std::fill(block_counts.begin(), block_counts.end(), 0.0);

        for (std::int64_t i = 0; i < m; ++i) {
            for (std::int64_t j = 0; j < n_clusters; ++j) {
                scores[static_cast<std::size_t>(i * n_clusters + j)] =
                    centre_norms[static_cast<std::size_t>(j)];
            }
        }
        // scores += -2 X_block C^T, a row-major product evaluated as the
        // column-major product C X_block^T of the same memory.
        const char trans_centers = 't';
        const char trans_samples = 'n';
        const int k_int = as_int(n_clusters);
        const int m_int = as_int(m);
        const int d_int = as_int(d);
        const double alpha = -2.0;
        const double beta = 1.0;
        dgemm_(&trans_centers, &trans_samples, &k_int, &m_int, &d_int, &alpha,
               centers_old.data(), &d_int, X.row(start), &d_int, &beta, scores.data(), &k_int);

        for (std::int64_t i = 0; i < m; ++i) {
            double least = scores[static_cast<std::size_t>(i * n_clusters)];
            std::int32_t label = 0;
            for (std::int64_t j = 1; j < n_clusters; ++j) {
                const double score = scores[static_cast<std::size_t>(i * n_clusters + j)];
                if (score < least) {
                    least = score;
                    label = static_cast<std::int32_t>(j);
                }
            }
            labels[static_cast<std::size_t>(start + i)] = label;
            if (update_centers) {
                block_counts[static_cast<std::size_t>(label)] += 1.0;
                const double* sample = X.row(start + i);
                for (std::int64_t k = 0; k < d; ++k) {
                    block_sums[static_cast<std::size_t>(label * d + k)] += sample[k] * 1.0;
                }
            }
        }
        if (update_centers) {
            for (std::int64_t j = 0; j < n_clusters; ++j) {
                counts[static_cast<std::size_t>(j)] += block_counts[static_cast<std::size_t>(j)];
                for (std::int64_t k = 0; k < d; ++k) {
                    centers_new[static_cast<std::size_t>(j * d + k)] +=
                        block_sums[static_cast<std::size_t>(j * d + k)];
                }
            }
        }
    }
    if (!update_centers) {
        return;
    }
    reseat_empty_clusters(X, centers_old, centers_new, counts, labels, n_clusters);
    // Sums to means: multiplication by the reciprocal of the count.
    for (std::int64_t j = 0; j < n_clusters; ++j) {
        const double count = counts[static_cast<std::size_t>(j)];
        if (count > 0) {
            const double reciprocal = 1.0 / count;
            for (std::int64_t k = 0; k < d; ++k) {
                centers_new[static_cast<std::size_t>(j * d + k)] *= reciprocal;
            }
        }
    }
    for (std::int64_t j = 0; j < n_clusters; ++j) {
        center_shift[static_cast<std::size_t>(j)] = std::sqrt(squared_distance(
            centers_new.data() + j * d, centers_old.data() + j * d, d));
    }
}

struct LloydResult {
    std::vector<std::int32_t> labels;
    std::vector<double> centers;
    double inertia = 0.0;
};

// Lloyd iterations from the given seeds until the labels stop changing or the
// summed squared centre movement is at most `tol` (the documented stopping
// rules of KMeans), at most kMaxIterations passes. When it stops on the
// tolerance, one more assignment pass makes the labels consistent with the
// final centres.
LloydResult lloyd_run(const Samples& X, std::vector<double> centers, std::int64_t n_clusters,
                      double tol) {
    const std::int64_t n = X.samples;
    const std::int64_t d = X.features;
    std::vector<double> centers_new(centers.size(), 0.0);
    std::vector<std::int32_t> labels(static_cast<std::size_t>(n), -1);
    std::vector<std::int32_t> labels_old = labels;
    std::vector<double> counts(static_cast<std::size_t>(n_clusters), 0.0);
    std::vector<double> center_shift(static_cast<std::size_t>(n_clusters), 0.0);

    bool labels_settled = false;
    for (int iteration = 0; iteration < kMaxIterations; ++iteration) {
        lloyd_pass(X, centers, centers_new, counts, labels, center_shift, n_clusters, true);
        std::swap(centers, centers_new);
        if (labels == labels_old) {
            labels_settled = true;
            break;
        }
        std::vector<double> squared(center_shift.size());
        for (std::size_t j = 0; j < center_shift.size(); ++j) {
            squared[j] = center_shift[j] * center_shift[j];
        }
        if (npy::pairwise_sum(squared.data(), squared.size()) <= tol) {
            break;
        }
        labels_old = labels;
    }
    if (!labels_settled) {
        std::vector<double> unused;
        lloyd_pass(X, centers, unused, counts, labels, center_shift, n_clusters, false);
    }
    LloydResult result;
    for (std::int64_t i = 0; i < n; ++i) {
        result.inertia += squared_distance(
            X.row(i), centers.data() + labels[static_cast<std::size_t>(i)] * d, d) * 1.0;
    }
    result.labels = std::move(labels);
    result.centers = std::move(centers);
    return result;
}

// True when the two label vectors describe the same partition up to renaming
// of the clusters (every cluster of the first maps onto one of the second).
bool same_partition(const std::vector<std::int32_t>& labels1,
                    const std::vector<std::int32_t>& labels2, std::int64_t n_clusters) {
    std::vector<std::int32_t> mapping(static_cast<std::size_t>(n_clusters), -1);
    for (std::size_t i = 0; i < labels1.size(); ++i) {
        std::int32_t& mapped = mapping[static_cast<std::size_t>(labels1[i])];
        if (mapped == -1) {
            mapped = labels2[i];
        } else if (mapped != labels2[i]) {
            return false;
        }
    }
    return true;
}

// Binary min-heap of integers in an array (Williams 1964, CACM 7(6); Cormen
// et al., Introduction to Algorithms, ch. 6): the children of position p are
// 2p + 1 and 2p + 2 and no child is smaller than its parent. hc_cut reads the
// cluster labels off the array positions, so the exact layout matters, not only
// the heap property; the layout is that of Python's heapq module:
//  * push appends and moves the item up while it is smaller than its parent;
//  * replacing the root moves the hole down along the smaller child (the right
//    child on equal keys) to a leaf and then moves the displaced item up from
//    there (the bottom-up variant of Floyd 1964, Knuth TAOCP vol. 3, 5.2.3).
// Checked against heapq black box through the label order of
// AgglomerativeClustering in test_sklearn_oracle.cpp.
void heap_move_up(std::vector<std::int64_t>& heap, std::size_t root, std::size_t pos) {
    const std::int64_t item = heap[pos];
    while (pos > root) {
        const std::size_t parent_pos = (pos - 1) >> 1;
        if (!(item < heap[parent_pos])) {
            break;
        }
        heap[pos] = heap[parent_pos];
        pos = parent_pos;
    }
    heap[pos] = item;
}

void heap_replace_at(std::vector<std::int64_t>& heap, std::size_t pos) {
    const std::size_t size = heap.size();
    const std::size_t origin = pos;
    const std::int64_t item = heap[pos];
    std::size_t child = 2 * pos + 1;
    while (child < size) {
        if (child + 1 < size && !(heap[child] < heap[child + 1])) {
            child += 1;
        }
        heap[pos] = heap[child];
        pos = child;
        child = 2 * pos + 1;
    }
    heap[pos] = item;
    heap_move_up(heap, origin, pos);
}

void heap_push(std::vector<std::int64_t>& heap, std::int64_t item) {
    heap.push_back(item);
    heap_move_up(heap, 0, heap.size() - 1);
}

// Offers `item`: when it is larger than the smallest element it takes the
// root's place; otherwise the heap is unchanged.
void heap_push_pop(std::vector<std::int64_t>& heap, std::int64_t item) {
    if (!heap.empty() && heap[0] < item) {
        std::swap(item, heap[0]);
        heap_replace_at(heap, 0);
    }
}

// Disjoint sets over the 2n - 1 nodes of a dendrogram (Mullner 2011, section
// 3.2). Leaves are 0..n-1; the cluster made by the i-th merge is node n + i and
// becomes the parent of the roots it joins. Finding a root compresses the path.
struct DendrogramSets {
    explicit DendrogramSets(std::int64_t leaves)
        : parent(static_cast<std::size_t>(2 * leaves - 1)),
          members(static_cast<std::size_t>(2 * leaves - 1), 1), next_node(leaves) {
        std::iota(parent.begin(), parent.end(), std::int64_t{0});
    }
    std::int64_t root(std::int64_t node) {
        std::int64_t top = node;
        while (parent[static_cast<std::size_t>(top)] != top) {
            top = parent[static_cast<std::size_t>(top)];
        }
        while (parent[static_cast<std::size_t>(node)] != top) {
            const std::int64_t up = parent[static_cast<std::size_t>(node)];
            parent[static_cast<std::size_t>(node)] = top;
            node = up;
        }
        return top;
    }
    // Joins two roots under a new node and returns the number of leaves below it.
    std::int64_t join(std::int64_t first, std::int64_t second) {
        const std::int64_t count =
            members[static_cast<std::size_t>(first)] + members[static_cast<std::size_t>(second)];
        parent[static_cast<std::size_t>(first)] = next_node;
        parent[static_cast<std::size_t>(second)] = next_node;
        members[static_cast<std::size_t>(next_node)] = count;
        next_node += 1;
        return count;
    }
    std::vector<std::int64_t> parent;
    std::vector<std::int64_t> members;
    std::int64_t next_node;
};

std::int64_t condensed_index(std::int64_t n, std::int64_t i, std::int64_t j) {
    if (i < j) {
        return n * i - (i * (i + 1) / 2) + (j - i - 1);
    }
    return n * j - (j * (j + 1) / 2) + (i - j - 1);
}

// Lance-Williams update for Ward linkage (Mullner 2011, table 1, Ward row), for
// distances that are not squared: with cluster k and the merged pair x, y,
//   d(k, x+y)^2 = ((n_k + n_x) d(k,x)^2 + (n_k + n_y) d(k,y)^2 - n_k d(x,y)^2)
//                 / (n_k + n_x + n_y).
// The common denominator is applied as one reciprocal factor to each term.
double ward_update(double d_xi, double d_yi, double d_xy, int size_x, int size_y, int size_i) {
    const double reciprocal = 1.0 / (size_x + size_y + size_i);
    return std::sqrt((size_i + size_x) * reciprocal * d_xi * d_xi +
                     (size_i + size_y) * reciprocal * d_yi * d_yi -
                     size_i * reciprocal * d_xy * d_xy);
}

void require_finite(const Samples& X, const char* estimator) {
    for (const double value : X.values) {
        if (!std::isfinite(value)) {
            throw ClusteringError(std::string("Input X contains ") +
                                  (std::isnan(value) ? "NaN." : "infinity or a value too large for dtype('float64').") +
                                  " (" + estimator + ")");
        }
    }
}

}  // namespace

namespace detail {

double einsum_row_norm_squared(const double* row, std::int64_t features) {
    // Sum of squares of one row as numpy.einsum('ij,ij->i') evaluates it on a
    // CPU with two float64 lanes (the SSE baseline numpy dispatches to). Even
    // positions accumulate in one lane and odd positions in the other. A
    // multiply-add is formed as x*x + lane (not an FMA; numpy's baseline has
    // none). Full blocks of eight elements add their four pairs from the
    // last pair to the first; the remaining elements are taken two at a time,
    // an odd last element paired with zero. The lanes are added at the end.
    // The summation order was found by comparison with einsum on random rows
    // of every length from 1 to 130 (test_sklearn_oracle.cpp and
    // test_aggregate_support.cpp).
    double lane0 = 0.0;
    double lane1 = 0.0;
    std::int64_t remaining = features;
    const double* x = row;
    for (; remaining >= 8; remaining -= 8, x += 8) {
        for (int pair = 3; pair >= 0; --pair) {
            lane0 = x[2 * pair] * x[2 * pair] + lane0;
            lane1 = x[2 * pair + 1] * x[2 * pair + 1] + lane1;
        }
    }
    for (; remaining > 0; remaining -= 2, x += 2) {
        lane0 = x[0] * x[0] + lane0;
        const double second = remaining >= 2 ? x[1] : 0.0;
        lane1 = second * second + lane1;
    }
    return 0.0 + (lane0 + lane1);
}

namespace {

enum class LinkageMethod { Ward, Complete };

std::vector<double> nn_chain_linkage(const Samples& X, LinkageMethod method);

}  // namespace

std::vector<double> ward_linkage(const Samples& X) {
    return nn_chain_linkage(X, LinkageMethod::Ward);
}

std::vector<double> complete_linkage(const Samples& X) {
    if (X.samples < 2) {
        throw ClusteringError("The number of observations cannot be determined on an empty "
                              "distance matrix. (scipy.cluster.hierarchy.linkage needs at least "
                              "two observations)");
    }
    require_finite(X, "linkage");
    return nn_chain_linkage(X, LinkageMethod::Complete);
}

std::vector<std::int64_t> dendrogram_leaves(const std::vector<double>& Z, std::int64_t n) {
    std::vector<std::int64_t> leaves;
    if (n <= 0) {
        return leaves;
    }
    if (n == 1) {
        leaves.push_back(0);
        return leaves;
    }
    std::vector<std::int64_t> stack = {2 * n - 2};
    while (!stack.empty()) {
        const std::int64_t id = stack.back();
        stack.pop_back();
        if (id < n) {
            leaves.push_back(id);
            continue;
        }
        const auto row = static_cast<std::size_t>(id - n);
        // The second child goes on the stack first, so the first is visited first.
        stack.push_back(static_cast<std::int64_t>(Z[row * 4 + 1]));
        stack.push_back(static_cast<std::int64_t>(Z[row * 4]));
    }
    return leaves;
}

namespace {

// Agglomerative clustering by the nearest-neighbour-chain algorithm (Murtagh
// 1983; Mullner 2011, "Modern hierarchical, agglomerative clustering
// algorithms", arXiv:1109.2378, section 3), for linkages that satisfy the
// reducibility property (Ward and complete).
//
// Phase 1 grows a chain of clusters in which every cluster is the nearest
// neighbour of its predecessor. When the last two clusters are nearest
// neighbours of each other they are merged by the Lance-Williams update, the
// two chain entries are dropped and the search continues from the shorter chain.
// An empty chain restarts at the lowest-numbered live cluster.
//
// Ties, in the order the tools observe them: a cluster's nearest neighbour is
// the live cluster at the smallest distance, a later cluster replacing an
// earlier one only by being strictly closer, so the lowest index wins; the
// predecessor in the chain wins over every other cluster at equal distance
// (this is what ends the chain); the merged cluster occupies the slot of the
// higher-numbered member and the lower slot is retired.
//
// Phase 2 (section 3.2): the merges come out in chain order, not by height, so
// they are sorted by height with a stable sort (equal heights keep their chain
// order) and the clusters are renamed in that order. The cluster made by the
// i-th sorted merge is called n + i, and each row lists its two children with
// the smaller name first.
std::vector<double> nn_chain_linkage(const Samples& X, LinkageMethod method) {
    const std::int64_t n = X.samples;
    const std::int64_t dims = X.features;
    const auto sz = [](std::int64_t v) { return static_cast<std::size_t>(v); };

    // Euclidean distances between all pairs of observations, row-major upper triangle.
    std::vector<double> dist(sz(n * (n - 1) / 2));
    {
        std::size_t at = 0;
        for (std::int64_t i = 0; i < n; ++i) {
            for (std::int64_t j = i + 1; j < n; ++j) {
                double squares = 0.0;
                const double* p = X.row(i);
                const double* q = X.row(j);
                for (std::int64_t k = 0; k < dims; ++k) {
                    const double gap = std::abs(p[k] - q[k]);
                    squares = squares + gap * gap;
                }
                dist[at++] = std::sqrt(squares);
            }
        }
    }
    const auto between = [&](std::int64_t a, std::int64_t b) -> double& {
        return dist[sz(condensed_index(n, a, b))];
    };

    // Phase 1. members[s] is the number of observations in the cluster held by
    // slot s, 0 once the slot is retired.
    struct Merge {
        std::int64_t low_slot;
        std::int64_t high_slot;
        double height;
    };
    std::vector<Merge> merges;
    merges.reserve(sz(n - 1));
    std::vector<int> members(sz(n), 1);
    std::vector<std::int64_t> chain;
    chain.reserve(sz(n));
    const double infinity = std::numeric_limits<double>::infinity();
    while (static_cast<std::int64_t>(merges.size()) < n - 1) {
        if (chain.empty()) {
            for (std::int64_t s = 0; s < n; ++s) {
                if (members[sz(s)] > 0) {
                    chain.push_back(s);
                    break;
                }
            }
        }
        std::int64_t tail = 0;
        std::int64_t neighbour = 0;
        double gap = infinity;
        while (true) {
            tail = chain.back();
            const bool has_predecessor = chain.size() > 1;
            if (has_predecessor) {
                neighbour = chain[chain.size() - 2];
                gap = between(tail, neighbour);
            } else {
                gap = infinity;
            }
            for (std::int64_t s = 0; s < n; ++s) {
                if (members[sz(s)] == 0 || s == tail) {
                    continue;
                }
                const double candidate = between(tail, s);
                if (candidate < gap) {
                    gap = candidate;
                    neighbour = s;
                }
            }
            if (has_predecessor && neighbour == chain[chain.size() - 2]) {
                break;  // tail and its predecessor are mutual nearest neighbours
            }
            chain.push_back(neighbour);
        }
        chain.pop_back();
        chain.pop_back();

        const std::int64_t low = std::min(tail, neighbour);
        const std::int64_t high = std::max(tail, neighbour);
        const int n_low = members[sz(low)];
        const int n_high = members[sz(high)];
        merges.push_back({low, high, gap});
        members[sz(low)] = 0;
        members[sz(high)] = n_low + n_high;
        for (std::int64_t k = 0; k < n; ++k) {
            if (members[sz(k)] == 0 || k == high) {
                continue;
            }
            const double to_low = between(k, low);
            const double to_high = between(k, high);
            between(k, high) = method == LinkageMethod::Complete
                                   ? std::max(to_low, to_high)
                                   : ward_update(to_low, to_high, gap, n_low, n_high, members[sz(k)]);
        }
    }

    // Phase 2.
    std::vector<std::size_t> by_height(merges.size());
    std::iota(by_height.begin(), by_height.end(), std::size_t{0});
    std::stable_sort(by_height.begin(), by_height.end(), [&merges](std::size_t a, std::size_t b) {
        return merges[a].height < merges[b].height;
    });
    std::vector<double> linkage(sz((n - 1) * 4));
    DendrogramSets sets(n);
    for (std::size_t row = 0; row < by_height.size(); ++row) {
        const Merge& merge = merges[by_height[row]];
        const std::int64_t first = sets.root(merge.low_slot);
        const std::int64_t second = sets.root(merge.high_slot);
        linkage[row * 4] = static_cast<double>(std::min(first, second));
        linkage[row * 4 + 1] = static_cast<double>(std::max(first, second));
        linkage[row * 4 + 2] = merge.height;
        linkage[row * 4 + 3] = static_cast<double>(sets.join(first, second));
    }
    return linkage;
}

}  // namespace

std::vector<std::int64_t> hc_cut(std::int64_t n_clusters, const std::vector<std::int64_t>& children,
                                 std::int64_t n_leaves) {
    if (n_clusters > n_leaves) {
        throw ClusteringError("Cannot extract more clusters than samples: " +
                              std::to_string(n_clusters) + " clusters where given for a tree with " +
                              std::to_string(n_leaves) + " leaves.");
    }
    // Cutting a dendrogram into k clusters means undoing its last k - 1
    // merges: starting from the root, the cluster created last (the largest
    // node id; node n_leaves + i is made by merge i) is replaced by its two
    // children. The open clusters are kept in a min-heap of negated node ids so
    // that the largest id is at the root; the label of a cluster is its final
    // position in the heap array (the labelling scikit-learn produces).
    const std::size_t root = children.size() - 2;
    std::vector<std::int64_t> open{-(std::max(children[root], children[root + 1]) + 1)};
    for (std::int64_t step = 0; step < n_clusters - 1; ++step) {
        const std::int64_t node = -open[0] - n_leaves;
        const std::int64_t first = children[static_cast<std::size_t>(node * 2)];
        const std::int64_t second = children[static_cast<std::size_t>(node * 2 + 1)];
        heap_push(open, -first);
        heap_push_pop(open, -second);
    }
    // Every leaf below an open cluster gets that cluster's label.
    std::vector<std::int64_t> labels(static_cast<std::size_t>(n_leaves), 0);
    for (std::size_t i = 0; i < open.size(); ++i) {
        std::vector<std::int64_t> pending{-open[i]};
        while (!pending.empty()) {
            const std::int64_t current = pending.back();
            pending.pop_back();
            if (current < n_leaves) {
                labels[static_cast<std::size_t>(current)] = static_cast<std::int64_t>(i);
            } else {
                pending.push_back(children[static_cast<std::size_t>((current - n_leaves) * 2)]);
                pending.push_back(children[static_cast<std::size_t>((current - n_leaves) * 2 + 1)]);
            }
        }
    }
    return labels;
}

}  // namespace detail

namespace detail {

KMeansFit kmeans_fit(const Samples& X_in, std::int64_t n_clusters, std::uint32_t random_state) {
    if (X_in.samples < n_clusters) {
        throw ClusteringError("n_samples=" + std::to_string(X_in.samples) +
                              " should be >= n_clusters=" + std::to_string(n_clusters) + ".");
    }
    require_finite(X_in, "KMeans");
    // Seeding and assignment call into BLAS on one thread; see blas().

    const double tol = convergence_threshold(X_in);
    npy::RandomState rng(random_state);

    // The data are centred on the feature means before clustering (the mean is
    // the sum over the samples as in convergence_threshold) and the means are
    // added back to the returned centres.
    Samples X = X_in;
    std::vector<double> feature_means(static_cast<std::size_t>(X.features));
    for (std::int64_t j = 0; j < X.features; ++j) {
        const double mean =
            sum_over_samples(feature_column(X_in, j), X.features) / static_cast<double>(X.samples);
        feature_means[static_cast<std::size_t>(j)] = mean;
        for (std::int64_t i = 0; i < X.samples; ++i) {
            X.values[static_cast<std::size_t>(i * X.features + j)] -= mean;
        }
    }
    const std::vector<double> x_squared_norms = row_norms(X.values.data(), X.samples, X.features);

    // Ten restarts; a restart replaces the best run only when its inertia is
    // lower and it partitions the samples differently.
    LloydResult best;
    bool have_best = false;
    for (int restart = 0; restart < kRestarts; ++restart) {
        std::vector<double> seeds = kmeans_plusplus_centers(X, x_squared_norms, n_clusters, rng);
        LloydResult run = lloyd_run(X, std::move(seeds), n_clusters, tol);
        if (!have_best ||
            (run.inertia < best.inertia && !same_partition(run.labels, best.labels, n_clusters))) {
            best = std::move(run);
            have_best = true;
        }
    }
    KMeansFit fit;
    fit.labels.assign(best.labels.begin(), best.labels.end());
    fit.centers = std::move(best.centers);
    for (std::int64_t j = 0; j < n_clusters; ++j) {
        for (std::int64_t k = 0; k < X.features; ++k) {
            fit.centers[static_cast<std::size_t>(j * X.features + k)] +=
                feature_means[static_cast<std::size_t>(k)];
        }
    }
    return fit;
}

}  // namespace detail

std::vector<std::int64_t> kmeans(const Samples& X, std::int64_t n_clusters,
                                 std::uint32_t random_state) {
    return detail::kmeans_fit(X, n_clusters, random_state).labels;
}

std::vector<std::int64_t> ward(const Samples& X, std::int64_t n_clusters) {
    if (X.samples < 2) {
        throw ClusteringError("Found array with " + std::to_string(X.samples) +
                              " sample(s) (shape=(" + std::to_string(X.samples) + ", " +
                              std::to_string(X.features) +
                              ")) while a minimum of 2 is required by AgglomerativeClustering.");
    }
    require_finite(X, "AgglomerativeClustering");
    const std::vector<double> linkage = detail::ward_linkage(X);
    std::vector<std::int64_t> children(static_cast<std::size_t>((X.samples - 1) * 2));
    for (std::int64_t i = 0; i < X.samples - 1; ++i) {
        children[static_cast<std::size_t>(i * 2)] =
            static_cast<std::int64_t>(linkage[static_cast<std::size_t>(i * 4)]);
        children[static_cast<std::size_t>(i * 2 + 1)] =
            static_cast<std::int64_t>(linkage[static_cast<std::size_t>(i * 4 + 1)]);
    }
    return detail::hc_cut(n_clusters, children, X.samples);
}

}  // namespace hicx::cluster
