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

constexpr std::int64_t kChunkSize = 256;  // sklearn.cluster._k_means_common.CHUNK_SIZE
constexpr int kMaxIter = 300;
constexpr int kNInit = 10;
constexpr double kTol = 1e-4;

int as_int(std::int64_t value) { return static_cast<int>(value); }

// _euclidean_dense_dense: four-way unrolled, squared.
double euclidean_squared(const double* a, const double* b, std::int64_t n_features) {
    const std::int64_t n = n_features / 4;
    const std::int64_t rem = n_features % 4;
    double result = 0.0;
    for (std::int64_t i = 0; i < n; ++i) {
        result += ((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                   (a[2] - b[2]) * (a[2] - b[2]) + (a[3] - b[3]) * (a[3] - b[3]));
        a += 4;
        b += 4;
    }
    for (std::int64_t i = 0; i < rem; ++i) {
        result += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return result;
}

std::vector<double> row_norms(const double* values, std::int64_t rows, std::int64_t features) {
    std::vector<double> norms(static_cast<std::size_t>(rows));
    for (std::int64_t i = 0; i < rows; ++i) {
        norms[static_cast<std::size_t>(i)] =
            detail::einsum_row_norm_squared(values + i * features, features);
    }
    return norms;
}

// _tolerance(X, tol) = np.mean(np.var(X, axis=0)) * tol. The column
// reductions of np.var run sequentially over the samples; np.mean over the
// variances is numpy's pairwise sum.
double tolerance(const Samples& X) {
    const std::int64_t n = X.samples;
    const std::int64_t d = X.features;
    std::vector<double> variances(static_cast<std::size_t>(d));
    for (std::int64_t j = 0; j < d; ++j) {
        double mean = 0.0;
        for (std::int64_t i = 0; i < n; ++i) {
            mean += X.row(i)[j];
        }
        mean /= static_cast<double>(n);
        double squares = 0.0;
        for (std::int64_t i = 0; i < n; ++i) {
            const double deviation = X.row(i)[j] - mean;
            squares += deviation * deviation;
        }
        variances[static_cast<std::size_t>(j)] = squares / static_cast<double>(n);
    }
    const double mean_variance =
        npy::pairwise_sum(variances.data(), variances.size()) / static_cast<double>(d);
    return mean_variance * kTol;
}

// _kmeans_plusplus with unit sample weights.
std::vector<double> kmeans_plusplus(const Samples& X, const std::vector<double>& x_squared_norms,
                                    std::int64_t n_clusters, npy::RandomState& random_state) {
    const std::int64_t n = X.samples;
    const std::int64_t d = X.features;
    std::vector<double> centers(static_cast<std::size_t>(n_clusters * d));
    const std::int64_t n_local_trials =
        2 + static_cast<std::int64_t>(std::log(static_cast<double>(n_clusters)));

    // sample_weight / sample_weight.sum(): the sum of n ones is exactly n.
    std::vector<double> p(static_cast<std::size_t>(n), 1.0 / static_cast<double>(n));
    const std::int64_t center_id = random_state.choice(p);
    std::copy(X.row(center_id), X.row(center_id) + d, centers.begin());

    // _euclidean_distances(centers[0, np.newaxis], X, Y_norm_squared=...):
    // numpy's matmul of a (1, d) row and X.T is the vector @ matrix case, a
    // gemv on X in column major order.
    const std::vector<double> ones(static_cast<std::size_t>(n), 1.0);
    std::vector<double> closest(static_cast<std::size_t>(n));
    cblas_dgemv(kCblasColMajor, kCblasTrans, as_int(d), as_int(n), 1.0, X.values.data(),
                as_int(d), centers.data(), 1, 0.0, closest.data(), 1);
    {
        const double xx = detail::einsum_row_norm_squared(centers.data(), d);
        for (std::int64_t j = 0; j < n; ++j) {
            double value = -2.0 * closest[static_cast<std::size_t>(j)];
            value += xx;
            value += x_squared_norms[static_cast<std::size_t>(j)];
            closest[static_cast<std::size_t>(j)] = std::max(value, 0.0);
        }
    }
    // closest_dist_sq @ sample_weight: a (1, n) row times a vector is numpy's
    // scalar-output case, DOUBLE_dot, which starts from 0.0.
    double current_pot = 0.0 + cblas_ddot(as_int(n), closest.data(), 1, ones.data(), 1);

    std::vector<double> candidates(static_cast<std::size_t>(n_local_trials * d));
    std::vector<double> distances(static_cast<std::size_t>(n_local_trials * n));
    std::vector<double> candidate_pots(static_cast<std::size_t>(n_local_trials));
    std::vector<double> cumulative(static_cast<std::size_t>(n));
    std::vector<std::int64_t> candidate_ids(static_cast<std::size_t>(n_local_trials));

    for (std::int64_t c = 1; c < n_clusters; ++c) {
        std::vector<double> rand_vals = random_state.uniform(static_cast<std::size_t>(n_local_trials));
        for (double& value : rand_vals) {
            value *= current_pot;
        }
        // stable_cumsum(sample_weight * closest_dist_sq): np.cumsum, sequential.
        double running = 0.0;
        for (std::int64_t j = 0; j < n; ++j) {
            running += 1.0 * closest[static_cast<std::size_t>(j)];
            cumulative[static_cast<std::size_t>(j)] = running;
        }
        for (std::int64_t t = 0; t < n_local_trials; ++t) {
            // np.searchsorted(side='left'), then np.clip(..., None, n - 1).
            auto position = static_cast<std::int64_t>(
                std::lower_bound(cumulative.begin(), cumulative.end(),
                                 rand_vals[static_cast<std::size_t>(t)]) -
                cumulative.begin());
            position = std::min(position, n - 1);
            candidate_ids[static_cast<std::size_t>(t)] = position;
            std::copy(X.row(position), X.row(position) + d,
                      candidates.begin() + t * d);
        }
        // X[candidate_ids] @ X.T: matrix @ matrix, a C contiguous left operand
        // and a transposed view on the right, so gemm(RowMajor, NoTrans, Trans).
        cblas_dgemm(kCblasRowMajor, kCblasNoTrans, kCblasTrans, as_int(n_local_trials),
                    as_int(n), as_int(d), 1.0, candidates.data(), as_int(d), X.values.data(),
                    as_int(d), 0.0, distances.data(), as_int(n));
        for (std::int64_t t = 0; t < n_local_trials; ++t) {
            const double xx = detail::einsum_row_norm_squared(candidates.data() + t * d, d);
            for (std::int64_t j = 0; j < n; ++j) {
                double& value = distances[static_cast<std::size_t>(t * n + j)];
                value = -2.0 * value;
                value += xx;
                value += x_squared_norms[static_cast<std::size_t>(j)];
                value = std::max(value, 0.0);
                // np.minimum(closest_dist_sq, distance_to_candidates)
                const double current = closest[static_cast<std::size_t>(j)];
                value = current <= value ? current : value;
            }
        }
        // distance_to_candidates @ sample_weight.reshape(-1, 1): matrix @
        // column vector, a gemv in column major order.
        cblas_dgemv(kCblasColMajor, kCblasTrans, as_int(n), as_int(n_local_trials), 1.0,
                    distances.data(), as_int(n), ones.data(), 1, 0.0, candidate_pots.data(), 1);
        std::int64_t best = 0;
        for (std::int64_t t = 1; t < n_local_trials; ++t) {
            if (candidate_pots[static_cast<std::size_t>(t)] <
                candidate_pots[static_cast<std::size_t>(best)]) {
                best = t;
            }
        }
        current_pot = candidate_pots[static_cast<std::size_t>(best)];
        std::copy(distances.begin() + best * n, distances.begin() + (best + 1) * n,
                  closest.begin());
        const std::int64_t chosen = candidate_ids[static_cast<std::size_t>(best)];
        std::copy(X.row(chosen), X.row(chosen) + d, centers.begin() + c * d);
    }
    return centers;
}

// _relocate_empty_clusters_dense with unit weights.
void relocate_empty_clusters(const Samples& X, const std::vector<double>& centers_old,
                             std::vector<double>& centers_new,
                             std::vector<double>& weight_in_clusters,
                             const std::vector<std::int32_t>& labels, std::int64_t n_clusters) {
    std::vector<std::int64_t> empty;
    for (std::int64_t j = 0; j < n_clusters; ++j) {
        if (weight_in_clusters[static_cast<std::size_t>(j)] == 0.0) {
            empty.push_back(j);
        }
    }
    if (empty.empty()) {
        return;
    }
    const std::int64_t n = X.samples;
    const std::int64_t d = X.features;
    const auto n_empty = static_cast<std::int64_t>(empty.size());
    // ((X - centers_old[labels]) ** 2).sum(axis=1): a pairwise sum per row.
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
    // np.argpartition(distances, -n_empty)[:-n_empty-1:-1]
    const std::vector<std::int64_t> partition = npy::argpartition(distances, -n_empty);
    for (std::int64_t idx = 0; idx < n_empty; ++idx) {
        const std::int64_t new_cluster = empty[static_cast<std::size_t>(idx)];
        const std::int64_t far = partition[static_cast<std::size_t>(n - 1 - idx)];
        const double weight = 1.0;
        const std::int64_t old_cluster = labels[static_cast<std::size_t>(far)];
        for (std::int64_t k = 0; k < d; ++k) {
            centers_new[static_cast<std::size_t>(old_cluster * d + k)] -= X.row(far)[k] * weight;
            centers_new[static_cast<std::size_t>(new_cluster * d + k)] = X.row(far)[k] * weight;
        }
        weight_in_clusters[static_cast<std::size_t>(new_cluster)] = weight;
        weight_in_clusters[static_cast<std::size_t>(old_cluster)] -= weight;
    }
}

// lloyd_iter_chunked_dense with unit weights.
void lloyd_iteration(const Samples& X, const std::vector<double>& centers_old,
                     std::vector<double>& centers_new, std::vector<double>& weight_in_clusters,
                     std::vector<std::int32_t>& labels, std::vector<double>& center_shift,
                     std::int64_t n_clusters, bool update_centers) {
    const std::int64_t n = X.samples;
    const std::int64_t d = X.features;
    const std::int64_t samples_per_chunk = n > kChunkSize ? kChunkSize : n;
    const std::int64_t remainder = n % samples_per_chunk;
    const std::int64_t chunks = n / samples_per_chunk + (remainder != 0 ? 1 : 0);
    const std::vector<double> centers_squared_norms = row_norms(centers_old.data(), n_clusters, d);

    if (update_centers) {
        std::fill(centers_new.begin(), centers_new.end(), 0.0);
        std::fill(weight_in_clusters.begin(), weight_in_clusters.end(), 0.0);
    }
    std::vector<double> chunk_centers(static_cast<std::size_t>(n_clusters * d));
    std::vector<double> chunk_weights(static_cast<std::size_t>(n_clusters));
    std::vector<double> pairwise(static_cast<std::size_t>(samples_per_chunk * n_clusters));

    for (std::int64_t chunk = 0; chunk < chunks; ++chunk) {
        const std::int64_t start = chunk * samples_per_chunk;
        const std::int64_t end = (chunk == chunks - 1 && remainder > 0)
                                     ? start + remainder
                                     : start + samples_per_chunk;
        const std::int64_t m = end - start;
        std::fill(chunk_centers.begin(), chunk_centers.end(), 0.0);
        std::fill(chunk_weights.begin(), chunk_weights.end(), 0.0);

        for (std::int64_t i = 0; i < m; ++i) {
            for (std::int64_t j = 0; j < n_clusters; ++j) {
                pairwise[static_cast<std::size_t>(i * n_clusters + j)] =
                    centers_squared_norms[static_cast<std::size_t>(j)];
            }
        }
        // _gemm(RowMajor, NoTrans, Trans, m, k, d, -2.0, X, d, centers, d,
        // 1.0, pairwise, k), which scikit-learn turns into this Fortran call.
        const char trans_b = 't';
        const char trans_a = 'n';
        const int k_int = as_int(n_clusters);
        const int m_int = as_int(m);
        const int d_int = as_int(d);
        const double alpha = -2.0;
        const double beta = 1.0;
        dgemm_(&trans_b, &trans_a, &k_int, &m_int, &d_int, &alpha, centers_old.data(), &d_int,
               X.row(start), &d_int, &beta, pairwise.data(), &k_int);

        for (std::int64_t i = 0; i < m; ++i) {
            double min_sq_dist = pairwise[static_cast<std::size_t>(i * n_clusters)];
            std::int32_t label = 0;
            for (std::int64_t j = 1; j < n_clusters; ++j) {
                const double sq_dist = pairwise[static_cast<std::size_t>(i * n_clusters + j)];
                if (sq_dist < min_sq_dist) {
                    min_sq_dist = sq_dist;
                    label = static_cast<std::int32_t>(j);
                }
            }
            labels[static_cast<std::size_t>(start + i)] = label;
            if (update_centers) {
                chunk_weights[static_cast<std::size_t>(label)] += 1.0;
                const double* sample = X.row(start + i);
                for (std::int64_t k = 0; k < d; ++k) {
                    chunk_centers[static_cast<std::size_t>(label * d + k)] += sample[k] * 1.0;
                }
            }
        }
        if (update_centers) {
            for (std::int64_t j = 0; j < n_clusters; ++j) {
                weight_in_clusters[static_cast<std::size_t>(j)] +=
                    chunk_weights[static_cast<std::size_t>(j)];
                for (std::int64_t k = 0; k < d; ++k) {
                    centers_new[static_cast<std::size_t>(j * d + k)] +=
                        chunk_centers[static_cast<std::size_t>(j * d + k)];
                }
            }
        }
    }
    if (!update_centers) {
        return;
    }
    relocate_empty_clusters(X, centers_old, centers_new, weight_in_clusters, labels, n_clusters);
    // _average_centers
    for (std::int64_t j = 0; j < n_clusters; ++j) {
        const double weight = weight_in_clusters[static_cast<std::size_t>(j)];
        if (weight > 0) {
            const double alpha_j = 1.0 / weight;
            for (std::int64_t k = 0; k < d; ++k) {
                centers_new[static_cast<std::size_t>(j * d + k)] *= alpha_j;
            }
        }
    }
    // _center_shift
    for (std::int64_t j = 0; j < n_clusters; ++j) {
        center_shift[static_cast<std::size_t>(j)] = std::sqrt(euclidean_squared(
            centers_new.data() + j * d, centers_old.data() + j * d, d));
    }
}

struct LloydResult {
    std::vector<std::int32_t> labels;
    double inertia = 0.0;
};

// _kmeans_single_lloyd with unit weights.
LloydResult kmeans_single_lloyd(const Samples& X, std::vector<double> centers,
                                std::int64_t n_clusters, double tol) {
    const std::int64_t n = X.samples;
    const std::int64_t d = X.features;
    std::vector<double> centers_new(centers.size(), 0.0);
    std::vector<std::int32_t> labels(static_cast<std::size_t>(n), -1);
    std::vector<std::int32_t> labels_old = labels;
    std::vector<double> weight_in_clusters(static_cast<std::size_t>(n_clusters), 0.0);
    std::vector<double> center_shift(static_cast<std::size_t>(n_clusters), 0.0);

    bool strict_convergence = false;
    for (int iteration = 0; iteration < kMaxIter; ++iteration) {
        lloyd_iteration(X, centers, centers_new, weight_in_clusters, labels, center_shift,
                        n_clusters, true);
        std::swap(centers, centers_new);
        if (labels == labels_old) {
            strict_convergence = true;
            break;
        }
        std::vector<double> squared(center_shift.size());
        for (std::size_t j = 0; j < center_shift.size(); ++j) {
            squared[j] = center_shift[j] * center_shift[j];
        }
        const double center_shift_tot = npy::pairwise_sum(squared.data(), squared.size());
        if (center_shift_tot <= tol) {
            break;
        }
        labels_old = labels;
    }
    if (!strict_convergence) {
        std::vector<double> unused;
        lloyd_iteration(X, centers, unused, weight_in_clusters, labels, center_shift, n_clusters,
                        false);
    }
    LloydResult result;
    for (std::int64_t i = 0; i < n; ++i) {
        result.inertia += euclidean_squared(
            X.row(i), centers.data() + labels[static_cast<std::size_t>(i)] * d, d) * 1.0;
    }
    result.labels = std::move(labels);
    return result;
}

bool is_same_clustering(const std::vector<std::int32_t>& labels1,
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

// Python's heapq on a list of ints.
void heap_siftdown(std::vector<std::int64_t>& heap, std::size_t startpos, std::size_t pos) {
    const std::int64_t newitem = heap[pos];
    while (pos > startpos) {
        const std::size_t parentpos = (pos - 1) >> 1;
        const std::int64_t parent = heap[parentpos];
        if (newitem < parent) {
            heap[pos] = parent;
            pos = parentpos;
            continue;
        }
        break;
    }
    heap[pos] = newitem;
}

void heap_siftup(std::vector<std::int64_t>& heap, std::size_t pos) {
    const std::size_t endpos = heap.size();
    const std::size_t startpos = pos;
    const std::int64_t newitem = heap[pos];
    std::size_t childpos = 2 * pos + 1;
    while (childpos < endpos) {
        const std::size_t rightpos = childpos + 1;
        if (rightpos < endpos && !(heap[childpos] < heap[rightpos])) {
            childpos = rightpos;
        }
        heap[pos] = heap[childpos];
        pos = childpos;
        childpos = 2 * pos + 1;
    }
    heap[pos] = newitem;
    heap_siftdown(heap, startpos, pos);
}

void heappush(std::vector<std::int64_t>& heap, std::int64_t item) {
    heap.push_back(item);
    heap_siftdown(heap, 0, heap.size() - 1);
}

void heappushpop(std::vector<std::int64_t>& heap, std::int64_t item) {
    if (!heap.empty() && heap[0] < item) {
        std::swap(item, heap[0]);
        heap_siftup(heap, 0);
    }
}

// scipy's LinkageUnionFind.
struct LinkageUnionFind {
    explicit LinkageUnionFind(std::int64_t n)
        : parent(static_cast<std::size_t>(2 * n - 1)), size(static_cast<std::size_t>(2 * n - 1), 1),
          next_label(n) {
        std::iota(parent.begin(), parent.end(), std::int64_t{0});
    }
    std::int64_t merge(std::int64_t x, std::int64_t y) {
        parent[static_cast<std::size_t>(x)] = next_label;
        parent[static_cast<std::size_t>(y)] = next_label;
        const std::int64_t merged =
            size[static_cast<std::size_t>(x)] + size[static_cast<std::size_t>(y)];
        size[static_cast<std::size_t>(next_label)] = merged;
        next_label += 1;
        return merged;
    }
    std::int64_t find(std::int64_t x) {
        std::int64_t p = x;
        while (parent[static_cast<std::size_t>(x)] != x) {
            x = parent[static_cast<std::size_t>(x)];
        }
        while (parent[static_cast<std::size_t>(p)] != x) {
            const std::int64_t next = parent[static_cast<std::size_t>(p)];
            parent[static_cast<std::size_t>(p)] = x;
            p = next;
        }
        return x;
    }
    std::vector<std::int64_t> parent;
    std::vector<std::int64_t> size;
    std::int64_t next_label;
};

std::int64_t condensed_index(std::int64_t n, std::int64_t i, std::int64_t j) {
    if (i < j) {
        return n * i - (i * (i + 1) / 2) + (j - i - 1);
    }
    return n * j - (j * (j + 1) / 2) + (i - j - 1);
}

// _ward in _hierarchy_distance_update.pxi, evaluated left to right as C does.
double ward_update(double d_xi, double d_yi, double d_xy, int size_x, int size_y, int size_i) {
    const double t = 1.0 / (size_x + size_y + size_i);
    return std::sqrt((size_i + size_x) * t * d_xi * d_xi + (size_i + size_y) * t * d_yi * d_yi -
                     size_i * t * d_xy * d_xy);
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
    // DOUBLE_sum_of_products_contig_contig_outstride0_two at an SSE baseline:
    // two lanes, blocks of four vectors combined from the last to the first,
    // the tail loaded with zero padding, the lanes added at the end.
    double lane0 = 0.0;
    double lane1 = 0.0;
    std::int64_t count = features;
    const double* x = row;
    for (; count >= 8; count -= 8, x += 8) {
        lane0 = x[6] * x[6] + lane0;
        lane1 = x[7] * x[7] + lane1;
        lane0 = x[4] * x[4] + lane0;
        lane1 = x[5] * x[5] + lane1;
        lane0 = x[2] * x[2] + lane0;
        lane1 = x[3] * x[3] + lane1;
        lane0 = x[0] * x[0] + lane0;
        lane1 = x[1] * x[1] + lane1;
    }
    for (; count > 0; count -= 2, x += 2) {
        lane0 = x[0] * x[0] + lane0;
        const double second = count >= 2 ? x[1] : 0.0;
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

std::vector<double> nn_chain_linkage(const Samples& X, LinkageMethod method) {
    const std::int64_t n = X.samples;
    const std::int64_t d = X.features;
    // pdist(X, 'euclidean'): sequential sum of squared differences, sqrt.
    std::vector<double> D(static_cast<std::size_t>(n * (n - 1) / 2));
    {
        std::size_t position = 0;
        for (std::int64_t i = 0; i < n; ++i) {
            for (std::int64_t j = i + 1; j < n; ++j) {
                double sum = 0.0;
                const double* a = X.row(i);
                const double* b = X.row(j);
                for (std::int64_t k = 0; k < d; ++k) {
                    const double difference = std::abs(a[k] - b[k]);
                    sum = sum + difference * difference;
                }
                D[position++] = std::sqrt(sum);
            }
        }
    }

    // nn_chain(dists, n, method=ward)
    std::vector<double> Z(static_cast<std::size_t>((n - 1) * 4));
    std::vector<int> size(static_cast<std::size_t>(n), 1);
    std::vector<std::int64_t> cluster_chain(static_cast<std::size_t>(n));
    std::int64_t chain_length = 0;
    std::int64_t x = 0;
    std::int64_t y = 0;
    for (std::int64_t k = 0; k < n - 1; ++k) {
        if (chain_length == 0) {
            chain_length = 1;
            for (std::int64_t i = 0; i < n; ++i) {
                if (size[static_cast<std::size_t>(i)] > 0) {
                    cluster_chain[0] = i;
                    break;
                }
            }
        }
        double current_min = 0.0;
        while (true) {
            x = cluster_chain[static_cast<std::size_t>(chain_length - 1)];
            if (chain_length > 1) {
                y = cluster_chain[static_cast<std::size_t>(chain_length - 2)];
                current_min = D[static_cast<std::size_t>(condensed_index(n, x, y))];
            } else {
                current_min = std::numeric_limits<double>::infinity();
            }
            for (std::int64_t i = 0; i < n; ++i) {
                if (size[static_cast<std::size_t>(i)] == 0 || x == i) {
                    continue;
                }
                const double dist = D[static_cast<std::size_t>(condensed_index(n, x, i))];
                if (dist < current_min) {
                    current_min = dist;
                    y = i;
                }
            }
            if (chain_length > 1 && y == cluster_chain[static_cast<std::size_t>(chain_length - 2)]) {
                break;
            }
            cluster_chain[static_cast<std::size_t>(chain_length)] = y;
            chain_length += 1;
        }
        chain_length -= 2;
        if (x > y) {
            std::swap(x, y);
        }
        const int nx = size[static_cast<std::size_t>(x)];
        const int ny = size[static_cast<std::size_t>(y)];
        Z[static_cast<std::size_t>(k * 4)] = static_cast<double>(x);
        Z[static_cast<std::size_t>(k * 4 + 1)] = static_cast<double>(y);
        Z[static_cast<std::size_t>(k * 4 + 2)] = current_min;
        Z[static_cast<std::size_t>(k * 4 + 3)] = static_cast<double>(nx + ny);
        size[static_cast<std::size_t>(x)] = 0;
        size[static_cast<std::size_t>(y)] = nx + ny;
        for (std::int64_t i = 0; i < n; ++i) {
            const int ni = size[static_cast<std::size_t>(i)];
            if (ni == 0 || i == y) {
                continue;
            }
            const double d_xi = D[static_cast<std::size_t>(condensed_index(n, i, x))];
            const double d_yi = D[static_cast<std::size_t>(condensed_index(n, i, y))];
            // _complete in _hierarchy_distance_update.pxi: max(d_xi, d_yi).
            D[static_cast<std::size_t>(condensed_index(n, i, y))] =
                method == LinkageMethod::Complete ? std::max(d_xi, d_yi)
                                                  : ward_update(d_xi, d_yi, current_min, nx, ny, ni);
        }
    }

    // np.argsort(Z[:, 2], kind='mergesort'), a stable sort.
    std::vector<std::int64_t> order(static_cast<std::size_t>(n - 1));
    std::iota(order.begin(), order.end(), std::int64_t{0});
    std::stable_sort(order.begin(), order.end(), [&Z](std::int64_t a, std::int64_t b) {
        return Z[static_cast<std::size_t>(a * 4 + 2)] < Z[static_cast<std::size_t>(b * 4 + 2)];
    });
    std::vector<double> sorted(Z.size());
    for (std::size_t row = 0; row < order.size(); ++row) {
        std::copy(Z.begin() + order[row] * 4, Z.begin() + order[row] * 4 + 4,
                  sorted.begin() + static_cast<std::ptrdiff_t>(row * 4));
    }

    // label(Z, n)
    LinkageUnionFind union_find(n);
    for (std::int64_t i = 0; i < n - 1; ++i) {
        const auto a = static_cast<std::int64_t>(sorted[static_cast<std::size_t>(i * 4)]);
        const auto b = static_cast<std::int64_t>(sorted[static_cast<std::size_t>(i * 4 + 1)]);
        const std::int64_t root_a = union_find.find(a);
        const std::int64_t root_b = union_find.find(b);
        if (root_a < root_b) {
            sorted[static_cast<std::size_t>(i * 4)] = static_cast<double>(root_a);
            sorted[static_cast<std::size_t>(i * 4 + 1)] = static_cast<double>(root_b);
        } else {
            sorted[static_cast<std::size_t>(i * 4)] = static_cast<double>(root_b);
            sorted[static_cast<std::size_t>(i * 4 + 1)] = static_cast<double>(root_a);
        }
        sorted[static_cast<std::size_t>(i * 4 + 3)] =
            static_cast<double>(union_find.merge(root_a, root_b));
    }
    return sorted;
}

}  // namespace

std::vector<std::int64_t> hc_cut(std::int64_t n_clusters, const std::vector<std::int64_t>& children,
                                 std::int64_t n_leaves) {
    if (n_clusters > n_leaves) {
        throw ClusteringError("Cannot extract more clusters than samples: " +
                              std::to_string(n_clusters) + " clusters where given for a tree with " +
                              std::to_string(n_leaves) + " leaves.");
    }
    const std::size_t root = children.size() - 2;
    std::vector<std::int64_t> nodes{-(std::max(children[root], children[root + 1]) + 1)};
    for (std::int64_t step = 0; step < n_clusters - 1; ++step) {
        const std::int64_t node = -nodes[0] - n_leaves;
        const std::int64_t first = children[static_cast<std::size_t>(node * 2)];
        const std::int64_t second = children[static_cast<std::size_t>(node * 2 + 1)];
        heappush(nodes, -first);
        heappushpop(nodes, -second);
    }
    std::vector<std::int64_t> labels(static_cast<std::size_t>(n_leaves), 0);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        // _hc_get_descendent(-node, children, n_leaves)
        std::vector<std::int64_t> pending{-nodes[i]};
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

std::vector<std::int64_t> kmeans(const Samples& X_in, std::int64_t n_clusters,
                                 std::uint32_t random_state) {
    if (X_in.samples < n_clusters) {
        throw ClusteringError("n_samples=" + std::to_string(X_in.samples) +
                              " should be >= n_clusters=" + std::to_string(n_clusters) + ".");
    }
    require_finite(X_in, "KMeans");
    // The E-step and the seeding call into BLAS on one thread; see blas().

    const double tol = tolerance(X_in);
    npy::RandomState rng(random_state);

    // X -= X.mean(axis=0)
    Samples X = X_in;
    for (std::int64_t j = 0; j < X.features; ++j) {
        double mean = 0.0;
        for (std::int64_t i = 0; i < X.samples; ++i) {
            mean += X_in.row(i)[j];
        }
        mean /= static_cast<double>(X.samples);
        for (std::int64_t i = 0; i < X.samples; ++i) {
            X.values[static_cast<std::size_t>(i * X.features + j)] -= mean;
        }
    }
    const std::vector<double> x_squared_norms = row_norms(X.values.data(), X.samples, X.features);

    std::vector<std::int32_t> best_labels;
    double best_inertia = 0.0;
    bool have_best = false;
    for (int init = 0; init < kNInit; ++init) {
        std::vector<double> centers = kmeans_plusplus(X, x_squared_norms, n_clusters, rng);
        LloydResult run = kmeans_single_lloyd(X, std::move(centers), n_clusters, tol);
        if (!have_best || (run.inertia < best_inertia &&
                           !is_same_clustering(run.labels, best_labels, n_clusters))) {
            best_labels = std::move(run.labels);
            best_inertia = run.inertia;
            have_best = true;
        }
    }
    return {best_labels.begin(), best_labels.end()};
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
