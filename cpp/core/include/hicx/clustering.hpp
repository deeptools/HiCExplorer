// Clustering of dense sample vectors, reproducing scikit-learn's partitions.
//
// Two methods are implemented, because they are the two a ported tool
// actually executes (hicAggregateContacts, --kmeans and --hclust):
//
//  * kmeans: sklearn.cluster.KMeans(n_clusters=k, random_state=seed).fit(X)
//    .labels_ with the scikit-learn 1.3.2 defaults, that is k-means++ seeding
//    with 2 + int(log(k)) local trials, n_init=10 (the "warn" default resolves
//    to 10 in 1.3), max_iter=300, tol=1e-4 scaled by the mean feature
//    variance, and the Lloyd algorithm on mean-centred data.
//  * ward: sklearn.cluster.AgglomerativeClustering(n_clusters=k).fit(X)
//    .labels_, which for Euclidean ward linkage without connectivity builds
//    the full tree with scipy.cluster.hierarchy.ward (pdist, then the
//    nearest-neighbour chain algorithm with the Lance-Williams ward update)
//    and cuts it with sklearn's _hc_cut.
//
// Spectral clustering is not implemented: hicAggregateContacts parses
// --spectral but never reads it (pinned by
// test_spectral_option_is_ignored), so no code path reaches it.
//
// "The same result" for a clustering is the partition, and because the output
// file of a cluster is named after its label the labels themselves are
// reproduced, not only the partition up to relabelling. That requires the same
// arithmetic where a comparison can flip, so:
//
//  * the random draws are numpy's RandomState (numpy_random.hpp);
//  * every dot product that numpy or scikit-learn hands to BLAS is handed to
//    the same OpenBLAS entry point with the same arguments here, loaded with
//    dlopen on first use (see clustering.cpp): cblas_dgemv
//    and cblas_dgemm for the matmul calls of the seeding (numpy's matmul
//    selects gemv for a row or column operand and gemm otherwise), cblas_ddot
//    for the potential, and Fortran dgemm with the transposition scikit-learn's
//    _gemm applies for the RowMajor E-step. The library is the one numpy is
//    linked against, run on one thread;
//  * row norms follow np.einsum('ij,ij->i'), whose inner loop at numpy's
//    baseline (SSE3, two float64 lanes) accumulates even and odd elements in
//    separate lanes, in blocks of eight taken from the back, and adds the two
//    lanes at the end;
//  * the Euclidean kernels are scikit-learn's four-way unrolled loop and
//    scipy's sequential pdist loop;
//  * the E-step works on chunks of 256 samples whose centre sums are
//    accumulated per chunk and combined in chunk order. scikit-learn combines
//    them in the order its OpenMP threads finish, which gives the same bits
//    for up to two chunks (512 samples) because a two-term sum commutes, and
//    may differ in the last bits beyond that; this port is deterministic.
//  * the inertia that picks the best of the ten initialisations is summed
//    sequentially; scikit-learn sums it with an OpenMP reduction. It only
//    decides between two different partitions whose inertias agree to the
//    last bits, which does not happen on real data.
//
// Measured agreement is recorded in the hicAggregateContacts port.

#ifndef HICX_CLUSTERING_HPP
#define HICX_CLUSTERING_HPP

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace hicx::cluster {

// A dense float64 matrix of samples, row major: sample i is
// values[i * features .. (i + 1) * features).
struct Samples {
    std::int64_t samples = 0;
    std::int64_t features = 0;
    std::vector<double> values;

    [[nodiscard]] const double* row(std::int64_t i) const {
        return values.data() + i * features;
    }
};

// The ValueError scikit-learn raises for input it refuses, with its message.
class ClusteringError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

// KMeans(n_clusters, random_state).fit(X).labels_
[[nodiscard]] std::vector<std::int64_t> kmeans(const Samples& X, std::int64_t n_clusters,
                                               std::uint32_t random_state = 0);

// AgglomerativeClustering(n_clusters).fit(X).labels_ (ward, Euclidean)
[[nodiscard]] std::vector<std::int64_t> ward(const Samples& X, std::int64_t n_clusters);

namespace detail {

// np.einsum('ij,ij->i', X, X) for one row, bit for bit.
[[nodiscard]] double einsum_row_norm_squared(const double* row, std::int64_t features);

// scipy.cluster.hierarchy.ward(X): the (n - 1) by 4 linkage matrix, row major.
[[nodiscard]] std::vector<double> ward_linkage(const Samples& X);

// scipy.cluster.hierarchy.linkage(X, method='complete') for hicCorrelate: the
// same pdist, nn_chain, stable sort and label() as ward_linkage, with the
// complete linkage update max(d_xi, d_yi). Throws ClusteringError for fewer
// than two observations or a non-finite value, as scipy raises.
[[nodiscard]] std::vector<double> complete_linkage(const Samples& X);

// scipy.cluster.hierarchy.dendrogram(Z)['leaves'] with its defaults
// (count_sort and distance_sort off): every merge lists its first child's
// leaves before its second child's.
[[nodiscard]] std::vector<std::int64_t> dendrogram_leaves(const std::vector<double>& Z,
                                                          std::int64_t n);

// sklearn.cluster._agglomerative._hc_cut(n_clusters, children, n_leaves), with
// children the first two columns of the linkage matrix.
[[nodiscard]] std::vector<std::int64_t> hc_cut(std::int64_t n_clusters,
                                               const std::vector<std::int64_t>& children,
                                               std::int64_t n_leaves);

}  // namespace detail

}  // namespace hicx::cluster

#endif  // HICX_CLUSTERING_HPP
