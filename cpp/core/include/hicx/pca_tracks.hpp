// The sign decisions of hicPCA's --extraTrack (correlateEigenvectorWithGeneTrack
// and correlateEigenvectorWithHistonMarkTrack of hicexplorer/hicPCA.py), for
// the Python API (hicexplorer._hicx.hicPCA).

#ifndef HICX_PCA_TRACKS_HPP
#define HICX_PCA_TRACKS_HPP

#include <vector>

namespace hicx::pca {

// scipy.stats.pearsonr(x, y).statistic with scipy 1.14's float64 arithmetic:
// for two values the product of the signs of the differences, otherwise the
// sum of (x - mean x) / |x - mean x| * (y - mean y) / |y - mean y| with numpy's
// pairwise sums, the norms scaled by the largest deviation, clipped to
// [-1, 1]; NaN when either input is constant.
struct Pearson {
    double statistic = 0.0;
    // scipy's ConstantInputWarning and NearConstantInputWarning
    bool constant = false;
    bool nearly_constant = false;
};

// Throws std::invalid_argument with scipy's message for inputs of different
// lengths or of fewer than two values.
[[nodiscard]] Pearson pearsonr(const std::vector<double>& x, const std::vector<double>& y);

// Whether correlateEigenvectorWithHistonMarkTrack flips an eigenvector: for an
// active mark when the mean signal of the positive bins is below that of the
// negative ones, for an inactive mark when it is above, and never when either
// mean is 0.
[[nodiscard]] bool histone_flip(double positive_mean, double negative_mean, bool active);

}  // namespace hicx::pca

#endif  // HICX_PCA_TRACKS_HPP
