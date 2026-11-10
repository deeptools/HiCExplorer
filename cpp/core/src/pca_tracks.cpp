#include "hicx/pca_tracks.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "hicx/numpy_compat.hpp"

namespace hicx::pca {

namespace {

double sign(double value) {
    if (std::isnan(value)) {
        return value;
    }
    return value > 0.0 ? 1.0 : (value < 0.0 ? -1.0 : 0.0);
}

bool all_equal_first(const std::vector<double>& v) {
    for (const double value : v) {
        if (!(value == v[0])) {
            return false;
        }
    }
    return true;
}

double mean(const std::vector<double>& v) {
    return npy::pairwise_sum(v) / static_cast<double>(v.size());
}

// xmax * vector_norm(xm / xmax): np.linalg.norm along an axis is the square
// root of the pairwise sum of the squares.
double scaled_norm(const std::vector<double>& centred) {
    double largest = 0.0;
    bool nan = false;
    for (const double value : centred) {
        if (std::isnan(value)) {
            nan = true;
        }
        largest = std::max(largest, std::fabs(value));
    }
    if (nan) {
        largest = std::numeric_limits<double>::quiet_NaN();
    }
    std::vector<double> squares(centred.size());
    for (std::size_t i = 0; i < centred.size(); ++i) {
        const double scaled = centred[i] / largest;
        squares[i] = scaled * scaled;
    }
    return largest * std::sqrt(npy::pairwise_sum(squares));
}

}  // namespace

Pearson pearsonr(const std::vector<double>& x, const std::vector<double>& y) {
    const std::size_t n = x.size();
    if (n != y.size()) {
        throw std::invalid_argument("`x` and `y` must have the same length along `axis`.");
    }
    if (n < 2) {
        throw std::invalid_argument("`x` and `y` must have length at least 2.");
    }
    Pearson result;
    result.constant = all_equal_first(x) || all_equal_first(y);
    if (n == 2) {
        result.statistic = sign(x[1] - x[0]) * sign(y[1] - y[0]);
        return result;
    }
    const double x_mean = mean(x);
    const double y_mean = mean(y);
    std::vector<double> xm(n);
    std::vector<double> ym(n);
    for (std::size_t i = 0; i < n; ++i) {
        xm[i] = x[i] - x_mean;
        ym[i] = y[i] - y_mean;
    }
    const double norm_x = scaled_norm(xm);
    const double norm_y = scaled_norm(ym);
    const double threshold = std::pow(std::numeric_limits<double>::epsilon(), 0.75);
    const bool near_x = norm_x < threshold * std::fabs(x_mean);
    const bool near_y = norm_y < threshold * std::fabs(y_mean);
    result.nearly_constant = (near_x || near_y) && !result.constant;
    std::vector<double> products(n);
    for (std::size_t i = 0; i < n; ++i) {
        products[i] = xm[i] / norm_x * ym[i] / norm_y;
    }
    double r = npy::pairwise_sum(products);
    if (!std::isnan(r)) {
        r = std::min(std::max(r, -1.0), 1.0);
    }
    result.statistic = result.constant ? std::numeric_limits<double>::quiet_NaN() : r;
    return result;
}

bool histone_flip(double positive_mean, double negative_mean, bool active) {
    if (negative_mean == 0.0 || positive_mean == 0.0) {
        return false;
    }
    return active ? positive_mean < negative_mean : positive_mean > negative_mean;
}

}  // namespace hicx::pca
