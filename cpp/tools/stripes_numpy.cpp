// See stripes_numpy.hpp. Each routine follows the library source named in the
// comment above it; the unit tests compare them with the libraries
// themselves, on fixtures produced by cpp/scripts/stripenn_oracle.

#include "stripes_numpy.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace hicx::stripes::np {

// ---------------------------------------------------------------------------
// numpy

namespace {

double pairwise(const double* a, std::size_t n) {
    if (n < 8) {
        double res = -0.0;
        for (std::size_t i = 0; i < n; ++i) {
            res += a[i];
        }
        return res;
    }
    if (n <= 128) {
        double r[8];
        for (int j = 0; j < 8; ++j) {
            r[j] = a[j];
        }
        std::size_t i = 8;
        for (; i < n - (n % 8); i += 8) {
            for (int j = 0; j < 8; ++j) {
                r[j] += a[i + j];
            }
        }
        double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; ++i) {
            res += a[i];
        }
        return res;
    }
    std::size_t n2 = n / 2;
    n2 -= n2 % 8;
    return pairwise(a, n2) + pairwise(a + n2, n - n2);
}

}  // namespace

double pairwise_sum(const double* values, std::size_t count) { return pairwise(values, count); }

double sum_compact(const double* compact, std::size_t rows, std::size_t cols,
                   bool contiguous_in_source) {
    if (rows == 0 || cols == 0) {
        return 0.0;
    }
    if (rows == 1 || cols == 1 || contiguous_in_source) {
        return pairwise(compact, rows * cols);
    }
    constexpr std::size_t kBuffer = 8192;
    const std::size_t total = rows * cols;
    double acc = 0.0;
    bool first = true;
    for (std::size_t done = 0; done < total; done += kBuffer) {
        const std::size_t take = std::min(kBuffer, total - done);
        const double part = pairwise(compact + done, take);
        if (first) {
            acc = part;
            first = false;
        } else {
            acc += part;
        }
    }
    return acc;
}

double mean_compact(const double* compact, std::size_t rows, std::size_t cols,
                    bool contiguous_in_source) {
    if (rows == 0 || cols == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return sum_compact(compact, rows, cols, contiguous_in_source) / static_cast<double>(rows * cols);
}

double sum_block(const double* base, std::size_t rows, std::size_t cols, std::size_t row_stride) {
    if (rows == 0 || cols == 0) {
        return 0.0;
    }
    if (rows == 1 || row_stride == cols) {
        return pairwise(base, rows * cols);
    }
    std::vector<double> compact(rows * cols);
    for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t k = 0; k < cols; ++k) {
            compact[r * cols + k] = base[r * row_stride + k];
        }
    }
    return sum_compact(compact.data(), rows, cols, false);
}

double mean_block(const double* base, std::size_t rows, std::size_t cols, std::size_t row_stride) {
    if (rows == 0 || cols == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return sum_block(base, rows, cols, row_stride) / static_cast<double>(rows * cols);
}

double mean_of_row(const double* row, std::size_t cols) {
    if (cols == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return pairwise(row, cols) / static_cast<double>(cols);
}

double mean_of(const std::vector<double>& values) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return pairwise(values.data(), values.size()) / static_cast<double>(values.size());
}

double round_half_even(double value) { return std::nearbyint(value); }

double median(std::vector<double> values) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    for (const double v : values) {
        if (std::isnan(v)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    std::sort(values.begin(), values.end());
    const std::size_t n = values.size();
    if (n % 2 == 1) {
        return values[n / 2];
    }
    const double pair[2] = {values[n / 2 - 1], values[n / 2]};
    return pairwise(pair, 2) / 2.0;
}

double nanmedian(std::vector<double> values) {
    values.erase(std::remove_if(values.begin(), values.end(), [](double v) { return std::isnan(v); }),
                 values.end());
    return median(std::move(values));
}

QuantileRanks quantile_ranks(std::int64_t n, double q) {
    QuantileRanks ranks;
    const double virtual_index = static_cast<double>(n - 1) * q;
    if (virtual_index >= static_cast<double>(n - 1)) {
        ranks.previous = n - 1;
        ranks.next = n - 1;
        ranks.gamma = virtual_index - static_cast<double>(n - 1);
        return ranks;
    }
    const double previous = std::floor(virtual_index);
    ranks.previous = static_cast<std::int64_t>(previous);
    ranks.next = ranks.previous + 1;
    ranks.gamma = virtual_index - previous;
    return ranks;
}

double quantile_lerp(double previous_value, double next_value, double gamma) {
    const double diff = next_value - previous_value;
    double result = previous_value + diff * gamma;
    if (gamma >= 0.5) {
        result = next_value - diff * (1.0 - gamma);
    }
    return result;
}

// ---------------------------------------------------------------------------
// Python random

PyRandom::PyRandom(std::uint64_t seed) {
    // init_genrand(19650218)
    state_[0] = 19650218U;
    for (int i = 1; i < 624; ++i) {
        state_[i] = 1812433253U * (state_[i - 1] ^ (state_[i - 1] >> 30)) + static_cast<std::uint32_t>(i);
    }
    // the seed's 32 bit words, least significant first (a seed of 0 is one word)
    std::uint32_t key[2];
    int key_length = 1;
    key[0] = static_cast<std::uint32_t>(seed & 0xFFFFFFFFULL);
    if ((seed >> 32) != 0) {
        key[1] = static_cast<std::uint32_t>(seed >> 32);
        key_length = 2;
    }
    int i = 1;
    int j = 0;
    for (int k = std::max(624, key_length); k > 0; --k) {
        state_[i] = (state_[i] ^ ((state_[i - 1] ^ (state_[i - 1] >> 30)) * 1664525U)) + key[j] +
                    static_cast<std::uint32_t>(j);
        ++i;
        ++j;
        if (i >= 624) {
            state_[0] = state_[623];
            i = 1;
        }
        if (j >= key_length) {
            j = 0;
        }
    }
    for (int k = 623; k > 0; --k) {
        state_[i] = (state_[i] ^ ((state_[i - 1] ^ (state_[i - 1] >> 30)) * 1566083941U)) -
                    static_cast<std::uint32_t>(i);
        ++i;
        if (i >= 624) {
            state_[0] = state_[623];
            i = 1;
        }
    }
    state_[0] = 0x80000000U;
    index_ = 624;
}

std::uint32_t PyRandom::next32() {
    if (index_ >= 624) {
        constexpr std::uint32_t kUpper = 0x80000000U;
        constexpr std::uint32_t kLower = 0x7fffffffU;
        constexpr std::uint32_t kMatrix = 0x9908b0dfU;
        int kk = 0;
        std::uint32_t y;
        for (; kk < 624 - 397; ++kk) {
            y = (state_[kk] & kUpper) | (state_[kk + 1] & kLower);
            state_[kk] = state_[kk + 397] ^ (y >> 1) ^ ((y & 1U) ? kMatrix : 0U);
        }
        for (; kk < 623; ++kk) {
            y = (state_[kk] & kUpper) | (state_[kk + 1] & kLower);
            state_[kk] = state_[kk + (397 - 624)] ^ (y >> 1) ^ ((y & 1U) ? kMatrix : 0U);
        }
        y = (state_[623] & kUpper) | (state_[0] & kLower);
        state_[623] = state_[396] ^ (y >> 1) ^ ((y & 1U) ? kMatrix : 0U);
        index_ = 0;
    }
    std::uint32_t y = state_[index_++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680U;
    y ^= (y << 15) & 0xefc60000U;
    y ^= (y >> 18);
    return y;
}

double PyRandom::random() {
    const std::uint32_t a = next32() >> 5;
    const std::uint32_t b = next32() >> 6;
    return (static_cast<double>(a) * 67108864.0 + static_cast<double>(b)) * (1.0 / 9007199254740992.0);
}

std::size_t PyRandom::choice_index(std::size_t size) {
    // random.choices: population[floor(random() * (n + 0.0))]
    return static_cast<std::size_t>(std::floor(random() * static_cast<double>(size)));
}

// ---------------------------------------------------------------------------
// scipy.ndimage

namespace {

enum class Boundary { Constant, Reflect };

// Index into a line of `length` samples for a position outside it, scipy's
// 'reflect' (d c b a | a b c d | d c b a).
inline std::ptrdiff_t reflect_index(std::ptrdiff_t k, std::ptrdiff_t length) {
    const std::ptrdiff_t period = 2 * length;
    std::ptrdiff_t m = k % period;
    if (m < 0) {
        m += period;
    }
    return m < length ? m : period - 1 - m;
}

// ndimage.correlate1d for a float32 image along `axis`: NI_Correlate1D
// (ni_filters.c) on a line buffer extended by size1 samples before and size2
// after, in double, with the symmetric / antisymmetric shortcuts.
ImageF correlate1d(const ImageF& in, const std::vector<double>& weights, int axis, Boundary mode,
                   double cval) {
    const std::ptrdiff_t filter_size = static_cast<std::ptrdiff_t>(weights.size());
    const std::ptrdiff_t size1 = filter_size / 2;
    const std::ptrdiff_t size2 = filter_size - size1 - 1;
    int symmetric = 0;
    if (filter_size & 1) {
        symmetric = 1;
        for (std::ptrdiff_t ii = 1; ii <= filter_size / 2; ++ii) {
            if (std::fabs(weights[ii + size1] - weights[size1 - ii]) > DBL_EPSILON) {
                symmetric = 0;
                break;
            }
        }
        if (symmetric == 0) {
            symmetric = -1;
            for (std::ptrdiff_t ii = 1; ii <= filter_size / 2; ++ii) {
                if (std::fabs(weights[size1 + ii] + weights[size1 - ii]) > DBL_EPSILON) {
                    symmetric = 0;
                    break;
                }
            }
        }
    }
    ImageF out(in.rows, in.cols);
    const std::ptrdiff_t length = axis == 0 ? in.rows : in.cols;
    const std::ptrdiff_t lines = axis == 0 ? in.cols : in.rows;
    std::vector<double> ext(static_cast<std::size_t>(length + size1 + size2));
    std::vector<double> oline(static_cast<std::size_t>(length));
    const double* fw = weights.data() + size1;  // fw[0] is the centre tap
    for (std::ptrdiff_t line = 0; line < lines; ++line) {
        auto sample = [&](std::ptrdiff_t k) -> double {
            const std::size_t r = axis == 0 ? static_cast<std::size_t>(k) : static_cast<std::size_t>(line);
            const std::size_t c = axis == 0 ? static_cast<std::size_t>(line) : static_cast<std::size_t>(k);
            return static_cast<double>(in.v[r * static_cast<std::size_t>(in.cols) + c]);
        };
        for (std::ptrdiff_t k = -size1; k < length + size2; ++k) {
            double value;
            if (k >= 0 && k < length) {
                value = sample(k);
            } else if (mode == Boundary::Constant) {
                value = cval;
            } else {
                value = sample(reflect_index(k, length));
            }
            ext[static_cast<std::size_t>(k + size1)] = value;
        }
        const double* iline = ext.data() + size1;
        for (std::ptrdiff_t ll = 0; ll < length; ++ll, ++iline) {
            double acc;
            if (symmetric > 0) {
                acc = iline[0] * fw[0];
                for (std::ptrdiff_t jj = -size1; jj < 0; ++jj) {
                    acc += (iline[jj] + iline[-jj]) * fw[jj];
                }
            } else if (symmetric < 0) {
                acc = iline[0] * fw[0];
                for (std::ptrdiff_t jj = -size1; jj < 0; ++jj) {
                    acc += (iline[jj] - iline[-jj]) * fw[jj];
                }
            } else {
                acc = iline[size2] * fw[size2];
                for (std::ptrdiff_t jj = -size1; jj < size2; ++jj) {
                    acc += iline[jj] * fw[jj];
                }
            }
            oline[static_cast<std::size_t>(ll)] = acc;
        }
        for (std::ptrdiff_t k = 0; k < length; ++k) {
            const std::size_t r = axis == 0 ? static_cast<std::size_t>(k) : static_cast<std::size_t>(line);
            const std::size_t c = axis == 0 ? static_cast<std::size_t>(line) : static_cast<std::size_t>(k);
            out.v[r * static_cast<std::size_t>(out.cols) + c] = static_cast<float>(oline[static_cast<std::size_t>(k)]);
        }
    }
    return out;
}

// ndimage.sobel(image, axis): [-1, 0, 1] along `axis`, then [1, 2, 1] along the
// other axis, mode 'reflect'.
ImageF sobel(const ImageF& in, int axis) {
    static const std::vector<double> derivative{-1.0, 0.0, 1.0};
    static const std::vector<double> smoothing{1.0, 2.0, 1.0};
    ImageF first = correlate1d(in, derivative, axis, Boundary::Reflect, 0.0);
    return correlate1d(first, smoothing, 1 - axis, Boundary::Reflect, 0.0);
}

}  // namespace

namespace {

// numpy.exp(-0.5 / 4 * x ** 2) for x = 0..8, the exponentials of Stripenn's
// default Canny sigma of 2 (radius 8). numpy 1.26 evaluates float64 exp on
// CPUs with AVX512 through Intel's SVML, which differs from the C library's
// exp in the last bit for five of these nine values; the table holds the
// values numpy returns on such a CPU (the one Stripenn was run on for the
// comparison), so that the default kernel is bit identical. Other sigmas use
// the C library's exp and can differ from numpy's kernel by one unit in the
// last place of a weight, which moves a float32 pixel of the smoothed image
// in about one image in a thousand and an edge of the Canny map in none of
// the 1200 frames images checked (cpp/scripts/stripenn_oracle/check_primitives.py).
constexpr double kExpSigma2[9] = {0x1.0000000000000p+0,  0x1.c3d6a24ed8221p-1, 0x1.368b2fc6f960ap-1,
                                  0x1.4c71b2477ab20p-2,  0x1.152aaa3bf81cbp-3, 0x1.67ee6d9ff847bp-5,
                                  0x1.6c0504695c418p-7,  0x1.1eb805a03296cp-9, 0x1.5fc21041027acp-12};

}  // namespace

std::vector<double> gaussian_kernel1d(double sigma, int radius) {
    const double sigma2 = sigma * sigma;
    std::vector<double> phi(static_cast<std::size_t>(2 * radius + 1));
    for (int i = -radius; i <= radius; ++i) {
        const double x2 = static_cast<double>(static_cast<std::int64_t>(i) * i);
        phi[static_cast<std::size_t>(i + radius)] =
            (sigma == 2.0 && radius == 8) ? kExpSigma2[std::abs(i)] : std::exp((-0.5 / sigma2) * x2);
    }
    const double total = pairwise(phi.data(), phi.size());
    for (double& v : phi) {
        v = v / total;
    }
    return phi;
}

ImageF gaussian_filter_constant(const ImageF& image, double sigma) {
    const int radius = static_cast<int>(4.0 * sigma + 0.5);
    std::vector<double> weights = gaussian_kernel1d(sigma, radius);
    std::reverse(weights.begin(), weights.end());
    ImageF pass = correlate1d(image, weights, 0, Boundary::Constant, 0.0);
    return correlate1d(pass, weights, 1, Boundary::Constant, 0.0);
}

// ---------------------------------------------------------------------------
// skimage.feature.canny

namespace {

// _canny_cy._nonmaximum_suppression_bilinear for float32 arrays, with the
// expression types of the Cython source: w and the neighbour products are
// float, (1.0 - w) is double.
ImageF nonmaximum_suppression(const ImageF& isobel, const ImageF& jsobel, const ImageF& magnitude,
                              float low_threshold) {
    const int rows = magnitude.rows;
    const int cols = magnitude.cols;
    ImageF out(rows, cols, 0.0F);
    for (int x = 1; x < rows - 1; ++x) {  // eroded_mask is 0 on the one pixel border
        for (int y = 1; y < cols - 1; ++y) {
            const float m = magnitude.at(x, y);
            if (!(m >= low_threshold)) {
                continue;
            }
            const float iso = isobel.at(x, y);
            const float jso = jsobel.at(x, y);
            const bool is_down = iso <= 0;
            const bool is_up = iso >= 0;
            const bool is_left = jso <= 0;
            const bool is_right = jso >= 0;
            const bool cond1 = (is_up && is_right) || (is_down && is_left);
            const bool cond2 = (is_down && is_right) || (is_up && is_left);
            if (!cond1 && !cond2) {
                continue;
            }
            const float abs_isobel = std::fabs(iso);
            const float abs_jsobel = std::fabs(jso);
            float w = 0.0F;
            float n11 = 0.0F;
            float n12 = 0.0F;
            float n21 = 0.0F;
            float n22 = 0.0F;
            if (cond1) {
                if (abs_isobel > abs_jsobel) {
                    w = abs_jsobel / abs_isobel;
                    n11 = magnitude.at(x + 1, y);
                    n12 = magnitude.at(x + 1, y + 1);
                    n21 = magnitude.at(x - 1, y);
                    n22 = magnitude.at(x - 1, y - 1);
                } else {
                    w = abs_isobel / abs_jsobel;
                    n11 = magnitude.at(x, y + 1);
                    n12 = magnitude.at(x + 1, y + 1);
                    n21 = magnitude.at(x, y - 1);
                    n22 = magnitude.at(x - 1, y - 1);
                }
            } else {
                if (abs_isobel < abs_jsobel) {
                    w = abs_isobel / abs_jsobel;
                    n11 = magnitude.at(x, y + 1);
                    n12 = magnitude.at(x - 1, y + 1);
                    n21 = magnitude.at(x, y - 1);
                    n22 = magnitude.at(x + 1, y - 1);
                } else {
                    w = abs_jsobel / abs_isobel;
                    n11 = magnitude.at(x - 1, y);
                    n12 = magnitude.at(x - 1, y + 1);
                    n21 = magnitude.at(x + 1, y);
                    n22 = magnitude.at(x + 1, y - 1);
                }
            }
            const float p1 = n12 * w;
            const double t1 = static_cast<double>(n11) * (1.0 - static_cast<double>(w));
            const bool c_plus = (static_cast<double>(p1) + t1) <= static_cast<double>(m);
            if (c_plus) {
                const float p2 = n22 * w;
                const double t2 = static_cast<double>(n21) * (1.0 - static_cast<double>(w));
                const bool c_minus = (static_cast<double>(p2) + t2) <= static_cast<double>(m);
                if (c_minus) {
                    out.at(x, y) = m;
                }
            }
        }
    }
    return out;
}

}  // namespace

CannyStages canny_stages(const ImageF& image, double sigma) {
    CannyStages st;
    const int rows = image.rows;
    const int cols = image.cols;
    // _preprocess: mode='constant' computes the bleed-over of the border.
    ImageF mask(rows, cols, 1.0F);
    ImageF bleed = gaussian_filter_constant(mask, sigma);
    const float eps = std::numeric_limits<float>::epsilon();
    for (float& v : bleed.v) {
        v = v + eps;
    }
    st.smoothed = gaussian_filter_constant(image, sigma);
    for (std::size_t i = 0; i < st.smoothed.v.size(); ++i) {
        st.smoothed.v[i] = st.smoothed.v[i] / bleed.v[i];
    }
    st.jsobel = sobel(st.smoothed, 1);
    st.isobel = sobel(st.smoothed, 0);
    st.magnitude = ImageF(rows, cols);
    for (std::size_t i = 0; i < st.magnitude.v.size(); ++i) {
        const float a = st.isobel.v[i] * st.isobel.v[i];
        const float b = st.jsobel.v[i] * st.jsobel.v[i];
        const float s = a + b;
        st.magnitude.v[i] = std::sqrt(s);
    }
    const float low_threshold = 0.1F;
    const float high_threshold = 0.2F;
    st.suppressed = nonmaximum_suppression(st.isobel, st.jsobel, st.magnitude, low_threshold);

    // ndi.label(low_mask, ones((3, 3))): 8-connected components; a component is
    // kept when it holds a pixel at or above the high threshold.
    const std::size_t total = st.suppressed.v.size();
    st.edges.assign(total, 0);
    std::vector<std::uint8_t> seen(total, 0);
    std::vector<std::size_t> stack;
    std::vector<std::size_t> component;
    for (std::size_t start = 0; start < total; ++start) {
        if (!(st.suppressed.v[start] > 0.0F) || seen[start]) {
            continue;
        }
        component.clear();
        stack.clear();
        stack.push_back(start);
        seen[start] = 1;
        bool high = false;
        while (!stack.empty()) {
            const std::size_t at = stack.back();
            stack.pop_back();
            component.push_back(at);
            if (st.suppressed.v[at] >= high_threshold) {
                high = true;
            }
            const int r = static_cast<int>(at / static_cast<std::size_t>(cols));
            const int c = static_cast<int>(at % static_cast<std::size_t>(cols));
            for (int dr = -1; dr <= 1; ++dr) {
                for (int dc = -1; dc <= 1; ++dc) {
                    const int rr = r + dr;
                    const int cc = c + dc;
                    if (rr < 0 || rr >= rows || cc < 0 || cc >= cols) {
                        continue;
                    }
                    const std::size_t nb = static_cast<std::size_t>(rr) * static_cast<std::size_t>(cols) +
                                           static_cast<std::size_t>(cc);
                    if (seen[nb] || !(st.suppressed.v[nb] > 0.0F)) {
                        continue;
                    }
                    seen[nb] = 1;
                    stack.push_back(nb);
                }
            }
        }
        if (high) {
            for (const std::size_t at : component) {
                st.edges[at] = 1;
            }
        }
    }
    return st;
}

std::vector<std::uint8_t> canny(const ImageF& image, double sigma) {
    return canny_stages(image, sigma).edges;
}

// ---------------------------------------------------------------------------
// OpenCV

std::vector<double> box_filter(const std::vector<double>& plane, int rows, int cols, int ksize) {
    const double kv = 1.0 / static_cast<double>(ksize * ksize);
    const int anchor = ksize / 2;
    std::vector<double> out(plane.size());
    auto reflect101 = [](int p, int len) {
        if (len == 1) {
            return 0;
        }
        while (p < 0 || p >= len) {
            if (p < 0) {
                p = -p;
            } else {
                p = 2 * (len - 1) - p;
            }
        }
        return p;
    };
    std::vector<int> row_index(static_cast<std::size_t>(ksize));
    std::vector<int> col_index(static_cast<std::size_t>(cols) * static_cast<std::size_t>(ksize));
    for (int c = 0; c < cols; ++c) {
        for (int dx = 0; dx < ksize; ++dx) {
            col_index[static_cast<std::size_t>(c) * ksize + dx] = reflect101(c + dx - anchor, cols);
        }
    }
    for (int r = 0; r < rows; ++r) {
        for (int dy = 0; dy < ksize; ++dy) {
            row_index[static_cast<std::size_t>(dy)] = reflect101(r + dy - anchor, rows);
        }
        for (int c = 0; c < cols; ++c) {
            double s = 0.0;
            for (int dy = 0; dy < ksize; ++dy) {
                const double* line = plane.data() + static_cast<std::size_t>(row_index[dy]) * cols;
                for (int dx = 0; dx < ksize; ++dx) {
                    s = std::fma(kv, line[col_index[static_cast<std::size_t>(c) * ksize + dx]], s);
                }
            }
            out[static_cast<std::size_t>(r) * cols + c] = s;
        }
    }
    return out;
}

float rgb_to_gray(float r, float g, float b) {
    const float c0 = 0.299F;
    const float c1 = 0.587F;
    const float c2 = 0.114F;
    const float inner = std::fma(r, c0, g * c1);
    return std::fma(b, c2, inner);
}

// ---------------------------------------------------------------------------
// scipy.signal

std::vector<int> convolve2d_same_3x3(const std::vector<std::uint8_t>& m, int rows, int cols,
                                     const int kernel[3][3]) {
    std::vector<int> out(static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols), 0);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int acc = 0;
            for (int u = 0; u < 3; ++u) {
                const int si = i + 1 - u;
                if (si < 0 || si >= rows) {
                    continue;
                }
                for (int v = 0; v < 3; ++v) {
                    const int sj = j + 1 - v;
                    if (sj < 0 || sj >= cols) {
                        continue;
                    }
                    acc += static_cast<int>(m[static_cast<std::size_t>(si) * cols + sj]) * kernel[u][v];
                }
            }
            out[static_cast<std::size_t>(i) * cols + j] = acc;
        }
    }
    return out;
}

}  // namespace hicx::stripes::np
