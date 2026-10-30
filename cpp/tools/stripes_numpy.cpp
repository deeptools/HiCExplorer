// See stripes_numpy.hpp. Every routine is written from the published algorithm
// or the documented behaviour named in the comment above it; the unit tests
// and cpp/scripts/stripenn_oracle/check_primitives.py compare each with the
// original library (as a black-box oracle) bit for bit.

#include "stripes_numpy.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace hicx::stripes::np {

// ---------------------------------------------------------------------------
// Summation, means and order statistics

namespace {

// Pairwise summation (Higham 1993) of n doubles with an unrolled base case,
// the grouping numpy.sum uses for a contiguous float64 array.
//  * n < 8: left to right.
//  * 8 <= n <= 128: eight partial sums over elements i, i + 8, i + 16, ...,
//    combined as ((s0 + s1) + (s2 + s3)) + ((s4 + s5) + (s6 + s7)), then the
//    n mod 8 trailing elements added in order.
//  * n > 128: the two halves (the first rounded down to a multiple of 8) are
//    summed recursively and added.
double unrolled_sum(const double* a, std::size_t n) {
    double partial[8];
    for (int j = 0; j < 8; ++j) {
        partial[j] = a[j];
    }
    const std::size_t body = n - n % 8;
    for (std::size_t i = 8; i < body; i += 8) {
        for (int j = 0; j < 8; ++j) {
            partial[j] += a[i + static_cast<std::size_t>(j)];
        }
    }
    double total = ((partial[0] + partial[1]) + (partial[2] + partial[3])) +
                   ((partial[4] + partial[5]) + (partial[6] + partial[7]));
    for (std::size_t i = body; i < n; ++i) {
        total += a[i];
    }
    return total;
}

double grouped_sum(const double* a, std::size_t n) {
    if (n < 8) {
        double total = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            total += a[i];
        }
        return total;
    }
    if (n <= 128) {
        return unrolled_sum(a, n);
    }
    std::size_t half = n / 2;
    half -= half % 8;
    return grouped_sum(a, half) + grouped_sum(a + half, n - half);
}

// numpy's reduction adds the grouped sum to the additive identity +0.0, so a
// sum of negative zeros is +0.0 (numpy 1.26; checked in check_primitives.py).
double pairwise(const double* a, std::size_t n) { return 0.0 + grouped_sum(a, n); }

// numpy reduces a non-contiguous view in chunks of its working buffer of 8192
// elements (numpy.getbufsize()).
constexpr std::size_t kReductionBuffer = 8192;

}  // namespace

double pairwise_sum(const double* values, std::size_t count) { return pairwise(values, count); }

double sum_compact(const double* compact, std::size_t rows, std::size_t cols,
                   bool contiguous_in_source) {
    if (rows == 0 || cols == 0) {
        return 0.0;
    }
    const std::size_t total = rows * cols;
    if (rows == 1 || cols == 1 || contiguous_in_source) {
        return pairwise(compact, total);
    }
    // One pairwise sum per full buffer, the first buffer's sum starting the
    // accumulation and every later one added to it.
    double accumulated = 0.0;
    for (std::size_t done = 0; done < total; done += kReductionBuffer) {
        accumulated += grouped_sum(compact + done, std::min(kReductionBuffer, total - done));
    }
    return accumulated;
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
        std::copy(base + r * row_stride, base + r * row_stride + cols, compact.begin() + static_cast<std::ptrdiff_t>(r * cols));
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
    return mean_of_row(values.data(), values.size());
}

// Round half to even (IEEE 754 roundToIntegralTiesToEven), what numpy.round
// and Python's round do for zero decimals; the default rounding mode of the
// floating point environment is to nearest even.
double round_half_even(double value) { return std::nearbyint(value); }

// The median of n values is the middle one for odd n and the mean of the two
// middle ones for even n; any NaN makes it NaN, and an empty input gives NaN.
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
    // numpy takes the mean of the one or two middle values, a sum from +0.0.
    if (n % 2 == 1) {
        return 0.0 + values[n / 2];
    }
    return (0.0 + values[n / 2 - 1] + values[n / 2]) / 2.0;
}

double nanmedian(std::vector<double> values) {
    values.erase(std::remove_if(values.begin(), values.end(), [](double v) { return std::isnan(v); }),
                 values.end());
    return median(std::move(values));
}

// Hyndman and Fan (1996) definition 7: virtual index h = (n - 1) q; the
// quantile lies between the order statistics floor(h) and floor(h) + 1 (both
// n - 1 when q reaches 1) with weight h - floor(h) on the upper one.
QuantileRanks quantile_ranks(std::int64_t n, double q) {
    QuantileRanks ranks;
    const double last = static_cast<double>(n - 1);
    const double virtual_index = last * q;
    if (virtual_index >= last) {
        ranks.previous = n - 1;
        ranks.next = n - 1;
        ranks.gamma = virtual_index - last;
        return ranks;
    }
    const double lower = std::floor(virtual_index);
    ranks.previous = static_cast<std::int64_t>(lower);
    ranks.next = ranks.previous + 1;
    ranks.gamma = virtual_index - lower;
    return ranks;
}

// Linear interpolation lo + gamma (hi - lo), evaluated from the nearer end
// (hi - (1 - gamma)(hi - lo) for gamma >= 1/2) so that it is exact at both
// ends and monotone in gamma, the usual stable form of a lerp.
double quantile_lerp(double previous_value, double next_value, double gamma) {
    const double span = next_value - previous_value;
    if (gamma >= 0.5) {
        return next_value - span * (1.0 - gamma);
    }
    return previous_value + span * gamma;
}

// ---------------------------------------------------------------------------
// Mersenne Twister MT19937 (Matsumoto and Nishimura 1998), seeded as Python's
// random.Random does for an integer: init_by_array of the authors' 2002
// reference implementation with the seed's 32 bit words as the key.

namespace {

constexpr int kStateSize = 624;
constexpr int kShift = 397;
constexpr std::uint32_t kMatrixA = 0x9908b0dfU;
constexpr std::uint32_t kUpperMask = 0x80000000U;
constexpr std::uint32_t kLowerMask = 0x7fffffffU;

}  // namespace

PyRandom::PyRandom(std::uint64_t seed) {
    // Linear congruential fill of the state (init_genrand with 19650218).
    state_[0] = 19650218U;
    for (int i = 1; i < kStateSize; ++i) {
        state_[i] = 1812433253U * (state_[i - 1] ^ (state_[i - 1] >> 30)) + static_cast<std::uint32_t>(i);
    }
    // The key: the seed's 32 bit words, least significant first (seed 0 is one word).
    std::uint32_t key[2] = {static_cast<std::uint32_t>(seed & 0xFFFFFFFFULL),
                            static_cast<std::uint32_t>(seed >> 32)};
    const int key_length = key[1] != 0 ? 2 : 1;
    // Mix the key into the state, then mix the state with itself; both passes
    // walk the state cyclically and wrap element 623 into element 0.
    int at = 1;
    int key_at = 0;
    for (int steps = std::max(kStateSize, key_length); steps > 0; --steps) {
        state_[at] = (state_[at] ^ ((state_[at - 1] ^ (state_[at - 1] >> 30)) * 1664525U)) + key[key_at] +
                     static_cast<std::uint32_t>(key_at);
        ++at;
        ++key_at;
        if (at >= kStateSize) {
            state_[0] = state_[kStateSize - 1];
            at = 1;
        }
        if (key_at >= key_length) {
            key_at = 0;
        }
    }
    for (int steps = kStateSize - 1; steps > 0; --steps) {
        state_[at] = (state_[at] ^ ((state_[at - 1] ^ (state_[at - 1] >> 30)) * 1566083941U)) -
                     static_cast<std::uint32_t>(at);
        ++at;
        if (at >= kStateSize) {
            state_[0] = state_[kStateSize - 1];
            at = 1;
        }
    }
    state_[0] = kUpperMask;  // the first element is set to 2^31 so the initial array is non-zero
    index_ = kStateSize;
}

std::uint32_t PyRandom::next32() {
    if (index_ >= kStateSize) {
        // Regenerate all 624 words in place: word i takes the upper bit of
        // itself and the lower 31 bits of the next word (cyclically), shifted
        // right by one and, when the combined low bit is set, xor-ed with the
        // matrix constant, then combined with the word kShift positions ahead.
        for (int i = 0; i < kStateSize; ++i) {
            const std::uint32_t mixed =
                (state_[i] & kUpperMask) | (state_[(i + 1) % kStateSize] & kLowerMask);
            state_[i] = state_[(i + kShift) % kStateSize] ^ (mixed >> 1) ^ ((mixed & 1U) ? kMatrixA : 0U);
        }
        index_ = 0;
    }
    // Tempering.
    std::uint32_t y = state_[index_++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680U;
    y ^= (y << 15) & 0xefc60000U;
    y ^= (y >> 18);
    return y;
}

double PyRandom::random() {
    const std::uint32_t upper = next32() >> 5;  // 27 bits
    const std::uint32_t lower = next32() >> 6;  // 26 bits
    return (static_cast<double>(upper) * 67108864.0 + static_cast<double>(lower)) * (1.0 / 9007199254740992.0);
}

std::size_t PyRandom::choice_index(std::size_t size) {
    return static_cast<std::size_t>(std::floor(random() * static_cast<double>(size)));
}

std::size_t PyRandom::weighted_choice_index(const std::vector<double>& cum_weights) {
    const double total = cum_weights.back() + 0.0;
    const double target = random() * total;
    // bisect_right over cum_weights[0 .. n - 2]: the first index whose
    // cumulative weight exceeds the target (n - 1 if none does).
    return static_cast<std::size_t>(std::upper_bound(cum_weights.begin(), cum_weights.end() - 1, target) -
                                    cum_weights.begin());
}

// ---------------------------------------------------------------------------
// Image filters (the definitions of 1-D correlation, Sobel and Gaussian smoothing)

namespace {

// Position in a line of `length` samples that stands for position k outside it
// under the 'reflect' boundary (d c b a | a b c d | d c b a): the line is
// mirrored about its ends with the end sample repeated, with period 2 * length.
inline std::ptrdiff_t reflect_index(std::ptrdiff_t k, std::ptrdiff_t length) {
    const std::ptrdiff_t period = 2 * length;
    std::ptrdiff_t m = k % period;
    if (m < 0) {
        m += period;
    }
    return m < length ? m : period - 1 - m;
}

// Mirror structure of a 1-D kernel about its centre tap (index size1 = n / 2).
enum class KernelSymmetry { None, Even, Odd };

KernelSymmetry kernel_symmetry(const std::vector<double>& w) {
    const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(w.size());
    if (n % 2 == 0) {
        return KernelSymmetry::None;
    }
    const std::ptrdiff_t centre = n / 2;
    bool even = true;
    bool odd = true;
    for (std::ptrdiff_t k = 1; k <= centre; ++k) {
        even = even && w[centre + k] == w[centre - k];
        odd = odd && w[centre + k] == -w[centre - k];
    }
    return even ? KernelSymmetry::Even : (odd ? KernelSymmetry::Odd : KernelSymmetry::None);
}

}  // namespace

// 1-D correlation of a float32 image with `weights` along `axis`. By definition
//   out[i] = sum_{t = 0}^{n-1} weights[t] * in[i + t - n / 2]
// (n = weights.size(), centre tap n / 2), samples beyond the line taken from the
// boundary mode (constant: cval; reflect: see reflect_index). Each line is
// computed in double and the result is stored as float32.
//
// The order in which the terms are added is part of the specification, since
// the output must equal the oracle's bit for bit (floating point addition is not
// associative):
//   even kernel (w[c+k] == w[c-k]): start with in[i] * w[c], then for k = c down
//     to 1 add (in[i-k] + in[i+k]) * w[c-k]; the mirrored taps are paired first;
//   odd kernel (w[c+k] == -w[c-k]): the same with (in[i-k] - in[i+k]) * w[c-k];
//   other kernels: start with the last tap in[i + n-1-c] * w[n-1], then add the
//     remaining taps in increasing offset, from t = 0 to n - 2.
ImageF correlate1d(const ImageF& in, const std::vector<double>& weights, int axis, Boundary mode,
                   double cval) {
    const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(weights.size());
    const std::ptrdiff_t before = n / 2;       // taps left of the centre
    const std::ptrdiff_t after = n - before - 1;  // taps right of the centre
    const KernelSymmetry symmetry = kernel_symmetry(weights);
    const std::ptrdiff_t length = axis == 0 ? in.rows : in.cols;
    const std::ptrdiff_t lines = axis == 0 ? in.cols : in.rows;
    const auto index_of = [&](std::ptrdiff_t line, std::ptrdiff_t k) {
        return axis == 0 ? static_cast<std::size_t>(k) * static_cast<std::size_t>(in.cols) +
                               static_cast<std::size_t>(line)
                         : static_cast<std::size_t>(line) * static_cast<std::size_t>(in.cols) +
                               static_cast<std::size_t>(k);
    };

    ImageF out(in.rows, in.cols);
    // padded[before + k] holds sample k of the current line, k = -before .. length + after - 1.
    std::vector<double> padded(static_cast<std::size_t>(length + before + after));
    for (std::ptrdiff_t line = 0; line < lines; ++line) {
        for (std::ptrdiff_t k = -before; k < length + after; ++k) {
            double value = cval;
            if (k >= 0 && k < length) {
                value = static_cast<double>(in.v[index_of(line, k)]);
            } else if (mode == Boundary::Reflect) {
                value = static_cast<double>(in.v[index_of(line, reflect_index(k, length))]);
            }
            padded[static_cast<std::size_t>(k + before)] = value;
        }
        for (std::ptrdiff_t i = 0; i < length; ++i) {
            const double* centre = padded.data() + before + i;  // centre[t] is sample i + t
            double sum;
            if (symmetry != KernelSymmetry::None) {
                sum = centre[0] * weights[static_cast<std::size_t>(before)];
                for (std::ptrdiff_t k = before; k >= 1; --k) {
                    const double pair = symmetry == KernelSymmetry::Even ? centre[-k] + centre[k]
                                                                         : centre[-k] - centre[k];
                    sum += pair * weights[static_cast<std::size_t>(before - k)];
                }
            } else {
                sum = centre[after] * weights[static_cast<std::size_t>(n - 1)];
                for (std::ptrdiff_t t = 0; t < n - 1; ++t) {
                    sum += centre[t - before] * weights[static_cast<std::size_t>(t)];
                }
            }
            out.v[index_of(line, i)] = static_cast<float>(sum);
        }
    }
    return out;
}

namespace {

// Sobel derivative along `axis`: correlation with [-1, 0, 1] along `axis`, then
// with the smoothing kernel [1, 2, 1] along the other axis, boundary 'reflect'.
ImageF sobel(const ImageF& in, int axis) {
    static const std::vector<double> derivative{-1.0, 0.0, 1.0};
    static const std::vector<double> smoothing{1.0, 2.0, 1.0};
    ImageF first = correlate1d(in, derivative, axis, Boundary::Reflect, 0.0);
    return correlate1d(first, smoothing, 1 - axis, Boundary::Reflect, 0.0);
}

}  // namespace

namespace {

// exp(-0.5 / sigma^2 * x^2) for sigma = 2 and x = 0..8 (the radius of the
// default Canny sigma of Stripenn), as numpy.exp returns it for the
// float64 array of exponents. These nine numbers are data, not code: numpy
// 1.26 evaluates float64 exp with a vectorised routine on CPUs with AVX512,
// and that routine differs from the C library's exp in the last bit for five
// of them, so the table fixes the values the comparison with Stripenn was made
// with (an AVX512 machine). cpp/scripts/stripenn_oracle/gen_exp_table.py
// regenerates the table with numpy and prints the hexadecimal floats below;
// on a CPU without AVX512 it reports which entries differ. Other sigmas use
// the C library's exp and can differ from numpy's kernel by one unit in the
// last place of a weight, which moves a float32 pixel of the smoothed image in
// about one image in a thousand and an edge of the Canny map in none of the
// 1200 frame images checked (cpp/scripts/stripenn_oracle/check_primitives.py).
constexpr double kExpSigma2[9] = {0x1.0000000000000p+0,  0x1.c3d6a24ed8221p-1, 0x1.368b2fc6f960ap-1,
                                  0x1.4c71b2477ab20p-2,  0x1.152aaa3bf81cbp-3, 0x1.67ee6d9ff847bp-5,
                                  0x1.6c0504695c418p-7,  0x1.1eb805a03296cp-9, 0x1.5fc21041027acp-12};

}  // namespace

// The Gaussian kernel exp(-x^2 / (2 sigma^2)), x = -radius..radius, scaled to
// sum 1 (the sum taken pairwise, as numpy.sum would).
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
// Canny edge detection (Canny 1986, IEEE TPAMI 8:679-698)

namespace {

// Non-maximum suppression of the gradient magnitude. A pixel survives when its
// magnitude is at least the low threshold and not smaller than the magnitude
// interpolated at the two points where the gradient line through the pixel
// crosses the ring of its eight neighbours (Canny 1986, section on
// non-maximum suppression). The gradient (g_row, g_col) lies in one of two
// diagonal quadrants: with equal signs (or a zero component) between the
// directions towards (+1, +1) and (-1, -1), otherwise between (-1, +1) and
// (+1, -1). Inside the quadrant the gradient is closer to one axis; the line
// then crosses between the axis neighbour and the diagonal neighbour on each
// side, at the fraction w = minor / major of the way to the diagonal one:
//   interpolated = w * diagonal + (1 - w) * axis neighbour.
// The precision of the arithmetic reproduces scikit-image 0.24 (the magnitude,
// w and the diagonal product are float32, the axis term is formed in double;
// a pixel with a NaN weight is dropped); ties |g_row| == |g_col| take the
// column-axis neighbour in the first quadrant and the row-axis neighbour in
// the second.
ImageF nonmaximum_suppression(const ImageF& g_row, const ImageF& g_col, const ImageF& magnitude,
                              float low_threshold) {
    const int rows = magnitude.rows;
    const int cols = magnitude.cols;
    ImageF out(rows, cols, 0.0F);
    for (int r = 1; r < rows - 1; ++r) {  // the one pixel border is never an edge
        for (int c = 1; c < cols - 1; ++c) {
            const float m = magnitude.at(r, c);
            if (!(m >= low_threshold)) {
                continue;
            }
            const float gr = g_row.at(r, c);
            const float gc = g_col.at(r, c);
            const bool toward_down_right = gr >= 0 && gc >= 0;   // gradient in (+1, +1) quadrant
            const bool toward_up_left = gr <= 0 && gc <= 0;      // gradient in (-1, -1) quadrant
            const bool toward_down_left = gr >= 0 && gc <= 0;
            const bool toward_up_right = gr <= 0 && gc >= 0;
            const bool main_diagonal = toward_down_right || toward_up_left;
            const bool anti_diagonal = toward_down_left || toward_up_right;
            if (!main_diagonal && !anti_diagonal) {
                continue;
            }
            const float abs_row = std::fabs(gr);
            const float abs_col = std::fabs(gc);
            // Offsets of the "plus" side; the "minus" side is the negation.
            // The diagonal neighbour is at (diag_row, +1). The axis neighbour
            // is at (diag_row, 0) when the gradient is closer to the row axis
            // and at (0, +1) when it is closer to the column axis.
            const int diag_row = main_diagonal ? 1 : -1;
            const bool near_row_axis = main_diagonal ? abs_row > abs_col : !(abs_row < abs_col);
            const float w = near_row_axis ? abs_col / abs_row : abs_row / abs_col;
            const int axis_row = near_row_axis ? diag_row : 0;
            const int axis_col = near_row_axis ? 0 : 1;
            const double one_minus_w = 1.0 - static_cast<double>(w);
            const auto interpolated = [&](int sign) {
                const float diagonal = magnitude.at(r + sign * diag_row, c + sign);
                const float axis = magnitude.at(r + sign * axis_row, c + sign * axis_col);
                const float weighted_diagonal = diagonal * w;
                return static_cast<double>(weighted_diagonal) + static_cast<double>(axis) * one_minus_w;
            };
            if (interpolated(1) <= static_cast<double>(m) && interpolated(-1) <= static_cast<double>(m)) {
                out.at(r, c) = m;
            }
        }
    }
    return out;
}

// Hysteresis thresholding: the edges are the pixels of `suppressed` above zero
// that are 8-connected, through such pixels, to a pixel at or above `high`.
// Region growing from those strong pixels visits every pixel once.
std::vector<std::uint8_t> hysteresis(const ImageF& suppressed, float high) {
    const int rows = suppressed.rows;
    const int cols = suppressed.cols;
    std::vector<std::uint8_t> edge(suppressed.v.size(), 0);
    std::vector<std::size_t> frontier;
    for (std::size_t i = 0; i < suppressed.v.size(); ++i) {
        if (suppressed.v[i] > 0.0F && suppressed.v[i] >= high) {
            edge[i] = 1;
            frontier.push_back(i);
        }
    }
    while (!frontier.empty()) {
        const std::size_t at = frontier.back();
        frontier.pop_back();
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
                if (!edge[nb] && suppressed.v[nb] > 0.0F) {
                    edge[nb] = 1;
                    frontier.push_back(nb);
                }
            }
        }
    }
    return edge;
}

}  // namespace

CannyStages canny_stages(const ImageF& image, double sigma, double low_threshold, double high_threshold) {
    CannyStages st;
    const int rows = image.rows;
    const int cols = image.cols;
    // Smoothing with zero fill outside the image darkens the border; dividing
    // by the smoothed all-ones image (plus one float32 epsilon) undoes it.
    const ImageF bleed = gaussian_filter_constant(ImageF(rows, cols, 1.0F), sigma);
    const float eps = std::numeric_limits<float>::epsilon();
    st.smoothed = gaussian_filter_constant(image, sigma);
    for (std::size_t i = 0; i < st.smoothed.v.size(); ++i) {
        st.smoothed.v[i] = st.smoothed.v[i] / (bleed.v[i] + eps);
    }
    // Sobel gradient along the columns (jsobel) and the rows (isobel), and its magnitude.
    st.jsobel = sobel(st.smoothed, 1);
    st.isobel = sobel(st.smoothed, 0);
    st.magnitude = ImageF(rows, cols);
    for (std::size_t i = 0; i < st.magnitude.v.size(); ++i) {
        const float row_sq = st.isobel.v[i] * st.isobel.v[i];
        const float col_sq = st.jsobel.v[i] * st.jsobel.v[i];
        st.magnitude.v[i] = std::sqrt(row_sq + col_sq);
    }
    // The thresholds are float32 in every comparison.
    st.suppressed = nonmaximum_suppression(st.isobel, st.jsobel, st.magnitude, static_cast<float>(low_threshold));
    st.edges = hysteresis(st.suppressed, static_cast<float>(high_threshold));
    return st;
}

std::vector<std::uint8_t> canny(const ImageF& image, double sigma, double low_threshold,
                                double high_threshold) {
    return canny_stages(image, sigma, low_threshold, high_threshold).edges;
}

// ---------------------------------------------------------------------------
// Box filter and grey conversion

namespace {

// Index of sample p of a line of `length` samples under the border
// reflect-101 (gfedcb|abcdefgh|gfedcba): mirrored about the end samples,
// which are not repeated. A line of one sample maps everything to it.
int reflect101(int p, int length) {
    if (length == 1) {
        return 0;
    }
    const int period = 2 * (length - 1);
    int m = p % period;
    if (m < 0) {
        m += period;
    }
    return m < length ? m : period - m;
}

}  // namespace

std::vector<double> box_filter(const std::vector<double>& plane, int rows, int cols, int ksize) {
    const double weight = 1.0 / static_cast<double>(ksize * ksize);
    const int anchor = ksize / 2;
    // The plane padded by the anchor on every side with the reflect-101 border.
    const int padded_cols = cols + ksize - 1;
    std::vector<double> padded(static_cast<std::size_t>(rows + ksize - 1) * static_cast<std::size_t>(padded_cols));
    for (int pr = 0; pr < rows + ksize - 1; ++pr) {
        const double* source = plane.data() + static_cast<std::size_t>(reflect101(pr - anchor, rows)) * cols;
        for (int pc = 0; pc < padded_cols; ++pc) {
            padded[static_cast<std::size_t>(pr) * padded_cols + pc] = source[reflect101(pc - anchor, cols)];
        }
    }
    std::vector<double> out(plane.size());
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            double sum = 0.0;
            for (int dy = 0; dy < ksize; ++dy) {
                const double* line = padded.data() + static_cast<std::size_t>(r + dy) * padded_cols + c;
                for (int dx = 0; dx < ksize; ++dx) {
                    sum = std::fma(weight, line[dx], sum);
                }
            }
            out[static_cast<std::size_t>(r) * cols + c] = sum;
        }
    }
    return out;
}

// Luma of ITU-R BT.601: 0.299 R + 0.587 G + 0.114 B.
float rgb_to_gray(float r, float g, float b) {
    const float inner = std::fma(r, 0.299F, g * 0.587F);
    return std::fma(b, 0.114F, inner);
}

// The weights as integers over 2^15: 0.299 * 32768 = 9798.1 and 0.587 * 32768
// = 19234.8 rounded to the nearest integer, and the blue weight the remainder
// that makes the three sum to exactly 2^15 (3735; 0.114 * 32768 = 3735.6 would
// round to 3736), so that white stays white. The weighted sum is rounded to
// nearest, ties upwards, by adding 2^14 before the shift.
std::uint8_t rgb_to_gray_u8(std::uint8_t r, std::uint8_t g, std::uint8_t b) {
    const std::uint32_t weighted = 9798U * r + 19235U * g + 3735U * b + (1U << 14);
    return static_cast<std::uint8_t>(weighted >> 15);
}

// ---------------------------------------------------------------------------
// 2-D convolution

// Full 2-D convolution cropped to the size of the input ('same'), zero fill:
//   out[i][j] = sum_{u,v} kernel[u][v] * m[i + 1 - u][j + 1 - v]
// (the kernel is flipped, unlike correlation). Integer arithmetic, so the
// order of the terms does not matter.
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
