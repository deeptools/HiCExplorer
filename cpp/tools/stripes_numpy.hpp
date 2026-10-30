// The numerical building blocks of the Stripenn method (Yoon et al. 2022,
// Stripenn 1.1.65.22) that hicDetectStripes needs: summation, order statistics
// and a seeded random generator, a Gaussian / Sobel / Canny edge detector, a
// box filter, a grey conversion and a 3 x 3 convolution. hicDetectStripes
// implements the Stripenn method and its acceptance gate is identical output
// to Stripenn on the same input, so rounding is part of the specification: a
// sum whose terms are added in another order, or a product kept in another
// precision, can move a p-value comparison or an edge decision.
//
// All routines are written here from the published algorithms and the
// documented behaviour of the routines Stripenn calls (Matsumoto and
// Nishimura 1998 for the generator, Canny 1986 for the edge detector, the
// definitions of correlation, convolution and linear interpolation of order
// statistics). The result of each is compared bit for bit with the original
// library as a black-box oracle, never read from its source:
//   numpy 1.26.4          sum / mean / median / quantile / round
//   Python 3.9            random.Random and random.choices
//   scipy                 ndimage.gaussian_filter, sobel; signal.convolve2d
//   scikit-image 0.24.0   feature.canny
//   OpenCV 4.9.0          filter2D (normalised box kernel), cvtColor RGB2GRAY
// None of those libraries is a dependency of HiCExplorer. The comparisons are
// cpp/scripts/stripenn_oracle/check_primitives.py (randomised, large) and the
// fixtures of tests/test_detect_stripes.cpp and tests/test_stripes_primitives.cpp.

#ifndef HICX_TOOLS_STRIPES_NUMPY_HPP
#define HICX_TOOLS_STRIPES_NUMPY_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace hicx::stripes::np {

// The sum numpy.sum returns for a contiguous double array. Floating point
// addition is not associative, so the grouping is part of the result: it is
// the pairwise summation of Higham (1993, "The accuracy of floating point
// summation", SIAM J. Sci. Comput. 14:783-799) with an unrolled base case:
// fewer than 8 terms are added in order, up to 128 terms with eight
// interleaved partial sums, more by halving at a multiple of 8.
[[nodiscard]] double pairwise_sum(const double* values, std::size_t count);

// numpy.sum of a rows x cols block cut out of a larger row-major array
// (stride between rows: row_stride). numpy reduces a view whose rows are not
// adjacent in memory through a working buffer of 8192 elements (its
// documented default buffer size, numpy.getbufsize): one pairwise sum per
// buffer, filled in row-major order, and the buffer sums added in sequence. A
// single row or a block whose rows are adjacent (row_stride == cols) is one
// contiguous reduction.
[[nodiscard]] double sum_block(const double* base, std::size_t rows, std::size_t cols,
                               std::size_t row_stride);

// The same for a block that has already been copied into a compact row-major
// buffer; `contiguous_in_source` tells whether the block was contiguous in the
// array it was cut from (a single row, or full-width rows).
[[nodiscard]] double sum_compact(const double* compact, std::size_t rows, std::size_t cols,
                                 bool contiguous_in_source);
[[nodiscard]] double mean_compact(const double* compact, std::size_t rows, std::size_t cols,
                                  bool contiguous_in_source);

// numpy.mean of such a block: sum_block / (rows * cols); NaN when empty.
[[nodiscard]] double mean_block(const double* base, std::size_t rows, std::size_t cols,
                                std::size_t row_stride);

// numpy.mean(array, axis=1) for one row of `cols` elements.
[[nodiscard]] double mean_of_row(const double* row, std::size_t cols);

// numpy.mean of a 1-D array (a Python list passed to np.mean).
[[nodiscard]] double mean_of(const std::vector<double>& values);

// numpy.round / Python round: half to even.
[[nodiscard]] double round_half_even(double value);

// numpy.median and numpy.nanmedian of a 1-D array (NaN for an empty input).
[[nodiscard]] double median(std::vector<double> values);
[[nodiscard]] double nanmedian(std::vector<double> values);

// numpy.quantile(values, q) with the default 'linear' method, which is
// definition 7 of Hyndman and Fan (1996, Am. Stat. 50:361-365): the virtual
// index is (n - 1) q, and the quantile is interpolated linearly between the
// order statistics sorted[previous] and sorted[next] with weight gamma on the
// upper one. quantile_ranks gives previous, next and gamma for n values;
// quantile_lerp interpolates.
struct QuantileRanks {
    std::int64_t previous = 0;
    std::int64_t next = 0;
    double gamma = 0.0;
};
[[nodiscard]] QuantileRanks quantile_ranks(std::int64_t n, double q);
[[nodiscard]] double quantile_lerp(double previous_value, double next_value, double gamma);

// Python's random.Random(seed) for a non-negative integer seed: the Mersenne
// Twister MT19937 (Matsumoto and Nishimura 1998, ACM TOMACS 8:3-30) seeded
// with the init_by_array procedure of the authors' 2002 reference
// implementation over the seed's 32 bit words (least significant first), and
// random() as the 53 bit value (a * 2^26 + b) / 2^53 with a the upper 27 bits
// of one output and b the upper 26 bits of the next (the authors' genrand_res53).
class PyRandom {
  public:
    explicit PyRandom(std::uint64_t seed);
    [[nodiscard]] double random();
    // random.choices(population, k=1)[0] over a population of `size`
    // elements: floor(random() * size), the index of the pick.
    [[nodiscard]] std::size_t choice_index(std::size_t size);
    // random.choices(population, cum_weights=cum_weights, k=1)[0]: the index
    // found by bisecting the cumulative weights (bisect_right over the first
    // n - 1 entries) at random() times the total weight.
    [[nodiscard]] std::size_t weighted_choice_index(const std::vector<double>& cum_weights);

  private:
    std::uint32_t next32();
    std::uint32_t state_[624];
    int index_ = 624;
};

// A single channel float32 image, row major.
struct ImageF {
    int rows = 0;
    int cols = 0;
    std::vector<float> v;
    ImageF() = default;
    ImageF(int r, int c, float fill = 0.0F)
        : rows(r), cols(c), v(static_cast<std::size_t>(r) * static_cast<std::size_t>(c), fill) {}
    [[nodiscard]] float& at(int r, int c) { return v[static_cast<std::size_t>(r) * cols + c]; }
    [[nodiscard]] float at(int r, int c) const { return v[static_cast<std::size_t>(r) * cols + c]; }
};

// Samples beyond the ends of a line: a constant value, or the line mirrored
// about its ends with the end sample repeated (d c b a | a b c d | d c b a).
enum class Boundary { Constant, Reflect };

// 1-D correlation of a float32 image with `weights` along `axis` (0: down the
// columns, 1: along the rows), out[i] = sum_t weights[t] * in[i + t - n / 2],
// computed in double and stored as float32; the order of the additions is
// fixed and documented at the definition.
[[nodiscard]] ImageF correlate1d(const ImageF& in, const std::vector<double>& weights, int axis,
                                 Boundary mode, double cval);

// Gaussian smoothing of a float32 image with zero fill outside the image and the
// kernel truncated at 4 sigma (what scipy.ndimage.gaussian_filter returns with
// mode='constant', cval=0.0, truncate=4.0): one 1-D correlation per axis
// (axis 0 then axis 1), each line converted to double, filtered with the
// symmetric kernel (summation order documented at correlate1d in the .cpp) and
// stored back as float32.
[[nodiscard]] ImageF gaussian_filter_constant(const ImageF& image, double sigma);

// The sampled, normalised Gaussian kernel exp(-x^2 / (2 sigma^2)) / sum, x = -radius..radius, exposed for tests.
[[nodiscard]] std::vector<double> gaussian_kernel1d(double sigma, int radius);

// The Canny (1986) edge detector as skimage.feature.canny(image, sigma=sigma,
// low_threshold, high_threshold) computes it on a float32 image with zero fill
// outside the image, no mask and thresholds as absolute gradient magnitudes
// (defaults 0.1 and 0.2): Gaussian smoothing, Sobel gradient, non-maximum
// suppression of the gradient magnitude along the gradient direction with
// bilinear interpolation between neighbours, and hysteresis thresholding
// (a suppressed pixel above the low threshold is an edge when it is 8-connected
// through such pixels to one at or above the high threshold). The one pixel
// border of the image is never an edge. Returns 0/1 per pixel.
[[nodiscard]] std::vector<std::uint8_t> canny(const ImageF& image, double sigma,
                                              double low_threshold = 0.1,
                                              double high_threshold = 0.2);

// The stages of canny, exposed so tests compare each against skimage.
struct CannyStages {
    ImageF smoothed;
    ImageF isobel;
    ImageF jsobel;
    ImageF magnitude;
    ImageF suppressed;
    std::vector<std::uint8_t> edges;
};
[[nodiscard]] CannyStages canny_stages(const ImageF& image, double sigma, double low_threshold = 0.1,
                                       double high_threshold = 0.2);

// The normalised box filter, cv2.filter2D(img, -1, ones((k, k)) / (k * k)), on
// a single double channel: the correlation of the plane with a k x k kernel of
// weights 1 / k^2, anchored at its centre, with the default border of OpenCV
// (BORDER_REFLECT_101: gfedcb|abcdefgh|gfedcba, the edge sample not repeated).
// The taps are accumulated in row-major order, each by a fused multiply-add of
// the weight, which is how the OpenCV 4.9 build used as oracle rounds for
// kernels up to 7 x 7 (identity checked for k = 1, 3, 5, 7; OpenCV computes
// larger kernels, area 50 and above, by a Fourier transform whose rounding
// is not reproduced, so the result then agrees to about 1e-15 relative only;
// Stripenn's default is k = 3).
[[nodiscard]] std::vector<double> box_filter(const std::vector<double>& plane, int rows, int cols,
                                             int ksize);

// Grey value of one RGB pixel, cv2.cvtColor(..., COLOR_RGB2GRAY): the weights
// 0.299, 0.587, 0.114 (ITU-R BT.601 luma). Float32 pixels are combined in
// float32 as fma(b, 0.114f, fma(r, 0.299f, g * 0.587f)); 8 bit pixels with the
// weights as 15 bit fixed point integers, (9798 r + 19235 g + 3735 b + 2^14) >> 15
// (0.299 and 0.587 times 2^15, rounded; the blue weight is the remainder to 2^15),
// which rounds half up.
[[nodiscard]] float rgb_to_gray(float r, float g, float b);
[[nodiscard]] std::uint8_t rgb_to_gray_u8(std::uint8_t r, std::uint8_t g, std::uint8_t b);

// 2-D convolution of a 0/1 matrix M with a 3 x 3 integer kernel, cropped to the
// size of M with zero fill (what scipy.signal.convolve2d(M, kernel, mode='same')
// returns), as Stripenn's vertical line detection uses it; returns an int
// matrix.
[[nodiscard]] std::vector<int> convolve2d_same_3x3(const std::vector<std::uint8_t>& m, int rows,
                                                   int cols, const int kernel[3][3]);

}  // namespace hicx::stripes::np

#endif  // HICX_TOOLS_STRIPES_NUMPY_HPP
