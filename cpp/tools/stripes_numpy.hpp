// The numerical building blocks of the Stripenn method, reproduced operation
// by operation: the numpy, scipy, OpenCV and scikit-image routines that
// Stripenn 1.1.65.22 (Yoon et al. 2022) calls, in the order of floating point
// operations those libraries perform them. hicDetectStripes implements the
// Stripenn approach and its acceptance gate is identical output to Stripenn on
// the same input, so rounding is part of the specification here: a sum whose
// terms are added in another order can move a p-value comparison or a Canny
// threshold test.
//
// Every routine states the library routine it follows and the library
// version the order was read from or measured against:
//   numpy 1.26.4     pairwise summation, quantile, median, round
//   Python 3.9       random.Random (MT19937) and random.choices
//   (scipy 1.13.1 ndimage.gaussian_filter, sobel and signal.convolve2d are the
//   oracle for the filters below, which are written from their definitions)
//   scikit-image 0.24.0   feature.canny, _nonmaximum_suppression_bilinear
//   OpenCV 4.9.0     filter2D with a normalised box kernel, cvtColor RGB2GRAY
// None of those libraries is a dependency of HiCExplorer; they are only the
// external oracle the unit tests and cpp/scripts/stripenn_oracle compare to.

#ifndef HICX_TOOLS_STRIPES_NUMPY_HPP
#define HICX_TOOLS_STRIPES_NUMPY_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace hicx::stripes::np {

// numpy's pairwise summation of a contiguous double array (loops_utils.h.src
// DOUBLE_pairwise_sum): fewer than 8 terms are added in order, up to 128 terms
// with eight interleaved accumulators, more by halving at a multiple of 8.
[[nodiscard]] double pairwise_sum(const double* values, std::size_t count);

// numpy.sum of a rows x cols block cut out of a larger row-major array
// (stride between rows: row_stride). A view whose rows are not adjacent in
// memory is reduced through numpy's 8192 element buffer, one pairwise sum per
// buffer filled in row-major order and the buffers added in sequence. A single
// row or a block whose rows are adjacent (row_stride == cols) is one
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

// numpy.quantile(values, q) with the default 'linear' method given the two
// order statistics around the virtual index: sorted[k] and sorted[k + 1].
// `n` is the number of values; the function returns the interpolated value
// and reports the indices it needs through the two out parameters.
struct QuantileRanks {
    std::int64_t previous = 0;
    std::int64_t next = 0;
    double gamma = 0.0;
};
[[nodiscard]] QuantileRanks quantile_ranks(std::int64_t n, double q);
[[nodiscard]] double quantile_lerp(double previous_value, double next_value, double gamma);

// Python's random.Random(seed) for a non-negative integer seed: the MT19937
// generator initialised by init_by_array over the seed's 32 bit words, and
// random() as (a * 67108864 + b) / 2^53 with a = genrand >> 5, b = genrand >> 6.
class PyRandom {
  public:
    explicit PyRandom(std::uint64_t seed);
    [[nodiscard]] double random();
    // random.choices(population, k=1)[0] over a population of `size`
    // elements: floor(random() * size), the index of the pick.
    [[nodiscard]] std::size_t choice_index(std::size_t size);

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

// Gaussian smoothing of a float32 image with zero fill outside the image and the
// kernel truncated at 4 sigma (the result of scipy.ndimage.gaussian_filter with
// mode='constant', cval=0.0, truncate=4.0, the external oracle): one 1-D
// correlation per axis (axis 0 then axis 1), each line converted to double,
// filtered with the symmetric kernel (summation order documented at
// correlate1d in the .cpp) and stored back as float32.
[[nodiscard]] ImageF gaussian_filter_constant(const ImageF& image, double sigma);

// The sampled, normalised Gaussian kernel exp(-x^2 / (2 sigma^2)) / sum, x = -radius..radius, exposed for tests.
[[nodiscard]] std::vector<double> gaussian_kernel1d(double sigma, int radius);

// skimage.feature.canny(image, sigma=sigma) with every other argument at its
// default (thresholds 0.1 and 0.2 of the dtype maximum 1, mode='constant',
// no mask), on a float32 image. Returns 0/1 per pixel.
[[nodiscard]] std::vector<std::uint8_t> canny(const ImageF& image, double sigma);

// The stages of canny, exposed so tests compare each against skimage.
struct CannyStages {
    ImageF smoothed;
    ImageF isobel;
    ImageF jsobel;
    ImageF magnitude;
    ImageF suppressed;
    std::vector<std::uint8_t> edges;
};
[[nodiscard]] CannyStages canny_stages(const ImageF& image, double sigma);

// cv2.filter2D(img, -1, ones((k, k)) / (k * k)) on a single double channel
// with OpenCV's default border (BORDER_REFLECT_101) and anchor (centre): the
// kernel taps are accumulated in row-major order by fused multiply-add, as the
// dispatched OpenCV build does.
[[nodiscard]] std::vector<double> box_filter(const std::vector<double>& plane, int rows, int cols,
                                             int ksize);

// cv2.cvtColor(float32 RGB image, COLOR_RGB2GRAY) for one pixel:
// fma(b, 0.114f, fma(r, 0.299f, g * 0.587f)) in float32.
[[nodiscard]] float rgb_to_gray(float r, float g, float b);

// 2-D convolution of a 0/1 matrix M with a 3 x 3 integer kernel, cropped to the
// size of M with zero fill (what scipy.signal.convolve2d(M, kernel, mode='same')
// returns, the oracle), as stripenn.ImageProcessing.verticalLine uses it;
// returns an int matrix.
[[nodiscard]] std::vector<int> convolve2d_same_3x3(const std::vector<std::uint8_t>& m, int rows,
                                                   int cols, const int kernel[3][3]);

}  // namespace hicx::stripes::np

#endif  // HICX_TOOLS_STRIPES_NUMPY_HPP
