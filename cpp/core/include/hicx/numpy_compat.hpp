// Bit exact reimplementations of the numpy and CPython behaviour that the
// textual output of the HiCExplorer tools depends on.
//
// Three things are needed:
//
//  * numpy's pairwise summation. numpy does not sum a float64 array
//    sequentially, it uses eight accumulators inside blocks of 128 elements
//    and splits recursively above that. Any other order gives a different
//    last digit, and hicInfo prints the full repr of the sum.
//  * CPython's float repr. The shortest round tripping decimal form, printed
//    in exponential notation only when the decimal point position is <= -4 or
//    > 16. This differs from printf %g and from std::to_chars(general).
//  * str() of a numpy array of strings, used by hicInfo for the list of
//    available cooler bin columns.

#ifndef HICX_NUMPY_COMPAT_HPP
#define HICX_NUMPY_COMPAT_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace hicx::npy {

// numpy's pairwise summation for contiguous float64 data, as implemented in
// numpy/core/src/umath/loops_arithm_fp.dispatch.c.src (PW_BLOCKSIZE = 128).
double pairwise_sum(const double* data, std::size_t n);

inline double pairwise_sum(const std::vector<double>& data) {
    return pairwise_sum(data.data(), data.size());
}

// The same reduction for float32 arrays. scipy keeps float32 matrices in
// float32 all the way through matrix.sum(), so the accumulation has to happen
// in single precision to reproduce the printed value.
float pairwise_sum(const float* data, std::size_t n);

inline float pairwise_sum(const std::vector<float>& data) {
    return pairwise_sum(data.data(), data.size());
}

// repr() of a Python float, which is also str() of a numpy.float64 scalar.
std::string float_repr(double value);

// "{:,}".format(i)
std::string int_with_thousands_separator(std::int64_t value);

// str() of a one dimensional numpy array of strings, for example
// ['chrom' 'start' 'end']. numpy wraps at a line width of 75 characters and
// indents continuation lines by one space.
std::string array_str(const std::vector<std::string>& items);

}  // namespace hicx::npy

#endif  // HICX_NUMPY_COMPAT_HPP
