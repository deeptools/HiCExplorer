// Iterative correction, a port of hicexplorer/iterativeCorrection.py:10-86.
//
// The Python holds the matrix as a COO and repeats, up to M times:
//
//     s = W.sum(axis=1)                  marginals of the symmetric matrix
//     s = s / mean(s[s != 0])
//     total_bias *= s
//     deviation = max|s - 1|
//     s = 1 / s
//     W.data *= s[row]; W.data *= s[col]
//     stop when deviation < 1e-5
//
// and finally rescales so that the mean of the nonzero biases is one:
//
//     corr = mean(total_bias[total_bias != 0])
//     total_bias /= corr
//     W.data = W.data * corr * corr
//
// Three properties of that loop are behaviour and are reproduced exactly:
//
//  1. The marginal is a coo_matvec against a vector of ones, so each row is
//     accumulated sequentially in ascending column order. sparse_kernels.hpp
//     explains how the upper triangle storage reproduces that order.
//  2. The two multiplications by s[row] and s[col] are separate and in that
//     order. Multiplying by the product instead rounds differently.
//  3. The scaling is applied *before* the convergence test, so the pass that
//     detects convergence still scales the matrix.
//
// What is not reproduced is the structure of the loop body: this port fuses the
// scaling of one iteration with the marginal accumulation of the next, so the
// 12 bytes per stored entry are streamed once per iteration rather than twice.
// The arithmetic per entry is unchanged, and the order in which each marginal
// is accumulated is unchanged.
//
// The two guards that call exit(1) in the Python (a scaled value above 1e100
// inside the loop, above 1e10 after the final rescale) are reported to the
// caller instead, which turns them into a message and a non zero exit status.

#ifndef HICX_MATH_ICE_HPP
#define HICX_MATH_ICE_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "hicx/math/sparse_kernels.hpp"

namespace hicx::ice {

struct Options {
    // hicCorrectMatrix --iterNum, which the Python defaults to 500 and passes
    // as iterativeCorrection's M.
    std::int64_t max_iterations = 500;
    double tolerance = 1e-5;
    int threads = 1;
};

struct Result {
    // total_bias, which becomes the correction factor vector.
    std::vector<double> correction_factors;
    // Number of passes executed. The Python logs iternum + 1 on convergence,
    // an off by one in the message only; this is the number of passes.
    std::int64_t iterations = 0;
    bool converged = false;
    // Set when the Python would have called exit(1). `message` is the line it
    // logs first.
    bool failed = false;
    std::string message;
};

// Corrects `block` in place and returns the bias vector. The block is the whole
// matrix for the default path and one chromosome for --perchr.
Result correct(const kernels::DiagonalBlock& block, const Options& options);

}  // namespace hicx::ice

#endif  // HICX_MATH_ICE_HPP
