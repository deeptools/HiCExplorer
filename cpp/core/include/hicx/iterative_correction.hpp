// hicexplorer.iterativeCorrection.iterativeCorrection on a general sparse
// matrix, for the public Python API.
//
// The tool path (hicx/math/ice.hpp) balances an upper triangle CSR and reads
// it as the symmetric matrix it represents. The Python function takes any
// scipy matrix, works on W = matrix.astype(float).tocoo() with every stored
// entry, duplicates and explicit zeros included, in the storage order of the
// COO, and only requires the matrix to be symmetric within a relative 1e-10.
// This kernel follows that COO literally:
//
//   s = W.sum(axis=1)          coo_matvec: s[row[k]] += data[k], k ascending
//   s = s / mean(s[s != 0])    numpy's pairwise mean
//   total_bias *= s
//   deviation = max|s - 1|     NaN propagates, as np.max does
//   s = 1 / s
//   W.data *= s[row]; W.data *= s[col]
//   stop when deviation < tolerance
//
// and the final rescale corr = mean(total_bias[total_bias != 0]),
// total_bias / corr, W.data * corr * corr. On the upper triangle storage the
// tool uses, the two compute the same values (cpp/tests/test_api_loops_correct.cpp).

#ifndef HICX_ITERATIVE_CORRECTION_HPP
#define HICX_ITERATIVE_CORRECTION_HPP

#include <cstdint>
#include <functional>
#include <vector>

namespace hicx::iterative_correction {

struct Coo {
    std::int64_t rows = 0;
    std::int64_t cols = 0;
    std::vector<std::int64_t> row;
    std::vector<std::int64_t> col;
    std::vector<double> data;
};

// np.abs(matrix - matrix.T).mean() / (1. * np.abs(matrix.mean())), the
// symmetry measure the Python compares against 1e-10.
[[nodiscard]] double asymmetry(const Coo& matrix);

enum class Outcome {
    Finished,
    // np.any(W.data > 1e100) inside the loop: the Python logs an error and
    // calls exit(1)
    LargeValuesInLoop,
    // np.any(W.data > 1e10) after the final rescale: the same
    LargeValuesAtEnd,
};

struct Result {
    Outcome outcome = Outcome::Finished;
    std::vector<double> total_bias;
    // The pass at which deviation < tolerance held (the Python's iternum), 0
    // when the loop ran to M without converging.
    std::int64_t converged_at = 0;
};

// Called after the scaling of pass `iternum` (1 based) with its deviation,
// before the convergence test, as the Python's verbose block is.
using Progress = std::function<void(std::int64_t iternum, double deviation)>;

// Balances `matrix` in place (its data) for at most `max_iterations` passes.
[[nodiscard]] Result correct(Coo& matrix, std::int64_t max_iterations, double tolerance,
                             const Progress& progress = {});

}  // namespace hicx::iterative_correction

#endif  // HICX_ITERATIVE_CORRECTION_HPP
