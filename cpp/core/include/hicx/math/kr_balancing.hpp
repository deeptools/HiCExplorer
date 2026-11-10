// Knight-Ruiz matrix balancing.
//
// This is a reimplementation of the algorithm krbalancing 0.0.5 implements
// (deeptools/Knight-Ruiz-Matrix-balancing-algorithm), written from that source
// as the specification of the iteration. It is deliberately not a vendoring of
// it. The upstream file is 363 lines of Eigen and carries five defects that a
// vendoring would import wholesale, all recorded as findings F1, F2, F3, F5 and
// F6 in cpp/STATUS.md:
//
//   * it downcasts the float64 input to float32 on ingestion, so it balances a
//     rounded matrix;
//   * it accumulates the two rescaling sums in float32 over every stored value,
//     inside an OpenMP loop whose body is a critical section. The critical
//     section serialises the additions without fixing their order, so the
//     result depends on thread scheduling and the tool is not reproducible run
//     to run;
//   * it stages the whole matrix through a triplet vector before
//     setFromTriplets, having already reserved the destination, copies again
//     inside setFromTriplets, again for A = A + I and again for the
//     triangularView extraction, and never releases the staging buffer;
//   * it calls exit(0) after 300 outer iterations, terminating the host process
//     with a success status and no output file;
//   * its OpenMP loops wrap their whole body in omp critical, so they are
//     serial with added contention, against a hardcoded thread count of 10.
//
// What is kept is the mathematics: the outer Newton iteration on the balancing
// equation with an inner conjugate gradient solve, the trust region tests
// against delta = 0.1 and Delta = 3, and the eta update. The matrix it balances
// is A = M + 1e-5 * I, where M is the symmetric contact matrix; the identity
// term is the upstream's way of keeping empty rows out of a division by zero
// and it is part of the result, because it puts an entry on every position of
// the diagonal of the balanced matrix.
//
// Representation. Nothing is staged and nothing is copied. The matrix stays in
// the caller's upper triangle CSR and is reached through the symmetric kernels
// of sparse_kernels.hpp, so the balancing costs the ten working vectors, that
// is eighty bytes per bin, on top of the matrix that was going to be resident
// anyway. The identity term is a scalar the matrix vector product adds, not a
// second matrix.
//
// Modes. `v4`, the default, is float64 throughout and accumulates the two
// rescaling sums per row partition and combines them in partition order, which
// makes it deterministic and thread count invariant. `v3` reproduces the
// upstream's float32 input rounding and float32 accumulators, but sums them in
// a fixed row major order rather than in an arbitrary one, because there is no
// reproducible target to match and determinism is a hard requirement
// (cpp/PLAN.md 4.1). v3 exists to measure what the float32 arithmetic costs,
// which is a number this project owes its users; it is not a default.

#ifndef HICX_MATH_KR_BALANCING_HPP
#define HICX_MATH_KR_BALANCING_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "hicx/math/sparse_kernels.hpp"

namespace hicx::kr {

struct Options {
    // v3: round every stored value to float32 before balancing
    // (krbalancing.cpp:12,27). Inert for raw integer counts below 2^24, real
    // for an already balanced or normalised float matrix.
    bool float32_input = false;
    // v3: accumulate the two rescaling sums in float32
    // (krbalancing.cpp:228-229).
    bool float32_rescale = false;
    int threads = 1;
    // krbalancing calls exit(0) here. This port stops and reports instead.
    std::int64_t max_outer_iterations = 300;
    // The identity term added to the diagonal, krbalancing.cpp:46.
    double diagonal_addend = 0.00001;
};

class Balancer {
  public:
    Balancer(kernels::DiagonalBlock& block, const Options& options);

    // computeKR. Returns false when the outer iteration did not converge within
    // max_outer_iterations, where the upstream would have exited the process
    // with status zero.
    bool compute();

    [[nodiscard]] const std::string& message() const noexcept { return message_; }
    [[nodiscard]] std::int64_t outer_iterations() const noexcept { return outer_; }

    // get_normalisation_vector(rescale). Rescaling is an in place side effect
    // on the vector, and whether it has happened by the time the caller reads
    // it is observable behaviour of hicCorrectMatrix (finding F4), so the same
    // "rescale once, on request" bookkeeping is reproduced here.
    const std::vector<double>& normalisation_vector(bool rescale);

    // get_normalised_matrix(rescale): scales every stored entry of the block by
    // x[row] * x[col] in place. The caller is responsible for the entries of
    // the diagonal that the identity term creates and the storage does not
    // hold; see hicx::correct::add_missing_diagonal.
    void normalise_matrix(bool rescale);

    // sqrt(norm_vector_sum / original_sum), the value krbalancing prints. NaN
    // until the vector has been rescaled.
    [[nodiscard]] double normalisation_factor() const noexcept { return factor_; }
    [[nodiscard]] bool rescaled() const noexcept { return rescaled_; }

  private:
    void rescale_normalisation_vector();
    void product(const std::vector<double>& in, std::vector<double>& out);

    kernels::DiagonalBlock& block_;
    Options options_;
    std::int64_t n_ = 0;
    std::vector<double> x_;
    std::vector<double> v_;
    std::vector<double> rk_;
    std::vector<double> y_;
    std::vector<double> z_;
    std::vector<double> p_;
    std::vector<double> w_;
    std::vector<double> scratch_;
    std::vector<double> ap_;
    std::vector<double> ynew_;
    double rho_km1_ = 0.0;
    double rho_km2_ = 0.0;
    double innertol_ = 0.0;
    std::int64_t outer_ = 0;
    double factor_ = 0.0;
    bool rescaled_ = false;
    std::string message_;
};

// Rounds every stored value of a matrix to float32 and back, which is what
// krbalancing does on ingestion. Exposed because hicCorrectMatrix has to undo
// it when the output keeps the raw counts.
void round_values_to_float32(kernels::DiagonalBlock& block);

}  // namespace hicx::kr

#endif  // HICX_MATH_KR_BALANCING_HPP
