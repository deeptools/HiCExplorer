#include "hicx/math/kr_balancing.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace hicx::kr {

namespace {

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

}  // namespace

void round_values_to_float32(kernels::DiagonalBlock& block) {
    double* data = block.data();
    for (std::int64_t local = 0; local < block.size(); ++local) {
        for (std::int64_t k = block.row_begin(local); k < block.row_end(local); ++k) {
            data[k] = static_cast<double>(static_cast<float>(data[k]));
        }
    }
}

Balancer::Balancer(kernels::DiagonalBlock& block, const Options& options)
    : block_(block), options_(options), n_(block.size()) {
    const std::size_t n = static_cast<std::size_t>(n_);
    x_.assign(n, 1.0);
    v_.assign(n, 0.0);
    rk_.assign(n, 0.0);
    y_.assign(n, 1.0);
    z_.assign(n, 0.0);
    p_.assign(n, 0.0);
    w_.assign(n, 0.0);
    scratch_.assign(n, 0.0);
    ap_.assign(n, 0.0);
    ynew_.assign(n, 0.0);
    factor_ = std::numeric_limits<double>::quiet_NaN();
    if (options_.float32_input) {
        round_values_to_float32(block_);
    }
}

void Balancer::product(const std::vector<double>& in, std::vector<double>& out) {
    kernels::symmetric_spmv(block_, in.data(), out.data(), options_.diagonal_addend,
                            options_.threads);
}

bool Balancer::compute() {
    const std::size_t n = static_cast<std::size_t>(n_);
    // v = x .* (A * x); rk = 1 - v; rho = rk . rk
    product(x_, scratch_);
    for (std::size_t i = 0; i < n; ++i) {
        v_[i] = x_[i] * scratch_[i];
        rk_[i] = 1.0 - v_[i];
    }
    rho_km1_ = dot(rk_, rk_);
    rho_km2_ = rho_km1_;

    const double tol = 1e-6;
    const double g = 0.9;
    const double etamax = 0.1;
    const double delta = 0.1;
    const double Delta = 3.0;
    const double stop_tol = tol * 0.5;
    const double rt = tol * tol;
    double eta = etamax;
    double rout = rho_km1_;
    double rold = rout;

    outer_ = 0;
    while (rout > rt) {
        ++outer_;
        if (outer_ > options_.max_outer_iterations) {
            // krbalancing.cpp:115-119 prints the whole x vector and calls
            // exit(0) here, which leaves the caller with a success status and
            // no output file. Deliberate deviation, cpp/STATUS.md.
            message_ = "Knight-Ruiz balancing did not converge in " +
                       std::to_string(options_.max_outer_iterations) +
                       " outer iterations. The upstream krbalancing library exits with "
                       "status 0 and writes nothing at this point; this port reports the "
                       "failure instead.";
            return false;
        }
        std::int64_t k = 0;
        std::fill(y_.begin(), y_.end(), 1.0);
        innertol_ = std::max(eta * eta * rout, rt);

        // Inner conjugate gradient solve.
        while (rho_km1_ > innertol_) {
            ++k;
            if (k == 1) {
                for (std::size_t i = 0; i < n; ++i) {
                    z_[i] = rk_[i] / v_[i];
                    p_[i] = z_[i];
                }
                rho_km1_ = dot(rk_, z_);
            } else {
                const double beta = rho_km1_ / rho_km2_;
                for (std::size_t i = 0; i < n; ++i) {
                    p_[i] = z_[i] + beta * p_[i];
                }
            }
            // w = x .* (A * (x .* p)) + v .* p
            for (std::size_t i = 0; i < n; ++i) {
                scratch_[i] = x_[i] * p_[i];
            }
            product(scratch_, w_);
            for (std::size_t i = 0; i < n; ++i) {
                w_[i] = x_[i] * w_[i] + v_[i] * p_[i];
            }
            const double alpha = rho_km1_ / dot(p_, w_);
            double smallest = std::numeric_limits<double>::infinity();
            double largest = -std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < n; ++i) {
                ap_[i] = alpha * p_[i];
                ynew_[i] = y_[i] + ap_[i];
                smallest = std::min(smallest, ynew_[i]);
                largest = std::max(largest, ynew_[i]);
            }
            if (smallest <= delta) {
                // Distance to the lower boundary of the cone, over the
                // components that move towards it.
                double gamma = std::numeric_limits<double>::infinity();
                for (std::size_t i = 0; i < n; ++i) {
                    if (ap_[i] < 0.0) {
                        gamma = std::min(gamma, (delta - y_[i]) / ap_[i]);
                    }
                }
                for (std::size_t i = 0; i < n; ++i) {
                    y_[i] += gamma * ap_[i];
                }
                break;
            }
            if (largest >= Delta) {
                double gamma = std::numeric_limits<double>::infinity();
                for (std::size_t i = 0; i < n; ++i) {
                    if (ynew_[i] > Delta) {
                        gamma = std::min(gamma, (Delta - y_[i]) / ap_[i]);
                    }
                }
                for (std::size_t i = 0; i < n; ++i) {
                    y_[i] += gamma * ap_[i];
                }
                break;
            }
            y_.swap(ynew_);
            rho_km2_ = rho_km1_;
            double rho = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                rk_[i] -= alpha * w_[i];
                z_[i] = rk_[i] / v_[i];
                rho += rk_[i] * z_[i];
            }
            rho_km1_ = rho;
        }

        for (std::size_t i = 0; i < n; ++i) {
            x_[i] *= y_[i];
        }
        product(x_, scratch_);
        for (std::size_t i = 0; i < n; ++i) {
            v_[i] = x_[i] * scratch_[i];
            rk_[i] = 1.0 - v_[i];
        }
        rho_km1_ = dot(rk_, rk_);
        rout = rho_km1_;

        const double rat = rout / rold;
        rold = rout;
        const double res_norm = std::sqrt(rout);
        const double eta_o = eta;
        eta = g * rat;
        if (g * eta_o * eta_o > 0.1) {
            eta = std::max(eta, g * eta_o * eta_o);
        }
        eta = std::max(std::min(eta, etamax), stop_tol / res_norm);
    }
    return true;
}

void Balancer::rescale_normalisation_vector() {
    const std::int64_t rows = block_.size();
    const std::int32_t* indices = block_.indices();
    const double* data = block_.data();
    const std::int64_t first = block_.first();
    const double addend = options_.diagonal_addend;

    if (options_.float32_rescale) {
        // v3: float accumulators, as krbalancing.cpp:228-229 declares them, but
        // in a fixed row major order so that the result is reproducible. The
        // upstream order is thread scheduling dependent and produces a
        // different answer on every run.
        float original_sum = 0.0F;
        float norm_sum = 0.0F;
        for (std::int64_t local = 0; local < rows; ++local) {
            const std::int64_t begin = block_.row_begin(local);
            const std::int64_t end = block_.row_end(local);
            // A's diagonal entry exists on every row, whether or not the
            // storage holds one, because of the identity term.
            double diagonal = addend;
            std::int64_t k = begin;
            if (k < end && indices[k] - first == local) {
                diagonal += data[k];
                ++k;
            }
            // The upstream accumulators are float but the terms are double, so
            // C++ evaluates each addition in double and rounds the result to
            // float. Written out so that the narrowing is visible.
            const double row_x = x_[static_cast<std::size_t>(local)];
            original_sum = static_cast<float>(original_sum + diagonal);
            norm_sum = static_cast<float>(norm_sum + diagonal * row_x * row_x);
            for (; k < end; ++k) {
                const std::int64_t column = indices[k] - first;
                original_sum = static_cast<float>(original_sum + data[k] * 2.0);
                norm_sum = static_cast<float>(
                    norm_sum + data[k] * row_x * x_[static_cast<std::size_t>(column)] * 2.0);
            }
        }
        factor_ = std::sqrt(static_cast<double>(norm_sum) / static_cast<double>(original_sum));
    } else {
        // v4: float64, accumulated once per row partition and combined in
        // partition order, so the sum is deterministic and independent of the
        // thread count while still being blocked rather than a single 61.8
        // million term sequential fold.
        const std::size_t partitions = block_.partitions();
        std::vector<double> originals(partitions, 0.0);
        std::vector<double> norms(partitions, 0.0);
        for (std::size_t p = 0; p < partitions; ++p) {
            double original_sum = 0.0;
            double norm_sum = 0.0;
            for (std::int64_t local = block_.partition()[p]; local < block_.partition()[p + 1];
                 ++local) {
                const std::int64_t begin = block_.row_begin(local);
                const std::int64_t end = block_.row_end(local);
                const double row_x = x_[static_cast<std::size_t>(local)];
                double diagonal = addend;
                std::int64_t k = begin;
                if (k < end && indices[k] - first == local) {
                    diagonal += data[k];
                    ++k;
                }
                original_sum += diagonal;
                norm_sum += diagonal * row_x * row_x;
                for (; k < end; ++k) {
                    const std::int64_t column = indices[k] - first;
                    original_sum += data[k] * 2.0;
                    norm_sum += data[k] * row_x * x_[static_cast<std::size_t>(column)] * 2.0;
                }
            }
            originals[p] = original_sum;
            norms[p] = norm_sum;
        }
        double original_sum = 0.0;
        double norm_sum = 0.0;
        for (std::size_t p = 0; p < partitions; ++p) {
            original_sum += originals[p];
            norm_sum += norms[p];
        }
        factor_ = std::sqrt(norm_sum / original_sum);
    }

    for (double& value : x_) {
        value /= factor_;
    }
    rescaled_ = true;
}

const std::vector<double>& Balancer::normalisation_vector(bool rescale) {
    if (rescale && !rescaled_) {
        rescale_normalisation_vector();
    }
    return x_;
}

void Balancer::normalise_matrix(bool rescale) {
    if (rescale && !rescaled_) {
        rescale_normalisation_vector();
    }
    kernels::scale_symmetric(block_, x_.data(), options_.threads);
}

}  // namespace hicx::kr
