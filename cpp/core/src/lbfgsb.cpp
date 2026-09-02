// Projected limited-memory BFGS, and the negative binomial fit built on it.
//
// This is the optimiser half of hicx::stats. It replaces
// scipy.optimize.fmin_l_bfgs_b as the fit_nbinom package calls it:
//
//     optim(log_likelihood, x0=initial_params, args=(X,), approx_grad=1,
//           bounds=[(eps, None), (eps, 1)])
//
// with the scipy defaults m=10, factr=1e7, pgtol=1e-5, epsilon=1e-8,
// maxiter=15000, maxfun=15000.
//
// What is faithful: the objective, the initial values, the bounds, the forward
// difference gradient and its step, the history length and both stopping
// tests. What is not: the search direction. L-BFGS-B computes a generalised
// Cauchy point along the projected steepest descent path, then minimises the
// quadratic model over the variables that are still free; this walks the two
// loop recursion over the currently free variables and takes a projected
// backtracking step. Both converge to the same minimiser of a smooth
// objective; they stop at different points inside the same tolerance. The size
// of that difference on real data is reported rather than assumed, which is
// what cpp/PLAN.md risk 6 asks for.

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "hicx/numpy_compat.hpp"
#include "hicx/stats_ops.hpp"

namespace hicx::stats {

namespace {

constexpr double kMachineEpsilon = 2.220446049250313e-16;

double clamp_to_bounds(double value, const Bound& bound) {
    return std::min(std::max(value, bound.lower), bound.upper);
}

bool at_lower(double x, const Bound& bound) {
    return x <= bound.lower;
}

bool at_upper(double x, const Bound& bound) {
    return x >= bound.upper;
}

// The infinity norm of the projected gradient, which is the pgtol test scipy
// applies: max |x - P(x - g)| over the box.
double projected_gradient_norm(const std::vector<double>& x,
                               const std::vector<double>& g,
                               std::span<const Bound> bounds) {
    double norm = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double projected = clamp_to_bounds(x[i] - g[i], bounds[i]);
        norm = std::max(norm, std::abs(x[i] - projected));
    }
    return norm;
}

}  // namespace

LbfgsbResult minimise_lbfgsb(
    const std::function<double(std::span<const double>)>& objective,
    std::span<const double> x0, std::span<const Bound> bounds,
    const LbfgsbOptions& options) {
    const std::size_t n = x0.size();
    LbfgsbResult result;
    result.x.assign(x0.begin(), x0.end());
    for (std::size_t i = 0; i < n; ++i) {
        result.x[i] = clamp_to_bounds(result.x[i], bounds[i]);
    }
    if (n == 0) {
        return result;
    }

    int evaluations = 0;
    const auto evaluate = [&](const std::vector<double>& point) {
        ++evaluations;
        return objective(point);
    };

    // approx_grad=1: forward differences with a fixed absolute step, exactly
    // as scipy's approx_fprime does it.
    std::vector<double> probe(n);
    const auto gradient = [&](const std::vector<double>& point, double value,
                              std::vector<double>& out) {
        out.assign(n, 0.0);
        probe = point;
        for (std::size_t i = 0; i < n; ++i) {
            const double saved = probe[i];
            probe[i] = saved + options.epsilon;
            out[i] = (evaluate(probe) - value) / options.epsilon;
            probe[i] = saved;
        }
    };

    double f = evaluate(result.x);
    std::vector<double> g;
    gradient(result.x, f, g);

    std::vector<std::vector<double>> s_history;
    std::vector<std::vector<double>> y_history;
    std::vector<double> rho_history;

    std::vector<double> direction(n, 0.0);
    std::vector<double> alpha(static_cast<std::size_t>(options.memory), 0.0);
    std::vector<double> candidate(n, 0.0);
    std::vector<double> next_gradient;

    for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
        if (projected_gradient_norm(result.x, g, bounds) <= options.pgtol) {
            result.status = 0;
            break;
        }
        if (evaluations >= options.max_function_evaluations) {
            result.status = 1;
            break;
        }

        // The free set: a variable pinned at a bound with the gradient pushing
        // it further out stays where it is, and its search direction is zero.
        std::vector<char> free_variable(n, 1);
        for (std::size_t i = 0; i < n; ++i) {
            if ((at_lower(result.x[i], bounds[i]) && g[i] > 0.0) ||
                (at_upper(result.x[i], bounds[i]) && g[i] < 0.0)) {
                free_variable[i] = 0;
            }
        }

        // Two loop recursion over the free variables only.
        for (std::size_t i = 0; i < n; ++i) {
            direction[i] = free_variable[i] != 0 ? -g[i] : 0.0;
        }
        const std::size_t history = s_history.size();
        for (std::size_t k = history; k-- > 0;) {
            double dot = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                if (free_variable[i] != 0) {
                    dot += s_history[k][i] * direction[i];
                }
            }
            alpha[k] = rho_history[k] * dot;
            for (std::size_t i = 0; i < n; ++i) {
                if (free_variable[i] != 0) {
                    direction[i] -= alpha[k] * y_history[k][i];
                }
            }
        }
        if (history > 0) {
            const std::size_t last = history - 1;
            double sy = 0.0;
            double yy = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                if (free_variable[i] != 0) {
                    sy += s_history[last][i] * y_history[last][i];
                    yy += y_history[last][i] * y_history[last][i];
                }
            }
            const double scale = yy > 0.0 ? sy / yy : 1.0;
            if (scale > 0.0) {
                for (std::size_t i = 0; i < n; ++i) {
                    direction[i] *= scale;
                }
            }
        }
        for (std::size_t k = 0; k < history; ++k) {
            double dot = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                if (free_variable[i] != 0) {
                    dot += y_history[k][i] * direction[i];
                }
            }
            const double beta = rho_history[k] * dot;
            for (std::size_t i = 0; i < n; ++i) {
                if (free_variable[i] != 0) {
                    direction[i] += (alpha[k] - beta) * s_history[k][i];
                }
            }
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (free_variable[i] == 0) {
                direction[i] = 0.0;
            }
        }

        double slope = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            slope += g[i] * direction[i];
        }
        if (!(slope < 0.0)) {
            // The limited-memory model is not a descent direction here, which
            // happens after a bound becomes active. Restart from steepest
            // descent, as the Fortran does when its subspace step fails.
            s_history.clear();
            y_history.clear();
            rho_history.clear();
            slope = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                direction[i] = free_variable[i] != 0 ? -g[i] : 0.0;
                slope += g[i] * direction[i];
            }
            if (!(slope < 0.0)) {
                result.status = 0;
                break;
            }
        }

        // Projected backtracking with the Armijo condition. The step is taken
        // along the projection of x + t * d onto the box, so a bound that the
        // full step would cross simply becomes active.
        constexpr double kArmijo = 1e-4;
        double step = 1.0;
        double next_f = f;
        bool improved = false;
        for (int trial = 0; trial < 60; ++trial) {
            for (std::size_t i = 0; i < n; ++i) {
                candidate[i] = clamp_to_bounds(result.x[i] + step * direction[i],
                                               bounds[i]);
            }
            next_f = evaluate(candidate);
            double moved = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                moved += g[i] * (candidate[i] - result.x[i]);
            }
            if (std::isfinite(next_f) && next_f <= f + kArmijo * moved) {
                improved = true;
                break;
            }
            if (evaluations >= options.max_function_evaluations) {
                break;
            }
            step *= 0.5;
        }
        if (!improved) {
            result.status = 2;
            break;
        }

        gradient(candidate, next_f, next_gradient);

        std::vector<double> s(n, 0.0);
        std::vector<double> y(n, 0.0);
        double sy = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            s[i] = candidate[i] - result.x[i];
            y[i] = next_gradient[i] - g[i];
            sy += s[i] * y[i];
        }
        double yy = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            yy += y[i] * y[i];
        }
        if (sy > kMachineEpsilon * std::sqrt(yy) * std::sqrt(yy)) {
            s_history.push_back(std::move(s));
            y_history.push_back(std::move(y));
            rho_history.push_back(1.0 / sy);
            if (s_history.size() > static_cast<std::size_t>(options.memory)) {
                s_history.erase(s_history.begin());
                y_history.erase(y_history.begin());
                rho_history.erase(rho_history.begin());
            }
        }

        const double previous_f = f;
        result.x = candidate;
        f = next_f;
        g = next_gradient;
        result.iterations = iteration + 1;

        // The factr test, in the form the Fortran uses:
        // (f_k - f_{k+1}) / max(|f_k|, |f_{k+1}|, 1) <= factr * eps.
        const double denominator =
            std::max({std::abs(previous_f), std::abs(f), 1.0});
        if ((previous_f - f) / denominator <= options.factr * kMachineEpsilon) {
            result.status = 0;
            break;
        }
        if (iteration + 1 == options.max_iterations) {
            result.status = 1;
        }
    }

    result.f = f;
    result.function_evaluations = evaluations;
    return result;
}

NBinomFit fit_nbinom(std::span<const double> data) {
    NBinomFit fit;
    const std::size_t n = data.size();
    if (n == 0) {
        fit.status = 2;
        return fit;
    }

    // log_likelihood, negated. The gammaln and log terms are summed with
    // numpy's pairwise order, because that is how np.sum reduces them and the
    // objective value decides where the optimiser stops.
    //
    // np.log(factorial(X)) is a constant in the parameters and scipy's
    // factorial overflows to inf above 170, which would make the objective
    // infinite for every parameter value; it is kept because the Python keeps
    // it, and the caller sees the same objective.
    std::vector<double> log_factorial(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        log_factorial[i] = std::log(std::tgamma(data[i] + 1.0));
    }
    const double log_factorial_sum = npy::pairwise_sum(log_factorial);

    std::vector<double> scratch(n, 0.0);
    const auto negative_log_likelihood =
        [&](std::span<const double> parameters) -> double {
        const double r = parameters[0];
        const double p = parameters[1];
        const double safe_p = p < 1.0 ? p : 1.0 - kMachineEpsilon;
        for (std::size_t i = 0; i < n; ++i) {
            scratch[i] = gammaln(data[i] + r);
        }
        const double gammaln_sum = npy::pairwise_sum(scratch);
        for (std::size_t i = 0; i < n; ++i) {
            scratch[i] = data[i] * std::log(1.0 - safe_p);
        }
        const double tail = npy::pairwise_sum(scratch);
        const double value = gammaln_sum - log_factorial_sum -
                             static_cast<double>(n) * gammaln(r) +
                             static_cast<double>(n) * r * std::log(p) + tail;
        return -value;
    };

    // The initial values of R's fitdistr, as the Python computes them.
    const double mean = npy::pairwise_sum(data.data(), n) / static_cast<double>(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double difference = data[i] - mean;
        scratch[i] = difference * difference;
    }
    const double variance = npy::pairwise_sum(scratch) / static_cast<double>(n);
    const double size = variance > mean ? (mean * mean) / (variance - mean) : 10.0;
    const double p0 = size + mean != 0.0 ? size / (size + mean) : size;

    const std::vector<double> start{size, p0};
    const std::vector<Bound> bounds{
        Bound{kMachineEpsilon, std::numeric_limits<double>::infinity()},
        Bound{kMachineEpsilon, 1.0}};

    const LbfgsbResult solution =
        minimise_lbfgsb(negative_log_likelihood, start, bounds, LbfgsbOptions());
    fit.size = solution.x[0];
    fit.prob = solution.x[1];
    fit.status = solution.status;
    fit.iterations = solution.iterations;
    return fit;
}

}  // namespace hicx::stats
