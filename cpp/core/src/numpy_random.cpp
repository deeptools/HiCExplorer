#include "hicx/numpy_random.hpp"

#include <algorithm>
#include <stdexcept>

namespace hicx::npy {

double RandomState::random_sample() {
    // mt19937_next_double (numpy/random/src/mt19937/mt19937.h).
    const std::uint32_t a = static_cast<std::uint32_t>(engine_()) >> 5;
    const std::uint32_t b = static_cast<std::uint32_t>(engine_()) >> 6;
    return (static_cast<double>(a) * 67108864.0 + static_cast<double>(b)) /
           9007199254740992.0;
}

std::vector<double> RandomState::uniform(std::size_t count) {
    std::vector<double> values(count);
    for (double& value : values) {
        // random_uniform(bitgen, low=0.0, range=1.0) = low + range * next_double
        value = 0.0 + 1.0 * random_sample();
    }
    return values;
}

std::int64_t RandomState::choice(const std::vector<double>& p) {
    if (p.empty()) {
        throw std::invalid_argument("'a' cannot be empty unless no samples are taken");
    }
    std::vector<double> cdf(p.size());
    double running = 0.0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        running += p[i];
        cdf[i] = running;
    }
    const double last = cdf.back();
    for (double& value : cdf) {
        value /= last;
    }
    const double u = random_sample();
    // searchsorted(side='right'): the number of cdf entries <= u.
    return static_cast<std::int64_t>(std::upper_bound(cdf.begin(), cdf.end(), u) -
                                     cdf.begin());
}

}  // namespace hicx::npy
