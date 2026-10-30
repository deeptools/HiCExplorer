#include "hicx/numpy_random.hpp"

#include <algorithm>
#include <stdexcept>

namespace hicx::npy {

RandomState::RandomState(std::uint32_t seed) {
    // Matsumoto and Nishimura 1998: s[0] = seed,
    // s[i] = 1812433253 * (s[i-1] xor (s[i-1] >> 30)) + i  (mod 2^32).
    state_[0] = seed;
    for (int i = 1; i < kStateWords; ++i) {
        const std::uint32_t previous = state_[static_cast<std::size_t>(i - 1)];
        state_[static_cast<std::size_t>(i)] =
            1812433253U * (previous ^ (previous >> 30)) + static_cast<std::uint32_t>(i);
    }
    position_ = kStateWords;
}

void RandomState::regenerate() {
    // Linear recurrence x[i+624] = x[i+397] xor ((upper bit of x[i], lower 31
    // bits of x[i+1]) A), with the twist matrix A of constant 0x9908b0df.
    constexpr int kShift = 397;
    constexpr std::uint32_t kUpper = 0x80000000U;
    constexpr std::uint32_t kLower = 0x7fffffffU;
    constexpr std::uint32_t kTwist = 0x9908b0dfU;
    for (int i = 0; i < kStateWords; ++i) {
        const std::uint32_t joined =
            (state_[static_cast<std::size_t>(i)] & kUpper) |
            (state_[static_cast<std::size_t>((i + 1) % kStateWords)] & kLower);
        std::uint32_t next = state_[static_cast<std::size_t>((i + kShift) % kStateWords)] ^ (joined >> 1);
        if ((joined & 1U) != 0) {
            next ^= kTwist;
        }
        state_[static_cast<std::size_t>(i)] = next;
    }
    position_ = 0;
}

std::uint32_t RandomState::next_word() {
    if (position_ >= kStateWords) {
        regenerate();
    }
    // Tempering of the paper.
    std::uint32_t y = state_[static_cast<std::size_t>(position_++)];
    y ^= y >> 11;
    y ^= (y << 7) & 0x9d2c5680U;
    y ^= (y << 15) & 0xefc60000U;
    y ^= y >> 18;
    return y;
}

double RandomState::random_sample() {
    const std::uint32_t a = next_word() >> 5;
    const std::uint32_t b = next_word() >> 6;
    return (static_cast<double>(a) * 67108864.0 + static_cast<double>(b)) /
           9007199254740992.0;
}

std::vector<double> RandomState::uniform(std::size_t count) {
    std::vector<double> values(count);
    for (double& value : values) {
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
