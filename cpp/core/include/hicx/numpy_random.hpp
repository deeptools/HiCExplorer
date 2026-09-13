// numpy.random.RandomState, the legacy generator, as far as scikit-learn's
// seeded estimators draw from it.
//
// `KMeans(random_state=0)` turns the integer into `np.random.RandomState(0)`
// (sklearn.utils.check_random_state) and the k-means++ seeding draws from it
// with `choice(n, p=...)` and `uniform(size=...)`. Reproducing its partition
// on real data starts with reproducing those draws bit for bit, which is what
// this header provides and nothing more.
//
// What numpy does, from numpy/random/_mt19937.pyx, mtrand.pyx and
// src/mt19937 at 1.26.4:
//
//  * An integer seed goes through `_legacy_seeding`, which calls
//    `mt19937_seed`: the reference init_genrand with multiplier 1812433253.
//    That is also what std::mt19937's integer constructor does, and the
//    tempered 32-bit outputs of the two are identical.
//  * `random_sample()` is `mt19937_next_double`: two 32-bit draws a and b,
//    `((a >> 5) * 67108864.0 + (b >> 6)) / 9007199254740992.0`. This is not
//    std::generate_canonical, which combines the words differently.
//  * `uniform(0, 1)` is `0.0 + 1.0 * random_sample()`, the same double.
//  * `choice(n, p=p)` with replacement and no size: `cdf = p.cumsum()`
//    (sequential), `cdf /= cdf[-1]`, one `random_sample()`, and
//    `cdf.searchsorted(u, side='right')`.

#ifndef HICX_NUMPY_RANDOM_HPP
#define HICX_NUMPY_RANDOM_HPP

#include <cstdint>
#include <random>
#include <vector>

namespace hicx::npy {

class RandomState {
  public:
    // np.random.RandomState(seed) for an integer seed in [0, 2**32).
    explicit RandomState(std::uint32_t seed) : engine_(seed) {}

    // RandomState.random_sample(): a double in [0, 1) with 53 random bits.
    double random_sample();

    // RandomState.uniform(size=n) with the default bounds.
    std::vector<double> uniform(std::size_t count);

    // RandomState.choice(len(p), p=p): one index drawn with probabilities p.
    // The probabilities are taken as given; numpy's validation (non-negative,
    // summing to one within sqrt(eps)) raises before anything is drawn, so a
    // caller that passes invalid p gets no draw in either implementation.
    std::int64_t choice(const std::vector<double>& p);

  private:
    std::mt19937 engine_;
};

}  // namespace hicx::npy

#endif  // HICX_NUMPY_RANDOM_HPP
