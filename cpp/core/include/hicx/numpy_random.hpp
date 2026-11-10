// numpy.random.RandomState, the legacy generator, as far as scikit-learn's
// seeded estimators draw from it.
//
// `KMeans(random_state=0)` turns the integer into `np.random.RandomState(0)`
// (sklearn.utils.check_random_state) and the k-means++ seeding draws from it
// with `choice(n, p=...)` and `uniform(size=...)`. Reproducing its partition
// on real data starts with reproducing those draws bit for bit, which is what
// this header provides and nothing more.
//
// Documented behaviour of the legacy generator (numpy.random reference,
// "Legacy random generation"; Matsumoto and Nishimura 1998, ACM TOMACS 8(1),
// 3-30):
//
//  * The generator is MT19937. An integer seed is expanded with the reference
//    initialisation recurrence of the paper (multiplier 1812433253).
//  * `random_sample()` combines two 32-bit outputs a and b into 53 bits,
//    `((a >> 5) * 67108864.0 + (b >> 6)) / 9007199254740992.0`, as the
//    documentation of the legacy generator describes.
//  * `uniform(0, 1)` is `0.0 + 1.0 * random_sample()`, the same double.
//  * `choice(n, p=p)` with replacement and no size: the cumulative sum of p,
//    normalised by its last element, one `random_sample()`, and the
//    insertion point on the right (`searchsorted(side='right')`).
//
// Identity with numpy 1.26.4 is checked by test_sklearn_oracle.cpp against
// draws written by cpp/scripts/sklearn_oracle_fixture.py for 53 seeds (1500
// doubles each, so the 624-word state is regenerated twice, plus uniform and
// choice draws).

#ifndef HICX_NUMPY_RANDOM_HPP
#define HICX_NUMPY_RANDOM_HPP

#include <cstdint>
#include <array>
#include <vector>

namespace hicx::npy {

class RandomState {
  public:
    // np.random.RandomState(seed) for an integer seed in [0, 2**32).
    explicit RandomState(std::uint32_t seed);

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
    static constexpr int kStateWords = 624;

    // The next tempered 32-bit output of MT19937.
    std::uint32_t next_word();
    // Regenerates all 624 state words (the "twist" of the paper).
    void regenerate();

    std::array<std::uint32_t, kStateWords> state_{};
    int position_ = kStateWords;
};

}  // namespace hicx::npy

#endif  // HICX_NUMPY_RANDOM_HPP
