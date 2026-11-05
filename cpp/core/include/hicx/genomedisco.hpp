// GenomeDISCO (DIfferences in Smoothed COntact maps), reimplemented from the
// real kundajelab/genomedisco source (comparison_types/disco_random_walks.py,
// data_operations.py), not from the paper alone. Reference: Ursu et al.,
// Bioinformatics 2018 Oct 1;34(16):2701-2707. doi: 10.1093/bioinformatics/bty164.
//
// Unlike hicx::scc (cpp/core/include/hicx/scc.hpp), there is no pip-installable
// package to validate against: the real tool is Python 2 only and does not
// install in this project's environment. This port is instead cross-checked
// against a hand-transcribed python3 reproduction of the real source (not a
// fresh reimplementation from memory), run on real GM12878 replicate data; see
// the commit message for the values. That reproduction is not itself an
// independent, previously-published artifact the way hicrep is, so treat the
// agreement as a check that this port matches that transcription, not as
// confirmation against the published tool's own output.
//
// Algorithm, the "standard config" the GenomeDISCO paper and its own
// documentation recommend (example_parameters.txt: subsampling=lowest, tmin=3,
// tmax=3, norm=sqrtvc, removeDiag=yes, transition=yes):
//  1. The diagonal is already absent (hicx::scc::ChromosomeEntries never
//     stores it), matching removeDiag=yes.
//  2. Both chromosome matrices are (optionally) binomially downsampled to the
//     lower of the two total contact counts (data_operations.subsample_to_depth).
//  3. sqrtvc normalization: each stored value is divided by the square root
//     of its row's and column's coverage (the sum of every entry touching
//     that row, counted once per entry) in the symmetrized matrix
//     (data_operations.sqrtvc).
//  4. The normalized upper triangle is mirrored into a full symmetric matrix,
//     then row-normalized into a transition matrix (each row sums to 1;
//     empty rows are left as all-zero, matching to_transition's divide-by-1
//     rule for a zero row).
//  5. The matrix is raised to the t-th power by repeated sparse matrix
//     products (not elementwise), for t = tmin..tmax.
//  6. At each t, the score is the L1 distance between the two transition
//     matrices' t-th powers (every entry, not just the stored ones; a
//     position stored in only one side counts via the other's implicit
//     zero), divided by the average number of nonzero rows between the two
//     inputs.
//  7. The reported reproducibility is 1 minus the trapezoidal-rule average of
//     those per-t scores over t = tmin..tmax (scipy/sklearn's
//     metrics.auc(range(len(ts)), scores) / (len(ts) - 1)); with tmin == tmax
//     (the standard config) there is only one t, and the reproducibility is
//     just 1 minus that single score.
//
// Deviation from the real tool, documented the same way hicx::scc documents
// its own: GenomeDISCO's own CLI reports one score per chromosome and leaves
// combining them to the caller; hicCorrelate's --method genomedisco reports
// the mean over chromosomes, its own choice. --genomediscoSubsample uses a
// std::mt19937_64 with per-entry binomial draws, not numpy's Mersenne
// Twister, so it matches the resampling procedure but not numpy's exact draws.

#ifndef HICX_GENOMEDISCO_HPP
#define HICX_GENOMEDISCO_HPP

#include <cstdint>
#include <vector>

#include "hicx/scc.hpp"

namespace hicx::genomedisco {

struct Options {
    int tmin = 3;
    int tmax = 3;
    bool subsample = true;  // "subsampling lowest" in the standard config
    std::uint64_t subsample_seed = 7;
};

// hicx::scc::ChromosomeEntries is reused as-is: the upper-triangle entries of
// one chromosome, diagonal already excluded, local to [0, nbins).
using hicx::scc::ChromosomeEntries;

// One chromosome's GenomeDISCO reproducibility score for the given options.
// Returns NaN if either input chromosome has no stored entries.
double chromosome_score(const ChromosomeEntries& e1, const ChromosomeEntries& e2,
                        const Options& options);

}  // namespace hicx::genomedisco

#endif  // HICX_GENOMEDISCO_HPP
