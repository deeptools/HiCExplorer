// QuASAR-Rep's replicate reproducibility score, reimplemented from the real
// hifive 1.5.7 Python package's own compiled core
// (hifive/libraries/_quasar.pyx, called directly through Python 2, not
// read and guessed at: hifive.Quasar._normalize, ._find_correlations and
// the module-level find_replicate_scores, each run here as this port's own
// validation oracle), not from the paper alone. Reference: Sauria & Taylor,
// Genome Biol. 2017 Aug 1;18(1):145, doi: 10.1186/s13059-017-1270-4.
//
// hifive is Python 2 only and does not install in this project's own
// pinned reference environment, but it does install (and its Cython
// extension does build and run) under a plain Python 2.7 + numpy/scipy
// environment, which is what was used to generate the reference values in
// the commit message: _quasar.filter_bins, .find_bg_dist_norm and
// .find_correlations were called directly on real GM12878/A549 data,
// bypassing hifive's own Project/Fend/HiC file-format machinery (which
// only matters for I/O, not for the statistic itself).
//
// Algorithm, per chromosome, width = 100 (hifive's own hardcoded value):
//  1. A bin is "valid" if, after iterating to a fixed point, it still has
//     at least one stored contact within `width` bins of another valid
//     bin (Quasar._normalize's filter_bins loop: a bin can lose its only
//     neighbor when that neighbor is itself dropped, so this repeats until
//     the valid set stops shrinking).
//  2. For every pair of valid bins within `width` of each other, a
//     genomic-distance background curve is built (the mean stored count at
//     each distance, across every valid pair at that distance, not just
//     the ones with data); every valid bin gets a background-divided
//     "profile": a 2 x `width` row of its contacts to the `width` nearest
//     valid bins on each side, each divided by that distance's background
//     value (Quasar._normalize's find_bg_dist_norm).
//  3. For every pair of valid bins (X, Y) within `width` of each other,
//     the correlation between X's and Y's profiles, computed only at
//     positions k shared by both (k != X, Y, each within `width` of
//     whichever of X/Y it's compared against), is stored at (X, Y-X-1);
//     pairs with fewer than 3 shared positions or a degenerate (zero)
//     variance on either side are left undefined (Quasar._find_correlations).
//  4. The replicate score between two samples is the Pearson correlation,
//     over every (row, distance) position valid in both samples, between
//     each sample's own (correlation x sqrt(1 + count)) values -- the
//     "trans" value in find_replicate_scores. A position is valid in a
//     sample only if both its own bin and the bin `distance` away are
//     themselves valid bins, and its correlation and distance weight are
//     both finite.
//
// Deviation from the real tool: QuASAR-Rep's own CLI computes this per
// chromosome (and, in its multi-resolution/multi-coverage sweep, combines
// chromosomes into a running "All" total via accumulated sums, not a mean
// of per-chromosome scores) and leaves any further combining to the
// caller; hicCorrelate's --method quasar reports the mean over chromosomes
// of each one's own score, this port's own choice, the same deviation
// hicx::scc, hicx::genomedisco and hicx::hicspector all document for
// themselves.

#ifndef HICX_QUASAR_HPP
#define HICX_QUASAR_HPP

#include <cstdint>
#include <vector>

#include "hicx/scc.hpp"

namespace hicx::quasar {

struct Options {
    int width = 100;  // hifive's own hardcoded width
};

// hicx::scc::ChromosomeEntries is reused as-is: the upper-triangle entries
// of one chromosome, local to [0, nbins), diagonal excluded (filter_bins
// and find_bg_dist_norm both skip span == 0, so a stored diagonal entry
// would never contribute to the statistic anyway). Values are treated as
// integer counts; non-integer input is rounded.
using hicx::scc::ChromosomeEntries;

// One chromosome's QuASAR-Rep replicate score between two samples. Returns
// NaN if either input is empty, the chromosomes differ in size, or no
// position is valid in both samples' correlation matrices.
double chromosome_score(const ChromosomeEntries& e1, const ChromosomeEntries& e2,
                        const Options& options);

}  // namespace hicx::quasar

#endif  // HICX_QUASAR_HPP
