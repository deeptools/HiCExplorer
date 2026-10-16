// HiCRep's stratum-adjusted correlation coefficient (SCC), reimplemented from
// the actual hicrep 0.2.6 Python package (hicrep/hicrep.py, hicrep/utils.py),
// not from the paper alone. Reference: Yang et al., Genome Res. 2017
// Nov;27(11):1939-1949. doi: 10.1101/gr.220640.117.
//
// There is no counterpart in HiCExplorer's own Python tools; this is a
// cpp_only addition to hicCorrelate (cpp/AGENTS_CONTRACT.md tier for C++-only
// features). It is validated against the real hicrep package on real GM12878
// replicate data (see cpp/tools/hicCorrelate.cpp and the scc CLI test), not
// against a HiCExplorer reference, because none exists.
//
// Algorithm, per chromosome:
//  1. Trim the upper-triangle contact matrix to the diagonals
//     [1, nDiags), i.e. drop the main diagonal and anything nDiags or more
//     bins away (hicrep.utils.trimDiags).
//  2. Normalize by the matrix's genome-wide total contact count (the sum of
//     every stored entry in the whole file, not just this chromosome), or
//     downsample both matrices to equal total counts if requested.
//  3. If h > 0, apply hicrep's edge-corrected 2D mean filter of half-size h
//     (hicrep.utils.meanFilterSparse): a (2h+1) x (2h+1) box filter, clipped
//     at the matrix edges, divided by the true number of neighbours a cell
//     has rather than (2h+1)^2, so edge cells are not spuriously
//     down-weighted.
//  4. Per diagonal (genomic-distance stratum), a weighted Pearson correlation
//     over the union of nonzero positions in either matrix
//     (hicrep.hicrep.sccByDiag), weighted by nSamplesD * varVstran(nSamplesD)
//     with varVstran(n) = (1 + 1/n) / 12 for n >= 2, else undefined (treated
//     as zero weight).
//  5. The chromosome's SCC is the weighted average of the per-diagonal
//     correlations.
//
// Deviation from hicrep: hicrep's CLI reports one SCC score per chromosome
// and leaves combining them to the caller. hicCorrelate reports a single
// number per matrix pair, the mean over chromosomes (excluding any that come
// back NaN, e.g. because every diagonal above had fewer than two samples).
// This reduction is this port's own choice, not hicrep's.
//
// --sccDownsample reimplements hicrep.utils.resample (multinomial
// downsampling of the matrix with more contacts to the other's total) with a
// std::mt19937_64, not numpy's Mersenne Twister, so it reproduces the
// resampling *procedure* but not numpy's bit-exact draws.

#ifndef HICX_SCC_HPP
#define HICX_SCC_HPP

#include <cstdint>
#include <optional>
#include <vector>

#include "hicx/sparse_matrix.hpp"

namespace hicx::scc {

struct Options {
    int h = 0;                     // mean filter half-size; 0 disables smoothing
    std::int64_t max_distance_bp = -1;  // dBPMax; -1 means the whole chromosome
    bool downsample = false;       // resample to equal totals instead of normalizing
    std::uint64_t downsample_seed = 10;  // not numpy-equivalent, see file header
};

// One chromosome's upper-triangle entries, already localized to [0, nbins)
// row/column indices (the chromosome's own bin offset subtracted out).
struct ChromosomeEntries {
    std::int64_t nbins = 0;
    std::vector<std::int32_t> rows;
    std::vector<std::int32_t> cols;  // cols[i] > rows[i], diagonal already dropped
    std::vector<double> values;
};

// hicrep.hicrep.hicrepSCC's per-chromosome score for one chromosome, given
// its own two already-extracted entry lists and the two matrices' genome-wide
// totals (n1, n2 in hicrep). bin_size is that chromosome's bin size in bp.
// Returns NaN if no diagonal had at least two samples with nonzero weight.
double chromosome_scc(const ChromosomeEntries& e1, const ChromosomeEntries& e2,
                      double total_contacts_1, double total_contacts_2, std::int64_t bin_size,
                      const Options& options);

}  // namespace hicx::scc

#endif  // HICX_SCC_HPP
