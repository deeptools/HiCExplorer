// HiC-Spector's reproducibility score, reimplemented from the real
// gersteinlab/HiC-spector Python reference (run_reproducibility_v2.py,
// read in full and run as the validation oracle for this port), not from
// the paper alone. Reference: Yan, Sarkar et al., Bioinformatics 2017 Jul
// 1;33(14):2199-2201, doi: 10.1093/bioinformatics/btx152.
//
// Unlike hicx::genomedisco, the real tool's own Python reference needs
// nothing but numpy/scipy (no Python 2, no external binary), so this port
// is validated the same way hicx::scc is: against the real tool's own
// code, run on real GM12878/A549 data. See the commit message for the
// exact values.
//
// Algorithm, get_reproducibility(M1, M2, num_evec) in the real source:
//  1. A bin is dropped from consideration if, in *either* input matrix, its
//     only nonzero entry is its own diagonal (self-contacts only, no real
//     inter-bin signal) -- unless the other matrix still has signal there.
//     What survives is the shared index set `iz`.
//  2. Both matrices are restricted to iz x iz, then each is independently
//     restricted again to its own nonzero-row-sum indices and turned into
//     a symmetric normalized graph Laplacian: L = I - D^-1/2 A D^-1/2,
//     symmetrized by (L + L^T) / 2 for numerical cleanliness.
//  3. The `num_evec` smallest eigenpairs of each Laplacian are computed
//     (the real tool uses scipy.sparse.linalg.eigsh's iterative Lanczos
//     solver; this port uses LAPACK's dsyevr with an index range instead,
//     an exact dense computation rather than an iterative approximation --
//     for the modest, already-filtered sizes this runs on, this is strictly
//     more precise, not a shortcut. It does mean the cost is the dense
//     O(n^3) of a full tridiagonalization, not the sparse iterative
//     solver's cost, which matters at whole-chromosome, fine-resolution
//     scale: a 25kb human chr1 (~10,000 bins) does not finish in any
//     reasonable time, while a smaller chromosome (e.g. chr21, ~1,900
//     bins) takes well under a minute. This mirrors hicx::genomedisco's
//     own documented scalability caveat, not something specific to this
//     port's choice of solver -- the real tool's own paper and practice
//     likewise favor coarser resolutions for whole-chromosome spectral
//     analysis).
//  4. Eigenvectors are zero-padded back out to the iz-restricted index
//     space, and their inverse participation ratio (IPR = 1 / sum(v^4), a
//     measure of how delocalized an eigenvector is) is computed. Only
//     eigenvectors with IPR above --ipr_cut are kept; the number of
//     survivors in each matrix is independently counted, and the smaller
//     count is how many pairs actually get compared.
//  5. Each surviving pair's distance is
//     sqrt(min(||v1-v2||^2, ||v1+v2||^2)) (eigenvector sign is arbitrary,
//     so both relative orientations are tried). The reproducibility score
//     is |sqrt(2) - mean(distances)| / sqrt(2): 1.0 for identical inputs,
//     decreasing as the two Laplacians' low-order spectra diverge.
//
// Deviation from the real tool: HiC-Spector's own CLI reports one score
// per chromosome and leaves combining them to the caller; hicCorrelate's
// --method hicspector reports the mean over chromosomes, this port's own
// choice, the same deviation hicx::scc and hicx::genomedisco document for
// themselves.

#ifndef HICX_HICSPECTOR_HPP
#define HICX_HICSPECTOR_HPP

#include <cstdint>
#include <vector>

namespace hicx::hicspector {

struct Options {
    int num_evec = 20;     // run_reproducibility_v2.py's own hardcoded num_evec
    double ipr_cut = 5.0;  // get_reproducibility's own hardcoded ipr_cut
};

// One chromosome's entries, local to [0, nbins) row/column indices. Unlike
// hicx::scc::ChromosomeEntries, the diagonal is *not* excluded: the
// degenerate-bin filter (step 1 above) explicitly looks at each bin's own
// diagonal value, so this type keeps whatever diagonal entries the matrix
// stores (rows[i] <= cols[i], diagonal entries have rows[i] == cols[i]).
struct ChromosomeEntries {
    std::int64_t nbins = 0;
    std::vector<std::int32_t> rows;
    std::vector<std::int32_t> cols;
    std::vector<double> values;
};

// One chromosome's HiC-Spector reproducibility score. Returns NaN if either
// input is empty, the matrices differ in size, or no eigenvector survives
// the IPR cut in either matrix.
double chromosome_score(const ChromosomeEntries& e1, const ChromosomeEntries& e2,
                        const Options& options);

}  // namespace hicx::hicspector

#endif  // HICX_HICSPECTOR_HPP
