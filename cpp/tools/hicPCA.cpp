// Port of hicexplorer/hicPCA.py.
//
// Computes PCA eigenvectors of a Hi-C matrix, one chromosome at a time: the
// submatrix is turned into an obs/exp matrix, its covariance is taken, and the
// requested eigenvectors of that covariance are written as bedgraph or bigWig.
//
// Two questions had to be settled empirically before this could be written,
// and both are pinned in hicexplorer/test/general/test_hicPCA.py.
//
// **Which eigenvector.** hicPCA.py:305 calls scipy.linalg.eig, the general
// non-symmetric solver, on a symmetric covariance matrix and then slices
// columns k-1 without sorting (:314-318). Measured over 26 chromosome blocks
// of the four validation configurations, dgeev's column order coincides with
// descending eigenvalue order for 23 of them and does not for three, all of
// them near degenerate or exactly degenerate spectra. The one that matters is
// chrX of small_test_matrix.h5 under --method lieberman, 4,485 bins: the two
// leading eigenvalues are 2795.01342998 and 2792.60820894, a relative gap of
// 8.6e-04, and dgeev returns them the other way round, so hicPCA's "first
// eigenvector" is the eigenvector of the *second* largest eigenvalue and the
// two are orthogonal. A symmetric solver, which returns sorted eigenpairs,
// gives a different vector there, not a slightly different one, and no
// tolerance covers that. Compatibility mode therefore calls dgeev from the
// same OpenBLAS build scipy is linked against.
//
// Calling the same dgeev is necessary and not sufficient: the matrix has to be
// the same to the bit. chrX of the same file under --method dist_norm
// --ligation_factor has a massively repeated largest eigenvalue, and eig's
// columns 0 and 1 hold eigenvalues of descending rank 169 and 179, so no rule
// on the eigenvalues describes the selection. A covariance that agreed with
// np.cov to 7.7e-14 relative made the same dgeev put the largest eigenvalue in
// column 0 and shift every column by one. See numpy_covariance, and
// test_pca_dist_norm_ligation_selects_eig_columns_not_the_largest_eigenvalues.
//
// **Which sign.** Arbitrary, and demonstrably so: running the Python on
// small_test_matrix_50kb_res.h5 with OPENBLAS_NUM_THREADS=1 and with 16 gives
// bedgraph files whose chr2L block is sign-inverted while every magnitude
// agrees to all twelve printed digits. With --extraTrack the sign is decided
// by the correlation with the track and is reproducible; without it, it is
// not. The port matches dgeev's sign by calling dgeev, and the harness
// comparison is on the magnitudes with a per chromosome global sign, which is
// the rule cpp/PLAN.md 5.4 sets.
//
// The memory rewrite. The Python densifies the chromosome block, calls
// np.corrcoef and np.cov on it and keeps five to six live copies; on
// mm9_reduced_chr1.cool, 9,760 bins from a 722 KB file, that is 4,070 MB and
// 471 s. The port holds at most two dense blocks at a time: the centred
// obs/exp block and the covariance while dsyrk runs, then the covariance and
// dgeev's eigenvector matrix. The Pearson matrix is read row by row out of
// the covariance instead of being a third block. The obs/exp block used to
// stay sparse, with the covariance taken through the identity in
// transform_ops.hpp; that was given up for the reason above, and it did not
// lower the peak, which the eigensolver phase sets.
//
// Flags the Python does not have (cpp/PLAN.md 5.8, tier 11 and
// cpp/OPTIMIZATION.md):
//   --compatMode v3  (default) dgeev, LAPACK's column order, LAPACK's signs.
//   --compatMode v4  dsyevr with range='I', only the requested eigenvectors,
//                    sorted by descending eigenvalue, sign fixed by making the
//                    largest magnitude component positive.
//   --eigenSolver dense    (default) the dense covariance and the solver
//                    --compatMode selects.
//   --eigenSolver lanczos  the requested eigenvectors of the same covariance
//                    from a Krylov solver that never forms it
//                    (hicx/sparse_pca.hpp), ordered by eigenvalue, with the
//                    sign rule of v4 before the track rule.
//   --threads N      workers. The output is byte-identical for any N.
//
// Where the threads go. Multithreaded LAPACK was measured to change the bits of
// every stage of dgeev and of dsyrk with the OpenBLAS thread count
// (hicx/dgeev_selected.hpp), and those bits decide which eigenvector the
// Python writes, so LAPACK runs on one OpenBLAS thread. What runs in parallel
// is what is exact by construction:
//   * chromosomes: each one's covariance and eigensolver is a serial
//     computation, so several chromosomes run at once and the results are
//     written in chromosome order. The dense blocks in flight are limited to
//     the size of the largest chromosome's, so the peak memory is the one a
//     serial run has;
//   * the covariance rows (numpy_covariance), as before;
//   * with lanczos, the sparse products of the Krylov solver.
// And dgeev computes only the requested eigenvectors, bit for bit the columns
// it would return (dgeev_selected_eigenvectors), which removes the
// back-transform of all the others; the columns were checked against dgeev_
// itself on mm9_reduced_chr1.cool and on both chrX blocks of the harness.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <string>
#include <thread>
#include <vector>

#include "hicx/adjust_ops.hpp"
#include "hicx/argparse.hpp"
#include "hicx/bins.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/dgeev_selected.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/matrix_ops.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/obsexp_ops.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/sparse_pca.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/transform_ops.hpp"
#include "hicx/version.hpp"

#include "bigWig.h"

extern "C" {
// OpenBLAS's CBLAS dsyrk, the routine numpy's cblas_matrixproduct calls for
// dot(X, X.T). The dependency prefix ships no cblas header, so it is declared
// here; the symbol resolves from the same libopenblasp-r0.3.28.so that
// liblapack.so points to, which is the build numpy and scipy use.
void cblas_dsyrk(int order, int uplo, int trans, int n, int k, double alpha,
                 const double* a, int lda, double beta, double* c, int ldc);
}

namespace {

// CBLAS enumerators, as numpy passes them.
constexpr int kCblasRowMajor = 101;
constexpr int kCblasUpper = 121;
constexpr int kCblasNoTrans = 111;

// Runs fn(first, last) over [0, n) in `threads` contiguous chunks. Every
// caller writes disjoint rows, so the result does not depend on the split.
template <typename Function>
void for_row_chunks(std::int64_t n, int threads, Function&& fn) {
    const int workers = std::max(1, threads);
    if (workers == 1 || n < 64) {
        fn(std::int64_t{0}, n);
        return;
    }
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(workers));
    const std::int64_t chunk = (n + workers - 1) / workers;
    for (int w = 0; w < workers; ++w) {
        const std::int64_t first = static_cast<std::int64_t>(w) * chunk;
        const std::int64_t last = std::min(n, first + chunk);
        if (first >= last) {
            break;
        }
        pool.emplace_back([&fn, first, last] { fn(first, last); });
    }
    for (std::thread& worker : pool) {
        worker.join();
    }
}

// np.cov(obs_exp_matrix_) evaluated the way numpy 1.26 evaluates it, to the
// bit, which is what hicPCA.py:301 hands to scipy.linalg.eig:
//
//     X = array(m, dtype=float64)          elementwise, exact
//     X -= X.mean(axis=1)[:, None]          pairwise row sums, then / n
//     c = dot(X, X.T)                       cblas_dsyrk(RowMajor, Upper,
//                                           NoTrans), then the upper
//                                           triangle copied to the lower
//     c *= true_divide(1, n - 1)            elementwise, exact
//
// Why bit-exact and not the sparse identity (A A^T - n m m^T) / (n - 1) of
// covariance_of_symmetric, which agrees to 7.7e-14 relative: dgeev's column
// order is a discontinuous function of its input's last bits on a degenerate
// spectrum, and hicPCA selects columns by position. Measured on
// small_test_matrix.h5 chrX under --method dist_norm --ligation_factor, whose
// largest eigenvalue 4079.5559768 is massively repeated: the same dgeev from
// the same OpenBLAS, pinned to one thread, puts 3977.6724988 and
// 3891.1597454 in columns 0 and 1 on numpy's covariance, which is what the
// Python writes, and puts 4079.5559768 in column 0 on the identity's
// covariance, shifting everything by one. Those two columns have descending
// ranks 169 and 179 of 4,485, so no ordering rule stated on the eigenvalues
// reproduces the selection; only the identical matrix does.
//
// numpy's dsyrk result depends on the BLAS thread count (measured: np.cov of
// that block differs in bits between OPENBLAS_NUM_THREADS=1 and 32), and this
// runs under pin_blas_to_one_thread, so the matrix is np.cov's at one BLAS
// thread. The Python reference runs at the default thread count; its written
// eigenvectors agree with its own one-thread run to 1e-12 on every hicPCA
// case, and on this chrX both covariances put the same eigenvalues in columns
// 0 and 1 under dgeev at 1 and at 32 threads.
//
// `block` is consumed: it is released once the dense centred copy exists, so
// the peak is two dense blocks (centred X and the covariance) during dsyrk,
// the same two the eigensolver phase holds (covariance and dgeev's vectors).
// The result is raw, NaN and infinity included, because the Pearson matrix is
// derived from the uncleaned covariance (hicPCA.py:296).
hicx::DenseSymmetric numpy_covariance(hicx::CsrMatrix& block, int threads) {
    const std::int64_t n = block.rows();
    if (n == 0) {
        block = hicx::CsrMatrix();
        return hicx::DenseSymmetric();
    }
    const std::size_t width = static_cast<std::size_t>(n);
    std::unique_ptr<double[]> centred(new double[width * width]);
    {
        const std::vector<std::int64_t>& indptr = block.indptr();
        const std::vector<std::int32_t>& indices = block.indices();
        const std::vector<double>& data = block.data();
        const double observations = static_cast<double>(n);
        for_row_chunks(n, threads, [&](std::int64_t first, std::int64_t last) {
            for (std::int64_t i = first; i < last; ++i) {
                double* row = centred.get() + static_cast<std::size_t>(i) * width;
                std::fill(row, row + width, 0.0);
                const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(i)]);
                const std::size_t end =
                    static_cast<std::size_t>(indptr[static_cast<std::size_t>(i) + 1]);
                for (std::size_t k = begin; k < end; ++k) {
                    row[static_cast<std::size_t>(indices[k])] = data[k];
                }
                // _methods._mean: umr_sum over the contiguous row, which is
                // the buffered pairwise reduction, then true_divide by n.
                const double mean = hicx::npy::pairwise_sum(row, width) / observations;
                for (std::size_t j = 0; j < width; ++j) {
                    row[j] -= mean;
                }
            }
        });
    }
    block = hicx::CsrMatrix();

    // dsyrk with beta = 0 never reads its output, and the copy below writes
    // every element of the lower triangle, so the block need not be zeroed.
    hicx::DenseSymmetric covariance = hicx::DenseSymmetric::uninitialized(n);
    const int order = static_cast<int>(n);
    cblas_dsyrk(kCblasRowMajor, kCblasUpper, kCblasNoTrans, order, order, 1.0,
                centred.get(), order, 0.0, covariance.data(), order);
    centred.reset();

    // numpy's syrk() helper copies R[i, j] into R[j, i] for j > i. Reads only
    // the upper triangle and writes only the lower one, so the rows can be
    // split freely.
    for_row_chunks(n, threads, [&](std::int64_t first, std::int64_t last) {
        for (std::int64_t i = first; i < last; ++i) {
            double* row = covariance.row(i);
            for (std::int64_t j = 0; j < i; ++j) {
                row[static_cast<std::size_t>(j)] = covariance.at(j, i);
            }
        }
    });
    const double factor = 1.0 / static_cast<double>(n - 1);
    for_row_chunks(n, threads, [&](std::int64_t first, std::int64_t last) {
        for (std::int64_t i = first; i < last; ++i) {
            double* row = covariance.row(i);
            for (std::size_t j = 0; j < width; ++j) {
                row[j] *= factor;
            }
        }
    });
    return covariance;
}

// convertNansToZeros(csr_matrix(c)).todense() followed by the same for
// infinities (hicPCA.py:303-304). csr_matrix of a dense block keeps only the
// entries that compare unequal to zero, so besides NaN and infinity becoming
// 0.0 a negative zero comes back as a positive one. Zero signs can reach
// dgeev's arithmetic, so they are normalised as well.
void clean_covariance_like_hicpca(hicx::DenseSymmetric& covariance, int threads) {
    const std::int64_t n = covariance.size();
    for_row_chunks(n, threads, [&](std::int64_t first, std::int64_t last) {
        for (std::int64_t i = first; i < last; ++i) {
            double* row = covariance.row(i);
            for (std::int64_t j = 0; j < n; ++j) {
                double& value = row[static_cast<std::size_t>(j)];
                if (!std::isfinite(value) || value == 0.0) {
                    value = 0.0;
                }
            }
        }
    });
}

const char* const kUsage =
    "usage: hicPCA --matrix MATRIX --outputFileName OUTPUTFILENAME\n"
    "              [OUTPUTFILENAME ...]\n"
    "              [--whichEigenvectors WHICHEIGENVECTORS [WHICHEIGENVECTORS ...]]\n"
    "              [--format {bedgraph,bigwig}]\n"
    "              [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]\n"
    "              [--method {dist_norm,lieberman}] [--ligation_factor]\n"
    "              [--extraTrack EXTRATRACK] [--histonMarkType HISTONMARKTYPE]\n"
    "              [--pearsonMatrix PEARSONMATRIX] [--obsexpMatrix OBSEXPMATRIX]\n"
    "              [--ignoreMaskedBins] [--compatMode {v3,v4}]\n"
    "              [--eigenSolver {dense,lanczos}] [--threads THREADS]\n"
    "              [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Computes PCA eigenvectors for a Hi-C matrix.\n"
    "\n"
    "    $ hicPCA --matrix hic_matrix.h5 -o pca1.bedgraph pca2.bedgraph\n"
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        HiCExplorer matrix in h5 format.\n"
    "  --outputFileName OUTPUTFILENAME [OUTPUTFILENAME ...], -o ...\n"
    "                        File names for the result of the pca. Number of output\n"
    "                        files must match the number of computed eigenvectors.\n"
    "\n"
    "Optional arguments:\n"
    "  --whichEigenvectors WHICHEIGENVECTORS [...], -we ...\n"
    "                        The list of eigenvectors that the PCA should compute\n"
    "                        e.g. 1 2 5 will return the first, second and fifth\n"
    "                        eigenvector. (Default: 1 2).\n"
    "  --format {bedgraph,bigwig}, -f {bedgraph,bigwig}\n"
    "                        Output format. (Default: bigwig).\n"
    "  --chromosomes CHROMOSOMES [CHROMOSOMES ...]\n"
    "                        List of chromosomes to be included in the correlation.\n"
    "  --method {dist_norm,lieberman}\n"
    "                        Method used to build the obs-exp matrix.\n"
    "                        (Default: dist_norm).\n"
    "  --ligation_factor     Multiply a scaling factor to each entry of the expected\n"
    "                        matrix, as the Homer software does. Only effective with\n"
    "                        the dist_norm method.\n"
    "  --extraTrack EXTRATRACK\n"
    "                        A gene track (bed) or a histone mark coverage file\n"
    "                        (bigwig) used to decide the sign of the eigenvector.\n"
    "  --histonMarkType HISTONMARKTYPE\n"
    "                        active or inactive. (Default: active).\n"
    "  --pearsonMatrix PEARSONMATRIX, -pm PEARSONMATRIX\n"
    "                        Write the intermediate Pearson matrix to this file.\n"
    "  --obsexpMatrix OBSEXPMATRIX, -oem OBSEXPMATRIX\n"
    "                        Write the intermediate observed/expected matrix here.\n"
    "  --ignoreMaskedBins    Remove the masked bins before the PCA is computed.\n"
    "  --compatMode {v3,v4}  v3 reproduces scipy.linalg.eig, the general solver, in\n"
    "                        LAPACK's own column order and sign, which is what the\n"
    "                        Python does. v4 uses the symmetric solver dsyevr for\n"
    "                        the requested eigenvectors only, sorted by descending\n"
    "                        eigenvalue with a deterministic sign. Not a Python\n"
    "                        option; see cpp/PLAN.md 5.8. (Default: v3).\n"
    "  --eigenSolver {dense,lanczos}\n"
    "                        dense forms each chromosome's covariance and solves it\n"
    "                        as --compatMode says. lanczos computes the requested\n"
    "                        eigenvectors of the same covariance from the sparse\n"
    "                        obs/exp matrix without forming it, ordered by\n"
    "                        eigenvalue, which is not always the Python's order.\n"
    "                        Not a Python option; see cpp/PLAN.md tier 11.\n"
    "                        (Default: dense).\n"
    "  --threads THREADS     Workers for chromosomes, covariance rows and sparse\n"
    "                        products. The output is byte-identical for any value.\n"
    "                        (Default: 4).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::string matrix;
    std::vector<std::string> output_file_names;
    // The Python default is the *string* '1 2' rather than a list, and
    // argparse with nargs='+' leaves it as that string. len() of it is 3, so
    // the count check at hicPCA.py:242 compares 3 against the number of output
    // files and iterating it yields the characters '1', ' ' and '2'. Kept as
    // the three one-character tokens so that both behaviours follow.
    std::vector<std::string> which_eigenvectors{"1", " ", "2"};
    std::string format = "bigwig";
    std::vector<std::string> chromosomes;
    std::string method = "dist_norm";
    bool ligation_factor = false;
    std::optional<std::string> extra_track;
    std::string histone_mark_type = "active";
    std::optional<std::string> pearson_matrix;
    std::optional<std::string> obsexp_matrix;
    bool ignore_masked_bins = false;
    bool compat_v4 = false;
    bool lanczos = false;
    int threads = 4;
};

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// hicPCA.py parse_arguments, plus the C++-only --compatMode and --threads.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicPCA", "Computes PCA eigenvectors for a Hi-C matrix.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrix", "-m"})
        .required()
        .input({"h5", "cool", "mcool"})
        .help("HiCExplorer matrix in h5 format.");
    required.add({"--outputFileName", "-o"})
        .nargs("+")
        .required()
        .output({"bedgraph", "bigwig"})
        .help("File names for the result of the pca. Number of output files must match the "
              "number of computed eigenvectors.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--whichEigenvectors", "-we"})
        .default_value("1 2")
        .nargs("+")
        .required(false)
        .help("The list of eigenvectors that the PCA should compute e.g. 1 2 5 will return the "
              "first, second and fifth eigenvector.");
    optional.add({"--format", "-f"})
        .choices({"bedgraph", "bigwig"})
        .default_value("bigwig")
        .required(false)
        .help("Output format. Either bedgraph or bigwig.");
    optional.add({"--chromosomes"})
        .nargs("+")
        .help("List of chromosomes to be included in the correlation.");
    optional.add({"--method"})
        .choices({"dist_norm", "lieberman"})
        .default_value("dist_norm")
        .required(false)
        .help("possible methods which can be used to build the obs-exp matrix are dist_norm "
              "and lieberman.");
    optional.add({"--ligation_factor"})
        .action(cli::Action::StoreTrue)
        .help("Multiply a scaling factor to each entry of the expected matrix, as the Homer "
              "software does. Only effective with the dist_norm method.");
    optional.add({"--extraTrack"})
        .input({"bed", "bigwig"})
        .help("Either a gene track (bed) or a histone mark coverage file (bigwig) used to decide "
              "the sign of the eigenvector.");
    optional.add({"--histonMarkType"})
        .default_value("active")
        .help("Set it to active or inactive.");
    optional.add({"--pearsonMatrix", "-pm"})
        .output({"h5", "cool"})
        .help("Write the intermediate Pearson matrix to this file.");
    optional.add({"--obsexpMatrix", "-oem"})
        .output({"h5", "cool"})
        .help("Write the intermediate observed/expected matrix to this file.");
    optional.add({"--ignoreMaskedBins"})
        .action(cli::Action::StoreTrue)
        .help("Remove the masked bins before the PCA is computed.");
    optional.add({"--compatMode"})
        .choices({"v3", "v4"})
        .default_value("v3")
        .cpp_only("v4 replaces the general eigensolver by the symmetric one for the requested "
                  "eigenvectors only (cpp/PLAN.md 5.8).")
        .help("v3 reproduces scipy.linalg.eig; v4 uses the symmetric solver dsyevr.");
    optional.add({"--eigenSolver"})
        .choices({"dense", "lanczos"})
        .default_value("dense")
        .cpp_only("lanczos computes the requested eigenvectors without forming the dense "
                  "covariance (cpp/PLAN.md tier 11).")
        .help("dense forms the covariance; lanczos uses a Krylov solver on the sparse matrix.");
    optional.add({"--threads", "-t"})
        .type("int")
        .default_value(4)
        .cpp_only("The Python computes the covariance single-threaded; the output is "
                  "byte-identical for any value.")
        .help("Workers for chromosomes, covariance rows and sparse products.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show the help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.matrix = ns.str("matrix");
    args.output_file_names = ns.strs("outputFileName");
    if (ns.given("whichEigenvectors")) {
        args.which_eigenvectors = ns.strs("whichEigenvectors");
    } else {
        // The default is the string '1 2'; see the comment on the member.
        args.which_eigenvectors.clear();
        for (const char character : ns.str("whichEigenvectors")) {
            args.which_eigenvectors.emplace_back(1, character);
        }
    }
    args.format = ns.str("format");
    args.chromosomes = ns.strs("chromosomes");
    args.method = ns.str("method");
    args.ligation_factor = ns.flag("ligation_factor");
    args.extra_track = ns.opt_str("extraTrack");
    args.histone_mark_type = ns.str("histonMarkType");
    args.pearson_matrix = ns.opt_str("pearsonMatrix");
    args.obsexp_matrix = ns.opt_str("obsexpMatrix");
    args.ignore_masked_bins = ns.flag("ignoreMaskedBins");
    args.compat_v4 = ns.str("compatMode") == "v4";
    args.lanczos = ns.str("eigenSolver") == "lanczos";
    args.threads = static_cast<int>(ns.integer("threads"));
    if (args.threads < 1) {
        std::fputs(kUsage, stderr);
        std::fputs("hicPCA: error: argument --threads: must be at least 1\n", stderr);
        std::exit(2);
    }
    if (args.lanczos && args.compat_v4) {
        std::fputs(kUsage, stderr);
        std::fputs("hicPCA: error: argument --eigenSolver: lanczos replaces the dense solver "
                   "that --compatMode selects; --compatMode v4 cannot be combined with it\n",
                   stderr);
        std::exit(2);
    }
    return args;
}

// hicPCA.py:367 and :401 write with '{:.12f}', twelve decimal places.
std::string format_twelve_decimals(double value) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.12f", value);
    return std::string(buffer);
}

// A CSR matrix assembled one complete row at a time, in increasing row order.
// This is what replaces the lil_matrix accumulators of hicPCA.py:269 and :272:
// no per-element object, no second copy at .tocsr() time, and the values are
// appended as they are produced (cpp/PLAN.md 4.4 rules 4 and 7).
class RowwiseCsrBuilder {
  public:
    explicit RowwiseCsrBuilder(std::int64_t n) : n_(n), indptr_(1, 0) {}

    void skip_to(std::int64_t row) {
        while (next_row_ < row) {
            indptr_.push_back(static_cast<std::int64_t>(data_.size()));
            ++next_row_;
        }
    }
    void push(std::int64_t column, double value) {
        indices_.push_back(static_cast<std::int32_t>(column));
        data_.push_back(value);
    }
    void end_row() {
        indptr_.push_back(static_cast<std::int64_t>(data_.size()));
        ++next_row_;
    }

    // Consecutive complete rows collected away from the builder, so that
    // chromosomes can be computed out of order and appended in order.
    struct Rows {
        std::vector<std::int64_t> counts;
        std::vector<std::int32_t> indices;
        std::vector<double> data;
        void push(std::int64_t column, double value) {
            indices.push_back(static_cast<std::int32_t>(column));
            data.push_back(value);
        }
        void end_row(std::size_t row_start) {
            counts.push_back(static_cast<std::int64_t>(data.size() - row_start));
        }
    };

    // Appends `rows` as rows first, first + 1, ... The first block appended to
    // an empty builder is moved in rather than copied.
    void append(std::int64_t first, Rows&& rows) {
        skip_to(first);
        std::int64_t cursor = static_cast<std::int64_t>(data_.size());
        if (data_.empty()) {
            indices_ = std::move(rows.indices);
            data_ = std::move(rows.data);
        } else {
            indices_.insert(indices_.end(), rows.indices.begin(), rows.indices.end());
            data_.insert(data_.end(), rows.data.begin(), rows.data.end());
        }
        for (const std::int64_t count : rows.counts) {
            cursor += count;
            indptr_.push_back(cursor);
            ++next_row_;
        }
        rows = Rows();
    }
    [[nodiscard]] hicx::CsrMatrix finish(std::string dtype) {
        skip_to(n_);
        return hicx::CsrMatrix(n_, n_, std::move(indptr_), std::move(indices_),
                               std::move(data_), std::move(dtype));
    }

  private:
    std::int64_t n_;
    std::int64_t next_row_ = 0;
    std::vector<std::int64_t> indptr_;
    std::vector<std::int32_t> indices_;
    std::vector<double> data_;
};

// hicPCA builds its own MatrixFileHandler for the intermediate matrices, from
// the *output* file name, so an h5 input written to a .cool name really does
// produce a cooler. That is different from hiCMatrix.save, which reuses the
// input handler.
void write_intermediate(const std::string& path, const hicx::MatrixData& source,
                        hicx::CsrMatrix matrix) {
    hicx::MatrixData out;
    out.matrix = std::move(matrix);
    out.cut_intervals = source.cut_intervals;
    out.nan_bins = source.nan_bins;
    out.correction_factors = source.correction_factors;
    out.distance_counts = source.distance_counts;
    out.correction_factors_are_column = source.correction_factors_are_column;
    if (ends_with(path, ".h5")) {
        hicx::H5SaveOptions options;
        options.symmetric = true;
        hicx::write_hicexplorer_h5(path, out, options);
        return;
    }
    hicx::CoolSaveOptions options;
    options.symmetric = true;
    options.apply_correction = false;
    hicx::write_cool(path, out, options);
}

// ---------------------------------------------------------------------------
// --extraTrack

// hicexplorer.readBed on a gene track: only the chromosome, start and end of
// each interval are used (hicPCA.py:153-173). gzip is detected from the magic
// bytes, as utilities.opener does.
struct BedInterval {
    std::string chrom;
    std::int64_t start = 0;
    std::int64_t end = 0;
};

std::vector<BedInterval> read_bed(const std::string& path) {
    std::vector<BedInterval> intervals;
    std::string content;
    {
        std::FILE* handle = std::fopen(path.c_str(), "rb");
        if (handle == nullptr) {
            throw std::runtime_error("cannot open extraTrack file " + path);
        }
        unsigned char magic[2] = {0, 0};
        const std::size_t read = std::fread(magic, 1, 2, handle);
        std::fclose(handle);
        if (read == 2 && magic[0] == 0x1f && magic[1] == 0x8b) {
            gzFile gz = gzopen(path.c_str(), "rb");
            if (gz == nullptr) {
                throw std::runtime_error("cannot open extraTrack file " + path);
            }
            char buffer[65536];
            int got = 0;
            while ((got = gzread(gz, buffer, sizeof(buffer))) > 0) {
                content.append(buffer, static_cast<std::size_t>(got));
            }
            gzclose(gz);
        } else {
            std::FILE* plain = std::fopen(path.c_str(), "rb");
            char buffer[65536];
            std::size_t got = 0;
            while ((got = std::fread(buffer, 1, sizeof(buffer), plain)) > 0) {
                content.append(buffer, got);
            }
            std::fclose(plain);
        }
    }

    std::size_t position = 0;
    while (position < content.size()) {
        std::size_t newline = content.find('\n', position);
        if (newline == std::string::npos) {
            newline = content.size();
        }
        std::string line = content.substr(position, newline - position);
        position = newline + 1;
        while (!line.empty() && (line.back() == '\r' || line.back() == ' ')) {
            line.pop_back();
        }
        if (line.empty() || line[0] == '#' || line.rfind("track", 0) == 0 ||
            line.rfind("browser", 0) == 0) {
            continue;
        }
        std::vector<std::string> fields;
        std::size_t start = 0;
        while (true) {
            const std::size_t tab = line.find('\t', start);
            if (tab == std::string::npos) {
                fields.push_back(line.substr(start));
                break;
            }
            fields.push_back(line.substr(start, tab - start));
            start = tab + 1;
        }
        if (fields.size() < 3) {
            continue;
        }
        BedInterval interval;
        interval.chrom = fields[0];
        try {
            interval.start = std::stoll(fields[1]);
            interval.end = std::stoll(fields[2]);
        } catch (const std::exception&) {
            continue;
        }
        intervals.push_back(std::move(interval));
    }
    return intervals;
}

// scipy.stats.pearsonr's correlation coefficient. Only its sign is used
// (hicPCA.py:192), so the classic centred form is enough.
double pearson_correlation(const std::vector<double>& a, const std::vector<double>& b) {
    const std::size_t n = a.size();
    if (n == 0 || n != b.size()) {
        return 0.0;
    }
    double mean_a = 0.0;
    double mean_b = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        mean_a += a[i];
        mean_b += b[i];
    }
    mean_a /= static_cast<double>(n);
    mean_b /= static_cast<double>(n);
    double numerator = 0.0;
    double sum_a = 0.0;
    double sum_b = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double da = a[i] - mean_a;
        const double db = b[i] - mean_b;
        numerator += da * db;
        sum_a += da * da;
        sum_b += db * db;
    }
    const double denominator = std::sqrt(sum_a) * std::sqrt(sum_b);
    if (denominator == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return numerator / denominator;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);

    // hicPCA.py:242-248.
    if (args.which_eigenvectors.size() != args.output_file_names.size()) {
        std::fprintf(stderr,
                     "ERROR:hicexplorer.hicPCA:Number of output file names and "
                     "number of eigenvectors does not match. Please provide the "
                     "name of each file.\n");
        return 1;
    }
    std::vector<int> which;
    which.reserve(args.which_eigenvectors.size());
    for (const std::string& token : args.which_eigenvectors) {
        try {
            std::size_t consumed = 0;
            const int value = std::stoi(token, &consumed);
            if (consumed != token.size()) {
                throw std::invalid_argument("trailing characters");
            }
            which.push_back(value);
        } catch (const std::exception&) {
            // int(' ') in hicPCA.py:316 raises ValueError and the run ends.
            std::fprintf(stderr,
                         "hicPCA: --whichEigenvectors: invalid literal for int(): "
                         "'%s'\n",
                         token.c_str());
            return 1;
        }
    }

    try {
        // Before anything calls into LAPACK. See transform_ops.hpp: the
        // eigenvectors dgeev returns depend on the BLAS thread count, so the
        // tool would otherwise not be reproducible.
        hicx::pin_blas_to_one_thread();

        hicx::ToolMatrix hic = hicx::ToolMatrix::load(args.matrix);

        if (args.ignore_masked_bins) {
            if (!args.chromosomes.empty() && !hic.data().nan_bins.empty()) {
                // Pinned reference defect, reproduced and not fixed.
                // hicPCA.py:253-259 masks the NaN bins, replaces the bin
                // table with enlarge_bins and only then calls
                // keepOnlyTheseChr, whose first action (HiCMatrix.py:616) is
                // restoreMaskedBins. That undoes both steps: the masked bins
                // come back as empty rows and columns, the bin table comes
                // back from orig_cut_intervals, which maskBins saved before
                // the enlargement, and nan_bins becomes the masked set, which
                // keepOnlyTheseChr then narrows to the kept chromosomes. So
                // --ignoreMaskedBins is silently a no-op whenever
                // --chromosomes is given. The net effect of that round trip
                // is mask_and_restore_bins: masked entries dropped, shape and
                // bin table unchanged, float64 matrix, NaN correction factors
                // at masked bins. Measured on small_test_matrix_50kb_res.h5
                // chrX, which has four NaN bins: the Python writes all 449
                // bins at their original 50 kb boundaries.
                // test_pca_ignore_masked_bins_is_a_no_op_with_chromosomes pins
                // it.
                const std::vector<std::int64_t> masked = hic.data().nan_bins;
                hicx::mask_and_restore_bins(hic.data(), masked);
            } else {
                // Without --chromosomes nothing restores the mask: maskBins
                // removes the bins for good and setCutIntervals installs the
                // enlarged bins. With no NaN bins maskBins returns early,
                // restoreMaskedBins has nothing to undo and the enlargement
                // persists on either path, which is this branch too.
                hicx::delete_bins(hic.data(), hic.data().nan_bins);
                hicx::enlarge_bins(hic.data().cut_intervals);
            }
            hic.refresh_boundaries();
        }
        if (!args.chromosomes.empty()) {
            hicx::keep_only_chromosomes(hic.data(), args.chromosomes);
            hic.refresh_boundaries();
        }

        const std::vector<std::pair<std::string, hicx::BinRange>>& boundaries =
            hic.boundaries();
        const std::int64_t chromosome_count = static_cast<std::int64_t>(boundaries.size());
        std::int64_t length_chromosome = 0;
        for (const std::pair<std::string, hicx::BinRange>& entry : boundaries) {
            length_chromosome += entry.second.last - entry.second.first;
        }

        const std::int64_t nbins = hic.matrix().rows();
        // The row-wise builders require the chromosome blocks to arrive in
        // increasing bin order, which is the order a bin table taken from a
        // file has. A matrix whose chromosomes interleave would need the
        // rows sorted afterwards; say so rather than write a wrong file.
        std::int64_t previous_end = 0;
        for (const std::pair<std::string, hicx::BinRange>& entry : boundaries) {
            if (entry.second.first < previous_end) {
                throw std::runtime_error(
                    "the chromosomes of this matrix are not in bin order, which the "
                    "intermediate matrix writers require");
            }
            previous_end = entry.second.last;
        }

        std::optional<RowwiseCsrBuilder> pearson_builder;
        std::optional<RowwiseCsrBuilder> obsexp_builder;
        if (args.pearson_matrix.has_value()) {
            pearson_builder.emplace(nbins);
        }
        if (args.obsexp_matrix.has_value()) {
            obsexp_builder.emplace(nbins);
        }

        // The bedgraph and bigwig payload, one row per bin in chromosome order.
        std::vector<std::string> chrom_list;
        std::vector<std::int64_t> start_list;
        std::vector<std::int64_t> end_list;
        // components[k] is the k-th requested eigenvector, over all bins.
        std::vector<std::vector<double>> components(which.size());
        // Whether the bin's row carries all requested components. numpy's
        // hstack silently drops an out-of-range column, which shortens every
        // row of that chromosome and makes hicPCA skip it when writing.
        std::vector<bool> complete;

        const hicx::EigenSolver solver =
            args.compat_v4 ? hicx::EigenSolver::Dsyevr : hicx::EigenSolver::Dgeev;

        // --extraTrack, opened once.
        bigWigFile_t* histone_track = nullptr;
        const bool extra_track_is_bigwig =
            args.extra_track.has_value() && (ends_with(*args.extra_track, ".bw") ||
                                             ends_with(*args.extra_track, ".bigwig"));
        if (extra_track_is_bigwig) {
            if (bwInit(1 << 17) != 0) {
                throw std::runtime_error("libBigWig initialisation failed");
            }
            histone_track = bwOpen(const_cast<char*>(args.extra_track->c_str()), nullptr,
                                   "r");
            if (histone_track == nullptr) {
                throw std::runtime_error("cannot open extraTrack " + *args.extra_track);
            }
        }

        // One chromosome's share of the work, computed independently of every
        // other chromosome and consumed in chromosome order further down.
        struct ChromosomeResult {
            RowwiseCsrBuilder::Rows obsexp;
            RowwiseCsrBuilder::Rows pearson;
            hicx::EigenResult eigen;
        };
        const std::size_t chromosome_total = boundaries.size();
        std::vector<ChromosomeResult> results(chromosome_total);
        const bool want_obsexp = obsexp_builder.has_value();
        const bool want_pearson = pearson_builder.has_value();

        const auto compute = [&](std::size_t c) {
            const std::int64_t first = boundaries[c].second.first;
            const std::int64_t n = boundaries[c].second.last - first;
            ChromosomeResult& result = results[c];

            std::vector<std::int64_t> order(static_cast<std::size_t>(n));
            for (std::int64_t i = 0; i < n; ++i) {
                order[static_cast<std::size_t>(i)] = first + i;
            }
            hicx::CsrMatrix block = hicx::select_bins(hic.matrix(), order);
            // hicPCA densifies the block, so it is one of the tools cpp/PLAN.md
            // 4.4 rule 2 exempts from upper-triangle-only storage. The sparse
            // form is still what the covariance is computed from; only the
            // symmetric completion is materialised.
            block.materialize_full();

            if (args.method == "lieberman") {
                hicx::obs_exp_lieberman_in_place(block, length_chromosome,
                                                 chromosome_count);
            } else {
                hicx::obs_exp_non_zero_in_place(block, args.ligation_factor);
            }

            if (want_obsexp) {
                // lil_matrix(dense) keeps only the non-zeros, which is what
                // hicPCA.py:293 assigns into the accumulator.
                for (std::int64_t i = 0; i < n; ++i) {
                    const std::size_t row_start = result.obsexp.data.size();
                    const std::int64_t begin = block.indptr()[static_cast<std::size_t>(i)];
                    const std::int64_t stop =
                        block.indptr()[static_cast<std::size_t>(i) + 1];
                    for (std::int64_t k = begin; k < stop; ++k) {
                        const double value = block.data()[static_cast<std::size_t>(k)];
                        if (value != 0.0) {
                            result.obsexp.push(
                                first + block.indices()[static_cast<std::size_t>(k)], value);
                        }
                    }
                    result.obsexp.end_row(row_start);
                }
            }

            if (args.lanczos) {
                if (want_pearson) {
                    // np.corrcoef of the obs/exp block one row at a time out of
                    // the sparse matrix (transform_ops' DenseCorrelationRows).
                    // A one bin chromosome has no correlation and stores no
                    // entry, which is what the dense path writes for it.
                    std::vector<double> row(static_cast<std::size_t>(n));
                    std::optional<hicx::DenseCorrelationRows> rows;
                    if (n > 1) {
                        rows.emplace(block, std::vector<hicx::BinRange>{hicx::BinRange{0, n}},
                                     hicx::DenseCorrelationRows::Kind::Pearson, args.threads);
                    }
                    for (std::int64_t i = 0; i < n; ++i) {
                        const std::size_t row_start = result.pearson.data.size();
                        if (rows.has_value()) {
                            rows->fill_row(i, row.data());
                            for (std::int64_t j = 0; j < n; ++j) {
                                if (row[static_cast<std::size_t>(j)] != 0.0) {
                                    result.pearson.push(first + j,
                                                        row[static_cast<std::size_t>(j)]);
                                }
                            }
                        }
                        result.pearson.end_row(row_start);
                    }
                }
                result.eigen = hicx::covariance_eigenvectors_lanczos(block, which, args.threads);
                return;
            }

            // Consumes the sparse block. Bit-identical to np.cov, which is
            // what decides the column order dgeev returns; see
            // numpy_covariance.
            hicx::DenseSymmetric covariance = numpy_covariance(block, args.threads);

            if (want_pearson) {
                // np.corrcoef of the same obs/exp block, which is this
                // covariance rescaled by the square roots of its diagonal.
                // Read row by row so that the covariance survives for the
                // eigensolver and no second dense block is allocated.
                const std::vector<double> scaling = hicx::pearson_scaling(covariance);
                std::vector<double> row(static_cast<std::size_t>(n));
                for (std::int64_t i = 0; i < n; ++i) {
                    const std::size_t row_start = result.pearson.data.size();
                    hicx::pearson_row(covariance, scaling, i, row.data());
                    for (std::int64_t j = 0; j < n; ++j) {
                        if (row[static_cast<std::size_t>(j)] != 0.0) {
                            result.pearson.push(first + j, row[static_cast<std::size_t>(j)]);
                        }
                    }
                    result.pearson.end_row(row_start);
                }
            }

            clean_covariance_like_hicpca(covariance, args.threads);
            result.eigen = solver == hicx::EigenSolver::Dgeev
                               ? hicx::dgeev_selected_eigenvectors(covariance, which)
                               : hicx::leading_eigenvectors(covariance, which, solver);
        };

        // Chromosomes in parallel for the dense solver. Each computation is
        // serial LAPACK on its own block, so the results do not depend on
        // which run together; the scheduler only bounds memory, by admitting a
        // chromosome while the dense blocks in flight stay within the largest
        // chromosome's, which is what a serial run holds at its peak. The
        // Krylov solver threads its sparse products instead and takes the
        // chromosomes one at a time.
        {
            std::vector<std::size_t> by_size(chromosome_total);
            std::iota(by_size.begin(), by_size.end(), std::size_t{0});
            const auto bins_of = [&](std::size_t c) {
                return static_cast<double>(boundaries[c].second.last - boundaries[c].second.first);
            };
            std::stable_sort(by_size.begin(), by_size.end(), [&](std::size_t a, std::size_t b) {
                return bins_of(a) > bins_of(b);
            });
            double budget = 0.0;
            for (std::size_t c = 0; c < chromosome_total; ++c) {
                budget = std::max(budget, bins_of(c) * bins_of(c));
            }
            const int workers =
                args.lanczos ? 1
                             : std::max(1, std::min(args.threads, static_cast<int>(chromosome_total)));
            std::mutex mutex;
            std::condition_variable changed;
            std::vector<char> started(chromosome_total, 0);
            double in_flight = 0.0;
            std::exception_ptr failure;
            const auto worker = [&]() {
                std::unique_lock<std::mutex> lock(mutex);
                while (!failure) {
                    bool any_left = false;
                    std::size_t pick = chromosome_total;
                    for (const std::size_t c : by_size) {
                        if (started[c] != 0) {
                            continue;
                        }
                        any_left = true;
                        const double cost = bins_of(c) * bins_of(c);
                        if (in_flight == 0.0 || in_flight + cost <= budget) {
                            pick = c;
                            break;
                        }
                    }
                    if (!any_left) {
                        return;
                    }
                    if (pick == chromosome_total) {
                        changed.wait(lock);
                        continue;
                    }
                    started[pick] = 1;
                    const double cost = bins_of(pick) * bins_of(pick);
                    in_flight += cost;
                    lock.unlock();
                    std::exception_ptr error;
                    try {
                        compute(pick);
                    } catch (...) {
                        error = std::current_exception();
                    }
                    lock.lock();
                    in_flight -= cost;
                    if (error && !failure) {
                        failure = error;
                    }
                    changed.notify_all();
                }
            };
            if (workers == 1) {
                worker();
            } else {
                std::vector<std::thread> pool;
                pool.reserve(static_cast<std::size_t>(workers));
                for (int w = 0; w < workers; ++w) {
                    pool.emplace_back(worker);
                }
                for (std::thread& thread : pool) {
                    thread.join();
                }
            }
            if (failure) {
                std::rethrow_exception(failure);
            }
        }

        for (std::size_t c = 0; c < chromosome_total; ++c) {
            const std::string& chrname = boundaries[c].first;
            const std::int64_t first = boundaries[c].second.first;
            const std::int64_t n = boundaries[c].second.last - first;
            ChromosomeResult& result = results[c];
            if (want_obsexp) {
                obsexp_builder->append(first, std::move(result.obsexp));
            }
            if (want_pearson) {
                pearson_builder->append(first, std::move(result.pearson));
            }
            hicx::EigenResult& eigen = result.eigen;

            bool all_present = true;
            for (const std::vector<double>& vector : eigen.vectors) {
                if (vector.empty()) {
                    all_present = false;
                }
            }

            // correlateEigenvectorWithHistonMarkTrack, applied per chromosome
            // and per eigenvector before the values are appended. The Python
            // guards the whole block with `if chromosome in
            // bwTrack.chroms().keys()` (hicPCA.py:213), so a chromosome the
            // track does not carry is left alone rather than queried bin by
            // bin.
            bool track_has_chromosome = false;
            if (histone_track != nullptr && histone_track->cl != nullptr) {
                for (std::int64_t key = 0; key < histone_track->cl->nKeys; ++key) {
                    if (chrname == histone_track->cl->chrom[key]) {
                        track_has_chromosome = true;
                        break;
                    }
                }
            }
            if (histone_track != nullptr && track_has_chromosome && all_present) {
                for (std::vector<double>& vector : eigen.vectors) {
                    double positive_sum = 0.0;
                    double negative_sum = 0.0;
                    std::size_t positive_count = 0;
                    std::size_t negative_count = 0;
                    for (std::int64_t i = 0; i < n; ++i) {
                        const double value = vector[static_cast<std::size_t>(i)];
                        if (value == 0.0) {
                            continue;
                        }
                        const hicx::CutInterval& bin =
                            hic.data().cut_intervals[static_cast<std::size_t>(first + i)];
                        double* statistics =
                            bwStats(histone_track, const_cast<char*>(chrname.c_str()),
                                    static_cast<uint32_t>(bin.start),
                                    static_cast<uint32_t>(bin.end), 1, mean);
                        double statistic = 0.0;
                        if (statistics != nullptr) {
                            if (!std::isnan(statistics[0])) {
                                statistic = statistics[0];
                            }
                            free(statistics);
                        }
                        if (value > 0.0) {
                            ++positive_count;
                            positive_sum += statistic;
                        } else {
                            ++negative_count;
                            negative_sum += statistic;
                        }
                    }
                    const double positive_mean =
                        positive_sum != 0.0
                            ? positive_sum / static_cast<double>(positive_count)
                            : 0.0;
                    const double negative_mean =
                        negative_sum != 0.0
                            ? negative_sum / static_cast<double>(negative_count)
                            : 0.0;
                    const bool flip =
                        args.histone_mark_type == "active"
                            ? (positive_mean < negative_mean && negative_mean != 0.0 &&
                               positive_mean != 0.0)
                            : (positive_mean > negative_mean && negative_mean != 0.0 &&
                               positive_mean != 0.0);
                    if (flip) {
                        for (double& value : vector) {
                            if (value != 0.0) {
                                value = -value;
                            }
                        }
                    }
                }
            }

            for (std::int64_t i = 0; i < n; ++i) {
                const hicx::CutInterval& bin =
                    hic.data().cut_intervals[static_cast<std::size_t>(first + i)];
                chrom_list.push_back(bin.chrom);
                start_list.push_back(bin.start);
                end_list.push_back(bin.end);
                complete.push_back(all_present);
                for (std::size_t k = 0; k < which.size(); ++k) {
                    components[k].push_back(
                        eigen.vectors[k].empty()
                            ? 0.0
                            : eigen.vectors[k][static_cast<std::size_t>(i)]);
                }
            }
        }

        if (histone_track != nullptr) {
            bwClose(histone_track);
            bwCleanup();
            histone_track = nullptr;
        }

        if (pearson_builder.has_value()) {
            write_intermediate(*args.pearson_matrix, hic.data(),
                               pearson_builder->finish("float64"));
        }
        if (obsexp_builder.has_value()) {
            write_intermediate(*args.obsexp_matrix, hic.data(),
                               obsexp_builder->finish("float64"));
        }

        // correlateEigenvectorWithGeneTrack: a bed extraTrack flips a whole
        // chromosome when its eigenvector correlates negatively with the gene
        // density. Applied after both intermediate matrices are written, as in
        // hicPCA.py:354.
        if (args.extra_track.has_value() && !extra_track_is_bigwig) {
            const hicx::BinTable bins(hic.data().cut_intervals);
            std::vector<double> gene_occurrence(
                static_cast<std::size_t>(hic.data().cut_intervals.size()), 0.0);
            const std::vector<std::pair<std::string, std::int64_t>> sizes =
                bins.chromosome_sizes();
            for (const BedInterval& interval : read_bed(*args.extra_track)) {
                const std::optional<hicx::BinRange> range =
                    bins.chrom_bin_range(interval.chrom);
                if (!range.has_value()) {
                    continue;  // chromosome not in the matrix
                }
                std::int64_t chromosome_size = 0;
                for (const std::pair<std::string, std::int64_t>& size : sizes) {
                    if (size.first == interval.chrom) {
                        chromosome_size = size.second;
                        break;
                    }
                }
                if (interval.start > chromosome_size) {
                    continue;  // hicPCA.py:157, warns and skips
                }
                const std::optional<std::pair<std::int64_t, std::int64_t>> bin_id =
                    bins.region_bin_range(interval.chrom, interval.start, interval.end);
                if (!bin_id.has_value()) {
                    continue;  // hicPCA.py:166
                }
                // The Python counts the *end* bin, bin_id[1].
                gene_occurrence[static_cast<std::size_t>(bin_id->second)] += 1.0;
            }
            for (const std::pair<std::string, hicx::BinRange>& entry : boundaries) {
                const std::int64_t first = entry.second.first;
                const std::int64_t last = entry.second.last;
                std::vector<double> density(
                    gene_occurrence.begin() + static_cast<std::ptrdiff_t>(first),
                    gene_occurrence.begin() + static_cast<std::ptrdiff_t>(last));
                for (std::size_t k = 0; k < which.size(); ++k) {
                    std::vector<double> slice(
                        components[k].begin() + static_cast<std::ptrdiff_t>(first),
                        components[k].begin() + static_cast<std::ptrdiff_t>(last));
                    const double correlation = pearson_correlation(slice, density);
                    if (correlation < 0.0) {
                        for (std::int64_t i = first; i < last; ++i) {
                            components[k][static_cast<std::size_t>(i)] =
                                -components[k][static_cast<std::size_t>(i)];
                        }
                    }
                }
            }
        }

        if (args.format == "bedgraph") {
            for (std::size_t k = 0; k < args.output_file_names.size(); ++k) {
                std::FILE* out = std::fopen(args.output_file_names[k].c_str(), "w");
                if (out == nullptr) {
                    throw std::runtime_error("cannot write " + args.output_file_names[k]);
                }
                for (std::size_t i = 0; i < chrom_list.size(); ++i) {
                    if (!complete[i]) {
                        continue;
                    }
                    std::fprintf(out, "%s\t%lld\t%lld\t%s\n", chrom_list[i].c_str(),
                                 static_cast<long long>(start_list[i]),
                                 static_cast<long long>(end_list[i]),
                                 format_twelve_decimals(components[k][i]).c_str());
                }
                std::fclose(out);
            }
        } else {
            if (chrom_list.empty()) {
                throw std::runtime_error("no bins to write");
            }
            // hicPCA.py:374-381 builds the header from the last bin end of
            // every run of equal chromosome names.
            std::vector<std::string> header_names;
            std::vector<std::uint32_t> header_lengths;
            std::string previous = chrom_list[0];
            for (std::size_t i = 0; i < chrom_list.size(); ++i) {
                if (previous != chrom_list[i]) {
                    header_names.push_back(previous);
                    header_lengths.push_back(static_cast<std::uint32_t>(end_list[i - 1]));
                }
                previous = chrom_list[i];
            }
            header_names.push_back(chrom_list.back());
            header_lengths.push_back(static_cast<std::uint32_t>(end_list.back()));

            if (bwInit(1 << 17) != 0) {
                throw std::runtime_error("libBigWig initialisation failed");
            }
            for (std::size_t k = 0; k < args.output_file_names.size(); ++k) {
                bigWigFile_t* out =
                    bwOpen(const_cast<char*>(args.output_file_names[k].c_str()), nullptr,
                           "w");
                if (out == nullptr) {
                    throw std::runtime_error("cannot write " + args.output_file_names[k]);
                }
                // pyBigWig's addHeader default is maxZooms=10.
                if (bwCreateHdr(out, 10) != 0) {
                    throw std::runtime_error("bigwig header creation failed");
                }
                std::vector<char*> names;
                names.reserve(header_names.size());
                for (std::string& name : header_names) {
                    names.push_back(name.data());
                }
                out->cl = bwCreateChromList(names.data(), header_lengths.data(),
                                            static_cast<std::int64_t>(names.size()));
                if (out->cl == nullptr) {
                    throw std::runtime_error("bigwig chromosome list creation failed");
                }
                if (bwWriteHdr(out) != 0) {
                    throw std::runtime_error("bigwig header write failed");
                }

                // One bwAddIntervals call per run of equal chromosome names,
                // which is what pyBigWig's addEntries does internally.
                std::size_t i = 0;
                while (i < chrom_list.size()) {
                    if (!complete[i]) {
                        ++i;
                        continue;
                    }
                    const std::string& chrom = chrom_list[i];
                    std::vector<char*> block_chroms;
                    std::vector<std::uint32_t> block_starts;
                    std::vector<std::uint32_t> block_ends;
                    std::vector<float> block_values;
                    while (i < chrom_list.size() && chrom_list[i] == chrom) {
                        if (complete[i]) {
                            block_chroms.push_back(
                                const_cast<char*>(chrom_list[i].c_str()));
                            block_starts.push_back(
                                static_cast<std::uint32_t>(start_list[i]));
                            block_ends.push_back(static_cast<std::uint32_t>(end_list[i]));
                            block_values.push_back(static_cast<float>(components[k][i]));
                        }
                        ++i;
                    }
                    if (block_chroms.empty()) {
                        continue;
                    }
                    const int status = bwAddIntervals(
                        out, block_chroms.data(), block_starts.data(), block_ends.data(),
                        block_values.data(),
                        static_cast<std::uint32_t>(block_chroms.size()));
                    if (status != 0) {
                        throw std::runtime_error("bigwig interval write failed");
                    }
                }
                bwClose(out);
            }
            bwCleanup();
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicPCA: %s\n", error.what());
        return 1;
    }

    hicx::report_resource_usage("hicPCA");
    return 0;
}
