// See hicx/dgeev_selected.hpp for why and what. The transcriptions follow
// lapack-netlib/SRC/{dgeev,dtrevc3}.f of OpenBLAS 0.3.28
// statement by statement, 1-based indices kept, so that they can be checked
// against the Fortran line by line.

#include "hicx/dgeev_selected.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>


extern "C" {
// Fortran character arguments carry a hidden trailing length.
void dgeev_(const char*, const char*, const int*, double*, const int*, double*, double*, double*,
            const int*, double*, const int*, double*, const int*, int*, std::size_t, std::size_t);
double dlamch_(const char*, std::size_t);
double dlange_(const char*, const int*, const int*, const double*, const int*, double*,
               std::size_t);
void dlascl_(const char*, const int*, const int*, const double*, const double*, const int*,
             const int*, double*, const int*, int*, std::size_t);
void dgebal_(const char*, const int*, double*, const int*, int*, int*, double*, int*, std::size_t);
void dgehrd_(const int*, const int*, const int*, double*, const int*, double*, double*,
             const int*, int*);
void dorghr_(const int*, const int*, const int*, double*, const int*, const double*, double*,
             const int*, int*);
void dhseqr_(const char*, const char*, const int*, const int*, const int*, double*, const int*,
             double*, double*, double*, const int*, double*, const int*, int*, std::size_t,
             std::size_t);
void dgebak_(const char*, const char*, const int*, const int*, const int*, const double*,
             const int*, double*, const int*, int*, std::size_t, std::size_t);
void dlaln2_(const int*, const int*, const int*, const double*, const double*, const double*,
             const int*, const double*, const double*, const double*, const int*, const double*,
             const double*, double*, const int*, double*, double*, int*);
void dlacpy_(const char*, const int*, const int*, const double*, const int*, double*, const int*,
             std::size_t);
double dnrm2_(const int*, const double*, const int*);
double dlapy2_(const double*, const double*);
int idamax_(const int*, const double*, const int*);
void dlartg_(const double*, const double*, double*, double*, double*);
void drot_(const int*, double*, const int*, double*, const int*, const double*, const double*);
void dscal_(const int*, const double*, double*, const int*);
void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
void dgemm_(const char*, const char*, const int*, const int*, const int*, const double*,
            const double*, const int*, const double*, const int*, const double*, double*,
            const int*, std::size_t, std::size_t);
}

namespace hicx {
namespace {

// Pointer to element (i, j), 1-based, of a column major array.
inline double* at(double* p, int ld, int i, int j) {
    return p + static_cast<std::size_t>(i - 1) +
           static_cast<std::size_t>(j - 1) * static_cast<std::size_t>(ld);
}
inline const double* at(const double* p, int ld, int i, int j) {
    return p + static_cast<std::size_t>(i - 1) +
           static_cast<std::size_t>(j - 1) * static_cast<std::size_t>(ld);
}

// dtrevc3's block size for this workspace; 1 means the unblocked path, which
// is not replayed.
constexpr int kTrevcNbMin = 8;
constexpr int kTrevcNbMax = 128;
int trevc3_block_size(int n, int lwork) {
    if (lwork >= n + 2 * n * kTrevcNbMin) {
        return std::min((lwork - n) / (2 * n), kTrevcNbMax);
    }
    return 1;
}

struct TrevcMember {
    int ki;  // KI of dtrevc3's loop, 1-based
    int ip;  // 0 for a real eigenvalue, -1 for the complex pair (KI-1, KI)
    int iv;  // work column; for a pair the imaginary part, the real part at iv-1
};
struct TrevcBlock {
    std::vector<TrevcMember> members;
    int iv_flush = 0;
    int ki2 = 0;
};

// dtrevc3's right eigenvector loop, bookkeeping only: which KI goes to which
// work column and where each block is flushed through dgemm.
std::vector<TrevcBlock> trevc3_blocks(int n, const double* t, int nb) {
    std::vector<TrevcBlock> blocks;
    TrevcBlock current;
    int iv = nb > 2 ? nb : 2;
    int ip = 0;
    for (int ki = n; ki >= 1; --ki) {
        if (ip == -1) {
            ip = 1;
            continue;  // GO TO 140 skips the flush test as well
        }
        if (ki == 1 || *at(t, n, ki, ki - 1) == 0.0) {
            ip = 0;
        } else {
            ip = -1;
        }
        current.members.push_back({ki, ip, iv});
        if (ip != 0) {
            iv = iv - 1;
        }
        const int ki2 = ip == 0 ? ki : ki - 1;
        if (iv <= 2 || ki2 == 1) {
            current.iv_flush = iv;
            current.ki2 = ki2;
            blocks.push_back(std::move(current));
            current = TrevcBlock();
            iv = nb;
        } else {
            iv = iv - 1;
        }
    }
    return blocks;
}

// DTREVC3(SIDE='R', HOWMNY='B') with NB > 1 for the blocks holding one of
// `wanted` (1-based columns). `vr` holds the Schur vectors and is only read.
// Returns column -> the vector DTREVC3 leaves in that column of VR.
std::map<int, std::vector<double>> trevc3_selected(int n, const double* t, const double* vr,
                                                   int nb, const std::set<int>& wanted) {
    const std::size_t nn = static_cast<std::size_t>(n);
    const char safe_minimum = 'S';
    const char precision = 'P';
    const double unfl = dlamch_(&safe_minimum, 1);
    const double ulp = dlamch_(&precision, 1);
    const double smlnum = unfl * (static_cast<double>(n) / ulp);
    const double bignum = (1.0 - ulp) / smlnum;

    const std::vector<TrevcBlock> blocks = trevc3_blocks(n, t, nb);
    std::vector<const TrevcBlock*> needed;
    int top = 0;
    for (const TrevcBlock& block : blocks) {
        const int first = block.ki2;
        const int last = block.ki2 + nb - block.iv_flush;
        bool hit = false;
        for (const int column : wanted) {
            hit = hit || (column >= first && column <= last);
        }
        if (hit) {
            needed.push_back(&block);
            top = std::max(top, last);
        }
    }

    std::vector<double> work(nn * static_cast<std::size_t>(1 + 2 * nb), 0.0);
    const auto W = [&](int i, int column) -> double& {
        return work[static_cast<std::size_t>(i - 1) + static_cast<std::size_t>(column) * nn];
    };
    const auto T = [&](int i, int j) -> const double& { return *at(t, n, i, j); };
    // WORK(J): 1-norm of the strictly upper part of column J of T.
    W(1, 0) = 0.0;
    for (int j = 2; j <= top; ++j) {
        W(j, 0) = 0.0;
        for (int i = 1; i <= j - 1; ++i) {
            W(j, 0) = W(j, 0) + std::fabs(T(i, j));
        }
    }

    const int ione = 1;
    const int itwo = 2;
    const int logical_false = 0;
    const double one = 1.0;
    const double zero = 0.0;
    std::vector<int> iscomplex(static_cast<std::size_t>(kTrevcNbMax) + 1, 0);
    double x[4] = {0.0, 0.0, 0.0, 0.0};  // X(2,2), column major
    const auto X = [&](int i, int j) -> double& { return x[(i - 1) + (j - 1) * 2]; };
    std::map<int, std::vector<double>> result;

    for (const TrevcBlock* block : needed) {
        for (const TrevcMember& member : block->members) {
            const int ki = member.ki;
            const int ip = member.ip;
            const int iv = member.iv;
            const double wr = T(ki, ki);
            double wi = 0.0;
            if (ip != 0) {
                wi = std::sqrt(std::fabs(T(ki, ki - 1))) * std::sqrt(std::fabs(T(ki - 1, ki)));
            }
            const double smin = std::max(ulp * (std::fabs(wr) + std::fabs(wi)), smlnum);
            double scale = 0.0;
            double xnorm = 0.0;
            int ierr = 0;
            if (ip == 0) {
                // Real right eigenvector.
                W(ki, iv) = 1.0;
                for (int k = 1; k <= ki - 1; ++k) {
                    W(k, iv) = -T(k, ki);
                }
                int jnxt = ki - 1;
                for (int j = ki - 1; j >= 1; --j) {
                    if (j > jnxt) {
                        continue;
                    }
                    int j1 = j;
                    const int j2 = j;
                    jnxt = j - 1;
                    if (j > 1 && T(j, j - 1) != 0.0) {
                        j1 = j - 1;
                        jnxt = j - 2;
                    }
                    if (j1 == j2) {
                        dlaln2_(&logical_false, &ione, &ione, &smin, &one, &T(j, j), &n, &one, &one,
                                &W(j, iv), &n, &wr, &zero, x, &itwo, &scale, &xnorm, &ierr);
                        if (xnorm > 1.0 && W(j, 0) > bignum / xnorm) {
                            X(1, 1) = X(1, 1) / xnorm;
                            scale = scale / xnorm;
                        }
                        if (scale != 1.0) {
                            dscal_(&ki, &scale, &W(1, iv), &ione);
                        }
                        W(j, iv) = X(1, 1);
                        const int jm1 = j - 1;
                        const double factor = -X(1, 1);
                        daxpy_(&jm1, &factor, &T(1, j), &ione, &W(1, iv), &ione);
                    } else {
                        dlaln2_(&logical_false, &itwo, &ione, &smin, &one, &T(j - 1, j - 1), &n,
                                &one, &one, &W(j - 1, iv), &n, &wr, &zero, x, &itwo, &scale,
                                &xnorm, &ierr);
                        if (xnorm > 1.0) {
                            const double beta = std::max(W(j - 1, 0), W(j, 0));
                            if (beta > bignum / xnorm) {
                                X(1, 1) = X(1, 1) / xnorm;
                                X(2, 1) = X(2, 1) / xnorm;
                                scale = scale / xnorm;
                            }
                        }
                        if (scale != 1.0) {
                            dscal_(&ki, &scale, &W(1, iv), &ione);
                        }
                        W(j - 1, iv) = X(1, 1);
                        W(j, iv) = X(2, 1);
                        const int jm2 = j - 2;
                        double factor = -X(1, 1);
                        daxpy_(&jm2, &factor, &T(1, j - 1), &ione, &W(1, iv), &ione);
                        factor = -X(2, 1);
                        daxpy_(&jm2, &factor, &T(1, j), &ione, &W(1, iv), &ione);
                    }
                }
                for (int k = ki + 1; k <= n; ++k) {
                    W(k, iv) = 0.0;
                }
                iscomplex[static_cast<std::size_t>(iv)] = ip;
            } else {
                // Complex right eigenvector, real part in iv-1, imaginary in iv.
                if (std::fabs(T(ki - 1, ki)) >= std::fabs(T(ki, ki - 1))) {
                    W(ki - 1, iv - 1) = 1.0;
                    W(ki, iv) = wi / T(ki - 1, ki);
                } else {
                    W(ki - 1, iv - 1) = -wi / T(ki, ki - 1);
                    W(ki, iv) = 1.0;
                }
                W(ki, iv - 1) = 0.0;
                W(ki - 1, iv) = 0.0;
                for (int k = 1; k <= ki - 2; ++k) {
                    W(k, iv - 1) = -W(ki - 1, iv - 1) * T(k, ki - 1);
                    W(k, iv) = -W(ki, iv) * T(k, ki);
                }
                int jnxt = ki - 2;
                for (int j = ki - 2; j >= 1; --j) {
                    if (j > jnxt) {
                        continue;
                    }
                    int j1 = j;
                    const int j2 = j;
                    jnxt = j - 1;
                    if (j > 1 && T(j, j - 1) != 0.0) {
                        j1 = j - 1;
                        jnxt = j - 2;
                    }
                    if (j1 == j2) {
                        dlaln2_(&logical_false, &ione, &itwo, &smin, &one, &T(j, j), &n, &one, &one,
                                &W(j, iv - 1), &n, &wr, &wi, x, &itwo, &scale, &xnorm, &ierr);
                        if (xnorm > 1.0 && W(j, 0) > bignum / xnorm) {
                            X(1, 1) = X(1, 1) / xnorm;
                            X(1, 2) = X(1, 2) / xnorm;
                            scale = scale / xnorm;
                        }
                        if (scale != 1.0) {
                            dscal_(&ki, &scale, &W(1, iv - 1), &ione);
                            dscal_(&ki, &scale, &W(1, iv), &ione);
                        }
                        W(j, iv - 1) = X(1, 1);
                        W(j, iv) = X(1, 2);
                        const int jm1 = j - 1;
                        double factor = -X(1, 1);
                        daxpy_(&jm1, &factor, &T(1, j), &ione, &W(1, iv - 1), &ione);
                        factor = -X(1, 2);
                        daxpy_(&jm1, &factor, &T(1, j), &ione, &W(1, iv), &ione);
                    } else {
                        dlaln2_(&logical_false, &itwo, &itwo, &smin, &one, &T(j - 1, j - 1), &n,
                                &one, &one, &W(j - 1, iv - 1), &n, &wr, &wi, x, &itwo, &scale,
                                &xnorm, &ierr);
                        if (xnorm > 1.0) {
                            const double beta = std::max(W(j - 1, 0), W(j, 0));
                            if (beta > bignum / xnorm) {
                                const double rec = 1.0 / xnorm;
                                X(1, 1) = X(1, 1) * rec;
                                X(1, 2) = X(1, 2) * rec;
                                X(2, 1) = X(2, 1) * rec;
                                X(2, 2) = X(2, 2) * rec;
                                scale = scale * rec;
                            }
                        }
                        if (scale != 1.0) {
                            dscal_(&ki, &scale, &W(1, iv - 1), &ione);
                            dscal_(&ki, &scale, &W(1, iv), &ione);
                        }
                        W(j - 1, iv - 1) = X(1, 1);
                        W(j, iv - 1) = X(2, 1);
                        W(j - 1, iv) = X(1, 2);
                        W(j, iv) = X(2, 2);
                        const int jm2 = j - 2;
                        double factor = -X(1, 1);
                        daxpy_(&jm2, &factor, &T(1, j - 1), &ione, &W(1, iv - 1), &ione);
                        factor = -X(2, 1);
                        daxpy_(&jm2, &factor, &T(1, j), &ione, &W(1, iv - 1), &ione);
                        factor = -X(1, 2);
                        daxpy_(&jm2, &factor, &T(1, j - 1), &ione, &W(1, iv), &ione);
                        factor = -X(2, 2);
                        daxpy_(&jm2, &factor, &T(1, j), &ione, &W(1, iv), &ione);
                    }
                }
                for (int k = ki + 1; k <= n; ++k) {
                    W(k, iv - 1) = 0.0;
                    W(k, iv) = 0.0;
                }
                iscomplex[static_cast<std::size_t>(iv - 1)] = -ip;
                iscomplex[static_cast<std::size_t>(iv)] = ip;
            }
        }

        // The flush: one dgemm with dtrevc3's dimensions, then normalisation.
        const int ivf = block->iv_flush;
        const int ki2 = block->ki2;
        const int columns = nb - ivf + 1;
        const int inner = ki2 + nb - ivf;
        const char cN = 'N';
        dgemm_(&cN, &cN, &n, &columns, &inner, &one, vr, &n,
               &work[static_cast<std::size_t>(ivf) * nn], &n, &zero,
               &work[static_cast<std::size_t>(nb + ivf) * nn], &n, 1, 1);
        double remax = 0.0;
        for (int k = ivf; k <= nb; ++k) {
            const double* column = &work[static_cast<std::size_t>(nb + k) * nn];
            if (iscomplex[static_cast<std::size_t>(k)] == 0) {
                const int ii = idamax_(&n, column, &ione);
                remax = 1.0 / std::fabs(W(ii, nb + k));
            } else if (iscomplex[static_cast<std::size_t>(k)] == 1) {
                double emax = 0.0;
                for (int ii = 1; ii <= n; ++ii) {
                    emax = std::max(emax, std::fabs(W(ii, nb + k)) + std::fabs(W(ii, nb + k + 1)));
                }
                remax = 1.0 / emax;
            }
            // ISCOMPLEX(K) = -1 keeps the REMAX of K-1, as the Fortran does.
            dscal_(&n, &remax, &work[static_cast<std::size_t>(nb + k) * nn], &ione);
        }
        for (int k = ivf; k <= nb; ++k) {
            const double* source = &work[static_cast<std::size_t>(nb + k) * nn];
            result[ki2 + (k - ivf)].assign(source, source + n);
        }
    }
    return result;
}

EigenResult to_result(int n, const std::vector<int>& which, const std::vector<double>& wr,
                      const std::vector<double>& wi,
                      const std::function<const double*(int)>& column,
                      std::vector<double>* values_imaginary) {
    EigenResult result;
    result.values.assign(which.size(), 0.0);
    result.vectors.assign(which.size(), {});
    if (values_imaginary != nullptr) {
        values_imaginary->assign(which.size(), 0.0);
    }
    for (std::size_t entry = 0; entry < which.size(); ++entry) {
        const int index = which[entry] - 1;
        if (index < 0 || index >= n) {
            continue;
        }
        // A conjugate pair occupies columns (j, j+1) as (real, imaginary);
        // scipy builds both complex vectors from them and hicPCA keeps the real
        // part, which is column j for both members.
        int source = index;
        if (wi[static_cast<std::size_t>(index)] < 0.0 && index > 0) {
            source = index - 1;
        }
        const double* values = column(source);
        result.values[entry] = wr[static_cast<std::size_t>(index)];
        result.vectors[entry].assign(values, values + n);
        if (values_imaginary != nullptr) {
            (*values_imaginary)[entry] = wi[static_cast<std::size_t>(index)];
        }
    }
    return result;
}

EigenResult selected_impl(int n, double* a, const std::vector<int>& which,
                          std::vector<double>* values_imaginary, const std::function<void()>& release) {
    if (n == 0) {
        release();
        EigenResult empty;
        empty.values.assign(which.size(), 0.0);
        empty.vectors.assign(which.size(), {});
        if (values_imaginary != nullptr) {
            values_imaginary->assign(which.size(), 0.0);
        }
        return empty;
    }
    const char cV = 'V';
    const char cN = 'N';
    const char cS = 'S';
    const char cB = 'B';
    const char cR = 'R';
    const char cL = 'L';
    const char cG = 'G';
    const char cM = 'M';
    const char cP = 'P';
    std::vector<double> wr(static_cast<std::size_t>(n));
    std::vector<double> wi(static_cast<std::size_t>(n));
    int info = 0;
    int lwork = -1;
    double query = 0.0;
    double dummy = 0.0;
    dgeev_(&cN, &cV, &n, a, &n, wr.data(), wi.data(), &dummy, &n, &dummy, &n, &query, &lwork,
           &info, 1, 1);
    if (info != 0) {
        throw std::runtime_error("dgeev workspace query failed");
    }
    lwork = static_cast<int>(query);
    std::vector<double> work(static_cast<std::size_t>(std::max(lwork, 1)));
    std::vector<double> vr(static_cast<std::size_t>(n) * static_cast<std::size_t>(n));

    // dgeev hands dtrevc3 WORK(IWRK) with IWRK = ITAU = N + 1.
    const int trevc_nb = trevc3_block_size(n, lwork - n);
    if (trevc_nb <= 1) {
        // The unblocked back-transform is not replayed: run dgeev itself.
        dgeev_(&cN, &cV, &n, a, &n, wr.data(), wi.data(), &dummy, &n, vr.data(), &n, work.data(),
               &lwork, &info, 1, 1);
        if (info != 0) {
            throw std::runtime_error("dgeev failed with info " + std::to_string(info));
        }
        release();
        return to_result(
            n, which, wr, wi,
            [&](int column) {
                return vr.data() + static_cast<std::size_t>(column) * static_cast<std::size_t>(n);
            },
            values_imaginary);
    }

    // dgeev, statement by statement, up to dhseqr.
    const double eps = dlamch_(&cP, 1);
    double smlnum = dlamch_(&cS, 1);
    double bignum = 1.0 / smlnum;
    smlnum = std::sqrt(smlnum) / eps;
    bignum = 1.0 / smlnum;
    double dum = 0.0;
    const double anrm = dlange_(&cM, &n, &n, a, &n, &dum, 1);
    bool scalea = false;
    double cscale = 0.0;
    if (anrm > 0.0 && anrm < smlnum) {
        scalea = true;
        cscale = smlnum;
    } else if (anrm > bignum) {
        scalea = true;
        cscale = bignum;
    }
    const int izero = 0;
    int ierr = 0;
    if (scalea) {
        dlascl_(&cG, &izero, &izero, &anrm, &cscale, &n, &n, a, &n, &ierr, 1);
    }
    int ilo = 0;
    int ihi = 0;
    dgebal_(&cB, &n, a, &n, &ilo, &ihi, work.data(), &ierr, 1);
    const int itau = n;  // 0-based WORK(ITAU)
    int iwrk = itau + n;
    int remaining = lwork - iwrk;
    dgehrd_(&n, &ilo, &ihi, a, &n, work.data() + itau, work.data() + iwrk, &remaining, &ierr);
    dlacpy_(&cL, &n, &n, a, &n, vr.data(), &n, 1);
    dorghr_(&n, &ilo, &ihi, vr.data(), &n, work.data() + itau, work.data() + iwrk, &remaining,
            &ierr);
    iwrk = itau;
    remaining = lwork - iwrk;
    dhseqr_(&cS, &cV, &n, &ilo, &ihi, a, &n, wr.data(), wi.data(), vr.data(), &n,
            work.data() + iwrk, &remaining, &info, 1, 1);
    if (info != 0) {
        throw std::runtime_error("dgeev failed with info " + std::to_string(info));
    }

    // The columns dgeev's result is read from, with both columns of every
    // complex pair, which the normalisation below mixes.
    std::set<int> columns;  // 1-based
    for (const int requested : which) {
        const int index = requested - 1;
        if (index < 0 || index >= n) {
            continue;
        }
        int first = index;
        if (wi[static_cast<std::size_t>(index)] < 0.0 && index > 0) {
            first = index - 1;
        }
        columns.insert(first + 1);
        if (wi[static_cast<std::size_t>(first)] != 0.0 && first + 1 < n) {
            columns.insert(first + 2);
        }
    }
    std::map<int, std::vector<double>> transformed =
        trevc3_selected(n, a, vr.data(), trevc_nb, columns);
    vr.clear();
    vr.shrink_to_fit();
    release();

    const int count = static_cast<int>(columns.size());
    const std::vector<int> order(columns.begin(), columns.end());
    std::vector<double> v(static_cast<std::size_t>(n) * static_cast<std::size_t>(std::max(count, 1)));
    for (int q = 0; q < count; ++q) {
        const std::vector<double>& source = transformed[order[static_cast<std::size_t>(q)]];
        std::copy(source.begin(), source.end(),
                  v.begin() + static_cast<std::ptrdiff_t>(q) * static_cast<std::ptrdiff_t>(n));
    }
    // DGEBAK('B', 'R') scales and swaps whole rows, so on a subset of the
    // columns it does what it does to those columns of the full matrix.
    if (count > 0) {
        dgebak_(&cB, &cR, &n, &ilo, &ihi, work.data(), &count, v.data(), &n, &ierr, 1, 1);
    }
    // dgeev's normalisation (the loop at label 40).
    const int ione = 1;
    std::vector<double> squares(static_cast<std::size_t>(n));
    for (int q = 0; q < count; ++q) {
        const int column = order[static_cast<std::size_t>(q)];
        double* values = v.data() + static_cast<std::size_t>(q) * static_cast<std::size_t>(n);
        const double imaginary = wi[static_cast<std::size_t>(column - 1)];
        if (imaginary == 0.0) {
            const double scl = 1.0 / dnrm2_(&n, values, &ione);
            dscal_(&n, &scl, values, &ione);
        } else if (imaginary > 0.0) {
            // The partner column + 1 is the next entry of the sorted set.
            double* partner = values + n;
            const double norm_real = dnrm2_(&n, values, &ione);
            const double norm_imaginary = dnrm2_(&n, partner, &ione);
            const double scl = 1.0 / dlapy2_(&norm_real, &norm_imaginary);
            dscal_(&n, &scl, values, &ione);
            dscal_(&n, &scl, partner, &ione);
            for (int k = 0; k < n; ++k) {
                squares[static_cast<std::size_t>(k)] = values[k] * values[k] + partner[k] * partner[k];
            }
            const int k = idamax_(&n, squares.data(), &ione) - 1;
            double cs = 0.0;
            double sn = 0.0;
            double r = 0.0;
            dlartg_(&values[k], &partner[k], &cs, &sn, &r);
            drot_(&n, values, &ione, partner, &ione, &cs, &sn);
            partner[k] = 0.0;
        }
    }
    if (scalea) {
        const int rows = n - info;
        const int lda = std::max(n - info, 1);
        const int one_column = 1;
        dlascl_(&cG, &izero, &izero, &cscale, &anrm, &rows, &one_column, wr.data() + info, &lda,
                &ierr, 1);
        dlascl_(&cG, &izero, &izero, &cscale, &anrm, &rows, &one_column, wi.data() + info, &lda,
                &ierr, 1);
    }
    return to_result(
        n, which, wr, wi,
        [&](int column) {
            const auto position =
                std::find(order.begin(), order.end(), column + 1) - order.begin();
            return v.data() + static_cast<std::size_t>(position) * static_cast<std::size_t>(n);
        },
        values_imaginary);
}

}  // namespace

EigenResult dgeev_selected_eigenvectors(DenseSymmetric& matrix, const std::vector<int>& which) {
    const int n = static_cast<int>(matrix.size());
    return selected_impl(n, matrix.data(), which, nullptr,
                         [&matrix] { matrix.release(); });
}

EigenResult dgeev_selected_general(int n, double* matrix, const std::vector<int>& which,
                                   std::vector<double>* values_imaginary) {
    return selected_impl(n, matrix, which, values_imaginary, [] {});
}

}  // namespace hicx
