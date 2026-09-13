// L-BFGS-B 3.0 as scipy 1.14.1 runs it. See hicx/lbfgsb_scipy.hpp for what is
// translated from where and why this exists beside the projected BFGS of
// lbfgsb.cpp.
//
// Translation conventions. The Fortran is followed statement by statement, and
// the arrays keep their Fortran shape: every two dimensional array is stored
// column major in a flat vector and addressed through a 1-based accessor with
// the leading dimension the Fortran declares, so `ws(i, j)` reads as
// `ws(i, j)` here too. The places where the Fortran passes a pointer into the
// middle of an array (`wn(1, js)`, `wa(2*m+1)`, `wn(col+1, col+1)`) pass the
// same offset. The reverse communication of setulb is replaced by direct calls
// of the objective at the points where mainlb returns 'FG_START', 'FG_LNSRCH'
// and 'NEW_X', which changes the control flow's shape and nothing it computes.
//
// Original code: L-BFGS-B (version 3.0), Ciyou Zhu, Richard Byrd, Peihuang Lu,
// Jorge Nocedal, Jose Luis Morales. Copyright of the Fortran by its authors,
// distributed under the 3-clause BSD licence as part of scipy.

#include "hicx/lbfgsb_scipy.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "hicx/numpy_compat.hpp"

namespace hicx::stats {

namespace {

// ---------------------------------------------------------------------------
// Fortran array helpers

// Element (i, j) of a column major array with leading dimension ld, 1-based.
inline double& el(double* a, int ld, int i, int j) {
    return a[(i - 1) + static_cast<std::ptrdiff_t>(j - 1) * ld];
}
// Pointer to element (i, j), for passing a sub-array the way Fortran does.
inline double* ptr(double* a, int ld, int i, int j) {
    return a + (i - 1) + static_cast<std::ptrdiff_t>(j - 1) * ld;
}

// The reference BLAS, whose loops reduce strictly left to right.
double ddot(int n, const double* x, const double* y) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum = sum + x[i] * y[i];
    }
    return sum;
}
void daxpy(int n, double a, const double* x, double* y) {
    if (n <= 0 || a == 0.0) {
        return;
    }
    for (int i = 0; i < n; ++i) {
        y[i] = y[i] + a * x[i];
    }
}
void dscal(int n, double a, double* x) {
    for (int i = 0; i < n; ++i) {
        x[i] = a * x[i];
    }
}
void dcopy(int n, const double* x, double* y) {
    for (int i = 0; i < n; ++i) {
        y[i] = x[i];
    }
}

// ---------------------------------------------------------------------------
// LINPACK

// Cholesky factorisation of a symmetric positive definite matrix; the upper
// triangle of `a` is replaced by R with A = R'R.
void dpofa(double* a, int lda, int n, int& info) {
    for (int j = 1; j <= n; ++j) {
        info = j;
        double s = 0.0;
        const int jm1 = j - 1;
        for (int k = 1; k <= jm1; ++k) {
            double t = el(a, lda, k, j) - ddot(k - 1, ptr(a, lda, 1, k), ptr(a, lda, 1, j));
            t = t / el(a, lda, k, k);
            el(a, lda, k, j) = t;
            s = s + t * t;
        }
        s = el(a, lda, j, j) - s;
        if (s <= 0.0) {
            return;
        }
        el(a, lda, j, j) = std::sqrt(s);
    }
    info = 0;
}

// Solves T x = b or T' x = b for triangular T. job 00: T lower, T x = b;
// 01: T upper, T x = b; 10: T lower, T' x = b; 11: T upper, T' x = b.
void dtrsl(double* t, int ldt, int n, double* b, int job, int& info) {
    for (info = 1; info <= n; ++info) {
        if (el(t, ldt, info, info) == 0.0) {
            return;
        }
    }
    info = 0;
    int kase = 1;
    if (job % 10 != 0) {
        kase = 2;
    }
    if ((job % 100) / 10 != 0) {
        kase = kase + 2;
    }
    const auto bb = [b](int i) -> double& { return b[i - 1]; };
    switch (kase) {
        case 1:
            bb(1) = bb(1) / el(t, ldt, 1, 1);
            for (int j = 2; j <= n; ++j) {
                const double temp = -bb(j - 1);
                daxpy(n - j + 1, temp, ptr(t, ldt, j, j - 1), &bb(j));
                bb(j) = bb(j) / el(t, ldt, j, j);
            }
            break;
        case 2:
            bb(n) = bb(n) / el(t, ldt, n, n);
            for (int jj = 2; jj <= n; ++jj) {
                const int j = n - jj + 1;
                const double temp = -bb(j + 1);
                daxpy(j, temp, ptr(t, ldt, 1, j + 1), &bb(1));
                bb(j) = bb(j) / el(t, ldt, j, j);
            }
            break;
        case 3:
            bb(n) = bb(n) / el(t, ldt, n, n);
            for (int jj = 2; jj <= n; ++jj) {
                const int j = n - jj + 1;
                bb(j) = bb(j) - ddot(jj - 1, ptr(t, ldt, j + 1, j), &bb(j + 1));
                bb(j) = bb(j) / el(t, ldt, j, j);
            }
            break;
        default:
            bb(1) = bb(1) / el(t, ldt, 1, 1);
            for (int j = 2; j <= n; ++j) {
                bb(j) = bb(j) - ddot(j - 1, ptr(t, ldt, 1, j), &bb(1));
                bb(j) = bb(j) / el(t, ldt, j, j);
            }
            break;
    }
}

// ---------------------------------------------------------------------------
// L-BFGS-B subroutines

void active(int n, const double* l, const double* u, const int* nbd, double* x,
            int* iwhere, bool& prjctd, bool& cnstnd, bool& boxed) {
    prjctd = false;
    cnstnd = false;
    boxed = true;
    for (int i = 0; i < n; ++i) {
        if (nbd[i] > 0) {
            if (nbd[i] <= 2 && x[i] <= l[i]) {
                if (x[i] < l[i]) {
                    prjctd = true;
                    x[i] = l[i];
                }
            } else if (nbd[i] >= 2 && x[i] >= u[i]) {
                if (x[i] > u[i]) {
                    prjctd = true;
                    x[i] = u[i];
                }
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        if (nbd[i] != 2) {
            boxed = false;
        }
        if (nbd[i] == 0) {
            iwhere[i] = -1;
        } else {
            cnstnd = true;
            if (nbd[i] == 2 && u[i] - l[i] <= 0.0) {
                iwhere[i] = 3;
            } else {
                iwhere[i] = 0;
            }
        }
    }
}

// Product of the 2m x 2m middle matrix with v, into p.
void bmv(int m, double* sy, double* wt, int col, const double* v, double* p, int& info) {
    if (col == 0) {
        return;
    }
    const auto vv = [v](int i) { return v[i - 1]; };
    const auto pp = [p](int i) -> double& { return p[i - 1]; };
    pp(col + 1) = vv(col + 1);
    for (int i = 2; i <= col; ++i) {
        const int i2 = col + i;
        double sum = 0.0;
        for (int k = 1; k <= i - 1; ++k) {
            sum = sum + el(sy, m, i, k) * vv(k) / el(sy, m, k, k);
        }
        pp(i2) = vv(i2) + sum;
    }
    dtrsl(wt, m, col, &pp(col + 1), 11, info);
    if (info != 0) {
        return;
    }
    for (int i = 1; i <= col; ++i) {
        pp(i) = vv(i) / std::sqrt(el(sy, m, i, i));
    }
    dtrsl(wt, m, col, &pp(col + 1), 1, info);
    if (info != 0) {
        return;
    }
    for (int i = 1; i <= col; ++i) {
        pp(i) = -pp(i) / std::sqrt(el(sy, m, i, i));
    }
    for (int i = 1; i <= col; ++i) {
        double sum = 0.0;
        for (int k = i + 1; k <= col; ++k) {
            sum = sum + el(sy, m, k, i) * pp(col + k) / el(sy, m, i, i);
        }
        pp(i) = pp(i) + sum;
    }
}

void hpsolb(int n, double* t, int* iorder, int iheap) {
    const auto tt = [t](int i) -> double& { return t[i - 1]; };
    const auto io = [iorder](int i) -> int& { return iorder[i - 1]; };
    if (iheap == 0) {
        for (int k = 2; k <= n; ++k) {
            const double ddum = tt(k);
            const int indxin = io(k);
            int i = k;
            while (i > 1) {
                const int j = i / 2;
                if (ddum < tt(j)) {
                    tt(i) = tt(j);
                    io(i) = io(j);
                    i = j;
                } else {
                    break;
                }
            }
            tt(i) = ddum;
            io(i) = indxin;
        }
    }
    if (n > 1) {
        int i = 1;
        const double out = tt(1);
        const int indxou = io(1);
        const double ddum = tt(n);
        const int indxin = io(n);
        while (true) {
            int j = i + i;
            if (j <= n - 1) {
                if (tt(j + 1) < tt(j)) {
                    j = j + 1;
                }
                if (tt(j) < ddum) {
                    tt(i) = tt(j);
                    io(i) = io(j);
                    i = j;
                    continue;
                }
            }
            break;
        }
        tt(i) = ddum;
        io(i) = indxin;
        tt(n) = out;
        io(n) = indxou;
    }
}

// The generalised Cauchy point.
void cauchy(int n, double* x, double* l, double* u, int* nbd, double* g, int* iorder,
            int* iwhere, double* t, double* d, double* xcp, int m, double* wy, double* ws,
            double* sy, double* wt, double theta, int col, int head, double* p, double* c,
            double* wbp, double* v, int& nseg, double sbgnrm, int& info, double epsmch) {
    const auto X = [x](int i) -> double& { return x[i - 1]; };
    const auto L = [l](int i) { return l[i - 1]; };
    const auto U = [u](int i) { return u[i - 1]; };
    const auto NBD = [nbd](int i) { return nbd[i - 1]; };
    const auto G = [g](int i) { return g[i - 1]; };
    const auto IORDER = [iorder](int i) -> int& { return iorder[i - 1]; };
    const auto IWHERE = [iwhere](int i) -> int& { return iwhere[i - 1]; };
    const auto T = [t](int i) -> double& { return t[i - 1]; };
    const auto D = [d](int i) -> double& { return d[i - 1]; };
    const auto XCP = [xcp](int i) -> double& { return xcp[i - 1]; };
    const auto P = [p](int i) -> double& { return p[i - 1]; };
    const auto C = [c](int i) -> double& { return c[i - 1]; };
    const auto WBP = [wbp](int i) -> double& { return wbp[i - 1]; };

    if (sbgnrm <= 0.0) {
        dcopy(n, x, xcp);
        return;
    }
    bool bnded = true;
    int nfree = n + 1;
    int nbreak = 0;
    int ibkmin = 0;
    double tl = 0.0;
    double tu = 0.0;
    double bkmin = 0.0;
    const int col2 = 2 * col;
    double f1 = 0.0;

    for (int i = 1; i <= col2; ++i) {
        P(i) = 0.0;
    }

    for (int i = 1; i <= n; ++i) {
        const double neggi = -G(i);
        if (IWHERE(i) != 3 && IWHERE(i) != -1) {
            if (NBD(i) <= 2) {
                tl = X(i) - L(i);
            }
            if (NBD(i) >= 2) {
                tu = U(i) - X(i);
            }
            const bool xlower = NBD(i) <= 2 && tl <= 0.0;
            const bool xupper = NBD(i) >= 2 && tu <= 0.0;
            IWHERE(i) = 0;
            if (xlower) {
                if (neggi <= 0.0) {
                    IWHERE(i) = 1;
                }
            } else if (xupper) {
                if (neggi >= 0.0) {
                    IWHERE(i) = 2;
                }
            } else {
                if (std::abs(neggi) <= 0.0) {
                    IWHERE(i) = -3;
                }
            }
        }
        int pointr = head;
        if (IWHERE(i) != 0 && IWHERE(i) != -1) {
            D(i) = 0.0;
        } else {
            D(i) = neggi;
            f1 = f1 - neggi * neggi;
            for (int j = 1; j <= col; ++j) {
                P(j) = P(j) + el(wy, n, i, pointr) * neggi;
                P(col + j) = P(col + j) + el(ws, n, i, pointr) * neggi;
                pointr = pointr % m + 1;
            }
            if (NBD(i) <= 2 && NBD(i) != 0 && neggi < 0.0) {
                nbreak = nbreak + 1;
                IORDER(nbreak) = i;
                T(nbreak) = tl / (-neggi);
                if (nbreak == 1 || T(nbreak) < bkmin) {
                    bkmin = T(nbreak);
                    ibkmin = nbreak;
                }
            } else if (NBD(i) >= 2 && neggi > 0.0) {
                nbreak = nbreak + 1;
                IORDER(nbreak) = i;
                T(nbreak) = tu / neggi;
                if (nbreak == 1 || T(nbreak) < bkmin) {
                    bkmin = T(nbreak);
                    ibkmin = nbreak;
                }
            } else {
                nfree = nfree - 1;
                IORDER(nfree) = i;
                if (std::abs(neggi) > 0.0) {
                    bnded = false;
                }
            }
        }
    }

    if (theta != 1.0) {
        dscal(col, theta, &P(col + 1));
    }

    dcopy(n, x, xcp);

    if (nbreak == 0 && nfree == n + 1) {
        return;
    }

    for (int j = 1; j <= col2; ++j) {
        C(j) = 0.0;
    }

    double f2 = -theta * f1;
    const double f2_org = f2;
    if (col > 0) {
        bmv(m, sy, wt, col, p, v, info);
        if (info != 0) {
            return;
        }
        f2 = f2 - ddot(col2, v, p);
    }
    double dtm = -f1 / f2;
    double tsum = 0.0;
    nseg = 1;

    bool skip_to_888 = nbreak == 0;
    bool reached_999 = false;
    if (!skip_to_888) {
        int nleft = nbreak;
        int iter = 1;
        double tj = 0.0;
        while (true) {
            const double tj0 = tj;
            int ibp = 0;
            if (iter == 1) {
                tj = bkmin;
                ibp = IORDER(ibkmin);
            } else {
                if (iter == 2) {
                    if (ibkmin != nbreak) {
                        T(ibkmin) = T(nbreak);
                        IORDER(ibkmin) = IORDER(nbreak);
                    }
                }
                hpsolb(nleft, t, iorder, iter - 2);
                tj = T(nleft);
                ibp = IORDER(nleft);
            }

            const double dt = tj - tj0;

            if (dtm < dt) {
                break;  // goto 888
            }

            tsum = tsum + dt;
            nleft = nleft - 1;
            iter = iter + 1;
            const double dibp = D(ibp);
            D(ibp) = 0.0;
            double zibp = 0.0;
            if (dibp > 0.0) {
                zibp = U(ibp) - X(ibp);
                XCP(ibp) = U(ibp);
                IWHERE(ibp) = 2;
            } else {
                zibp = L(ibp) - X(ibp);
                XCP(ibp) = L(ibp);
                IWHERE(ibp) = 1;
            }
            if (nleft == 0 && nbreak == n) {
                dtm = dt;
                reached_999 = true;
                break;
            }

            nseg = nseg + 1;
            const double dibp2 = dibp * dibp;

            f1 = f1 + dt * f2 + dibp2 - theta * dibp * zibp;
            f2 = f2 - theta * dibp2;

            if (col > 0) {
                daxpy(col2, dt, p, c);
                int pointr = head;
                for (int j = 1; j <= col; ++j) {
                    WBP(j) = el(wy, n, ibp, pointr);
                    WBP(col + j) = theta * el(ws, n, ibp, pointr);
                    pointr = pointr % m + 1;
                }
                bmv(m, sy, wt, col, wbp, v, info);
                if (info != 0) {
                    return;
                }
                const double wmc = ddot(col2, c, v);
                const double wmp = ddot(col2, p, v);
                const double wmw = ddot(col2, wbp, v);
                daxpy(col2, -dibp, wbp, p);
                f1 = f1 + dibp * wmc;
                f2 = f2 + 2.0 * dibp * wmp - dibp2 * wmw;
            }
            f2 = std::max(epsmch * f2_org, f2);
            if (nleft > 0) {
                dtm = -f1 / f2;
                continue;  // goto 777
            } else if (bnded) {
                f1 = 0.0;
                f2 = 0.0;
                dtm = 0.0;
            } else {
                dtm = -f1 / f2;
            }
            break;  // falls through to 888
        }
    }

    if (!reached_999) {
        // 888
        if (dtm <= 0.0) {
            dtm = 0.0;
        }
        tsum = tsum + dtm;
        daxpy(n, tsum, d, xcp);
    }
    // 999
    if (col > 0) {
        daxpy(col2, dtm, p, c);
    }
}

void cmprlb(int n, int m, double* x, double* g, double* ws, double* wy, double* sy,
            double* wt, double* z, double* r, double* wa, int* index, double theta, int col,
            int head, int nfree, bool cnstnd, int& info) {
    if (!cnstnd && col > 0) {
        for (int i = 0; i < n; ++i) {
            r[i] = -g[i];
        }
    } else {
        for (int i = 1; i <= nfree; ++i) {
            const int k = index[i - 1];
            r[i - 1] = -theta * (z[k - 1] - x[k - 1]) - g[k - 1];
        }
        bmv(m, sy, wt, col, &wa[2 * m], &wa[0], info);
        if (info != 0) {
            info = -8;
            return;
        }
        int pointr = head;
        for (int j = 1; j <= col; ++j) {
            const double a1 = wa[j - 1];
            const double a2 = theta * wa[col + j - 1];
            for (int i = 1; i <= nfree; ++i) {
                const int k = index[i - 1];
                r[i - 1] = r[i - 1] + el(wy, n, k, pointr) * a1 + el(ws, n, k, pointr) * a2;
            }
            pointr = pointr % m + 1;
        }
    }
}

void formk(int n, int nsub, int* ind, int nenter, int ileave, int* indx2, int iupdat,
           bool updatd, double* wn, double* wn1, int m, double* ws, double* wy, double* sy,
           double theta, int col, int head, int& info) {
    const int m2 = 2 * m;
    const auto IND = [ind](int i) { return ind[i - 1]; };
    const auto INDX2 = [indx2](int i) { return indx2[i - 1]; };
    int upcl = 0;

    if (updatd) {
        if (iupdat > m) {
            for (int jy = 1; jy <= m - 1; ++jy) {
                const int js = m + jy;
                dcopy(m - jy, ptr(wn1, m2, jy + 1, jy + 1), ptr(wn1, m2, jy, jy));
                dcopy(m - jy, ptr(wn1, m2, js + 1, js + 1), ptr(wn1, m2, js, js));
                dcopy(m - 1, ptr(wn1, m2, m + 2, jy + 1), ptr(wn1, m2, m + 1, jy));
            }
        }

        const int pbegin = 1;
        const int pend = nsub;
        const int dbegin = nsub + 1;
        const int dend = n;
        const int iy = col;
        const int is = m + col;
        int ipntr = head + col - 1;
        if (ipntr > m) {
            ipntr = ipntr - m;
        }
        int jpntr = head;
        for (int jy = 1; jy <= col; ++jy) {
            const int js = m + jy;
            double temp1 = 0.0;
            double temp2 = 0.0;
            double temp3 = 0.0;
            for (int k = pbegin; k <= pend; ++k) {
                const int k1 = IND(k);
                temp1 = temp1 + el(wy, n, k1, ipntr) * el(wy, n, k1, jpntr);
            }
            for (int k = dbegin; k <= dend; ++k) {
                const int k1 = IND(k);
                temp2 = temp2 + el(ws, n, k1, ipntr) * el(ws, n, k1, jpntr);
                temp3 = temp3 + el(ws, n, k1, ipntr) * el(wy, n, k1, jpntr);
            }
            el(wn1, m2, iy, jy) = temp1;
            el(wn1, m2, is, js) = temp2;
            el(wn1, m2, is, jy) = temp3;
            jpntr = jpntr % m + 1;
        }

        const int jy = col;
        jpntr = head + col - 1;
        if (jpntr > m) {
            jpntr = jpntr - m;
        }
        ipntr = head;
        for (int i = 1; i <= col; ++i) {
            const int is2 = m + i;
            double temp3 = 0.0;
            for (int k = pbegin; k <= pend; ++k) {
                const int k1 = IND(k);
                temp3 = temp3 + el(ws, n, k1, ipntr) * el(wy, n, k1, jpntr);
            }
            ipntr = ipntr % m + 1;
            el(wn1, m2, is2, jy) = temp3;
        }
        upcl = col - 1;
    } else {
        upcl = col;
    }

    int ipntr = head;
    for (int iy = 1; iy <= upcl; ++iy) {
        const int is = m + iy;
        int jpntr = head;
        for (int jy = 1; jy <= iy; ++jy) {
            const int js = m + jy;
            double temp1 = 0.0;
            double temp2 = 0.0;
            double temp3 = 0.0;
            double temp4 = 0.0;
            for (int k = 1; k <= nenter; ++k) {
                const int k1 = INDX2(k);
                temp1 = temp1 + el(wy, n, k1, ipntr) * el(wy, n, k1, jpntr);
                temp2 = temp2 + el(ws, n, k1, ipntr) * el(ws, n, k1, jpntr);
            }
            for (int k = ileave; k <= n; ++k) {
                const int k1 = INDX2(k);
                temp3 = temp3 + el(wy, n, k1, ipntr) * el(wy, n, k1, jpntr);
                temp4 = temp4 + el(ws, n, k1, ipntr) * el(ws, n, k1, jpntr);
            }
            el(wn1, m2, iy, jy) = el(wn1, m2, iy, jy) + temp1 - temp3;
            el(wn1, m2, is, js) = el(wn1, m2, is, js) - temp2 + temp4;
            jpntr = jpntr % m + 1;
        }
        ipntr = ipntr % m + 1;
    }

    ipntr = head;
    for (int is = m + 1; is <= m + upcl; ++is) {
        int jpntr = head;
        for (int jy = 1; jy <= upcl; ++jy) {
            double temp1 = 0.0;
            double temp3 = 0.0;
            for (int k = 1; k <= nenter; ++k) {
                const int k1 = INDX2(k);
                temp1 = temp1 + el(ws, n, k1, ipntr) * el(wy, n, k1, jpntr);
            }
            for (int k = ileave; k <= n; ++k) {
                const int k1 = INDX2(k);
                temp3 = temp3 + el(ws, n, k1, ipntr) * el(wy, n, k1, jpntr);
            }
            if (is <= jy + m) {
                el(wn1, m2, is, jy) = el(wn1, m2, is, jy) + temp1 - temp3;
            } else {
                el(wn1, m2, is, jy) = el(wn1, m2, is, jy) - temp1 + temp3;
            }
            jpntr = jpntr % m + 1;
        }
        ipntr = ipntr % m + 1;
    }

    for (int iy = 1; iy <= col; ++iy) {
        const int is = col + iy;
        const int is1 = m + iy;
        for (int jy = 1; jy <= iy; ++jy) {
            const int js = col + jy;
            const int js1 = m + jy;
            el(wn, m2, jy, iy) = el(wn1, m2, iy, jy) / theta;
            el(wn, m2, js, is) = el(wn1, m2, is1, js1) * theta;
        }
        for (int jy = 1; jy <= iy - 1; ++jy) {
            el(wn, m2, jy, is) = -el(wn1, m2, is1, jy);
        }
        for (int jy = iy; jy <= col; ++jy) {
            el(wn, m2, jy, is) = el(wn1, m2, is1, jy);
        }
        el(wn, m2, iy, iy) = el(wn, m2, iy, iy) + el(sy, m, iy, iy);
    }

    dpofa(wn, m2, col, info);
    if (info != 0) {
        info = -1;
        return;
    }
    const int col2 = 2 * col;
    for (int js = col + 1; js <= col2; ++js) {
        dtrsl(wn, m2, col, ptr(wn, m2, 1, js), 11, info);
    }
    for (int is = col + 1; is <= col2; ++is) {
        for (int js = is; js <= col2; ++js) {
            el(wn, m2, is, js) =
                el(wn, m2, is, js) + ddot(col, ptr(wn, m2, 1, is), ptr(wn, m2, 1, js));
        }
    }
    dpofa(ptr(wn, m2, col + 1, col + 1), m2, col, info);
    if (info != 0) {
        info = -2;
        return;
    }
}

void formt(int m, double* wt, double* sy, double* ss, int col, double theta, int& info) {
    for (int j = 1; j <= col; ++j) {
        el(wt, m, 1, j) = theta * el(ss, m, 1, j);
    }
    for (int i = 2; i <= col; ++i) {
        for (int j = i; j <= col; ++j) {
            const int k1 = std::min(i, j) - 1;
            double ddum = 0.0;
            for (int k = 1; k <= k1; ++k) {
                ddum = ddum + el(sy, m, i, k) * el(sy, m, j, k) / el(sy, m, k, k);
            }
            el(wt, m, i, j) = ddum + theta * el(ss, m, i, j);
        }
    }
    dpofa(wt, m, col, info);
    if (info != 0) {
        info = -3;
    }
}

void freev(int n, int& nfree, int* index, int& nenter, int& ileave, int* indx2, int* iwhere,
           bool& wrk, bool updatd, bool cnstnd, int iter) {
    nenter = 0;
    ileave = n + 1;
    if (iter > 0 && cnstnd) {
        for (int i = 1; i <= nfree; ++i) {
            const int k = index[i - 1];
            if (iwhere[k - 1] > 0) {
                ileave = ileave - 1;
                indx2[ileave - 1] = k;
            }
        }
        for (int i = 1 + nfree; i <= n; ++i) {
            const int k = index[i - 1];
            if (iwhere[k - 1] <= 0) {
                nenter = nenter + 1;
                indx2[nenter - 1] = k;
            }
        }
    }
    wrk = (ileave < n + 1) || (nenter > 0) || updatd;

    nfree = 0;
    int iact = n + 1;
    for (int i = 1; i <= n; ++i) {
        if (iwhere[i - 1] <= 0) {
            nfree = nfree + 1;
            index[nfree - 1] = i;
        } else {
            iact = iact - 1;
            index[iact - 1] = i;
        }
    }
}

void matupd(int n, int m, double* ws, double* wy, double* sy, double* ss, double* d, double* r,
            int& itail, int iupdat, int& col, int& head, double& theta, double rr, double dr,
            double stp, double dtd) {
    if (iupdat <= m) {
        col = iupdat;
        itail = (head + iupdat - 2) % m + 1;
    } else {
        itail = itail % m + 1;
        head = head % m + 1;
    }

    dcopy(n, d, ptr(ws, n, 1, itail));
    dcopy(n, r, ptr(wy, n, 1, itail));

    theta = rr / dr;

    if (iupdat > m) {
        for (int j = 1; j <= col - 1; ++j) {
            dcopy(j, ptr(ss, m, 2, j + 1), ptr(ss, m, 1, j));
            dcopy(col - j, ptr(sy, m, j + 1, j + 1), ptr(sy, m, j, j));
        }
    }
    int pointr = head;
    for (int j = 1; j <= col - 1; ++j) {
        el(sy, m, col, j) = ddot(n, d, ptr(wy, n, 1, pointr));
        el(ss, m, j, col) = ddot(n, ptr(ws, n, 1, pointr), d);
        pointr = pointr % m + 1;
    }
    if (stp == 1.0) {
        el(ss, m, col, col) = dtd;
    } else {
        el(ss, m, col, col) = stp * stp * dtd;
    }
    el(sy, m, col, col) = dr;
}

void projgr(int n, const double* l, const double* u, const int* nbd, const double* x,
            const double* g, double& sbgnrm) {
    sbgnrm = 0.0;
    for (int i = 0; i < n; ++i) {
        double gi = g[i];
        if (gi != gi) {
            sbgnrm = gi;
            return;
        }
        if (nbd[i] != 0) {
            if (gi < 0.0) {
                if (nbd[i] >= 2) {
                    gi = std::max(x[i] - u[i], gi);
                }
            } else {
                if (nbd[i] <= 2) {
                    gi = std::min(x[i] - l[i], gi);
                }
            }
        }
        sbgnrm = std::max(sbgnrm, std::abs(gi));
    }
}

void subsm(int n, int m, int nsub, int* ind, double* l, double* u, int* nbd, double* x,
           double* d, double* xp, double* ws, double* wy, double theta, double* xx, double* gg,
           int col, int head, int& iword, double* wv, double* wn, int& info) {
    if (nsub <= 0) {
        return;
    }
    const auto IND = [ind](int i) { return ind[i - 1]; };
    const auto D = [d](int i) -> double& { return d[i - 1]; };
    const auto X = [x](int i) -> double& { return x[i - 1]; };
    const auto L = [l](int i) { return l[i - 1]; };
    const auto U = [u](int i) { return u[i - 1]; };
    const auto NBD = [nbd](int i) { return nbd[i - 1]; };
    const auto WV = [wv](int i) -> double& { return wv[i - 1]; };

    int pointr = head;
    for (int i = 1; i <= col; ++i) {
        double temp1 = 0.0;
        double temp2 = 0.0;
        for (int j = 1; j <= nsub; ++j) {
            const int k = IND(j);
            temp1 = temp1 + el(wy, n, k, pointr) * D(j);
            temp2 = temp2 + el(ws, n, k, pointr) * D(j);
        }
        WV(i) = temp1;
        WV(col + i) = theta * temp2;
        pointr = pointr % m + 1;
    }

    const int m2 = 2 * m;
    const int col2 = 2 * col;
    dtrsl(wn, m2, col2, wv, 11, info);
    if (info != 0) {
        return;
    }
    for (int i = 1; i <= col; ++i) {
        WV(i) = -WV(i);
    }
    dtrsl(wn, m2, col2, wv, 1, info);
    if (info != 0) {
        return;
    }

    pointr = head;
    for (int jy = 1; jy <= col; ++jy) {
        const int js = col + jy;
        for (int i = 1; i <= nsub; ++i) {
            const int k = IND(i);
            D(i) = D(i) + el(wy, n, k, pointr) * WV(jy) / theta + el(ws, n, k, pointr) * WV(js);
        }
        pointr = pointr % m + 1;
    }

    dscal(nsub, 1.0 / theta, d);

    iword = 0;
    dcopy(n, x, xp);

    for (int i = 1; i <= nsub; ++i) {
        const int k = IND(i);
        const double dk = D(i);
        double xk = X(k);
        if (NBD(k) != 0) {
            if (NBD(k) == 1) {
                X(k) = std::max(L(k), xk + dk);
                if (X(k) == L(k)) {
                    iword = 1;
                }
            } else {
                if (NBD(k) == 2) {
                    xk = std::max(L(k), xk + dk);
                    X(k) = std::min(U(k), xk);
                    if (X(k) == L(k) || X(k) == U(k)) {
                        iword = 1;
                    }
                } else {
                    if (NBD(k) == 3) {
                        X(k) = std::min(U(k), xk + dk);
                        if (X(k) == U(k)) {
                            iword = 1;
                        }
                    }
                }
            }
        } else {
            X(k) = xk + dk;
        }
    }

    if (iword == 0) {
        return;
    }

    double dd_p = 0.0;
    for (int i = 0; i < n; ++i) {
        dd_p = dd_p + (x[i] - xx[i]) * gg[i];
    }
    if (dd_p > 0.0) {
        dcopy(n, xp, x);
    } else {
        return;
    }

    double alpha = 1.0;
    double temp1 = alpha;
    int ibd = 0;
    for (int i = 1; i <= nsub; ++i) {
        const int k = IND(i);
        const double dk = D(i);
        if (NBD(k) != 0) {
            if (dk < 0.0 && NBD(k) <= 2) {
                const double temp2 = L(k) - X(k);
                if (temp2 >= 0.0) {
                    temp1 = 0.0;
                } else if (dk * alpha < temp2) {
                    temp1 = temp2 / dk;
                }
            } else if (dk > 0.0 && NBD(k) >= 2) {
                const double temp2 = U(k) - X(k);
                if (temp2 <= 0.0) {
                    temp1 = 0.0;
                } else if (dk * alpha > temp2) {
                    temp1 = temp2 / dk;
                }
            }
            if (temp1 < alpha) {
                alpha = temp1;
                ibd = i;
            }
        }
    }

    if (alpha < 1.0) {
        const double dk = D(ibd);
        const int k = IND(ibd);
        if (dk > 0.0) {
            X(k) = U(k);
            D(ibd) = 0.0;
        } else if (dk < 0.0) {
            X(k) = L(k);
            D(ibd) = 0.0;
        }
    }
    for (int i = 1; i <= nsub; ++i) {
        const int k = IND(i);
        X(k) = X(k) + alpha * D(i);
    }
}

// dcsrch and dcstep, the More-Thuente line search of MINPACK-2.

enum class SearchTask { Start, Fg, Convergence, Warning, Error };

struct SearchState {
    bool brackt = false;
    int stage = 0;
    double ginit = 0.0;
    double gtest = 0.0;
    double gx = 0.0;
    double gy = 0.0;
    double finit = 0.0;
    double fx = 0.0;
    double fy = 0.0;
    double stx = 0.0;
    double sty = 0.0;
    double stmin = 0.0;
    double stmax = 0.0;
    double width = 0.0;
    double width1 = 0.0;
};

void dcstep(double& stx, double& fx, double& dx, double& sty, double& fy, double& dy,
            double& stp, double fp, double dp, bool& brackt, double stpmin, double stpmax) {
    const double p66 = 0.66;
    const double sgnd = dp * (dx / std::abs(dx));
    double stpf = 0.0;

    if (fp > fx) {
        const double theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;
        const double s = std::max({std::abs(theta), std::abs(dx), std::abs(dp)});
        double gamma = s * std::sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
        if (stp < stx) {
            gamma = -gamma;
        }
        const double p = (gamma - dx) + theta;
        const double q = ((gamma - dx) + gamma) + dp;
        const double r = p / q;
        const double stpc = stx + r * (stp - stx);
        const double stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2.0) * (stp - stx);
        if (std::abs(stpc - stx) < std::abs(stpq - stx)) {
            stpf = stpc;
        } else {
            stpf = stpc + (stpq - stpc) / 2.0;
        }
        brackt = true;
    } else if (sgnd < 0.0) {
        const double theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;
        const double s = std::max({std::abs(theta), std::abs(dx), std::abs(dp)});
        double gamma = s * std::sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
        if (stp > stx) {
            gamma = -gamma;
        }
        const double p = (gamma - dp) + theta;
        const double q = ((gamma - dp) + gamma) + dx;
        const double r = p / q;
        const double stpc = stp + r * (stx - stp);
        const double stpq = stp + (dp / (dp - dx)) * (stx - stp);
        if (std::abs(stpc - stp) > std::abs(stpq - stp)) {
            stpf = stpc;
        } else {
            stpf = stpq;
        }
        brackt = true;
    } else if (std::abs(dp) < std::abs(dx)) {
        const double theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;
        const double s = std::max({std::abs(theta), std::abs(dx), std::abs(dp)});
        double gamma =
            s * std::sqrt(std::max(0.0, (theta / s) * (theta / s) - (dx / s) * (dp / s)));
        if (stp > stx) {
            gamma = -gamma;
        }
        const double p = (gamma - dp) + theta;
        const double q = (gamma + (dx - dp)) + gamma;
        const double r = p / q;
        double stpc = 0.0;
        if (r < 0.0 && gamma != 0.0) {
            stpc = stp + r * (stx - stp);
        } else if (stp > stx) {
            stpc = stpmax;
        } else {
            stpc = stpmin;
        }
        const double stpq = stp + (dp / (dp - dx)) * (stx - stp);

        if (brackt) {
            if (std::abs(stpc - stp) < std::abs(stpq - stp)) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
            if (stp > stx) {
                stpf = std::min(stp + p66 * (sty - stp), stpf);
            } else {
                stpf = std::max(stp + p66 * (sty - stp), stpf);
            }
        } else {
            if (std::abs(stpc - stp) > std::abs(stpq - stp)) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
            stpf = std::min(stpmax, stpf);
            stpf = std::max(stpmin, stpf);
        }
    } else {
        if (brackt) {
            const double theta = 3.0 * (fp - fy) / (sty - stp) + dy + dp;
            const double s = std::max({std::abs(theta), std::abs(dy), std::abs(dp)});
            double gamma = s * std::sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
            if (stp > sty) {
                gamma = -gamma;
            }
            const double p = (gamma - dp) + theta;
            const double q = ((gamma - dp) + gamma) + dy;
            const double r = p / q;
            const double stpc = stp + r * (sty - stp);
            stpf = stpc;
        } else if (stp > stx) {
            stpf = stpmax;
        } else {
            stpf = stpmin;
        }
    }

    if (fp > fx) {
        sty = stp;
        fy = fp;
        dy = dp;
    } else {
        if (sgnd < 0.0) {
            sty = stx;
            fy = fx;
            dy = dx;
        }
        stx = stp;
        fx = fp;
        dx = dp;
    }
    stp = stpf;
}

// `task` is both input and output, as in the Fortran: Start begins a search,
// Fg asks for f and g at the returned stp.
void dcsrch(double f, double g, double& stp, double ftol, double gtol, double xtol,
            double stpmin, double stpmax, SearchTask& task, SearchState& s) {
    const double p5 = 0.5;
    const double p66 = 0.66;
    const double xtrapl = 1.1;
    const double xtrapu = 4.0;

    if (task == SearchTask::Start) {
        if (stp < stpmin || stp > stpmax || g >= 0.0 || ftol < 0.0 || gtol < 0.0 ||
            xtol < 0.0 || stpmin < 0.0 || stpmax < stpmin) {
            task = SearchTask::Error;
            return;
        }
        s.brackt = false;
        s.stage = 1;
        s.finit = f;
        s.ginit = g;
        s.gtest = ftol * s.ginit;
        s.width = stpmax - stpmin;
        s.width1 = s.width / p5;
        s.stx = 0.0;
        s.fx = s.finit;
        s.gx = s.ginit;
        s.sty = 0.0;
        s.fy = s.finit;
        s.gy = s.ginit;
        s.stmin = 0.0;
        s.stmax = stp + xtrapu * stp;
        task = SearchTask::Fg;
        return;
    }

    const double ftest = s.finit + stp * s.gtest;
    if (s.stage == 1 && f <= ftest && g >= 0.0) {
        s.stage = 2;
    }

    if (s.brackt && (stp <= s.stmin || stp >= s.stmax)) {
        task = SearchTask::Warning;
    }
    if (s.brackt && s.stmax - s.stmin <= xtol * s.stmax) {
        task = SearchTask::Warning;
    }
    if (stp == stpmax && f <= ftest && g <= s.gtest) {
        task = SearchTask::Warning;
    }
    if (stp == stpmin && (f > ftest || g >= s.gtest)) {
        task = SearchTask::Warning;
    }
    if (f <= ftest && std::abs(g) <= gtol * (-s.ginit)) {
        task = SearchTask::Convergence;
    }
    if (task == SearchTask::Warning || task == SearchTask::Convergence) {
        return;
    }

    if (s.stage == 1 && f <= s.fx && f > ftest) {
        const double fm = f - stp * s.gtest;
        double fxm = s.fx - s.stx * s.gtest;
        double fym = s.fy - s.sty * s.gtest;
        const double gm = g - s.gtest;
        double gxm = s.gx - s.gtest;
        double gym = s.gy - s.gtest;
        dcstep(s.stx, fxm, gxm, s.sty, fym, gym, stp, fm, gm, s.brackt, s.stmin, s.stmax);
        s.fx = fxm + s.stx * s.gtest;
        s.fy = fym + s.sty * s.gtest;
        s.gx = gxm + s.gtest;
        s.gy = gym + s.gtest;
    } else {
        dcstep(s.stx, s.fx, s.gx, s.sty, s.fy, s.gy, stp, f, g, s.brackt, s.stmin, s.stmax);
    }

    if (s.brackt) {
        if (std::abs(s.sty - s.stx) >= p66 * s.width1) {
            stp = s.stx + p5 * (s.sty - s.stx);
        }
        s.width1 = s.width;
        s.width = std::abs(s.sty - s.stx);
    }

    if (s.brackt) {
        s.stmin = std::min(s.stx, s.sty);
        s.stmax = std::max(s.stx, s.sty);
    } else {
        s.stmin = stp + xtrapl * (stp - s.stx);
        s.stmax = stp + xtrapu * (stp - s.stx);
    }

    stp = std::max(stp, stpmin);
    stp = std::min(stp, stpmax);

    if ((s.brackt && (stp <= s.stmin || stp >= s.stmax)) ||
        (s.brackt && s.stmax - s.stmin <= xtol * s.stmax)) {
        stp = s.stx;
    }
    task = SearchTask::Fg;
}

// ---------------------------------------------------------------------------
// scipy's ScalarFunction with a '2-point' absolute step gradient.

class ScalarFunction {
  public:
    ScalarFunction(const std::function<double(std::span<const double>)>& objective,
                   std::vector<double> x0, std::vector<double> lower,
                   std::vector<double> upper, double epsilon)
        : objective_(objective),
          x_(std::move(x0)),
          lower_(std::move(lower)),
          upper_(std::move(upper)),
          epsilon_(epsilon) {
        f_ = evaluate(x_);
        g_ = gradient(x_, f_);
    }

    // fun_and_grad: memoised on exact equality of x, as np.array_equal does.
    void fun_and_grad(const std::vector<double>& x, double& f, std::vector<double>& g) {
        bool same = x.size() == x_.size();
        for (std::size_t i = 0; same && i < x.size(); ++i) {
            same = x[i] == x_[i];
        }
        if (!same) {
            x_ = x;
            f_ = evaluate(x_);
            g_ = gradient(x_, f_);
        }
        f = f_;
        g = g_;
    }

    [[nodiscard]] int nfev() const noexcept { return nfev_; }

  private:
    double evaluate(const std::vector<double>& x) {
        ++nfev_;
        return objective_(std::span<const double>(x.data(), x.size()));
    }

    std::vector<double> gradient(const std::vector<double>& x0, double f0) {
        const std::size_t n = x0.size();
        std::vector<double> h(n, epsilon_);
        // _eps_for_method('2-point') for float64: EPS ** 0.5.
        const double relative_step = std::sqrt(std::numeric_limits<double>::epsilon());
        for (std::size_t i = 0; i < n; ++i) {
            const double sign = x0[i] >= 0.0 ? 1.0 : -1.0;
            if ((x0[i] + h[i]) - x0[i] == 0.0) {
                h[i] = relative_step * sign * std::max(1.0, std::abs(x0[i]));
            }
        }
        // _adjust_scheme_to_bounds(x0, h, 1, '1-sided', lb, ub)
        bool unbounded = true;
        for (std::size_t i = 0; i < n; ++i) {
            unbounded = unbounded && lower_[i] == -std::numeric_limits<double>::infinity() &&
                        upper_[i] == std::numeric_limits<double>::infinity();
        }
        if (!unbounded) {
            for (std::size_t i = 0; i < n; ++i) {
                const double lower_dist = x0[i] - lower_[i];
                const double upper_dist = upper_[i] - x0[i];
                const double stepped = x0[i] + h[i];
                const bool violated = stepped < lower_[i] || stepped > upper_[i];
                const bool fitting = std::abs(h[i]) <= std::max(lower_dist, upper_dist);
                if (violated && fitting) {
                    h[i] = -h[i];
                }
                if (!fitting) {
                    h[i] = upper_dist >= lower_dist ? upper_dist : -lower_dist;
                }
            }
        }
        // _dense_difference, method '2-point'
        std::vector<double> g(n, 0.0);
        std::vector<double> x1 = x0;
        for (std::size_t i = 0; i < n; ++i) {
            x1[i] += h[i];
            const double dx = x1[i] - x0[i];
            const double df = evaluate(x1) - f0;
            g[i] = df / dx;
            x1[i] = x0[i];
        }
        return g;
    }

    const std::function<double(std::span<const double>)>& objective_;
    std::vector<double> x_;
    std::vector<double> lower_;
    std::vector<double> upper_;
    double epsilon_;
    double f_ = 0.0;
    std::vector<double> g_;
    int nfev_ = 0;
};

}  // namespace

LbfgsbResult minimise_lbfgsb_scipy(
    const std::function<double(std::span<const double>)>& objective,
    std::span<const double> x0_in, std::span<const Bound> bounds,
    const LbfgsbOptions& options, int max_line_search_steps) {
    const int n = static_cast<int>(x0_in.size());
    const int m = options.memory;
    LbfgsbResult result;
    if (n == 0) {
        return result;
    }

    // _minimize_lbfgsb: x0 = np.clip(x0, lb, ub); np.clip propagates NaN.
    std::vector<double> lower(static_cast<std::size_t>(n));
    std::vector<double> upper(static_cast<std::size_t>(n));
    std::vector<double> l(static_cast<std::size_t>(n), 0.0);
    std::vector<double> u(static_cast<std::size_t>(n), 0.0);
    std::vector<int> nbd(static_cast<std::size_t>(n), 0);
    std::vector<double> x(x0_in.begin(), x0_in.end());
    for (int i = 0; i < n; ++i) {
        const Bound& bound = bounds[static_cast<std::size_t>(i)];
        lower[static_cast<std::size_t>(i)] = bound.lower;
        upper[static_cast<std::size_t>(i)] = bound.upper;
        double& xi = x[static_cast<std::size_t>(i)];
        if (!std::isnan(xi)) {
            xi = std::min(std::max(xi, bound.lower), bound.upper);
        }
        const bool has_lower = !std::isinf(bound.lower);
        const bool has_upper = !std::isinf(bound.upper);
        if (has_lower) {
            l[static_cast<std::size_t>(i)] = bound.lower;
        }
        if (has_upper) {
            u[static_cast<std::size_t>(i)] = bound.upper;
        }
        nbd[static_cast<std::size_t>(i)] =
            has_lower ? (has_upper ? 2 : 1) : (has_upper ? 3 : 0);
    }

    // The factr scipy hands to setulb is ftol / eps with ftol = factr * eps.
    const double machine_eps = std::numeric_limits<double>::epsilon();
    const double factr = (options.factr * machine_eps) / machine_eps;
    const double pgtol = options.pgtol;

    ScalarFunction function(objective, x, lower, upper, options.epsilon);

    // Workspace, in the shapes mainlb declares.
    std::vector<double> ws(static_cast<std::size_t>(n * m), 0.0);
    std::vector<double> wy(static_cast<std::size_t>(n * m), 0.0);
    std::vector<double> sy(static_cast<std::size_t>(m * m), 0.0);
    std::vector<double> ss(static_cast<std::size_t>(m * m), 0.0);
    std::vector<double> wt(static_cast<std::size_t>(m * m), 0.0);
    std::vector<double> wn(static_cast<std::size_t>(4 * m * m), 0.0);
    std::vector<double> snd(static_cast<std::size_t>(4 * m * m), 0.0);
    std::vector<double> z(static_cast<std::size_t>(n), 0.0);
    std::vector<double> r(static_cast<std::size_t>(n), 0.0);
    std::vector<double> d(static_cast<std::size_t>(n), 0.0);
    std::vector<double> t(static_cast<std::size_t>(n), 0.0);
    std::vector<double> xp(static_cast<std::size_t>(n), 0.0);
    std::vector<double> wa(static_cast<std::size_t>(8 * m), 0.0);
    std::vector<int> index(static_cast<std::size_t>(n), 0);
    std::vector<int> iwhere(static_cast<std::size_t>(n), 0);
    std::vector<int> indx2(static_cast<std::size_t>(n), 0);
    std::vector<double> g(static_cast<std::size_t>(n), 0.0);
    double f = 0.0;

    // mainlb, task = 'START'
    const double epsmch = machine_eps;  // 2 * dlamch('e')
    int col = 0;
    int head = 1;
    double theta = 1.0;
    int iupdat = 0;
    bool updatd = false;
    int iback = 0;
    int itail = 0;
    int iword = 0;
    int nact = 0;
    int ileave = 0;
    int nenter = 0;
    double fold = 0.0;
    double dnorm = 0.0;
    double gd = 0.0;
    double stpmx = 0.0;
    double sbgnrm = 0.0;
    double stp = 0.0;
    double gdold = 0.0;
    double dtd = 0.0;
    int iter = 0;
    int nfgv = 0;
    int nseg = 0;
    int nintol = 0;
    int nskip = 0;
    int nfree = n;
    int ifun = 0;
    const double tol = factr * epsmch;
    int info = 0;
    bool prjctd = false;
    bool cnstnd = false;
    bool boxed = false;
    bool wrk = false;
    SearchState search;
    SearchTask search_task = SearchTask::Start;
    (void)nact;
    (void)nintol;
    (void)nskip;
    (void)nfgv;
    (void)iword;
    (void)prjctd;

    enum class Outcome { Converged, Limit, Abnormal };
    Outcome outcome = Outcome::Abnormal;
    int iterations = 0;  // the driver's n_iterations

    active(n, l.data(), u.data(), nbd.data(), x.data(), iwhere.data(), prjctd, cnstnd, boxed);

    // FG_START
    function.fun_and_grad(x, f, g);
    nfgv = 1;
    projgr(n, l.data(), u.data(), nbd.data(), x.data(), g.data(), sbgnrm);

    const auto refresh_memory = [&]() {
        info = 0;
        col = 0;
        head = 1;
        theta = 1.0;
        iupdat = 0;
        updatd = false;
    };

    if (sbgnrm <= pgtol) {
        outcome = Outcome::Converged;
    } else {
        bool finished = false;
        while (!finished) {
            // 222
            iword = -1;
            bool skip_cauchy = false;
            if (!cnstnd && col > 0) {
                dcopy(n, x.data(), z.data());
                wrk = updatd;
                nseg = 0;
                skip_cauchy = true;
            }
            if (!skip_cauchy) {
                cauchy(n, x.data(), l.data(), u.data(), nbd.data(), g.data(), indx2.data(),
                       iwhere.data(), t.data(), d.data(), z.data(), m, wy.data(), ws.data(),
                       sy.data(), wt.data(), theta, col, head, &wa[0], &wa[2 * m], &wa[4 * m],
                       &wa[6 * m], nseg, sbgnrm, info, epsmch);
                if (info != 0) {
                    refresh_memory();
                    continue;
                }
                nintol = nintol + nseg;
                freev(n, nfree, index.data(), nenter, ileave, indx2.data(), iwhere.data(), wrk,
                      updatd, cnstnd, iter);
                nact = n - nfree;
            }

            // 333
            if (!(nfree == 0 || col == 0)) {
                if (wrk) {
                    formk(n, nfree, index.data(), nenter, ileave, indx2.data(), iupdat, updatd,
                          wn.data(), snd.data(), m, ws.data(), wy.data(), sy.data(), theta, col,
                          head, info);
                }
                if (info != 0) {
                    refresh_memory();
                    continue;
                }
                cmprlb(n, m, x.data(), g.data(), ws.data(), wy.data(), sy.data(), wt.data(),
                       z.data(), r.data(), wa.data(), index.data(), theta, col, head, nfree,
                       cnstnd, info);
                if (info == 0) {
                    subsm(n, m, nfree, index.data(), l.data(), u.data(), nbd.data(), z.data(),
                          r.data(), xp.data(), ws.data(), wy.data(), theta, x.data(), g.data(),
                          col, head, iword, wa.data(), wn.data(), info);
                }
                // 444
                if (info != 0) {
                    refresh_memory();
                    continue;
                }
            }

            // 555: the search direction d = z - x, then the line search.
            for (int i = 0; i < n; ++i) {
                d[static_cast<std::size_t>(i)] =
                    z[static_cast<std::size_t>(i)] - x[static_cast<std::size_t>(i)];
            }

            bool line_search_start = true;
            bool restart = false;
            while (true) {
                // lnsrlb
                if (line_search_start) {
                    dtd = ddot(n, d.data(), d.data());
                    dnorm = std::sqrt(dtd);
                    stpmx = 1.0e10;
                    if (cnstnd) {
                        if (iter == 0) {
                            stpmx = 1.0;
                        } else {
                            for (int i = 0; i < n; ++i) {
                                const double a1 = d[static_cast<std::size_t>(i)];
                                if (nbd[static_cast<std::size_t>(i)] != 0) {
                                    if (a1 < 0.0 && nbd[static_cast<std::size_t>(i)] <= 2) {
                                        const double a2 = l[static_cast<std::size_t>(i)] -
                                                          x[static_cast<std::size_t>(i)];
                                        if (a2 >= 0.0) {
                                            stpmx = 0.0;
                                        } else if (a1 * stpmx < a2) {
                                            stpmx = a2 / a1;
                                        }
                                    } else if (a1 > 0.0 && nbd[static_cast<std::size_t>(i)] >= 2) {
                                        const double a2 = u[static_cast<std::size_t>(i)] -
                                                          x[static_cast<std::size_t>(i)];
                                        if (a2 <= 0.0) {
                                            stpmx = 0.0;
                                        } else if (a1 * stpmx > a2) {
                                            stpmx = a2 / a1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (iter == 0 && !boxed) {
                        stp = std::min(1.0 / dnorm, stpmx);
                    } else {
                        stp = 1.0;
                    }
                    dcopy(n, x.data(), t.data());
                    dcopy(n, g.data(), r.data());
                    fold = f;
                    ifun = 0;
                    iback = 0;
                    search_task = SearchTask::Start;
                    line_search_start = false;
                }
                // 556
                gd = ddot(n, g.data(), d.data());
                bool new_x = false;
                if (ifun == 0) {
                    gdold = gd;
                    if (gd >= 0.0) {
                        // ascent direction in projection
                        info = -4;
                    }
                }
                if (info == 0) {
                    dcsrch(f, gd, stp, 1.0e-3, 0.9, 0.1, 0.0, stpmx, search_task, search);
                    if (search_task != SearchTask::Convergence &&
                        search_task != SearchTask::Warning) {
                        ifun = ifun + 1;
                        nfgv = nfgv + 1;
                        iback = ifun - 1;
                        if (stp == 1.0) {
                            dcopy(n, z.data(), x.data());
                        } else {
                            for (int i = 0; i < n; ++i) {
                                const std::size_t k = static_cast<std::size_t>(i);
                                x[k] = stp * d[k] + t[k];
                                if (nbd[k] == 1 || nbd[k] == 2) {
                                    x[k] = std::max(x[k], l[k]);
                                }
                                if (nbd[k] == 2 || nbd[k] == 3) {
                                    x[k] = std::min(x[k], u[k]);
                                }
                            }
                        }
                        // An 'ERROR' from dcsrch is not CONV or WARN either,
                        // so the Fortran asks for f and g here too.
                        if (search_task == SearchTask::Error) {
                            search_task = SearchTask::Fg;
                        }
                    } else {
                        new_x = true;
                    }
                }
                // back in mainlb
                if (info != 0 || iback >= max_line_search_steps) {
                    dcopy(n, t.data(), x.data());
                    dcopy(n, r.data(), g.data());
                    f = fold;
                    if (col == 0) {
                        if (info == 0) {
                            info = -9;
                            nfgv = nfgv - 1;
                            ifun = ifun - 1;
                            iback = iback - 1;
                        }
                        iter = iter + 1;
                        outcome = Outcome::Abnormal;
                        finished = true;
                    } else {
                        if (info == 0) {
                            nfgv = nfgv - 1;
                        }
                        refresh_memory();
                        restart = true;
                    }
                    break;
                }
                if (!new_x) {
                    function.fun_and_grad(x, f, g);
                    continue;  // 666
                }
                iter = iter + 1;
                projgr(n, l.data(), u.data(), nbd.data(), x.data(), g.data(), sbgnrm);
                break;
            }
            if (finished) {
                break;
            }
            if (restart) {
                continue;
            }

            // NEW_X returns to the driver, which counts the iteration and
            // applies its own limits before calling setulb again.
            ++iterations;
            if (iterations >= options.max_iterations ||
                function.nfev() > options.max_function_evaluations) {
                outcome = Outcome::Limit;
                break;
            }

            // 777
            if (sbgnrm <= pgtol) {
                outcome = Outcome::Converged;
                break;
            }
            {
                const double ddum = std::max({std::abs(fold), std::abs(f), 1.0});
                if ((fold - f) <= tol * ddum) {
                    outcome = Outcome::Converged;
                    break;
                }
            }

            for (int i = 0; i < n; ++i) {
                const std::size_t k = static_cast<std::size_t>(i);
                r[k] = g[k] - r[k];
            }
            const double rr = ddot(n, r.data(), r.data());
            double dr = 0.0;
            double ddum = 0.0;
            if (stp == 1.0) {
                dr = gd - gdold;
                ddum = -gdold;
            } else {
                dr = (gd - gdold) * stp;
                dscal(n, stp, d.data());
                ddum = -gdold * stp;
            }

            if (dr <= epsmch * ddum) {
                nskip = nskip + 1;
                updatd = false;
                continue;  // 888
            }

            updatd = true;
            iupdat = iupdat + 1;
            matupd(n, m, ws.data(), wy.data(), sy.data(), ss.data(), d.data(), r.data(), itail,
                   iupdat, col, head, theta, rr, dr, stp, dtd);
            formt(m, wt.data(), sy.data(), ss.data(), col, theta, info);
            if (info != 0) {
                refresh_memory();
            }
        }
    }

    result.x = x;
    result.f = f;
    result.iterations = iterations;
    result.function_evaluations = function.nfev();
    switch (outcome) {
        case Outcome::Converged:
            result.status = 0;
            break;
        case Outcome::Limit:
            result.status = 1;
            break;
        case Outcome::Abnormal:
            // warnflag: nfev > maxfun or n_iterations >= maxiter count as 1
            // even when the Fortran stopped for another reason.
            result.status = (function.nfev() > options.max_function_evaluations ||
                             iterations >= options.max_iterations)
                                ? 1
                                : 2;
            break;
    }
    return result;
}

NBinomFit fit_nbinom_scipy(std::span<const double> data) {
    constexpr double infinitesimal = std::numeric_limits<double>::epsilon();
    NBinomFit fit;
    const std::size_t n = data.size();
    if (n == 0) {
        fit.size = 10.0;
        fit.prob = std::numeric_limits<double>::quiet_NaN();
        fit.status = 2;
        return fit;
    }

    // np.sum(np.log(factorial(X))). scipy's factorial is gamma(x + 1) for a
    // float array and 0 for x < 0.
    std::vector<double> log_factorial(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double gamma = data[i] < 0.0 ? 0.0 : std::tgamma(data[i] + 1.0);
        log_factorial[i] = std::log(gamma);
    }
    const double log_factorial_sum = npy::pairwise_sum(log_factorial);

    std::vector<double> scratch(n, 0.0);
    const std::function<double(std::span<const double>)> negative_log_likelihood =
        [&](std::span<const double> parameters) -> double {
        const double r = parameters[0];
        const double p = parameters[1];
        for (std::size_t i = 0; i < n; ++i) {
            scratch[i] = gammaln(data[i] + r);
        }
        const double gammaln_sum = npy::pairwise_sum(scratch);
        const double log_one_minus_p = std::log(1.0 - (p < 1.0 ? p : 1.0 - infinitesimal));
        for (std::size_t i = 0; i < n; ++i) {
            scratch[i] = data[i] * log_one_minus_p;
        }
        const double tail = npy::pairwise_sum(scratch);
        const double count = static_cast<double>(n);
        const double result =
            gammaln_sum - log_factorial_sum - count * gammaln(r) + count * r * std::log(p) + tail;
        return -result;
    };

    // R's fitdistr moment estimator, as fit_nbinom computes it.
    const double mean = npy::pairwise_sum(data.data(), n) / static_cast<double>(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double difference = data[i] - mean;
        scratch[i] = difference * difference;
    }
    const double variance = npy::pairwise_sum(scratch) / static_cast<double>(n);
    const double size = variance > mean ? (mean * mean) / (variance - mean) : 10.0;
    const double prob = size / ((size + mean) != 0.0 ? (size + mean) : 1.0);

    const std::vector<double> start{size, prob};
    const std::vector<Bound> bounds{
        Bound{infinitesimal, std::numeric_limits<double>::infinity()},
        Bound{infinitesimal, 1.0}};

    if (std::isinf(log_factorial_sum)) {
        // The objective is +inf at every parameter pair. scipy's forward
        // difference is then inf - inf = NaN, and the NaN propagates through
        // the Cauchy point and the line search until the search is abandoned
        // and the starting point restored: measured, fmin_l_bfgs_b returns x0
        // bit for bit on both such distributions of the cHi-C test data. The
        // NaN path itself depends on how gfortran's MIN and MAX treat a NaN
        // argument, which C++ does not share, so the result is returned
        // directly rather than by walking that path.
        fit.size = std::min(std::max(size, bounds[0].lower), bounds[0].upper);
        fit.prob = std::min(std::max(prob, bounds[1].lower), bounds[1].upper);
        fit.status = 2;
        return fit;
    }

    const LbfgsbResult solution =
        minimise_lbfgsb_scipy(negative_log_likelihood, start, bounds, LbfgsbOptions());
    fit.size = solution.x[0];
    fit.prob = solution.x[1];
    fit.status = solution.status;
    fit.iterations = solution.iterations;
    return fit;
}

}  // namespace hicx::stats
