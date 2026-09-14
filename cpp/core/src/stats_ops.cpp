#include "hicx/stats_ops.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "hicx/numpy_compat.hpp"

namespace hicx::stats {

namespace {

// --------------------------------------------------------------------------
// Cephes, translated from scipy/special/special/cephes/{ndtr,psi}.h. The
// coefficient tables and the branch structure are copied so that a value
// printed with twelve decimals matches scipy's. Cephes Math Library Release
// 2.9, Copyright 1984 to 2000 Stephen L. Moshier; redistributed by SciPy under
// the BSD licence.

constexpr double kMaxLog = 7.09782712893383996732E2;
constexpr double kEuler = 0.577215664901532860606512090082402431;
// np.finfo(np.float64).eps, the bound and the guard fit_nbinom uses.
constexpr double kMachineEpsilonStats = 2.220446049250313e-16;

double polevl(double x, const double* coefficients, int degree) {
    double answer = coefficients[0];
    for (int i = 1; i <= degree; ++i) {
        answer = answer * x + coefficients[i];
    }
    return answer;
}

// polevl for a polynomial whose leading coefficient is 1 and therefore absent
// from the table.
double p1evl(double x, const double* coefficients, int degree) {
    double answer = x + coefficients[0];
    for (int i = 1; i < degree; ++i) {
        answer = answer * x + coefficients[i];
    }
    return answer;
}

constexpr double kNdtrP[] = {
    2.46196981473530512524E-10, 5.64189564831068821977E-1, 7.46321056442269912687E0,
    4.86371970985681366614E1,   1.96520832956077098242E2,  5.26445194995477358631E2,
    9.34528527171957607540E2,   1.02755188689515710272E3,  5.57535335369399327526E2};

constexpr double kNdtrQ[] = {
    1.32281951154744992508E1, 8.67072140885989742329E1, 3.54937778887819891062E2,
    9.75708501743205489753E2, 1.82390916687909736289E3, 2.24633760818710981792E3,
    1.65666309194161350182E3, 5.57535340817727675546E2};

constexpr double kNdtrR[] = {5.64189583547755073984E-1, 1.27536670759978104416E0,
                             5.01905042251180477414E0,  6.16021097993053585195E0,
                             7.40974269950448939160E0,  2.97886665372100240670E0};

constexpr double kNdtrS[] = {2.26052863220117276590E0, 9.39603524938001434673E0,
                             1.20489539808096656605E1, 1.70814450747565897222E1,
                             9.60896809063285878198E0, 3.36907645100081516050E0};

constexpr double kNdtrT[] = {9.60497373987051638749E0, 9.00260197203842689217E1,
                             2.23200534594684319226E3, 7.00332514112805075473E3,
                             5.55923013010394962768E4};

constexpr double kNdtrU[] = {3.35617141647503099647E1, 5.21357949780152679795E2,
                             4.59432382970980127987E3, 2.26290000613890934246E4,
                             4.92673942608635921086E4};

constexpr double kPsiA[] = {
    8.33333333333333333333E-2,  -2.10927960927960927961E-2, 7.57575757575757575758E-3,
    -4.16666666666666666667E-3, 3.96825396825396825397E-3,  -8.33333333333333333333E-3,
    8.33333333333333333333E-2};

constexpr float kPsiY = 0.99558162689208984f;
constexpr double kPsiRoot1 = 1569415565.0 / 1073741824.0;
constexpr double kPsiRoot2 = (381566830.0 / 1073741824.0) / 1073741824.0;
constexpr double kPsiRoot3 = 0.9016312093258695918615325266959189453125e-19;

constexpr double kPsiP[] = {-0.0020713321167745952, -0.045251321448739056,
                            -0.28919126444774784,   -0.65031853770896507,
                            -0.32555031186804491,   0.25479851061131551};
constexpr double kPsiQ[] = {-0.55789841321675513e-6, 0.0021284987017821144,
                            0.054151797245674225,   0.43593529692665969,
                            1.4606242909763515,     2.0767117023730469,
                            1.0};

double digamma_imp_1_2(double x) {
    double g = x - kPsiRoot1;
    g -= kPsiRoot2;
    g -= kPsiRoot3;
    const double r = polevl(x - 1.0, kPsiP, 5) / polevl(x - 1.0, kPsiQ, 6);
    return g * static_cast<double>(kPsiY) + g * r;
}

double psi_asy(double x) {
    double y = 0.0;
    if (x < 1.0e17) {
        const double z = 1.0 / (x * x);
        y = z * polevl(z, kPsiA, 6);
    }
    return std::log(x) - (0.5 / x) - y;
}

// --------------------------------------------------------------------------
// incbet, translated from scipy/special/special/cephes/incbet.h. The two
// continued fractions, the power series and the thresholds that pick between
// them are the cephes ones unchanged; `Gamma` and `lgam` are libm's tgamma and
// lgamma_r rather than a second vendored table, which is where the two
// implementations may differ, by at most an ulp.

constexpr double kMachEp = 1.11022302462515654042E-16;
constexpr double kMaxGam = 171.624376956302725;
constexpr double kMinLog = -7.08396418532264106224E2;
constexpr double kIncbetBig = 4.503599627370496e15;
constexpr double kIncbetBigInv = 2.22044604925031308085e-16;

// Power series for the incomplete beta integral, cephes pseries. Used when
// b * x is small, where the continued fractions converge slowly.
double incbet_pseries(double a, double b, double x) {
    const double ai = 1.0 / a;
    double u = (1.0 - b) * x;
    double v = u / (a + 1.0);
    const double t1 = v;
    double t = u;
    double n = 2.0;
    double s = 0.0;
    const double z = kMachEp * ai;

    while (std::abs(v) > z) {
        u = (n - b) * x / n;
        t *= u;
        v = t / (a + n);
        s += v;
        n += 1.0;
    }
    s += t1;
    s += ai;

    u = a * std::log(x);
    if ((a + b) < kMaxGam && std::abs(u) < kMaxLog) {
        t = std::tgamma(a + b) / (std::tgamma(a) * std::tgamma(b));
        s = s * t * std::pow(x, a);
    } else {
        t = gammaln(a + b) - gammaln(a) - gammaln(b) + u + std::log(s);
        s = t < kMinLog ? 0.0 : std::exp(t);
    }
    return s;
}

// Continued fraction expansion #1 for the incomplete beta integral.
double incbet_cf1(double a, double b, double x) {
    double k1 = a;
    double k2 = a + b;
    double k3 = a;
    double k4 = a + 1.0;
    double k5 = 1.0;
    double k6 = b - 1.0;
    double k7 = k4;
    double k8 = a + 2.0;

    double pkm2 = 0.0;
    double qkm2 = 1.0;
    double pkm1 = 1.0;
    double qkm1 = 1.0;
    double ans = 1.0;
    double r = 1.0;
    const double thresh = 3.0 * kMachEp;

    for (int n = 0; n < 300; ++n) {
        double xk = -(x * k1 * k2) / (k3 * k4);
        double pk = pkm1 + pkm2 * xk;
        double qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        xk = (x * k5 * k6) / (k7 * k8);
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        if (qk != 0.0) {
            r = pk / qk;
        }
        double t = 1.0;
        if (r != 0.0) {
            t = std::abs((ans - r) / r);
            ans = r;
        }
        if (t < thresh) {
            return ans;
        }

        k1 += 1.0;
        k2 += 1.0;
        k3 += 2.0;
        k4 += 2.0;
        k5 += 1.0;
        k6 -= 1.0;
        k7 += 2.0;
        k8 += 2.0;

        if ((std::abs(qk) + std::abs(pk)) > kIncbetBig) {
            pkm2 *= kIncbetBigInv;
            pkm1 *= kIncbetBigInv;
            qkm2 *= kIncbetBigInv;
            qkm1 *= kIncbetBigInv;
        }
        if ((std::abs(qk) < kIncbetBigInv) || (std::abs(pk) < kIncbetBigInv)) {
            pkm2 *= kIncbetBig;
            pkm1 *= kIncbetBig;
            qkm2 *= kIncbetBig;
            qkm1 *= kIncbetBig;
        }
    }
    return ans;
}

// Continued fraction expansion #2 for the incomplete beta integral.
double incbet_cf2(double a, double b, double x) {
    double k1 = a;
    double k2 = b - 1.0;
    double k3 = a;
    double k4 = a + 1.0;
    double k5 = 1.0;
    double k6 = a + b;
    double k7 = a + 1.0;
    double k8 = a + 2.0;

    double pkm2 = 0.0;
    double qkm2 = 1.0;
    double pkm1 = 1.0;
    double qkm1 = 1.0;
    const double z = x / (1.0 - x);
    double ans = 1.0;
    double r = 1.0;
    const double thresh = 3.0 * kMachEp;

    for (int n = 0; n < 300; ++n) {
        double xk = -(z * k1 * k2) / (k3 * k4);
        double pk = pkm1 + pkm2 * xk;
        double qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        xk = (z * k5 * k6) / (k7 * k8);
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        if (qk != 0.0) {
            r = pk / qk;
        }
        double t = 1.0;
        if (r != 0.0) {
            t = std::abs((ans - r) / r);
            ans = r;
        }
        if (t < thresh) {
            return ans;
        }

        k1 += 1.0;
        k2 -= 1.0;
        k3 += 2.0;
        k4 += 2.0;
        k5 += 1.0;
        k6 += 1.0;
        k7 += 2.0;
        k8 += 2.0;

        if ((std::abs(qk) + std::abs(pk)) > kIncbetBig) {
            pkm2 *= kIncbetBigInv;
            pkm1 *= kIncbetBigInv;
            qkm2 *= kIncbetBigInv;
            qkm1 *= kIncbetBigInv;
        }
        if ((std::abs(qk) < kIncbetBigInv) || (std::abs(pk) < kIncbetBigInv)) {
            pkm2 *= kIncbetBig;
            pkm1 *= kIncbetBig;
            qkm2 *= kIncbetBig;
            qkm1 *= kIncbetBig;
        }
    }
    return ans;
}

}  // namespace

double betainc(double a, double b, double x) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (std::isnan(a) || std::isnan(b) || std::isnan(x)) {
        return nan;
    }
    if (a <= 0.0 || b <= 0.0) {
        return nan;
    }
    if (x <= 0.0 || x >= 1.0) {
        if (x == 0.0) {
            return 0.0;
        }
        if (x == 1.0) {
            return 1.0;
        }
        return nan;
    }

    double t = 0.0;
    bool reflected = false;

    if (b * x <= 1.0 && x <= 0.95) {
        t = incbet_pseries(a, b, x);
    } else {
        double w = 1.0 - x;
        double aa = a;
        double bb = b;
        double xx = x;
        double xc = w;
        // Reverse a and b when x is greater than the mean, so that the
        // expansion is always evaluated on the convergent side.
        if (x > (a / (a + b))) {
            reflected = true;
            aa = b;
            bb = a;
            xc = x;
            xx = w;
        }

        if (reflected && (bb * xx) <= 1.0 && xx <= 0.95) {
            t = incbet_pseries(aa, bb, xx);
        } else {
            const double y0 = xx * (aa + bb - 2.0) - (aa - 1.0);
            w = y0 < 0.0 ? incbet_cf1(aa, bb, xx) : incbet_cf2(aa, bb, xx) / xc;

            // Multiply w by x^a (1-x)^b Gamma(a+b) / (a Gamma(a) Gamma(b)).
            double y = aa * std::log(xx);
            double s = bb * std::log(xc);
            if ((aa + bb) < kMaxGam && std::abs(y) < kMaxLog &&
                std::abs(s) < kMaxLog) {
                t = std::pow(xc, bb);
                t *= std::pow(xx, aa);
                t /= aa;
                t *= w;
                t *= std::tgamma(aa + bb) / (std::tgamma(aa) * std::tgamma(bb));
            } else {
                y += s + gammaln(aa + bb) - gammaln(aa) - gammaln(bb);
                y += std::log(w / aa);
                t = y < kMinLog ? 0.0 : std::exp(y);
            }
        }
    }

    if (reflected) {
        t = t <= kMachEp ? 1.0 - kMachEp : 1.0 - t;
    }
    return t;
}

double nbinom_sf(double x, double r, double p) {
    // hicexplorer/lib/cnb.py:29 plus hicDetectLoops.py:163, literally.
    return 1.0 - betainc(r, x + 1.0, p);
}

double erfc(double a) {
    if (std::isnan(a)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double x = a < 0.0 ? -a : a;
    if (x < 1.0) {
        return 1.0 - erf(a);
    }
    double z = -a * a;
    if (z < -kMaxLog) {
        return a < 0 ? 2.0 : 0.0;
    }
    z = std::exp(z);
    double p = 0.0;
    double q = 0.0;
    if (x < 8.0) {
        p = polevl(x, kNdtrP, 8);
        q = p1evl(x, kNdtrQ, 8);
    } else {
        p = polevl(x, kNdtrR, 5);
        q = p1evl(x, kNdtrS, 6);
    }
    double y = (z * p) / q;
    if (a < 0) {
        y = 2.0 - y;
    }
    if (y != 0.0) {
        return y;
    }
    return a < 0 ? 2.0 : 0.0;
}

double erf(double x) {
    if (std::isnan(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (x < 0.0) {
        return -erf(-x);
    }
    if (std::abs(x) > 1.0) {
        return 1.0 - erfc(x);
    }
    const double z = x * x;
    return x * polevl(z, kNdtrT, 4) / p1evl(z, kNdtrU, 5);
}

double ndtr(double a) {
    if (std::isnan(a)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double x = a * M_SQRT1_2;
    const double z = std::abs(x);
    double y = 0.0;
    if (z < M_SQRT1_2) {
        y = 0.5 + 0.5 * erf(x);
    } else {
        y = 0.5 * erfc(z);
        if (x > 0) {
            y = 1.0 - y;
        }
    }
    return y;
}

double normal_sf(double x) { return ndtr(-x); }

double gammaln(double x) {
    // scipy.special.gammaln for a positive real argument is cephes lgam, which
    // agrees with glibc's lgamma to within one unit in the last place. The
    // difference is 1e-16 relative, twelve orders of magnitude below the
    // optimiser difference measured in fit_nbinom, so the libm routine is used
    // rather than a second vendored table.
    //
    // lgamma_r rather than lgamma: glibc's lgamma writes the sign of the
    // gamma function into the global `signgam`, which is a data race as soon
    // as two threads fit two distributions at once. hicDetectLoops will do
    // exactly that.
#if defined(__GLIBC__) || defined(_POSIX_VERSION)
    int sign = 0;
    return ::lgamma_r(x, &sign);
#else
    return std::lgamma(x);
#endif
}

double digamma(double x) {
    double y = 0.0;
    if (std::isnan(x)) {
        return x;
    }
    if (x == std::numeric_limits<double>::infinity()) {
        return x;
    }
    if (x == -std::numeric_limits<double>::infinity()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (x == 0.0) {
        return std::copysign(std::numeric_limits<double>::infinity(), -x);
    }
    if (x < 0.0) {
        double q = 0.0;
        const double r = std::modf(x, &q);
        if (r == 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        y = -M_PI / std::tan(M_PI * r);
        x = 1.0 - x;
    }
    if (x <= 10.0 && x == std::floor(x)) {
        const int n = static_cast<int>(x);
        for (int i = 1; i < n; ++i) {
            y += 1.0 / i;
        }
        y -= kEuler;
        return y;
    }
    if (x < 1.0) {
        y -= 1.0 / x;
        x += 1.0;
    } else if (x < 10.0) {
        while (x > 2.0) {
            x -= 1.0;
            y += 1.0 / x;
        }
    }
    if (x >= 1.0 && x <= 2.0) {
        return y + digamma_imp_1_2(x);
    }
    return y + psi_asy(x);
}

// --------------------------------------------------------------------------

std::vector<double> rankdata_average(std::span<const double> values) {
    const std::size_t n = values.size();
    std::vector<double> ranks(n, 0.0);
    if (n == 0) {
        return ranks;
    }
    // np.argsort(kind='quicksort') is not stable, but rankdata's result does
    // not depend on the order within a tie group, so a stable sort gives the
    // same ranks and is deterministic.
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(),
                     [&values](std::size_t a, std::size_t b) {
                         return values[a] < values[b];
                     });

    std::size_t i = 0;
    while (i < n) {
        std::size_t j = i + 1;
        while (j < n && values[order[j]] == values[order[i]]) {
            ++j;
        }
        // Ranks i+1 .. j average to (i + j + 1) / 2.
        const double rank = 0.5 * (static_cast<double>(i) + static_cast<double>(j) + 1.0);
        for (std::size_t k = i; k < j; ++k) {
            ranks[order[k]] = rank;
        }
        i = j;
    }
    return ranks;
}

RanksumsResult ranksums(std::span<const double> x, std::span<const double> y) {
    RanksumsResult result;
    const std::size_t n1 = x.size();
    const std::size_t n2 = y.size();
    if (n1 == 0 && n2 == 0) {
        result.defined = false;
        result.statistic = std::numeric_limits<double>::quiet_NaN();
        result.pvalue = std::numeric_limits<double>::quiet_NaN();
        return result;
    }

    std::vector<double> all;
    all.reserve(n1 + n2);
    all.insert(all.end(), x.begin(), x.end());
    all.insert(all.end(), y.begin(), y.end());
    const std::vector<double> ranked = rankdata_average(all);

    // np.sum over the first n1 ranks, in numpy's reduction order.
    const double s = npy::pairwise_sum(ranked.data(), n1);
    const double total = static_cast<double>(n1) + static_cast<double>(n2);
    const double expected = static_cast<double>(n1) * (total + 1.0) / 2.0;
    const double variance =
        static_cast<double>(n1) * static_cast<double>(n2) * (total + 1.0) / 12.0;
    const double z = (s - expected) / std::sqrt(variance);
    result.statistic = z;
    result.pvalue = 2.0 * normal_sf(std::abs(z));
    return result;
}

double benjamini_hochberg_cutoff(std::vector<double> pvalues, double q) {
    if (pvalues.empty()) {
        return 0.0;
    }
    std::sort(pvalues.begin(), pvalues.end());
    const double n = static_cast<double>(pvalues.size());
    double largest = 0.0;
    for (std::size_t i = 0; i < pvalues.size(); ++i) {
        const double p = pvalues[i];
        if (p <= q * (static_cast<double>(i) + 1.0) / n) {
            if (p >= largest) {
                largest = p;
            }
        }
    }
    return largest;
}

std::vector<double> benjamini_hochberg_adjusted(std::span<const double> pvalues) {
    std::vector<double> out(pvalues.begin(), pvalues.end());
    std::vector<std::size_t> order;
    order.reserve(pvalues.size());
    for (std::size_t i = 0; i < pvalues.size(); ++i) {
        if (!std::isnan(pvalues[i])) {
            order.push_back(i);
        }
    }
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return pvalues[a] < pvalues[b]; });
    const double m = static_cast<double>(order.size());
    std::vector<double> adjusted(order.size());
    for (std::size_t k = 0; k < order.size(); ++k) {
        adjusted[k] = pvalues[order[k]] * m / static_cast<double>(k + 1);
    }
    for (std::size_t k = adjusted.size(); k > 1; --k) {
        adjusted[k - 2] = std::min(adjusted[k - 2], adjusted[k - 1]);
    }
    for (std::size_t k = 0; k < order.size(); ++k) {
        out[order[k]] = std::min(adjusted[k], 1.0);
    }
    return out;
}

std::vector<double> bonferroni_adjusted(std::span<const double> pvalues) {
    std::vector<double> out(pvalues.begin(), pvalues.end());
    const double m = static_cast<double>(
        std::count_if(pvalues.begin(), pvalues.end(), [](double p) { return !std::isnan(p); }));
    for (double& p : out) {
        if (!std::isnan(p)) {
            p = std::min(p * m, 1.0);
        }
    }
    return out;
}

// --------------------------------------------------------------------------
// The float32 flavour of fit_nbinom. See the comment on NBinomPrecision for
// why this exists; the float64 flavour lives beside the optimiser in
// lbfgsb.cpp and this one is deliberately a separate function rather than a
// branch inside it, so that the tools that do not need it are untouched.

NBinomFit fit_nbinom(std::span<const double> data, NBinomPrecision precision) {
    if (precision == NBinomPrecision::Float64) {
        return fit_nbinom(data);
    }
    NBinomFit fit;
    const std::size_t n = data.size();
    if (n == 0) {
        fit.status = 2;
        return fit;
    }

    // The array as scipy holds it. The caller has already rounded the values
    // to float32; this is the array the ufuncs actually see.
    std::vector<float> x(n, 0.0F);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = static_cast<float>(data[i]);
    }

    // np.sum(np.log(factorial(X))): scipy's factorial returns a float64 array
    // even for a float32 argument, but it evaluates gamma(x + 1) in the
    // float32 loop and widens the result, so it **overflows above about 34.6**
    // rather than above 170.6. That is not a detail: as soon as one obs/exp
    // value at a genomic distance exceeds it, this term is infinite, the whole
    // objective is infinite for every parameter pair, the forward difference
    // gradient is NaN and fmin_l_bfgs_b returns its starting point. Measured
    // on the corpus: 3 of the 199 distances of gm12878_chr1.cool are in that
    // state, and computing the term in float64 instead lets the optimiser run
    // on those three and call a loop the reference does not, which is exactly
    // how this was found.
    std::vector<double> log_factorial(n, 0.0);
    bool factorial_overflowed = false;
    for (std::size_t i = 0; i < n; ++i) {
        const float gamma =
            static_cast<float>(std::tgamma(static_cast<double>(x[i] + 1.0F)));
        if (std::isinf(gamma)) {
            factorial_overflowed = true;
            break;
        }
        log_factorial[i] = std::log(static_cast<double>(gamma));
    }

    const double log_factorial_sum =
        factorial_overflowed ? std::numeric_limits<double>::infinity()
                             : npy::pairwise_sum(log_factorial);

    // gammaln is by far the most expensive thing in the objective and it is
    // evaluated over the same array on every iteration, so the distinct values
    // are collected once and each objective evaluation looks the result up
    // instead of recomputing it. This is exact, not an approximation: gammaln
    // is a function of its argument, and the values are scattered back into
    // their original positions before the reduction, so the summation order is
    // untouched and the result is bit identical to the direct loop.
    //
    // It pays because an obs/exp value is a small integer count divided by a
    // per distance constant, so the same value recurs constantly. Measured on
    // gm12878_chr1.cool, the six sampled distance distributions hold 130,256
    // values of which 4,469 are distinct, 3.4 percent. The table is only built
    // when it at least halves the work; on a distribution of genuinely
    // distinct values, such as the 19 distances of the GSE63525 cool where the
    // ratio is 99.9 percent, it is skipped and the direct loop runs.
    std::vector<float> unique;
    std::vector<std::int32_t> slot;
    bool use_table = false;
    {
        bool finite = true;
        for (std::size_t i = 0; i < n && finite; ++i) {
            finite = std::isfinite(x[i]);
        }
        if (finite) {
            unique = x;
            std::sort(unique.begin(), unique.end());
            unique.erase(std::unique(unique.begin(), unique.end()), unique.end());
            if (unique.size() * 2 <= n) {
                slot.resize(n);
                for (std::size_t i = 0; i < n; ++i) {
                    slot[i] = static_cast<std::int32_t>(std::distance(
                        unique.begin(),
                        std::lower_bound(unique.begin(), unique.end(), x[i])));
                }
                use_table = true;
            }
        }
        if (!use_table) {
            unique.clear();
            unique.shrink_to_fit();
        }
    }
    std::vector<float> table(unique.size(), 0.0F);

    std::vector<float> scratch(n, 0.0F);
    const auto negative_log_likelihood =
        [&](std::span<const double> parameters) -> double {
        const double r = parameters[0];
        const double p = parameters[1];
        const double safe_p = p < 1.0 ? p : 1.0 - kMachineEpsilonStats;

        // gammaln(X + r): the addition runs in the float32 loop, with the
        // scalar cast to float32 first, and so does the reduction.
        const float r32 = static_cast<float>(r);
        if (use_table) {
            for (std::size_t u = 0; u < unique.size(); ++u) {
                table[u] =
                    static_cast<float>(gammaln(static_cast<double>(unique[u] + r32)));
            }
            for (std::size_t i = 0; i < n; ++i) {
                scratch[i] = table[static_cast<std::size_t>(slot[i])];
            }
        } else {
            for (std::size_t i = 0; i < n; ++i) {
                scratch[i] =
                    static_cast<float>(gammaln(static_cast<double>(x[i] + r32)));
            }
        }
        const float gammaln_sum = npy::pairwise_sum(scratch.data(), n);

        // X * log(1 - p), likewise float32 throughout.
        const float log1p32 = static_cast<float>(std::log(1.0 - safe_p));
        for (std::size_t i = 0; i < n; ++i) {
            scratch[i] = x[i] * log1p32;
        }
        const float tail = npy::pairwise_sum(scratch.data(), n);

        // The five terms are combined left to right, and the float32 partial
        // sums widen to float64 as they meet a float64 term.
        double value = static_cast<double>(gammaln_sum);
        value -= log_factorial_sum;
        value -= static_cast<double>(n) * gammaln(r);
        value += static_cast<double>(n) * r * std::log(p);
        value += static_cast<double>(tail);
        return -value;
    };

    // np.mean and np.var of a float32 array reduce in float32 as well, so the
    // moment estimator that seeds the optimiser is a float32 quantity.
    const float mean32 =
        npy::pairwise_sum(x.data(), n) / static_cast<float>(n);
    for (std::size_t i = 0; i < n; ++i) {
        const float difference = x[i] - mean32;
        scratch[i] = difference * difference;
    }
    const float variance32 =
        npy::pairwise_sum(scratch.data(), n) / static_cast<float>(n);

    double start_size = 0.0;
    double start_prob = 0.0;
    if (variance32 > mean32) {
        // `size = (m ** 2) / (v - m)`. `m ** 2` with a Python int exponent
        // promotes a numpy float32 scalar to float64 in numpy 1.26, while
        // `v - m` stays float32, so the quotient and everything after it is
        // float64 even though the data is not. `m * m` would have stayed
        // float32; the difference is the eighth significant digit of the
        // starting point, which is enough to move the fit by 1e-08.
        const double mean = static_cast<double>(mean32);
        start_size = (mean * mean) / static_cast<double>(variance32 - mean32);
        const double denominator = start_size + mean;
        start_prob = denominator != 0.0 ? start_size / denominator : start_size;
    } else {
        // `size = 10` is a Python int here, and numpy 1.26 promotes
        // int + float32 to float64 for scalars, so this branch finishes in
        // double precision even though the data is float32.
        start_size = 10.0;
        const double denominator = 10.0 + static_cast<double>(mean32);
        start_prob = denominator != 0.0 ? 10.0 / denominator : 10.0;
    }

    const std::vector<double> begin{start_size, start_prob};
    const std::vector<Bound> bounds{
        Bound{kMachineEpsilonStats, std::numeric_limits<double>::infinity()},
        Bound{kMachineEpsilonStats, 1.0}};

    if (factorial_overflowed) {
        // The objective is +inf everywhere, so the forward difference gradient
        // is inf - inf = NaN, the projected gradient norm test compares
        // against NaN and std::max leaves it at zero, and minimise_lbfgsb
        // breaks out at iteration zero with the clamped starting point.
        // Returning it directly is the same answer, asserted bit for bit by
        // the unit test in cpp/tests/test_detect_loops.cpp.
        //
        // This is not an optimisation and is not claimed as one. Measured
        // interleaved on gm12878_chr1.cool, three runs each: 7.12 s of CPU
        // with it against 7.12 s without, because only 3 of the 199 distances
        // reach it. It is here because the alternative is to depend on
        // std::max(0.0, NaN) evaluating to 0.0 for the tool to terminate at
        // the right point, which is an accident of the comparison order in
        // projected_gradient_norm rather than something this code states.
        fit.size = std::min(std::max(start_size, bounds[0].lower), bounds[0].upper);
        fit.prob = std::min(std::max(start_prob, bounds[1].lower), bounds[1].upper);
        fit.status = 0;
        fit.iterations = 0;
        return fit;
    }

    const LbfgsbResult solution =
        minimise_lbfgsb(negative_log_likelihood, begin, bounds, LbfgsbOptions());
    fit.size = solution.x[0];
    fit.prob = solution.x[1];
    fit.status = solution.status;
    fit.iterations = solution.iterations;
    return fit;
}

void bonferroni_in_place(std::vector<double>& pvalues) {
    const double n = static_cast<double>(pvalues.size());
    for (double& p : pvalues) {
        p *= n;
    }
    for (double& p : pvalues) {
        if (!std::isnan(p) && p > 1.0) {
            p = 1.0;
        }
    }
}

}  // namespace hicx::stats
