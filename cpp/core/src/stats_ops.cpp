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

}  // namespace

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
