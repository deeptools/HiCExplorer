// See hicx/cnb.hpp.

#include "hicx/cnb.hpp"

#include <cmath>

#include "hicx/scipy_special.hpp"
#include "special/cephes/gamma.h"

namespace hicx::cnb {

double gammaln(double x) {
    return special::cephes::lgam(x);
}

double xlog1py(double x, double y) {
    if (x == 0.0 && !std::isnan(y)) {
        return 0.0;
    }
    return x * std::log1p(y);
}

double pdf(double x, double r, double p) {
    const double gamma_part = gammaln(r + x) - gammaln(x + 1.0) - gammaln(r);
    return std::exp(gamma_part + (r * std::log(p)) + xlog1py(x, -p));
}

double cdf(double x, double r, double p) {
    return hicx::scipy::betainc(r, x + 1.0, p);
}

}  // namespace hicx::cnb
