// See hicx/scipy_special.hpp.

#include "hicx/scipy_special.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

#include <boost/math/special_functions/beta.hpp>

namespace hicx::scipy {

namespace {

// scipy/special/boost_special_functions.h, SpecialPolicy.
using SpecialPolicy = boost::math::policies::policy<
    boost::math::policies::promote_float<false>,
    boost::math::policies::promote_double<false>,
    boost::math::policies::max_root_iterations<400>>;

}  // namespace

double betainc(double a, double b, double x) {
    // ibeta_wrap
    if (std::isnan(a) || std::isnan(b) || std::isnan(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if ((a <= 0) || (b <= 0) || (x < 0) || (x > 1)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    try {
        return boost::math::ibeta(a, b, x, SpecialPolicy());
    } catch (const std::domain_error&) {
        return std::numeric_limits<double>::quiet_NaN();
    } catch (const std::overflow_error&) {
        return std::numeric_limits<double>::infinity();
    } catch (const std::underflow_error&) {
        return 0.0;
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

}  // namespace hicx::scipy
