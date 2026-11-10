#include "norms.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

namespace hicfilecpp::detail {

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();

}  // namespace

double percentileMath2(std::vector<double> values, double p) {
    const size_t length = values.size();
    if (length == 0) {
        return kNaN;
    }
    if (length == 1) {
        return values[0];
    }
    const auto n = static_cast<double>(length);
    const double pos = p * (n + 1) / 100;
    const double fpos = std::floor(pos);
    const auto intPos = static_cast<int64_t>(fpos);
    const double dif = pos - fpos;
    std::sort(values.begin(), values.end());
    if (pos < 1) {
        return values[0];
    }
    if (pos >= n) {
        return values[length - 1];
    }
    const double lower = values[static_cast<size_t>(intPos - 1)];
    const double upper = values[static_cast<size_t>(intPos)];
    return lower + dif * (upper - lower);
}

double percentileMath3(std::vector<double> values, double p) {
    const size_t length = values.size();
    if (length == 0) {
        return kNaN;
    }
    if (length == 1) {
        return values[0];
    }
    const auto n = static_cast<double>(length);
    const double q = p / 100;
    const double pos = q == 0 ? 0 : (q == 1 ? n : q * (n + 1));
    std::sort(values.begin(), values.end());
    if (pos < 1) {
        return values[0];
    }
    if (pos >= n) {
        return values[length - 1];
    }
    const double fpos = std::floor(pos);
    const auto intPos = static_cast<int64_t>(fpos);
    const double dif = pos - fpos;
    const double lower = values[static_cast<size_t>(intPos - 1)];
    const double upper = values[static_cast<size_t>(intPos)];
    return lower + dif * (upper - lower);
}

// ------------------------------------------------ ExpectedValueCalculation

ExpectedValueCalculation::ExpectedValueCalculation(const std::vector<int64_t>& lengths, int32_t gridSize,
                                                   std::string type)
    : lengths_(lengths), gridSize_(gridSize), type_(std::move(type)) {
    int64_t maxLen = 0;
    for (size_t i = 1; i < lengths_.size(); ++i) {
        maxLen = std::max(maxLen, lengths_[i]);
    }
    numberOfBins_ = maxLen / gridSize_ + 1;
    actualDistances_.assign(static_cast<size_t>(numberOfBins_), 0.0);
}

void ExpectedValueCalculation::addDistance(int32_t chrIndex, int32_t bin1, int32_t bin2, double weight) {
    if (std::isnan(weight)) {
        return;
    }
    if (chrIndex <= 0 || static_cast<size_t>(chrIndex) >= lengths_.size()) {
        return;
    }
    const auto found = chromosomeCounts_.find(chrIndex);
    if (found == chromosomeCounts_.end()) {
        chromosomeCounts_.emplace(chrIndex, weight);
    } else {
        found->second = found->second + weight;
    }
    const int64_t dist = std::abs(static_cast<int64_t>(bin1) - bin2);
    actualDistances_[static_cast<size_t>(dist)] += weight;
}

void ExpectedValueCalculation::computeDensity() {
    int64_t maxNumBins = 0;
    std::vector<double> possibleDistances(static_cast<size_t>(numberOfBins_), 0.0);
    for (const auto& entry : chromosomeCounts_) {
        const int64_t nChrBins = lengths_[static_cast<size_t>(entry.first)] / gridSize_;
        maxNumBins = std::max(maxNumBins, nChrBins);
        for (int64_t i = 0; i < nChrBins; ++i) {
            possibleDistances[static_cast<size_t>(i)] += static_cast<double>(nChrBins - i);
        }
    }
    densityAvg_.assign(static_cast<size_t>(maxNumBins), 0.0);
    const auto& actual = actualDistances_;
    const auto at = [](const std::vector<double>& v, int64_t i) { return v[static_cast<size_t>(i)]; };
    double numSum = actual.empty() ? 0.0 : actual[0];
    double denSum = possibleDistances.empty() ? 0.0 : possibleDistances[0];
    int64_t bound1 = 0;
    int64_t bound2 = 0;
    for (int64_t ii = 0; ii < maxNumBins; ++ii) {
        if (numSum < 400) {
            while (numSum < 400 && bound2 < maxNumBins) {
                bound2++;
                numSum += at(actual, bound2);
                denSum += at(possibleDistances, bound2);
            }
        } else if (numSum >= 400 && bound2 - bound1 > 0) {
            while (bound2 - bound1 > 0 && bound2 < numberOfBins_ && bound1 < numberOfBins_ &&
                   numSum - at(actual, bound1) - at(actual, bound2) >= 400) {
                numSum = numSum - at(actual, bound1) - at(actual, bound2);
                denSum = denSum - at(possibleDistances, bound1) - at(possibleDistances, bound2);
                bound1++;
                bound2--;
            }
        }
        densityAvg_[static_cast<size_t>(ii)] = numSum / denSum;
        if (bound2 + 2 < maxNumBins) {
            numSum += at(actual, bound2 + 1) + at(actual, bound2 + 2);
            denSum += at(possibleDistances, bound2 + 1) + at(possibleDistances, bound2 + 2);
            bound2 += 2;
        } else if (bound2 + 1 < maxNumBins) {
            numSum += at(actual, bound2 + 1);
            denSum += at(possibleDistances, bound2 + 1);
            bound2++;
        }
    }
    chrScaleFactors_.clear();
    for (const auto& entry : chromosomeCounts_) {
        const int64_t nChrBins = lengths_[static_cast<size_t>(entry.first)] / gridSize_;
        double expectedCount = 0;
        for (int64_t n = 0; n < nChrBins; ++n) {
            if (n < maxNumBins) {
                const double v = densityAvg_[static_cast<size_t>(n)];
                expectedCount += static_cast<double>(nChrBins - n) * v;
            }
        }
        chrScaleFactors_[entry.first] = expectedCount / entry.second;
    }
}

// ------------------------------------------------------------------- VC

template <class T>
std::vector<T> computeVC(const std::vector<NormRecord>& records, int64_t size) {
    std::vector<T> rowsums(static_cast<size_t>(size), T(0));
    for (const auto& r : records) {
        rowsums[static_cast<size_t>(r.x)] += r.counts;
        if (r.x != r.y) {
            rowsums[static_cast<size_t>(r.y)] += r.counts;
        }
    }
    return rowsums;
}

template <class T>
double sumFactor(const std::vector<NormRecord>& records, const std::vector<T>& norm) {
    double matrixSum = 0;
    double normSum = 0;
    for (const auto& r : records) {
        const double vx = norm[static_cast<size_t>(r.x)];
        const double vy = norm[static_cast<size_t>(r.y)];
        if (!std::isnan(vx) && !std::isnan(vy) && vx > 0 && vy > 0) {
            if (r.x == r.y) {
                normSum += r.counts / (vx * vy);
                matrixSum += r.counts;
            } else {
                normSum += 2.0f * r.counts / (vx * vy);
                matrixSum += 2.0f * r.counts;
            }
        }
    }
    return std::sqrt(normSum / matrixSum);
}

// ------------------------------------------------------------------- KR

namespace {

std::vector<int32_t> getOffset(const std::vector<NormRecord>& records, int64_t size, double percent, bool math3) {
    std::vector<double> rowSums(static_cast<size_t>(size), 0.0);
    for (const auto& r : records) {
        rowSums[static_cast<size_t>(r.x)] += r.counts;
        if (r.x != r.y) {
            rowSums[static_cast<size_t>(r.y)] += r.counts;
        }
    }
    double thresh = 0;
    if (percent > 0) {
        std::vector<double> positive;
        for (const double sum : rowSums) {
            if (sum != 0) {
                positive.push_back(sum);
            }
        }
        thresh = math3 ? percentileMath3(std::move(positive), percent) : percentileMath2(std::move(positive), percent);
    }
    std::vector<int32_t> offset(rowSums.size());
    int32_t index = 0;
    for (size_t i = 0; i < rowSums.size(); ++i) {
        offset[i] = rowSums[i] <= thresh ? -1 : index++;
    }
    return offset;
}

std::vector<double> krMultiply(const std::vector<NormRecord>& records, const std::vector<int32_t>& offset,
                               const std::vector<double>& vector) {
    std::vector<double> result(vector.size(), 0.0);
    for (const auto& r : records) {
        const int32_t row = offset[static_cast<size_t>(r.x)];
        const int32_t col = offset[static_cast<size_t>(r.y)];
        if (row != -1 && col != -1) {
            result[static_cast<size_t>(row)] += vector[static_cast<size_t>(col)] * r.counts;
            if (row != col) {
                result[static_cast<size_t>(col)] += vector[static_cast<size_t>(row)] * r.counts;
            }
        }
    }
    return result;
}

// computeKRNormVector (the bnewt balancing of Knight and Ruiz)
std::optional<std::vector<double>> computeKRNormVector(const std::vector<NormRecord>& records,
                                                       const std::vector<int32_t>& offset, double tol,
                                                       std::vector<double> x0, double delta) {
    const size_t n = x0.size();
    const std::vector<double> e(n, 1.0);
    const double g = 0.9;
    const double etamax = 0.1;
    double eta = etamax;
    const double rt = std::pow(tol, 2);
    std::vector<double> v = krMultiply(records, offset, x0);
    std::vector<double> rk(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        v[i] = v[i] * x0[i];
        rk[i] = 1 - v[i];
    }
    double rho_km1 = 0;
    for (const double aRk : rk) {
        rho_km1 += aRk * aRk;
    }
    double rout = rho_km1;
    double rold = rout;
    int notChanging = 0;
    while (rout > rt && notChanging < 100) {
        int k = 0;
        std::vector<double> y = e;
        std::vector<double> ynew(n, 0.0);
        std::vector<double> Z(n, 0.0);
        std::vector<double> p(n, 0.0);
        std::vector<double> w(n, 0.0);
        double rho_km2 = rho_km1;
        const double innertol = std::max(std::pow(eta, 2) * rout, rt);
        while (rho_km1 > innertol) {
            k++;
            if (k == 1) {
                rho_km1 = 0;
                for (size_t i = 0; i < n; ++i) {
                    Z[i] = rk[i] / v[i];
                    p[i] = Z[i];
                    rho_km1 += rk[i] * Z[i];
                }
            } else {
                const double beta = rho_km1 / rho_km2;
                for (size_t i = 0; i < n; ++i) {
                    p[i] = Z[i] + beta * p[i];
                }
            }
            std::vector<double> tmp(n);
            for (size_t i = 0; i < n; ++i) {
                tmp[i] = x0[i] * p[i];
            }
            tmp = krMultiply(records, offset, tmp);
            double alpha = 0;
            for (size_t i = 0; i < n; ++i) {
                w[i] = x0[i] * tmp[i] + v[i] * p[i];
                alpha += p[i] * w[i];
            }
            alpha = rho_km1 / alpha;
            double minynew = std::numeric_limits<double>::max();
            for (size_t i = 0; i < n; ++i) {
                ynew[i] = y[i] + alpha * p[i];
                if (ynew[i] < minynew) {
                    minynew = ynew[i];
                }
            }
            if (minynew <= delta) {
                if (delta == 0) {
                    break;
                }
                double gamma = std::numeric_limits<double>::max();
                for (size_t i = 0; i < n; ++i) {
                    if (alpha * p[i] < 0) {
                        if ((delta - y[i]) / (alpha * p[i]) < gamma) {
                            gamma = (delta - y[i]) / (alpha * p[i]);
                        }
                    }
                }
                for (size_t i = 0; i < n; ++i) {
                    y[i] = y[i] + gamma * alpha * p[i];
                }
                break;
            }
            rho_km2 = rho_km1;
            rho_km1 = 0;
            for (size_t i = 0; i < n; ++i) {
                y[i] = ynew[i];
                rk[i] = rk[i] - alpha * w[i];
                Z[i] = rk[i] / v[i];
                rho_km1 += rk[i] * Z[i];
            }
        }
        for (size_t i = 0; i < n; ++i) {
            x0[i] = x0[i] * y[i];
        }
        v = krMultiply(records, offset, x0);
        rho_km1 = 0;
        for (size_t i = 0; i < v.size(); ++i) {
            v[i] = v[i] * x0[i];
            rk[i] = 1 - v[i];
            rho_km1 += rk[i] * rk[i];
        }
        if (std::abs(rho_km1 - rout) < 0.000001 || std::isinf(rho_km1)) {
            notChanging++;
        }
        rout = rho_km1;
        const double rat = rout / rold;
        rold = rout;
        const double rNorm = std::sqrt(rout);
        const double etaO = eta;
        eta = g * rat;
        if (g * std::pow(etaO, 2) > 0.1) {
            eta = std::max(eta, g * std::pow(etaO, 2));
        }
        eta = std::max(std::min(eta, etamax), 0.5 * tol / rNorm);
    }
    if (notChanging >= 100) {
        return std::nullopt;
    }
    return x0;
}

}  // namespace

template <class T>
std::vector<T> computeKR(const std::vector<NormRecord>& records, int64_t size, bool math3) {
    bool recalculate = true;
    std::vector<int32_t> offset = getOffset(records, size, 0, math3);
    std::vector<T> kr;
    int iteration = 1;
    while (recalculate && iteration <= 6) {
        size_t newSize = 0;
        for (const int32_t o : offset) {
            if (o != -1) {
                newSize++;
            }
        }
        auto x0 = computeKRNormVector(records, offset, 0.000001, std::vector<double>(newSize, 1.0), 0.1);
        recalculate = false;
        if (!x0 || iteration == 5) {
            recalculate = true;
            offset = getOffset(records, size, iteration < 5 ? iteration : 10, math3);
        } else {
            kr.assign(static_cast<size_t>(size), T(0));
            size_t krIndex = 0;
            for (const int32_t o : offset) {
                kr[krIndex++] = o == -1 ? std::numeric_limits<T>::quiet_NaN()
                                        : static_cast<T>(1.0 / (*x0)[static_cast<size_t>(o)]);
            }
            const double mySum = sumFactor(records, kr);
            int32_t index = 0;
            for (size_t i = 0; i < kr.size(); ++i) {
                if (static_cast<double>(kr[i]) * mySum < 0.01) {
                    offset[i] = -1;
                    recalculate = true;
                } else if (offset[i] != -1) {
                    offset[i] = index++;
                }
            }
        }
        iteration++;
    }
    if (iteration > 6 && recalculate) {
        kr.assign(static_cast<size_t>(size), std::numeric_limits<T>::quiet_NaN());
    }
    return kr;
}

// ---------------------------------------------------------------- SCALE

namespace {

template <class T>
std::vector<T> rowSumsTimes(const std::vector<NormRecord>& records, const std::vector<T>& vector) {
    std::vector<double> sums(vector.size(), 0.0);
    for (const auto& r : records) {
        double counts = r.counts;
        if (r.x == r.y) {
            counts *= .5;
        }
        sums[static_cast<size_t>(r.x)] += counts * vector[static_cast<size_t>(r.y)];
        sums[static_cast<size_t>(r.y)] += counts * vector[static_cast<size_t>(r.x)];
    }
    return std::vector<T>(sums.begin(), sums.end());
}

// FinalScale.scaleToTargetVector
template <class T>
std::vector<T> scaleToTargetVector(const std::vector<NormRecord>& records, const std::vector<T>& targetVectorInitial) {
    constexpr float tol = .0005f;
    constexpr float percentLowRowSumExcluded = 0.0001f;
    constexpr float dp = percentLowRowSumExcluded / 2;
    constexpr float percentZValsToIgnore = 0;
    constexpr float dp1 = 0;
    constexpr float tolerance = .0005f;
    constexpr int maxIter = 100;
    constexpr int totalIterations = 3 * maxIter;
    constexpr float minErrorThreshold = .02f;
    constexpr float OFFSET = .5f;
    const T nan = std::numeric_limits<T>::quiet_NaN();

    float localPercentLowRowSumExcluded = percentLowRowSumExcluded;
    float localPercentZValsToIgnore = percentZValsToIgnore;
    const size_t k = targetVectorInitial.size();

    std::vector<T> current(k, T(0));
    std::vector<T> row;
    std::vector<T> col;
    std::vector<T> dr(k, T(0));
    std::vector<T> dc(k, T(0));
    std::vector<int> bad(k, 0);
    std::vector<int> bad1(k, 0);
    std::vector<T> s(k, T(0));
    std::vector<double> zz;
    std::vector<double> r0;
    std::vector<T> zTargetVector = targetVectorInitial;
    std::vector<T> calculatedVectorB(k, T(0));
    std::vector<T> one(k, T(1));
    std::vector<double> reportErrorForIteration(static_cast<size_t>(totalIterations + 3), 0.0);
    std::vector<int> numNonZero(k, 0);

    for (size_t p = 0; p < k; ++p) {
        if (std::isnan(zTargetVector[p])) {
            continue;
        }
        if (zTargetVector[p] > 0) {
            zz.push_back(zTargetVector[p]);
        }
    }
    std::sort(zz.begin(), zz.end());
    const auto l = static_cast<int64_t>(zz.size());
    auto zlind = static_cast<int64_t>(std::max(0.0f, static_cast<float>(l) * localPercentZValsToIgnore + OFFSET));
    auto zhind = static_cast<int64_t>(
        std::min(static_cast<double>(l - 1), static_cast<double>(l) * (1.0 - localPercentZValsToIgnore) + OFFSET));
    double zLow = zz.empty() ? kNaN : zz[static_cast<size_t>(zlind)];
    double zHigh = zz.empty() ? kNaN : zz[static_cast<size_t>(zhind)];
    for (size_t p = 0; p < k; ++p) {
        if (zTargetVector[p] > 0 && (zTargetVector[p] < zLow || zTargetVector[p] > zHigh)) {
            zTargetVector[p] = nan;
        }
    }
    for (size_t p = 0; p < k; ++p) {
        if (zTargetVector[p] == 0) {
            one[p] = 0;
        }
    }
    for (const auto& r : records) {
        numNonZero[static_cast<size_t>(r.x)]++;
        if (r.x != r.y) {
            numNonZero[static_cast<size_t>(r.y)]++;
        }
    }
    for (size_t p = 0; p < k; ++p) {
        if (numNonZero[p] > 0) {
            r0.push_back(numNonZero[p]);
        }
    }
    std::sort(r0.begin(), r0.end());
    const auto n0 = static_cast<int64_t>(r0.size());
    auto rlind = static_cast<int64_t>(std::max(0.0f, static_cast<float>(n0) * localPercentLowRowSumExcluded + OFFSET));
    double low = r0.empty() ? kNaN : r0[static_cast<size_t>(std::min(rlind, n0 - 1))];

    for (size_t p = 0; p < k; ++p) {
        if ((numNonZero[p] < low && zTargetVector[p] > 0) || std::isnan(zTargetVector[p])) {
            bad[p] = 1;
            zTargetVector[p] = T(1);
        }
    }
    row = rowSumsTimes(records, one);
    const std::vector<T> rowBackup = row;
    for (size_t p = 0; p < k; ++p) {
        dr[p] = static_cast<T>(1 - bad[p]);
    }
    dc = dr;
    one = dr;
    for (size_t p = 0; p < k; ++p) {
        if (zTargetVector[p] == 0) {
            one[p] = 0;
        }
    }
    for (size_t p = 0; p < k; ++p) {
        bad1[p] = static_cast<int>(1 - one[p]);
    }
    current = dr;

    double ber = 10.0 * (1.0 + tolerance);
    double err = ber;
    int iter = 0;
    int nerr = 0;
    std::vector<double> errors(10000, 0.0);
    int allItersI = 0;

    while ((ber > tolerance || err > 5.0 * tolerance) && iter < maxIter && allItersI < totalIterations &&
           localPercentLowRowSumExcluded <= 0.2 && localPercentZValsToIgnore <= 0.1) {
        iter++;
        allItersI++;
        int fail = 1;
        for (size_t p = 0; p < k; ++p) {
            if (bad1[p] == 1) {
                row[p] = T(1);
            }
        }
        for (size_t p = 0; p < k; ++p) {
            s[p] = zTargetVector[p] / row[p];
        }
        for (size_t p = 0; p < k; ++p) {
            dr[p] *= s[p];
        }
        col = rowSumsTimes(records, dr);
        for (size_t p = 0; p < k; ++p) {
            col[p] *= dc[p];
        }
        for (size_t p = 0; p < k; ++p) {
            if (bad1[p] == 1) {
                col[p] = T(1);
            }
        }
        for (size_t p = 0; p < k; ++p) {
            s[p] = zTargetVector[p] / col[p];
        }
        for (size_t p = 0; p < k; ++p) {
            dc[p] *= s[p];
        }
        row = rowSumsTimes(records, dc);
        for (size_t p = 0; p < k; ++p) {
            row[p] *= dr[p];
        }
        for (size_t p = 0; p < k; ++p) {
            calculatedVectorB[p] = static_cast<T>(std::sqrt(static_cast<double>(dr[p] * dc[p])));
        }
        ber = 0;
        for (size_t p = 0; p < k; ++p) {
            if (bad1[p] == 1) {
                continue;
            }
            const double tempErr = std::abs(calculatedVectorB[p] - current[p]);
            if (tempErr > ber) {
                ber = tempErr;
            }
        }
        reportErrorForIteration[static_cast<size_t>(allItersI - 1)] = ber;
        if (iter % 10 == 0) {
            col = rowSumsTimes(records, calculatedVectorB);
            err = 0;
            for (size_t p = 0; p < k; ++p) {
                if (bad1[p] == 1) {
                    continue;
                }
                const double tempErr = std::abs(col[p] * calculatedVectorB[p] - zTargetVector[p]);
                if (err < tempErr) {
                    err = tempErr;
                }
            }
            errors[static_cast<size_t>(nerr++)] = err;
        }
        current = calculatedVectorB;
        if ((ber < tolerance) &&
            (nerr < 2 || (nerr >= 2 && errors[static_cast<size_t>(nerr - 1)] < 0.5 * errors[static_cast<size_t>(nerr - 2)]))) {
            continue;
        }
        if (iter > 5) {
            for (int q = 1; q <= 5; ++q) {
                if (reportErrorForIteration[static_cast<size_t>(allItersI - q)] * (1.0 + minErrorThreshold) <
                    reportErrorForIteration[static_cast<size_t>(allItersI - q - 1)]) {
                    fail = 0;
                }
            }
            if (nerr >= 2 && errors[static_cast<size_t>(nerr - 1)] > 0.75 * errors[static_cast<size_t>(nerr - 2)]) {
                fail = 1;
            }
            if (iter >= maxIter) {
                fail = 1;
            }
            if (fail == 1) {
                localPercentLowRowSumExcluded += dp;
                localPercentZValsToIgnore += dp1;
                nerr = 0;
                rlind = static_cast<int64_t>(
                    std::max(0.0f, static_cast<float>(n0) * localPercentLowRowSumExcluded + OFFSET));
                low = r0.empty() ? kNaN : r0[static_cast<size_t>(std::min(rlind, n0 - 1))];
                zlind = static_cast<int64_t>(
                    std::max(0.0f, static_cast<float>(l) * localPercentZValsToIgnore + OFFSET));
                zhind = static_cast<int64_t>(std::min(static_cast<double>(l - 1),
                                                      static_cast<double>(l) * (1.0 - localPercentZValsToIgnore) + OFFSET));
                zLow = zz.empty() ? kNaN : zz[static_cast<size_t>(zlind)];
                zHigh = zz.empty() ? kNaN : zz[static_cast<size_t>(zhind)];
                for (size_t p = 0; p < k; ++p) {
                    if (zTargetVector[p] > 0 && (zTargetVector[p] < zLow || zTargetVector[p] > zHigh)) {
                        zTargetVector[p] = nan;
                    }
                }
                for (size_t p = 0; p < k; ++p) {
                    if ((numNonZero[p] < low && zTargetVector[p] > 0) || std::isnan(zTargetVector[p])) {
                        bad[p] = 1;
                        bad1[p] = 1;
                        one[p] = 0;
                        zTargetVector[p] = T(1);
                    }
                }
                ber = 10.0 * (1.0 + tol);
                err = 10.0 * (1.0 + tol);
                if (reportErrorForIteration[static_cast<size_t>(allItersI - 1)] >
                    reportErrorForIteration[static_cast<size_t>(allItersI - 6)]) {
                    for (size_t p = 0; p < k; ++p) {
                        dr[p] = static_cast<T>(1 - bad[p]);
                    }
                    dc = dr;
                    one = dr;
                    current = dr;
                    row = rowBackup;
                } else {
                    for (size_t p = 0; p < k; ++p) {
                        dr[p] *= static_cast<T>(1 - bad[p]);
                    }
                    for (size_t p = 0; p < k; ++p) {
                        dc[p] *= static_cast<T>(1 - bad[p]);
                    }
                }
                iter = 0;
            }
        }
    }
    for (size_t p = 0; p < k; ++p) {
        if (bad[p] == 1) {
            calculatedVectorB[p] = nan;
        }
    }
    return calculatedVectorB;
}

}  // namespace

template <class T>
std::vector<T> computeScale(const std::vector<NormRecord>& records, int64_t size) {
    std::vector<T> vector = scaleToTargetVector(records, std::vector<T>(static_cast<size_t>(size), T(1)));
    // ZeroScale.normalizeVectorByScaleFactor
    for (auto& value : vector) {
        if (value <= 0 || std::isnan(value)) {
            value = std::numeric_limits<T>::quiet_NaN();
        } else {
            value = static_cast<T>(T(1) / value);
        }
    }
    double normalizedSumTotal = 0;
    double sumTotal = 0;
    for (const auto& r : records) {
        const double valX = vector[static_cast<size_t>(r.x)];
        const double valY = vector[static_cast<size_t>(r.y)];
        if (!std::isnan(valX) && !std::isnan(valY)) {
            const double normalizedValue = r.counts / (valX * valY);
            normalizedSumTotal += normalizedValue;
            sumTotal += r.counts;
            if (r.x != r.y) {
                normalizedSumTotal += normalizedValue;
                sumTotal += r.counts;
            }
        }
    }
    const double scaleFactor = std::sqrt(normalizedSumTotal / sumTotal);
    for (auto& value : vector) {
        if (!std::isnan(value)) {
            value = static_cast<T>(scaleFactor * value);
        }
    }
    return vector;
}

template <class T>
void addDistancesFromRecords(ExpectedValueCalculation& ev, int32_t chrIndex, const std::vector<NormRecord>& records,
                             const std::vector<T>& vector) {
    for (const auto& r : records) {
        const T vx = vector[static_cast<size_t>(r.x)];
        const T vy = vector[static_cast<size_t>(r.y)];
        if (vx > 0 && !std::isnan(vx) && vy > 0 && !std::isnan(vy)) {
            const T product = vx * vy;
            const double value = r.counts / product;
            ev.addDistance(chrIndex, r.x, r.y, value);
        }
    }
}

template std::vector<double> computeVC(const std::vector<NormRecord>&, int64_t);
template std::vector<float> computeVC(const std::vector<NormRecord>&, int64_t);
template double sumFactor(const std::vector<NormRecord>&, const std::vector<double>&);
template double sumFactor(const std::vector<NormRecord>&, const std::vector<float>&);
template std::vector<double> computeKR(const std::vector<NormRecord>&, int64_t, bool);
template std::vector<float> computeKR(const std::vector<NormRecord>&, int64_t, bool);
template std::vector<double> computeScale(const std::vector<NormRecord>&, int64_t);
template std::vector<float> computeScale(const std::vector<NormRecord>&, int64_t);
template void addDistancesFromRecords(ExpectedValueCalculation&, int32_t, const std::vector<NormRecord>&,
                                      const std::vector<double>&);
template void addDistancesFromRecords(ExpectedValueCalculation&, int32_t, const std::vector<NormRecord>&,
                                      const std::vector<float>&);

}  // namespace hicfilecpp::detail
