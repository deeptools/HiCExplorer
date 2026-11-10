// Contact input of hicDifferentialAnalysis; see diff_contacts_impl.hpp.

#include "diff_contacts_impl.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "hicx/stats_ops.hpp"

namespace hicx::diffc {

std::uint64_t mix64(std::uint64_t z) noexcept {
    z = (z ^ (z >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27U)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31U);
}

CounterRng::CounterRng(std::initializer_list<std::uint64_t> keys) noexcept {
    std::uint64_t state = 0x6a09e667f3bcc909ULL;
    for (std::uint64_t key : keys) {
        state = mix64(state ^ key) + 0x9e3779b97f4a7c15ULL;
    }
    state_ = state;
}

std::uint64_t CounterRng::next() noexcept {
    state_ += 0x9e3779b97f4a7c15ULL;
    return mix64(state_);
}

double CounterRng::uniform() noexcept {
    return static_cast<double>(next() >> 11U) * 0x1.0p-53;
}

std::int64_t binomial(std::int64_t n, double p, CounterRng& rng) {
    if (n <= 0 || p <= 0.0) {
        return 0;
    }
    if (p >= 1.0) {
        return n;
    }
    if (n <= 32) {
        std::int64_t k = 0;
        for (std::int64_t t = 0; t < n; ++t) {
            if (rng.uniform() < p) {
                ++k;
            }
        }
        return k;
    }
    const double q = 1.0 - p;
    std::int64_t mode = static_cast<std::int64_t>(std::floor(static_cast<double>(n + 1) * p));
    mode = std::clamp<std::int64_t>(mode, 0, n);
    const double dn = static_cast<double>(n);
    const double dm = static_cast<double>(mode);
    const double log_pm = hicx::stats::gammaln(dn + 1.0) - hicx::stats::gammaln(dm + 1.0) -
                          hicx::stats::gammaln(dn - dm + 1.0) + dm * std::log(p) +
                          (dn - dm) * std::log1p(-p);
    const double pm = std::exp(log_pm);
    double u = rng.uniform() - pm;
    if (u <= 0.0) {
        return mode;
    }
    const double ratio = p / q;
    double up = pm;
    double down = pm;
    std::int64_t k_up = mode;
    std::int64_t k_down = mode;
    while (k_up < n || k_down > 0) {
        if (k_up < n) {
            up *= static_cast<double>(n - k_up) / static_cast<double>(k_up + 1) * ratio;
            ++k_up;
            u -= up;
            if (u <= 0.0) {
                return k_up;
            }
        }
        if (k_down > 0) {
            down *= static_cast<double>(k_down) / static_cast<double>(n - k_down + 1) / ratio;
            --k_down;
            u -= down;
            if (u <= 0.0) {
                return k_down;
            }
        }
        if (up < 1e-300 && down < 1e-300) {
            break;
        }
    }
    // Only the rounding of the accumulated mass remains.
    return mode;
}

// --------------------------------------------------------------------------

std::pair<std::size_t, std::size_t> ChromosomeData::row_range(std::int64_t i,
                                                              std::int64_t col_first,
                                                              std::int64_t col_last) const {
    const auto begin = col.begin() + row_start[static_cast<std::size_t>(i)];
    const auto end = col.begin() + row_start[static_cast<std::size_t>(i) + 1];
    const auto lo = std::lower_bound(begin, end, static_cast<std::int32_t>(col_first));
    const auto hi = std::lower_bound(lo, end, static_cast<std::int32_t>(col_last));
    return {static_cast<std::size_t>(lo - col.begin()), static_cast<std::size_t>(hi - col.begin())};
}

std::vector<ChromosomeData> load_chromosome(const hicx::CoolFile& file, std::int64_t first,
                                            std::int64_t last, std::int64_t band,
                                            const LoadRequest& request) {
    std::vector<ChromosomeData> out(request.split ? 2 : 1);
    const std::int64_t bins = std::max<std::int64_t>(last - first, 0);
    for (ChromosomeData& data : out) {
        data.first_bin = first;
        data.bins = bins;
        data.coverage.assign(static_cast<std::size_t>(bins), 0.0);
        data.row_start.assign(static_cast<std::size_t>(bins) + 1, 0);
    }
    if (bins == 0) {
        return out;
    }
    std::int64_t previous_row = -1;
    std::int64_t previous_col = -1;
    file.for_each_pixel_chunk(first, last, first, last, [&](const hicx::PixelChunk& chunk) {
        for (std::size_t k = 0; k < chunk.bin1.size(); ++k) {
            std::int64_t b1 = chunk.bin1[k];
            std::int64_t b2 = chunk.bin2[k];
            if (b1 > b2) {
                std::swap(b1, b2);
            }
            const double c = chunk.count[k];
            if (!(c >= 0.0) || c != std::floor(c) || c > 4294967295.0) {
                throw std::runtime_error(
                    "the model works on raw contact counts, but '" + request.path +
                    "' holds the count " + std::to_string(c) +
                    ", which is not a non-negative integer below 2^32 (a balanced or "
                    "normalised matrix?); give the matrix of raw counts");
            }
            const std::int64_t i = b1 - first;
            const std::int64_t j = b2 - first;
            if (i < previous_row || (i == previous_row && j <= previous_col)) {
                throw std::runtime_error("the pixels of '" + request.path +
                                         "' are not sorted by (bin1, bin2), as the cooler "
                                         "format requires");
            }
            previous_row = i;
            previous_col = j;
            const auto whole = static_cast<std::int64_t>(c);
            std::int64_t values[2] = {whole, 0};
            if (request.split) {
                CounterRng rng({request.split_seed, kSplitStream,
                                static_cast<std::uint64_t>(request.index),
                                static_cast<std::uint64_t>(b1), static_cast<std::uint64_t>(b2)});
                values[0] = binomial(whole, 0.5, rng);
                values[1] = whole - values[0];
            }
            for (std::size_t h = 0; h < out.size(); ++h) {
                if (values[h] == 0) {
                    continue;
                }
                ChromosomeData& data = out[h];
                const auto v = static_cast<double>(values[h]);
                if (i != j) {
                    data.coverage[static_cast<std::size_t>(i)] += v;
                    data.coverage[static_cast<std::size_t>(j)] += v;
                }
                if (j - i <= band) {
                    ++data.row_start[static_cast<std::size_t>(i) + 1];
                    data.col.push_back(static_cast<std::int32_t>(j));
                    data.count.push_back(static_cast<std::uint32_t>(values[h]));
                }
            }
        }
    });
    for (ChromosomeData& data : out) {
        for (std::size_t i = 1; i < data.row_start.size(); ++i) {
            data.row_start[i] += data.row_start[i - 1];
        }
        data.col.shrink_to_fit();
        data.count.shrink_to_fit();
    }
    return out;
}

double median(std::vector<double> values) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::sort(values.begin(), values.end());
    const std::size_t n = values.size();
    return n % 2 == 1 ? values[n / 2] : 0.5 * (values[n / 2 - 1] + values[n / 2]);
}

std::vector<char> coverage_outliers(const std::vector<double>& coverage, double lower,
                                    double upper) {
    std::vector<char> invalid(coverage.size(), 0);
    std::vector<double> positive;
    for (double c : coverage) {
        if (c > 0.0) {
            positive.push_back(c);
        }
    }
    if (positive.empty()) {
        std::fill(invalid.begin(), invalid.end(), 1);
        return invalid;
    }
    const double centre = median(positive);
    std::vector<double> deviation(coverage.size());
    for (std::size_t i = 0; i < coverage.size(); ++i) {
        deviation[i] = std::abs(coverage[i] - centre);
    }
    const double mad = median(deviation);
    for (std::size_t i = 0; i < coverage.size(); ++i) {
        if (coverage[i] == 0.0) {
            invalid[i] = 1;
        } else if (mad > 0.0) {
            const double z = 0.6745 * (coverage[i] - centre) / mad;
            invalid[i] = (z < lower || z > upper) ? 1 : 0;
        }
    }
    return invalid;
}

std::vector<double> valid_pairs(const std::vector<char>& invalid, std::int64_t band) {
    const auto n = static_cast<std::int64_t>(invalid.size());
    const std::int64_t top = std::max<std::int64_t>(0, std::min(band, n - 1));
    std::vector<double> pairs(static_cast<std::size_t>(band) + 1, 0.0);
    for (std::int64_t d = 0; d <= top; ++d) {
        std::int64_t count = 0;
        for (std::int64_t i = 0; i + d < n; ++i) {
            count += (invalid[static_cast<std::size_t>(i)] == 0 &&
                      invalid[static_cast<std::size_t>(i + d)] == 0)
                         ? 1
                         : 0;
        }
        pairs[static_cast<std::size_t>(d)] = static_cast<double>(count);
    }
    return pairs;
}

std::vector<double> distance_decay(const ChromosomeData& data, const std::vector<char>& invalid,
                                   std::int64_t band) {
    std::vector<double> sums(static_cast<std::size_t>(band) + 1, 0.0);
    for (std::int64_t i = 0; i < data.bins; ++i) {
        if (invalid[static_cast<std::size_t>(i)] != 0) {
            continue;
        }
        for (std::int64_t k = data.row_start[static_cast<std::size_t>(i)];
             k < data.row_start[static_cast<std::size_t>(i) + 1]; ++k) {
            const std::int64_t j = data.col[static_cast<std::size_t>(k)];
            const std::int64_t d = j - i;
            if (d <= band && invalid[static_cast<std::size_t>(j)] == 0) {
                sums[static_cast<std::size_t>(d)] += data.count[static_cast<std::size_t>(k)];
            }
        }
    }
    const std::vector<double> pairs = valid_pairs(invalid, band);
    for (std::size_t d = 0; d < sums.size(); ++d) {
        sums[d] = pairs[d] > 0.0 ? sums[d] / pairs[d] : 0.0;
    }
    return sums;
}

void thin_pixels(ChromosomeData& data,
                 const std::function<bool(std::int64_t, std::int64_t)>& select, double keep,
                 std::uint64_t seed, std::uint64_t stream) {
    for (std::int64_t i = 0; i < data.bins; ++i) {
        for (std::int64_t k = data.row_start[static_cast<std::size_t>(i)];
             k < data.row_start[static_cast<std::size_t>(i) + 1]; ++k) {
            const std::int64_t j = data.col[static_cast<std::size_t>(k)];
            if (!select(i, j)) {
                continue;
            }
            std::uint32_t& c = data.count[static_cast<std::size_t>(k)];
            CounterRng rng({seed, kPlantStream, stream,
                            static_cast<std::uint64_t>(data.first_bin + i),
                            static_cast<std::uint64_t>(data.first_bin + j)});
            c = static_cast<std::uint32_t>(binomial(static_cast<std::int64_t>(c), keep, rng));
        }
    }
}

std::vector<std::vector<std::string>> read_fields(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot read '" + path + "'");
    }
    std::vector<std::vector<std::string>> rows;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty() || line[0] == '#' || line.rfind("track", 0) == 0 ||
            line.rfind("browser", 0) == 0) {
            continue;
        }
        std::istringstream fields(line);
        std::vector<std::string> row;
        std::string field;
        while (fields >> field) {
            row.push_back(field);
        }
        if (!row.empty()) {
            rows.push_back(std::move(row));
        }
    }
    return rows;
}

}  // namespace hicx::diffc
