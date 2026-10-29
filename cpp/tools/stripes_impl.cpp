// See stripes_impl.hpp. Each function names the Stripenn 1.1.65.22 function it
// ports (stripenn/getStripe.py unless stated). The order of the arithmetic
// follows the Python so that results are identical, not merely close.

#include "stripes_impl.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstring>
#include <map>
#include <stdexcept>

#include "hicx/parallel.hpp"

namespace hicx::stripes {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr int kFrameDistance = 400;  // Stripenn's table height, in bins

inline bool is_nan(double v) { return v != v; }

// Python's slice normalisation for a step 1 slice on a dimension of length L.
inline void normalise_slice(std::int64_t& start, std::int64_t& stop, std::int64_t length) {
    if (start < 0) {
        start += length;
        if (start < 0) {
            start = 0;
        }
    } else if (start > length) {
        start = length;
    }
    if (stop < 0) {
        stop += length;
        if (stop < 0) {
            stop = 0;
        }
    } else if (stop > length) {
        stop = length;
    }
}

// Python's list index with negative wrap; -1 when out of range.
inline std::int64_t py_index(std::int64_t index, std::int64_t length) {
    if (index < 0) {
        index += length;
    }
    return (index < 0 || index >= length) ? -1 : index;
}

}  // namespace

// ---------------------------------------------------------------------------
// Geometry

NullGeometry null_geometry(std::int64_t chrom_size, std::int64_t resolution) {
    NullGeometry g;
    // itera = min(chrsize / resol / 500, 25), unitsize = floor(chrsize / resol / itera)
    const double raw = std::min(static_cast<double>(chrom_size) / static_cast<double>(resolution) / 500.0, 25.0);
    g.itera = static_cast<int>(raw);
    if (g.itera > 0) {
        g.unitsize = static_cast<int>(std::floor(static_cast<double>(chrom_size) /
                                                 static_cast<double>(resolution) / static_cast<double>(g.itera)));
    }
    return g;
}

std::int64_t required_depth(std::int64_t resolution, int max_width) {
    const std::int64_t background = static_cast<std::int64_t>(50000.0 / static_cast<double>(resolution));
    // 399 bins of the stripe table, a further 400 for the offset of the
    // null distribution's stale unit index (see compute_null), the windows.
    return 2 * kFrameDistance + 2 * background + max_width + 16;
}

// ---------------------------------------------------------------------------
// BandMatrix

BandMatrix::BandMatrix(std::int64_t bins, std::int64_t depth, std::int64_t resolution,
                       std::int64_t chrom_size)
    : bins_(bins), depth_(depth), resolution_(resolution), chrom_size_(chrom_size) {
    values_.assign(static_cast<std::size_t>(bins) * static_cast<std::size_t>(depth + 1), 0.0);
    geometry_ = null_geometry(chrom_size, resolution);
    unit_flags_.assign(static_cast<std::size_t>(bins), 0);
    window_flags_.assign(static_cast<std::size_t>(bins), 0);
    unit_row_end_.assign(static_cast<std::size_t>(std::max(geometry_.itera, 0)), 0);
    for (int it = 0; it < geometry_.itera; ++it) {
        std::int64_t end1 = static_cast<std::int64_t>(geometry_.unitsize) * resolution * (it + 1);
        if (end1 > chrom_size - kFrameDistance * resolution) {
            end1 = chrom_size - kFrameDistance * resolution;
        }
        unit_row_end_[static_cast<std::size_t>(it)] = ceil_div(end1, resolution);
    }
}

void BandMatrix::add(std::int64_t a, std::int64_t b, double value) {
    const std::int64_t d = b - a;
    if (d <= depth_) {
        values_[static_cast<std::size_t>(a) * static_cast<std::size_t>(depth_ + 1) + static_cast<std::size_t>(d)] += value;
    }
    if (geometry_.itera <= 0 || geometry_.unitsize <= 0) {
        return;
    }
    const std::int64_t unit = geometry_.unitsize;
    const std::int64_t itera = geometry_.itera;
    const std::int64_t ua = a / unit;
    const std::int64_t ub = b / unit;
    if (ua < itera && ua == ub) {
        unit_flags_[static_cast<std::size_t>(a)] = 1;
        unit_flags_[static_cast<std::size_t>(b)] = 1;
    }
    // The wide window of row a (unit ua): columns [max(unit*ua - 400, 0), min(unit*(ua+1) + 400, ...)).
    auto in_window = [&](std::int64_t row_unit, std::int64_t column) {
        const std::int64_t c0 = std::max<std::int64_t>(unit * row_unit - kFrameDistance, 0);
        std::int64_t end2 = unit * resolution_ * (row_unit + 1) + kFrameDistance * resolution_;
        if (end2 > chrom_size_) {
            end2 = chrom_size_ - 1;
        }
        const std::int64_t c1 = ceil_div(end2, resolution_);
        return column >= c0 && column < c1;
    };
    if (ua < itera && a < unit_row_end_[static_cast<std::size_t>(ua)] && in_window(ua, b)) {
        window_flags_[static_cast<std::size_t>(a)] = 1;
    }
    if (b != a && ub < itera && b < unit_row_end_[static_cast<std::size_t>(ub)] && in_window(ub, a)) {
        window_flags_[static_cast<std::size_t>(b)] = 1;
    }
}

double BandMatrix::get(std::int64_t row, std::int64_t column) const {
    if (row > column) {
        std::swap(row, column);
    }
    if (bad(row) || bad(column)) {
        return kNaN;
    }
    const std::int64_t d = column - row;
    if (d > depth_) {
        throw std::runtime_error("hicDetectStripes: internal error, matrix entry " + std::to_string(row) + "," +
                                 std::to_string(column) + " is outside the stored band of " +
                                 std::to_string(depth_) + " bins");
    }
    return values_[static_cast<std::size_t>(row) * static_cast<std::size_t>(depth_ + 1) + static_cast<std::size_t>(d)];
}

std::vector<double> BandMatrix::dense(std::int64_t row0, std::int64_t row1, std::int64_t col0,
                                      std::int64_t col1) const {
    const std::int64_t rows = std::max<std::int64_t>(row1 - row0, 0);
    const std::int64_t cols = std::max<std::int64_t>(col1 - col0, 0);
    std::vector<double> out(static_cast<std::size_t>(rows * cols));
    for (std::int64_t r = 0; r < rows; ++r) {
        for (std::int64_t c = 0; c < cols; ++c) {
            out[static_cast<std::size_t>(r * cols + c)] = get(row0 + r, col0 + c);
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// getQuantile_original

QuantileSelector::QuantileSelector(std::vector<double> quantiles)
    : quantiles_(std::move(quantiles)), histogram_(32768, 0), bucket_slot_(32768, -1) {}

namespace {
inline std::size_t bucket_of(double value) {
    std::uint64_t bits;
    std::memcpy(&bits, &value, sizeof bits);
    return static_cast<std::size_t>(bits >> 48);
}
}  // namespace

void QuantileSelector::count(double value, bool diagonal) {
    if (!(value > 0.0)) {
        return;
    }
    const std::uint64_t multiplicity = diagonal ? 1 : 2;
    histogram_[bucket_of(value)] += multiplicity;
    total_ += multiplicity;
}

void QuantileSelector::prepare() {
    ranks_.clear();
    for (const double q : quantiles_) {
        ranks_.push_back(np::quantile_ranks(static_cast<std::int64_t>(total_), q));
    }
    std::vector<std::uint64_t> before(32768, 0);
    std::uint64_t running = 0;
    for (std::size_t k = 0; k < 32768; ++k) {
        before[k] = running;
        running += histogram_[k];
    }
    auto locate = [&](std::int64_t rank) {
        std::size_t lo = 0;
        for (std::size_t k = 0; k < 32768; ++k) {
            if (histogram_[k] != 0 && before[k] <= static_cast<std::uint64_t>(rank) &&
                static_cast<std::uint64_t>(rank) < before[k] + histogram_[k]) {
                lo = k;
                break;
            }
        }
        return lo;
    };
    bucket_values_.clear();
    bucket_start_.clear();
    std::fill(bucket_slot_.begin(), bucket_slot_.end(), -1);
    for (const np::QuantileRanks& r : ranks_) {
        for (const std::int64_t rank : {r.previous, r.next}) {
            const std::size_t b = locate(rank);
            if (bucket_slot_[b] < 0) {
                bucket_slot_[b] = static_cast<std::int64_t>(bucket_values_.size());
                bucket_values_.emplace_back();
                bucket_start_.push_back(before[b]);
            }
        }
    }
}

void QuantileSelector::gather(double value, bool diagonal) {
    if (!(value > 0.0)) {
        return;
    }
    const std::int64_t slot = bucket_slot_[bucket_of(value)];
    if (slot < 0) {
        return;
    }
    bucket_values_[static_cast<std::size_t>(slot)].push_back(value);
    if (!diagonal) {
        bucket_values_[static_cast<std::size_t>(slot)].push_back(value);
    }
}

std::vector<double> QuantileSelector::result() {
    for (auto& values : bucket_values_) {
        std::sort(values.begin(), values.end());
    }
    std::vector<double> out;
    if (total_ == 0) {
        out.assign(quantiles_.size(), kNaN);
        return out;
    }
    auto order_statistic = [&](std::int64_t rank) {
        // find the bucket holding this rank
        for (std::size_t k = 0; k < 32768; ++k) {
            const std::int64_t slot = bucket_slot_[k];
            if (slot < 0) {
                continue;
            }
            const std::uint64_t start = bucket_start_[static_cast<std::size_t>(slot)];
            const auto& values = bucket_values_[static_cast<std::size_t>(slot)];
            if (static_cast<std::uint64_t>(rank) >= start && static_cast<std::uint64_t>(rank) < start + values.size()) {
                return values[static_cast<std::size_t>(static_cast<std::uint64_t>(rank) - start)];
            }
        }
        throw std::runtime_error("hicDetectStripes: quantile rank not found, the two passes saw different pixels");
    };
    for (const np::QuantileRanks& r : ranks_) {
        out.push_back(np::quantile_lerp(order_statistic(r.previous), order_statistic(r.next), r.gamma));
    }
    return out;
}

// ---------------------------------------------------------------------------
// mpmean

std::vector<double> expected_values(const BandMatrix& band) {
    const std::int64_t n = band.bins();
    const std::int64_t nframes = ceil_div(n, kFrameDistance);
    std::vector<std::vector<double>> pixels(static_cast<std::size_t>(nframes), std::vector<double>(kFrameDistance, 0.0));
    std::vector<std::vector<std::int64_t>> counts(static_cast<std::size_t>(nframes),
                                                  std::vector<std::int64_t>(kFrameDistance, 0));
    for (std::int64_t f = 0; f < nframes; ++f) {
        const std::int64_t r0 = f * kFrameDistance;
        const std::int64_t r1 = std::min((f + 1) * kFrameDistance, n);
        const std::int64_t cols_end = std::min((f + 2) * kFrameDistance, n);
        const std::int64_t cfm_cols = cols_end - r0;
        for (std::int64_t i = 0; i < r1 - r0; ++i) {
            for (std::int64_t j = 0; j < std::min<std::int64_t>(kFrameDistance, cfm_cols); ++j) {
                if (i + j >= cfm_cols) {
                    continue;
                }
                double count = band.get(r0 + i, r0 + i + j);
                if (is_nan(count)) {
                    count = 0.0;
                }
                pixels[static_cast<std::size_t>(f)][static_cast<std::size_t>(j)] += count;
                counts[static_cast<std::size_t>(f)][static_cast<std::size_t>(j)] += 1;
            }
        }
    }
    std::vector<double> means(kFrameDistance);
    for (int j = 0; j < kFrameDistance; ++j) {
        double sum = 0.0;
        std::int64_t count = 0;
        for (std::int64_t f = 0; f < nframes; ++f) {
            sum += pixels[static_cast<std::size_t>(f)][static_cast<std::size_t>(j)];
            count += counts[static_cast<std::size_t>(f)][static_cast<std::size_t>(j)];
        }
        means[static_cast<std::size_t>(j)] = count == 0 ? kNaN : sum / static_cast<double>(count);
    }
    return means;
}

// ---------------------------------------------------------------------------
// nulldist

std::int64_t available_columns(const BandMatrix& band) {
    const NullGeometry g = null_geometry(band.chrom_size(), band.resolution());
    std::int64_t pool = 0;
    for (int it = 0; it < g.itera; ++it) {
        for (std::int64_t r = static_cast<std::int64_t>(g.unitsize) * it;
             r < static_cast<std::int64_t>(g.unitsize) * (it + 1); ++r) {
            if (band.unit_flags()[static_cast<std::size_t>(r)] != 0) {
                ++pool;
            }
        }
    }
    return pool;
}

std::vector<std::int64_t> allocate_samples(const std::vector<ChromInfo>& chroms,
                                           const std::vector<std::int64_t>& available) {
    const std::size_t count = chroms.size();
    std::vector<std::int64_t> result(count, -1);
    if (count == 0) {
        return result;
    }
    // First stage, from the chromosome sizes.
    double total_size = 0.0;
    for (const ChromInfo& c : chroms) {
        total_size += static_cast<double>(c.size);
    }
    std::int64_t sum_size = 0;
    for (const ChromInfo& c : chroms) {
        sum_size += c.size;
    }
    std::vector<std::uint64_t> first(count);
    std::uint64_t first_sum = 0;
    for (std::size_t i = 0; i < count; ++i) {
        const double share = (static_cast<double>(chroms[i].size) / static_cast<double>(sum_size)) * 1000.0;
        first[i] = static_cast<std::uint64_t>(share);
        first_sum += first[i];
    }
    (void)total_size;
    std::vector<std::size_t> names2;  // indices into chroms
    for (std::size_t i = 0; i < count; ++i) {
        if (first[i] != 0) {
            names2.push_back(i);
        }
    }
    // Second stage, from the available columns of those chromosomes.
    std::int64_t sum_available = 0;
    for (const std::size_t i : names2) {
        sum_available += available[i];
    }
    std::vector<std::uint64_t> second(names2.size());
    std::uint64_t second_sum = 0;
    for (std::size_t k = 0; k < names2.size(); ++k) {
        const double share = (static_cast<double>(available[names2[k]]) / static_cast<double>(sum_available)) * 1000.0;
        second[k] = share >= 0.0 ? static_cast<std::uint64_t>(share) : 0;
        second_sum += second[k];
    }
    const std::uint64_t dif = 1000 - second_sum;
    std::vector<std::size_t> names3;  // chromnames2 after the second notzero filter
    std::vector<std::size_t> names3_position;
    for (std::size_t k = 0; k < names2.size(); ++k) {
        if (second[k] != 0) {
            names3.push_back(names2[k]);
        }
    }
    if (!second.empty()) {
        second[0] += dif;  // samplesize[0] += dif, on the unfiltered array
    }
    // main_null_calc: c = chromnames2.index(chr), ss = samplesize[c]
    for (std::size_t p = 0; p < names3.size(); ++p) {
        result[names3[p]] = p < second.size() ? static_cast<std::int64_t>(second[p]) : 0;
    }
    (void)first_sum;
    return result;
}

namespace {

struct NullSettings {
    int background_up = 0;
    int background_down = 0;
    int background_size = 0;
};

NullSettings null_settings(std::int64_t resolution) {
    NullSettings s;
    const double background_size = 50000.0 / static_cast<double>(resolution);
    const double up = std::floor(background_size / 2.0);
    const double down = background_size - up;
    s.background_up = static_cast<int>(up);
    s.background_down = static_cast<int>(down);
    s.background_size = static_cast<int>(background_size);
    return s;
}

}  // namespace

NullColumns compute_null(const BandMatrix& band, std::int64_t ss, np::PyRandom& rng) {
    const std::int64_t res = band.resolution();
    const std::int64_t chrsize = band.chrom_size();
    const NullSettings cfg = null_settings(res);
    const NullGeometry geometry = null_geometry(chrsize, res);
    const std::int64_t itera = geometry.itera;
    const std::int64_t unitsize = geometry.unitsize;

    std::vector<std::vector<double>> cols_lu;
    std::vector<std::vector<double>> cols_ru;
    std::vector<std::vector<double>> cols_ld;
    std::vector<std::vector<double>> cols_rd;
    std::vector<std::int64_t> n_pool;

    struct UnitWindow {
        std::int64_t rows0 = 0;
        std::int64_t nrow = 0;
        std::int64_t cols0 = 0;
        std::int64_t ncol = 0;
    };
    auto window_of = [&](std::int64_t unit) {
        UnitWindow w;
        const std::int64_t t_start1 = unitsize * res * unit + 1;
        std::int64_t t_start0 = t_start1 - 400 * res;
        std::int64_t t_end1 = unitsize * res * (unit + 1);
        std::int64_t t_end2 = unitsize * res * (unit + 1) + 400 * res;
        if (t_end2 > chrsize) {
            t_end2 = chrsize - 1;
        }
        if (t_end1 > chrsize - 400 * res) {
            t_end1 = chrsize - 400 * res;
        }
        if (t_start0 <= 1) {
            t_start0 = 1;
        }
        if (t_end1 < t_start1 || t_end2 < t_start0) {
            throw std::runtime_error("hicDetectStripes: chromosome too short for the null distribution");
        }
        w.rows0 = floor_div(t_start1, res);
        w.nrow = ceil_div(t_end1, res) - w.rows0;
        w.cols0 = floor_div(t_start0, res);
        w.ncol = ceil_div(t_end2, res) - w.cols0;
        return w;
    };

    // Mean of mat[r0:r1, c0:c1] as numpy computes it: NaN when empty, the
    // values of masked bins read as zero.
    std::vector<double> buffer;
    auto block_mean = [&](const UnitWindow& w, std::int64_t r0, std::int64_t r1, std::int64_t c0,
                          std::int64_t c1) {
        normalise_slice(r0, r1, w.nrow);
        normalise_slice(c0, c1, w.ncol);
        if (r0 >= r1 || c0 >= c1) {
            return kNaN;
        }
        const std::size_t rows = static_cast<std::size_t>(r1 - r0);
        const std::size_t cols = static_cast<std::size_t>(c1 - c0);
        buffer.resize(rows * cols);
        for (std::size_t r = 0; r < rows; ++r) {
            for (std::size_t c = 0; c < cols; ++c) {
                buffer[r * cols + c] = band.get_zero(w.rows0 + r0 + static_cast<std::int64_t>(r),
                                                     w.cols0 + c0 + static_cast<std::int64_t>(c));
            }
        }
        const bool contiguous = rows == 1 || (c0 == 0 && c1 == w.ncol);
        return np::mean_compact(buffer.data(), rows, cols, contiguous);
    };

    auto draw_columns = [&](const UnitWindow& w, const std::vector<std::int64_t>& pool, std::int64_t k,
                            bool offset) {
        std::vector<std::int64_t> picks(static_cast<std::size_t>(k));
        for (std::int64_t i = 0; i < k; ++i) {
            picks[static_cast<std::size_t>(i)] = pool[rng.choice_index(pool.size())];
        }
        for (std::int64_t i = 0; i < k; ++i) {
            const std::int64_t x = picks[static_cast<std::size_t>(i)];
            std::vector<double> lu(kFrameDistance);
            std::vector<double> ru(kFrameDistance);
            std::vector<double> ld(kFrameDistance);
            std::vector<double> rd(kFrameDistance);
            for (int j = 0; j < kFrameDistance; ++j) {
                std::int64_t y_down = x + j;
                std::int64_t y_up = x - j;
                if (offset) {
                    y_down += 400;
                    y_up += 400;
                }
                const std::int64_t bu = cfg.background_up;
                const std::int64_t bd = cfg.background_down;
                const std::int64_t bs = cfg.background_size;
                const double left_up = block_mean(w, x - bu - bs, x - bu, y_up - bu, y_up + bd);
                const double center_up = block_mean(w, x - bu, x + bd, y_up - bu, y_up + bd);
                const double right_up = block_mean(w, x + bd, x + bd + bs, y_up - bu, y_up + bd);
                const double left_down = block_mean(w, x - bu - bs, x - bu, y_down - bu, y_down + bd);
                const double center_down = block_mean(w, x - bu, x + bd, y_down - bu, y_down + bd);
                const double right_down = block_mean(w, x + bd, x + bd + bs, y_down - bu, y_down + bd);
                lu[static_cast<std::size_t>(j)] = center_up - left_up;
                ru[static_cast<std::size_t>(j)] = center_up - right_up;
                ld[static_cast<std::size_t>(j)] = center_down - left_down;
                rd[static_cast<std::size_t>(j)] = center_down - right_down;
            }
            cols_lu.push_back(std::move(lu));
            cols_ru.push_back(std::move(ru));
            cols_ld.push_back(std::move(ld));
            cols_rd.push_back(std::move(rd));
        }
    };

    for (std::int64_t it = 0; it < itera; ++it) {
        const std::int64_t sss = static_cast<std::int64_t>(static_cast<double>(ss) / static_cast<double>(itera));
        const UnitWindow w = window_of(it);
        std::vector<std::int64_t> pool;
        for (std::int64_t x = 0; x < w.nrow; ++x) {
            if (band.window_flags()[static_cast<std::size_t>(w.rows0 + x)] == 0) {
                continue;
            }
            if (x > 20 && x < unitsize - 20) {
                if (it == 0 && !(x > 410 && x < w.ncol)) {
                    continue;
                }
                pool.push_back(x);
            }
        }
        n_pool.push_back(static_cast<std::int64_t>(pool.size()));
        if (pool.empty()) {
            continue;
        }
        const std::int64_t k = static_cast<std::int64_t>(pool.size()) < sss ? static_cast<std::int64_t>(pool.size()) : sss;
        draw_columns(w, pool, k, it > 0);
    }

    const std::int64_t depl = ss - static_cast<std::int64_t>(cols_lu.size());
    if (depl > 0 && !n_pool.empty()) {
        const std::int64_t rich = static_cast<std::int64_t>(
            std::max_element(n_pool.begin(), n_pool.end()) - n_pool.begin());
        const UnitWindow w = window_of(rich);
        std::vector<std::int64_t> pool;
        for (std::int64_t x = 0; x < w.nrow; ++x) {
            if (band.window_flags()[static_cast<std::size_t>(w.rows0 + x)] != 0 && x > 20 && x < unitsize - 20) {
                pool.push_back(x);
            }
        }
        if (!pool.empty()) {
            // `it` here is the loop variable left over from the loop above,
            // not `rich`, as in the Python.
            draw_columns(w, pool, depl, itera - 1 > 0);
        }
    }

    NullColumns out;
    out.columns = static_cast<std::int64_t>(cols_lu.size());
    const std::size_t n = static_cast<std::size_t>(out.columns);
    out.leftup.assign(kFrameDistance * n, 0.0);
    out.rightup.assign(kFrameDistance * n, 0.0);
    out.leftdown.assign(kFrameDistance * n, 0.0);
    out.rightdown.assign(kFrameDistance * n, 0.0);
    for (std::size_t k = 0; k < n; ++k) {
        for (int d = 0; d < kFrameDistance; ++d) {
            out.leftup[static_cast<std::size_t>(d) * n + k] = cols_lu[k][static_cast<std::size_t>(d)];
            out.rightup[static_cast<std::size_t>(d) * n + k] = cols_ru[k][static_cast<std::size_t>(d)];
            out.leftdown[static_cast<std::size_t>(d) * n + k] = cols_ld[k][static_cast<std::size_t>(d)];
            out.rightdown[static_cast<std::size_t>(d) * n + k] = cols_rd[k][static_cast<std::size_t>(d)];
        }
    }
    return out;
}

NullTables join_null(const std::vector<NullColumns>& parts) {
    NullTables t;
    for (const NullColumns& p : parts) {
        t.columns += p.columns;
    }
    const std::size_t n = static_cast<std::size_t>(t.columns);
    t.leftup.assign(kFrameDistance * n, 0.0);
    t.rightup.assign(kFrameDistance * n, 0.0);
    t.leftdown.assign(kFrameDistance * n, 0.0);
    t.rightdown.assign(kFrameDistance * n, 0.0);
    std::size_t offset = 0;
    for (const NullColumns& p : parts) {
        const std::size_t pn = static_cast<std::size_t>(p.columns);
        for (int d = 0; d < kFrameDistance; ++d) {
            const std::size_t row = static_cast<std::size_t>(d);
            std::copy_n(p.leftup.begin() + static_cast<std::ptrdiff_t>(row * pn), pn, t.leftup.begin() + static_cast<std::ptrdiff_t>(row * n + offset));
            std::copy_n(p.rightup.begin() + static_cast<std::ptrdiff_t>(row * pn), pn, t.rightup.begin() + static_cast<std::ptrdiff_t>(row * n + offset));
            std::copy_n(p.leftdown.begin() + static_cast<std::ptrdiff_t>(row * pn), pn, t.leftdown.begin() + static_cast<std::ptrdiff_t>(row * n + offset));
            std::copy_n(p.rightdown.begin() + static_cast<std::ptrdiff_t>(row * pn), pn, t.rightdown.begin() + static_cast<std::ptrdiff_t>(row * n + offset));
        }
        offset += pn;
    }
    t.build_sorted();
    return t;
}

void NullTables::build_sorted() {
    const std::size_t n = static_cast<std::size_t>(columns);
    auto build = [&](const std::vector<double>& table, std::vector<std::vector<double>>& sorted) {
        sorted.assign(kFrameDistance, {});
        for (int d = 0; d < kFrameDistance; ++d) {
            auto& row = sorted[static_cast<std::size_t>(d)];
            for (std::size_t k = 0; k < n; ++k) {
                const double v = table[static_cast<std::size_t>(d) * n + k];
                if (!is_nan(v)) {
                    row.push_back(v);
                }
            }
            std::sort(row.begin(), row.end());
        }
    };
    build(leftup, sorted_leftup);
    build(rightup, sorted_rightup);
    build(leftdown, sorted_leftdown);
    build(rightdown, sorted_rightdown);
}

// ---------------------------------------------------------------------------
// StripeSearch

namespace {

// ImageProcessing.verticalLine(edges, L=60, H=120)
std::vector<std::uint8_t> vertical_line(const std::vector<std::uint8_t>& edges, int S) {
    static const int gx[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    static const int gy[3][3] = {{1, 2, 1}, {0, 0, 0}, {-1, -2, -1}};
    const std::vector<int> fx = np::convolve2d_same_3x3(edges, S, S, gx);
    const std::vector<int> fy = np::convolve2d_same_3x3(edges, S, S, gy);
    std::vector<std::uint8_t> out(static_cast<std::size_t>(S) * static_cast<std::size_t>(S), 0);
    constexpr double kPi = 3.14159265358979323846;
    for (int r = 0; r < S; ++r) {
        for (int c = 0; c < S; ++c) {
            const std::size_t at = static_cast<std::size_t>(r) * S + c;
            double ori = std::atan2(static_cast<double>(fx[at]), static_cast<double>(fy[at]));
            ori = ori * 180 / kPi;
            if (ori < 0) {
                ori += 360;
            }
            if (ori > 60 && ori < 120) {
                int y = c - 1;
                if (y < 0) {
                    y += S;  // numpy wraps the index -1 to the last column
                }
                out[static_cast<std::size_t>(r) * S + y] = 1;
            }
        }
    }
    return out;
}

// ImageProcessing.block(mat, c)
void block_scan(const std::vector<std::uint8_t>& vert, int S, int c, int& length, int& end) {
    const int lo = std::max(0, c - 1);
    const int hi = std::min(S, c + 2);
    int count = 0;
    int maximum = 0;
    int j = 0;
    int end_row = 0;
    int buffer = 0;
    for (int i = 0; i < S; ++i) {
        bool one = false;
        for (int k = lo; k < hi; ++k) {
            if (vert[static_cast<std::size_t>(i) * S + k] == 1) {
                one = true;
                break;
            }
        }
        if (one) {
            ++count;
            j = i;
        } else if (buffer < 5) {
            ++buffer;
        } else {
            if (count > maximum) {
                maximum = count;
                end_row = j;
            }
            count = 0;
            buffer = 0;
        }
    }
    if (count > maximum) {
        maximum = count;
        end_row = j;
    }
    length = maximum;
    end = end_row;
    if (end < c) {
        end = end - length + 1;
    }
}

}  // namespace

std::vector<RawStripe> stripe_search(const std::vector<double>& submat, int S, double M,
                                     const Parameters& params, FrameTrace* trace) {
    std::vector<RawStripe> found;
    const std::size_t pixels = static_cast<std::size_t>(S) * static_cast<std::size_t>(S);

    // blue = 255 * (M - submat) / M, clipped below at 0, green = blue; the
    // image channels are red 1, green blue / 255, clipped to [0, 1].
    std::vector<double> green(pixels);
    for (std::size_t i = 0; i < pixels; ++i) {
        double blue = 255.0 * (M - submat[i]);
        blue = blue / M;
        if (blue < 0) {
            blue = 0;
        }
        double g = blue / 255.0;
        g = std::min(std::max(g, 0.0), 1.0);
        green[i] = g;
    }
    if (trace != nullptr) {
        trace->size = S;
        trace->gray.assign(6 * pixels, 0.0F);
        trace->edges.assign(6 * pixels, 0);
        trace->vert.assign(6 * pixels, 0);
    }

    const double delta = (0.5 + 0.1) - 0.5;  // numpy.arange(0.5, 1.01, 0.1)
    // The red channel is 1 everywhere; its blur is the same constant at every pixel.
    std::vector<double> ones(9, 1.0);
    const double red_blur = np::box_filter(ones, 3, 3, params.bfilter)[0];
    for (int bi = 0; bi < 6; ++bi) {
        const double b = 0.5 + bi * delta;
        // imBrightness3D(In=([0,0,0],[1,b,b]), Out=([0,0,0],[1,1,1]))
        const double k = (1.0 - 0.0) / (b - 0.0);
        std::vector<double> adj(pixels);
        for (std::size_t i = 0; i < pixels; ++i) {
            const double x = green[i];
            double out;
            if (x <= 0.0) {
                out = 0.0;
            } else if (x > b) {
                out = 1.0;
            } else {
                out = k * (x - 0.0) + 0.0;
            }
            adj[i] = out;
        }
        std::vector<double> blur = np::box_filter(adj, S, S, params.bfilter);
        np::ImageF gray(S, S);
        {
            const double red_clipped = std::min(std::max(red_blur, 0.0), 1.0);
            const float r32 = static_cast<float>(red_clipped);
            for (std::size_t i = 0; i < pixels; ++i) {
                const double gb = std::min(std::max(blur[i], 0.0), 1.0);
                const float g32 = static_cast<float>(gb);
                gray.v[i] = np::rgb_to_gray(r32, g32, g32);
            }
        }
        const std::vector<std::uint8_t> edges = np::canny(gray, params.canny);
        const std::vector<std::uint8_t> vert = vertical_line(edges, S);
        if (trace != nullptr) {
            std::copy(gray.v.begin(), gray.v.end(), trace->gray.begin() + static_cast<std::ptrdiff_t>(bi * pixels));
            std::copy(edges.begin(), edges.end(), trace->edges.begin() + static_cast<std::ptrdiff_t>(bi * pixels));
            std::copy(vert.begin(), vert.end(), trace->vert.begin() + static_cast<std::ptrdiff_t>(bi * pixels));
        }

        std::vector<int> test_column;
        std::vector<int> end_points;
        std::vector<int> start_points;
        std::vector<int> updown;
        for (int c = 0; c < S; ++c) {
            int line_length = 0;
            int end = 0;
            block_scan(vert, S, c, line_length, end);
            const int above = std::min(c, end);
            const int bottom = std::max(c, end);
            bool any = false;
            for (int r = above; r <= bottom; ++r) {
                if (vert[static_cast<std::size_t>(r) * S + c] != 0) {
                    any = true;
                    break;
                }
            }
            if (line_length > params.min_length && any) {
                test_column.push_back(c);
                end_points.push_back(end);
                start_points.push_back(c);
                updown.push_back(end > c ? 2 : 1);
            }
        }

        struct Pair {
            int first;
            int second;
        };
        std::vector<Pair> pairs;
        std::vector<int> min_vec;
        std::vector<int> max_vec;
        for (int ud = 1; ud <= 2; ++ud) {
            std::vector<std::uint8_t> testmat(pixels, 0);
            std::vector<std::size_t> udidx;
            for (std::size_t i = 0; i < updown.size(); ++i) {
                if (updown[i] == ud) {
                    udidx.push_back(i);
                }
            }
            for (const std::size_t c : udidx) {
                int st = test_column[c];
                int en = end_points[c];
                if (ud == 1) {
                    std::swap(st, en);
                }
                // testmat[st:en, test_column[c]] = 1
                for (int r = std::max(st, 0); r < std::min(en, S); ++r) {
                    testmat[static_cast<std::size_t>(r) * S + test_column[c]] = 1;
                }
            }
            // line refinement
            for (int r = 0; r < S; ++r) {
                std::uint8_t* vec = testmat.data() + static_cast<std::size_t>(r) * S;
                std::vector<int> st;
                std::vector<int> en;
                for (int i = 0; i + 1 < S; ++i) {
                    if (vec[i + 1] > vec[i]) {
                        st.push_back(i + 1);
                    }
                    if (vec[i + 1] < vec[i]) {
                        en.push_back(i);
                    }
                }
                if (vec[0] == 1) {
                    st.insert(st.begin(), 0);
                }
                if (vec[S - 1] == 1) {
                    en.push_back(S - 1);
                }
                const std::size_t n_lines = st.size();
                for (std::size_t l = 0; l < n_lines; ++l) {
                    int sum = 0;
                    for (int c = st[l]; c <= en[l] && c < S; ++c) {
                        sum += edges[static_cast<std::size_t>(r) * S + c];
                    }
                    if (sum > 0) {
                        for (int c = st[l]; c < en[l]; ++c) {
                            vec[c] = vert[static_cast<std::size_t>(r) * S + c];
                        }
                    } else {
                        const int med = static_cast<int>(np::round_half_even(static_cast<double>(st[l] + en[l]) / 2.0));
                        for (int c = st[l]; c < en[l]; ++c) {
                            vec[c] = 0;
                        }
                        vec[med] = 1;
                    }
                }
            }
            // value_counts of the columns holding a 1
            std::vector<int> col_index;
            std::vector<int> col_length;
            for (int c = 0; c < S; ++c) {
                int n = 0;
                for (int r = 0; r < S; ++r) {
                    n += testmat[static_cast<std::size_t>(r) * S + c] == 1;
                }
                if (n >= 3) {
                    col_index.push_back(c);
                    col_length.push_back(n);
                }
            }
            const std::size_t nrow = col_index.size();
            std::vector<double> mean_x;
            std::vector<int> continuous;
            std::vector<double> len;
            bool is_continue = false;
            auto finish_run = [&]() {
                double total = 0.0;
                for (const double a : len) {
                    total += a;
                }
                // Len = [a / sum(Len) for a in Len]; the Python sums numpy ints
                // here, so the sum is exact.
                std::vector<double> shares;
                for (const double a : len) {
                    shares.push_back(a / total);
                }
                double temp = 0.0;
                for (std::size_t i = 0; i < continuous.size(); ++i) {
                    temp = temp + static_cast<double>(continuous[i]) * shares[i];
                }
                mean_x.push_back(np::round_half_even(temp));
            };
            for (std::size_t c = 0; c + 1 < nrow; ++c) {
                const int current = col_index[c];
                const int next = col_index[c + 1];
                if (next - current == 1 && is_continue) {
                    continuous.push_back(next);
                    len.push_back(col_length[c + 1]);
                } else if (next - current == 1 && !is_continue) {
                    continuous = {current, next};
                    len = {static_cast<double>(col_length[c]), static_cast<double>(col_length[c + 1])};
                    is_continue = true;
                } else if (next - current != 1 && !is_continue) {
                    continuous = {current};
                    len = {static_cast<double>(col_length[c])};
                    is_continue = false;
                    finish_run();
                } else {
                    finish_run();
                    continuous = {current};
                    len = {static_cast<double>(col_length[c])};
                    is_continue = false;
                }
            }
            finish_run();

            std::vector<double> X = mean_x;
            std::sort(X.begin(), X.end());
            X.erase(std::unique(X.begin(), X.end()), X.end());
            const std::size_t x_size = X.size();
            for (std::size_t c = 0; c + 1 < x_size; ++c) {
                const int n = static_cast<int>(X[c]);
                const int m = static_cast<int>(X[c + 1]);
                const int span = std::abs(m - n);
                if (span > 1 && span <= params.max_width) {
                    if (span > 4) {
                        pairs.push_back({n, m - 2});
                    } else {
                        pairs.push_back({n, m});
                    }
                    int min1 = S;
                    int max1 = -1;
                    for (int cc = std::max(0, n - 1); cc < std::min(n + 2, S); ++cc) {
                        for (int r = 0; r < S; ++r) {
                            if (testmat[static_cast<std::size_t>(r) * S + cc] == 1) {
                                min1 = std::min(min1, r);
                                max1 = std::max(max1, r);
                            }
                        }
                    }
                    int min2 = S;
                    int max2 = -1;
                    for (int cc = std::max(0, m - 1); cc < std::min(m + 2, S); ++cc) {
                        for (int r = 0; r < S; ++r) {
                            if (testmat[static_cast<std::size_t>(r) * S + cc] == 1) {
                                min2 = std::min(min2, r);
                                max2 = std::max(max2, r);
                            }
                        }
                    }
                    if (max1 < 0 || max2 < 0) {
                        throw std::runtime_error("hicDetectStripes: empty column in the stripe pairing (Stripenn would fail here)");
                    }
                    int mn = std::min(min1, min2);
                    int mx = std::max(max1, max2);
                    if (ud == 1) {
                        mx = span > 4 ? m - 2 : m;
                    } else {
                        mn = static_cast<int>(X[c]);
                    }
                    min_vec.push_back(mn);
                    max_vec.push_back(mx);
                }
            }
        }
        for (std::size_t c = 0; c < pairs.size(); ++c) {
            RawStripe s;
            s.x = pairs[c].first;
            s.y = min_vec[c];
            s.w = pairs[c].second - pairs[c].first + 1;
            s.h = max_vec[c] - min_vec[c] + 1;
            if (s.h > 0 && s.w > 0) {
                s.total = np::sum_block(submat.data() + static_cast<std::size_t>(s.y) * S + s.x,
                                        static_cast<std::size_t>(s.h), static_cast<std::size_t>(s.w),
                                        static_cast<std::size_t>(S));
            } else {
                s.total = 0.0;
            }
            s.mean = s.total / static_cast<double>(s.h) / static_cast<double>(s.w);
            found.push_back(s);
        }
    }
    return found;
}

// ---------------------------------------------------------------------------
// search_frame

FrameResult search_frame(const BandMatrix& band, const ChromInfo& chrom, std::int64_t idx, double M,
                         const Parameters& params, FrameTrace* trace) {
    FrameResult out;
    const std::int64_t res = band.resolution();
    const std::int64_t rowsize = chrom.nbins;
    std::int64_t start = idx * 200 - 100;
    std::int64_t end = (idx + 1) * 200 + 99;
    if (end >= rowsize) {
        end = rowsize - 1;
    }
    if (idx == 0) {
        start = 0;
    }
    const std::int64_t framesize = end - start + 1;
    if (framesize <= 0) {
        return out;
    }
    std::vector<double> D = band.dense(start, end + 1, start, end + 1);
    for (double& v : D) {
        if (is_nan(v)) {
            v = 0.0;
        }
    }
    std::vector<std::int64_t> nonzero;
    for (std::int64_t c = 0; c < framesize; ++c) {
        bool any = false;
        for (std::int64_t r = 0; r < framesize; ++r) {
            if (D[static_cast<std::size_t>(r * framesize + c)] != 0.0) {
                any = true;
                break;
            }
        }
        if (any) {
            nonzero.push_back(c);
        }
    }
    const std::int64_t S = static_cast<std::int64_t>(nonzero.size());
    if (!(S > 10)) {
        return out;
    }
    std::vector<std::int64_t> start_array(static_cast<std::size_t>(framesize));
    std::vector<std::int64_t> end_array(static_cast<std::size_t>(framesize));
    for (std::int64_t j = 0; j < framesize; ++j) {
        start_array[static_cast<std::size_t>(j)] = (start + j) * res + 1;
        end_array[static_cast<std::size_t>(j)] = start_array[static_cast<std::size_t>(j)] + res - 1;
    }
    if (end_array.back() >= chrom.size) {
        end_array.back() = chrom.size;
    }
    out.submat.resize(static_cast<std::size_t>(S * S));
    out.start_array.resize(static_cast<std::size_t>(S));
    out.end_array.resize(static_cast<std::size_t>(S));
    for (std::int64_t r = 0; r < S; ++r) {
        for (std::int64_t c = 0; c < S; ++c) {
            out.submat[static_cast<std::size_t>(r * S + c)] =
                D[static_cast<std::size_t>(nonzero[static_cast<std::size_t>(r)] * framesize + nonzero[static_cast<std::size_t>(c)])];
        }
        out.start_array[static_cast<std::size_t>(r)] = start_array[static_cast<std::size_t>(nonzero[static_cast<std::size_t>(r)])];
        out.end_array[static_cast<std::size_t>(r)] = end_array[static_cast<std::size_t>(nonzero[static_cast<std::size_t>(r)])];
    }
    out.searched = true;
    const std::vector<RawStripe> raw = stripe_search(out.submat, static_cast<int>(S), M, params, trace);
    for (const RawStripe& s : raw) {
        Stripe st;
        const std::int64_t x = s.x;
        const std::int64_t y = s.y;
        const std::int64_t w = s.w;
        const std::int64_t h = s.h;
        auto at = [&](const std::vector<std::int64_t>& arr, std::int64_t i) {
            const std::int64_t k = py_index(i, S);
            if (k < 0) {
                throw std::runtime_error("hicDetectStripes: stripe outside its frame (Stripenn would fail here)");
            }
            return arr[static_cast<std::size_t>(k)];
        };
        st.pos1 = at(out.start_array, x);
        st.pos2 = at(out.end_array, x + w - 1);
        st.pos3 = at(out.start_array, y);
        st.pos4 = at(out.end_array, y + h - 1);
        st.length = st.pos4 - st.pos3 + 1;
        st.width = st.pos2 - st.pos1 + 1;
        st.total = s.total;
        st.mean = s.mean;
        st.num = idx;
        st.start = start;
        st.end = end;
        st.x = s.x;
        st.y = s.y;
        st.h = s.h;
        st.w = s.w;
        out.stripes.push_back(st);
    }
    // StripeSearch ends with RemoveRedundant(result, 'size') on the frame's own boxes.
    out.stripes = remove_redundant(out.stripes, RedundantBy::Size);
    return out;
}

// ---------------------------------------------------------------------------
// RemoveRedundant

std::vector<Stripe> remove_redundant(const std::vector<Stripe>& stripes, RedundantBy by) {
    const std::size_t total = stripes.size();
    if (total == 0) {
        return stripes;
    }
    std::vector<char> keep(total, 1);
    std::map<std::int64_t, std::vector<std::size_t>> by_num;
    for (std::size_t i = 0; i < total; ++i) {
        by_num[stripes[i].num].push_back(i);
    }
    for (const auto& [n, own] : by_num) {
        std::vector<std::size_t> merged = own;
        const auto next = by_num.find(n + 1);
        if (next != by_num.end()) {
            merged.insert(merged.end(), next->second.begin(), next->second.end());
            std::sort(merged.begin(), merged.end());
        }
        const std::size_t length = merged.size();
        for (std::size_t i = 0; i + 1 < length; ++i) {
            for (std::size_t j = i + 1; j < length; ++j) {
                const Stripe& a = stripes[merged[i]];
                const Stripe& b = stripes[merged[j]];
                const std::int64_t int_x = std::max<std::int64_t>(0, std::min(a.pos2, b.pos2) + 1 - std::max(a.pos1, b.pos1));
                const std::int64_t int_y = std::max<std::int64_t>(0, std::min(a.pos4, b.pos4) + 1 - std::max(a.pos3, b.pos3));
                const double s_x = static_cast<double>(int_x) / static_cast<double>(std::min(a.pos2 - a.pos1, b.pos2 - b.pos1));
                const double s_y = static_cast<double>(int_y) / static_cast<double>(std::min(a.pos4 - a.pos3, b.pos4 - b.pos3));
                if (s_x > 0.2 && s_y > 0.2) {
                    bool delete_first;
                    if (by == RedundantBy::Size) {
                        delete_first = static_cast<double>(a.h) / static_cast<double>(a.w) <=
                                       static_cast<double>(b.h) / static_cast<double>(b.w);
                    } else {
                        delete_first = a.pvalue > b.pvalue;
                    }
                    keep[delete_first ? merged[i] : merged[j]] = 0;
                }
            }
        }
    }
    std::vector<Stripe> out;
    for (std::size_t i = 0; i < total; ++i) {
        if (keep[i]) {
            out.push_back(stripes[i]);
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// pvalue

void compute_pvalues(const BandMatrix& band, const ChromInfo& chrom, const NullTables& null,
                     std::vector<Stripe>& stripes, unsigned int threads) {
    const std::int64_t res = band.resolution();
    const std::int64_t bs = static_cast<std::int64_t>(50000.0 / static_cast<double>(res));
    const double all_columns = static_cast<double>(null.columns);
    hicx::parallel_for(stripes.size(), threads, [&](std::size_t index) {
        Stripe& s = stripes[index];
        std::int64_t leftmost = s.pos1 - bs * res;
        std::int64_t rightmost = s.pos2 + bs * res;
        if (leftmost < 1) {
            leftmost = 1;
        }
        if (rightmost > chrom.size) {
            rightmost = chrom.size;
        }
        const std::int64_t col0 = floor_div(leftmost, res);
        const std::int64_t col1 = ceil_div(rightmost, res);
        const std::int64_t row0 = floor_div(s.pos3, res);
        const std::int64_t row1 = ceil_div(s.pos4, res);
        const std::int64_t R = row1 - row0;
        const std::int64_t C = col1 - col0;
        std::vector<double> mat = band.dense(row0, row1, col0, col1);
        for (double& v : mat) {
            if (is_nan(v)) {
                v = 0.0;
            }
        }
        // mat[:, bs:-bs], mat[:, :bs], mat[:, -bs:]
        auto column_means = [&](std::int64_t c0, std::int64_t c1) {
            normalise_slice(c0, c1, C);
            std::vector<double> means(static_cast<std::size_t>(R));
            std::vector<double> row;
            for (std::int64_t r = 0; r < R; ++r) {
                row.assign(mat.begin() + static_cast<std::ptrdiff_t>(r * C + c0),
                           mat.begin() + static_cast<std::ptrdiff_t>(r * C + std::max(c0, c1)));
                means[static_cast<std::size_t>(r)] =
                    c1 > c0 ? np::mean_of_row(row.data(), static_cast<std::size_t>(c1 - c0)) : kNaN;
            }
            return means;
        };
        const std::vector<double> center = column_means(bs, -bs);
        const std::vector<double> left = column_means(0, bs);
        const std::vector<double> right = column_means(-bs, C);
        const std::int64_t x1 = (s.pos1 - 1) / res;
        const std::int64_t x2 = s.pos2 / res;
        const std::int64_t y1 = (s.pos3 - 1) / res;
        const std::int64_t y2 = s.pos4 / res;

        std::vector<double> pvalues;
        pvalues.reserve(static_cast<std::size_t>(R));
        for (std::int64_t j = 0; j < R; ++j) {
            const double left_diff = center[static_cast<std::size_t>(j)] - left[static_cast<std::size_t>(j)];
            const double right_diff = center[static_cast<std::size_t>(j)] - right[static_cast<std::size_t>(j)];
            const std::vector<std::vector<double>>* sorted_left = nullptr;
            const std::vector<std::vector<double>>* sorted_right = nullptr;
            std::int64_t difference;
            if (x1 == y1) {
                difference = j;
                if (difference >= 400) {
                    difference = 399;
                }
                sorted_left = &null.sorted_leftdown;
                sorted_right = &null.sorted_rightdown;
            } else if (x2 == y2) {
                difference = y2 - y1 - j - 1;
                if (difference >= 400) {
                    difference = 399;
                }
                if (difference < 0) {
                    difference += 400;  // numpy's negative row index
                }
                sorted_left = &null.sorted_leftup;
                sorted_right = &null.sorted_rightup;
            } else {
                throw std::runtime_error("hicDetectStripes: stripe touches the diagonal at neither end (Stripenn would fail here)");
            }
            const std::vector<double>& bl = (*sorted_left)[static_cast<std::size_t>(difference)];
            const std::vector<double>& br = (*sorted_right)[static_cast<std::size_t>(difference)];
            auto count_ge = [](const std::vector<double>& sorted, double v) {
                if (is_nan(v)) {
                    return static_cast<std::size_t>(0);
                }
                return static_cast<std::size_t>(sorted.end() - std::lower_bound(sorted.begin(), sorted.end(), v));
            };
            const double p1 = static_cast<double>(count_ge(bl, left_diff)) / static_cast<double>(bl.size());
            const double p2 = static_cast<double>(count_ge(br, right_diff)) / static_cast<double>(br.size());
            double pval = std::max(p1, p2);
            if (pval == 0) {
                pval = 1 / all_columns;
            }
            pvalues.push_back(pval);
        }
        s.pvalue = np::median(pvalues);
    });
}

// ---------------------------------------------------------------------------
// scoringstripes

namespace {

struct Matrix {
    std::int64_t rows = 0;
    std::int64_t cols = 0;
    std::vector<double> v;
    Matrix() = default;
    Matrix(std::int64_t r, std::int64_t c, double fill = 0.0)
        : rows(r), cols(c), v(static_cast<std::size_t>(std::max<std::int64_t>(r, 0) * std::max<std::int64_t>(c, 0)), fill) {}
    double& at(std::int64_t r, std::int64_t c) { return v[static_cast<std::size_t>(r * cols + c)]; }
    double at(std::int64_t r, std::int64_t c) const { return v[static_cast<std::size_t>(r * cols + c)]; }
};

// expecMatrix: (len(y_seq) x len(x_seq)), exval[min(|x - y|, 399)]
Matrix expec_matrix(const std::vector<double>& exval, std::int64_t xs, std::int64_t xe, std::int64_t ys,
                    std::int64_t ye) {
    Matrix m(std::max<std::int64_t>(ye - ys, 0), std::max<std::int64_t>(xe - xs, 0));
    for (std::int64_t y = 0; y < m.rows; ++y) {
        for (std::int64_t x = 0; x < m.cols; ++x) {
            std::int64_t idx = std::llabs((xs + x) - (ys + y));
            if (idx >= 400) {
                idx = 399;
            }
            m.at(y, x) = exval[static_cast<std::size_t>(idx)];
        }
    }
    return m;
}

// np.divide(observed, expected + 1e-8), observed cut to the expected shape
// (Stripenn raises where the two differ, at an unaligned chromosome end).
Matrix observed_over_expected(const BandMatrix& band, std::int64_t row0, std::int64_t col0, Matrix expected) {
    for (double& e : expected.v) {
        e += .00000001;
    }
    Matrix out(expected.rows, expected.cols);
    for (std::int64_t r = 0; r < out.rows; ++r) {
        for (std::int64_t c = 0; c < out.cols; ++c) {
            const std::int64_t br = row0 + r;
            const std::int64_t bc = col0 + c;
            const double obs = (br < band.bins() && bc < band.bins()) ? band.get(br, bc) : kNaN;
            out.at(r, c) = obs / expected.at(r, c);
        }
    }
    return out;
}

// masking(matrix, mask_index_start, mask_index_end, start_index, end_index, direc)
void apply_mask(Matrix& m, std::int64_t mask_start, std::int64_t mask_end, std::int64_t start_index,
                std::int64_t end_index, int direc) {
    const std::int64_t rel_start = mask_start - start_index;
    const std::int64_t mask_size = mask_end - mask_start;
    const std::int64_t rel_end = rel_start + mask_size;
    const std::int64_t L = end_index - start_index + 1;
    std::vector<std::int64_t> idx;
    for (std::int64_t x = rel_start; x <= rel_end; ++x) {
        if (x >= 0 && x < L) {
            idx.push_back(x);
        }
    }
    if (idx.empty()) {
        return;
    }
    if (direc == 1) {
        for (const std::int64_t y : idx) {
            if (y < m.cols) {
                for (std::int64_t x = 0; x < m.rows; ++x) {
                    m.at(x, y) = kNaN;
                }
            }
        }
    } else {
        for (const std::int64_t x : idx) {
            if (x < m.rows) {
                for (std::int64_t y = 0; y < m.cols; ++y) {
                    m.at(x, y) = kNaN;
                }
            }
        }
    }
}

// stats.elementwise_product_sum(K, *COL)
std::vector<double> elementwise_product_sum(const std::vector<std::array<int, 3>>& K_columns,
                                            const std::vector<const std::vector<double>*>& cols) {
    // K_columns[j] holds column j of K: {K[0][j], K[1][j], K[2][j]}
    const std::size_t ncol_k = K_columns.size();
    std::vector<double> result;
    const std::size_t col_size = cols[0]->size();
    for (std::size_t i = 1; i + 1 < col_size; ++i) {
        double temp = 0.0;
        for (std::size_t j = 0; j < ncol_k; ++j) {
            const std::vector<double>& col = *cols[j];
            const double a = static_cast<double>(K_columns[j][0]) * col[i - 1];
            const double b = static_cast<double>(K_columns[j][1]) * col[i];
            const double c = static_cast<double>(K_columns[j][2]) * col[i + 1];
            temp += (a + b) + c;
        }
        if (ncol_k == 3 && temp < 0) {
            temp *= -1;
        }
        result.push_back(temp);
    }
    return result;
}

}  // namespace

void compute_scores(const BandMatrix& band, const ChromInfo& chrom, const std::vector<double>& exval,
                    const Parameters& params, std::vector<Stripe>& stripes, unsigned int threads) {
    const std::int64_t res = band.resolution();
    const std::int64_t bs = static_cast<std::int64_t>(50000.0 / static_cast<double>(res));
    const std::int64_t chrom_bin_size = chrom.nbins;

    bool use_mask = false;
    std::int64_t mask_start = 0;
    std::int64_t mask_end = 0;
    std::string mask_chr;
    std::int64_t mask_x_start = 0;
    std::int64_t mask_x_end = 0;
    if (params.mask != "0") {
        const std::size_t colon = params.mask.find(':');
        const std::size_t dash = params.mask.find('-', colon == std::string::npos ? 0 : colon);
        if (colon == std::string::npos || dash == std::string::npos) {
            throw std::runtime_error("hicDetectStripes: --mask must look like chr9:12345678-12345789");
        }
        mask_chr = params.mask.substr(0, colon);
        mask_start = std::stoll(params.mask.substr(colon + 1, dash - colon - 1));
        mask_end = std::stoll(params.mask.substr(dash + 1));
        mask_x_start = static_cast<std::int64_t>(static_cast<double>(mask_start) / static_cast<double>(res));
        mask_x_end = static_cast<std::int64_t>(static_cast<double>(mask_end) / static_cast<double>(res));
        use_mask = mask_chr == chrom.name;
    }

    static const std::vector<std::array<int, 3>> K_xl{{-1, -2, -1}, {1, 2, 1}};
    static const std::vector<std::array<int, 3>> K_xr{{1, 2, 1}, {-1, -2, -1}};
    static const std::vector<std::array<int, 3>> K_y{{1, 0, -1}, {2, 0, -2}, {1, 0, -1}};

    hicx::parallel_for(stripes.size(), threads, [&](std::size_t index) {
        Stripe& s = stripes[index];
        const double xs = static_cast<double>(s.pos1);
        const double xe = static_cast<double>(s.pos2);
        const double ys = static_cast<double>(s.pos3);
        const double ye = static_cast<double>(s.pos4);
        const std::int64_t x_start_index = static_cast<std::int64_t>(xs / static_cast<double>(res));
        const std::int64_t x_end_index = static_cast<std::int64_t>(xe / static_cast<double>(res));
        const std::int64_t y_start_index = static_cast<std::int64_t>(ys / static_cast<double>(res));
        const std::int64_t y_end_index = static_cast<std::int64_t>(ye / static_cast<double>(res));
        std::int64_t leftmost = x_start_index - bs;
        std::int64_t rightmost = x_end_index + bs;
        if (leftmost < 1) {
            leftmost = 1;
        }
        if (rightmost >= chrom_bin_size) {
            rightmost = chrom_bin_size - 1;
        }
        const std::int64_t row0 = floor_div(s.pos3, res);
        Matrix center = observed_over_expected(band, row0, floor_div(s.pos1, res),
                                               expec_matrix(exval, x_start_index, x_end_index, y_start_index, y_end_index));
        Matrix left = observed_over_expected(band, row0, leftmost,
                                             expec_matrix(exval, leftmost, x_start_index, y_start_index, y_end_index));
        Matrix right = observed_over_expected(band, row0, x_end_index,
                                              expec_matrix(exval, x_end_index + 1, rightmost + 1, y_start_index, y_end_index));

        if (use_mask && static_cast<double>(mask_start) > std::min(xs - 50000, ys) &&
            static_cast<double>(mask_start) < std::max(xe + 50000, ye)) {
            apply_mask(center, mask_x_start, mask_x_end, x_start_index, x_end_index, 1);
            apply_mask(center, mask_x_start, mask_x_end, y_start_index, y_end_index, 2);
            apply_mask(left, mask_x_start, mask_x_end, leftmost, x_start_index, 1);
            apply_mask(left, mask_x_start, mask_x_end, y_start_index, y_end_index, 2);
            apply_mask(right, mask_x_start, mask_x_end, x_end_index + 1, rightmost, 1);
            apply_mask(right, mask_x_start, mask_x_end, y_start_index, y_end_index, 2);
        }

        auto all_nan_columns = [](const Matrix& m) {
            std::vector<std::int64_t> cols;
            for (std::int64_t c = 0; c < m.cols; ++c) {
                bool all = true;
                for (std::int64_t r = 0; r < m.rows; ++r) {
                    if (!is_nan(m.at(r, c))) {
                        all = false;
                        break;
                    }
                }
                if (all) {
                    cols.push_back(c);
                }
            }
            return cols;
        };
        auto delete_columns = [](Matrix& m, const std::vector<std::int64_t>& cols) {
            Matrix out(m.rows, m.cols - static_cast<std::int64_t>(cols.size()));
            std::int64_t oc = 0;
            for (std::int64_t c = 0; c < m.cols; ++c) {
                if (std::find(cols.begin(), cols.end(), c) != cols.end()) {
                    continue;
                }
                for (std::int64_t r = 0; r < m.rows; ++r) {
                    out.at(r, oc) = m.at(r, c);
                }
                ++oc;
            }
            m = std::move(out);
        };
        auto delete_rows = [](Matrix& m, const std::vector<std::int64_t>& rows) {
            Matrix out(m.rows - static_cast<std::int64_t>(rows.size()), m.cols);
            std::int64_t orow = 0;
            for (std::int64_t r = 0; r < m.rows; ++r) {
                if (std::find(rows.begin(), rows.end(), r) != rows.end()) {
                    continue;
                }
                for (std::int64_t c = 0; c < m.cols; ++c) {
                    out.at(orow, c) = m.at(r, c);
                }
                ++orow;
            }
            m = std::move(out);
        };

        const std::vector<std::int64_t> center_del = all_nan_columns(center);
        const std::vector<std::int64_t> left_del = all_nan_columns(left);
        const std::vector<std::int64_t> right_del = all_nan_columns(right);
        std::vector<std::int64_t> rowdel;
        const bool same_start = s.pos1 == s.pos3;
        if (!center_del.empty()) {
            delete_columns(center, center_del);
            for (const std::int64_t x : center_del) {
                rowdel.push_back(same_start ? x : center.rows - 1 - x);
            }
        }
        if (!left_del.empty()) {
            delete_columns(left, left_del);
            for (const std::int64_t x : left_del) {
                rowdel.push_back(same_start ? x : left.rows - 1 - x);
            }
        }
        if (!right_del.empty()) {
            delete_columns(right, right_del);
            for (const std::int64_t x : right_del) {
                rowdel.push_back(same_start ? x : right.rows - 1 - x);
            }
        }
        if (!rowdel.empty()) {
            const std::int64_t n_rows = center.rows;
            std::vector<std::int64_t> unique_rows;
            for (std::int64_t r : rowdel) {
                const std::int64_t k = py_index(r, n_rows);
                if (k >= 0) {
                    unique_rows.push_back(k);
                }
            }
            std::sort(unique_rows.begin(), unique_rows.end());
            unique_rows.erase(std::unique(unique_rows.begin(), unique_rows.end()), unique_rows.end());
            if (!unique_rows.empty()) {
                delete_rows(center, unique_rows);
                delete_rows(left, unique_rows);
                delete_rows(right, unique_rows);
            }
        }
        for (double& v : center.v) {
            if (is_nan(v)) {
                v = 0.0;
            }
        }
        for (double& v : left.v) {
            if (is_nan(v)) {
                v = 0.0;
            }
        }
        for (double& v : right.v) {
            if (is_nan(v)) {
                v = 0.0;
            }
        }
        auto row_means = [](const Matrix& m) {
            std::vector<double> means(static_cast<std::size_t>(std::max<std::int64_t>(m.rows, 0)));
            for (std::int64_t r = 0; r < m.rows; ++r) {
                means[static_cast<std::size_t>(r)] =
                    np::mean_of_row(m.v.data() + static_cast<std::size_t>(r * m.cols), static_cast<std::size_t>(m.cols));
            }
            return means;
        };
        const std::vector<double> centerm = row_means(center);
        const std::vector<double> leftm = row_means(left);
        const std::vector<double> rightm = row_means(right);
        if (!center.v.empty()) {
            s.oe_mean = np::mean_of(center.v);
        }

        const std::vector<double> g_xl = elementwise_product_sum(K_xl, {&leftm, &centerm});
        const std::vector<double> g_xr = elementwise_product_sum(K_xr, {&centerm, &rightm});
        const std::vector<double> g_y = elementwise_product_sum(K_y, {&leftm, &centerm, &rightm});
        std::vector<double> diff;
        const std::size_t n = std::min({g_xl.size(), g_xr.size(), g_y.size()});
        for (std::size_t i = 0; i < n; ++i) {
            double g_x;
            if (is_nan(g_xl[i]) || is_nan(g_xr[i])) {
                g_x = kNaN;
            } else {
                g_x = std::min(g_xl[i], g_xr[i]);
            }
            const double d = g_x - g_y[i];
            if (d >= 0 || d < 0) {
                diff.push_back(d);
            }
        }
        const double avgdiff = np::mean_of(diff);
        s.stripiness = np::nanmedian(centerm) * avgdiff;
    });
}

// ---------------------------------------------------------------------------
// Python repr(float)

std::string py_repr(double value) {
    if (is_nan(value)) {
        return "nan";
    }
    if (std::isinf(value)) {
        return value < 0 ? "-inf" : "inf";
    }
    if (value == 0.0) {
        return std::signbit(value) ? "-0.0" : "0.0";
    }
    char buf[64];
    const auto result = std::to_chars(buf, buf + sizeof buf, value, std::chars_format::scientific);
    std::string text(buf, result.ptr);
    bool negative = false;
    if (text[0] == '-') {
        negative = true;
        text.erase(0, 1);
    }
    const std::size_t e = text.find('e');
    std::string mantissa = text.substr(0, e);
    const int exponent = std::stoi(text.substr(e + 1));
    std::string digits;
    for (const char c : mantissa) {
        if (c != '.') {
            digits += c;
        }
    }
    std::string out;
    if (exponent >= -4 && exponent < 16) {
        const int decpt = exponent + 1;
        if (decpt <= 0) {
            out = "0." + std::string(static_cast<std::size_t>(-decpt), '0') + digits;
        } else if (static_cast<std::size_t>(decpt) >= digits.size()) {
            out = digits + std::string(static_cast<std::size_t>(decpt) - digits.size(), '0') + ".0";
        } else {
            out = digits.substr(0, static_cast<std::size_t>(decpt)) + "." + digits.substr(static_cast<std::size_t>(decpt));
        }
    } else {
        out = digits.substr(0, 1);
        if (digits.size() > 1) {
            out += "." + digits.substr(1);
        }
        char ebuf[16];
        std::snprintf(ebuf, sizeof ebuf, "e%c%02d", exponent < 0 ? '-' : '+', std::abs(exponent));
        out += ebuf;
    }
    return negative ? "-" + out : out;
}

}  // namespace hicx::stripes
