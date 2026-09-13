#include "detect_loops_impl.hpp"

#include <stdexcept>

namespace hicx::loops {

namespace {

// The epsilon obs_exp_matrix substitutes for NaN and infinity
// (utilities.py:137,144 with pToEpsilon=True).
constexpr double kObsExpEpsilon = 0.000001;
// obs_exp_matrix_non_zero uses a different one (utilities.py:543).
constexpr double kObsExpNonZeroEpsilon = 0.000000001;

}  // namespace

DenseBlock dense_block(const CsrMatrix& matrix, std::int64_t row0, std::int64_t row1,
                       std::int64_t col0, std::int64_t col1) {
    DenseBlock block;
    block.rows = row1 > row0 ? row1 - row0 : 0;
    block.cols = col1 > col0 ? col1 - col0 : 0;
    block.values.assign(static_cast<std::size_t>(block.rows * block.cols), 0.0);
    if (block.rows == 0 || block.cols == 0) {
        return block;
    }
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& values = matrix.data();
    for (std::int64_t row = row0; row < row1; ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        // The column indices of a row are sorted, so the window is one range.
        const auto first = std::lower_bound(indices.begin() + static_cast<std::ptrdiff_t>(begin),
                                            indices.begin() + static_cast<std::ptrdiff_t>(end),
                                            static_cast<std::int32_t>(col0));
        for (auto it = first; it != indices.begin() + static_cast<std::ptrdiff_t>(end); ++it) {
            if (*it >= col1) {
                break;
            }
            const std::size_t k = static_cast<std::size_t>(std::distance(indices.begin(), it));
            const std::int64_t local_row = row - row0;
            const std::int64_t local_col = static_cast<std::int64_t>(*it) - col0;
            block.values[static_cast<std::size_t>(local_row * block.cols + local_col)] =
                values[k];
        }
    }
    return block;
}

std::vector<double> flatten_slice(const DenseBlock& block, std::int64_t row_start,
                                  std::int64_t row_stop, std::int64_t col_start,
                                  std::int64_t col_stop) {
    const AxisRange rows = numpy_slice(block.rows, row_start, row_stop);
    const AxisRange cols = numpy_slice(block.cols, col_start, col_stop);
    std::vector<double> flat;
    flat.reserve(static_cast<std::size_t>(rows.size() * cols.size()));
    for (std::int64_t r = rows.begin; r < rows.end; ++r) {
        for (std::int64_t c = cols.begin; c < cols.end; ++c) {
            flat.push_back(block.at(r, c));
        }
    }
    return flat;
}

double element(const CsrMatrix& matrix, std::int64_t row, std::int64_t col) {
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
    const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
    const auto first = std::lower_bound(indices.begin() + static_cast<std::ptrdiff_t>(begin),
                                        indices.begin() + static_cast<std::ptrdiff_t>(end),
                                        static_cast<std::int32_t>(col));
    if (first == indices.begin() + static_cast<std::ptrdiff_t>(end) || *first != col) {
        return 0.0;
    }
    return matrix.data()[static_cast<std::size_t>(std::distance(indices.begin(), first))];
}

ValueKind value_kind(const std::string& dtype) {
    if (dtype.rfind("int", 0) == 0 || dtype.rfind("uint", 0) == 0) {
        return ValueKind::Integer;
    }
    if (dtype == "float32") {
        return ValueKind::Float32;
    }
    return ValueKind::Float64;
}

DistanceGroups group_by_distance(const CsrMatrix& matrix) {
    DistanceGroups groups;
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::size_t nnz = matrix.stored_nnz();
    if (nnz == 0) {
        groups.offsets.assign(1, 0);
        return groups;
    }

    std::int64_t minimum = std::numeric_limits<std::int64_t>::max();
    std::int64_t maximum = std::numeric_limits<std::int64_t>::min();
    std::vector<std::int32_t> distances(nnz, 0);
    for (std::int64_t row = 0; row < matrix.rows(); ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t distance =
                std::abs(static_cast<std::int64_t>(indices[k]) - row);
            distances[k] = static_cast<std::int32_t>(distance);
            minimum = std::min(minimum, distance);
            maximum = std::max(maximum, distance);
        }
    }
    groups.min_distance = minimum;
    groups.max_distance = maximum;

    const std::size_t buckets = static_cast<std::size_t>(maximum - minimum + 1);
    groups.offsets.assign(buckets + 1, 0);
    for (const std::int32_t distance : distances) {
        ++groups.offsets[static_cast<std::size_t>(distance - minimum) + 1];
    }
    for (std::size_t i = 1; i < groups.offsets.size(); ++i) {
        groups.offsets[i] += groups.offsets[i - 1];
    }
    groups.positions.assign(nnz, 0);
    std::vector<std::int64_t> cursor(groups.offsets.begin(), groups.offsets.end() - 1);
    // A single forward pass over the CSR positions, so each bucket ends up in
    // CSR storage order, which is the order numpy's boolean mask preserves.
    for (std::size_t k = 0; k < nnz; ++k) {
        const std::size_t bucket = static_cast<std::size_t>(distances[k] - minimum);
        groups.positions[static_cast<std::size_t>(cursor[bucket]++)] =
            static_cast<std::int64_t>(k);
    }
    return groups;
}

std::vector<double> expected_interactions(const CsrMatrix& matrix,
                                          const DistanceGroups& groups,
                                          unsigned int threads) {
    const std::int64_t shape = matrix.rows();
    std::vector<double> expected(static_cast<std::size_t>(shape), 0.0);
    if (groups.empty() || shape == 0) {
        return expected;
    }
    const std::vector<double>& values = matrix.data();
    const std::int64_t buckets = groups.bucket_count();

    parallel_for(static_cast<std::size_t>(buckets), threads, [&](std::size_t index) {
        const std::int64_t distance =
            groups.min_distance + static_cast<std::int64_t>(index);
        if (distance >= shape) {
            return;  // outside the expected array, as np.zeros(shape[0]) is
        }
        const std::span<const std::int64_t> bucket =
            groups.bucket(static_cast<std::int64_t>(index));
        if (bucket.empty()) {
            return;
        }
        std::vector<double> gathered(bucket.size(), 0.0);
        for (std::size_t i = 0; i < bucket.size(); ++i) {
            gathered[i] = values[static_cast<std::size_t>(bucket[i])];
        }
        expected[static_cast<std::size_t>(distance)] =
            simd::pairwise_sum(gathered.data(), gathered.size());
    });

    // occurrences = np.arange(shape + 1, 1, -1), so occurrences[d] = shape + 1 - d.
    for (std::int64_t distance = 0; distance < shape; ++distance) {
        const double occurrences = static_cast<double>(shape + 1 - distance);
        double value = expected[static_cast<std::size_t>(distance)] / occurrences;
        if (std::isnan(value) || std::isinf(value)) {
            value = 0.0;
        }
        expected[static_cast<std::size_t>(distance)] = value;
    }
    return expected;
}

std::vector<double> expected_interactions_non_zero(const CsrMatrix& matrix) {
    const std::int64_t shape = matrix.rows();
    std::vector<double> expected(static_cast<std::size_t>(shape), 0.0);
    std::vector<double> occurrences(static_cast<std::size_t>(shape), 0.0);
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& values = matrix.data();
    // Sequential accumulation in nonzero() order, which is what the Python
    // loop at utilities.py:327-329 does.
    for (std::int64_t row = 0; row < matrix.rows(); ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t distance =
                std::abs(static_cast<std::int64_t>(indices[k]) - row);
            if (distance >= shape) {
                continue;
            }
            expected[static_cast<std::size_t>(distance)] += values[k];
            occurrences[static_cast<std::size_t>(distance)] += 1.0;
        }
    }
    for (std::int64_t distance = 0; distance < shape; ++distance) {
        double value = expected[static_cast<std::size_t>(distance)] /
                       occurrences[static_cast<std::size_t>(distance)];
        if (std::isnan(value) || std::isinf(value)) {
            value = 0.0;
        }
        expected[static_cast<std::size_t>(distance)] = value;
    }
    return expected;
}

CsrMatrix obs_exp_matrix(const CsrMatrix& matrix, const DistanceGroups& groups,
                         unsigned int threads) {
    const std::vector<double> expected =
        expected_interactions(matrix, groups, threads);
    const ValueKind kind = value_kind(matrix.dtype());
    const std::int64_t shape = matrix.rows();

    CsrMatrix::Arrays arrays;
    arrays.rows = matrix.rows();
    arrays.cols = matrix.cols();
    arrays.indptr = matrix.indptr();
    arrays.indices = matrix.indices();
    arrays.data = matrix.data();
    arrays.dtype = matrix.dtype();
    arrays.symmetry = Symmetry::Full;

    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    for (std::int64_t row = 0; row < matrix.rows(); ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t distance =
                std::abs(static_cast<std::int64_t>(indices[k]) - row);
            // np.ceil(distance / 2).astype(np.int32): the halved lookup.
            const std::int64_t lookup =
                static_cast<std::int64_t>(std::ceil(static_cast<double>(distance) / 2.0));
            const double divisor =
                lookup < shape ? expected[static_cast<std::size_t>(lookup)] : 0.0;
            // data.astype(np.float32) then np.divide against the float64
            // expected value, which promotes the quotient back to float64.
            double value =
                static_cast<double>(static_cast<float>(arrays.data[k])) / divisor;
            if (std::isnan(value)) {
                value = kObsExpEpsilon;
            } else if (std::isinf(value)) {
                value = kObsExpEpsilon;
            }
            arrays.data[k] = cast_to_kind(value, kind);
        }
    }
    return CsrMatrix::adopt(std::move(arrays));
}

CsrMatrix obs_exp_matrix_non_zero(const CsrMatrix& matrix, bool ligation_factor) {
    const std::vector<double> expected = expected_interactions_non_zero(matrix);
    const std::int64_t shape = matrix.rows();

    // row_sums = np.array(submatrix.sum(axis=1).T).flatten() and
    // total_interactions = submatrix.sum(), both taken before the float32
    // cast, so they are reduced in the dtype the matrix arrived in. scipy's
    // axis=1 sum is a matrix-vector product against a vector of ones, which
    // accumulates each row sequentially in that dtype; the scalar sum then
    // reduces the resulting dense vector with numpy's pairwise order.
    const bool single_precision = value_kind(matrix.dtype()) == ValueKind::Float32;
    std::vector<double> row_sums(static_cast<std::size_t>(shape), 0.0);
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& values = matrix.data();
    for (std::int64_t row = 0; row < matrix.rows(); ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        if (single_precision) {
            float sum = 0.0F;
            for (std::size_t k = begin; k < end; ++k) {
                sum += static_cast<float>(values[k]);
            }
            row_sums[static_cast<std::size_t>(row)] = static_cast<double>(sum);
        } else {
            double sum = 0.0;
            for (std::size_t k = begin; k < end; ++k) {
                sum += values[k];
            }
            row_sums[static_cast<std::size_t>(row)] = sum;
        }
    }
    double total = 0.0;
    if (single_precision) {
        std::vector<float> narrow(row_sums.size(), 0.0F);
        for (std::size_t i = 0; i < row_sums.size(); ++i) {
            narrow[i] = static_cast<float>(row_sums[i]);
        }
        total = static_cast<double>(npy::pairwise_sum(narrow.data(), narrow.size()));
    } else {
        total = npy::pairwise_sum(row_sums.data(), row_sums.size());
    }

    CsrMatrix::Arrays arrays;
    arrays.rows = matrix.rows();
    arrays.cols = matrix.cols();
    arrays.indptr = matrix.indptr();
    arrays.indices = matrix.indices();
    arrays.data = matrix.data();
    // The Python never casts back, so the result stays float32.
    arrays.dtype = "float32";
    arrays.symmetry = Symmetry::Full;

    for (std::int64_t row = 0; row < matrix.rows(); ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t column = static_cast<std::int64_t>(indices[k]);
            const std::int64_t distance = std::abs(column - row);
            double divisor =
                distance < shape ? expected[static_cast<std::size_t>(distance)] : 0.0;
            if (ligation_factor) {
                // The factor is formed in the dtype of the row sums, which is
                // the dtype the matrix arrived in, and only then multiplied
                // into the float64 expected value.
                const double factor =
                    single_precision
                        ? static_cast<double>(
                              static_cast<float>(row_sums[static_cast<std::size_t>(row)]) *
                              static_cast<float>(row_sums[static_cast<std::size_t>(column)]) /
                              static_cast<float>(total))
                        : row_sums[static_cast<std::size_t>(row)] *
                              row_sums[static_cast<std::size_t>(column)] / total;
                divisor *= factor;
            }
            double value =
                static_cast<double>(static_cast<float>(arrays.data[k])) / divisor;
            // The assignment target is a float32 array, so the quotient is
            // rounded to single precision on the way in.
            float stored = static_cast<float>(value);
            if (std::isnan(stored)) {
                stored = static_cast<float>(kObsExpNonZeroEpsilon);
            } else if (std::isinf(stored)) {
                stored = static_cast<float>(kObsExpNonZeroEpsilon);
            }
            arrays.data[k] = static_cast<double>(stored);
        }
    }
    CsrMatrix result = CsrMatrix::adopt(std::move(arrays));
    result.eliminate_zeros();
    return result;
}

PreselectionResult preselect_candidates(
    const CsrMatrix& obs_exp, const DistanceGroups& groups, double threshold,
    const std::map<std::int64_t, double>& threshold_by_genomic_distance,
    std::int64_t resolution, double obs_exp_threshold, unsigned int threads) {
    PreselectionResult result;
    result.mask.assign(obs_exp.stored_nnz(), 0);
    if (groups.empty()) {
        return result;
    }
    const std::int64_t buckets = groups.bucket_count();
    result.fits.assign(static_cast<std::size_t>(buckets), DistanceFit{});
    const std::vector<double>& values = obs_exp.data();
    const bool use_dictionary = !threshold_by_genomic_distance.empty();
    // fit_nbinom sees the .data array of this matrix, so its objective is
    // evaluated in the matrix's own dtype. See hicx::stats::NBinomPrecision.
    const stats::NBinomPrecision precision =
        value_kind(obs_exp.dtype()) == ValueKind::Float32
            ? stats::NBinomPrecision::Float32
            : stats::NBinomPrecision::Float64;

    std::vector<std::string> errors(static_cast<std::size_t>(buckets));

    parallel_for(static_cast<std::size_t>(buckets), threads, [&](std::size_t index) {
        const std::int64_t distance =
            groups.min_distance + static_cast<std::int64_t>(index);
        const std::span<const std::int64_t> bucket =
            groups.bucket(static_cast<std::int64_t>(index));
        DistanceFit& fit = result.fits[index];
        fit.distance = distance;
        fit.count = bucket.size();
        if (bucket.empty()) {
            // The Python still fits an empty array here, which yields NaN
            // parameters and then selects nothing because the masked value
            // list is empty as well. Skipping it produces the same mask.
            return;
        }

        double limit = threshold;
        if (use_dictionary) {
            const std::int64_t key = distance * resolution;
            const auto found = threshold_by_genomic_distance.find(key);
            if (found == threshold_by_genomic_distance.end()) {
                errors[index] = "no threshold for genomic distance " +
                                std::to_string(key) +
                                " in the --pValuePreselection file";
                return;
            }
            limit = found->second;
        }

        std::vector<double> data(bucket.size(), 0.0);
        for (std::size_t i = 0; i < bucket.size(); ++i) {
            data[i] = values[static_cast<std::size_t>(bucket[i])];
        }
        const stats::NBinomFit parameters = stats::fit_nbinom(data, precision);
        fit.size = parameters.size;
        fit.prob = parameters.prob;
        fit.degenerate = parameters.iterations == 0;

        // The tail probability and the two comparisons around it run in the
        // dtype of the obs/exp matrix, for the same reason the fit does.
        // scipy.special.betainc resolves to its 'fff->f' loop when its array
        // argument is float32, whatever the dtype of the two scalars, so
        // cnb.cdf and the `1 - cdf` that follows it are single precision, and
        // so is the comparison against the threshold, which rounds the
        // threshold to float32 as well.
        //
        // The one thing not reproduced exactly is scipy's single precision
        // cephes incbet: this evaluates the double precision one and rounds.
        // Measured against scipy over the corpus the two agree to within one
        // float32 ulp of the cdf, which the cancellation in `1 - cdf`
        // amplifies to about 5e-07 of a tail probability near 0.1. A pixel
        // whose tail sits closer than that to --pValuePreselection can be
        // selected differently, which is exactly the kind of disagreement
        // class E5 is defined to tolerate.
        if (precision == stats::NBinomPrecision::Float32) {
            const float threshold32 = static_cast<float>(obs_exp_threshold);
            const float limit32 = static_cast<float>(limit);
            const float size32 = static_cast<float>(parameters.size);
            const float prob32 = static_cast<float>(parameters.prob);
            for (std::size_t i = 0; i < bucket.size(); ++i) {
                const float value = static_cast<float>(data[i]);
                if (!(value >= threshold32)) {
                    continue;
                }
                const float cdf = static_cast<float>(
                    stats::betainc(static_cast<double>(size32),
                                   static_cast<double>(value + 1.0F),
                                   static_cast<double>(prob32)));
                const float p = 1.0F - cdf;
                if (p <= limit32) {
                    result.mask[static_cast<std::size_t>(bucket[i])] = 1;
                }
            }
            return;
        }
        for (std::size_t i = 0; i < bucket.size(); ++i) {
            if (!(data[i] >= obs_exp_threshold)) {
                continue;
            }
            const double p =
                stats::nbinom_sf(data[i], parameters.size, parameters.prob);
            if (p <= limit) {
                result.mask[static_cast<std::size_t>(bucket[i])] = 1;
            }
        }
    });

    for (const std::string& message : errors) {
        if (!message.empty()) {
            throw std::runtime_error(message);
        }
    }
    return result;
}

std::vector<Candidate> neighborhood_merge(const std::vector<Candidate>& candidates,
                                          std::int64_t window_size,
                                          const CsrMatrix& obs_exp,
                                          unsigned int threads) {
    std::vector<char> keep(candidates.size(), 0);
    const std::int64_t x_max = obs_exp.rows();
    const std::int64_t y_max = obs_exp.cols();

    parallel_for(candidates.size(), threads, [&](std::size_t index) {
        const Candidate& candidate = candidates[index];
        const std::int64_t start_x =
            (candidate.row - window_size) > 0 ? candidate.row - window_size : 0;
        const std::int64_t start_y =
            (candidate.col - window_size) > 0 ? candidate.col - window_size : 0;
        const std::int64_t end_x = candidate.row + window_size + 1 < x_max
                                       ? candidate.row + window_size + 1
                                       : x_max;
        const std::int64_t end_y = candidate.col + window_size + 1 < y_max
                                       ? candidate.col + window_size + 1
                                       : y_max;
        const DenseBlock block = dense_block(obs_exp, start_x, end_x, start_y, end_y);
        if (block.size() == 0) {
            return;
        }
        const double maximum = *std::max_element(block.values.begin(), block.values.end());
        if (maximum == element(obs_exp, candidate.row, candidate.col)) {
            keep[index] = 1;
        }
    });

    std::vector<Candidate> selected;
    selected.reserve(candidates.size());
    for (std::size_t index = 0; index < candidates.size(); ++index) {
        if (keep[index] != 0) {
            selected.push_back(candidates[index]);
        }
    }
    return selected;
}

RegionTestResult candidate_region_test(const CsrMatrix& obs_exp,
                                       const std::vector<Candidate>& candidates,
                                       std::int64_t window_size, double p_value,
                                       std::int64_t peak_window_size,
                                       unsigned int threads) {
    RegionTestResult result;
    std::vector<char> accepted(candidates.size(), 0);
    std::vector<double> pvalues(candidates.size(), 0.0);
    const std::int64_t x_max = obs_exp.rows();
    const std::int64_t y_max = obs_exp.cols();
    const std::int64_t pw = peak_window_size;

    parallel_for(candidates.size(), threads, [&](std::size_t index) {
        const Candidate& candidate = candidates[index];
        const std::int64_t start_x =
            (candidate.row - window_size) > 0 ? candidate.row - window_size : 0;
        const std::int64_t start_y =
            (candidate.col - window_size) > 0 ? candidate.col - window_size : 0;
        const std::int64_t end_x = candidate.row + window_size + 1 < x_max
                                       ? candidate.row + window_size + 1
                                       : x_max;
        const std::int64_t end_y = candidate.col + window_size + 1 < y_max
                                       ? candidate.col + window_size + 1
                                       : y_max;
        const DenseBlock block = dense_block(obs_exp, start_x, end_x, start_y, end_y);
        // `if len(neighborhood) == 0` tests the number of rows, not the number
        // of elements, because len() of a 2D array is its first dimension.
        if (block.rows == 0) {
            return;
        }
        const double value = element(obs_exp, candidate.row, candidate.col);

        // np.array(np.where(block == value)).flatten() is
        // [rows..., columns...]; the code then reads element 0 and element 1.
        std::int64_t first_row = -1;
        std::int64_t second_entry = -1;
        std::int64_t first_column = -1;
        std::int64_t matches = 0;
        for (std::int64_t r = 0; r < block.rows && matches < 2; ++r) {
            for (std::int64_t c = 0; c < block.cols; ++c) {
                if (block.at(r, c) != value) {
                    continue;
                }
                if (matches == 0) {
                    first_row = r;
                    first_column = c;
                } else {
                    second_entry = r;
                }
                ++matches;
                if (matches >= 2) {
                    break;
                }
            }
        }
        if (matches == 0) {
            return;  // the Python raises IndexError here; it cannot happen
        }
        const std::int64_t peak_row = first_row;
        const std::int64_t peak_col = matches > 1 ? second_entry : first_column;

        const std::vector<double> peak = flatten_slice(
            block, peak_row - pw, peak_row + pw + 1, peak_col - pw, peak_col + pw + 1);

        std::vector<double> background;
        const auto extend = [&background](const std::vector<double>& part) {
            background.insert(background.end(), part.begin(), part.end());
        };
        // top to peak, then peak to bottom, then right middle, then left
        // middle, in exactly that order (:546-557).
        extend(flatten_slice(block, 0, peak_row - pw, 0, block.cols));
        extend(flatten_slice(block, peak_row + pw + 1, block.rows, 0, block.cols));
        extend(flatten_slice(block, peak_row - pw, peak_row + pw + 1, peak_col + pw + 1,
                             block.cols));
        extend(flatten_slice(block, peak_row - pw, peak_row + pw + 1, 0, peak_col - pw));

        if (static_cast<std::int64_t>(background.size()) < window_size) {
            return;
        }
        if (static_cast<std::int64_t>(peak.size()) < window_size) {
            return;
        }
        if (numpy_mean(peak) < numpy_mean(background)) {
            return;
        }
        const double peak_max = *std::max_element(peak.begin(), peak.end());
        const double background_max = *std::max_element(background.begin(), background.end());
        if (peak_max < background_max) {
            return;
        }

        std::vector<double> horizontal;
        {
            const std::vector<double> top =
                flatten_slice(block, 0, peak_row - pw, peak_col - pw, peak_col + pw + 1);
            const std::vector<double> bottom = flatten_slice(
                block, peak_row + pw + 1, block.rows, peak_col - pw, peak_col + pw + 1);
            horizontal.insert(horizontal.end(), top.begin(), top.end());
            horizontal.insert(horizontal.end(), bottom.begin(), bottom.end());
        }
        std::vector<double> vertical;
        {
            const std::vector<double> left =
                flatten_slice(block, peak_row - pw, peak_row + pw + 1, 0, peak_col - pw);
            const std::vector<double> right = flatten_slice(
                block, peak_row - pw, peak_row + pw + 1, peak_col + pw + 1, block.cols);
            vertical.insert(vertical.end(), left.begin(), left.end());
            vertical.insert(vertical.end(), right.begin(), right.end());
        }
        std::vector<double> bottom_left;
        {
            const std::vector<double> a =
                flatten_slice(block, peak_row, block.rows, 0, peak_col - pw);
            const std::vector<double> b = flatten_slice(block, peak_row + pw + 1,
                                                        block.rows, peak_col - pw,
                                                        peak_col + 1);
            bottom_left.insert(bottom_left.end(), a.begin(), a.end());
            bottom_left.insert(bottom_left.end(), b.begin(), b.end());
        }

        std::vector<double> sorted_peak = peak;
        std::sort(sorted_peak.begin(), sorted_peak.end());

        int accept_count = 0;
        for (std::vector<double>* data : {&bottom_left, &horizontal, &vertical}) {
            std::vector<double> sorted_data = *data;
            std::sort(sorted_data.begin(), sorted_data.end());
            const stats::RanksumsResult test = stats::ranksums(sorted_peak, sorted_data);
            if (test.pvalue <= p_value) {
                ++accept_count;
            }
        }
        if (accept_count < 3) {
            return;
        }
        std::vector<double> sorted_background = background;
        std::sort(sorted_background.begin(), sorted_background.end());
        const stats::RanksumsResult test =
            stats::ranksums(sorted_peak, sorted_background);
        if (test.pvalue <= p_value) {
            accepted[index] = 1;
            pvalues[index] = test.pvalue;
        }
    });

    for (std::size_t index = 0; index < candidates.size(); ++index) {
        if (accepted[index] != 0) {
            result.candidates.push_back(candidates[index]);
            result.pvalues.push_back(pvalues[index]);
        }
    }
    return result;
}

}  // namespace hicx::loops
