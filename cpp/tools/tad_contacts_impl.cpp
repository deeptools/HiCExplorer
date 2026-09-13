// See tad_contacts_impl.hpp for what is reproduced and why.

#include "tad_contacts_impl.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "hicx/cool_adapter.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/json_lite.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/stats_ops.hpp"
#include "hicx/text_formats.hpp"
#include "hicx/text_table.hpp"

namespace hicx::tads {

namespace {

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string region_text(const std::string& chrom, std::int64_t start, std::int64_t end) {
    return chrom + ":" + std::to_string(start) + "-" + std::to_string(end);
}

// A hiCMatrix loaded for one region of a cool file: the bins it holds and
// whether hicmatrix applied the weights to it.
struct Region {
    std::int64_t origin = 0;
    std::int64_t length = 0;
    bool apply_correction = false;
};

// Stored entries of cooler's square region fetch, both triangles, after
// eliminate_zeros. Only "more than one" matters, so it stops at two.
bool more_than_one_stored(const ContactMatrix& matrix, std::int64_t lo, std::int64_t hi) {
    const CsrMatrix& m = matrix.matrix;
    const std::vector<std::int64_t>& indptr = m.indptr();
    const std::vector<std::int32_t>& indices = m.indices();
    const std::vector<double>& data = m.data();
    const bool upper = m.symmetry() == Symmetry::UpperTriangle;
    std::int64_t count = 0;
    for (std::int64_t row = lo; row < hi; ++row) {
        const auto begin = indices.begin() + indptr[static_cast<std::size_t>(row)];
        const auto end = indices.begin() + indptr[static_cast<std::size_t>(row) + 1];
        for (auto it = std::lower_bound(begin, end, lo); it != end && *it < hi; ++it) {
            const std::size_t k = static_cast<std::size_t>(it - indices.begin());
            if (data[k] == 0.0) {
                continue;
            }
            count += (upper && *it != row) ? 2 : 1;
            if (count > 1) {
                return true;
            }
        }
    }
    return false;
}

// hiCMatrix(pMatrixFile=cool, pChrnameList=['chrom:start-end']): cooler's
// parse_region checks, then region_to_extent (cooler/core/_rangequery.py:11-30),
// then the correction test of hicmatrix/lib/cool.py:209-217.
Region cool_region(const ContactMatrix& matrix, const std::string& chrom,
                   std::int64_t start, std::int64_t end) {
    const auto name = std::find(matrix.chrom_names.begin(), matrix.chrom_names.end(), chrom);
    if (name == matrix.chrom_names.end()) {
        throw std::runtime_error("Unknown sequence label: " + chrom);
    }
    const std::int64_t length =
        matrix.chrom_lengths[static_cast<std::size_t>(name - matrix.chrom_names.begin())];
    if (end < start) {
        throw std::runtime_error("End cannot be less than start: " +
                                 region_text(chrom, start, end));
    }
    if (start < 0 || end > length) {
        throw std::runtime_error("Genomic region out of bounds: " +
                                 region_text(chrom, start, end));
    }
    const std::optional<BinRange> range = matrix.bins.chrom_bin_range(chrom);
    if (!range.has_value()) {
        throw std::runtime_error("chromosome " + chrom + " has no bins in the matrix");
    }

    std::int64_t lo = 0;
    std::int64_t hi = 0;
    if (matrix.bin_size > 0) {
        // int(np.floor(start / binsize)) and int(np.ceil(end / binsize)); both
        // coordinates are non negative here.
        lo = range->first + start / matrix.bin_size;
        hi = range->first + (end + matrix.bin_size - 1) / matrix.bin_size;
    } else {
        const std::vector<CutInterval>& bins = matrix.bins.intervals();
        const auto first = bins.begin() + range->first;
        const auto last = bins.begin() + range->last;
        const auto by_start = [](const CutInterval& bin, std::int64_t value) {
            return bin.start < value;
        };
        const auto upper = std::upper_bound(
            first, last, start,
            [](std::int64_t value, const CutInterval& bin) { return value < bin.start; });
        lo = range->first + (upper - first) - 1;
        hi = range->first + (std::lower_bound(first, last, end, by_start) - first);
        if (lo < 0) {
            throw std::runtime_error("region starts before the first bin: " +
                                     region_text(chrom, start, end));
        }
    }

    Region region;
    region.origin = lo;
    region.length = std::max<std::int64_t>(0, hi - lo);
    if (matrix.weights.has_value() && more_than_one_stored(matrix, lo, lo + region.length)) {
        const std::vector<double>& weights = *matrix.weights;
        const bool all_nan = std::all_of(
            weights.begin() + lo, weights.begin() + lo + region.length,
            [](double weight) { return std::isnan(weight); });
        region.apply_correction = !all_nan;
    }
    return region;
}

// getRegionBinRange(chrom, pos, pos)[0] on a region loaded from a cool file:
// the region's own interval tree, so the index is local to the region.
std::int64_t local_bin(const ContactMatrix& matrix, const Region& region,
                       const std::string& chrom, std::int64_t position) {
    const std::vector<CutInterval>& bins = matrix.bins.intervals();
    if (region.length <= 0) {
        throw std::runtime_error("the region holds no bins, so no position can be looked up");
    }
    const auto first = bins.begin() + region.origin;
    const auto last = first + region.length;
    if (first->chrom != chrom) {
        throw std::runtime_error("chromosome: " + chrom + " name not found in matrix");
    }
    auto it = std::upper_bound(
        first, last, position,
        [](std::int64_t value, const CutInterval& bin) { return value < bin.start; });
    if (it != first) {
        --it;
        if (position < it->end) {
            return static_cast<std::int64_t>(it - first);
        }
    }
    throw std::runtime_error("no bin of " + chrom + " covers position " +
                             std::to_string(position) + " in the loaded region");
}

std::int64_t whole_bin(const ContactMatrix& matrix, const std::string& chrom,
                       std::int64_t position) {
    const std::optional<std::int64_t> bin = matrix.bins.bin_at(chrom, position);
    if (!bin.has_value()) {
        throw std::runtime_error("no bin of " + chrom + " covers position " +
                                 std::to_string(position));
    }
    return *bin;
}

std::pair<std::int64_t, std::int64_t> whole_range(const ContactMatrix& matrix,
                                                  const std::string& chrom,
                                                  std::int64_t start, std::int64_t end) {
    const auto range = matrix.bins.region_bin_range(chrom, start, end);
    if (!range.has_value()) {
        throw std::runtime_error("getRegionBinRange found no bins for " +
                                 region_text(chrom, start, end));
    }
    return *range;
}

// The inter-TAD hiCMatrix and the four boundary indices of one iteration of
// the loop, before any slicing.
struct Frame {
    std::int64_t origin = 0;
    std::int64_t length = 0;
    bool apply_correction = false;
    std::int64_t left = 0;
    std::int64_t outer_left = 0;
    std::int64_t outer_right = -1;
    std::optional<std::int64_t> right;
};

Frame inter_frame(const ContactMatrix& matrix, const std::vector<Domain>& chromosome,
                  std::size_t index) {
    const std::size_t n = chromosome.size();
    const Domain& row = chromosome[index];
    const Domain& previous = index >= 1 ? chromosome[index - 1] : row;
    const std::string& chrom = previous.chrom;
    const std::int64_t start = previous.start;
    const std::int64_t end = index + 1 < n ? chromosome[index + 1].end : row.end;
    const bool has_next = index + 1 < n;

    Frame frame;
    if (matrix.format == Format::Cool) {
        const Region region = cool_region(matrix, chrom, start, end);
        frame.origin = region.origin;
        frame.length = region.length;
        frame.apply_correction = region.apply_correction;
        frame.left = local_bin(matrix, region, chrom, row.start);
        frame.outer_left = 0;
        frame.outer_right = -1;
        if (has_next) {
            frame.right = local_bin(matrix, region, chrom, row.end);
        }
    } else {
        frame.origin = 0;
        frame.length = matrix.matrix.rows();
        frame.left = whole_bin(matrix, chrom, row.start);
        const auto outer = whole_range(matrix, chrom, start, end);
        frame.outer_left = outer.first;
        frame.outer_right = outer.second;
        if (has_next) {
            frame.right = whole_bin(matrix, chrom, row.end);
        }
    }
    return frame;
}

}  // namespace

// --------------------------------------------------------------------------

std::vector<Domain> read_domains(const std::string& path) {
    const TextTable table = TextTable::read_tsv(path);
    if (table.cols() < 6) {
        throw std::runtime_error("the TAD domains file '" + path +
                                 "' needs six columns, it has " +
                                 std::to_string(table.cols()));
    }
    std::vector<Domain> domains;
    domains.reserve(table.rows());
    for (std::size_t row = 0; row < table.rows(); ++row) {
        Domain domain;
        for (std::size_t column = 0; column < 6; ++column) {
            domain.text[column] = table.column(column).as_str(row);
        }
        domain.chrom = domain.text[0];
        domain.start = table.column(1).as_int(row);
        domain.end = table.column(2).as_int(row);
        domains.push_back(std::move(domain));
    }
    return domains;
}

std::vector<std::vector<Domain>> group_by_chromosome(const std::vector<Domain>& domains) {
    std::vector<std::vector<Domain>> groups;
    for (const Domain& domain : domains) {
        if (groups.empty() || groups.back().front().chrom != domain.chrom) {
            groups.emplace_back();
        }
        groups.back().push_back(domain);
    }
    return groups;
}

ContactMatrix load_contact_matrix(const std::string& path, bool is_cooler) {
    ContactMatrix result;
    if (is_cooler) {
        // The whole pixel table, raw. Every region the Python loads is a view of
        // it, with the weights applied per region at extraction time.
        CoolFile cool(path);
        result.format = Format::Cool;
        result.chrom_names = cool.chrom_names();
        result.chrom_lengths = cool.chrom_lengths();
        result.bins = BinTable(cool.read_bins());
        if (cool.has_column("weight")) {
            result.weights = cool.read_column("weight");
        }
        if (const json::Value* bin_size = cool.info_value("bin-size");
            bin_size != nullptr && bin_size->is_number()) {
            result.bin_size = bin_size->as_int();
        }
        result.matrix = cool.read_matrix();
        if (!result.matrix.has_lower_entries()) {
            result.matrix.set_symmetry(Symmetry::UpperTriangle);
        }
    } else {
        // check_cooler is false, so the tools take the whole-matrix path, and
        // hm.hiCMatrix picks its loader from the file name alone.
        result.format = Format::H5;
        MatrixData data;
        if (ends_with(path, ".h5")) {
            data = read_hicexplorer_h5(path);
        } else {
            data = read_cool(path).data;
        }
        data.matrix.symmetrize_in_place();
        result.bins = BinTable(std::move(data.cut_intervals));
        result.matrix = std::move(data.matrix);
    }
    sort_row_indices(result.matrix);
    return result;
}

void sort_row_indices(CsrMatrix& matrix) {
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    bool sorted = true;
    for (std::int64_t row = 0; row < matrix.rows() && sorted; ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin + 1; k < end; ++k) {
            if (indices[k] < indices[k - 1]) {
                sorted = false;
                break;
            }
        }
    }
    if (sorted) {
        return;
    }
    CsrMatrix::Arrays arrays = matrix.release();
    std::vector<std::pair<std::int32_t, double>> buffer;
    for (std::int64_t row = 0; row < arrays.rows; ++row) {
        const std::size_t begin =
            static_cast<std::size_t>(arrays.indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(arrays.indptr[static_cast<std::size_t>(row) + 1]);
        buffer.clear();
        for (std::size_t k = begin; k < end; ++k) {
            buffer.emplace_back(arrays.indices[k], arrays.data[k]);
        }
        std::stable_sort(buffer.begin(), buffer.end(),
                         [](const auto& a, const auto& b) { return a.first < b.first; });
        for (std::size_t k = begin; k < end; ++k) {
            arrays.indices[k] = buffer[k - begin].first;
            arrays.data[k] = buffer[k - begin].second;
        }
    }
    matrix = CsrMatrix::adopt(std::move(arrays));
}

std::pair<std::int64_t, std::int64_t> normalise_slice(std::int64_t start, std::int64_t stop,
                                                      std::int64_t length) {
    const auto clamp = [length](std::int64_t value) {
        if (value < 0) {
            value += length;
            return value < 0 ? std::int64_t{0} : value;
        }
        return value > length ? length : value;
    };
    const std::int64_t first = clamp(start);
    const std::int64_t last = std::max(first, clamp(stop));
    return {first, last};
}

TadGeometry tad_geometry(const ContactMatrix& matrix, const std::vector<Domain>& chromosome,
                         std::size_t index) {
    const std::size_t n = chromosome.size();
    const Domain& row = chromosome[index];
    TadGeometry geometry;

    if (matrix.format == Format::Cool) {
        const Region region = cool_region(matrix, row.chrom, row.start, row.end);
        const std::int64_t end = region.origin + region.length;
        geometry.intra = Block{region.origin, end, region.origin, end, region.apply_correction};
    } else {
        const auto range = whole_range(matrix, row.chrom, row.start, row.end);
        const auto [first, last] = normalise_slice(range.first, range.second, matrix.matrix.rows());
        geometry.intra = Block{first, last, first, last, false};
    }

    const Frame frame = inter_frame(matrix, chromosome, index);
    const auto block = [&frame](std::int64_t row_start, std::int64_t row_stop,
                                std::int64_t col_start, std::int64_t col_stop) {
        const auto rows = normalise_slice(row_start, row_stop, frame.length);
        const auto cols = normalise_slice(col_start, col_stop, frame.length);
        return Block{frame.origin + rows.first, frame.origin + rows.second,
                     frame.origin + cols.first, frame.origin + cols.second,
                     frame.apply_correction};
    };

    const bool has_previous = index >= 1;
    const bool has_next = index + 1 < n;
    if (has_next) {
        geometry.right = block(frame.left, *frame.right, *frame.right, frame.outer_right);
    }
    if (has_previous && has_next) {
        geometry.left = block(frame.outer_left, frame.left, frame.left, *frame.right);
    } else if (index >= 2 && !has_next) {
        // `i - 1 > 0 and i + 1 >= len(chromosome_list)`: right_boundary_index
        // is the value the previous iteration left behind, in that
        // iteration's coordinates.
        const Frame previous = inter_frame(matrix, chromosome, index - 1);
        geometry.left = block(frame.outer_left, frame.left, frame.left, *previous.right);
    }
    return geometry;
}

DType block_dtype(const ContactMatrix& matrix, const Block& block) {
    if (block.apply_correction) {
        return DType::Float64;  // matrix.data.astype(float)
    }
    return matrix.matrix.dtype_kind();
}

DenseBlock extract_block(const ContactMatrix& matrix, const Block& block) {
    DenseBlock out;
    out.rows = std::max<std::int64_t>(0, block.rows());
    out.cols = std::max<std::int64_t>(0, block.cols());
    out.values.assign(static_cast<std::size_t>(out.rows * out.cols), 0.0);
    if (out.rows == 0 || out.cols == 0) {
        return out;
    }

    const CsrMatrix& m = matrix.matrix;
    const std::vector<std::int64_t>& indptr = m.indptr();
    const std::vector<std::int32_t>& indices = m.indices();
    const std::vector<double>& data = m.data();
    const bool cool = matrix.format == Format::Cool;
    const bool upper = m.symmetry() == Symmetry::UpperTriangle;
    // A cool region went through eliminate_zeros, and fillLowerTriangle's CSR
    // addition drops zero results; only an h5 matrix that was not filled keeps
    // explicitly stored zeros, and scipy counts those in .nnz.
    const bool count_stored = !cool && !upper;
    const bool correct = cool && block.apply_correction && matrix.weights.has_value();
    const std::vector<double>* weights = correct ? &*matrix.weights : nullptr;

    const std::int64_t r0 = block.row_begin;
    const std::int64_t r1 = block.row_end;
    const std::int64_t c0 = block.col_begin;
    const std::int64_t c1 = block.col_end;
    const std::int64_t width = out.cols;

    const auto put = [&](std::int64_t i, std::int64_t j, double value) {
        if (cool) {
            if (weights != nullptr) {
                // instances_factors *= features_factors; data *= instances_factors
                value *= (*weights)[static_cast<std::size_t>(i)] *
                         (*weights)[static_cast<std::size_t>(j)];
            }
            if (std::isnan(value)) {
                value = 0.0;
            }
        }
        out.values[static_cast<std::size_t>((i - r0) * width + (j - c0))] = value;
        if (count_stored || value != 0.0) {
            ++out.nnz;
        }
    };
    const auto visit_row = [&](std::int64_t row, std::int64_t lo, std::int64_t hi,
                               const auto& emit) {
        const auto begin = indices.begin() + indptr[static_cast<std::size_t>(row)];
        const auto end = indices.begin() + indptr[static_cast<std::size_t>(row) + 1];
        for (auto it = std::lower_bound(begin, end, lo); it != end && *it < hi; ++it) {
            emit(static_cast<std::int64_t>(*it),
                 data[static_cast<std::size_t>(it - indices.begin())]);
        }
    };

    if (!upper) {
        for (std::int64_t i = r0; i < r1; ++i) {
            visit_row(i, c0, c1, [&](std::int64_t j, double value) { put(i, j, value); });
        }
        return out;
    }
    // The stored entries of rows r0..r1 that lie on or above the diagonal.
    for (std::int64_t i = r0; i < r1; ++i) {
        const std::int64_t lo = std::max(c0, i);
        if (lo < c1) {
            visit_row(i, lo, c1, [&](std::int64_t j, double value) { put(i, j, value); });
        }
    }
    // The mirrored entries below the diagonal: (k, j) with k > j is stored as
    // (j, k), so walk the rows j of the block's column range.
    for (std::int64_t j = c0; j < std::min(c1, r1); ++j) {
        const std::int64_t lo = std::max(r0, j + 1);
        if (lo < r1) {
            visit_row(j, lo, r1, [&](std::int64_t k, double value) { put(k, j, value); });
        }
    }
    return out;
}

std::string PyNumber::str() const {
    switch (kind) {
        case Kind::PyInt:
        case Kind::Int64:
            return std::to_string(integer);
        case Kind::Float32:
            return value_repr(real, "float32");
        case Kind::Float64:
        default:
            return npy::float_repr(real);
    }
}

namespace {

using Kind = PyNumber::Kind;

// numpy 1.26 value based casting: a Python int operand takes the type of the
// numpy operand.
Kind effective(Kind self, Kind other) { return self == Kind::PyInt ? other : self; }

bool both_float32(Kind a, Kind b) { return a == Kind::Float32 && b == Kind::Float32; }

}  // namespace

PyNumber py_add(const PyNumber& a, const PyNumber& b) {
    if (a.kind == Kind::PyInt && b.kind == Kind::PyInt) {
        return PyNumber::py_int(a.integer + b.integer);
    }
    const Kind ka = effective(a.kind, b.kind);
    const Kind kb = effective(b.kind, a.kind);
    if (ka == Kind::Int64 && kb == Kind::Int64) {
        return PyNumber::int64(a.integer + b.integer);
    }
    if (both_float32(ka, kb)) {
        return PyNumber::float32(static_cast<double>(static_cast<float>(a.as_double()) +
                                                     static_cast<float>(b.as_double())));
    }
    return PyNumber::float64(a.as_double() + b.as_double());
}

PyNumber py_divide(const PyNumber& a, const PyNumber& b) {
    if (a.kind == Kind::PyInt && b.kind == Kind::PyInt && b.integer == 0) {
        throw std::runtime_error("division by zero");
    }
    const Kind ka = effective(a.kind, b.kind);
    const Kind kb = effective(b.kind, a.kind);
    if (both_float32(ka, kb)) {
        return PyNumber::float32(static_cast<double>(static_cast<float>(a.as_double()) /
                                                     static_cast<float>(b.as_double())));
    }
    return PyNumber::float64(a.as_double() / b.as_double());
}

PyNumber block_sum(const DenseBlock& block, DType dtype) {
    const std::size_t rows = static_cast<std::size_t>(block.rows);
    const std::size_t cols = static_cast<std::size_t>(block.cols);
    if (dtype == DType::Integer) {
        std::int64_t total = 0;
        for (const double value : block.values) {
            total += static_cast<std::int64_t>(value);
        }
        return PyNumber::int64(total);
    }
    if (dtype == DType::Float32) {
        std::vector<float> row_sums(rows, 0.0F);
        for (std::size_t r = 0; r < rows; ++r) {
            float sum = 0.0F;
            for (std::size_t c = 0; c < cols; ++c) {
                sum += static_cast<float>(block.values[r * cols + c]);
            }
            row_sums[r] = sum;
        }
        return PyNumber::float32(
            static_cast<double>(npy::pairwise_sum(row_sums.data(), row_sums.size())));
    }
    std::vector<double> row_sums(rows, 0.0);
    for (std::size_t r = 0; r < rows; ++r) {
        double sum = 0.0;
        for (std::size_t c = 0; c < cols; ++c) {
            sum += block.values[r * cols + c];
        }
        row_sums[r] = sum;
    }
    return PyNumber::float64(npy::pairwise_sum(row_sums.data(), row_sums.size()));
}

RankTest rank_sum_test(std::span<const double> x, std::span<const double> y) {
    const auto has_nan = [](std::span<const double> values) {
        return std::any_of(values.begin(), values.end(),
                           [](double value) { return std::isnan(value); });
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (has_nan(x) || has_nan(y)) {
        return RankTest{nan, nan};
    }
    const stats::RanksumsResult result = stats::ranksums(x, y);
    return RankTest{result.statistic, result.pvalue};
}

}  // namespace hicx::tads
