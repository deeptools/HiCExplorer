// hicexplorer._hicx.hicAggregateContacts: the public functions of
// hicexplorer/hicAggregateContacts.py and the data of its two drawing
// functions.
//
// read_bed_per_chrom, aggregate_contacts, aggregate_contacts_per_row,
// count_contacts, get_outlier_indices, compute_clusters, cluster_matrices and
// compute_avg are answered here. They keep the 3.7.6 data structures: agg_info
// is the caller's dict of dicts and is filled in place (submatrices as
// numpy.matrix, diagonals as 1 x w numpy.matrix, centre values as numpy
// scalars of the matrix dtype, positions as given), seen_chrs and
// agg_info["seen"] are the caller's lists. The arithmetic is
// hicx::aggregate (cpp/core/include/hicx/aggregate_contacts.hpp) and the
// clustering hicx::cluster (KMeans and ward AgglomerativeClustering with
// scikit-learn 1.3.2's labels), the matrix queries hicmatrix_queries.hpp.
//
// plot_aggregated_contacts and plot_diagnostic_heatmaps stay matplotlib
// drawing code (cpp/scripts/public_api_allowlist.tsv); the helpers whose names
// start with an underscore compute what they draw and write: the cluster
// averages and their log records, the --outFilePrefixMatrix and
// --outFileContactPairs files, and the sorted diagonal rows with their means.
//
// Deviations from 3.7.6, documented in cases_hicAggregateContacts.py:
//  * method='spectral' raises NotImplementedError (main() never reaches
//    SpectralClustering, cpp/core/include/hicx/clustering.hpp);
//  * largeRegionsOperation other than 'first', 'last' and 'center' raises
//    TypeError naming HiCExplorer 4 (argparse allows only those);
//  * the matrix of ma must be csr, compute_avg and the clustering take float64
//    submatrices of one shape, and the BED lines must be str;
//  * a window that reaches past the matrix raises ValueError (it cannot for a
//    hiCMatrix whose matrix matches its bins).

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicmatrix_queries.hpp"
#include "hicx/aggregate_contacts.hpp"
#include "hicx/clustering.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/numpy_sort.hpp"
#include "hicx/py_text.hpp"
#include "hicx/utilities.hpp"
#include "utilities.hpp"

namespace hicx::python {

namespace {

namespace ag = hicx::aggregate;

constexpr const char* kLogger = "hicexplorer.hicAggregateContacts";

py::object numpy() {
    return py::module_::import("numpy");
}

[[noreturn]] void system_exit(const std::string& message) {
    raise_error("SystemExit", message);
}

[[noreturn]] void unbound(const char* variable) {
    raise_error("UnboundLocalError", std::string("cannot access local variable '") + variable +
                                         "' where it is not associated with a value");
}

bool assertions_enabled() {
    return py::module_::import("sys").attr("flags").attr("optimize").cast<int>() == 0;
}

bool equals(py::handle a, const char* text) {
    return hicobj::py_equal(a, py::str(text));
}

void add_one(py::handle dict, const char* key) {
    dict[py::str(key)] = py::int_(hicobj::int64_of(hicobj::item(dict, py::str(key))) + 1);
}

py::object utilities_change(const std::string& name) {
    return to_str(hicx::utilities::change_chrom_name(name));
}

std::vector<std::int64_t> index_vector_of(py::handle indices) {
    std::vector<std::int64_t> out;
    for (const py::handle index : indices) {
        out.push_back(to_int64(index));
    }
    return out;
}

// --------------------------------------------------------------------------
// the csr matrix of ma, read in place

class CsrView {
  public:
    explicit CsrView(py::handle matrix_object) : sparse_(sparse(matrix_object, "count_contacts", "r")) {
        indptr_ = py::array::ensure(sparse_.indptr, py::array::c_style);
        indices_ = py::array::ensure(sparse_.indices, py::array::c_style);
        data_ = py::array::ensure(sparse_.data, py::array::c_style);
        if (!indptr_ || !indices_ || !data_) {
            throw_python_error();
        }
        kind_ = dtype_of(data_);
        if (kind_ != DType::Float64 && kind_ != DType::Float32 && kind_ != DType::Int32 && kind_ != DType::Int64) {
            raise_unsupported_dtype(data_, "count_contacts");
        }
        dtype_ = py::reinterpret_borrow<py::object>(data_.attr("dtype"));
    }

    [[nodiscard]] std::int64_t rows() const noexcept { return sparse_.rows; }
    [[nodiscard]] std::int64_t cols() const noexcept { return sparse_.cols; }

    // matrix[row, col], summing duplicates as scipy's indexing does.
    [[nodiscard]] double at(std::int64_t row, std::int64_t col) const {
        double value = 0.0;
        for (std::int64_t p = pointer(row); p < pointer(row + 1); ++p) {
            if (index(p) == col) {
                value += datum(p);
            }
        }
        return value;
    }

    // matrix[r0:r0 + w, :][:, c0:c0 + w].todense().astype(float), row major
    [[nodiscard]] std::vector<double> window(std::int64_t r0, std::int64_t c0, std::int64_t w) const {
        if (r0 < 0 || c0 < 0 || r0 + w > rows() || c0 + w > cols()) {
            raise_error("ValueError", "count_contacts: HiCExplorer 4 reads windows inside the matrix only");
        }
        const auto width = static_cast<std::size_t>(w);
        std::vector<double> out(width * width, 0.0);
        for (std::int64_t i = 0; i < w; ++i) {
            for (std::int64_t p = pointer(r0 + i); p < pointer(r0 + i + 1); ++p) {
                const std::int64_t c = index(p);
                if (c >= c0 && c < c0 + w) {
                    out[static_cast<std::size_t>(i) * width + static_cast<std::size_t>(c - c0)] += datum(p);
                }
            }
        }
        return out;
    }

    // the numpy scalar matrix[row, col] has
    [[nodiscard]] py::object scalar(double value) const { return dtype_.attr("type")(value); }

  private:
    [[nodiscard]] std::int64_t pointer(std::int64_t row) const {
        if (dtype_of(indptr_) == DType::Int32) {
            return static_cast<const std::int32_t*>(indptr_.data())[row];
        }
        return static_cast<const std::int64_t*>(indptr_.data())[row];
    }
    [[nodiscard]] std::int64_t index(std::int64_t p) const {
        if (dtype_of(indices_) == DType::Int32) {
            return static_cast<const std::int32_t*>(indices_.data())[p];
        }
        return static_cast<const std::int64_t*>(indices_.data())[p];
    }
    [[nodiscard]] double datum(std::int64_t p) const {
        switch (kind_) {
            case DType::Float32: return static_cast<const float*>(data_.data())[p];
            case DType::Int32: return static_cast<const std::int32_t*>(data_.data())[p];
            case DType::Int64: return static_cast<double>(static_cast<const std::int64_t*>(data_.data())[p]);
            default: return static_cast<const double*>(data_.data())[p];
        }
    }

    Sparse sparse_;
    py::array indptr_;
    py::array indices_;
    py::array data_;
    DType kind_ = DType::Float64;
    py::object dtype_;
};

// numpy.matrix of a rows by cols float64 block
py::object np_matrix(const std::vector<double>& values, std::size_t rows, std::size_t cols) {
    py::array_t<double> array({static_cast<py::ssize_t>(rows), static_cast<py::ssize_t>(cols)});
    if (!values.empty()) {
        std::memcpy(array.mutable_data(), values.data(), values.size() * sizeof(double));
    }
    return numpy().attr("asmatrix")(array);
}

py::array float_array(const std::vector<double>& values, const std::vector<py::ssize_t>& shape) {
    py::array_t<double> out(shape);
    if (!values.empty()) {
        std::memcpy(out.mutable_data(), values.data(), values.size() * sizeof(double));
    }
    return out;
}

// --------------------------------------------------------------------------
// count_contacts

struct Context {
    py::handle ma;
    std::optional<CsrView> matrix;
    // the bins of the interval tree of each chromosome, read once
    std::map<std::string, std::vector<hicobj::TreeBin>> bins;

    const CsrView& view() {
        if (!matrix) {
            matrix.emplace(ma.attr("matrix"));
        }
        return *matrix;
    }

    // ma.getRegionBinRange(chrname, start, end)
    std::optional<std::pair<std::int64_t, std::int64_t>> region(py::handle chrname, py::handle start, py::handle end) {
        const py::object chrom = hicobj::region_chrom(ma, chrname);
        const std::int64_t s = hicobj::int64_of(start);
        const std::int64_t e = hicobj::int64_of(end);
        const std::string key = repr_of(chrom);
        auto it = bins.find(key);
        if (it == bins.end()) {
            it = bins.emplace(key, hicobj::tree_bins(ma, chrom)).first;
        }
        const auto first = hicobj::first_overlap(it->second, s, s + 1);
        const auto last = first ? hicobj::first_overlap(it->second, e, e + 1) : std::nullopt;
        if (!first || !last) {
            log(hicobj::kHicMatrixLogger, "error", "Index error");
            return std::nullopt;
        }
        return std::make_pair(*first, *last);
    }
};

ag::LargeRegions large_regions(py::handle operation) {
    if (equals(operation, "first")) {
        return ag::LargeRegions::First;
    }
    if (equals(operation, "last")) {
        return ag::LargeRegions::Last;
    }
    if (equals(operation, "center")) {
        return ag::LargeRegions::Center;
    }
    raise_error("TypeError", "count_contacts: HiCExplorer 4 takes largeRegionsOperation 'first', 'last' or 'center'");
}

std::optional<std::string> strand(py::handle value) {
    if (value.is_none()) {
        return std::nullopt;
    }
    if (is_str(value)) {
        return utf8(value);
    }
    return std::string("\x01not a strand");
}

std::int64_t floor_div(std::int64_t a, std::int64_t b) {
    if (b == 0) {
        raise_error("ZeroDivisionError", "integer division or modulo by zero");
    }
    std::int64_t q = a / b;
    if ((a % b != 0) && ((a < 0) != (b < 0))) {
        --q;
    }
    return q;
}

void count_contacts(Context& context, py::handle interval, py::handle M_half_object, py::handle mode,
                    py::handle agg_info, py::handle largeRegionsOperation, py::handle seen_chrs, py::handle range,
                    py::handle transform) {
    py::handle ma = context.ma;
    log(kLogger, "debug", "interval " + str_of(interval));
    const py::object first = hicobj::item(interval, 0);
    const py::object second = hicobj::item(interval, 1);
    const std::size_t n1 = hicobj::length(first);
    const std::size_t n2 = hicobj::length(second);
    py::object orientation1 = py::none();
    py::object orientation2 = py::none();
    if (!((n1 == 3 && n2 == 3) || (n1 == 4 && n2 == 4))) {
        unbound("chrom1");
    }
    py::object chrom1 = hicobj::item(first, 0);
    py::object start1 = hicobj::item(first, 1);
    py::object end1 = hicobj::item(first, 2);
    py::object chrom2 = hicobj::item(second, 0);
    py::object start2 = hicobj::item(second, 1);
    py::object end2 = hicobj::item(second, 2);
    if (n1 == 4) {
        orientation1 = hicobj::item(first, 3);
        orientation2 = hicobj::item(second, 3);
    }
    const py::object coord = hicobj::item(agg_info, py::str("chrom_coord"));
    if (!hicobj::py_contains(coord, chrom1) || !hicobj::py_contains(coord, chrom2)) {
        return;
    }
    const py::object coord1 = hicobj::item(coord, chrom1);
    const py::object coord2 = hicobj::item(coord, chrom2);
    if (hicobj::py_less(hicobj::item(coord1, 1), hicobj::py_int(end1)) ||
        hicobj::py_less(hicobj::item(coord2, 1), hicobj::py_int(end2))) {
        return;
    }
    if (hicobj::py_less(hicobj::py_int(start1), hicobj::item(coord1, 0)) ||
        hicobj::py_less(hicobj::py_int(start2), hicobj::item(coord2, 0))) {
        return;
    }
    const auto range1 = context.region(to_string(chrom1), start1, end1);
    const auto range2 = context.region(to_string(chrom2), start2, end2);
    if (range1 == range2) {
        return;  // includes None == None
    }
    if (!range1 || !range2) {
        return;
    }
    const ag::LargeRegions operation = large_regions(largeRegionsOperation);
    std::int64_t bin1 = ag::pick_bin(range1->first, range1->second, operation);
    std::int64_t bin2 = ag::pick_bin(range2->first, range2->second, operation);
    if (bin1 > bin2) {
        if (hicobj::py_equal(chrom1, chrom2)) {
            std::swap(bin1, bin2);
        } else {
            // chromosome, start, end and bin move; the orientations do not
            std::swap(chrom1, chrom2);
            std::swap(start1, start2);
            std::swap(end1, end2);
            std::swap(bin1, bin2);
        }
    }
    add_one(agg_info, "counter");
    const std::int64_t counter = hicobj::int64_of(hicobj::item(agg_info, py::str("counter")));
    if (counter % 50000 == 0) {
        log(kLogger, "info", "Number of contacts considered: " + hicx::npy::int_with_thousands_separator(counter));
    }
    const std::int64_t bin_size = hicobj::bin_size_value(ma);
    const auto chrom1_range = hicobj::range_values(hicobj::chr_bin_range(ma, to_string(chrom1)));
    const auto chrom2_range = hicobj::range_values(hicobj::chr_bin_range(ma, to_string(chrom2)));

    py::list pair;
    pair.append(chrom1);
    pair.append(chrom2);
    if (!hicobj::py_contains(seen_chrs, pair)) {
        hicobj::item(hicobj::item(agg_info, py::str("agg_total")), chrom1)[chrom2] = py::int_(0);
        for (const char* key : {"agg_matrix", "agg_diagonals", "agg_contact_position", "agg_center_values"}) {
            hicobj::item(hicobj::item(agg_info, py::str(key)), chrom1)[chrom2] = py::list();
        }
        seen_chrs.attr("append")(pair);
    }
    if (equals(mode, "intra-chr")) {
        if (!is_str(range)) {
            raise_error("AttributeError", "'" + type_name(range) + "' object has no attribute 'split'");
        }
        const std::vector<std::string_view> parts = hicx::pytext::split(utf8(range), ":");
        if (parts.size() != 2) {
            raise_error("ValueError", parts.size() < 2 ? "not enough values to unpack (expected 2, got " +
                                                             std::to_string(parts.size()) + ")"
                                                       : std::string("too many values to unpack (expected 2)"));
        }
        const std::int64_t min_bins = floor_div(to_int64(parse_int(to_str(parts[0]))), bin_size);
        const std::int64_t max_bins = floor_div(to_int64(parse_int(to_str(parts[1]))), bin_size);
        const std::int64_t distance = std::llabs(bin2 - bin1);
        if (min_bins > distance || distance > max_bins) {
            return;
        }
    }
    const py::object seen = hicobj::item(agg_info, py::str("seen"));
    const py::tuple key = py::make_tuple(bin1, bin2);
    if (hicobj::py_contains(seen, key)) {
        return;
    }
    seen.attr("append")(key);
    const std::int64_t M_half = to_int64(M_half_object);
    if (bin1 - M_half < chrom1_range.first || bin1 + M_half >= chrom1_range.second) {
        log(kLogger, "info", "The given interval exceeds the chromosome range on " + str_of(chrom1) + ". It is skipped.");
        return;
    }
    if (bin2 - M_half < chrom2_range.first || bin2 + M_half >= chrom2_range.second) {
        log(kLogger, "info", "The given interval exceeds the chromosome range on " + str_of(chrom2) + ". It is skipped.");
        return;
    }
    const std::int64_t width = 2 * M_half + 1;
    if (width <= 0) {
        add_one(agg_info, "empty_mat");  // an empty slice sums to zero
        return;
    }
    const CsrView& view = context.view();
    const std::vector<double> window = view.window(bin1 - M_half, bin2 - M_half, width);
    const auto w = static_cast<std::size_t>(width);
    ag::Oriented oriented = ag::orient(window, w, ag::orientation_of(strand(orientation1), strand(orientation2)));
    if (oriented.sum == 0.0) {
        add_one(agg_info, "empty_mat");
        return;
    }
    add_one(agg_info, "used_counter");
    const std::int64_t used = hicobj::int64_of(hicobj::item(agg_info, py::str("used_counter")));
    if (used % 50000 == 0) {
        log(kLogger, "info",
            "Number of used contacts within the given range: " + hicx::npy::int_with_thousands_separator(used));
    }
    if (equals(transform, "total-counts") && oriented.sum > 0.0) {
        for (double& value : oriented.values) {
            value = value / oriented.sum;
        }
    }
    const py::object totals = hicobj::item(hicobj::item(agg_info, py::str("agg_total")), chrom1);
    totals[chrom2] = py::int_(hicobj::int64_of(hicobj::item(totals, chrom2)) + 1);
    std::vector<double> diagonal(w);
    for (std::size_t i = 0; i < w; ++i) {
        diagonal[i] = oriented.values[i * w + i];
    }
    const auto group = [&](const char* name) {
        return hicobj::item(hicobj::item(hicobj::item(agg_info, py::str(name)), chrom1), chrom2);
    };
    group("agg_matrix").attr("append")(np_matrix(oriented.values, w, w));
    group("agg_diagonals").attr("append")(np_matrix(diagonal, 1, w));
    group("agg_center_values").attr("append")(view.scalar(view.at(bin1, bin2)));
    group("agg_contact_position").attr("append")(py::make_tuple(start1, end1, start2, end2));
}

// the cleanup at the end of aggregate_contacts
void drop_empty(py::handle agg_info) {
    py::list to_del;
    py::object agg_matrix = hicobj::item(agg_info, py::str("agg_matrix"));
    const py::list items(agg_matrix.attr("items")());
    for (const py::handle entry : items) {
        const py::object k1 = hicobj::item(entry, 0);
        const py::object v1 = hicobj::item(entry, 1);
        if (hicobj::py_equal(v1, py::dict())) {
            log(kLogger, "info", "no matrix on " + str_of(k1));
            to_del.append(k1);
            continue;
        }
        py::dict kept;
        for (const py::handle inner : py::list(v1.attr("items")())) {
            const py::object value = hicobj::item(inner, 1);
            if (!hicobj::py_equal(value, py::list())) {
                kept[hicobj::item(inner, 0)] = value;
            }
        }
        if (!hicobj::py_equal(kept, py::dict())) {
            agg_matrix[k1] = kept;
        } else {
            to_del.append(k1);
        }
    }
    for (const char* name : {"agg_matrix", "agg_contact_position"}) {
        py::dict filtered;
        for (const py::handle entry : py::list(hicobj::item(agg_info, py::str(name)).attr("items")())) {
            if (!hicobj::py_contains(to_del, hicobj::item(entry, 0))) {
                filtered[hicobj::item(entry, 0)] = hicobj::item(entry, 1);
            }
        }
        agg_info[py::str(name)] = filtered;
    }
}

std::vector<std::string> fields_of(py::handle line) {
    if (!is_str(line)) {
        raise_error("TypeError", "hicAggregateContacts: HiCExplorer 4 reads str lines");
    }
    const std::string text = utf8(line);
    std::vector<std::string> fields;
    for (const std::string_view field : hicx::pytext::split_whitespace(hicx::pytext::strip(text))) {
        fields.emplace_back(field);
    }
    return fields;
}

const std::string& field(const std::vector<std::string>& fields, std::size_t index) {
    if (index >= fields.size()) {
        raise_error("IndexError", "list index out of range");
    }
    return fields[index];
}

// --------------------------------------------------------------------------
// compute_avg and the clustering

struct Stack {
    std::vector<double> values;
    std::vector<py::ssize_t> shape;
    [[nodiscard]] std::size_t count() const { return shape.empty() ? 0 : static_cast<std::size_t>(shape[0]); }
    [[nodiscard]] std::size_t cell() const {
        std::size_t c = 1;
        for (std::size_t d = 1; d < shape.size(); ++d) {
            c *= static_cast<std::size_t>(shape[d]);
        }
        return c;
    }
};

Stack stack_of(py::handle submatrices, const char* function) {
    const py::array array = py::array::ensure(submatrices);
    if (!array) {
        throw_python_error();
    }
    if (dtype_of(array) != DType::Float64) {
        raise_unsupported_dtype(array, function);
    }
    Stack out;
    out.values = to_vector<double>(array);
    for (py::ssize_t d = 0; d < array.ndim(); ++d) {
        out.shape.push_back(array.shape(d));
    }
    return out;
}

py::object compute_avg(const Stack& stack, py::handle operationType) {
    if (stack.shape.empty()) {
        raise_error("numpy.exceptions.AxisError", "axis 0 is out of bounds for array of dimension 0");
    }
    const std::vector<py::ssize_t> out_shape(stack.shape.begin() + 1, stack.shape.end());
    const auto result = [&](const std::vector<double>& values) -> py::object {
        if (out_shape.empty()) {
            return numpy_scalar<double>(values.empty() ? std::nan("") : values[0]);
        }
        return float_array(values, out_shape);
    };
    if (equals(operationType, "median")) {
        const std::vector<double> median = ag::axis0(stack.values, stack.count(), stack.cell(), ag::Average::Median);
        const double sum = npy::pairwise_sum(median);
        if (sum == 0.0 || std::isnan(sum)) {
            const std::vector<double> mean = ag::axis0(stack.values, stack.count(), stack.cell(), ag::Average::Mean);
            if (npy::pairwise_sum(mean) != 0.0) {
                log(kLogger, "info", "The median of the matrices is zero. Consider using the mean instead.");
            } else {
                log(kLogger, "info", "Apparently no matrices could be computed. All of them are zeros or nans.");
            }
        }
        return result(median);
    }
    if (equals(operationType, "mean")) {
        return result(ag::axis0(stack.values, stack.count(), stack.cell(), ag::Average::Mean));
    }
    return result(ag::axis0(stack.values, stack.count(), stack.cell(), ag::Average::Sum));
}

// get_outlier_indices: None or a boolean array per row
py::object outlier_indices(const py::array& data, double max_deviation) {
    const py::array_t<double> values = py::array_t<double>::ensure(data);
    if (!values) {
        throw_python_error();
    }
    const std::vector<double> flat = to_vector<double>(values);
    const std::size_t rows = values.ndim() >= 1 ? static_cast<std::size_t>(values.shape(0)) : 1;
    const std::size_t cols = rows == 0 ? 0 : flat.size() / rows;
    const auto outliers = ag::outlier_rows(flat, rows, cols, max_deviation);
    if (!outliers) {
        return py::none();
    }
    if (values.ndim() != 2) {
        raise_error("numpy.exceptions.AxisError",
                    "axis 1 is out of bounds for array of dimension " + std::to_string(values.ndim()));
    }
    py::array_t<bool> out(static_cast<py::ssize_t>(rows));
    for (std::size_t r = 0; r < rows; ++r) {
        out.mutable_at(static_cast<py::ssize_t>(r)) = (*outliers)[r];
    }
    return out;
}

py::list compress(py::handle list, const std::vector<std::size_t>& kept) {
    py::list out;
    const std::size_t n = hicobj::length(list);
    for (const std::size_t k : kept) {
        if (k < n) {
            out.append(hicobj::item(list, static_cast<std::int64_t>(k)));
        }
    }
    return out;
}

py::object compute_clusters(py::handle updated_info, std::int64_t k, py::handle method, py::handle how,
                            py::handle max_deviation, py::handle keep_outlier) {
    const py::object submatrices = hicobj::item(updated_info, py::str("submatrices"));
    const py::array first = py::array::ensure(hicobj::item(submatrices, 0));
    if (!first) {
        throw_python_error();
    }
    if (first.ndim() != 2) {
        raise_error("TypeError", "compute_clusters: HiCExplorer 4 takes two dimensional submatrices");
    }
    const auto r = static_cast<std::size_t>(first.shape(0));
    const auto c = static_cast<std::size_t>(first.shape(1));
    const std::size_t center_bin = (r + 1) / 2;
    const bool diagonal = equals(how, "diagonal");
    const bool center = equals(how, "center");
    const std::size_t features = diagonal ? std::min(r, c) : (center ? 1 : r * c);
    const std::size_t n = hicobj::length(submatrices);
    std::vector<double> matrix(n * features);
    for (std::size_t s = 0; s < n; ++s) {
        const py::array sub = py::array::ensure(hicobj::item(submatrices, static_cast<std::int64_t>(s)));
        if (!sub || dtype_of(sub) != DType::Float64 || sub.ndim() != 2 ||
            static_cast<std::size_t>(sub.shape(0)) != r || static_cast<std::size_t>(sub.shape(1)) != c) {
            raise_error("ValueError", "compute_clusters: HiCExplorer 4 takes float64 submatrices of one shape");
        }
        const std::vector<double> m = to_vector<double>(sub);
        if (diagonal) {
            for (std::size_t i = 0; i < features; ++i) {
                matrix[s * features + i] = m[i * c + i];
            }
        } else if (center) {
            // submatrix[cb - 2:cb + 1, cb - 2:cb + 1].reshape((1, 9)).mean()
            const auto lo = static_cast<std::int64_t>(center_bin) - 2;
            const auto norm = [](std::int64_t index, std::size_t size) {
                const auto n_ = static_cast<std::int64_t>(size);
                if (index < 0) {
                    index += n_;
                }
                return std::clamp<std::int64_t>(index, 0, n_);
            };
            const std::int64_t row_lo = norm(lo, r);
            const std::int64_t row_hi = std::max(row_lo, norm(static_cast<std::int64_t>(center_bin) + 1, r));
            const std::int64_t col_lo = norm(lo, c);
            const std::int64_t col_hi = std::max(col_lo, norm(static_cast<std::int64_t>(center_bin) + 1, c));
            std::vector<double> block;
            for (std::int64_t i = row_lo; i < row_hi; ++i) {
                for (std::int64_t j = col_lo; j < col_hi; ++j) {
                    block.push_back(m[static_cast<std::size_t>(i) * c + static_cast<std::size_t>(j)]);
                }
            }
            if (block.size() != 9) {
                raise_error("ValueError", "cannot reshape array of size " + std::to_string(block.size()) +
                                              " into shape (1,9)");
            }
            matrix[s] = npy::pairwise_sum(block) / 9.0;
        } else {
            std::copy(m.begin(), m.end(), matrix.begin() + static_cast<std::ptrdiff_t>(s * features));
        }
    }

    py::object updated_coords = hicobj::item(updated_info, py::str("coords"));
    py::object updated_centers = hicobj::item(updated_info, py::str("centers"));
    py::object updated_submatrices = submatrices;
    py::object updated_diags = hicobj::item(updated_info, py::str("diagonal"));
    std::vector<std::size_t> kept(n);
    for (std::size_t s = 0; s < n; ++s) {
        kept[s] = s;
    }
    if (!hicobj::py_truthy(keep_outlier)) {
        const auto outliers = ag::outlier_rows(matrix, n, features, to_double(max_deviation));
        if (outliers) {
            std::vector<std::size_t> survivors;
            for (std::size_t s = 0; s < n; ++s) {
                if (!(*outliers)[s]) {
                    survivors.push_back(s);
                }
            }
            if (survivors.size() != n) {
                if (survivors.empty()) {
                    system_exit("ERROR: all submatrices have been detected as outliers. You can consider changing "
                                "the threshold");
                }
                kept = survivors;
                updated_coords = compress(hicobj::item(updated_info, py::str("coords")), kept);
                updated_centers = compress(hicobj::item(updated_info, py::str("centers")), kept);
                updated_submatrices = compress(submatrices, kept);
                updated_diags = compress(hicobj::item(updated_info, py::str("diagonal")), kept);
            }
        }
    }
    hicx::cluster::Samples samples;
    samples.samples = static_cast<std::int64_t>(kept.size());
    samples.features = static_cast<std::int64_t>(features);
    bool any_nan = false;
    for (const std::size_t s : kept) {
        for (std::size_t f = 0; f < features; ++f) {
            const double value = matrix[s * features + f];
            any_nan = any_nan || std::isnan(value);
            samples.values.push_back(std::isnan(value) ? 0.0 : value);
        }
    }
    if (any_nan) {
        log(kLogger, "warning", "For clustering nan values have to be replaced by zeros.");
    }
    std::optional<std::vector<std::int64_t>> labels;
    if (k == 1 && equals(method, "no_clust")) {
        labels = std::vector<std::int64_t>(kept.size(), 0);
    }
    try {
        if (equals(method, "kmeans")) {
            labels = hicx::cluster::kmeans(samples, k, 0);
        }
        if (equals(method, "hierarchical")) {
            labels = hicx::cluster::ward(samples, k);
        }
    } catch (const hicx::cluster::ClusteringError& error) {
        raise_error("ValueError", error.what());
    }
    if (equals(method, "spectral")) {
        raise_error("NotImplementedError",
                    "compute_clusters: HiCExplorer 4 has no spectral clustering (hicAggregateContacts never uses it)");
    }
    if (!labels) {
        unbound("cluster_labels");
    }
    py::list clustered_dict;
    for (std::int64_t cluster = 0; cluster < k; ++cluster) {
        std::vector<std::int64_t> members;
        for (std::size_t i = 0; i < labels->size(); ++i) {
            if ((*labels)[i] == cluster) {
                members.push_back(static_cast<std::int64_t>(i));
            }
        }
        clustered_dict.append(from_vector(members));
    }
    py::dict out;
    out["coords"] = updated_coords;
    out["centers"] = updated_centers;
    out["submatrices"] = updated_submatrices;
    out["clustered_dict"] = clustered_dict;
    out["diagonal"] = updated_diags;
    return out;
}

// '%0.5f' % value, as np.savetxt formats it
std::string fixed5(double value) {
    if (std::isnan(value)) {
        return "nan";
    }
    if (std::isinf(value)) {
        return value > 0 ? "inf" : "-inf";
    }
    char buffer[512];
    std::snprintf(buffer, sizeof(buffer), "%0.5f", value);
    return buffer;
}

void write_text(py::handle name, const std::string& content) {
    const std::string path = fs_path(name);
    std::FILE* handle = std::fopen(path.c_str(), "w");
    if (handle == nullptr) {
        raise_os_error(errno, name);
    }
    const std::size_t written = std::fwrite(content.data(), 1, content.size(), handle);
    const int error = errno;
    std::fclose(handle);
    if (written != content.size()) {
        raise_os_error(error, name);
    }
}

}  // namespace

void bind_hicAggregateContacts(py::module_& m) {
    m.def("read_bed_per_chrom", [](py::handle fh, py::handle chrom_list, py::handle pConsiderStrandDirection) {
        py::dict interval;
        const bool strand_wanted = hicobj::py_truthy(pConsiderStrandDirection);
        for (const py::handle line : fh) {
            if (!is_str(line)) {
                raise_error("TypeError", "read_bed_per_chrom: HiCExplorer 4 reads str lines");
            }
            const std::string text = utf8(line);
            if (text.empty()) {
                raise_error("IndexError", "string index out of range");
            }
            if (text[0] == '#') {
                continue;
            }
            std::vector<std::string> fields = fields_of(line);
            if (strand_wanted && fields.size() < 6) {
                log(kLogger, "error",
                    "Strand information should be considered, but BED file has not at least six columns. Exiting!");
                PyErr_SetObject(PyExc_SystemExit, py::int_(1).ptr());
                throw_python_error();
            }
            if (!hicobj::py_contains(chrom_list, to_str(field(fields, 0)))) {
                const std::string changed = hicx::utilities::change_chrom_name(fields[0]);
                if (hicobj::py_contains(chrom_list, to_str(changed))) {
                    fields[0] = changed;
                } else {
                    continue;
                }
            }
            const py::str chrom = to_str(fields[0]);
            if (!hicobj::py_contains(interval, chrom)) {
                interval[chrom] = py::list();
            }
            const py::object start = parse_int(to_str(field(fields, 1)));
            const py::object end = parse_int(to_str(field(fields, 2)));
            if (strand_wanted) {
                interval[chrom].attr("append")(py::make_tuple(start, end, to_str(field(fields, 5))));
            } else {
                interval[chrom].attr("append")(py::make_tuple(start, end));
            }
        }
        return interval;
    }, py::arg("fh"), py::arg("chrom_list"), py::arg("pConsiderStrandDirection") = py::none());

    m.def("aggregate_contacts", [](py::handle bed1, py::handle bed2, py::handle agg_info, py::handle ma,
                                   py::handle M_half, py::handle largeRegionsOperation, py::handle range,
                                   py::handle transform, py::handle mode, py::handle pConsiderStrandDirection) {
        Context context{ma, std::nullopt, {}};
        py::list seen_chrs;
        const bool strand_wanted = hicobj::py_truthy(pConsiderStrandDirection);
        const bool inter = equals(mode, "inter-chr");
        const bool intra = equals(mode, "intra-chr");
        for (const py::handle entry1 : py::list(bed1.attr("items")())) {
            const py::object k1 = hicobj::item(entry1, 0);
            const py::object v1 = hicobj::item(entry1, 1);
            for (const py::handle entry2 : py::list(bed2.attr("items")())) {
                const py::object k2 = hicobj::item(entry2, 0);
                const py::object v2 = hicobj::item(entry2, 1);
                const bool same = hicobj::py_equal(k1, k2);
                if (inter && same) {
                    if (hicobj::length(bed1) == 1 && hicobj::length(bed2) == 1) {
                        system_exit("Error: 'inter-chr' mode needs at least a pair of coordinates with different "
                                    "chromoses to be available in the bed files.");
                    }
                    continue;
                }
                if (intra && !same) {
                    continue;
                }
                for (const py::handle coord1 : v1) {
                    for (const py::handle coord2 : v2) {
                        if (same && hicobj::py_equal(coord1, coord2)) {
                            continue;
                        }
                        py::list interval;
                        if (strand_wanted) {
                            interval.append(py::make_tuple(k1, hicobj::item(coord1, 0), hicobj::item(coord1, 1),
                                                           hicobj::item(coord1, 2)));
                            interval.append(py::make_tuple(k2, hicobj::item(coord2, 0), hicobj::item(coord2, 1),
                                                           hicobj::item(coord2, 2)));
                        } else {
                            interval.append(py::make_tuple(k1, hicobj::item(coord1, 0), hicobj::item(coord1, 1)));
                            interval.append(py::make_tuple(k2, hicobj::item(coord2, 0), hicobj::item(coord2, 1)));
                        }
                        count_contacts(context, interval, M_half, mode, agg_info, largeRegionsOperation, seen_chrs,
                                       range, transform);
                    }
                }
            }
        }
        drop_empty(agg_info);
        return py::none();
    }, py::arg("bed1"), py::arg("bed2"), py::arg("agg_info"), py::arg("ma"), py::arg("M_half"),
       py::arg("largeRegionsOperation"), py::arg("range") = py::none(), py::arg("transform") = py::none(),
       py::arg("mode") = "", py::arg("pConsiderStrandDirection") = py::none());

    m.def("aggregate_contacts_per_row", [](py::handle bed1, py::handle bed2, py::handle agg_info, py::handle ma,
                                           py::handle chrom_list, py::handle M_half, py::handle largeRegionsOperation,
                                           py::handle range, py::handle transform, py::handle mode, py::handle perChr,
                                           py::handle pConsiderStrandDirection) {
        (void)perChr;
        Context context{ma, std::nullopt, {}};
        py::list seen_chrs;
        const bool strand_wanted = hicobj::py_truthy(pConsiderStrandDirection);
        const bool inter = equals(mode, "inter-chr");
        const bool intra = equals(mode, "intra-chr");
        py::object it1 = py::reinterpret_steal<py::object>(PyObject_GetIter(bed1.ptr()));
        py::object it2 = py::reinterpret_steal<py::object>(PyObject_GetIter(bed2.ptr()));
        if (!it1 || !it2) {
            throw_python_error();
        }
        while (true) {
            PyObject* raw1 = PyIter_Next(it1.ptr());
            if (raw1 == nullptr) {
                if (PyErr_Occurred()) {
                    throw_python_error();
                }
                break;
            }
            const py::object line1_object = py::reinterpret_steal<py::object>(raw1);
            PyObject* raw2 = PyIter_Next(it2.ptr());
            if (raw2 == nullptr) {
                if (PyErr_Occurred()) {
                    throw_python_error();
                }
                break;
            }
            const py::object line2_object = py::reinterpret_steal<py::object>(raw2);
            std::vector<std::string> line1 = fields_of(line1_object);
            std::vector<std::string> line2 = fields_of(line2_object);
            const auto known = [&](std::vector<std::string>& fields) {
                if (!hicobj::py_contains(chrom_list, to_str(field(fields, 0)))) {
                    fields[0] = utf8(utilities_change(fields[0]));
                    if (!hicobj::py_contains(chrom_list, to_str(fields[0]))) {
                        return false;
                    }
                }
                return true;
            };
            if (!known(line1) || !known(line2)) {
                continue;
            }
            if (inter && line1[0] == line2[0]) {
                continue;
            }
            if (!inter && intra && line1[0] != line2[0]) {
                continue;
            }
            py::list interval;
            if (strand_wanted) {
                interval.append(py::make_tuple(to_str(line1[0]), to_str(field(line1, 1)), to_str(field(line1, 2)),
                                               to_str(field(line1, 5))));
                interval.append(py::make_tuple(to_str(line2[0]), to_str(field(line2, 1)), to_str(field(line2, 2)),
                                               to_str(field(line2, 5))));
            } else {
                interval.append(py::make_tuple(to_str(line1[0]), to_str(field(line1, 1)), to_str(field(line1, 2))));
                interval.append(py::make_tuple(to_str(line2[0]), to_str(field(line2, 1)), to_str(field(line2, 2))));
            }
            count_contacts(context, interval, M_half, mode, agg_info, largeRegionsOperation, seen_chrs, range,
                           transform);
        }
        py::list to_del;
        for (const py::handle entry : py::list(hicobj::item(agg_info, py::str("agg_matrix")).attr("items")())) {
            const py::object k1 = hicobj::item(entry, 0);
            const py::object v1 = hicobj::item(entry, 1);
            if (hicobj::py_equal(v1, py::dict())) {
                log(kLogger, "info", "no matrix on " + str_of(k1));
                to_del.append(k1);
            } else if (intra && hicobj::py_equal(hicobj::item(v1, k1), py::list())) {
                log(kLogger, "info", "no matrix between " + str_of(k1) + " and " + str_of(k1));
                to_del.append(k1);
            }
        }
        for (const char* name : {"agg_matrix", "agg_contact_position"}) {
            py::dict filtered;
            for (const py::handle entry : py::list(hicobj::item(agg_info, py::str(name)).attr("items")())) {
                if (!hicobj::py_contains(to_del, hicobj::item(entry, 0))) {
                    filtered[hicobj::item(entry, 0)] = hicobj::item(entry, 1);
                }
            }
            agg_info[py::str(name)] = filtered;
        }
        return py::none();
    }, py::arg("bed1"), py::arg("bed2"), py::arg("agg_info"), py::arg("ma"), py::arg("chrom_list"), py::arg("M_half"),
       py::arg("largeRegionsOperation"), py::arg("range") = py::none(), py::arg("transform") = py::none(),
       py::arg("mode") = "", py::arg("perChr") = false, py::arg("pConsiderStrandDirection") = py::none());

    m.def("count_contacts", [](py::handle interval, py::handle ma, py::handle M_half, py::handle mode,
                               py::handle agg_info, py::handle largeRegionsOperation, py::handle seen_chrs,
                               py::handle range, py::handle transform) {
        Context context{ma, std::nullopt, {}};
        count_contacts(context, interval, M_half, mode, agg_info, largeRegionsOperation, seen_chrs, range, transform);
        return py::none();
    }, py::arg("interval"), py::arg("ma"), py::arg("M_half"), py::arg("mode"), py::arg("agg_info"),
       py::arg("largeRegionsOperation"), py::arg("seen_chrs"), py::arg("range") = py::none(),
       py::arg("transform") = py::none());

    m.def("get_outlier_indices", [](py::handle data, py::handle max_deviation) {
        const py::array array = py::array::ensure(data);
        if (!array) {
            throw_python_error();
        }
        return outlier_indices(array, to_double(max_deviation));
    }, py::arg("data"), py::arg("max_deviation") = 200);

    m.def("compute_clusters", [](py::handle updated_info, py::handle k, py::handle method, py::handle how,
                                 py::handle max_deviation, py::handle keep_outlier) {
        return compute_clusters(updated_info, to_int64(k), method, how, max_deviation, keep_outlier);
    }, py::arg("updated_info"), py::arg("k"), py::arg("method") = "kmeans", py::arg("how") = "full",
       py::arg("max_deviation") = 2, py::arg("keep_outlier") = false);

    m.def("cluster_matrices", [](py::handle agg_info, py::handle k_object, py::handle method, py::handle how,
                                 py::handle perChr, py::handle max_deviation, py::handle keep_outlier) {
        std::int64_t k = to_int64(k_object);
        py::dict updated_info;
        const bool per_chr = hicobj::py_truthy(perChr);
        if (!per_chr) {
            py::dict genome;
            for (const char* key : {"coords", "centers", "submatrices", "clustered_dict", "diagonal"}) {
                genome[key] = py::list();
            }
            updated_info["genome"] = genome;
        }
        const py::object agg_matrix = hicobj::item(agg_info, py::str("agg_matrix"));
        for (const py::handle chrom1 : py::list(agg_matrix.attr("keys")())) {
            const py::object per_chrom = hicobj::item(agg_matrix, chrom1);
            for (const py::handle entry : py::list(per_chrom.attr("items")())) {
                const py::object chrom2 = hicobj::item(entry, 0);
                const py::object matrices = hicobj::item(entry, 1);
                const auto group = [&](const char* name) {
                    return hicobj::item(hicobj::item(hicobj::item(agg_info, py::str(name)), chrom1), chrom2);
                };
                if (per_chr) {
                    py::list full_coords;
                    for (const py::handle i : group("agg_contact_position")) {
                        full_coords.append(py::make_tuple(chrom1, hicobj::item(i, 0), hicobj::item(i, 1), chrom2,
                                                          hicobj::item(i, 2), hicobj::item(i, 3)));
                    }
                    py::dict info;
                    info["coords"] = full_coords;
                    info["centers"] = group("agg_center_values");
                    info["submatrices"] = group("agg_matrix");
                    info["clustered_dict"] = py::list();
                    info["diagonal"] = group("agg_diagonals");
                    updated_info[chrom1] = info;
                    if (assertions_enabled() && !hicobj::py_equal(chrom1, chrom2)) {
                        raise_error("AssertionError", "");
                    }
                    const std::size_t count = hicobj::length(group("agg_matrix"));
                    log(kLogger, "info", "Length of entry on chr " + str_of(chrom1) + ": " + std::to_string(count));
                    if (static_cast<std::int64_t>(count) < k) {
                        log(kLogger, "info", "number of the submatrices on chromosome " + str_of(chrom1) +
                                                 " is less than " + std::to_string(k) + ". Clustering is skipped.");
                        k = 1;
                    }
                    updated_info[chrom1] =
                        compute_clusters(hicobj::item(updated_info, chrom1), k, method, how, max_deviation, keep_outlier);
                } else {
                    const py::object genome = hicobj::item(updated_info, py::str("genome"));
                    hicobj::item(genome, py::str("submatrices")).attr("extend")(group("agg_matrix"));
                    const py::object positions = group("agg_contact_position");
                    const py::object centers = group("agg_center_values");
                    const py::object diagonals = group("agg_diagonals");
                    const std::size_t count = hicobj::length(matrices);
                    for (std::size_t idx = 0; idx < count; ++idx) {
                        const py::object position = hicobj::item(positions, static_cast<std::int64_t>(idx));
                        if (hicobj::length(position) != 4) {
                            raise_error("ValueError", hicobj::length(position) < 4
                                                          ? "not enough values to unpack (expected 4, got " +
                                                                std::to_string(hicobj::length(position)) + ")"
                                                          : std::string("too many values to unpack (expected 4)"));
                        }
                        hicobj::item(genome, py::str("coords"))
                            .attr("append")(py::make_tuple(chrom1, hicobj::item(position, 0), hicobj::item(position, 1),
                                                           chrom2, hicobj::item(position, 2),
                                                           hicobj::item(position, 3)));
                        hicobj::item(genome, py::str("centers"))
                            .attr("append")(hicobj::item(centers, static_cast<std::int64_t>(idx)));
                        hicobj::item(genome, py::str("diagonal"))
                            .attr("append")(hicobj::item(diagonals, static_cast<std::int64_t>(idx)));
                    }
                }
            }
        }
        if (!per_chr) {
            updated_info["genome"] = compute_clusters(hicobj::item(updated_info, py::str("genome")), k, method, how,
                                                      max_deviation, keep_outlier);
        }
        return updated_info;
    }, py::arg("agg_info"), py::arg("k"), py::arg("method") = "kmeans", py::arg("how") = "full",
       py::arg("perChr") = false, py::arg("max_deviation") = 2, py::arg("keep_outlier") = false);

    m.def("compute_avg", [](py::handle submatrices, py::handle operationType) {
        return compute_avg(stack_of(submatrices, "compute_avg"), operationType);
    }, py::arg("submatrices"), py::arg("operationType"));

    // plot_aggregated_contacts: compute_avg of a cluster with its log records;
    // (average, empty) where empty means the Python's `continue`.
    m.def("_cluster_average", [](py::handle chrom_info, py::handle cluster_indices, py::handle operationType,
                                 py::handle cluster_number) {
        const py::object submatrices = hicobj::item(chrom_info, py::str("submatrices"));
        py::list members;
        for (const py::handle x : cluster_indices) {
            members.append(hicobj::item(submatrices, x));
        }
        const Stack stack = stack_of(numpy().attr("array")(members), "plot_aggregated_contacts");
        const py::object average = compute_avg(stack, operationType);
        const py::array values = py::array_t<double>::ensure(average);
        const std::vector<double> flat = to_vector<double>(values);
        const double mean = npy::pairwise_sum(flat) / static_cast<double>(flat.size());
        log(kLogger, "info", "Mean aggregate matrix values: " + npy::float_repr(mean));
        const std::int64_t number = to_int64(cluster_number) + 1;
        log(kLogger, "info", "total pairs considered on cluster_" + std::to_string(number) + ": " +
                                 std::to_string(hicobj::length(cluster_indices)));
        if (values.ndim() == 0) {
            raise_error("IndexError", "tuple index out of range");
        }
        bool empty = false;
        if (values.shape(0) == 0) {
            log(kLogger, "debug", "matrix for contacts on cluster_" + std::to_string(number) + " is empty");
            empty = true;
        }
        return py::make_tuple(average, py::bool_(empty));
    }, py::arg("chrom_info"), py::arg("cluster_indices"), py::arg("operationType"), py::arg("cluster_number"));

    // plot_aggregated_contacts: the --outFilePrefixMatrix and
    // --outFileContactPairs files of a cluster.
    m.def("_cluster_outputs", [](py::handle chrom_info, py::handle chrom1, py::handle cluster_number,
                                 py::handle cluster_indices, py::handle num_clusters, py::handle average,
                                 py::handle args) {
        const std::int64_t clusters = to_int64(num_clusters);
        const std::string id = std::to_string(to_int64(cluster_number) + 1);
        const auto name_of = [&](py::handle prefix) {
            if (clusters == 1) {
                return str_of(prefix) + "_" + str_of(chrom1) + ".tab";
            }
            return str_of(prefix) + "_" + str_of(chrom1) + "_cluster_" + id + ".tab";
        };
        const py::object prefix_matrix = args.attr("outFilePrefixMatrix");
        if (hicobj::py_truthy(prefix_matrix)) {
            const py::array values = py::array_t<double>::ensure(average);
            const std::vector<double> flat = to_vector<double>(values);
            const std::size_t cols = values.ndim() == 2 ? static_cast<std::size_t>(values.shape(1)) : 1;
            std::string text;
            for (std::size_t k = 0; k < flat.size(); ++k) {
                text += fixed5(flat[k]);
                text += (k + 1) % cols == 0 ? "\n" : "\t";
            }
            write_text(to_str(name_of(prefix_matrix)), text);
        }
        const py::object pairs = args.attr("outFileContactPairs");
        if (hicobj::py_truthy(pairs)) {
            const py::array centers_all = py::array::ensure(hicobj::item(chrom_info, py::str("centers")));
            if (!centers_all) {
                throw_python_error();
            }
            const std::vector<std::int64_t> indices = index_vector_of(cluster_indices);
            const bool integer = dtype_of(centers_all) != DType::Float64 && dtype_of(centers_all) != DType::Float32;
            const py::array_t<double> centers_double = py::array_t<double>::ensure(centers_all);
            const std::vector<double> all = to_vector<double>(centers_double);
            std::vector<double> chosen;
            std::vector<std::int64_t> chosen_int;
            for (const std::int64_t index : indices) {
                const std::int64_t n = static_cast<std::int64_t>(all.size());
                const std::int64_t i = index < 0 ? index + n : index;
                if (i < 0 || i >= n) {
                    raise_error("IndexError", "index " + std::to_string(index) + " is out of bounds for axis 0 with size " +
                                                  std::to_string(n));
                }
                chosen.push_back(all[static_cast<std::size_t>(i)]);
                chosen_int.push_back(static_cast<std::int64_t>(all[static_cast<std::size_t>(i)]));
            }
            const std::vector<std::int64_t> order = integer ? npy::argsort(chosen_int) : npy::argsort(chosen);
            const py::object coords = hicobj::item(chrom_info, py::str("coords"));
            std::string text;
            for (auto it = order.rbegin(); it != order.rend(); ++it) {
                const auto cl_idx = static_cast<std::size_t>(*it);
                const py::object coord = hicobj::item(coords, static_cast<std::int64_t>(cl_idx));
                if (hicobj::length(coord) != 6) {
                    raise_error("ValueError", "not enough values to unpack (expected 6, got " +
                                                  std::to_string(hicobj::length(coord)) + ")");
                }
                for (std::int64_t f = 0; f < 6; ++f) {
                    text += str_of(hicobj::item(coord, f));
                    text += "\t";
                }
                text += integer ? std::to_string(chosen_int[cl_idx]) : npy::float_repr(chosen[cl_idx]);
                text += "\n";
            }
            write_text(to_str(name_of(pairs)), text);
        }
        return py::none();
    }, py::arg("chrom_info"), py::arg("chrom1"), py::arg("cluster_number"), py::arg("cluster_indices"),
       py::arg("num_clusters"), py::arg("average"), py::arg("args"));

    // plot_diagnostic_heatmaps: None (after its log record) when the
    // diagonals do not stack, else (cluster sizes, [(rows sorted by their
    // centre value, descending, their column means, x positions, rows)]).
    m.def("_diagnostic_data", [](py::handle chrom_info, py::handle chrom_name, py::handle M_half) -> py::object {
        const py::object diagonals = hicobj::item(chrom_info, py::str("diagonal"));
        std::vector<double> heatmap;
        std::size_t width = 0;
        const std::size_t rows = hicobj::length(diagonals);
        bool stackable = rows > 0;
        for (std::size_t r = 0; r < rows && stackable; ++r) {
            const py::array row = py::array_t<double>::ensure(hicobj::item(diagonals, static_cast<std::int64_t>(r)));
            if (!row) {
                throw_python_error();
            }
            const std::size_t size = static_cast<std::size_t>(row.size());
            if (r == 0) {
                width = size;
            } else if (size != width) {
                stackable = false;
                break;
            }
            const std::vector<double> values = to_vector<double>(row);
            heatmap.insert(heatmap.end(), values.begin(), values.end());
        }
        if (!stackable) {
            log(kLogger, "error", "Error computing diagnostic heatmap for chrom: " + str_of(chrom_name));
            return py::none();
        }
        const py::object clustered = hicobj::item(chrom_info, py::str("clustered_dict"));
        py::list sizes;
        py::list clusters;
        const std::int64_t half = to_int64(M_half);
        for (const py::handle cluster_indices : clustered) {
            sizes.append(py::int_(hicobj::length(cluster_indices)));
            const std::vector<std::int64_t> members = index_vector_of(cluster_indices);
            if (half < 0 || static_cast<std::size_t>(half) >= width) {
                raise_error("IndexError", "index " + std::to_string(half) + " is out of bounds for axis 1 with size " +
                                              std::to_string(width));
            }
            for (const std::int64_t member : members) {
                if (member < 0 || static_cast<std::size_t>(member) >= rows) {
                    raise_error("IndexError", "index " + std::to_string(member) +
                                                  " is out of bounds for axis 0 with size " + std::to_string(rows));
                }
            }
            const ag::DiagnosticCluster cluster =
                ag::diagnostic_cluster(heatmap, width, members, static_cast<std::size_t>(half));
            std::vector<double> sorted;
            for (const std::size_t row : cluster.rows) {
                sorted.insert(sorted.end(), heatmap.begin() + static_cast<std::ptrdiff_t>(row * width),
                              heatmap.begin() + static_cast<std::ptrdiff_t>((row + 1) * width));
            }
            const std::size_t count = cluster.rows.size();
            const std::vector<double>& mean = cluster.mean;
            std::vector<std::int64_t> x(width);
            for (std::size_t c = 0; c < width; ++c) {
                x[c] = static_cast<std::int64_t>(c) - half;
            }
            clusters.append(py::make_tuple(
                float_array(sorted, {static_cast<py::ssize_t>(count), static_cast<py::ssize_t>(width)}),
                float_array(mean, {static_cast<py::ssize_t>(width)}), from_vector(x), py::int_(count)));
        }
        return py::make_tuple(sizes, clusters);
    }, py::arg("chrom_info"), py::arg("chrom_name"), py::arg("M_half"));
}

}  // namespace hicx::python
