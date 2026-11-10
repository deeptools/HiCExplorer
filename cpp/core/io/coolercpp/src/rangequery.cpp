// Port of cooler/core/_rangequery.py (cooler 0.10.2, BSD-3-Clause).

#include "rangequery.hpp"

#include <algorithm>
#include <limits>

#include "numpy_compat.hpp"

namespace coolercpp::detail {

namespace {

constexpr std::int64_t kInt32Min = std::numeric_limits<std::int32_t>::min();
constexpr std::int64_t kInt32Max = std::numeric_limits<std::int32_t>::max();

std::string join(const std::string& group, const std::string& name) {
    return group == "/" ? "/" + name : group + "/" + name;
}

// numpy floor division of int64 by a Python int (division by zero gives 0
// with a RuntimeWarning in numpy).
std::int64_t floordiv(std::int64_t a, std::int64_t b) {
    if (b == 0) {
        return 0;
    }
    std::int64_t q = a / b;
    if ((a % b != 0) && ((a < 0) != (b < 0))) {
        --q;
    }
    return q;
}

bool comes_before(std::int64_t a0, std::int64_t a1, std::int64_t b0, std::int64_t b1,
                  bool strict) {
    if (a0 < b0) {
        return strict ? a1 <= b0 : a1 <= b1;
    }
    return false;
}

bool contains(std::int64_t a0, std::int64_t a1, std::int64_t b0, std::int64_t b1) {
    if (a0 > b0 || a1 < b1) {
        return false;
    }
    return a0 <= b0 && a1 >= b1;
}

bool fits_int32(std::int64_t v) { return v >= kInt32Min && v <= kInt32Max; }

template <typename Id>
std::vector<Id>& ids(IdArray& array) {
    if constexpr (std::is_same_v<Id, std::int32_t>) {
        return array.narrow_values();
    } else {
        return array.wide_values();
    }
}

}  // namespace

PixelDict run_tasks(const CSRReader& reader, const std::string& field, const Bbox& query,
                    const std::vector<Task>& tasks, bool reflect, bool return_index) {
    if (tasks.empty()) {
        return reader.meta(field, return_index);
    }
    bool narrow = reader.narrow_ids(query);
    std::size_t capacity = 0;
    for (const Task& task : tasks) {
        narrow = narrow && reader.narrow_ids(task.bbox);
        capacity += reader.span_pixels(task.span) * (reflect ? 2 : 1);
    }
    PixelDict out;
    out.bin1 = IdArray(narrow);
    out.bin2 = IdArray(narrow);
    out.bin1.reserve(capacity);
    out.bin2.reserve(capacity);
    out.value = Column::empty(reader.dtype_of(field));
    std::visit(
        [capacity](auto& v) {
            using V = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<V, CategoricalData>) {
                v.codes.reserve(capacity);
            } else {
                v.reserve(capacity);
            }
        },
        out.value.data());
    if (return_index) {
        out.index.emplace();
        out.index->reserve(capacity);
    }
    for (const Task& task : tasks) {
        reader.read_into(out, field, task.bbox, task.span, reflect, return_index, task.transpose);
    }
    return out;
}

CSRReader::CSRReader(const h5::File& file, std::string pixels_group,
                     std::vector<std::int64_t> bin1_offsets)
    : file_(file), pixels_(std::move(pixels_group)), offsets_(std::move(bin1_offsets)) {
    for (const std::string& column : file_.keys(pixels_)) {
        const std::string path = join(pixels_, column);
        if (file_.is_dataset(path)) {
            dtypes_.emplace_back(column, file_.open_dataset(path).type().dtype);
        }
    }
}

DType CSRReader::dtype_of(const std::string& column) const {
    for (const auto& [name, dtype] : dtypes_) {
        if (name == column) {
            return dtype;
        }
    }
    throw KeyError(column);
}

std::vector<std::pair<std::int64_t, std::int64_t>> CSRReader::get_spans(
    const Bbox& bbox, std::int64_t chunksize) const {
    std::vector<std::pair<std::int64_t, std::int64_t>> spans;
    if (bbox.i1 - bbox.i0 < 1 || bbox.j1 - bbox.j0 < 1) {
        return spans;
    }
    // edges = i0 + arg_prune_partition(bin1_offsets[i0 : i1 + 1], chunksize)
    const auto n = static_cast<std::int64_t>(offsets_.size());
    const npy::SliceRange range = npy::slice_indices(bbox.i0, bbox.i1 + 1, n);
    if (range.size() == 0) {
        throw IndexError("index 0 is out of bounds for axis 0 with size 0");
    }
    const auto first = offsets_.begin() + range.start;
    const auto last = offsets_.begin() + range.stop;
    const std::int64_t lo = *first;
    const std::int64_t hi = *(last - 1);
    const std::int64_t num = 2 + floordiv(hi - lo, chunksize);
    const std::vector<std::int64_t> cuts = npy::linspace_int(lo, hi, num);
    std::vector<std::int64_t> positions;
    positions.reserve(cuts.size());
    for (const std::int64_t cut : cuts) {
        positions.push_back(static_cast<std::int64_t>(std::lower_bound(first, last, cut) - first));
    }
    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
    for (std::size_t k = 1; k < positions.size(); ++k) {
        spans.emplace_back(bbox.i0 + positions[k - 1], bbox.i0 + positions[k]);
    }
    return spans;
}

PixelDict CSRReader::meta(const std::string& field, bool return_index) const {
    PixelDict dict;
    dict.bin1_dtype = dtype_of("bin1_id");
    dict.bin2_dtype = dtype_of("bin2_id");
    dict.value = Column::empty(dtype_of(field));
    if (return_index) {
        dict.index.emplace();
    }
    return dict;
}

std::size_t CSRReader::span_pixels(std::pair<std::int64_t, std::int64_t> span) const {
    const auto n = static_cast<std::int64_t>(offsets_.size());
    if (span.first < -n || span.first >= n || span.second < -n || span.second >= n) {
        return 0;
    }
    const std::int64_t lo = offsets_[static_cast<std::size_t>(npy::wrap_index(span.first, n))];
    const std::int64_t hi = offsets_[static_cast<std::size_t>(npy::wrap_index(span.second, n))];
    return hi > lo ? static_cast<std::size_t>(hi - lo) : 0;
}

bool CSRReader::narrow_ids(const Bbox& bbox) const {
    return offsets_.size() <= static_cast<std::size_t>(kInt32Max) && fits_int32(bbox.i0) &&
           fits_int32(bbox.i1) && fits_int32(bbox.j0) && fits_int32(bbox.j1);
}

void CSRReader::read_into(PixelDict& out, const std::string& field, const Bbox& bbox,
                          std::pair<std::int64_t, std::int64_t> span, bool reflect,
                          bool return_index, bool transpose) const {
    const std::int64_t s0 = span.first;
    const std::int64_t s1 = span.second;
    const auto n_offsets = static_cast<std::int64_t>(offsets_.size());
    const std::int64_t offset_lo = offsets_[static_cast<std::size_t>(npy::wrap_index(s0, n_offsets))];
    const std::int64_t offset_hi = offsets_[static_cast<std::size_t>(npy::wrap_index(s1, n_offsets))];

    const DType bin1_dtype = dtype_of("bin1_id");
    const h5::Dataset bin2_ds = file_.open_dataset(join(pixels_, "bin2_id"));
    const h5::Dataset data_ds = file_.open_dataset(join(pixels_, field));
    const DType bin2_dtype = bin2_ds.type().dtype;

    // The dtypes this part contributes to the concatenation.
    DType part1 = bin1_dtype;
    DType part2 = bin2_dtype;
    if (reflect && s1 > s0) {
        part1 = part2 = result_type(bin1_dtype, bin2_dtype);
    }
    if (transpose) {
        std::swap(part1, part2);
    }
    if (out.has_parts) {
        out.bin1_dtype = result_type(out.bin1_dtype, part1);
        out.bin2_dtype = result_type(out.bin2_dtype, part2);
    } else {
        out.bin1_dtype = part1;
        out.bin2_dtype = part2;
        out.has_parts = true;
    }
    if (s1 <= s0) {
        return;
    }

    const npy::SliceRange r2 =
        npy::slice_indices(offset_lo, offset_hi, static_cast<std::int64_t>(bin2_ds.length()));
    const npy::SliceRange rd =
        npy::slice_indices(offset_lo, offset_hi, static_cast<std::int64_t>(data_ds.length()));
    const Column data_extracted = data_ds.read_column(static_cast<std::size_t>(rd.start),
                                                      static_cast<std::size_t>(rd.stop));

    const auto run = [&]<typename Id>(Id /*tag*/) {
        const std::vector<Id> bin2_extracted = bin2_ds.read<Id>(
            static_cast<std::size_t>(r2.start), static_cast<std::size_t>(r2.stop));
        const auto extracted = static_cast<std::int64_t>(bin2_extracted.size());
        std::vector<Id>& rows = ids<Id>(transpose ? out.bin2 : out.bin1);
        std::vector<Id>& cols = ids<Id>(transpose ? out.bin1 : out.bin2);
        std::visit(
            [&](auto& dst) {
                using V = std::decay_t<decltype(dst)>;
                if constexpr (std::is_same_v<V, CategoricalData> ||
                              std::is_same_v<V, std::vector<std::string>>) {
                    throw TypeError("pixel value columns must be numeric");
                } else {
                    const V& src = std::get<V>(data_extracted.data());
                    const std::size_t first = rows.size();
                    for (std::int64_t i = s0; i < s1; ++i) {
                        const std::int64_t lo =
                            offsets_[static_cast<std::size_t>(npy::wrap_index(i, n_offsets))] - offset_lo;
                        const std::int64_t hi =
                            offsets_[static_cast<std::size_t>(npy::wrap_index(i + 1, n_offsets))] -
                            offset_lo;
                        const npy::SliceRange row = npy::slice_indices(lo, hi, extracted);
                        for (std::int64_t k = row.start; k < row.stop; ++k) {
                            const auto col = static_cast<std::int64_t>(bin2_extracted[static_cast<std::size_t>(k)]);
                            if (col >= bbox.j0 && col < bbox.j1) {
                                rows.push_back(static_cast<Id>(i));
                                cols.push_back(static_cast<Id>(col));
                                dst.push_back(src[static_cast<std::size_t>(k)]);
                                if (return_index) {
                                    out.index->push_back(offset_lo + k);
                                }
                            }
                        }
                    }
                    if (reflect) {
                        const std::size_t last = rows.size();
                        for (std::size_t k = first; k < last; ++k) {
                            const auto r = static_cast<std::int64_t>(rows[k]);
                            const auto c = static_cast<std::int64_t>(cols[k]);
                            if (r != c && c < bbox.i1) {
                                rows.push_back(static_cast<Id>(c));
                                cols.push_back(static_cast<Id>(r));
                                dst.push_back(dst[k]);
                                if (return_index) {
                                    out.index->push_back((*out.index)[k]);
                                }
                            }
                        }
                    }
                }
            },
            out.value.data());
    };
    if (out.bin1.narrow()) {
        run(std::int32_t{});
    } else {
        run(std::int64_t{});
    }
}

std::vector<Task> direct_tasks(const CSRReader& reader, const Bbox& bbox,
                               std::int64_t chunksize) {
    std::vector<Task> tasks;
    for (const auto& span : reader.get_spans(bbox, chunksize)) {
        tasks.push_back(Task{bbox, span, false});
    }
    return tasks;
}

std::vector<Task> fill_lower_tasks(const CSRReader& reader, const Bbox& bbox,
                                   std::int64_t chunksize) {
    std::int64_t i0 = bbox.i0;
    std::int64_t i1 = bbox.i1;
    std::int64_t j0 = bbox.j0;
    std::int64_t j1 = bbox.j1;
    const bool use_transpose = i1 > j1;
    if (use_transpose) {
        std::swap(i0, j0);
        std::swap(i1, j1);
    }
    std::vector<Bbox> boxes;
    std::vector<bool> transposed;
    if (i0 == j0 || comes_before(i0, i1, j0, j1, true)) {
        boxes = {Bbox{i0, i1, j0, j1}};
        transposed = {use_transpose};
    } else if (comes_before(i0, i1, j0, j1, false)) {
        boxes = {Bbox{i0, j0, j0, j1}, Bbox{j0, i1, j0, j1}};
        transposed = {use_transpose, use_transpose};
    } else if (contains(j0, j1, i0, i1)) {
        boxes = {Bbox{j0, i0, i0, i1}, Bbox{i0, i1, i0, j1}};
        transposed = {!use_transpose, use_transpose};
    } else {
        throw ValueError("This shouldn't happen.");
    }
    std::vector<Task> tasks;
    for (std::size_t b = 0; b < boxes.size(); ++b) {
        for (const auto& span : reader.get_spans(boxes[b], chunksize)) {
            tasks.push_back(Task{boxes[b], span, transposed[b]});
        }
    }
    return tasks;
}

PixelDict direct_query(const CSRReader& reader, const std::string& field, const Bbox& bbox,
                       std::int64_t chunksize, bool return_index) {
    return run_tasks(reader, field, bbox, direct_tasks(reader, bbox, chunksize), false,
                     return_index);
}

PixelDict fill_lower_query(const CSRReader& reader, const std::string& field, const Bbox& bbox,
                           std::int64_t chunksize, bool return_index) {
    return run_tasks(reader, field, bbox, fill_lower_tasks(reader, bbox, chunksize), true,
                     return_index);
}

namespace {

// Moves an ID array into a Column of the requested numpy dtype, reusing the
// buffer when the storage already has that type.
Column take_ids(IdArray& array, DType dtype) {
    if (array.narrow()) {
        if (dtype == DType::Int32) {
            return Column(std::move(array.narrow_values()));
        }
        Column out = Column(std::vector<std::int32_t>(std::move(array.narrow_values()))).astype(dtype);
        array.release();
        return out;
    }
    if (dtype == DType::Int64) {
        return Column(std::move(array.wide_values()));
    }
    Column out = Column(std::vector<std::int64_t>(std::move(array.wide_values()))).astype(dtype);
    array.release();
    return out;
}

// scipy coo_matrix((data, (bin1 - i0, bin2 - j0)), shape): the index dtype
// scipy picks and the checks it runs. Returns whether int64 indices are used.
bool check_coo(const PixelDict& dict, const Bbox& bbox) {
    const std::int64_t rows = bbox.i1 - bbox.i0;
    const std::int64_t cols = bbox.j1 - bbox.j0;
    if (rows < 0 || cols < 0) {
        throw ValueError("'shape' elements cannot be negative");
    }
    const auto minmax = [](const IdArray& ids, std::int64_t shift) {
        std::int64_t mn = std::numeric_limits<std::int64_t>::max();
        std::int64_t mx = std::numeric_limits<std::int64_t>::min();
        for (std::size_t k = 0; k < ids.size(); ++k) {
            const std::int64_t v = ids[k] - shift;
            mn = std::min(mn, v);
            mx = std::max(mx, v);
        }
        return std::pair<std::int64_t, std::int64_t>{mn, mx};
    };
    const auto [row_min, row_max] = minmax(dict.bin1, bbox.i0);
    const auto [col_min, col_max] = minmax(dict.bin2, bbox.j0);
    bool int64 = std::max(rows, cols) > kInt32Max;
    if (dict.bin1.size() > 0 &&
        (!fits_int32(row_min) || !fits_int32(row_max) || !fits_int32(col_min) || !fits_int32(col_max))) {
        int64 = true;
    }
    if (dict.value.size() > 0) {
        const std::array<std::int64_t, 2> shape{rows, cols};
        const std::array<std::pair<std::int64_t, std::int64_t>, 2> extremes{
            std::pair{row_min, row_max}, std::pair{col_min, col_max}};
        for (std::size_t axis = 0; axis < 2; ++axis) {
            if (extremes[axis].second >= shape[axis]) {
                throw ValueError("axis " + std::to_string(axis) + " index " +
                                 std::to_string(extremes[axis].second) + " exceeds matrix dimension " +
                                 std::to_string(shape[axis]));
            }
            if (extremes[axis].first < 0) {
                throw ValueError("negative axis " + std::to_string(axis) +
                                 " index: " + std::to_string(extremes[axis].first));
            }
        }
    }
    return int64;
}

Column shifted_indices(IdArray& ids, std::int64_t shift, bool int64) {
    if (!int64) {
        if (ids.narrow()) {
            std::vector<std::int32_t>& v = ids.narrow_values();
            for (std::int32_t& x : v) {
                x = static_cast<std::int32_t>(static_cast<std::int64_t>(x) - shift);
            }
            return Column(std::move(v));
        }
        const std::vector<std::int64_t>& v = ids.wide_values();
        std::vector<std::int32_t> out(v.size());
        for (std::size_t k = 0; k < v.size(); ++k) {
            out[k] = static_cast<std::int32_t>(v[k] - shift);
        }
        ids.release();
        return Column(std::move(out));
    }
    if (!ids.narrow()) {
        std::vector<std::int64_t>& v = ids.wide_values();
        for (std::int64_t& x : v) {
            x -= shift;
        }
        return Column(std::move(v));
    }
    const std::vector<std::int32_t>& v = ids.narrow_values();
    std::vector<std::int64_t> out(v.size());
    for (std::size_t k = 0; k < v.size(); ++k) {
        out[k] = static_cast<std::int64_t>(v[k]) - shift;
    }
    ids.release();
    return Column(std::move(out));
}

}  // namespace

Table frame_from_dict(PixelDict&& dict, const std::string& field) {
    Table table;
    const std::size_t n = dict.bin1.size();
    table.set("bin1_id", take_ids(dict.bin1, dict.bin1_dtype));
    table.set("bin2_id", take_ids(dict.bin2, dict.bin2_dtype));
    table.set(field, std::move(dict.value));
    if (dict.index.has_value()) {
        table.set_index(Index::values(std::move(*dict.index)));
    } else {
        table.set_index(Index::range(0, n));
    }
    return table;
}

SparseMatrix sparse_from_dict(PixelDict&& dict, const Bbox& bbox) {
    const bool int64 = check_coo(dict, bbox);
    SparseMatrix out;
    out.rows = bbox.i1 - bbox.i0;
    out.cols = bbox.j1 - bbox.j0;
    out.row = shifted_indices(dict.bin1, bbox.i0, int64);
    out.col = shifted_indices(dict.bin2, bbox.j0, int64);
    out.data = std::move(dict.value);
    return out;
}

DenseMatrix dense_from_dict(const PixelDict& dict, const Bbox& bbox) {
    (void)check_coo(dict, bbox);
    DenseMatrix out;
    out.rows = bbox.i1 - bbox.i0;
    out.cols = bbox.j1 - bbox.j0;
    const auto total = static_cast<std::size_t>(out.rows) * static_cast<std::size_t>(out.cols);
    out.values = Column::empty(dict.value.dtype(), total);
    const auto cols = static_cast<std::size_t>(out.cols);
    std::visit(
        [&](auto& dst) {
            using V = std::decay_t<decltype(dst)>;
            if constexpr (std::is_same_v<V, CategoricalData> ||
                          std::is_same_v<V, std::vector<std::string>>) {
                throw TypeError("a string column cannot fill a matrix");
            } else {
                using T = typename V::value_type;
                const auto& src = std::get<V>(dict.value.data());
                for (std::size_t k = 0; k < src.size(); ++k) {
                    const auto r = static_cast<std::size_t>(dict.bin1[k] - bbox.i0);
                    const auto c = static_cast<std::size_t>(dict.bin2[k] - bbox.j0);
                    auto& cell = dst[r * cols + c];
                    if constexpr (std::is_same_v<T, Bool8>) {
                        cell = Bool8(cell || src[k]);
                    } else if constexpr (std::is_integral_v<T>) {
                        // Wrapping addition, as numpy does for integer overflow.
                        using U = std::make_unsigned_t<T>;
                        cell = static_cast<T>(static_cast<U>(cell) + static_cast<U>(src[k]));
                    } else {
                        cell += src[k];
                    }
                }
            }
        },
        out.values.data());
    return out;
}

}  // namespace coolercpp::detail
