#include "hicx/h5_file.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "hicx/hdf5_util.hpp"

namespace hicx {

namespace {

// One block of a streamed column. 65,536 elements is 512 kB for a float64
// column, small enough to stay a rounding error against the budget of
// cpp/PLAN.md 4.5 and large enough that the per call HDF5 overhead disappears.
constexpr std::size_t kBlock = 65536;

// scipy picks the narrowest index type that can address the matrix
// (scipy.sparse.sputils.get_index_dtype), and the h5 file stores whatever
// scipy chose, so the writer has to make the same choice.
bool index_fits_in_int32(std::size_t nnz, std::int64_t rows) {
    const std::int64_t limit = 2147483647;
    return static_cast<std::int64_t>(nnz) <= limit && rows <= limit;
}

}  // namespace

bool is_hicexplorer_h5(const std::string& path) {
    if (!h5::is_hdf5(path)) {
        return false;
    }
    try {
        const h5::File file(path);
        return file.exists("/matrix/data") && file.exists("/matrix/indices") &&
               file.exists("/matrix/indptr") && file.exists("/intervals/chr_list");
    } catch (const h5::Error&) {
        return false;
    }
}

H5MatrixData read_hicexplorer_h5(const std::string& path) {
    const h5::File file(path);

    const std::vector<std::int64_t> shape = file.read_int64("/matrix/shape");
    if (shape.size() != 2) {
        throw h5::Error("unexpected /matrix/shape in " + path);
    }
    const std::string dtype = file.dataset_dtype("/matrix/data");
    std::vector<double> data = file.read_doubles("/matrix/data");
    std::vector<std::int32_t> indices = file.read_int32("/matrix/indices");
    std::vector<std::int64_t> indptr = file.read_int64("/matrix/indptr");

    H5MatrixData result;
    result.matrix = CsrMatrix(shape[0], shape[1], std::move(indptr), std::move(indices),
                              std::move(data), dtype);

    const std::vector<std::string> chroms = file.read_strings("/intervals/chr_list");
    const std::vector<std::int64_t> starts = file.read_int64("/intervals/start_list");
    const std::vector<std::int64_t> ends = file.read_int64("/intervals/end_list");
    // extra_list is normally float64, but matrices written by hicFindTADs
    // store text there (for example the z-score matrices under
    // test_data/find_TADs). Both have to be readable.
    const bool extra_is_text = file.dataset_dtype("/intervals/extra_list") == "string";
    std::vector<double> extra;
    std::vector<std::string> extra_text;
    if (extra_is_text) {
        extra_text = file.read_strings("/intervals/extra_list");
    } else {
        extra = file.read_doubles("/intervals/extra_list");
    }
    const std::size_t extra_size = extra_is_text ? extra_text.size() : extra.size();
    if (chroms.size() != starts.size() || chroms.size() != ends.size() ||
        chroms.size() != extra_size) {
        throw h5::Error("the interval lists of " + path + " have different lengths");
    }
    result.cut_intervals.reserve(chroms.size());
    for (std::size_t i = 0; i < chroms.size(); ++i) {
        CutInterval interval{chroms[i], starts[i], ends[i], 0.0, {}};
        if (extra_is_text) {
            interval.extra = std::numeric_limits<double>::quiet_NaN();
            interval.extra_text = extra_text[i];
        } else {
            interval.extra = extra[i];
        }
        result.cut_intervals.push_back(std::move(interval));
    }
    if (static_cast<std::int64_t>(result.cut_intervals.size()) != shape[0]) {
        throw h5::Error("Error loading matrix. Length of bin intervals (" +
                        std::to_string(result.cut_intervals.size()) +
                        ") is different than the size of the matrix (" +
                        std::to_string(shape[0]) + ")");
    }

    if (file.exists("/nan_bins")) {
        result.nan_bins = file.read_int64("/nan_bins");
    }

    if (file.exists("/correction_factors")) {
        result.correction_factors_are_column =
            file.dataset_is_column("/correction_factors");
        std::vector<double> factors = file.read_doubles("/correction_factors");
        if (static_cast<std::int64_t>(factors.size()) != shape[0]) {
            throw h5::Error(
                "Error loading matrix. Length of correction factors does not"
                "match size of matrix");
        }
        for (double& value : factors) {
            if (std::isnan(value) || std::isinf(value)) {
                value = 0.0;
            }
        }
        result.correction_factors = std::move(factors);
    }

    // The Python reader has a copy and paste defect here: when
    // /distance_counts exists it reads /correction_factors instead. It is
    // reproduced so that the two implementations agree bit for bit.
    if (file.exists("/distance_counts") && file.exists("/correction_factors")) {
        result.distance_counts = file.read_doubles("/correction_factors");
    }

    return result;
}

namespace {

// Writes a numeric column that is produced element by element. `produce` is
// called with a callback it invokes once per value; the values are buffered in
// blocks of kBlock and never assembled into a full array.
template <typename T, class Produce>
void write_streamed(h5::FileWriter& file, const std::string& path, hid_t file_type,
                    hid_t mem_type, std::size_t length, Produce&& produce,
                    std::size_t minor = 0) {
    const h5::Handle dataset = file.create_dataset(path, file_type, length, length,
                                                   h5::Filter::PyTablesBlosc, 0, minor);
    std::vector<T> buffer;
    buffer.reserve(std::min(length, kBlock));
    std::size_t offset = 0;
    const auto flush = [&]() {
        h5::FileWriter::write_block(dataset.get(), mem_type, offset, buffer.size(),
                                    buffer.data());
        offset += buffer.size();
        buffer.clear();
    };
    produce([&](T value) {
        buffer.push_back(value);
        if (buffer.size() == kBlock) {
            flush();
        }
    });
    flush();
    if (offset != length) {
        throw h5::Error(path + ": produced " + std::to_string(offset) +
                        " values but the dataset holds " + std::to_string(length));
    }
}

template <typename T>
void write_vector(h5::FileWriter& file, const std::string& path, hid_t file_type,
                  hid_t mem_type, const std::vector<T>& values) {
    write_streamed<T>(file, path, file_type, mem_type, values.size(),
                      [&](auto emit) {
                          for (const T value : values) {
                              emit(value);
                          }
                      });
}

// The /matrix/data column, written in the dtype the matrix carries. The values
// live as double in memory (see sparse_matrix.hpp), so an integer matrix is
// converted back on the way out; every value came from an integer dataset in
// the first place, so the conversion is exact.
void write_matrix_data(h5::FileWriter& file, const CsrMatrix& matrix, bool upper_only,
                       std::size_t length) {
    const auto produce = [&](auto emit) {
        const auto visit = [&](std::int64_t, std::int64_t, double value) { emit(value); };
        if (upper_only) {
            matrix.for_each_upper(visit);
        } else {
            matrix.for_each_stored(visit);
        }
    };
    const std::string& dtype = matrix.dtype();
    if (dtype == "float32") {
        write_streamed<float>(file, "/matrix/data", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT,
                              length, [&](auto emit) {
                                  produce([&](double value) {
                                      emit(static_cast<float>(value));
                                  });
                              });
    } else if (dtype == "int64") {
        write_streamed<std::int64_t>(file, "/matrix/data", H5T_STD_I64LE,
                                     H5T_NATIVE_INT64, length, [&](auto emit) {
                                         produce([&](double value) {
                                             emit(static_cast<std::int64_t>(value));
                                         });
                                     });
    } else if (dtype == "int32") {
        write_streamed<std::int32_t>(file, "/matrix/data", H5T_STD_I32LE,
                                     H5T_NATIVE_INT32, length, [&](auto emit) {
                                         produce([&](double value) {
                                             emit(static_cast<std::int32_t>(value));
                                         });
                                     });
    } else {
        write_streamed<double>(file, "/matrix/data", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
                               length, [&](auto emit) {
                                   produce([&](double value) { emit(value); });
                               });
    }
}

// A fixed width byte string column, the layout numpy gives an array of Python
// strings and therefore the layout PyTables writes for /intervals/chr_list.
void write_string_column(h5::FileWriter& file, const std::string& path,
                         const std::vector<CutInterval>& intervals,
                         bool use_extra_text) {
    std::size_t width = 0;
    for (const CutInterval& interval : intervals) {
        width = std::max(width,
                         (use_extra_text ? interval.extra_text : interval.chrom).size());
    }
    const h5::Handle type = h5::fixed_string_type(width);
    const h5::Handle dataset = file.create_dataset(path, type.get(), intervals.size(),
                                                   intervals.size(),
                                                   h5::Filter::PyTablesBlosc);
    const std::size_t item = std::max<std::size_t>(width, 1);
    std::vector<char> buffer;
    buffer.reserve(std::min(intervals.size(), kBlock) * item);
    std::size_t offset = 0;
    const auto flush = [&]() {
        h5::FileWriter::write_block(dataset.get(), type.get(), offset,
                                    buffer.size() / item, buffer.data());
        offset += buffer.size() / item;
        buffer.clear();
    };
    for (const CutInterval& interval : intervals) {
        const std::string& text = use_extra_text ? interval.extra_text : interval.chrom;
        const std::size_t previous = buffer.size();
        buffer.resize(previous + item, '\0');
        std::copy_n(text.data(), std::min(text.size(), item), buffer.begin() + static_cast<std::ptrdiff_t>(previous));
        if (buffer.size() == kBlock * item) {
            flush();
        }
    }
    flush();
}

// PyTables marks every node it creates with CLASS, VERSION and TITLE, and the
// root additionally with PYTABLES_FORMAT_VERSION and the file title. The
// markers are not needed to read the file back (PyTables infers a chunked
// dataset as a CArray without them), but writing them keeps a C++ written file
// indistinguishable from a Python written one at the node level.
void mark_pytables_nodes(h5::FileWriter& file, const std::string& title,
                         const std::vector<std::string>& groups,
                         const std::vector<std::string>& leaves) {
    file.set_bytes_attribute("/", "CLASS", "GROUP");
    file.set_bytes_attribute("/", "PYTABLES_FORMAT_VERSION", "2.1");
    file.set_bytes_attribute("/", "TITLE", title);
    file.set_bytes_attribute("/", "VERSION", "1.0");
    for (const std::string& group : groups) {
        file.set_bytes_attribute(group, "CLASS", "GROUP");
        file.set_bytes_attribute(group, "TITLE", "");
        file.set_bytes_attribute(group, "VERSION", "1.0");
    }
    for (const std::string& leaf : leaves) {
        file.set_bytes_attribute(leaf, "CLASS", "CARRAY");
        file.set_bytes_attribute(leaf, "TITLE", "", true);
        file.set_bytes_attribute(leaf, "VERSION", "1.1");
    }
}

// Defined below write_hicexplorer_h5, next to the rest of the metadata
// writing, and declared here because both writers call it.
void write_h5_metadata(h5::FileWriter& file, const MatrixData& data, std::int64_t rows,
                       std::int64_t cols);

}  // namespace

void write_hicexplorer_h5(const std::string& path, const MatrixData& data,
                          const H5SaveOptions& options) {
    // hicmatrix appends the suffix and unlinks an existing file
    // (hicmatrix/lib/h5.py:100-109). H5Fcreate with H5F_ACC_TRUNC is the same
    // thing for a regular file and also handles the case where the path exists
    // but is not an HDF5 file.
    std::string filename = path;
    if (filename.size() < 3 || filename.compare(filename.size() - 3, 3, ".h5") != 0) {
        filename += ".h5";
    }
    std::remove(filename.c_str());

    const CsrMatrix& matrix = data.matrix;
    if (data.cut_intervals.size() != static_cast<std::size_t>(matrix.rows())) {
        throw h5::Error("the bin table has " + std::to_string(data.cut_intervals.size()) +
                        " entries but the matrix has " + std::to_string(matrix.rows()) +
                        " rows");
    }
    const bool upper_only = options.symmetric;
    const std::vector<std::int64_t> indptr = upper_only
                                                 ? matrix.upper_triangle_indptr()
                                                 : matrix.stored_indptr_without_zeros();
    const std::size_t nnz = static_cast<std::size_t>(indptr.back());

    h5::FileWriter file(filename);
    file.create_group("/matrix");
    file.create_group("/intervals");

    write_matrix_data(file, matrix, upper_only, nnz);

    // scipy holds the column indices and the row offsets in the same index
    // type, so both follow the same choice.
    if (index_fits_in_int32(nnz, matrix.rows())) {
        write_streamed<std::int32_t>(
            file, "/matrix/indices", H5T_STD_I32LE, H5T_NATIVE_INT32, nnz,
            [&](auto emit) {
                const auto visit = [&](std::int64_t, std::int64_t column, double) {
                    emit(static_cast<std::int32_t>(column));
                };
                if (upper_only) {
                    matrix.for_each_upper(visit);
                } else {
                    matrix.for_each_stored(visit);
                }
            });
        write_streamed<std::int32_t>(
            file, "/matrix/indptr", H5T_STD_I32LE, H5T_NATIVE_INT32, indptr.size(),
            [&](auto emit) {
                for (const std::int64_t offset : indptr) {
                    emit(static_cast<std::int32_t>(offset));
                }
            });
    } else {
        write_streamed<std::int64_t>(
            file, "/matrix/indices", H5T_STD_I64LE, H5T_NATIVE_INT64, nnz,
            [&](auto emit) {
                const auto visit = [&](std::int64_t, std::int64_t column, double) {
                    emit(column);
                };
                if (upper_only) {
                    matrix.for_each_upper(visit);
                } else {
                    matrix.for_each_stored(visit);
                }
            });
        write_vector(file, "/matrix/indptr", H5T_STD_I64LE, H5T_NATIVE_INT64, indptr);
    }

    write_h5_metadata(file, data, matrix.rows(), matrix.cols());
}

namespace {

// Everything of the file except /matrix/{data,indices,indptr}: the shape, the
// bin table, the three optional nodes and the PyTables node markers. Shared by
// the CsrMatrix writer above and the DenseRowSource one below, which differ
// only in how the CSR arrays are produced.
void write_h5_metadata(h5::FileWriter& file, const MatrixData& data, std::int64_t rows,
                       std::int64_t cols) {
    const std::vector<std::int64_t> shape{rows, cols};
    write_vector(file, "/matrix/shape", H5T_STD_I64LE, H5T_NATIVE_INT64, shape);

    write_string_column(file, "/intervals/chr_list", data.cut_intervals, false);
    write_streamed<std::int64_t>(file, "/intervals/start_list", H5T_STD_I64LE,
                                 H5T_NATIVE_INT64, data.cut_intervals.size(),
                                 [&](auto emit) {
                                     for (const CutInterval& bin : data.cut_intervals) {
                                         emit(bin.start);
                                     }
                                 });
    write_streamed<std::int64_t>(file, "/intervals/end_list", H5T_STD_I64LE,
                                 H5T_NATIVE_INT64, data.cut_intervals.size(),
                                 [&](auto emit) {
                                     for (const CutInterval& bin : data.cut_intervals) {
                                         emit(bin.end);
                                     }
                                 });
    // extra_list is float64 for every matrix in the corpus except the z score
    // matrices of hicFindTADs, which store text there.
    const bool extra_is_text =
        !data.cut_intervals.empty() &&
        std::all_of(data.cut_intervals.begin(), data.cut_intervals.end(),
                    [](const CutInterval& bin) { return !bin.extra_text.empty(); });
    if (extra_is_text) {
        write_string_column(file, "/intervals/extra_list", data.cut_intervals, true);
    } else {
        write_streamed<double>(file, "/intervals/extra_list", H5T_IEEE_F64LE,
                               H5T_NATIVE_DOUBLE, data.cut_intervals.size(),
                               [&](auto emit) {
                                   for (const CutInterval& bin : data.cut_intervals) {
                                       emit(bin.extra);
                                   }
                               });
    }

    // The three optional nodes are omitted when they are empty, which is what
    // the `if len(...)` guards in the Python do.
    if (!data.nan_bins.empty()) {
        write_vector(file, "/nan_bins", H5T_STD_I64LE, H5T_NATIVE_INT64, data.nan_bins);
    }
    if (data.correction_factors.has_value() && !data.correction_factors->empty()) {
        // h5.py:157-158 replaces NaN with zero before writing; Inf is left as
        // it is, unlike on the read side.
        write_streamed<double>(file, "/correction_factors", H5T_IEEE_F64LE,
                               H5T_NATIVE_DOUBLE, data.correction_factors->size(),
                               [&](auto emit) {
                                   for (const double value : *data.correction_factors) {
                                       emit(std::isnan(value) ? 0.0 : value);
                                   }
                               },
                               data.correction_factors_are_column ? 1 : 0);
    }
    if (data.distance_counts.has_value() && !data.distance_counts->empty()) {
        write_vector(file, "/distance_counts", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
                     *data.distance_counts);
    }

    std::vector<std::string> leaves{"/matrix/data",         "/matrix/indices",
                                    "/matrix/indptr",       "/matrix/shape",
                                    "/intervals/chr_list",  "/intervals/start_list",
                                    "/intervals/end_list",  "/intervals/extra_list"};
    if (!data.nan_bins.empty()) {
        leaves.emplace_back("/nan_bins");
    }
    if (data.correction_factors.has_value() && !data.correction_factors->empty()) {
        leaves.emplace_back("/correction_factors");
    }
    if (data.distance_counts.has_value() && !data.distance_counts->empty()) {
        leaves.emplace_back("/distance_counts");
    }
    mark_pytables_nodes(file, "HiCExplorer matrix", {"/matrix", "/intervals"}, leaves);
}

}  // namespace

namespace {

// How many rows the streaming writer computes before it emits them. About
// 32 MB of staging: a rounding error against every budget in cpp/PLAN.md 4.5,
// and enough rows that sixteen workers each get a contiguous range.
std::int64_t staging_row_count(std::int64_t n) {
    if (n <= 0) {
        return 1;
    }
    const std::int64_t rows = 32000000 / (n * 8);
    return std::min<std::int64_t>(256, std::max<std::int64_t>(1, rows));
}

// Fills `count` rows from `first` into `staging`, row r at
// staging + (r - first) * n. Every row is produced by exactly one worker and
// the workers take fixed contiguous ranges, so the buffer is the same for any
// thread count (cpp/OPTIMIZATION.md section 3).
void fill_rows(const DenseRowSource& source, std::int64_t first, std::int64_t count,
               std::int64_t n, double* staging, int threads) {
    const auto work = [&](std::int64_t begin, std::int64_t end) {
        for (std::int64_t row = begin; row < end; ++row) {
            source.fill_row(row, staging + (row - first) * n);
        }
    };
    const int workers = std::max(1, threads);
    if (workers == 1 || count < 2) {
        work(first, first + count);
        return;
    }
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(workers));
    const std::int64_t chunk = (count + workers - 1) / workers;
    for (int w = 0; w < workers; ++w) {
        const std::int64_t begin = first + static_cast<std::int64_t>(w) * chunk;
        const std::int64_t end = std::min(first + count, begin + chunk);
        if (begin >= end) {
            break;
        }
        pool.emplace_back(work, begin, end);
    }
    for (std::thread& worker : pool) {
        worker.join();
    }
}

// Walks the source in chunks, handing each produced row to `consume` in
// increasing row order.
template <class F>
void for_each_produced_row(const DenseRowSource& source, std::int64_t n, int threads,
                           std::vector<double>& staging, F&& consume) {
    const std::int64_t block_rows = static_cast<std::int64_t>(staging.size() /
                                                              static_cast<std::size_t>(n));
    for (std::int64_t first = 0; first < n; first += block_rows) {
        const std::int64_t count = std::min(block_rows, n - first);
        fill_rows(source, first, count, n, staging.data(), threads);
        for (std::int64_t r = 0; r < count; ++r) {
            consume(first + r, staging.data() + r * n);
        }
    }
}

}  // namespace

void write_hicexplorer_h5(const std::string& path, const MatrixData& metadata,
                          const DenseRowSource& source, int threads,
                          const H5SaveOptions& options) {
    std::string filename = path;
    if (filename.size() < 3 || filename.compare(filename.size() - 3, 3, ".h5") != 0) {
        filename += ".h5";
    }
    std::remove(filename.c_str());

    const std::int64_t n = source.rows();
    if (metadata.cut_intervals.size() != static_cast<std::size_t>(n)) {
        throw h5::Error("the bin table has " +
                        std::to_string(metadata.cut_intervals.size()) +
                        " entries but the result has " + std::to_string(n) + " rows");
    }
    if (source.dtype() != "float64") {
        // The only callers are the dense transforms, which are float64. A
        // typed staging buffer per dtype would be dead code with no case
        // behind it, and cpp/AGENTS_CONTRACT.md rule 4 asks for the gap to be
        // named rather than filled with something untested.
        throw h5::Error("the streaming h5 writer only handles float64, not " +
                        source.dtype());
    }
    const bool upper_only = options.symmetric;

    std::vector<double> staging(static_cast<std::size_t>(staging_row_count(n)) *
                                static_cast<std::size_t>(n));

    // Pass one: the row lengths. A PyTables CArray has a fixed length that has
    // to be known before the dataset is created, so the entries have to be
    // counted before any of them can be written.
    std::vector<std::int64_t> indptr(static_cast<std::size_t>(n) + 1, 0);
    for_each_produced_row(source, n, threads, staging,
                          [&](std::int64_t row, const double* values) {
                              std::int64_t stored = 0;
                              for (std::int64_t j = upper_only ? row : 0; j < n; ++j) {
                                  // eliminate_zeros: an exact zero is not
                                  // stored, a NaN is.
                                  if (values[j] != 0.0) {
                                      ++stored;
                                  }
                              }
                              indptr[static_cast<std::size_t>(row) + 1] = stored;
                          });
    for (std::int64_t row = 0; row < n; ++row) {
        indptr[static_cast<std::size_t>(row) + 1] +=
            indptr[static_cast<std::size_t>(row)];
    }
    const std::size_t nnz = static_cast<std::size_t>(indptr.back());

    h5::FileWriter file(filename);
    file.create_group("/matrix");
    file.create_group("/intervals");

    const bool narrow = index_fits_in_int32(nnz, n);
    const h5::Handle data_set =
        file.create_dataset("/matrix/data", H5T_IEEE_F64LE, nnz, nnz,
                            h5::Filter::PyTablesBlosc, 0, 0);
    const h5::Handle index_set = file.create_dataset(
        "/matrix/indices", narrow ? H5T_STD_I32LE : H5T_STD_I64LE, nnz, nnz,
        h5::Filter::PyTablesBlosc, 0, 0);

    // Pass two: the values and the column indices, written into both datasets
    // in the same block so that the source is walked once more, not twice.
    std::vector<double> value_buffer;
    std::vector<std::int32_t> index_buffer32;
    std::vector<std::int64_t> index_buffer64;
    value_buffer.reserve(kBlock);
    if (narrow) {
        index_buffer32.reserve(kBlock);
    } else {
        index_buffer64.reserve(kBlock);
    }
    std::size_t written = 0;
    const auto flush = [&]() {
        const std::size_t count = value_buffer.size();
        if (count == 0) {
            return;
        }
        h5::FileWriter::write_block(data_set.get(), H5T_NATIVE_DOUBLE, written, count,
                                    value_buffer.data());
        if (narrow) {
            h5::FileWriter::write_block(index_set.get(), H5T_NATIVE_INT32, written, count,
                                        index_buffer32.data());
            index_buffer32.clear();
        } else {
            h5::FileWriter::write_block(index_set.get(), H5T_NATIVE_INT64, written, count,
                                        index_buffer64.data());
            index_buffer64.clear();
        }
        written += count;
        value_buffer.clear();
    };
    for_each_produced_row(
        source, n, threads, staging, [&](std::int64_t row, const double* values) {
            for (std::int64_t j = upper_only ? row : 0; j < n; ++j) {
                if (values[j] == 0.0) {
                    continue;
                }
                value_buffer.push_back(values[j]);
                if (narrow) {
                    index_buffer32.push_back(static_cast<std::int32_t>(j));
                } else {
                    index_buffer64.push_back(j);
                }
                if (value_buffer.size() == kBlock) {
                    flush();
                }
            }
        });
    flush();
    if (written != nnz) {
        throw h5::Error("/matrix/data: produced " + std::to_string(written) +
                        " values but the two passes counted " + std::to_string(nnz));
    }

    if (narrow) {
        write_streamed<std::int32_t>(file, "/matrix/indptr", H5T_STD_I32LE,
                                     H5T_NATIVE_INT32, indptr.size(), [&](auto emit) {
                                         for (const std::int64_t offset : indptr) {
                                             emit(static_cast<std::int32_t>(offset));
                                         }
                                     });
    } else {
        write_vector(file, "/matrix/indptr", H5T_STD_I64LE, H5T_NATIVE_INT64, indptr);
    }

    write_h5_metadata(file, metadata, n, n);
}

CsrMatrix materialize_dense_row_source(const DenseRowSource& source, int threads) {
    const std::int64_t n = source.rows();
    std::vector<std::int64_t> indptr(1, 0);
    indptr.reserve(static_cast<std::size_t>(n) + 1);
    std::vector<std::int32_t> indices;
    std::vector<double> data;
    std::vector<double> staging(static_cast<std::size_t>(staging_row_count(n)) *
                                static_cast<std::size_t>(std::max<std::int64_t>(n, 1)));
    if (n > 0) {
        for_each_produced_row(source, n, threads, staging,
                              [&](std::int64_t, const double* values) {
                                  for (std::int64_t j = 0; j < n; ++j) {
                                      if (values[j] != 0.0) {
                                          indices.push_back(static_cast<std::int32_t>(j));
                                          data.push_back(values[j]);
                                      }
                                  }
                                  indptr.push_back(
                                      static_cast<std::int64_t>(data.size()));
                              });
    }
    return CsrMatrix(n, n, std::move(indptr), std::move(indices), std::move(data),
                     source.dtype());
}

}  // namespace hicx
