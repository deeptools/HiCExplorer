#include "hicx/correlate_api.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

#include "hicx/clustering.hpp"
#include "hicx/python_error.hpp"

namespace hicx::correlate {

namespace {

// ---------------------------------------------------------------------------
// the arithmetic of scipy's sparsetools in the type T: integers wrap

template <class T>
T add(T a, T b) {
    if constexpr (std::is_integral_v<T>) {
        using U = std::make_unsigned_t<T>;
        return static_cast<T>(static_cast<U>(static_cast<U>(a) + static_cast<U>(b)));
    } else {
        return a + b;
    }
}

template <class T>
T subtract(T a, T b) {
    if constexpr (std::is_integral_v<T>) {
        using U = std::make_unsigned_t<T>;
        return static_cast<T>(static_cast<U>(static_cast<U>(a) - static_cast<U>(b)));
    } else {
        return a - b;
    }
}

// A compressed matrix: indptr over the major axis (rows for csr, columns for
// csc), `minor` the indices, `data` the values in storage order.
template <class T>
struct Compressed {
    std::int64_t major_size = 0;
    std::int64_t minor_size = 0;
    std::vector<std::int64_t> indptr;
    std::vector<std::int64_t> minor;
    std::vector<T> data;
};

// csr_has_canonical_format: sorted, duplicate free indices in every row.
template <class T>
bool has_canonical_format(const Compressed<T>& m) {
    for (std::int64_t i = 0; i < m.major_size; ++i) {
        const auto begin = static_cast<std::size_t>(m.indptr[static_cast<std::size_t>(i)]);
        const auto end = static_cast<std::size_t>(m.indptr[static_cast<std::size_t>(i) + 1]);
        if (begin > end) {
            return false;
        }
        for (std::size_t k = begin + 1; k < end; ++k) {
            if (!(m.minor[k - 1] < m.minor[k])) {
                return false;
            }
        }
    }
    return true;
}

enum class Operation { Plus, Minus };

template <class T>
T apply(Operation operation, T a, T b) {
    return operation == Operation::Plus ? add(a, b) : subtract(a, b);
}

// csr_plus_csr / csr_minus_csr (and the csc versions, which swap the axes):
// the canonical merge when both operands have canonical indices, the general
// scan otherwise; a result equal to 0 is not stored.
template <class T>
Compressed<T> binop(const Compressed<T>& a, const Compressed<T>& b, Operation operation) {
    Compressed<T> c;
    c.major_size = a.major_size;
    c.minor_size = a.minor_size;
    c.indptr.assign(static_cast<std::size_t>(a.major_size) + 1, 0);
    const T zero = T{};
    const auto at = [](std::int64_t value) { return static_cast<std::size_t>(value); };
    if (has_canonical_format(a) && has_canonical_format(b)) {
        for (std::int64_t i = 0; i < a.major_size; ++i) {
            std::size_t pa = at(a.indptr[at(i)]);
            const std::size_t ea = at(a.indptr[at(i) + 1]);
            std::size_t pb = at(b.indptr[at(i)]);
            const std::size_t eb = at(b.indptr[at(i) + 1]);
            const auto emit = [&](std::int64_t index, T value) {
                if (value != zero) {
                    c.minor.push_back(index);
                    c.data.push_back(value);
                }
            };
            while (pa < ea && pb < eb) {
                const std::int64_t ja = a.minor[pa];
                const std::int64_t jb = b.minor[pb];
                if (ja == jb) {
                    emit(ja, apply(operation, a.data[pa], b.data[pb]));
                    ++pa;
                    ++pb;
                } else if (ja < jb) {
                    emit(ja, apply(operation, a.data[pa], zero));
                    ++pa;
                } else {
                    emit(jb, apply(operation, zero, b.data[pb]));
                    ++pb;
                }
            }
            for (; pa < ea; ++pa) {
                emit(a.minor[pa], apply(operation, a.data[pa], zero));
            }
            for (; pb < eb; ++pb) {
                emit(b.minor[pb], apply(operation, zero, b.data[pb]));
            }
            c.indptr[at(i) + 1] = static_cast<std::int64_t>(c.data.size());
        }
    } else {
        std::vector<std::int64_t> next(at(a.minor_size), -1);
        std::vector<T> a_row(at(a.minor_size), zero);
        std::vector<T> b_row(at(a.minor_size), zero);
        for (std::int64_t i = 0; i < a.major_size; ++i) {
            std::int64_t head = -2;
            std::int64_t length = 0;
            for (std::size_t k = at(a.indptr[at(i)]); k < at(a.indptr[at(i) + 1]); ++k) {
                const std::int64_t j = a.minor[k];
                a_row[at(j)] = add(a_row[at(j)], a.data[k]);
                if (next[at(j)] == -1) {
                    next[at(j)] = head;
                    head = j;
                    ++length;
                }
            }
            for (std::size_t k = at(b.indptr[at(i)]); k < at(b.indptr[at(i) + 1]); ++k) {
                const std::int64_t j = b.minor[k];
                b_row[at(j)] = add(b_row[at(j)], b.data[k]);
                if (next[at(j)] == -1) {
                    next[at(j)] = head;
                    head = j;
                    ++length;
                }
            }
            for (std::int64_t k = 0; k < length; ++k) {
                const T result = apply(operation, a_row[at(head)], b_row[at(head)]);
                if (result != zero) {
                    c.minor.push_back(head);
                    c.data.push_back(result);
                }
                const std::int64_t temp = head;
                head = next[at(head)];
                next[at(temp)] = -1;
                a_row[at(temp)] = zero;
                b_row[at(temp)] = zero;
            }
            c.indptr[at(i) + 1] = static_cast<std::int64_t>(c.data.size());
        }
    }
    return c;
}

// csr_tocsc: the entries regrouped by their minor index, in the order they
// are stored (so each group is sorted by the old major index).
template <class T>
Compressed<T> transpose_storage(const Compressed<T>& m) {
    Compressed<T> t;
    t.major_size = m.minor_size;
    t.minor_size = m.major_size;
    const auto at = [](std::int64_t value) { return static_cast<std::size_t>(value); };
    t.indptr.assign(at(t.major_size) + 1, 0);
    const std::size_t stored = at(m.indptr[at(m.major_size)]);
    for (std::size_t k = 0; k < stored; ++k) {
        ++t.indptr[at(m.minor[k]) + 1];
    }
    for (std::size_t j = 0; j < at(t.major_size); ++j) {
        t.indptr[j + 1] += t.indptr[j];
    }
    t.minor.resize(stored);
    t.data.resize(stored);
    std::vector<std::int64_t> fill(t.indptr.begin(), t.indptr.end() - 1);
    for (std::int64_t i = 0; i < m.major_size; ++i) {
        for (std::size_t k = at(m.indptr[at(i)]); k < at(m.indptr[at(i) + 1]); ++k) {
            const std::size_t destination = at(fill[at(m.minor[k])]++);
            t.minor[destination] = i;
            t.data[destination] = m.data[k];
        }
    }
    return t;
}

// coo_tocsr (entries in input order within a row) followed by sum_duplicates:
// sorted indices, duplicates added.
template <class T>
Compressed<T> coo_to_compressed(const SparseInput<T>& in, bool to_csr) {
    Compressed<T> m;
    m.major_size = to_csr ? in.rows : in.cols;
    m.minor_size = to_csr ? in.cols : in.rows;
    const auto at = [](std::int64_t value) { return static_cast<std::size_t>(value); };
    const std::size_t stored = in.data.size();
    const std::vector<std::int64_t>& major = to_csr ? in.major : in.minor;
    const std::vector<std::int64_t>& minor = to_csr ? in.minor : in.major;
    std::vector<std::int64_t> pointer(at(m.major_size) + 1, 0);
    for (std::size_t k = 0; k < stored; ++k) {
        ++pointer[at(major[k]) + 1];
    }
    for (std::size_t i = 0; i < at(m.major_size); ++i) {
        pointer[i + 1] += pointer[i];
    }
    std::vector<std::int64_t> fill(pointer.begin(), pointer.end() - 1);
    std::vector<std::int64_t> indices(stored);
    std::vector<T> values(stored);
    for (std::size_t k = 0; k < stored; ++k) {
        const std::size_t destination = at(fill[at(major[k])]++);
        indices[destination] = minor[k];
        values[destination] = in.data[k];
    }
    m.indptr.assign(at(m.major_size) + 1, 0);
    std::vector<std::size_t> order;
    for (std::int64_t i = 0; i < m.major_size; ++i) {
        const std::size_t begin = at(pointer[at(i)]);
        const std::size_t end = at(pointer[at(i) + 1]);
        order.resize(end - begin);
        for (std::size_t k = 0; k < order.size(); ++k) {
            order[k] = begin + k;
        }
        std::stable_sort(order.begin(), order.end(),
                         [&](std::size_t x, std::size_t y) { return indices[x] < indices[y]; });
        for (std::size_t k = 0; k < order.size();) {
            const std::int64_t j = indices[order[k]];
            T sum = values[order[k]];
            std::size_t l = k + 1;
            for (; l < order.size() && indices[order[l]] == j; ++l) {
                sum = add(sum, values[order[l]]);
            }
            m.minor.push_back(j);
            m.data.push_back(sum);
            k = l;
        }
        m.indptr[at(i) + 1] = static_cast<std::int64_t>(m.data.size());
    }
    return m;
}

// self.__class__(other): the matrix in csr (to_csr) or csc form.
template <class T>
Compressed<T> as_format(const SparseInput<T>& in, bool to_csr) {
    if (in.format == SparseFormat::Coo) {
        return coo_to_compressed(in, to_csr);
    }
    Compressed<T> m;
    const bool is_csr = in.format == SparseFormat::Csr;
    m.major_size = is_csr ? in.rows : in.cols;
    m.minor_size = is_csr ? in.cols : in.rows;
    m.indptr = in.indptr;
    m.minor = in.minor;
    m.data = in.data;
    return is_csr == to_csr ? m : transpose_storage(m);
}

}  // namespace

template <class T>
std::pair<std::vector<T>, std::vector<T>> get_vectors(const SparseInput<T>& mat1,
                                                      const SparseInput<T>& mat2) {
    const bool csr = mat1.format != SparseFormat::Csc;
    const Compressed<T> first = as_format(mat1, csr);
    const Compressed<T> second = as_format(mat2, csr);
    Compressed<T> sum = binop(first, second, Operation::Plus);
    for (T& value : sum.data) {
        value = add(value, T{1});
    }
    std::pair<std::vector<T>, std::vector<T>> out;
    out.first = binop(sum, first, Operation::Minus).data;
    out.second = binop(sum, second, Operation::Minus).data;
    for (T& value : out.first) {
        value = subtract(value, T{1});
    }
    for (T& value : out.second) {
        value = subtract(value, T{1});
    }
    return out;
}

template std::pair<std::vector<std::int32_t>, std::vector<std::int32_t>> get_vectors(
    const SparseInput<std::int32_t>&, const SparseInput<std::int32_t>&);
template std::pair<std::vector<std::int64_t>, std::vector<std::int64_t>> get_vectors(
    const SparseInput<std::int64_t>&, const SparseInput<std::int64_t>&);
template std::pair<std::vector<float>, std::vector<float>> get_vectors(const SparseInput<float>&,
                                                                       const SparseInput<float>&);
template std::pair<std::vector<double>, std::vector<double>> get_vectors(
    const SparseInput<double>&, const SparseInput<double>&);

double matrix_minimum(const std::vector<double>& values) {
    if (values.empty()) {
        throw PythonError("ValueError",
                          "zero-size array to reduction operation minimum which has no identity");
    }
    double minimum = values[0];
    for (const double value : values) {
        if (std::isnan(value)) {
            return value;
        }
        minimum = std::min(minimum, value);
    }
    return minimum;
}

bool looks_like_distance_matrix(const std::vector<double>& values, std::int64_t rows,
                                std::int64_t cols) {
    if (rows != cols) {
        return false;
    }
    const auto n = static_cast<std::size_t>(rows);
    const auto close = [](double a, double b) {
        if (std::isnan(a) || std::isnan(b)) {
            return false;
        }
        if (std::isinf(a) || std::isinf(b)) {
            return a == b;
        }
        return std::abs(a - b) <= 1e-8 + 1e-5 * std::abs(b);
    };
    for (std::size_t i = 0; i < n; ++i) {
        if (!close(values[i * n + i], 0.0)) {
            return false;
        }
    }
    for (const double value : values) {
        if (!(value >= 0.0)) {
            return false;
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (!close(values[i * n + j], values[j * n + i])) {
                return false;
            }
        }
    }
    return true;
}

std::vector<double> complete_linkage_of_rows(const std::vector<double>& values, std::int64_t rows,
                                             std::int64_t cols) {
    if (rows < 2) {
        throw PythonError("ValueError",
                          "The number of observations cannot be determined on an empty distance "
                          "matrix.");
    }
    const char* const message = "The condensed distance matrix must contain only finite values.";
    for (const double value : values) {
        if (!std::isfinite(value)) {
            throw PythonError("ValueError", message);
        }
    }
    // pdist(rows, 'euclidean') overflows to inf for huge differences, which
    // linkage refuses.
    for (std::int64_t i = 0; i < rows; ++i) {
        for (std::int64_t j = i + 1; j < rows; ++j) {
            double sum = 0.0;
            for (std::int64_t k = 0; k < cols; ++k) {
                const double difference =
                    std::abs(values[static_cast<std::size_t>(i * cols + k)] -
                             values[static_cast<std::size_t>(j * cols + k)]);
                sum = sum + difference * difference;
            }
            if (!std::isfinite(std::sqrt(sum))) {
                throw PythonError("ValueError", message);
            }
        }
    }
    cluster::Samples observations;
    observations.samples = rows;
    observations.features = cols;
    observations.values = values;
    return cluster::detail::complete_linkage(observations);
}

}  // namespace hicx::correlate
