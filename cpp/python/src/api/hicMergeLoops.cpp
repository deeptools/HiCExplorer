// hicexplorer._hicx.hicMergeLoops: the C++ bodies of mergeLoops and readFile.
//
// mergeLoops runs hicx::merge_loops (cpp/core/include/hicx/merge_loops.hpp)
// on CPython-exact interval trees; this file performs the Python's dictionary
// lookups with Python's key equality, the int() conversions of the tuple and
// loop fields, and raises the 3.7.6 errors (UnboundLocalError for a first
// loop on an unknown chromosome, KeyError for a y tree looked up under
// loop[0], ValueError for an interval removed from a tree that does not hold
// it). readFile reads the file with hicx::TextTable, pandas.read_csv's
// dtype inference and float parser, and builds the DataFrame.
//
// Deviations: loop coordinates that are neither integers nor integral floats
// raise TypeError naming HiCExplorer 4 (3.7.6 computed with them only for a
// loop whose chromosome is in the trees); readFile infers int64, float64 and
// object columns, not bool, and raises pandas.errors.ParserError with the
// tool's message for a ragged file.

#include <cmath>
#include <exception>
#include <optional>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicx/merge_loops.hpp"
#include "hicx/text_table.hpp"
#include "loops_correct_common.hpp"

namespace hicx::python {

namespace ml = hicx::merge_loops;

namespace {

// int(object)
std::int64_t python_int(py::handle object) {
    PyObject* value = PyNumber_Long(object.ptr());
    if (value == nullptr) {
        throw_python_error();
    }
    const py::object number = py::reinterpret_steal<py::object>(value);
    return to_int64(number);
}

// A coordinate the overlap query computes with.
std::int64_t coordinate(py::handle object) {
    if (is_integer(object) && !PyBool_Check(object.ptr())) {
        return to_int64(object);
    }
    if (PyFloat_Check(object.ptr()) || py::isinstance(object, py::module_::import("numpy").attr("floating"))) {
        const double value = to_double(object);
        if (std::isfinite(value) && value == std::trunc(value)) {
            return static_cast<std::int64_t>(value);
        }
    }
    raise_error("TypeError", "mergeLoops: HiCExplorer 4 takes integer loop coordinates, not " + repr_of(object));
}

// hiCMatrix.intervalListToIntervalTree, the trees as the merge kernel takes
// them and the dictionary from chromosome to tree.
ml::Intervals interval_trees(py::handle interval_list, py::dict& dictionary) {
    ml::Intervals out;
    py::object previous = py::none();
    for (const py::handle interval : interval_list) {
        const py::object head = lc::item(interval, py::slice(0, 3, 1));
        if (lc::length(head) != 3) {
            const std::size_t got = lc::length(head);
            raise_error("ValueError", got < 3 ? "not enough values to unpack (expected 3, got " +
                                                    std::to_string(got) + ")"
                                              : "too many values to unpack (expected 3)");
        }
        const py::object chrom = lc::item(head, 0);
        const std::int64_t start = python_int(lc::item(head, 1));
        const std::int64_t end = python_int(lc::item(head, 2));
        if (lc::compare(previous, chrom, Py_NE)) {
            dictionary[chrom] = py::int_(out.tree_count);
            ++out.tree_count;
            previous = chrom;
        }
        out.tree.push_back(dictionary[chrom].cast<std::int64_t>());
        out.begin.push_back(start);
        out.end.push_back(end);
    }
    for (const auto& entry : dictionary) {
        out.dictionary.push_back(py::reinterpret_borrow<py::object>(entry.second).cast<std::int64_t>());
    }
    return out;
}

std::int64_t lookup(const py::dict& dictionary, py::handle key) {
    PyObject* value = PyDict_GetItemWithError(dictionary.ptr(), key.ptr());
    if (value == nullptr) {
        if (PyErr_Occurred()) {
            throw_python_error();
        }
        return -1;
    }
    return py::reinterpret_borrow<py::object>(value).cast<std::int64_t>();
}

bool contains(const py::dict& dictionary, py::handle key) {
    const int result = PyDict_Contains(dictionary.ptr(), key.ptr());
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

[[noreturn]] void raise_key_error(py::handle key) {
    PyErr_SetObject(PyExc_KeyError, key.ptr());
    throw py::error_already_set();
}

py::object merge_loops(py::handle pDataFrame, py::handle pLowestResolution, py::handle pTupleX,
                       py::handle pTupleY) {
    py::dict dictionary_x;
    py::dict dictionary_y;
    const ml::Intervals x = interval_trees(pTupleX, dictionary_x);
    const ml::Intervals y = interval_trees(pTupleY, dictionary_y);
    const py::object values = pDataFrame.attr("values");
    std::vector<ml::Loop> loops;
    std::vector<py::object> chromosomes;
    std::exception_ptr pending;
    try {
        for (const py::handle loop : values) {
            const std::int64_t lowest = python_int(pLowestResolution);
            const std::int64_t factor_x =
                lowest - std::llabs(python_int(lc::item(loop, 2)) - python_int(lc::item(loop, 1)));
            const std::int64_t factor_y =
                python_int(pLowestResolution) - std::llabs(python_int(lc::item(loop, 5)) - python_int(lc::item(loop, 4)));
            const py::object chrom = lc::item(loop, 0);
            ml::Loop row;
            row.x_tree = lookup(dictionary_x, chrom);
            row.y_guard = contains(dictionary_y, lc::item(loop, 3));
            row.y_tree = lookup(dictionary_y, chrom);
            row.x_begin = coordinate(lc::item(loop, 1)) - factor_x - 1;
            row.x_end = coordinate(lc::item(loop, 2)) + factor_x + 1;
            row.y_begin = coordinate(lc::item(loop, 4)) - factor_y - 1;
            row.y_end = coordinate(lc::item(loop, 5)) + factor_y + 1;
            loops.push_back(row);
            chromosomes.push_back(chrom);
        }
    } catch (...) {
        // raised when the Python reaches this loop, unless an earlier one fails
        pending = std::current_exception();
    }
    const ml::Result result = ml::merge(x, y, loops);
    switch (result.error) {
        case ml::Error::None:
            break;
        case ml::Error::UnboundX:
            raise_error("UnboundLocalError",
                        "cannot access local variable 'x_interval' where it is not associated with a value");
        case ml::Error::UnboundY:
            raise_error("UnboundLocalError",
                        "cannot access local variable 'y_interval' where it is not associated with a value");
        case ml::Error::MissingYTree:
        case ml::Error::MissingXRemove:
        case ml::Error::MissingYRemove:
            raise_key_error(chromosomes[result.loop]);
        case ml::Error::NotInXTree:
        case ml::Error::NotInYTree:
            raise_error("ValueError", "");
    }
    if (pending) {
        std::rethrow_exception(pending);
    }
    py::list ids;
    for (const std::int64_t id : result.ids) {
        ids.append(py::int_(id));
    }
    return ids;
}

py::object read_loop_file(py::handle pFile) {
    // open()'s errors, which are pandas' too
    (void)hicx::python::read_file(pFile);
    const std::string path = fs_path(pFile);
    hicx::TextTable table;
    try {
        table = hicx::TextTable::read_tsv(path);
    } catch (const std::exception& error) {
        raise_error("pandas.errors.ParserError", error.what());
    }
    if (table.cols() == 0) {
        raise_error("pandas.errors.EmptyDataError", "No columns to parse from file");
    }
    const py::module_ np = py::module_::import("numpy");
    py::dict columns;
    for (std::size_t c = 0; c < table.cols(); ++c) {
        const hicx::TableColumn& column = table.column(c);
        py::object array;
        switch (column.type) {
            case hicx::ColumnType::Int64:
                array = from_vector(column.ints);
                break;
            case hicx::ColumnType::Float64:
                array = from_vector(column.floats);
                break;
            default: {
                py::list cells;
                for (const std::string& text : column.strings) {
                    if (hicx::is_pandas_na(text)) {
                        cells.append(py::float_(std::nan("")));
                    } else {
                        cells.append(to_str(text));
                    }
                }
                array = np.attr("array")(cells, py::arg("dtype") = py::module_::import("builtins").attr("object"));
            }
        }
        columns[py::int_(c)] = array;
    }
    return py::module_::import("pandas").attr("DataFrame")(columns);
}

}  // namespace

void bind_hicMergeLoops(py::module_& m) {
    m.def("mergeLoops", &merge_loops, py::arg("pDataFrame"), py::arg("pLowestResolution"), py::arg("pTupleX"),
          py::arg("pTupleY"));
    m.def("readFile", &read_loop_file, py::arg("pFile"));
}

}  // namespace hicx::python
