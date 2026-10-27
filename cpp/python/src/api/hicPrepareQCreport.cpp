// hicexplorer._hicx.hicPrepareQCreport: the C++ bodies of the public
// functions of hicexplorer.hicPrepareQCreport.
//
// save_html renders the five DataFrames with hicx::qc_styler_html (the code
// the hicPrepareQCreport executable uses for hicQC.html); the make_figure_*
// functions keep their matplotlib/pandas drawing in Python and take the
// DataFrames they plot, and the ones they return, from the *_tables functions
// here (hicx/qc_report_api.hpp for the arithmetic). make_sure_path_exists is
// os.makedirs.
//
// Differences from 3.7.6 (cpp/python/API_BINDINGS.md):
//  * DataFrames have a single level index and single level columns; labels are
//    str, int or float; cells are Python/numpy ints, floats or str (others
//    raise TypeError naming HiCExplorer 4).
//  * make_figure_* with duplicate row labels or duplicate column labels in
//    the selected columns raise ValueError; pandas would return a cartesian
//    join.
//  * The messages of the TypeErrors for text in a count column follow numpy's
//    operand message for an object array; pandas' own wording for
//    nonstandard combinations of text cells may differ.
//  * save_html takes a path (str, bytes or PathLike), not an open file
//    descriptor.

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>
#include <cstring>
#include <set>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicx/qc_report_api.hpp"
#include "hicx/qc_html.hpp"

namespace hicx::python {

namespace {

constexpr const char* kLogger = "hicexplorer.hicPrepareQCreport";

using hicx::QcStylerCell;
using hicx::qcapi::Column;

// ---------------------------------------------------------------------------
// reading DataFrames

[[noreturn]] void refuse_value(py::handle value) {
    raise_error("TypeError", "hicPrepareQCreport: HiCExplorer 4 reads int, float and str cells of "
                             "the QC tables, got " + type_name(value));
}

QcStylerCell classify(py::handle value) {
    QcStylerCell cell;
    if (PyBool_Check(value.ptr())) {
        refuse_value(value);
    }
    if (PyLong_Check(value.ptr())) {
        int overflow = 0;
        const long long v = PyLong_AsLongLongAndOverflow(value.ptr(), &overflow);
        if (overflow != 0) {
            refuse_value(value);
        }
        cell.is_int = true;
        cell.integer = v;
    } else if (PyFloat_Check(value.ptr())) {
        cell.is_float = true;
        cell.real = PyFloat_AsDouble(value.ptr());
    } else if (is_str(value)) {
        cell.text = utf8(value);
    } else {
        refuse_value(value);
    }
    return cell;
}

// The default formatter of a label: str as it is, ints plainly, floats with
// six decimals.
std::string label_text(py::handle label) {
    if (is_str(label)) {
        return utf8(label);
    }
    if (PyBool_Check(label.ptr())) {
        refuse_value(label);
    }
    if (PyLong_Check(label.ptr())) {
        return str_of(label);
    }
    if (PyFloat_Check(label.ptr())) {
        return utf8(py::str(py::str("{:.6f}").attr("format")(label)));
    }
    refuse_value(label);
}

py::list to_list(py::handle index_or_series) {
    return py::list(index_or_series.attr("tolist")());
}

// A DataFrame with unique column labels addressed by str name.
struct Frame {
    py::object table;
    py::list columns;  // the column labels as Python objects
    std::vector<std::string> names;  // str labels, "" with a flag for other labels
    std::vector<bool> is_name;

    explicit Frame(py::handle t) : table(py::reinterpret_borrow<py::object>(t)) {
        columns = to_list(table.attr("columns"));
        for (py::handle label : columns) {
            is_name.push_back(is_str(label));
            names.push_back(is_str(label) ? utf8(label) : std::string());
        }
    }

    [[nodiscard]] std::ptrdiff_t find(const std::string& name) const {
        std::ptrdiff_t found = -1;
        for (std::size_t i = 0; i < names.size(); ++i) {
            if (is_name[i] && names[i] == name) {
                if (found >= 0) {
                    raise_error("ValueError", "hicPrepareQCreport: HiCExplorer 4 does not take "
                                              "duplicate column labels ('" + name + "')");
                }
                found = static_cast<std::ptrdiff_t>(i);
            }
        }
        return found;
    }

    [[nodiscard]] Column cells(std::size_t position) const {
        const py::object series = table.attr("iloc")[py::make_tuple(
            py::slice(py::none(), py::none(), py::none()), position)];
        Column out;
        for (py::handle value : to_list(series)) {
            out.push_back(classify(value));
        }
        return out;
    }

    [[nodiscard]] py::array values(std::size_t position) const {
        return py::array(table.attr("iloc")[py::make_tuple(
            py::slice(py::none(), py::none(), py::none()), position)].attr("to_numpy")());
    }

    // table[[names]]: the positions of the columns, with pandas' KeyErrors.
    [[nodiscard]] std::vector<std::size_t> select(const std::vector<std::string>& wanted) const {
        std::vector<std::size_t> positions;
        std::vector<std::string> missing;
        for (const std::string& name : wanted) {
            const std::ptrdiff_t at = find(name);
            if (at < 0) {
                if (std::find(missing.begin(), missing.end(), name) == missing.end()) {
                    missing.push_back(name);
                }
            } else {
                positions.push_back(static_cast<std::size_t>(at));
            }
        }
        if (!missing.empty()) {
            py::list names;
            for (const std::string& name : wanted) {
                names.append(to_str(name));
            }
            if (positions.empty()) {
                // the repr of an object Index, which wraps long lists as pandas does
                const py::object index = py::module_::import("pandas").attr("Index")(
                    names, py::arg("dtype") = py::module_::import("builtins").attr("object"));
                raise_error("KeyError", "None of [" + repr_of(index) + "] are in the [columns]");
            }
            py::list lost;
            for (const std::string& name : missing) {
                lost.append(to_str(name));
            }
            raise_error("KeyError", repr_of(lost) + " not in index");
        }
        return positions;
    }

    // table[name]
    [[nodiscard]] std::size_t series(const std::string& name) const {
        const std::ptrdiff_t at = find(name);
        if (at < 0) {
            raise_error("KeyError", name);
        }
        return static_cast<std::size_t>(at);
    }
};

void require_unique_index(py::handle table) {
    if (!table.attr("index").attr("is_unique").cast<bool>()) {
        raise_error("ValueError",
                    "hicPrepareQCreport: HiCExplorer 4 does not take duplicate row labels");
    }
}

// ---------------------------------------------------------------------------
// building DataFrames

py::array matrix_of(const std::vector<double>& values, std::size_t rows, std::size_t cols) {
    py::array out = from_vector(values);
    return out.reshape(std::vector<py::ssize_t>{static_cast<py::ssize_t>(rows),
                                                static_cast<py::ssize_t>(cols)});
}

py::object names_index(const std::vector<std::string>& names) {
    py::list labels;
    for (const std::string& name : names) {
        labels.append(to_str(name));
    }
    return py::module_::import("pandas").attr("Index")(labels);
}

py::object make_frame(py::handle data, py::handle index, py::handle columns) {
    return py::module_::import("pandas").attr("DataFrame")(data, py::arg("index") = index,
                                                           py::arg("columns") = columns);
}

// table[names] / scalar
py::object divided_by_scalar(const Frame& frame, const std::vector<std::string>& names,
                             double scalar) {
    const std::vector<std::size_t> positions = frame.select(names);
    const py::object index = frame.table.attr("index");
    const std::size_t samples = static_cast<std::size_t>(py::len(index));
    std::vector<std::vector<double>> divided;
    for (const std::size_t position : positions) {
        divided.push_back(hicx::qcapi::divide_by_scalar(frame.cells(position), scalar));
    }
    std::vector<double> values;
    values.reserve(samples * positions.size());
    for (std::size_t s = 0; s < samples; ++s) {
        for (const auto& column : divided) {
            values.push_back(column[s]);
        }
    }
    return make_frame(matrix_of(values, samples, positions.size()), index, names_index(names));
}

struct OutputColumn {
    std::string label;
    std::size_t name;  // index into the selected names
    bool count;        // the table's own column, else its fraction
};

// table[names].T / denominators and table[names].join(fractions.T, rsuffix)
// restricted to `output`, as (fractions, joined).
py::tuple ratio_tables(const Frame& frame, const std::vector<std::string>& names,
                       const std::vector<std::size_t>& positions, const Column& denominators,
                       const std::vector<OutputColumn>& output) {
    if (std::set<std::string>(names.begin(), names.end()).size() != names.size()) {
        raise_error("ValueError", "hicPrepareQCreport: HiCExplorer 4 does not take duplicate "
                                  "column labels");
    }
    const py::object index = frame.table.attr("index");
    const std::size_t samples = static_cast<std::size_t>(py::len(index));
    std::vector<std::vector<double>> fractions;
    for (const std::size_t position : positions) {
        fractions.push_back(hicx::qcapi::divide_columns(frame.cells(position), denominators));
    }
    std::vector<double> flat;
    for (const auto& row : fractions) {
        flat.insert(flat.end(), row.begin(), row.end());
    }
    py::object fraction_table =
        make_frame(matrix_of(flat, names.size(), samples), names_index(names), index);
    py::dict data;
    for (const OutputColumn& column : output) {
        data[to_str(column.label)] = column.count ? frame.values(positions[column.name])
                                                  : from_vector(fractions[column.name]);
    }
    py::object joined = py::module_::import("pandas").attr("DataFrame")(data, py::arg("index") = index);
    return py::make_tuple(fraction_table, joined);
}

std::vector<OutputColumn> paired_output(const std::vector<std::string>& names,
                                        const std::string& suffix) {
    std::vector<OutputColumn> output;
    for (std::size_t i = 0; i < names.size(); ++i) {
        output.push_back({names[i], i, true});
        output.push_back({names[i] + suffix, i, false});
    }
    return output;
}

// The four functions that divide the counts of several columns by one column.
py::tuple by_column(py::handle table, const std::vector<std::string>& names,
                    const std::string& denominator, const std::string& suffix) {
    const Frame frame(table);
    const std::vector<std::size_t> positions = frame.select(names);
    const std::size_t at = frame.series(denominator);
    require_unique_index(table);
    return ratio_tables(frame, names, positions, frame.cells(at), paired_output(names, suffix));
}

}  // namespace

// ---------------------------------------------------------------------------
// the public functions

namespace {

const std::vector<std::string> kPairsUsed = {"Hi-C contacts",
                                             "Pairs mappable, unique and high quality",
                                             "Sequenced reads"};
const std::vector<std::string> kUnmappable = {"Hi-C contacts", "Low mapping quality",
                                              "One mate not unique", "One mate unmapped"};
const std::vector<std::string> kDistance = {"inter chromosomal", "Intra short range (< 20kb)",
                                            "Intra long range (>= 20kb)"};
const std::vector<std::string> kOrientation = {
    "Read pair type: inward pairs", "Read pair type: outward pairs",
    "Read pair type: left pairs", "Read pair type: right pairs"};

py::tuple pairs_discarded_data(py::handle table) {
    static const char* const prefixes[] = {"One mate not close to rest site", "dangling end",
                                           "duplicated pairs", "same fragment", "self circle",
                                           "self ligation (removed)"};
    const Frame frame(table);
    std::vector<std::string> column_names;
    std::vector<OutputColumn> output;
    for (const char* prefix : prefixes) {
        const std::string name(prefix);
        const std::ptrdiff_t exact = frame.find(name);
        if (exact >= 0) {
            column_names.push_back(name);
            output.push_back({name, column_names.size() - 1, true});
            output.push_back({name + " %", column_names.size() - 1, false});
        } else {
            for (std::size_t i = 0; i < frame.names.size(); ++i) {
                if (!frame.is_name[i]) {
                    raise_error("TypeError", "argument of type '" + type_name(frame.columns[i]) +
                                                 "' is not iterable");
                }
                if (frame.names[i].find(name) != std::string::npos) {
                    column_names.push_back(frame.names[i]);
                    output.push_back({frame.names[i] + " %", column_names.size() - 1, false});
                }
            }
        }
    }
    py::list names;
    for (const std::string& name : column_names) {
        names.append(to_str(name));
    }
    log(kLogger, "debug", "column_names " + str_of(names));
    const std::vector<std::size_t> positions = frame.select(column_names);
    const std::size_t at = frame.series("Pairs mappable, unique and high quality");
    require_unique_index(table);
    return ratio_tables(frame, column_names, positions, frame.cells(at), output);
}

py::tuple read_orientation_data(py::handle table) {
    const Frame frame(table);
    const std::vector<std::size_t> positions = frame.select(kOrientation);
    require_unique_index(table);
    std::vector<Column> columns;
    for (const std::size_t position : positions) {
        columns.push_back(frame.cells(position));
    }
    return ratio_tables(frame, kOrientation, positions, hicx::qcapi::sum_columns(columns),
                        paired_output(kOrientation, " %"));
}

// ---------------------------------------------------------------------------
// save_html

hicx::QcDisplayTable display_table(py::handle frame) {
    hicx::QcDisplayTable out;
    const py::object index = frame.attr("index");
    const py::object columns = frame.attr("columns");
    if (index.attr("nlevels").cast<int>() != 1 || columns.attr("nlevels").cast<int>() != 1) {
        raise_error("TypeError", "save_html: HiCExplorer 4 takes single level indexes");
    }
    const py::object index_name = index.attr("name");
    if (!index_name.is_none()) {
        out.index_name = label_text(index_name);
    }
    const py::object columns_name = columns.attr("name");
    if (!columns_name.is_none()) {
        out.columns_name = label_text(columns_name);
    }
    for (py::handle label : to_list(index)) {
        out.index.push_back(label_text(label));
    }
    const py::list labels = to_list(columns);
    for (py::handle label : labels) {
        out.columns.push_back(label_text(label));
    }
    for (std::size_t c = 0; c < labels.size(); ++c) {
        std::vector<QcStylerCell> cells;
        const py::object series = frame.attr("iloc")[py::make_tuple(
            py::slice(py::none(), py::none(), py::none()), c)];
        for (py::handle value : to_list(series)) {
            cells.push_back(classify(value));
        }
        out.cells.push_back(std::move(cells));
    }
    return out;
}

void save_html(py::handle filename, py::handle unmap_table, py::handle discard_table,
               py::handle distance_table, py::handle orientation_table, py::handle all_table) {
    const hicx::QcDisplayTable unmapped = display_table(unmap_table);
    const hicx::QcDisplayTable discarded = display_table(discard_table);
    const hicx::QcDisplayTable distance = display_table(distance_table);
    const hicx::QcDisplayTable orientation = display_table(orientation_table);
    const hicx::QcDisplayTable all = display_table(all_table);
    const std::string html =
        hicx::qc_html_document(unmapped, discarded, distance, orientation, all);

    const py::object resolved = py::module_::import("os").attr("fspath")(filename);
    const std::string name = fs_path(filename);
    if (name.find('\0') != std::string::npos) {
        raise_error("ValueError", "embedded null byte");
    }
    const int fd = ::open(name.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0666);
    if (fd < 0) {
        raise_os_error(errno, resolved);
    }
    struct stat info {};
    if (::fstat(fd, &info) == 0 && S_ISDIR(info.st_mode)) {
        ::close(fd);
        raise_os_error(EISDIR, resolved);
    }
    std::string bytes;
    try {
        const py::object encoding = py::module_::import("locale").attr("getpreferredencoding")(false);
        bytes = py::bytes(to_str(html).attr("encode")(encoding)).cast<std::string>();
    } catch (...) {
        ::close(fd);
        throw;
    }
    std::size_t written = 0;
    while (written < bytes.size()) {
        const ssize_t n = ::write(fd, bytes.data() + written, bytes.size() - written);
        if (n < 0) {
            if (errno == EINTR) {
                continue;
            }
            const int error = errno;
            ::close(fd);
            raise_os_error(error, py::none());
        }
        written += static_cast<std::size_t>(n);
    }
    ::close(fd);
}

// ---------------------------------------------------------------------------
// make_sure_path_exists: os.makedirs(path) except for EEXIST

struct Split {
    std::string head;
    std::string tail;
};

// posixpath.split
Split split_path(const std::string& p) {
    const std::size_t i = p.rfind('/') == std::string::npos ? 0 : p.rfind('/') + 1;
    Split out{p.substr(0, i), p.substr(i)};
    if (!out.head.empty() && out.head.find_first_not_of('/') != std::string::npos) {
        const std::size_t end = out.head.find_last_not_of('/');
        out.head.resize(end + 1);
    }
    return out;
}

bool path_exists(const std::string& p) {
    struct stat info {};
    return ::stat(p.c_str(), &info) == 0;
}

py::object name_object(const std::string& text, bool bytes_path) {
    return bytes_path ? py::object(py::bytes(text)) : py::object(to_str(text));
}

// os.makedirs(name, 0o777, exist_ok=False); returns 0, or the errno and the
// path of the failing mkdir.
int make_directories(const std::string& name, std::string& failed) {
    Split parts = split_path(name);
    if (parts.tail.empty()) {
        parts = split_path(parts.head);
    }
    if (!parts.head.empty() && !parts.tail.empty() && !path_exists(parts.head)) {
        std::string inner;
        const int error = make_directories(parts.head, inner);
        if (error != 0 && error != EEXIST) {
            failed = inner;
            return error;
        }
        if (parts.tail == ".") {
            return 0;
        }
    }
    if (::mkdir(name.c_str(), 0777) != 0) {
        failed = name;
        return errno;
    }
    return 0;
}

void make_sure_path_exists(py::handle path) {
    const py::object resolved = py::module_::import("os").attr("fspath")(path);
    const std::string name = fs_path(path);
    if (name.find('\0') != std::string::npos) {
        raise_error("ValueError", "embedded null byte");
    }
    const bool bytes_path = PyBytes_Check(resolved.ptr());
    std::string failed;
    const int error = make_directories(name, failed);
    if (error != 0 && error != EEXIST) {
        raise_os_error(error, name_object(failed, bytes_path));
    }
}

}  // namespace

void bind_hicPrepareQCreport(py::module_& m) {
    m.def("save_html", &save_html, py::arg("filename"), py::arg("unmap_table"),
          py::arg("discard_table"), py::arg("distance_table"), py::arg("orientation_table"),
          py::arg("all_table"));
    m.def("make_sure_path_exists", &make_sure_path_exists, py::arg("path"));
    m.def("pairs_used_data", [](py::handle table) {
        return divided_by_scalar(Frame(table), kPairsUsed, 1e6);
    }, py::arg("table"));
    m.def("umappable_non_unique_reads_data", [](py::handle table) {
        return by_column(table, kUnmappable, "Sequenced reads", "_%");
    }, py::arg("table"));
    m.def("pairs_discarded_data", &pairs_discarded_data, py::arg("table"));
    m.def("distance_data", [](py::handle table) {
        return by_column(table, kDistance, "Hi-C contacts", " %");
    }, py::arg("table"));
    m.def("read_orientation_data", &read_orientation_data, py::arg("table"));
}

}  // namespace hicx::python
