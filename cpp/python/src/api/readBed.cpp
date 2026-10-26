// hicexplorer._hicx.readBed: the C++ bodies of the methods of
// hicexplorer.readBed.ReadBed.
//
// The Python class keeps its attributes (file_type, file_handle, line_number,
// prev_chrom, prev_start, prev_line, fields, BedInterval) and its methods;
// each method calls the function here with `self`. Calls from one method to
// another go through `self`, as in the Python, so a subclass that overrides a
// method still sees it called. The parsing rules are hicx::readbed
// (cpp/core/include/hicx/read_bed.hpp).

#include <optional>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicx/py_text.hpp"
#include "hicx/read_bed.hpp"
#include "utilities.hpp"

namespace hicx::python {

namespace rb = hicx::readbed;

namespace {

constexpr const char* kLogger = "hicexplorer.readBed";

bool compare(py::handle a, py::handle b, int op) {
    const int result = PyObject_RichCompareBool(a.ptr(), b.ptr(), op);
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

bool equal(py::handle a, py::handle b) {
    return compare(a, b, Py_EQ);
}

// `assert` statements run unless the interpreter runs with -O.
bool assertions_enabled() {
    return py::module_::import("sys").attr("flags").attr("optimize").cast<int>() == 0;
}

// next(iterator) with next()'s errors.
py::object next_of(py::handle iterator) {
    if (!PyIter_Check(iterator.ptr())) {
        raise_error("TypeError", "'" + type_name(iterator) + "' object is not an iterator");
    }
    PyObject* item = PyIter_Next(iterator.ptr());
    if (item == nullptr) {
        if (!PyErr_Occurred()) {
            PyErr_SetNone(PyExc_StopIteration);
        }
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(item);
}

// line.split('\t') of a str, or through the object's own method.
py::list split_tabs(const py::object& line) {
    if (is_str(line)) {
        py::list out;
        for (const std::string& field : rb::split_line(utf8(line), false)) {
            out.append(to_str(field));
        }
        return out;
    }
    return line.attr("split")("\t");
}

bool is_comment(const py::object& line) {
    if (is_str(line)) {
        return rb::is_comment_line(utf8(line));
    }
    auto starts = [&](const char* prefix) { return line.attr("startswith")(prefix).cast<bool>(); };
    return starts("#") || starts("track") || starts("browser") ||
           equal(line.attr("strip")(), py::str(""));
}

// int(text); nullopt, with the ValueError's message, when it does not parse.
std::optional<py::object> try_int(const py::object& text, std::string* error_message) {
    PyObject* value = PyLong_FromUnicodeObject(text.ptr(), 10);
    if (value == nullptr) {
        if (!PyErr_ExceptionMatches(PyExc_ValueError)) {
            throw_python_error();
        }
        py::error_already_set error;
        if (error_message != nullptr) {
            *error_message = str_of(error.value());
        }
        return std::nullopt;
    }
    return py::reinterpret_steal<py::object>(value);
}

// One field of get_bed_interval; `failed` for an integer column that does
// not parse, where the Python returns an empty dict.
py::object convert_field(std::size_t index, const py::object& field, const py::object& self,
                         const py::object& bed_line, bool& failed) {
    failed = false;
    switch (rb::column_kind(index)) {
        case rb::Column::Text:
            return field;
        case rb::Column::Strand: {
            bool invalid = false;
            const std::string strand = rb::normalize_strand(utf8(field), invalid);
            if (invalid) {
                log(kLogger, "warning",
                    rb::invalid_strand_message(utf8(field), str_of(bed_line), str_of(self.attr("line_number"))));
            }
            return to_str(strand);
        }
        case rb::Column::Integer: {
            const std::optional<py::object> value = try_int(field, nullptr);
            if (!value) {
                log(kLogger, "warning",
                    rb::not_integer_message(utf8(field), index + 1, str_of(self.attr("line_number"))));
                failed = true;
                return py::none();
            }
            return *value;
        }
        case rb::Column::Rgb: {
            const py::object r = to_string(field);
            const std::vector<std::string> parts = rb::split_commas(utf8(r), false);
            if (parts.size() != 3) {
                return r;
            }
            py::list numbers;
            for (const std::string& part : parts) {
                std::string detail;
                const std::optional<py::object> value = try_int(to_str(part), &detail);
                if (!value) {
                    log(kLogger, "debug", "Error reading line: #" + str_of(self.attr("line_number")) +
                                              ". The rgb field " + utf8(r) +
                                              " is not valid.\nError message: " + detail + "\n");
                    return r;
                }
                numbers.append(*value);
            }
            return numbers;
        }
        case rb::Column::Blocks: {
            const py::object r = to_string(field);
            py::list numbers;
            for (const std::string& part : rb::split_commas(utf8(r), true)) {
                std::string detail;
                const std::optional<py::object> value = try_int(to_str(part), &detail);
                if (!value) {
                    log(kLogger, "debug", "Error reading line #" + str_of(self.attr("line_number")) +
                                              ". The block field " + utf8(r) +
                                              " is not valid.\nError message: " + detail + "\n");
                    return r;
                }
                numbers.append(*value);
            }
            return numbers;
        }
        case rb::Column::Float:
        default: {
            PyObject* value = PyFloat_FromString(field.ptr());
            if (value == nullptr) {
                if (!PyErr_ExceptionMatches(PyExc_ValueError) && !PyErr_ExceptionMatches(PyExc_TypeError)) {
                    throw_python_error();
                }
                PyErr_Clear();
                return field;
            }
            return py::reinterpret_steal<py::object>(value);
        }
    }
}

}  // namespace

void bind_readBed(py::module_& m) {
    m.doc() = "C++ bodies of the methods of hicexplorer.readBed.ReadBed.";

    m.def("ReadBed___init__", [](py::object self, py::object file_handle) {
        self.attr("file_type") = py::none();
        self.attr("file_handle") = file_handle;
        self.attr("line_number") = py::int_(0);
        py::object fields = self.attr("get_no_comment_line")();
        fields = to_string(fields);
        fields = split_tabs(fields);
        self.attr("guess_file_type")(fields);
        self.attr("file_handle").attr("seek")(0);
        self.attr("prev_chrom") = py::none();
        self.attr("prev_start") = py::int_(-1);
        self.attr("prev_line") = py::none();
        py::list names;
        for (const std::string& name : rb::interval_fields("bed12")) {
            names.append(py::str(name));
        }
        self.attr("fields") = names;
        const py::object file_type = self.attr("file_type");
        py::ssize_t count = 6;
        if (equal(file_type, py::str("bed12"))) {
            count = 12;
        } else if (equal(file_type, py::str("bed9"))) {
            count = 9;
        }
        const py::object selected = names[py::slice(0, count, 1)];
        self.attr("BedInterval") = py::module_::import("collections").attr("namedtuple")("BedInterval", selected);
    }, py::arg("self"), py::arg("file_handle"));

    m.def("ReadBed_get_no_comment_line", [](py::object self) {
        std::size_t consumed = 0;
        py::object line;
        while (true) {
            line = to_string(next_of(self.attr("file_handle")));
            ++consumed;
            if (!is_comment(line)) {
                break;
            }
        }
        // one `self.line_number += 1` per line read, after the recursion returned
        for (std::size_t k = 0; k < consumed; ++k) {
            const py::object current = self.attr("line_number");
            PyObject* incremented = PyNumber_InPlaceAdd(current.ptr(), py::int_(1).ptr());
            if (incremented == nullptr) {
                throw_python_error();
            }
            self.attr("line_number") = py::reinterpret_steal<py::object>(incremented);
        }
        return line;
    }, py::arg("self"));

    m.def("ReadBed_guess_file_type", [](py::object self, py::handle line_values) {
        const Py_ssize_t n = PyObject_Length(line_values.ptr());
        if (n < 0) {
            throw_python_error();
        }
        std::string assumed;
        const std::string type = rb::guess_file_type(static_cast<std::size_t>(n), &assumed);
        if (!assumed.empty()) {
            log(kLogger, "debug", "Number of fields in BED file is not standard. Assuming " + assumed + ".");
        }
        self.attr("file_type") = py::str(type);
        return self.attr("file_type");
    }, py::arg("self"), py::arg("line_values"));

    m.def("ReadBed_next", [](py::object self) {
        const py::object line = self.attr("get_no_comment_line")();
        const py::object bed = self.attr("get_bed_interval")(line);
        if (equal(self.attr("prev_chrom"), bed.attr("chromosome")) && assertions_enabled()) {
            if (!compare(self.attr("prev_start"), bed.attr("start"), Py_LE)) {
                raise_error("AssertionError",
                            rb::unsorted_message(str_of(self.attr("file_handle").attr("name")),
                                                 str_of(self.attr("prev_line")), str_of(line)));
            }
        }
        self.attr("prev_chrom") = bed.attr("chromosome");
        self.attr("prev_start") = bed.attr("start");
        self.attr("prev_line") = line;
        return bed;
    }, py::arg("self"));

    m.def("ReadBed_get_bed_interval", [](py::object self, py::object bed_line) -> py::object {
        py::object stripped = is_str(bed_line) ? py::object(to_str(hicx::pytext::strip(utf8(bed_line))))
                                               : py::object(bed_line.attr("strip")());
        stripped = to_string(stripped);
        const py::list line_data = split_tabs(stripped);
        const std::size_t n = line_data.size();
        const bool checks = assertions_enabled();
        if (checks) {
            // `self.file_handle == 'bed12'`: the Python compares the handle
            std::optional<std::size_t> required;
            std::string file_type_text;
            if (equal(self.attr("file_handle"), py::str("bed12"))) {
                required = 12;
                file_type_text = "bed12";
            } else {
                const py::object file_type = self.attr("file_type");
                for (const char* candidate : {"bed3", "bed6"}) {
                    if (equal(file_type, py::str(candidate))) {
                        required = rb::required_field_count(candidate);
                        file_type_text = candidate;
                        break;
                    }
                }
            }
            if (required && n != *required) {
                raise_error("AssertionError",
                            rb::field_count_message(file_type_text, *required,
                                                    str_of(self.attr("line_number")), str_of(bed_line)));
            }
        }
        py::list line_values;
        for (std::size_t index = 0; index < n; ++index) {
            bool failed = false;
            const py::object field = line_data[index];
            py::object value = convert_field(index, field, self, bed_line, failed);
            if (failed) {
                return py::dict();
            }
            line_values.append(value);
        }
        if (checks) {
            if (line_values.size() < 3) {
                raise_error("IndexError", "list index out of range");
            }
            if (!compare(line_values[2], line_values[1], Py_GT)) {
                raise_error("AssertionError",
                            rb::start_end_message(str_of(self.attr("line_number")), str_of(bed_line)));
            }
        }
        const py::object file_type = self.attr("file_type");
        py::list fields = line_values;
        if (equal(file_type, py::str("bed3"))) {
            fields = line_values[py::slice(0, 3, 1)];
            fields.append(py::str("."));
            fields.append(py::int_(0));
            fields.append(py::str("."));
        } else if (equal(file_type, py::str("bed6"))) {
            fields = line_values[py::slice(0, 6, 1)];
        }
        return self.attr("BedInterval").attr("_make")(fields);
    }, py::arg("self"), py::arg("bed_line"));
}

}  // namespace hicx::python
