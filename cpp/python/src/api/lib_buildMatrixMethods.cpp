// hicexplorer._hicx.lib_buildMatrixMethods: the C++ bodies of the public
// names of hicexplorer.lib.buildMatrixMethods (cpp/python/API_BINDINGS.md).
//
// The computations are hicx::build_matrix (the run of the build matrix
// tools, core/include/hicx/build_matrix_run.hpp), the binning functions of
// hicx/build_matrix.hpp and hicx::build_matrix_api. The functions that take
// pysam objects (get_chrom_sizes, check_dangling_end,
// get_supplementary_alignment, get_correct_map, readBamFiles, process_data)
// keep taking them: the binding reads the attributes the 3.7.6 code read
// (references, lengths, flag, qname, mapq, rname, reference_id, pos, qlen,
// seq, is_reverse, cigartuples, has_tag/get_tag) from whatever object it is
// given and iterates the given iterators with next(), so pysam objects, and
// any object offering the same attributes, work as before; the decisions are
// C++. createMatrix opens the BAM files itself with htslib, by the .name of
// the file objects it is given, as 3.7.6 opened them with pysam by name.
//
// Deviations from 3.7.6 are listed in hicexplorer/test/api/cases_lib_buildMatrixMethods.py.

#include <unistd.h>

#include <cerrno>
#include <cstring>
#include <optional>
#include <string>
#include <vector>

#include "api.hpp"
#include "build_matrix_support.hpp"
#include "hicx/build_matrix.hpp"
#include "hicx/build_matrix_api.hpp"
#include "hicx/build_matrix_run.hpp"
#include "hicx/py_text.hpp"

namespace hicx::python {

namespace api = hicx::build_matrix_api;
namespace bm = hicx::build_matrix;

namespace {

constexpr const char* kLogger = "hicexplorer.lib.buildMatrixMethods";

static_assert(sizeof(hicx::SearchInterval) == 12, "C_Interval is three c_uint");

// --- Python protocol helpers ---------------------------------------------------

bool compare(py::handle a, py::handle b, int op) {
    const int result = PyObject_RichCompareBool(a.ptr(), b.ptr(), op);
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

bool truthy(py::handle value) {
    const int result = PyObject_IsTrue(value.ptr());
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

bool assertions_enabled() {
    return py::module_::import("sys").attr("flags").attr("optimize").cast<int>() == 0;
}

// format(value, ''), what "{}".format(value) calls.
std::string formatted(py::handle value) {
    PyObject* empty = PyUnicode_FromString("");
    if (empty == nullptr) {
        throw_python_error();
    }
    PyObject* text = PyObject_Format(value.ptr(), empty);
    Py_DECREF(empty);
    if (text == nullptr) {
        throw_python_error();
    }
    return utf8(py::reinterpret_steal<py::object>(text));
}

// next(iterator); StopIteration when it is exhausted.
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

// next(iterator), or nullopt for StopIteration.
std::optional<py::object> try_next(py::handle iterator) {
    try {
        return next_of(iterator);
    } catch (py::error_already_set& error) {
        if (error.matches(PyExc_StopIteration)) {
            return std::nullopt;
        }
        throw;
    }
}

// a, b, ... = value with the errors of the unpacking statement.
std::vector<py::object> unpack(py::handle value, std::size_t n) {
    PyObject* iterator = PyObject_GetIter(value.ptr());
    if (iterator == nullptr) {
        if (PyErr_ExceptionMatches(PyExc_TypeError)) {
            PyErr_Clear();
            raise_error("TypeError", "cannot unpack non-iterable " + type_name(value) + " object");
        }
        throw_python_error();
    }
    const py::object it = py::reinterpret_steal<py::object>(iterator);
    std::vector<py::object> items;
    while (items.size() <= n) {
        PyObject* item = PyIter_Next(it.ptr());
        if (item == nullptr) {
            if (PyErr_Occurred()) {
                throw_python_error();
            }
            break;
        }
        items.push_back(py::reinterpret_steal<py::object>(item));
    }
    if (items.size() < n) {
        raise_error("ValueError", "not enough values to unpack (expected " + std::to_string(n) + ", got " +
                                      std::to_string(items.size()) + ")");
    }
    if (items.size() > n) {
        raise_error("ValueError", "too many values to unpack (expected " + std::to_string(n) + ")");
    }
    return items;
}

py::object get_item(py::handle container, py::handle key) {
    PyObject* item = PyObject_GetItem(container.ptr(), key.ptr());
    if (item == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(item);
}

// int(value), Python's own conversion (str and bytes text included)
py::object python_int(py::handle value) {
    PyObject* number = PyNumber_Long(value.ptr());
    if (number == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(number);
}

py::object get_index(py::handle container, std::int64_t index) {
    return get_item(container, py::int_(index));
}

void set_index(py::handle container, std::int64_t index, py::handle value) {
    if (PyObject_SetItem(container.ptr(), py::int_(index).ptr(), value.ptr()) != 0) {
        throw_python_error();
    }
}

std::size_t length(py::handle value) {
    const Py_ssize_t n = PyObject_Length(value.ptr());
    if (n < 0) {
        throw_python_error();
    }
    return static_cast<std::size_t>(n);
}

bool is_numpy_int64(py::handle value) {
    return py::isinstance(value, py::module_::import("numpy").attr("int64"));
}

// An integer coordinate: a Python int or numpy.int64, the types the 3.7.6
// arithmetic keeps; anything else raises TypeError naming HiCExplorer 4.
std::int64_t coordinate(py::handle value, const char* function) {
    if (!(PyLong_Check(value.ptr()) || is_numpy_int64(value))) {
        raise_error("TypeError", std::string(function) +
                                     ": HiCExplorer 4 takes the coordinates as int or numpy.int64, not " +
                                     type_name(value));
    }
    return to_int64(value);
}

// An integer argument (Python int, bool or numpy integer).
std::int64_t integer_argument(py::handle value, const char* function, const char* name) {
    if (!is_integer(value)) {
        raise_error("TypeError", std::string(function) + ": HiCExplorer 4 takes " + name + " as an int, not " +
                                     type_name(value));
    }
    return to_int64(value);
}

std::string str_argument(py::handle value, const char* function, const char* what) {
    if (!is_str(value)) {
        raise_error("TypeError",
                    std::string(function) + ": HiCExplorer 4 takes " + what + " as a str, not " + type_name(value));
    }
    return utf8(value);
}

py::object int64_scalar(std::int64_t value) {
    return numpy_scalar<std::int64_t>(value);
}

py::object user_region_function() {
    return py::module_::import("hicexplorer._hicx.utilities").attr("getUserRegion");
}

[[noreturn]] void raise_system_exit(py::handle code) {
    PyErr_SetObject(PyExc_SystemExit, code.ptr());
    throw_python_error();
}

// --- the functions --------------------------------------------------------------

// check_dangling_end(read, dangling_sequences)
bool check_dangling_end(py::handle read, py::handle ds) {
    const py::str forward_key("pat_forw");
    const py::str reverse_key("pat_rev");
    const auto contains = [&](const py::str& key) {
        const int result = PySequence_Contains(ds.ptr(), key.ptr());
        if (result < 0) {
            throw_python_error();
        }
        return result == 1;
    };
    if (!contains(forward_key) || !contains(reverse_key)) {
        return false;
    }
    // read.seq.upper().startswith(pattern) / .endswith(pattern)
    const auto test = [&](bool prefix, const py::str& key) {
        const py::object seq = read.attr("seq");
        std::string upper;
        const std::string text = is_str(seq) ? utf8(seq) : std::string();
        if (is_str(seq) && std::all_of(text.begin(), text.end(),
                                       [](char c) { return static_cast<unsigned char>(c) < 0x80; })) {
            upper = api::ascii_upper(text);
        } else {
            // str.upper() beyond ASCII, or the AttributeError of None.upper
            upper = str_argument(seq.attr("upper")(), "check_dangling_end", "read.seq");
        }
        const py::object pattern = get_item(ds, key);
        const char* method = prefix ? "startswith" : "endswith";
        const auto matches = [&](const std::string& p) {
            return prefix ? api::starts_with(upper, p) : api::ends_with(upper, p);
        };
        if (is_str(pattern)) {
            return matches(utf8(pattern));
        }
        if (PyTuple_Check(pattern.ptr())) {
            for (py::handle item : pattern) {
                if (!is_str(item)) {
                    raise_error("TypeError", std::string("tuple for ") + method +
                                                 " must only contain str, not " + type_name(item));
                }
            }
            for (py::handle item : pattern) {
                if (matches(utf8(item))) {
                    return true;
                }
            }
            return false;
        }
        raise_error("TypeError", std::string(method) + " first arg must be str or a tuple of str, not " +
                                     type_name(pattern));
    };
    if (!truthy(read.attr("is_reverse")) && test(true, forward_key)) {
        return true;
    }
    if (truthy(read.attr("is_reverse")) && test(false, reverse_key)) {
        return true;
    }
    return false;
}

// get_supplementary_alignment(read, pysam_obj)
py::object get_supplementary_alignment(py::handle read, py::handle pysam_obj) {
    if (!truthy(read.attr("has_tag")("SA"))) {
        return py::none();
    }
    const py::object tag = read.attr("get_tag")("SA");
    if (!is_str(tag)) {
        (void)py::getattr(tag, "split");  // the AttributeError of tag.split
        raise_error("TypeError", "get_supplementary_alignment: HiCExplorer 4 takes the SA tag as a str, not " +
                                     type_name(tag));
    }
    const std::size_t count = api::supplementary_count(utf8(tag));
    py::list supplementary;
    for (std::size_t i = 0; i < count; ++i) {
        const py::object candidate = next_of(pysam_obj);
        if (compare(candidate.attr("qname"), read.attr("qname"), Py_EQ)) {
            supplementary.append(candidate);
        }
    }
    return supplementary;
}

// get_correct_map(primary, supplement_list)
py::object get_correct_map(py::handle primary, py::handle supplement_list) {
    if (assertions_enabled()) {
        for (py::handle supplement : py::iter(supplement_list)) {
            const py::object primary_name = primary.attr("qname");
            const py::object supplement_name = supplement.attr("qname");
            if (!compare(primary_name, supplement_name, Py_EQ)) {
                raise_error("AssertionError",
                            "ERROR, primary and supplementary reads do not have the same id. The ids are as "
                            "follows\n" + formatted(primary_name) + "\n" + formatted(supplement_name));
            }
        }
    }
    py::list head;
    head.append(primary);
    PyObject* joined = PyNumber_Add(head.ptr(), supplement_list.ptr());
    if (joined == nullptr) {
        throw_python_error();
    }
    const py::object read_list = py::reinterpret_steal<py::object>(joined);
    std::vector<std::int64_t> first_mapped;
    std::vector<py::object> reads;
    for (py::handle read : py::iter(read_list)) {
        reads.push_back(py::reinterpret_borrow<py::object>(read));
        const py::object slice = truthy(read.attr("is_reverse"))
                                     ? py::object(py::slice(py::none(), py::none(), py::int_(-1)))
                                     : py::object(py::slice(py::none(), py::none(), py::none()));
        const py::object cigartuples = get_item(read.attr("cigartuples"), slice);
        std::vector<api::CigarOperation> cigar;
        for (py::handle operation : py::iter(cigartuples)) {
            const std::vector<py::object> parts = unpack(operation, 2);
            const std::int64_t op = integer_argument(parts[0], "get_correct_map", "a CIGAR operation");
            if (op == 0) {
                cigar.emplace_back(op, 0);
                break;
            }
            cigar.emplace_back(op, integer_argument(parts[1], "get_correct_map", "a CIGAR length"));
        }
        first_mapped.push_back(api::bases_before_match(cigar));
    }
    return reads[api::first_index_of_minimum(first_mapped)];
}

// intervalListToIntervalTree(interval_list)
py::object interval_list_to_interval_tree(py::handle interval_list) {
    const py::module_ intervaltree = py::module_::import("intervaltree");
    const py::object Interval = intervaltree.attr("Interval");
    const py::object IntervalTree = intervaltree.attr("IntervalTree");
    py::dict trees;
    std::int64_t id = 0;
    for (py::handle interval : py::iter(interval_list)) {
        PyObject* head = PySequence_GetSlice(interval.ptr(), 0, 3);
        if (head == nullptr) {
            throw_python_error();
        }
        const std::vector<py::object> fields = unpack(py::reinterpret_steal<py::object>(head), 3);
        const py::object chrom = fields[0];
        if (!trees.contains(chrom)) {
            trees[chrom] = IntervalTree();
        }
        trees[chrom].attr("add")(Interval(fields[1], fields[2], py::int_(id)));
        ++id;
    }
    return trees;
}

// get_bins(bin_size, chrom_size, region)
py::list get_bins(py::handle bin_size, py::handle chrom_size, py::handle region) {
    py::object sizes = py::reinterpret_borrow<py::object>(chrom_size);
    std::int64_t start = 0;
    if (truthy(region)) {
        const std::vector<py::object> user = unpack(user_region_function()(chrom_size, region), 4);
        sizes = user[0];
        start = integer_argument(user[1], "get_bins", "the region start");
    }
    PyObject* step_index = PyNumber_Index(bin_size.ptr());
    if (step_index == nullptr) {
        throw_python_error();
    }
    Py_DECREF(step_index);
    const bool step_is_int = PyLong_Check(bin_size.ptr());
    if (!step_is_int && !is_numpy_int64(bin_size)) {
        raise_error("TypeError", "get_bins: HiCExplorer 4 takes bin_size as an int or numpy.int64, not " +
                                     type_name(bin_size));
    }
    const std::int64_t step = to_int64(bin_size);
    py::list bins;
    for (py::handle entry : py::iter(sizes)) {
        const std::vector<py::object> fields = unpack(entry, 2);
        PyObject* size_index = PyNumber_Index(fields[1].ptr());
        if (size_index == nullptr) {
            throw_python_error();
        }
        Py_DECREF(size_index);
        const std::int64_t size = coordinate(fields[1], "get_bins");
        const std::optional<std::vector<std::int64_t>> range = api::python_range(start, size, step);
        if (!range) {
            raise_error("ValueError", "range() arg 3 must not be zero");
        }
        for (const std::int64_t interval : *range) {
            const std::int64_t next = interval + step;
            // min(size, interval + bin_size) is size unless the sum is smaller
            const py::object end = next < size ? (step_is_int ? py::object(py::int_(next)) : int64_scalar(next))
                                               : fields[1];
            bins.append(py::make_tuple(fields[0], py::int_(interval), end));
        }
    }
    return bins;
}

// bed2interval_list(bed_file_handler, pChromosomeSize, pRegion)
py::list bed2interval_list(py::handle bed_file_handler, py::handle chromosome_size, py::handle region) {
    py::object region_sizes;
    py::object start_region;
    py::object end_region;
    const bool filtered = !region.is_none();
    if (filtered) {
        const std::vector<py::object> user = unpack(user_region_function()(chromosome_size, region), 4);
        region_sizes = user[0];
        start_region = user[1];
        end_region = user[2];
    }
    std::int64_t count = 0;
    py::list intervals;
    std::optional<std::array<py::object, 3>> fields_bound;  // chrom, start, end of the last parsed line
    for (py::handle line : py::iter(bed_file_handler)) {
        ++count;
        // line.strip().split()
        std::vector<py::object> fields;
        if (is_str(line)) {
            const std::string text = utf8(line);
            for (const std::string_view field : hicx::pytext::split_whitespace(hicx::pytext::strip(text))) {
                fields.push_back(to_str(field));
            }
        } else if (PyBytes_Check(line.ptr())) {
            const std::string text(PyBytes_AS_STRING(line.ptr()), static_cast<std::size_t>(PyBytes_GET_SIZE(line.ptr())));
            std::size_t i = 0;
            const auto space = [](char c) {
                return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\x0b' || c == '\x0c';
            };
            while (i < text.size()) {
                while (i < text.size() && space(text[i])) {
                    ++i;
                }
                std::size_t j = i;
                while (j < text.size() && !space(text[j])) {
                    ++j;
                }
                if (j > i) {
                    fields.push_back(py::bytes(text.substr(i, j - i)));
                }
                i = j;
            }
        } else {
            raise_error("TypeError", "bed2interval_list: HiCExplorer 4 reads lines of str or bytes, not " +
                                         type_name(line));
        }
        // chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        bool index_error = false;
        std::array<py::object, 3> parsed;
        if (fields.empty()) {
            index_error = true;
        } else {
            parsed[0] = fields[0];
            if (fields.size() < 2) {
                index_error = true;
            } else {
                parsed[1] = python_int(fields[1]);
                if (fields.size() < 3) {
                    index_error = true;
                } else {
                    parsed[2] = python_int(fields[2]);
                }
            }
        }
        if (index_error) {
            log(kLogger, "error", "error reading BED file at line " + std::to_string(count));
        } else {
            fields_bound = parsed;
        }
        if (!fields_bound) {
            raise_error("UnboundLocalError",
                        "cannot access local variable 'chrom' where it is not associated with a value");
        }
        const auto& [chrom, start, end] = *fields_bound;
        if (filtered) {
            // chrom == chrom_size[0] compares the str with a (chrom, size) tuple
            if (compare(chrom, get_index(region_sizes, 0), Py_EQ) && compare(start_region, start, Py_LE) &&
                compare(end_region, end, Py_LE)) {
                intervals.append(py::make_tuple(chrom, start, end));
            }
        } else {
            intervals.append(py::make_tuple(chrom, start, end));
        }
    }
    return intervals;
}

// get_rf_bins(rf_cut_intervals, min_distance, max_distance)
py::list get_rf_bins(py::handle rf_cut_intervals, py::handle min_distance, py::handle max_distance) {
    log(kLogger, "info", "Minimum distance considered between restriction sites is " + formatted(min_distance) +
                             "\nMax distance: " + formatted(max_distance) + "\n");
    // chrom, start, end = list(zip(*rf_cut_intervals))
    std::vector<py::object> items;
    for (py::handle item : py::iter(rf_cut_intervals)) {
        items.push_back(py::reinterpret_borrow<py::object>(item));
    }
    std::size_t columns = items.empty() ? 0 : std::numeric_limits<std::size_t>::max();
    std::vector<std::vector<py::object>> rows;
    for (std::size_t k = 0; k < items.size(); ++k) {
        PyObject* iterator = PyObject_GetIter(items[k].ptr());
        if (iterator == nullptr) {
            PyErr_Clear();
            raise_error("TypeError", "zip argument #" + std::to_string(k + 1) + " must support iteration");
        }
        Py_DECREF(iterator);
        std::vector<py::object> row;
        for (py::handle value : py::iter(items[k])) {
            row.push_back(py::reinterpret_borrow<py::object>(value));
        }
        columns = std::min(columns, row.size());
        rows.push_back(std::move(row));
    }
    if (columns < 3) {
        raise_error("ValueError", "not enough values to unpack (expected 3, got " + std::to_string(columns) + ")");
    }
    if (columns > 3) {
        raise_error("ValueError", "too many values to unpack (expected 3)");
    }
    if (min_distance.is_none()) {
        // np.diff(start) - rest_site_len <= None compares element-wise; with a
        // single site the difference is empty and the bin filter fails instead
        raise_error("TypeError", rows.size() > 1 ? "'<=' not supported between instances of 'int' and 'NoneType'"
                                                 : "'>=' not supported between instances of 'int' and 'NoneType'");
    }
    const std::int64_t min_d = integer_argument(min_distance, "get_rf_bins", "min_distance");
    const std::int64_t max_d = integer_argument(max_distance, "get_rf_bins", "max_distance");
    // the chromosomes interned by Python equality, as `chrom[idx] != chrom[idx - 1]` compares them
    py::dict ids;
    std::vector<py::object> chrom_objects;
    std::vector<hicx::GenomeInterval> sites;
    sites.reserve(rows.size());
    for (const auto& row : rows) {
        std::uint32_t id = 0;
        if (ids.contains(row[0])) {
            id = ids[row[0]].cast<std::uint32_t>();
        } else {
            id = static_cast<std::uint32_t>(chrom_objects.size());
            ids[row[0]] = py::int_(id);
            chrom_objects.push_back(row[0]);
        }
        sites.push_back({id, coordinate(row[1], "get_rf_bins"), coordinate(row[2], "get_rf_bins")});
    }
    std::vector<bool> start_is_zero;
    const std::vector<hicx::GenomeInterval> bins = hicx::get_rf_bins(sites, min_d, max_d, &start_is_zero);
    py::list out;
    for (std::size_t i = 0; i < bins.size(); ++i) {
        out.append(py::make_tuple(chrom_objects[bins[i].chrom],
                                  start_is_zero[i] ? py::object(py::int_(0)) : int64_scalar(bins[i].start),
                                  int64_scalar(bins[i].end)));
    }
    return out;
}

// get_chrom_sizes(bam_handle)
py::list get_chrom_sizes(py::handle bam_handle) {
    const py::object references = bam_handle.attr("references");
    const py::object lengths = bam_handle.attr("lengths");
    std::vector<py::object> names;
    std::vector<py::object> values;
    std::vector<std::string> keys;
    const py::object values_iterator = py::reinterpret_steal<py::object>(PyObject_GetIter(lengths.ptr()));
    if (!values_iterator) {
        throw_python_error();
    }
    for (py::handle name : py::iter(references)) {
        const std::optional<py::object> value = try_next(values_iterator);
        if (!value) {
            break;
        }
        keys.push_back(str_argument(name, "get_chrom_sizes", "the reference names"));
        names.push_back(py::reinterpret_borrow<py::object>(name));
        values.push_back(*value);
    }
    py::list out;
    for (const auto& [first, last] : api::ordered_dict_items(keys)) {
        out.append(py::make_tuple(names[first], values[last]));
    }
    return out;
}

// enlarge_bins(bin_intervals, chrom_sizes)
py::object enlarge_bins(py::handle bin_intervals, py::handle chrom_sizes) {
    PyObject* sizes_raw = PyObject_CallOneArg(reinterpret_cast<PyObject*>(&PyDict_Type), chrom_sizes.ptr());
    if (sizes_raw == nullptr) {
        throw_python_error();
    }
    const py::dict sizes = py::reinterpret_steal<py::dict>(sizes_raw);
    const std::size_t n = length(bin_intervals);
    std::vector<std::array<py::object, 3>> items(n);
    std::vector<std::int64_t> starts(n);
    std::vector<std::int64_t> ends(n);
    for (std::size_t i = 0; i < n; ++i) {
        const std::vector<py::object> fields = unpack(get_index(bin_intervals, static_cast<std::int64_t>(i)), 3);
        items[i] = {fields[0], fields[1], fields[2]};
        starts[i] = coordinate(fields[1], "enlarge_bins");
        ends[i] = coordinate(fields[2], "enlarge_bins");
    }
    std::vector<bool> same_chrom(n > 0 ? n - 1 : 0);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        same_chrom[i] = compare(items[i][0], items[i + 1][0], Py_EQ);
    }
    using Kind = api::EnlargeField::Kind;
    const auto field_object = [&](const api::EnlargeField& field, std::size_t chrom_from) -> py::object {
        switch (field.kind) {
            case Kind::OriginalStart: return items[field.from][1];
            case Kind::OriginalEnd: return items[field.from][2];
            case Kind::Zero: return py::int_(0);
            case Kind::Middle:
                // start_next - int(...) has the type of start_next
                return PyLong_Check(items[field.from][1].ptr()) ? py::object(py::int_(field.value))
                                                                 : int64_scalar(field.value);
            case Kind::Size: return get_item(sizes, items[chrom_from][0]);  // KeyError as the Python
        }
        return py::none();
    };
    for (const api::EnlargeAssignment& assignment : api::enlarge_bins(starts, ends, same_chrom)) {
        const py::object start = field_object(assignment.start, assignment.chrom_from);
        const py::object end = field_object(assignment.end, assignment.chrom_from);
        set_index(bin_intervals, static_cast<std::int64_t>(assignment.index),
                  py::make_tuple(items[assignment.chrom_from][0], start, end));
    }
    if (n == 0) {
        get_index(bin_intervals, -1);  // IndexError of bin_intervals[-1]
    }
    return py::reinterpret_borrow<py::object>(bin_intervals);
}

// readBamFiles(...)
py::tuple read_bam_files(py::handle file_one, py::handle file_two, py::handle items_per_buffer,
                         py::handle skip_duplication_check, py::handle read_pos_matrix, py::handle ref_id2name,
                         py::handle min_mapping_quality) {
    py::list buffer_mate1;
    py::list buffer_mate2;
    std::int64_t duplicated_pairs = 0;
    std::int64_t one_mate_unmapped = 0;
    std::int64_t one_mate_not_unique = 0;
    std::int64_t one_mate_low_quality = 0;
    bool all_data_read = false;
    std::int64_t accepted = 0;
    std::int64_t iter_num = 0;
    const auto flag_of = [](py::handle read) {
        return integer_argument(read.attr("flag"), "readBamFiles", "the flag");
    };
    while (compare(py::int_(accepted), items_per_buffer, Py_LT)) {
        std::optional<py::object> mate1 = try_next(file_one);
        std::optional<py::object> mate2;
        if (mate1) {
            mate2 = try_next(file_two);
        }
        if (!mate1 || !mate2) {
            all_data_read = true;
            break;
        }
        ++iter_num;
        // skip 'not primary' alignments
        while ((flag_of(*mate1) & 256) == 256) {
            std::optional<py::object> next = try_next(file_one);
            if (!next) {
                all_data_read = true;
                break;
            }
            mate1 = next;
        }
        while ((flag_of(*mate2) & 256) == 256) {
            std::optional<py::object> next = try_next(file_two);
            if (!next) {
                all_data_read = true;
                break;
            }
            mate2 = next;
        }
        if (assertions_enabled()) {
            const py::object name1 = mate1->attr("qname");
            const py::object name2 = mate2->attr("qname");
            if (!compare(name1, name2, Py_EQ)) {
                raise_error("AssertionError", "FATAL ERROR " + formatted(name1) + " " + formatted(name2) +
                                                  " Be sure that the sam files have the same read order If "
                                                  "using Bowtie2 or Hisat2 add the --reorder option");
            }
        }
        const py::object supplementary1 = get_supplementary_alignment(*mate1, file_one);
        const py::object supplementary2 = get_supplementary_alignment(*mate2, file_two);
        if (truthy(supplementary1)) {
            mate1 = get_correct_map(*mate1, supplementary1);
        }
        if (truthy(supplementary2)) {
            mate2 = get_correct_map(*mate2, supplementary2);
        }
        const std::int64_t flag1 = flag_of(*mate1);
        const std::int64_t flag2 = flag_of(*mate2);
        std::int64_t mapq1 = 0;
        std::int64_t mapq2 = 0;
        std::int64_t minimum = 0;
        if ((flag1 & 0x4) != 4 && (flag2 & 0x4) != 4) {
            mapq1 = integer_argument(mate1->attr("mapq"), "readBamFiles", "the mapping quality");
            mapq2 = integer_argument(mate2->attr("mapq"), "readBamFiles", "the mapping quality");
            minimum = integer_argument(min_mapping_quality, "readBamFiles", "pMinMappingQuality");
        }
        switch (bm::pair_verdict(flag1, flag2, mapq1, mapq2, minimum)) {
            case bm::PairVerdict::OneMateUnmapped: ++one_mate_unmapped; continue;
            case bm::PairVerdict::OneMateNotUnique: ++one_mate_not_unique; continue;
            case bm::PairVerdict::LowQuality: ++one_mate_low_quality; continue;
            case bm::PairVerdict::Keep: break;
        }
        // `if pSkipDuplicationCheck is False`: an identity test
        if (skip_duplication_check.ptr() == Py_False) {
            const py::object duplicated = read_pos_matrix.attr("is_duplicated")(
                get_item(ref_id2name, mate1->attr("rname")), mate1->attr("pos"),
                get_item(ref_id2name, mate2->attr("rname")), mate2->attr("pos"));
            if (truthy(duplicated)) {
                ++duplicated_pairs;
                continue;
            }
        }
        buffer_mate1.append(*mate1);
        buffer_mate2.append(*mate2);
        ++accepted;
    }
    const py::int_ discarded(iter_num - static_cast<std::int64_t>(buffer_mate1.size()));
    const auto counts = [&](py::object first, py::object second, bool done) {
        return py::make_tuple(first, second, py::bool_(done), py::int_(duplicated_pairs), py::int_(one_mate_unmapped),
                              py::int_(one_mate_not_unique), py::int_(one_mate_low_quality), discarded);
    };
    if (all_data_read && buffer_mate1.size() != 0 && buffer_mate2.size() != 0) {
        return counts(buffer_mate1, buffer_mate2, true);
    }
    if ((all_data_read && buffer_mate1.size() == 0) || buffer_mate2.size() == 0) {
        return counts(py::none(), py::none(), true);
    }
    return counts(buffer_mate1, buffer_mate2, false);
}

}  // namespace

}  // namespace hicx::python

// process_data and createMatrix
#include "lib_buildMatrixMethods_run.hpp"

namespace hicx::python {

void bind_lib_buildMatrixMethods(py::module_& m) {
    m.doc() = "C++ bodies of hicexplorer.lib.buildMatrixMethods.";

    // The ctypes layouts of C_Interval (hicx::SearchInterval) and C_Coverage
    // (the coverage begin and end of a bin), the bases of the two classes.
    {
        const py::module_ ctypes = py::module_::import("ctypes");
        const py::object structure_type = py::module_::import("builtins").attr("type");
        const auto layout = [&](const char* name, const std::vector<const char*>& fields, const char* doc) {
            py::list field_list;
            for (const char* field : fields) {
                field_list.append(py::make_tuple(field, ctypes.attr("c_uint")));
            }
            py::dict body;
            body["_fields_"] = field_list;
            body["__module__"] = py::str("hicexplorer._hicx.lib_buildMatrixMethods");
            body["__doc__"] = py::str(doc);
            return structure_type(name, py::make_tuple(ctypes.attr("Structure")), body);
        };
        m.attr("C_Interval") = layout("C_Interval", {"begin", "end", "data"},
                                      "The C++ layout of an interval of the shared bin array: three c_uint.");
        m.attr("C_Coverage") = layout("C_Coverage", {"begin", "end"},
                                      "The C++ layout of a coverage range: two c_uint.");
    }

    m.def("ReadPositionMatrix___init__", [](py::handle self) {
        self.attr("pos_matrix") = py::set();
    }, py::arg("self"));

    m.def("ReadPositionMatrix_is_duplicated", [](py::handle self, py::handle chrom1, py::handle start1,
                                                 py::handle chrom2, py::handle start2) {
        const bool chrom1_first = compare(chrom1, chrom2, Py_LT);
        const std::string c1 = formatted(chrom1);
        const std::string c2 = formatted(chrom2);
        const bool start1_first = compare(start1, start2, Py_LT);
        const std::string s1 = formatted(start1);
        const std::string s2 = formatted(start2);
        const py::str key = to_str(api::read_position_key(c1, c2, chrom1_first, s1, s2, start1_first));
        const py::object positions = self.attr("pos_matrix");
        const int found = PySequence_Contains(positions.ptr(), key.ptr());
        if (found < 0) {
            throw_python_error();
        }
        if (found == 1) {
            return true;
        }
        positions.attr("add")(key);
        return false;
    }, py::arg("self"), py::arg("chrom1"), py::arg("start1"), py::arg("chrom2"), py::arg("start2"));

    m.def("intervalListToIntervalTree", &interval_list_to_interval_tree, py::arg("interval_list"));
    m.def("get_bins", &get_bins, py::arg("bin_size"), py::arg("chrom_size"), py::arg("region") = py::none());
    m.def("bed2interval_list", &bed2interval_list, py::arg("bed_file_handler"), py::arg("pChromosomeSize"),
          py::arg("pRegion"));
    m.def("get_rf_bins", &get_rf_bins, py::arg("rf_cut_intervals"), py::arg("min_distance") = 200,
          py::arg("max_distance") = 800);
    m.def("get_chrom_sizes", &get_chrom_sizes, py::arg("bam_handle"));
    m.def("check_dangling_end", &check_dangling_end, py::arg("read"), py::arg("dangling_sequences"));
    m.def("get_supplementary_alignment", &get_supplementary_alignment, py::arg("read"), py::arg("pysam_obj"));
    m.def("get_correct_map", &get_correct_map, py::arg("primary"), py::arg("supplement_list"));
    m.def("enlarge_bins", &enlarge_bins, py::arg("bin_intervals"), py::arg("chrom_sizes"));
    m.def("readBamFiles", &read_bam_files, py::arg("pFileOneIterator"), py::arg("pFileTwoIterator"),
          py::arg("pNumberOfItemsPerBuffer"), py::arg("pSkipDuplicationCheck"), py::arg("pReadPosMatrix"),
          py::arg("pRefId2name"), py::arg("pMinMappingQuality"));

    // root is dirname(dirname(abspath(__file__))) + "/test/test_data/": two
    // levels up from hicexplorer/lib/buildMatrixMethods.py.
    m.def("Tester___init__", [](py::handle self, py::handle module_file) {
        tester_init(self, module_file, 2);
    }, py::arg("self"), py::arg("module_file"));

    bind_lib_buildMatrixMethods_process(m);
}

}  // namespace hicx::python
