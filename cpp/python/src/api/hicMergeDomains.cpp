// hicexplorer._hicx.hicMergeDomains: the C++ bodies of hicexplorer.hicMergeDomains.
//
// The functions work on the lists of lists of str the 3.7.6 functions took
// and returned. The algorithms are hicx::merge_domains
// (cpp/core/include/hicx/merge_domains.hpp), the code of the hicMergeDomains
// tool, which keeps the rows of a list as references to row objects; this
// file maps the Python row objects to those references by identity, so a row
// that appears twice in a list stays one object, and maps the results back
// to the same objects. Lists the 3.7.6 functions changed in place (add_id,
// add_relation_to_list, compare_boundaries_protein) are changed in place.
//
// Deviations from 3.7.6, recorded in hicexplorer/test/api/cases_hicMergeDomains.py:
//  * create_tree runs graphviz's dot itself with the command graphviz 0.20.3
//    runs (dot -Kdot -T<format> -O <file>), on a DOT source byte identical to
//    Digraph.source; it does not go through the graphviz package, so
//    patching graphviz.Digraph does not reach it.
//  * merge_protein raises ValueError for a bin size that is not positive,
//    where 3.7.6 loops forever.
//  * Row fields that are not str, and numbers that are not whole where the
//    Python compares them with int(), raise TypeError naming HiCExplorer 4.

#include <algorithm>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "api.hpp"
#include "hicx/merge_domains.hpp"
#include "hicx/py_text.hpp"

namespace hicx::python {

namespace md = hicx::merge_domains;

namespace {

std::size_t length(py::handle object) {
    const Py_ssize_t n = PyObject_Length(object.ptr());
    if (n < 0) {
        throw_python_error();
    }
    return static_cast<std::size_t>(n);
}

py::object item(py::handle sequence, std::int64_t index) {
    PyObject* value = PySequence_GetItem(sequence.ptr(), static_cast<Py_ssize_t>(index));
    if (value == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(value);
}

bool equal(py::handle a, py::handle b) {
    const int result = PyObject_RichCompareBool(a.ptr(), b.ptr(), Py_EQ);
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

bool truthy(py::handle object) {
    const int result = PyObject_IsTrue(object.ptr());
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

// md::PythonError carries "Type: message" for the builtin exceptions the
// Python raises; raise them as that type. Anything else (a failing dot run)
// is a RuntimeError with the whole message.
[[noreturn]] void raise_core_error(const std::exception& error) {
    const std::string what = error.what();
    const std::size_t colon = what.find(": ");
    if (colon != std::string::npos && colon >= 5) {
        const std::string type = what.substr(0, colon);
        const bool name = std::all_of(type.begin(), type.end(), [](char c) {
            return std::isalnum(static_cast<unsigned char>(c)) != 0;
        });
        if (name && type.size() > 5 && type.compare(type.size() - 5, 5, "Error") == 0) {
            raise_error(type, what.substr(colon + 2));
        }
    }
    raise_error("RuntimeError", what);
}

template <class F>
decltype(auto) guarded(F&& f) {
    try {
        return f();
    } catch (const md::ReferenceNeverTerminates& error) {
        raise_error("ValueError", std::string("HiCExplorer 4 refuses an input on which 3.7.6 never returns: ") +
                                      error.what());
    } catch (const md::PythonError& error) {
        raise_core_error(error);
    }
}

std::string text_of(py::handle value, const char* what) {
    if (!is_str(value)) {
        raise_error("TypeError", std::string("HiCExplorer 4 needs str ") + what + ", got " + type_name(value));
    }
    return utf8(value);
}

// str.rstrip() over Python's whitespace.
std::string rstrip(const std::string& text) {
    std::size_t end = text.size();
    while (end > 0) {
        const std::size_t width = hicx::pytext::previous_code_point_length(text, end);
        if (hicx::pytext::whitespace_length(text, end - width) == 0) {
            break;
        }
        end -= width;
    }
    return text.substr(0, end);
}

// [line.rstrip() for line in open(file)]
std::vector<std::string> stripped_lines(py::handle file) {
    std::vector<std::string> lines = read_text_lines(file);
    for (std::string& line : lines) {
        line = rstrip(line);
    }
    return lines;
}

py::list str_list(const std::vector<std::string>& values) {
    py::list out;
    for (const std::string& v : values) {
        out.append(to_str(v));
    }
    return out;
}

// The rows of Python lists, as references into one pool. A row object that
// occurs in several places is one pool row.
struct Rows {
    md::RowPool pool;
    std::vector<py::object> objects;  // the Python row of each pool row
    std::unordered_map<PyObject*, std::size_t> index;

    std::size_t add(py::handle row) {
        const auto found = index.find(row.ptr());
        if (found != index.end()) {
            return found->second;
        }
        md::Row fields;
        for (const py::handle field : py::reinterpret_borrow<py::iterable>(row)) {
            fields.push_back(text_of(field, "fields in a domain row"));
        }
        pool.rows.push_back(std::move(fields));
        objects.push_back(py::reinterpret_borrow<py::object>(row));
        index.emplace(row.ptr(), pool.rows.size() - 1);
        return pool.rows.size() - 1;
    }

    std::vector<std::size_t> list(py::handle rows) {
        std::vector<std::size_t> out;
        for (const py::handle row : py::reinterpret_borrow<py::iterable>(rows)) {
            out.push_back(add(row));
        }
        return out;
    }

    py::list objects_of(const std::vector<std::size_t>& rows) const {
        py::list out;
        for (const std::size_t r : rows) {
            out.append(objects[r]);
        }
        return out;
    }
};

// cList as merge_protein returns it: per chromosome [chrom, left, right, count].
md::MergedProtein merged_protein(py::handle c_list) {
    md::MergedProtein out;
    for (const py::handle chromosome : py::reinterpret_borrow<py::iterable>(c_list)) {
        out.emplace_back();
        for (const py::handle peak : py::reinterpret_borrow<py::iterable>(chromosome)) {
            md::ProteinBin bin;
            bin.chrom = text_of(item(peak, 0), "chromosome names");
            bin.left = to_int64(item(peak, 1));
            bin.right = to_int64(item(peak, 2));
            bin.count = to_int64(item(peak, 3));
            out.back().push_back(std::move(bin));
        }
    }
    return out;
}

// --------------------------------------------------------------------------

py::object create_list_of_file(py::handle file, py::handle /*binSizeList*/) {
    const std::vector<std::string> lines = stripped_lines(file);
    md::RowPool pool;
    const md::DomainList list = guarded([&] { return md::create_list_of_lines(pool, lines); });
    py::list rows;
    for (const std::size_t r : list.rows) {
        rows.append(str_list(pool.rows[r]));
    }
    // log.debug('binSize {}'.format(binSize))
    log("hicexplorer.hicMergeDomains", "debug", "binSize " + std::to_string(list.bin_size));
    return py::make_tuple(rows, list.bin_size);
}

py::object merge_list(py::handle d1, py::handle d2, py::handle pValue) {
    Rows rows;
    const std::vector<std::size_t> first = rows.list(d1);
    const std::vector<std::size_t> second = rows.list(d2);
    const std::int64_t value = to_int64(pValue);
    const std::vector<std::size_t> merged =
        guarded([&] { return md::merge_list(rows.pool, first, second, value); });
    return rows.objects_of(merged);
}

py::object add_id(py::handle domainList) {
    std::int64_t id = 1;
    for (const py::handle tad : py::reinterpret_borrow<py::iterable>(domainList)) {
        const py::str value = to_str("ID_" + std::to_string(id));
        if (PySequence_SetItem(tad.ptr(), 3, value.ptr()) != 0) {
            throw_python_error();
        }
        ++id;
    }
    return py::reinterpret_borrow<py::object>(domainList);
}

py::object create_relationsship_list(py::handle domainList, py::handle pPercent) {
    Rows rows;
    const std::vector<std::size_t> list = rows.list(domainList);
    const double percent = to_double(pPercent);
    const std::vector<md::Relation> relations =
        guarded([&] { return md::create_relationship_list(rows.pool, list, percent); });
    py::list out;
    for (const md::Relation& relation : relations) {
        py::list entry;
        entry.append(to_str(relation.chrom));
        entry.append(to_str(relation.parent));
        entry.append(str_list(relation.children));
        out.append(entry);
    }
    return out;
}

std::int64_t id_number(py::handle id) {
    const std::string text = text_of(id, "TAD ids");
    // int(id[3:]), [3:] counted in code points
    std::size_t at = 0;
    for (int k = 0; k < 3 && at < text.size(); ++k) {
        const unsigned char c = static_cast<unsigned char>(text[at]);
        at += c < 0x80 ? 1 : c < 0xE0 ? 2 : c < 0xF0 ? 3 : 4;
    }
    const std::string rest = at < text.size() ? text.substr(at) : std::string();
    return guarded([&] { return md::py_int(rest); });
}

py::object add_relation_to_list(py::handle rList, py::handle chromosom, py::handle parent, py::handle child) {
    py::list relations = py::reinterpret_borrow<py::list>(rList);
    const auto new_entry = [&]() {
        py::list children;
        children.append(child);
        py::list entry;
        entry.append(chromosom);
        entry.append(parent);
        entry.append(children);
        return entry;
    };
    if (length(rList) == 0) {
        relations.append(new_entry());
        return relations;
    }
    std::int64_t pos = static_cast<std::int64_t>(length(rList)) - 1;
    while (true) {
        const py::object entry = item(rList, pos);
        if (equal(item(entry, 1), parent)) {
            item(entry, 2).attr("append")(child);
            return relations;
        }
        if (id_number(item(entry, 1)) < id_number(parent)) {
            relations.append(new_entry());
            return relations;
        }
        --pos;
    }
}

void write_in_file(py::handle pList, py::handle name, py::handle relation) {
    const py::object handle = py::module_::import("builtins").attr("open")(name, "w");
    std::string text;
    const auto finish = [&]() {
        handle.attr("write")(to_str(text));
        handle.attr("close")();
    };
    try {
        const bool relations = truthy(relation);
        const std::size_t count = length(pList);
        for (std::size_t i = 0; i < count; ++i) {
            const py::object row = item(pList, static_cast<std::int64_t>(i));
            if (relations) {
                const py::object children = item(row, 2);
                for (std::size_t e = 0; e < length(children); ++e) {
                    text += text_of(item(row, 0), "chromosome names") + "\t" + text_of(item(row, 1), "TAD ids") +
                            "\t" + text_of(item(children, static_cast<std::int64_t>(e)), "TAD ids") + "\n";
                }
            } else {
                const std::int64_t fields = static_cast<std::int64_t>(length(row));
                std::string line;
                std::int64_t element = 0;
                while (element < fields - 1) {
                    line += text_of(item(row, element), "fields") + "\t";
                    ++element;
                }
                line += text_of(item(row, element), "fields") + "\n";
                text += line;
            }
        }
    } catch (...) {
        finish();
        throw;
    }
    finish();
}

py::object create_small_list(py::handle dList) {
    py::list small;
    const py::object first = item(dList, 0);
    py::list head;
    head.append(item(first, 0));
    head.append(py::list());
    small.append(head);
    py::object chrom = item(first, 0);
    const std::size_t count = length(dList);
    for (std::size_t pos = 0; pos < count; ++pos) {
        const py::object row = item(dList, static_cast<std::int64_t>(pos));
        if (equal(item(row, 0), chrom)) {
            item(item(small, static_cast<std::int64_t>(length(small)) - 1), 1).attr("append")(item(row, 3));
        } else {
            chrom = item(row, 0);
            py::list ids;
            ids.append(item(row, 3));
            py::list entry;
            entry.append(item(row, 0));
            entry.append(ids);
            small.append(entry);
        }
    }
    return small;
}

void create_tree(py::handle rList, py::handle dList, py::handle pFileNameSuffix, py::handle pFormat) {
    // name = pFileNameSuffix + '_' + rList[0][0]; Digraph(..., format=pFormat)
    const std::string prefix = text_of(pFileNameSuffix, "file name prefixes");
    const py::object first_relation = item(rList, 0);
    (void)text_of(item(first_relation, 0), "chromosome names");
    const std::string format = text_of(pFormat, "formats");
    if (!md::is_graphviz_format(format)) {
        std::string formats;
        for (const std::string& f : md::graphviz_formats()) {
            formats += (formats.empty() ? "'" : ", '") + f + "'";
        }
        raise_error("ValueError", "unknown format: " + repr_of(pFormat.attr("lower")()) + " (must be one of [" +
                                      formats + "])");
    }
    std::vector<md::Relation> relations;
    for (const py::handle entry : py::reinterpret_borrow<py::iterable>(rList)) {
        md::Relation relation;
        relation.chrom = text_of(item(entry, 0), "chromosome names");
        relation.parent = text_of(item(entry, 1), "TAD ids");
        for (const py::handle child : py::reinterpret_borrow<py::iterable>(item(entry, 2))) {
            relation.children.push_back(text_of(child, "TAD ids"));
        }
        relations.push_back(std::move(relation));
    }
    Rows rows;
    const std::vector<std::size_t> list = rows.list(dList);
    const py::object print = py::module_::import("builtins").attr("print");
    guarded([&] {
        md::create_tree(relations, rows.pool, list, prefix, [&](const md::TreeGraph& graph) {
            print(to_str("Saved relation tree of " + graph.chrom));
            if (!md::find_on_path("dot").has_value()) {
                // graphviz raises ExecutableNotFound([PosixPath('dot'), ...])
                const py::object error = py::module_::import("graphviz.backend.execute").attr("ExecutableNotFound")(
                    py::make_tuple(py::module_::import("pathlib").attr("PosixPath")("dot")));
                PyErr_SetObject(reinterpret_cast<PyObject*>(Py_TYPE(error.ptr())), error.ptr());
                throw_python_error();
            }
            md::render_graph(graph, format);
        });
        return 0;
    });
}

py::object read_protein(py::handle file) {
    const std::vector<std::string> lines = stripped_lines(file);
    const md::ProteinList proteins = guarded([&] { return md::read_protein_lines(lines); });
    py::list out;
    for (const auto& chromosome : proteins) {
        py::list rows;
        for (const md::Row& row : chromosome) {
            rows.append(str_list(row));
        }
        out.append(rows);
    }
    return out;
}

py::object merge_protein(py::handle pProteinList, py::handle binSize, py::handle minPeak) {
    md::ProteinList proteins;
    for (const py::handle chromosome : py::reinterpret_borrow<py::iterable>(pProteinList)) {
        proteins.emplace_back();
        for (const py::handle peak : py::reinterpret_borrow<py::iterable>(chromosome)) {
            md::Row row;
            for (const py::handle field : py::reinterpret_borrow<py::iterable>(peak)) {
                row.push_back(text_of(field, "fields of a protein peak"));
            }
            proteins.back().push_back(std::move(row));
        }
    }
    const std::int64_t bin_size = to_int64(binSize);
    const std::int64_t min_peak = to_int64(minPeak);
    const md::MergedProtein merged = guarded([&] { return md::merge_protein(proteins, bin_size, min_peak); });
    py::list out;
    for (const auto& chromosome : merged) {
        py::list bins;
        for (const md::ProteinBin& bin : chromosome) {
            py::list entry;
            entry.append(to_str(bin.chrom));
            entry.append(bin.left);
            entry.append(bin.right);
            entry.append(bin.count);
            bins.append(entry);
        }
        out.append(bins);
    }
    return out;
}

py::object compare_boundaries_protein(py::handle bList, py::handle cList, py::handle paraScore) {
    Rows rows;
    std::vector<std::size_t> list = rows.list(bList);
    const md::MergedProtein merged = merged_protein(cList);
    const double score = to_double(paraScore);
    const auto write_back = [&]() {
        const py::list kept = rows.objects_of(list);
        if (PySequence_SetSlice(bList.ptr(), 0, PY_SSIZE_T_MAX, kept.ptr()) != 0) {
            throw_python_error();
        }
    };
    try {
        md::compare_boundaries_protein(rows.pool, list, merged, score);
    } catch (const md::PythonError& error) {
        write_back();
        raise_core_error(error);
    }
    write_back();
    return py::reinterpret_borrow<py::object>(bList);
}

py::object create_list_with_protein(py::handle bList, py::handle minPeak, py::handle cList) {
    const py::module_ module = py::module_::import("hicexplorer.hicMergeDomains");
    const py::tuple loaded = module.attr("create_list_of_file")(bList, 0);
    py::object rows = loaded[0];
    const py::object bin_size = loaded[1];
    if (!cList.is_none()) {
        const py::object merged = module.attr("merge_protein")(cList, bin_size, minPeak);
        rows = module.attr("compare_boundaries_protein")(rows, merged);
    }
    return rows;
}

}  // namespace

void bind_hicMergeDomains(py::module_& m) {
    m.def("create_list_of_file", &create_list_of_file, py::arg("file"), py::arg("binSizeList"));
    m.def("merge_list", &merge_list, py::arg("d1"), py::arg("d2"), py::arg("pValue") = 5000);
    m.def("add_id", &add_id, py::arg("domainList"));
    m.def("create_relationsship_list", &create_relationsship_list, py::arg("domainList"),
          py::arg("pPercent") = 0.5);
    m.def("add_relation_to_list", &add_relation_to_list, py::arg("rList"), py::arg("chromosom"), py::arg("parent"),
          py::arg("child"));
    m.def("write_in_file", &write_in_file, py::arg("pList"), py::arg("name"), py::arg("relation") = false);
    m.def("create_tree", &create_tree, py::arg("rList"), py::arg("dList"), py::arg("pFileNameSuffix"),
          py::arg("pFormat"));
    m.def("create_small_list", &create_small_list, py::arg("dList"));
    m.def("read_protein", &read_protein, py::arg("file"));
    m.def("merge_protein", &merge_protein, py::arg("pProteinList"), py::arg("binSize"), py::arg("minPeak"));
    m.def("compare_boundaries_protein", &compare_boundaries_protein, py::arg("bList"), py::arg("cList"),
          py::arg("paraScore") = 0.2);
    m.def("create_list_with_protein", &create_list_with_protein, py::arg("bList"), py::arg("minPeak"),
          py::arg("cList") = py::none());
}

}  // namespace hicx::python
