// hicexplorer._hicx.list_tools: the C++ bodies of the public functions of
// hicexplorer.list_tools. The rules (which installed commands are tools, the
// one line summary, the aligned listing) are hicx/tool_list.hpp, shared with
// the `hicexplorer` executable. The package metadata and the tools'
// --help-json specifications are read through importlib.metadata and
// hicexplorer._cpp, which locates and runs the C++ executables.
//
// The functions call each other through the module (hicexplorer.list_tools),
// as the Python did, so a replaced tool_names or tool_list is honoured.

#include <string>
#include <vector>

#include "api.hpp"
#include "hicx/tool_list.hpp"

namespace hicx::python {

namespace tl = hicx::toollist;

namespace {

py::module_ module() {
    return py::module_::import("hicexplorer.list_tools");
}

py::list tool_names() {
    py::object distribution;
    try {
        distribution = py::module_::import("importlib.metadata").attr("distribution")("HiCExplorer");
    } catch (py::error_already_set& error) {
        if (error.matches(py::module_::import("importlib.metadata").attr("PackageNotFoundError").ptr())) {
            return py::list();
        }
        throw;
    }
    std::vector<std::string> scripts;
    for (py::handle entry : distribution.attr("entry_points")) {
        if (py::str(entry.attr("group")).equal(py::str("console_scripts"))) {
            scripts.push_back(utf8(entry.attr("name")));
        }
    }
    std::vector<tl::DistributionFile> files;
    const py::object listed = distribution.attr("files");
    if (!listed.is_none() && PyObject_IsTrue(listed.ptr()) == 1) {
        for (py::handle path : listed) {
            tl::DistributionFile file;
            for (py::handle part : path.attr("parts")) {
                file.parts.push_back(utf8(part));
            }
            file.name = utf8(path.attr("name"));
            file.stem = utf8(path.attr("stem"));
            file.suffix = utf8(path.attr("suffix"));
            files.push_back(std::move(file));
        }
    }
    py::list out;
    for (const std::string& name : tl::command_names(scripts, files)) {
        out.append(to_str(name));
    }
    return out;
}

std::string summary_of(const std::string& tool) {
    std::string description;
    try {
        const py::object spec = py::module_::import("hicexplorer._cpp").attr("tool_spec")(tool);
        const py::object value = spec.attr("get")("description");
        if (!value.is_none() && PyObject_IsTrue(value.ptr()) == 1) {
            description = utf8(value);
        }
    } catch (py::error_already_set& error) {
        if (error.matches(PyExc_OSError) || error.matches(PyExc_RuntimeError) ||
            error.matches(PyExc_ValueError)) {
            return "";
        }
        throw;
    }
    return tl::summary(description);
}

std::string tool_list() {
    const py::list names = module().attr("tool_names")();
    if (names.size() == 0) {
        return "The tool list is read from the package metadata; HiCExplorer is not installed.\n";
    }
    std::vector<std::pair<std::string, std::string>> entries;
    for (py::handle name : names) {
        const std::string tool = utf8(name);
        entries.emplace_back(tool, summary_of(tool));
    }
    return tl::format_tool_list(entries);
}

py::object parse_arguments(py::handle citation, py::handle version) {
    const py::module_ argparse = py::module_::import("argparse");
    const std::string listing = utf8(module().attr("tool_list")());
    const std::string description = utf8(citation) + "\nThe following is the list of tools:\n\n" +
                                    listing + "\nFor more information visit: http://hicexplorer.readthedocs.org\n";
    py::object parser = argparse.attr("ArgumentParser")(
        py::arg("formatter_class") = argparse.attr("RawDescriptionHelpFormatter"),
        py::arg("description") = to_str(description));
    parser.attr("add_argument")("--version", py::arg("action") = "version",
                                py::arg("version") = to_str("%(prog)s " + str_of(version)));
    return parser;
}

py::object process_args(py::handle args) {
    return module().attr("parse_arguments")().attr("parse_args")(args);
}

}  // namespace

void bind_list_tools(py::module_& m) {
    m.def("tool_names", &tool_names);
    m.def("tool_list", &tool_list);
    m.def("parse_arguments", &parse_arguments, py::arg("citation"), py::arg("version"));
    m.def("process_args", &process_args, py::arg("args") = py::none());
}

}  // namespace hicx::python
