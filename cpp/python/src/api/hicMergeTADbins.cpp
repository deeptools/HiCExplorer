// hicexplorer._hicx.hicMergeTADbins: the C++ bodies of
// hicexplorer.hicMergeTADbins.
//
// Both functions take the hicmatrix.HiCMatrix.hiCMatrix object of 3.7.6.
// get_boundary_bin_id reads the domains from the file handle and looks their
// first and last bins up in the object's bin table; merge_tad_bins restores
// the masked bins, merges the bins of every domain with hicx::plan_tad_merge
// and hicx::reduce_matrix (the code of the hicMergeTADbins tool), performs
// hiCMatrix.update_matrix on the object and saves it through the C++ writers
// with the settings of the object's file handler, as hiCMatrix.save does.
//
// Deviations from 3.7.6 (hicexplorer/test/api/cases_hicMergeTADbins.py): the
// file is written by the C++ cool and h5 writers of the tools (cool_adapter,
// h5_file), whose output the CLI harness compares with hicmatrix's, so the
// files hold the same matrix, bins and attributes but are not byte identical
// (creation date, HDF5 chunking).

#include <algorithm>
#include <cstdint>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicFindTADs.hpp"
#include "hicx/adjust_ops.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/py_text.hpp"
#include "hicx/reduce_matrix.hpp"
#include "utilities.hpp"

namespace hicx::python {

namespace {

constexpr const char* kLogger = "hicexplorer.hicMergeTADbins";

bool truthy(py::handle object) {
    const int result = PyObject_IsTrue(object.ptr());
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() && text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

[[noreturn]] void sys_exit() {
    py::module_::import("sys").attr("exit")();
    raise_error("SystemExit", "None");
}

// hiCMatrix.save(filename) with its defaults pSymmetric=True,
// pApplyCorrection=False and pHiCInfo=None.
void save(py::handle hic, py::handle filename) {
    const std::string name = utf8(py::str(filename));
    if (hic.attr("matrixFileHandler").is_none()) {
        const std::string type = ends_with(name, "h5") ? "h5" : "cool";
        hic.attr("matrixFileHandler") = py::module_::import("hicmatrix.lib").attr("MatrixFileHandler")(
            py::arg("pFileType") = type, py::arg("pHiCInfo") = py::none());
    }
    hic::restore_masked_bins(hic);
    if (!ends_with(name, "cool") && !ends_with(name, "h5")) {
        return;  // hiCMatrix.save writes nothing for any other name
    }
    hicx::MatrixData data = matrix_data(hic);
    const py::object file = hic.attr("matrixFileHandler").attr("matrixFile");
    const std::string kind = str_of(py::type::of(file).attr("__name__"));
    if (kind == "H5") {
        hicx::H5SaveOptions options;
        options.symmetric = true;
        hicx::write_hicexplorer_h5(name, data, options);
        return;
    }
    if (kind != "Cool") {
        raise_error("TypeError", "HiCExplorer 4 saves cool and h5 matrices, not " + kind);
    }
    hicx::CoolSaveOptions options;
    options.symmetric = true;
    options.apply_correction = false;
    options.enforce_integer = truthy(file.attr("enforceInteger"));
    options.append = truthy(file.attr("appendData"));
    options.file_was_h5 = truthy(file.attr("fileWasH5"));
    const py::object hic2cool = file.attr("hic2cool_version");
    if (!hic2cool.is_none()) {
        options.hic2cool_version = str_of(hic2cool);
    }
    const py::object op = file.attr("correctionOperator");
    if (!op.is_none()) {
        options.correction_operator = str_of(op)[0];
    }
    if (ends_with(name, "cool")) {
        file.attr("hic_metadata") = py::none();  // pHiCInfo
    }
    const py::object metadata = file.attr("hic_metadata");
    if (!metadata.is_none()) {
        options.has_hic_metadata = true;
        for (const char* key : {"genome-assembly", "matrix-generated-by", "matrix-generated-by-url"}) {
            if (metadata.contains(key)) {
                options.hic_metadata[key] = str_of(metadata[key]);
                PyObject_DelItemString(metadata.ptr(), key);
            }
        }
    }
    hicx::write_cool(name, data, options);
    file.attr("appendData") = py::str(options.append ? "a" : "w");
}

void merge_tad_bins(py::handle hic, py::handle boundary_id_list, py::handle filename) {
    hic::restore_masked_bins(hic);
    const py::object intervals = hic.attr("cut_intervals");
    const std::vector<hicx::CutInterval> table = cut_intervals(intervals);
    std::set<std::int64_t> ids;
    for (const py::handle id : py::reinterpret_borrow<py::iterable>(boundary_id_list)) {
        ids.insert(to_int64(id));
    }
    const hicx::BinMergePlan plan =
        hicx::plan_tad_merge(table, std::vector<std::int64_t>(ids.begin(), ids.end()));
    if (plan.bins_to_merge.empty()) {
        log(kLogger, "info", "Nothing to merge.");
        return;
    }
    py::list new_bins;
    for (std::size_t k = 0; k < plan.bins_to_merge.size(); ++k) {
        const std::vector<std::int64_t>& group = plan.bins_to_merge[k];
        for (const std::int64_t b : group) {
            if (!table[static_cast<std::size_t>(b)].extra_text.empty()) {
                raise_error("TypeError", "cannot perform reduce with flexible type");
            }
        }
        const py::object first = py::reinterpret_borrow<py::object>(intervals)[py::int_(group.front())];
        const py::object last = py::reinterpret_borrow<py::object>(intervals)[py::int_(group.back())];
        new_bins.append(py::make_tuple(first[py::int_(0)], first[py::int_(1)], last[py::int_(2)],
                                       numpy_scalar<double>(plan.intervals[k].extra)));
    }
    hic.attr("correction_factors") = py::none();
    // hiCMatrix.update_matrix(reduce_matrix(matrix, bins_to_merge, diagonal=True), new_bins)
    hicx::CsrMatrix merged = hicx::reduce_matrix(hic::matrix(hic), plan.bins_to_merge, true, true);
    if (PyObject_Length(hic.attr("orig_bin_ids").ptr()) > 0) {
        raise_error("ValueError", "matrix contains masked bins. Restore masked bins first");
    }
    const std::vector<std::int64_t> nan_bins = hicx::empty_column_bins(merged);
    hic.attr("matrix") = hic::csr_object(std::move(merged));
    hic::set_intervals(hic, new_bins);
    hic.attr("nan_bins") = from_vector(nan_bins);
    save(hic, filename);
}

py::object get_boundary_bin_id(py::handle hic, py::handle bed_fh) {
    const hicx::BinTable table = hic::bins(hic);
    std::set<std::int64_t> boundaries;
    std::int64_t line_number = 0;
    for (const py::handle raw : py::reinterpret_borrow<py::iterable>(bed_fh.attr("readlines")())) {
        ++line_number;
        const py::object line = to_string(raw);
        const std::string text = is_str(line) ? utf8(line) : str_of(line);
        if (text.rfind("browser", 0) == 0 || text.rfind("track", 0) == 0 || text.rfind("#", 0) == 0) {
            continue;
        }
        const std::string stripped(hicx::pytext::strip(text));
        const std::vector<std::string_view> fields = hicx::pytext::split(stripped, "\t");
        if (fields.size() < 3) {
            const std::string detail =
                "not enough values to unpack (expected 3, got " + std::to_string(fields.size()) + ")";
            log(kLogger, "error", "Could not read line\n" + text + "\n. " + detail);
            sys_exit();
        }
        const py::str chrom = to_str(std::string(fields[0]));
        py::object start;
        py::object end;
        try {
            start = parse_int(to_str(std::string(fields[1])));
            end = parse_int(to_str(std::string(fields[2])));
        } catch (py::error_already_set& error) {
            if (!error.matches(PyExc_ValueError)) {
                throw;
            }
            log(kLogger, "error",
                "Error reading line: " + std::to_string(line_number) +
                    ". One of the fields is not an integer.\nError message: " + str_of(error.value()));
            sys_exit();
        }
        if (PyObject_RichCompareBool(start.ptr(), end.ptr(), Py_LE) != 1 &&
            py::module_::import("sys").attr("flags").attr("optimize").cast<int>() == 0) {
            raise_error("AssertionError", "Error in line #" + std::to_string(line_number) +
                                              ", end1 larger than start1 in " + text);
        }
        const auto range = hic::region_bin_range(hic, table, chrom, start, end);
        if (!range.has_value()) {
            raise_error("TypeError", "cannot unpack non-iterable NoneType object");
        }
        boundaries.insert(range->first);
        boundaries.insert(range->second);
    }
    if (boundaries.empty()) {
        return from_vector(std::vector<double>{});  // np.sort([]) is float64
    }
    return from_vector(std::vector<std::int64_t>(boundaries.begin(), boundaries.end()));
}

}  // namespace

void bind_hicMergeTADbins(py::module_& m) {
    m.def("merge_tad_bins", &merge_tad_bins, py::arg("hic"), py::arg("boundary_id_list"), py::arg("filename"));
    m.def("get_boundary_bin_id", &get_boundary_bin_id, py::arg("hic"), py::arg("bed_fh"));
}

}  // namespace hicx::python
