// hicexplorer._hicx.hicInterIntraTAD: the C++ bodies of
// hicexplorer.hicInterIntraTAD.
//
// computeInterIntraTADs is the worker main() of 3.7.6 started once per share
// of a chromosome's TADs. Its blocks and sums are those of the
// hicInterIntraTAD tool (cpp/core/include/hicx/tad_contacts.hpp); this file
// reads the worker's arguments and puts into pQueue the fifteen lists and the
// row objects the Python put there, with the Python types of every value: the
// literal 0 where a neighbour is missing, numpy scalars of the matrix dtype
// for the sums and their ratios, Python ints and floats for the counts and
// densities. A failure puts 'Fail: <message>' and a traceback.
//
// Deviations from 3.7.6 (hicexplorer/test/api/cases_hicInterIntraTAD.py):
//  * the traceback after 'Fail: <message>' names the C++ worker instead of
//    Python frames;
//  * readDomainBoundaries reads the file with hicx::TextTable, which raises
//    pandas' ParserError also for a row with fewer fields, where pandas
//    fills them with NaN.

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicFindTADs.hpp"
#include "hicx/tad_contacts.hpp"

namespace hicx::python {

namespace tads = hicx::tads;

namespace {

constexpr const char* kLogger = "hicexplorer.hicInterIntraTAD";

py::object readDomainBoundaries(py::handle pFile) {
    return hic::read_domain_boundaries(pFile);
}

std::size_t length(py::handle object) {
    const Py_ssize_t n = PyObject_Length(object.ptr());
    if (n < 0) {
        throw_python_error();
    }
    return static_cast<std::size_t>(n);
}

py::object row_at(py::handle rows, std::size_t index) {
    PyObject* row = PySequence_GetItem(rows.ptr(), static_cast<Py_ssize_t>(index));
    if (row == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(row);
}

// A value as the Python holds it.
py::object number(const tads::PyNumber& n) {
    switch (n.kind) {
        case tads::PyNumber::Kind::PyInt: return py::int_(n.integer);
        case tads::PyNumber::Kind::Int64: return numpy_scalar<std::int64_t>(n.integer);
        case tads::PyNumber::Kind::Float32: return numpy_scalar<float>(static_cast<float>(n.real));
        case tads::PyNumber::Kind::Float64: return numpy_scalar<double>(n.real);
    }
    return py::none();
}

struct Side {
    tads::PyNumber sum = tads::PyNumber::py_int(0);
    py::object density = py::int_(0);
    std::int64_t contacts = 0;
    std::int64_t nnz = 0;
};

// sum, shape[0] * shape[1], nnz and nnz / (shape[0] * shape[1]) of one block.
Side measure(const tads::ContactMatrix& matrix, const tads::Block& block) {
    const tads::DenseBlock values = tads::extract_block(matrix, block);
    Side side;
    side.sum = tads::block_sum(values, tads::block_dtype(matrix, block));
    side.contacts = values.rows * values.cols;
    side.nnz = values.nnz;
    if (side.contacts == 0) {
        raise_error("ZeroDivisionError", "division by zero");
    }
    side.density = py::float_(static_cast<double>(side.nnz) / static_cast<double>(side.contacts));
    return side;
}

void computeInterIntraTADs(py::handle pMatrix, py::handle pDomainList, py::handle pCoolOrH5, py::handle pThreadId,
                           py::handle pQueue) {
    constexpr std::size_t kLists = 15;
    std::vector<py::list> lists(kLists);
    py::list rows;
    const auto fail = [&](const std::string& type, const std::string& message) {
        pQueue.attr("put")(to_str(hic::worker_failure("hicInterIntraTAD", type, message)));
    };
    try {
        const std::size_t count = length(pDomainList);
        const bool cool = PyObject_IsTrue(pCoolOrH5.ptr()) == 1;
        std::optional<tads::ContactMatrix> matrix;
        std::vector<tads::Domain> list;
        for (std::size_t i = 0; i < count; ++i) {
            if (pThreadId.is_none()) {
                log(kLogger, "debug", "first thread");
            } else if (PyObject_RichCompareBool(pThreadId.ptr(), Py_True, Py_EQ) == 1) {
                log(kLogger, "debug", "middle thread");
                if (i == 0 || i + 1 == count) {
                    log(kLogger, "debug", "i: " + std::to_string(i));
                    log(kLogger, "debug", "len(chromosome_list): " + std::to_string(count));
                }
            } else if (PyObject_RichCompareBool(pThreadId.ptr(), Py_False, Py_EQ) == 1) {
                log(kLogger, "debug", "last thread");
            }
            if (hic::skips(pThreadId, i, count)) {
                continue;
            }
            if (!matrix.has_value()) {
                list = hic::domains(pDomainList);
                matrix = hic::contact_matrix(pMatrix, cool);
            }
            const tads::TadGeometry geometry = tads::tad_geometry(*matrix, list, i);
            const Side intra = measure(*matrix, geometry.intra);
            Side left;
            Side right;
            if (geometry.left.has_value()) {
                left = measure(*matrix, *geometry.left);
            }
            if (geometry.right.has_value()) {
                right = measure(*matrix, *geometry.right);
            }
            const tads::PyNumber left_ratio = tads::py_divide(left.sum, intra.sum);
            const tads::PyNumber right_ratio = tads::py_divide(right.sum, intra.sum);
            const tads::PyNumber both_ratio = tads::py_divide(tads::py_add(left.sum, right.sum), intra.sum);
            const py::object values[kLists] = {
                number(left.sum),     number(right.sum),       left.density,
                right.density,        py::int_(left.contacts), py::int_(right.contacts),
                py::int_(left.nnz),   py::int_(right.nnz),     number(intra.sum),
                py::int_(intra.contacts), py::int_(intra.nnz), intra.density,
                number(left_ratio),   number(right_ratio),     number(both_ratio)};
            for (std::size_t k = 0; k < kLists; ++k) {
                lists[k].append(values[k]);
            }
            rows.append(row_at(pDomainList, i));
        }
    } catch (py::error_already_set& error) {
        fail(str_of(error.type().attr("__name__")), str_of(error.value()));
        return;
    } catch (const hicx::PythonError& error) {
        fail(error.type(), error.what());
        return;
    } catch (const std::exception& error) {
        fail("ValueError", hic::worker_message(error));
        return;
    }
    py::list result;
    for (const py::list& values : lists) {
        result.append(values);
    }
    result.append(rows);
    pQueue.attr("put")(result);
}

}  // namespace

void bind_hicInterIntraTAD(py::module_& m) {
    m.def("readDomainBoundaries", &readDomainBoundaries, py::arg("pFile"));
    m.def("computeInterIntraTADs", &computeInterIntraTADs, py::arg("pMatrix"), py::arg("pDomainList"),
          py::arg("pCoolOrH5"), py::arg("pThreadId"), py::arg("pQueue"));
}

}  // namespace hicx::python
