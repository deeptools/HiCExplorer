// hicexplorer._hicx.hicDifferentialTAD: the C++ bodies of
// hicexplorer.hicDifferentialTAD.
//
// computeDifferentialTADs is the worker main() of 3.7.6 started once per
// share of a chromosome's TADs. Its blocks and tests are those of the
// hicDifferentialTAD tool (cpp/core/include/hicx/tad_contacts.hpp, which
// documents every quirk of the geometry it reproduces); this file reads the
// worker's arguments (a cool file name or a hiCMatrix object per matrix, the
// TAD rows, the share flag pThreadId) and puts into pQueue what the Python put
// there: the statistics, the p-values and the three accepted flags per TAD
// and the row objects, or 'Fail: <message>' and a traceback.
//
// Deviations from 3.7.6 (hicexplorer/test/api/cases_hicDifferentialTAD.py):
//  * the traceback after 'Fail: <message>' names the C++ worker instead of
//    Python frames;
//  * readDomainBoundaries reads the file with hicx::TextTable, which raises
//    pandas' ParserError also for a row with fewer fields, where pandas
//    fills them with NaN.

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicFindTADs.hpp"
#include "hicx/tad_contacts.hpp"

namespace hicx::python {

namespace tads = hicx::tads;

namespace {

constexpr const char* kLogger = "hicexplorer.hicDifferentialTAD";

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

void computeDifferentialTADs(py::handle pMatrixTarget, py::handle pMatrixControl, py::handle pDomainList,
                             py::handle pCoolOrH5, py::handle pPValue, py::handle pThreadId, py::handle pQueue) {
    py::list accepted_left;
    py::list accepted_right;
    py::list accepted_intra;
    py::list p_values_list;
    py::list stats_list;
    py::list rows;
    const auto fail = [&](const std::string& type, const std::string& message) {
        pQueue.attr("put")(to_str(hic::worker_failure("hicDifferentialTAD", type, message)));
    };
    try {
        const std::size_t count = length(pDomainList);
        const bool cool = PyObject_IsTrue(pCoolOrH5.ptr()) == 1;
        std::optional<tads::ContactMatrix> target;
        std::optional<tads::ContactMatrix> control;
        std::vector<tads::Domain> list;
        const double nan = std::numeric_limits<double>::quiet_NaN();
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
            if (!target.has_value()) {
                list = hic::domains(pDomainList);
                target = hic::contact_matrix(pMatrixTarget, cool);
                control = hic::contact_matrix(pMatrixControl, cool);
            }
            const tads::TadGeometry geometry_target = tads::tad_geometry(*target, list, i);
            const tads::TadGeometry geometry_control = tads::tad_geometry(*control, list, i);
            const auto test = [&](const tads::Block& a, const tads::Block& b) {
                const tads::DenseBlock values_target = tads::extract_block(*target, a);
                const tads::DenseBlock values_control = tads::extract_block(*control, b);
                return tads::rank_sum_test(values_target.values, values_control.values);
            };
            std::array<std::optional<tads::RankTest>, 3> tests;
            if (geometry_target.left.has_value()) {
                tests[0] = test(*geometry_target.left, *geometry_control.left);
            }
            if (geometry_target.right.has_value()) {
                tests[1] = test(*geometry_target.right, *geometry_control.right);
            }
            tests[2] = test(geometry_target.intra, geometry_control.intra);

            const double threshold = to_double(pPValue);
            py::list p_values;
            py::list stats;
            const std::array<py::list*, 3> accepted = {&accepted_left, &accepted_right, &accepted_intra};
            for (std::size_t t = 0; t < 3; ++t) {
                if (!tests[t].has_value() || std::isnan(tests[t]->pvalue)) {
                    accepted[t]->append(0);
                    p_values.append(nan);
                    stats.append(nan);
                } else {
                    accepted[t]->append(tests[t]->pvalue <= threshold ? 1 : 0);
                    p_values.append(numpy_scalar<double>(tests[t]->pvalue));
                    stats.append(numpy_scalar<double>(tests[t]->statistic));
                }
            }
            p_values_list.append(p_values);
            stats_list.append(stats);
            rows.append(row_at(pDomainList, i));
        }
    } catch (py::error_already_set& error) {
        const std::string type = str_of(error.type().attr("__name__"));
        const std::string message = str_of(error.value());
        fail(type, message);
        return;
    } catch (const hicx::PythonError& error) {
        fail(error.type(), error.what());
        return;
    } catch (const std::exception& error) {
        fail("ValueError", hic::worker_message(error));
        return;
    }
    py::list result;
    result.append(stats_list);
    result.append(p_values_list);
    result.append(accepted_left);
    result.append(accepted_right);
    result.append(accepted_intra);
    result.append(rows);
    pQueue.attr("put")(result);
}

}  // namespace

void bind_hicDifferentialTAD(py::module_& m) {
    m.def("readDomainBoundaries", &readDomainBoundaries, py::arg("pFile"));
    m.def("computeDifferentialTADs", &computeDifferentialTADs, py::arg("pMatrixTarget"), py::arg("pMatrixControl"),
          py::arg("pDomainList"), py::arg("pCoolOrH5"), py::arg("pPValue"), py::arg("pThreadId"),
          py::arg("pQueue"));
}

}  // namespace hicx::python
