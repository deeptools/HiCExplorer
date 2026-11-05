// hicexplorer._hicx.chicDifferentialTest: the C++ bodies of the functions of
// hicexplorer.chicDifferentialTest besides main() and parse_arguments()
// (which run the executable cpp/tools/chicDifferentialTest.cpp).
//
// The tests are Fisher's exact test and the chi-squared test of independence
// without correction on 2x2 tables, and the chi-squared quantile, as
// hicx::contingency computes them (cpp/core/include/hicx/contingency_tests.hpp),
// with the types the Python (scipy) returns the p-values in: a Python float for
// a constant 1.0, numpy.float64 for a computed value. The result file is written with
// hicx::chic::Hdf5Writer in the layout h5py gave it.
//
// Not reproduced: chi2_contingency on a table of other than two columns
// (3.7.6 tests the two counts of a position, [sum of interactions, raw], so
// no caller builds one) raises TypeError naming HiCExplorer 4.

#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "api.hpp"
#include "chic_common.hpp"
#include "hicx/contingency_tests.hpp"

namespace hicx::python {

using namespace chic_api;

namespace {

constexpr const char* kLogger = "hicexplorer.chicDifferentialTest";

// The values of a group as numbers; nullopt when numpy cannot build the
// table (it raises ValueError, which the tests catch).
std::vector<double> group_values(py::handle group) {
    std::vector<double> out;
    for (const py::handle value : group) {
        out.push_back(to_double(value));
    }
    return out;
}

// np.asarray(np.ceil(x), dtype=np.int64)
std::int64_t ceil_int64(double value) {
    const double c = std::ceil(value);
    if (!(c >= -9223372036854775808.0 && c < 9223372036854775808.0)) {
        return std::numeric_limits<std::int64_t>::min();
    }
    return static_cast<std::int64_t>(c);
}

struct TestResult {
    py::list test_result;
    py::list accepted;
    py::list rejected;
};

py::list pair(py::ssize_t i, py::handle value) {
    py::list out;
    out.append(py::int_(i));
    out.append(value);
    return out;
}

TestResult chisquare(py::handle pDataFile1, py::handle pDataFile2, py::handle pAlpha) {
    TestResult out;
    const double critical_value =
        hicx::contingency::chi2_quantile(to_double(sub(py::int_(1), pAlpha)), 1.0);
    std::int64_t zero_values_counter = 0;
    py::ssize_t i = 0;
    const py::object pairs = builtin("zip")(pDataFile1, pDataFile2);
    for (const py::handle item : pairs) {
        const std::vector<py::object> groups = unpack(item, 2);
        const std::vector<double> a = group_values(groups[0]);
        const std::vector<double> b = group_values(groups[1]);
        std::optional<hicx::contingency::Chi2Result> result;
        bool value_error = a.size() != b.size() || a.empty();
        if (!value_error) {
            if (a.size() != 2) {
                raise_error("TypeError",
                            "chisquare_test: HiCExplorer 4 tests tables of two columns, not " + std::to_string(a.size()));
            }
            result = hicx::contingency::chi2_2x2(a[0], a[1], b[0], b[1]);
            value_error = !result.has_value();
        }
        if (value_error) {
            ++zero_values_counter;
            out.test_result.append(py::float_(std::numeric_limits<double>::quiet_NaN()));
            out.accepted.append(pair(i, py::float_(1.0)));
        } else {
            const py::object p_value = np_float64(result->pvalue);
            out.test_result.append(p_value);
            if (result->statistic >= critical_value) {
                out.rejected.append(pair(i, p_value));
            } else {
                out.accepted.append(pair(i, p_value));
            }
        }
        ++i;
    }
    if (zero_values_counter > 0) {
        log(kLogger, "info", std::to_string(zero_values_counter) +
                                 " samples were not tested because at least one condition contained no data in "
                                 "both groups.");
    }
    return out;
}

TestResult fisher(py::handle pDataFile1, py::handle pDataFile2, py::handle pAlpha) {
    TestResult out;
    py::ssize_t i = 0;
    const py::object pairs = builtin("zip")(pDataFile1, pDataFile2);
    for (const py::handle item : pairs) {
        const std::vector<py::object> groups = unpack(item, 2);
        const std::vector<double> a = group_values(groups[0]);
        const std::vector<double> b = group_values(groups[1]);
        std::optional<hicx::contingency::FisherExact> result;
        if (a.size() == 2 && b.size() == 2) {
            result = hicx::contingency::fisher_exact_2x2(ceil_int64(a[0]), ceil_int64(a[1]), ceil_int64(b[0]), ceil_int64(b[1]));
        }
        if (!result) {
            out.test_result.append(py::float_(std::numeric_limits<double>::quiet_NaN()));
            out.accepted.append(pair(i, py::float_(1.0)));
        } else {
            const py::object p_value =
                result->constant_one ? py::object(py::float_(result->pvalue)) : np_float64(result->pvalue);
            if (le(p_value, pAlpha)) {
                out.test_result.append(p_value);
                out.rejected.append(pair(i, p_value));
            } else {
                out.test_result.append(p_value);
                out.accepted.append(pair(i, p_value));
            }
        }
        ++i;
    }
    return out;
}

py::tuple as_tuple(const TestResult& result) {
    return py::make_tuple(result.test_result, result.accepted, result.rejected);
}

}  // namespace

void bind_chicDifferentialTest(py::module_& m) {
    m.doc() = "C++ bodies of hicexplorer.chicDifferentialTest.";

    m.def("readInteractionFile", [](py::handle pInteractionFile) {
        const std::vector<std::string> lines = read_text_lines(pInteractionFile);
        std::size_t next = 0;
        const auto readline = [&]() -> py::object {
            return next < lines.size() ? to_str(lines[next++]) : py::str("");
        };
        (void)readline();
        py::object header = readline();
        // float(header.strip().split('\t')[-1].split(' ')[-1])
        const py::object last_field = getitem(header.attr("strip")().attr("split")("\t"), -1);
        const py::object sum_of_all_interactions = float_of(getitem(last_field.attr("split")(" "), -1));
        header = add(header, readline());
        py::list line_content;
        py::list data;
        for (; next < lines.size(); ++next) {
            const py::str line = to_str(lines[next]);
            if (truthy(line.attr("startswith")("#"))) {
                continue;
            }
            const py::object fields = line.attr("strip")().attr("split")("\t");
            if (length(fields) <= 1) {
                continue;
            }
            line_content.append(fields);
            py::list datum;
            datum.append(sum_of_all_interactions);
            datum.append(float_of(getitem(fields, -1)));
            data.append(datum);
        }
        return py::make_tuple(header, line_content, data);
    }, py::arg("pInteractionFile"));

    m.def("chisquare_test", [](py::handle pDataFile1, py::handle pDataFile2, py::handle pAlpha) {
        return as_tuple(chisquare(pDataFile1, pDataFile2, pAlpha));
    }, py::arg("pDataFile1"), py::arg("pDataFile2"), py::arg("pAlpha"));

    m.def("fisher_exact_test", [](py::handle pDataFile1, py::handle pDataFile2, py::handle pAlpha) {
        return as_tuple(fisher(pDataFile1, pDataFile2, pAlpha));
    }, py::arg("pDataFile1"), py::arg("pDataFile2"), py::arg("pAlpha"));

    m.def("writeResult", [](py::handle pOutFileName, py::handle pData, py::handle pHeaderOld, py::handle pHeaderNew,
                            py::handle pAlpha, py::handle pTest) {
        with_open(pOutFileName, "w", [&](py::object& file) {
            std::string header = "# Differential analysis result file of HiCExplorer's chicDifferentialTest version ";
            header += utf8(str_object(package_version()));
            header += "\n";
            header += "# This file contains the p-values computed by " + utf8(format(pTest, "")) + " test\n";
            header += "# To test the smoothed (float) values they were rounded up to the next integer\n";
            header += "#\n";
            header += "# Alpha level " + utf8(str_object(pAlpha));
            header += "\n";
            header += "# Degrees of freedom 1";
            header += "\n#\n";
            write(file, header);
            file.attr("write")(add(getitem(pHeaderOld.attr("split")("\n"), 0), py::str("\n")));
            file.attr("write")(add(getitem(pHeaderNew.attr("split")("\n"), 0), py::str("\n")));
            write(file,
                  "#Chromosome\tStart\tEnd\tGene\tRelative distance\tsum of interactions 1\ttarget_1 raw\tsum of "
                  "interactions 2\ttarget_2 raw\tp-value\n");
            if (truthy(pData)) {
                for (const py::handle data : pData) {
                    const py::object first = getitem(data, 0);
                    py::object line = join("\t", getitem(first, py::slice(py::none(), py::int_(4), py::none())));
                    std::string text = utf8(line) + "\t" + utf8(format(getitem(first, 5), "")) + "\t";
                    const auto five = [](py::handle values) {
                        py::list parts;
                        for (const py::handle x : values) {
                            parts.append(format(x, ".5f"));
                        }
                        return utf8(join("\t", parts));
                    };
                    text += five(getitem(data, 3)) + "\t";
                    text += five(getitem(data, 4)) + "\t";
                    text += "\t" + utf8(format(format(getitem(data, 2), ".5f"), "")) + "\n";
                    write(file, text);
                }
            }
        });
    }, py::arg("pOutFileName"), py::arg("pData"), py::arg("pHeaderOld"), py::arg("pHeaderNew"), py::arg("pAlpha"),
       py::arg("pTest"));

    m.def("writeResultHDF", [](py::handle pOutFileName, py::handle pAcceptedData, py::handle pRejectedData,
                               py::handle pAllResultData, py::handle pInputData, py::handle pAlpha, py::handle pTest) {
        hicx::chic::Hdf5Writer writer = create_h5(pOutFileName);
        write_attribute(writer, "/", "type", py::str("differential"));
        write_attribute(writer, "/", "version", package_version());
        write_attribute(writer, "/", "alpha", pAlpha);
        write_attribute(writer, "/", "test", pTest);
        const char* const categories[] = {"accepted", "rejected", "all"};
        const py::handle category_data[] = {pAcceptedData, pRejectedData, pAllResultData};
        py::ssize_t i = 0;
        for (const py::handle input : pInputData) {
            const std::string matrix1 = utf8(getitem(getitem(input, 0), 1));
            const std::string matrix2 = utf8(getitem(getitem(input, 1), 1));
            const std::string chromosome = utf8(getitem(getitem(input, 0), 2));
            const py::object gene_name = getitem(getitem(input, 0), 3);
            const std::string gene = utf8(gene_name);
            ensure_group(writer, h5_path({matrix1}));
            ensure_group(writer, h5_path({matrix1, matrix2}));
            ensure_group(writer, h5_path({matrix1, matrix2, "genes"}));
            ensure_group(writer, h5_path({matrix1, matrix2, chromosome}));
            const std::string gene_path = h5_path({matrix1, matrix2, chromosome, gene});
            create_group(writer, gene_path);
            create_group(writer, gene_path + "/accepted");
            create_group(writer, gene_path + "/rejected");
            create_group(writer, gene_path + "/all");
            for (int c = 0; c < 3; ++c) {
                const std::string base = gene_path + "/" + categories[c];
                const py::object data_object = getitem(category_data[c], i);
                if (length(data_object) == 0) {
                    continue;
                }
                py::object row_chromosome = py::none();
                py::list start_list;
                py::list end_list;
                py::object sum_1 = py::none();
                py::object sum_2 = py::none();
                py::list relative_distance_list;
                py::list pvalue_list;
                py::list raw_1;
                py::list raw_2;
                for (const py::handle data : data_object) {
                    const py::object line = getitem(data, 0);
                    row_chromosome = getitem(line, 0);
                    start_list.append(getitem(line, 1));
                    end_list.append(getitem(line, 2));
                    relative_distance_list.append(getitem(line, 5));
                    sum_1 = float_of(getitem(getitem(data, 3), 0));
                    sum_2 = float_of(getitem(getitem(data, 4), 0));
                    raw_1.append(getitem(getitem(data, 3), 1));
                    raw_2.append(getitem(getitem(data, 4), 1));
                    pvalue_list.append(getitem(data, 2));
                }
                write_scalar(writer, base + "/chromosome", str_object(row_chromosome));
                write_values(writer, base + "/start_list", start_list, 9);
                write_values(writer, base + "/end_list", end_list, 9);
                write_scalar(writer, base + "/gene", str_object(gene_name));
                write_values(writer, base + "/relative_distance_list", relative_distance_list, 9);
                write_scalar(writer, base + "/sum_of_interactions_1", float_of(sum_1));
                write_scalar(writer, base + "/sum_of_interactions_2", float_of(sum_2));
                write_values(writer, base + "/raw_target_list_1", raw_1, 9);
                write_values(writer, base + "/raw_target_list_2", raw_2, 9);
                write_values(writer, base + "/pvalue_list", pvalue_list, 9);
            }
            if (!writer.hard_link(gene_path, h5_path({matrix1, matrix2, "genes", gene}))) {
                log(kLogger, "debug", "exception Unable to create link (name already exists)");
            }
            ++i;
        }
        writer.close();
    }, py::arg("pOutFileName"), py::arg("pAcceptedData"), py::arg("pRejectedData"), py::arg("pAllResultData"),
       py::arg("pInputData"), py::arg("pAlpha"), py::arg("pTest"));

    m.def("run_statistical_tests", [](py::handle pInteractionFilesList, py::handle pArgs, py::handle pViewpointObject,
                                      py::handle pQueue) -> py::object {
        py::list accepted_list;
        py::list rejected_list;
        py::list all_results_list;
        try {
            for (const py::handle interaction_file : pInteractionFilesList) {
                const py::object aggregated = pArgs.attr("aggregatedFile");
                const std::vector<py::object> first =
                    unpack(pViewpointObject.attr("readAggregatedFileHDF")(aggregated, getitem(interaction_file, 0)), 2);
                const std::vector<py::object> second =
                    unpack(pViewpointObject.attr("readAggregatedFileHDF")(aggregated, getitem(interaction_file, 1)), 2);
                const py::object& line_content1 = first[0];
                const py::object& data1 = first[1];
                const py::object& line_content2 = second[0];
                const py::object& data2 = second[1];
                if (length(line_content1) == 0 || length(line_content2) == 0) {
                    continue;
                }
                const py::object test = pArgs.attr("statisticTest");
                TestResult result;
                if (eq(test, py::str("chi2"))) {
                    result = chisquare(data1, data2, pArgs.attr("alpha"));
                } else if (eq(test, py::str("fisher"))) {
                    result = fisher(data1, data2, pArgs.attr("alpha"));
                } else {
                    raise_error("UnboundLocalError",
                                "cannot access local variable 'test_result' where it is not associated with a value");
                }
                const auto line = [&](py::handle index, py::handle value) {
                    py::list out;
                    out.append(getitem(line_content1, index));
                    out.append(getitem(line_content2, index));
                    out.append(value);
                    out.append(getitem(data1, index));
                    out.append(getitem(data2, index));
                    return out;
                };
                py::list write_out_lines;
                py::ssize_t i = 0;
                for (const py::handle value : result.test_result) {
                    write_out_lines.append(line(py::int_(i), value));
                    ++i;
                }
                py::list write_out_lines_accepted;
                for (const py::handle entry : result.accepted) {
                    write_out_lines_accepted.append(line(getitem(entry, 0), getitem(entry, 1)));
                }
                py::list write_out_lines_rejected;
                for (const py::handle entry : result.rejected) {
                    write_out_lines_rejected.append(line(getitem(entry, 0), getitem(entry, 1)));
                }
                accepted_list.append(write_out_lines_accepted);
                rejected_list.append(write_out_lines_rejected);
                all_results_list.append(write_out_lines);
            }
        } catch (py::error_already_set& error) {
            if (!is_exception(error)) {
                throw;
            }
            pQueue.attr("put")(to_str("Fail: " + str_of(error.value())));
            return py::none();
        } catch (const std::exception& error) {
            pQueue.attr("put")(to_str("Fail: " + fail_text(error)));
            return py::none();
        }
        if (pQueue.is_none()) {
            return py::none();
        }
        py::list result;
        result.append(accepted_list);
        result.append(rejected_list);
        result.append(all_results_list);
        pQueue.attr("put")(result);
        return py::none();
    }, py::arg("pInteractionFilesList"), py::arg("pArgs"), py::arg("pViewpointObject"), py::arg("pQueue") = py::none());
}

}  // namespace hicx::python
