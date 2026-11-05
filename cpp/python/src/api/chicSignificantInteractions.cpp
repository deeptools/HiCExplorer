// hicexplorer._hicx.chicSignificantInteractions: the C++ bodies of the
// functions of hicexplorer.chicSignificantInteractions besides main() and
// parse_arguments() (which run the executable
// cpp/tools/chicSignificantInteractions.cpp).
//
// The p-value of a merged candidate is 1 - betainc(size, raw + 1, prob)
// (hicx::stats::betainc, what hicexplorer.lib.cnb.cdf evaluates), the target
// regions are sorted and merged as pybedtools' sort().merge(d=resolution)
// does (hicx::chic::sort_merge_regions), the unique sample keys are
// np.unique(axis=0, return_index=True) with its stable lexicographic order,
// and the two output files are written with hicx::chic::Hdf5Writer in the
// layout h5py gave them. The interaction data are read and merged with the
// Viewpoint object's methods, as in 3.7.6.
//
// Differences to 3.7.6:
//  * call_multi_core ran compute_interaction_file in one process per thread
//    and polled the queues once per second; here it runs in this process,
//    once per thread, and the results are collected in thread order. A
//    worker that raises instead of putting its result made 3.7.6 wait
//    forever; here its exception propagates.
//  * writeSignificantHDF and writeTargetHDF write each viewpoint's group
//    themselves, in the layout of Viewpoint.writeInteractionFileHDF5 and
//    Viewpoint.createUniqueHDFGroup, instead of calling those methods with
//    an h5py group (the files are written by C++, not h5py). Where
//    createUniqueHDFGroup would loop forever (a group name that exists
//    already), a ValueError is raised, as the tool stops.
//  * A failure's traceback text (put into the queue after the message) names
//    the frames the error passed in this implementation.

#include <algorithm>
#include <cstdint>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>

#include "api.hpp"
#include "chic_common.hpp"
#include "hicx/chic_hdf5.hpp"
#include "hicx/stats_ops.hpp"

namespace hicx::python {

using namespace chic_api;
namespace hc = hicx::chic;

namespace {

constexpr const char* kLogger = "hicexplorer.chicSignificantInteractions";

bool is_float(py::handle value) { return PyFloat_Check(value.ptr()) != 0; }
bool is_dict(py::handle value) { return PyDict_Check(value.ptr()) != 0; }

// filter_by_pvalue
py::tuple filter_by_pvalue(py::handle pData, py::handle pPValue, py::handle pMergedLinesDict,
                           py::handle pPeakInteractionsThreshold) {
    py::dict accepted;
    py::list accepted_lines;
    const bool single = is_float(pPValue);
    if (single || is_dict(pPValue)) {
        for (const py::handle key : pData) {
            const py::object line = getitem(pData, key);
            const py::object limit = single ? py::reinterpret_borrow<py::object>(pPValue) : getitem(pPValue, key);
            if (le(getitem(line, -3), limit)) {
                if (ge(float_of(getitem(line, -1)), pPeakInteractionsThreshold)) {
                    const py::object merged = getitem(pMergedLinesDict, key);
                    accepted[key] = merged;
                    accepted_lines.append(getitem(merged, py::slice(py::none(), py::int_(3), py::none())));
                }
            }
        }
    }
    return py::make_tuple(accepted, accepted_lines);
}

// compute_new_p_values
py::tuple compute_new_p_values(py::handle pData, py::handle pBackgroundModel, py::handle pPValue,
                               py::handle pMergedLinesDict, py::handle pPeakInteractionsThreshold) {
    py::dict accepted;
    py::list accepted_lines;
    const bool single = is_float(pPValue);
    if (single || is_dict(pPValue)) {
        for (const py::handle key : pData) {
            if (!contains(pBackgroundModel, key)) {
                log(kLogger, "debug", "key not in background " + utf8(format(key, "")));
                continue;
            }
            const py::object line = getitem(pData, key);
            const double raw = to_double(float_of(getitem(line, -1)));
            const py::object model = getitem(pBackgroundModel, key);
            const double size = to_double(float_of(getitem(model, 0)));
            const double probability = to_double(float_of(getitem(model, 1)));
            setitem(line, -3, np_float64(1.0 - hicx::stats::betainc(size, raw + 1.0, probability)));
            const py::object limit = single ? py::reinterpret_borrow<py::object>(pPValue) : getitem(pPValue, key);
            if (le(getitem(line, -3), limit)) {
                if (ge(float_of(getitem(line, -1)), pPeakInteractionsThreshold)) {
                    accepted[key] = line;
                    const py::object merged = getitem(pMergedLinesDict, key);
                    const py::object target_content =
                        getitem(getitem(merged, 0), py::slice(py::none(), py::int_(3), py::none()));
                    setitem(target_content, 2, getitem(getitem(merged, -1), 2));
                    accepted_lines.append(target_content);
                }
            }
        }
    }
    return py::make_tuple(accepted, accepted_lines);
}

// merge_neighbors_x_fold and merge_neighbors_loose_p_value
py::object merge_selected(const py::dict& accepted_line, py::handle pViewpointObj, py::handle pResolution) {
    if (truthy(accepted_line)) {
        return pViewpointObj.attr("merge_neighbors")(accepted_line, py::arg("pMergeThreshold") = pResolution);
    }
    return py::make_tuple(accepted_line, py::none());
}

py::object merge_x_fold(py::handle pXfold, py::handle pData, py::handle pViewpointObj, py::handle pResolution) {
    py::dict accepted_line;
    const bool single = is_float(pXfold);
    if (single || is_dict(pXfold)) {
        const py::object values = getitem(pData, 0);
        const py::object lines = getitem(pData, 1);
        for (const py::handle key : values) {
            const py::object limit = single ? py::reinterpret_borrow<py::object>(pXfold) : getitem(pXfold, key);
            if (lt(getitem(getitem(values, key), -1), limit)) {
                continue;
            }
            accepted_line[key] = getitem(lines, key);
        }
    }
    return merge_selected(accepted_line, pViewpointObj, pResolution);
}

py::object merge_loose(py::handle pLoosePValue, py::handle pData, py::handle pViewpointObj, py::handle pResolution,
                       py::handle pTruncateZeroPvalues) {
    py::dict accepted_line;
    const bool single = is_float(pLoosePValue);
    if (single || is_dict(pLoosePValue)) {
        const py::object values = getitem(pData, 0);
        const py::object lines = getitem(pData, 1);
        for (const py::handle key : values) {
            const py::object limit =
                single ? py::reinterpret_borrow<py::object>(pLoosePValue) : getitem(pLoosePValue, key);
            const py::object pvalue = getitem(getitem(values, key), 1);
            if (truthy(pTruncateZeroPvalues)) {
                if (eq(pvalue, py::int_(0)) || gt(pvalue, limit)) {
                    continue;
                }
            } else if (gt(pvalue, limit)) {
                continue;
            }
            accepted_line[key] = getitem(lines, key);
        }
    }
    return merge_selected(accepted_line, pViewpointObj, pResolution);
}

py::object compute_interaction_file(py::handle pInteractionFilesList, py::handle pArgs, py::handle pViewpointObj,
                                    py::handle pBackground, py::handle pFilePath, py::handle pResolution,
                                    py::handle pQueue) {
    py::list significant_data_list;
    py::list significant_key_list;
    py::list target_data_list;
    py::list target_key_list;
    py::list reference_points_list_target;
    py::list reference_points_list_significant;
    const auto failure = [&](const std::string& message, const std::exception& error) -> py::object {
        const std::string text = "Fail: " + message + traceback_text(error);
        if (pQueue.is_none()) {
            return to_str(text);
        }
        pQueue.attr("put")(to_str(text));
        return py::none();
    };
    try {
        for (const py::handle interaction_file : pInteractionFilesList) {
            py::list target_list;
            py::list sample_prefix;
            py::object data;
            for (const py::handle sample : interaction_file) {
                data = pViewpointObj.attr("readInteractionFile")(pFilePath, sample);
                sample_prefix.append(getitem(sample, 0));
                bool compute_new = false;
                py::object accepted_scores;
                py::object merged_lines;
                py::object target_lines;
                const py::object x_fold = pArgs.attr("xFoldBackground");
                const py::object loose = pArgs.attr("loosePValue");
                if (!x_fold.is_none()) {
                    const std::vector<py::object> merged =
                        unpack(merge_x_fold(x_fold, data, pViewpointObj, pResolution), 2);
                    accepted_scores = merged[0];
                    merged_lines = merged[1];
                    compute_new = true;
                } else if (!loose.is_none()) {
                    const std::vector<py::object> merged = unpack(
                        merge_loose(loose, data, pViewpointObj, pResolution, pArgs.attr("truncateZeroPvalues")), 2);
                    accepted_scores = merged[0];
                    merged_lines = merged[1];
                    compute_new = true;
                }
                py::tuple selected;
                if (compute_new) {
                    selected = compute_new_p_values(accepted_scores, pBackground, pArgs.attr("pValue"), merged_lines,
                                                    pArgs.attr("peakInteractionsThreshold"));
                } else {
                    selected = filter_by_pvalue(getitem(data, 0), pArgs.attr("pValue"), getitem(data, 1),
                                                pArgs.attr("peakInteractionsThreshold"));
                }
                accepted_scores = selected[0];
                target_lines = selected[1];
                if (length(accepted_scores) == 0) {
                    with_open(py::str("errorLog.txt"), "a+", [&](py::object& errorlog) {
                        write(errorlog, "Failed for: " + utf8(format(interaction_file, "")) + ".\n");
                    });
                }
                target_list.append(target_lines);
                significant_data_list.append(accepted_scores);
                significant_key_list.append(sample);
                reference_points_list_significant.append(getitem(data, 2));
            }
            py::list flat;
            for (const py::handle sublist : target_list) {
                for (const py::handle item : sublist) {
                    flat.append(item);
                }
            }
            sample_prefix.append(py::str("::"));
            sample_prefix.append(getitem(getitem(interaction_file, 0), 1));
            sample_prefix.append(getitem(getitem(interaction_file, 0), 2));
            target_data_list.append(flat);
            target_key_list.append(sample_prefix);
            if (!data) {
                raise_error("UnboundLocalError",
                            "cannot access local variable 'data' where it is not associated with a value");
            }
            reference_points_list_target.append(getitem(data, 2));
        }
        if (pQueue.is_none()) {
            return py::list();
        }
    } catch (py::error_already_set& error) {
        if (!is_exception(error)) {
            throw;
        }
        return failure(str_of(error.value()), error);
    } catch (const std::exception& error) {
        return failure(fail_text(error), error);
    }
    py::list result;
    result.append(significant_data_list);
    result.append(significant_key_list);
    result.append(target_data_list);
    result.append(target_key_list);
    result.append(reference_points_list_target);
    result.append(reference_points_list_significant);
    pQueue.attr("put")(result);
    return py::none();
}

py::list flatten(const py::list& lists) {
    py::list out;
    for (const py::handle sublist : lists) {
        for (const py::handle item : sublist) {
            out.append(item);
        }
    }
    return out;
}

// np.unique(keys, axis=0, return_index=True) of the 2D str array np.array(keys)
// makes from the sample triplets: the rows in sorted order, each with the
// index of its first occurrence.
py::tuple unique_rows(const py::object& keys_array, const py::list& keys) {
    const py::object np = numpy();
    const py::array array = py::reinterpret_borrow<py::array>(keys_array);
    if (array.ndim() != 2 || array.dtype().kind() != 'U') {
        // no keys (a float64 array of shape (0,)) or other shapes: numpy's own
        const py::object result = np.attr("unique")(keys_array, py::arg("axis") = 0, py::arg("return_index") = true);
        return py::reinterpret_borrow<py::tuple>(result);
    }
    std::vector<std::vector<std::string>> rows;
    for (const py::handle row : array.attr("tolist")()) {
        std::vector<std::string> fields;
        for (const py::handle field : row) {
            fields.push_back(utf8(field));
        }
        rows.push_back(std::move(fields));
    }
    // UCS4 code point order, which for UTF-8 is the byte order
    std::vector<std::size_t> order(rows.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return rows[a] < rows[b]; });
    py::list unique;
    std::vector<std::int64_t> indices;
    for (std::size_t k = 0; k < order.size(); ++k) {
        if (k > 0 && rows[order[k]] == rows[order[k - 1]]) {
            continue;
        }
        unique.append(array.attr("__getitem__")(py::int_(order[k])).attr("tolist")());
        indices.push_back(static_cast<std::int64_t>(order[k]));
    }
    (void)keys;
    const py::object unique_array = np.attr("array")(unique, py::arg("dtype") = array.dtype());
    return py::make_tuple(unique_array, from_vector(indices));
}

// The root attributes of both output files.
void write_root_attributes(hc::Hdf5Writer& writer, const char* type, py::handle pArgs) {
    write_attribute(writer, "/", "type", py::str(type));
    write_attribute(writer, "/", "version", package_version());
    write_attribute(writer, "/", "pvalue", pArgs.attr("pValue"));
    if (!pArgs.attr("xFoldBackground").is_none()) {
        write_attribute(writer, "/", "mode_preselection", py::str("xfold"));
        write_attribute(writer, "/", "mode_preselection_value", pArgs.attr("xFoldBackground"));
    } else if (!pArgs.attr("loosePValue").is_none()) {
        write_attribute(writer, "/", "mode_preselection", py::str("loosePValue"));
        write_attribute(writer, "/", "mode_preselection_value", pArgs.attr("loosePValue"));
    } else {
        write_attribute(writer, "/", "mode_preselection", py::str("None"));
        // the target file spells it mode_preselection_calue
        write_attribute(writer, "/", std::string(type) == "target" ? "mode_preselection_calue" : "mode_preselection_value",
                        py::str("None"));
    }
    write_attribute(writer, "/", "range", pArgs.attr("range"));
    write_attribute(writer, "/", "combinationMode", pArgs.attr("combinationMode"));
    write_attribute(writer, "/", "truncateZeroPvalues", pArgs.attr("truncateZeroPvalues"));
    write_attribute(writer, "/", "fixateRange", pArgs.attr("fixateRange"));
    write_attribute(writer, "/", "peakInteractionsThreshold", pArgs.attr("peakInteractionsThreshold"));
}

// The name a gene's group gets: key, key_1, key_2, ... not yet seen.
std::string unique_name(std::set<std::string>& seen, const std::string& gene) {
    for (std::int64_t counter = 0;; ++counter) {
        const std::string name = counter == 0 ? gene : gene + "_" + std::to_string(counter);
        if (seen.insert(name).second) {
            return name;
        }
    }
}

// Viewpoint.createUniqueHDFGroup(group, name), which loops forever in 3.7.6
// when the name exists.
void create_unique_group(hc::Hdf5Writer& writer, const std::string& path) {
    if (writer.exists(path)) {
        raise_error("ValueError", "the group " + path +
                                      " exists already: Viewpoint.createUniqueHDFGroup does not terminate on it in "
                                      "HiCExplorer 3.7.6");
    }
    writer.create_group(path);
}

}  // namespace

void bind_chicSignificantInteractions(py::module_& m) {
    m.doc() = "C++ bodies of hicexplorer.chicSignificantInteractions.";

    m.def("compute_interaction_file", &compute_interaction_file, py::arg("pInteractionFilesList"), py::arg("pArgs"),
          py::arg("pViewpointObj"), py::arg("pBackground"), py::arg("pFilePath"), py::arg("pResolution"),
          py::arg("pQueue") = py::none());

    m.def("filter_by_pvalue", &filter_by_pvalue, py::arg("pData"), py::arg("pPValue"), py::arg("pMergedLinesDict"),
          py::arg("pPeakInteractionsThreshold"));

    m.def("compute_new_p_values", &compute_new_p_values, py::arg("pData"), py::arg("pBackgroundModel"),
          py::arg("pPValue"), py::arg("pMergedLinesDict"), py::arg("pPeakInteractionsThreshold"));

    m.def("merge_neighbors_x_fold", &merge_x_fold, py::arg("pXfold"), py::arg("pData"), py::arg("pViewpointObj"),
          py::arg("pResolution"));

    m.def("merge_neighbors_loose_p_value", &merge_loose, py::arg("pLoosePValue"), py::arg("pData"),
          py::arg("pViewpointObj"), py::arg("pResolution"), py::arg("pTruncateZeroPvalues"));

    m.def("call_multi_core", [](py::handle pInteractionFilesList, py::handle pArgs, py::handle pViewpointObj,
                                py::handle pBackground, py::handle pFilePath, py::handle pResolution) {
        const py::object threads_object = pArgs.attr("threads");
        const std::int64_t threads = to_int64(threads_object);
        const std::size_t slots = static_cast<std::size_t>(std::max<std::int64_t>(threads, 0));
        std::vector<py::object> parts[6];
        for (auto& part : parts) {
            part.assign(slots, py::none());
        }
        const py::object per_thread = floordiv(py::int_(length(pInteractionFilesList)), threads_object);
        bool fail_flag = false;
        py::object fail_message = py::str("");
        const py::object queue_type = py::module_::import("queue").attr("Queue");
        for (std::int64_t i = 0; i < threads; ++i) {
            const py::object begin = steal(PyNumber_Multiply(py::int_(i).ptr(), per_thread.ptr()));
            py::object files_thread;
            if (i < threads - 1) {
                const py::object stop = steal(PyNumber_Multiply(py::int_(i + 1).ptr(), per_thread.ptr()));
                files_thread = getitem(pInteractionFilesList, py::slice(begin, stop, py::none()));
            } else {
                files_thread = getitem(pInteractionFilesList, py::slice(begin, py::none(), py::none()));
            }
            const py::object queue = queue_type();
            compute_interaction_file(files_thread, pArgs, pViewpointObj, pBackground, pFilePath, pResolution, queue);
            const py::object result = queue.attr("get")();
            if (contains(result, py::str("Fail:"))) {
                fail_flag = true;
                fail_message = getitem(result, py::slice(py::int_(6), py::none(), py::none()));
            } else {
                const std::vector<py::object> values = unpack(result, 6);
                for (int k = 0; k < 6; ++k) {
                    parts[k][static_cast<std::size_t>(i)] = values[static_cast<std::size_t>(k)];
                }
            }
        }
        if (fail_flag) {
            log(kLogger, "error", utf8(str_object(fail_message)));
            python_exit(1);
        }
        const auto as_py_list = [](const std::vector<py::object>& objects) {
            py::list out;
            for (const py::object& object : objects) {
                out.append(object);
            }
            return out;
        };
        const py::object np = numpy();
        const py::list significant_data = flatten(as_py_list(parts[0]));
        const py::list significant_keys = flatten(as_py_list(parts[1]));
        const py::object keys_array = np.attr("array")(significant_keys);
        const py::object reference_points_target = np.attr("array")(flatten(as_py_list(parts[4])));
        const py::object reference_points_significant = np.attr("array")(flatten(as_py_list(parts[5])));
        const py::tuple unique = unique_rows(keys_array, significant_keys);
        py::list significant_data_new;
        for (const py::handle index : unique[1]) {
            significant_data_new.append(getitem(significant_data, index));
        }
        return py::make_tuple(significant_data_new, unique[0], flatten(as_py_list(parts[2])),
                              flatten(as_py_list(parts[3])), reference_points_target, reference_points_significant);
    }, py::arg("pInteractionFilesList"), py::arg("pArgs"), py::arg("pViewpointObj"), py::arg("pBackground"),
       py::arg("pFilePath"), py::arg("pResolution"));

    m.def("read_threshold_file", [](py::handle pFile) {
        py::dict distance_value_dict;
        const std::vector<std::string> lines = read_text_lines(pFile);
        for (const std::string& raw : lines) {
            const py::object line = to_str(raw).attr("strip")();
            if (truthy(line.attr("startswith")("#"))) {
                continue;
            }
            if (eq(line, py::str(""))) {
                break;
            }
            const std::vector<py::object> fields = unpack(line.attr("split")("\t"), 2);
            distance_value_dict[parse_int(fields[0])] = float_of(fields[1]);
        }
        return distance_value_dict;
    }, py::arg("pFile"));

    m.def("writeSignificantHDF", [](py::handle pOutFileName, py::handle pSignificantDataList,
                                    py::handle pSignificantKeyList, py::handle pViewpointObj,
                                    py::handle pReferencePointsList, py::handle pArgs) {
        (void)pViewpointObj;
        hc::Hdf5Writer writer = create_h5(pOutFileName);
        write_root_attributes(writer, "significant", pArgs);
        std::map<std::string, std::set<std::string>> keys_seen;
        py::ssize_t i = 0;
        for (const py::handle pair : builtin("zip")(pSignificantKeyList, pSignificantDataList)) {
            const std::vector<py::object> key_data = unpack(pair, 2);
            const py::object& key = key_data[0];
            const py::object& data = key_data[1];
            if (length(data) == 0) {
                ++i;
                continue;
            }
            py::object chromosome = py::none();
            py::list start_list;
            py::list end_list;
            py::object gene_name = py::none();
            py::object sum_of_interactions = py::none();
            py::list relative_distance_list;
            py::list relative_interactions_list;
            py::list pvalue_list;
            py::list xfold_list;
            py::list raw_target_list;
            for (const py::handle datum : data.attr("values")()) {
                chromosome = getitem(datum, 0);
                start_list.append(getitem(datum, 1));
                end_list.append(getitem(datum, 2));
                gene_name = getitem(datum, 3);
                sum_of_interactions = getitem(datum, 4);
                relative_distance_list.append(getitem(datum, 5));
                relative_interactions_list.append(getitem(datum, 6));
                pvalue_list.append(getitem(datum, 7));
                xfold_list.append(getitem(datum, 8));
                raw_target_list.append(getitem(datum, 9));
            }
            const std::string matrix = utf8(str_object(getitem(key, 0)));
            const std::string matrix_path = h5_path({matrix});
            if (!writer.exists(matrix_path)) {
                writer.create_group(matrix_path);
                keys_seen[matrix] = std::set<std::string>();
            }
            const std::string chromosome_path = matrix_path + "/" + utf8(str_object(chromosome));
            ensure_group(writer, chromosome_path);
            const std::string genes_path = matrix_path + "/genes";
            ensure_group(writer, genes_path);
            const auto seen = keys_seen.find(matrix);
            if (seen == keys_seen.end()) {
                raise_error("KeyError", repr_of(getitem(key, 0)));
            }
            const std::string gene_key = unique_name(seen->second, utf8(str_object(getitem(key, 2))));
            // Viewpoint.writeInteractionFileHDF5(chromosomeObject, gene_name_key, ...)
            const std::string group_path = chromosome_path + "/" + gene_key;
            create_unique_group(writer, group_path);
            write_scalar(writer, group_path + "/chromosome", str_object(chromosome));
            write_values(writer, group_path + "/start_list", start_list, 9);
            write_values(writer, group_path + "/end_list", end_list, 9);
            write_scalar(writer, group_path + "/gene", gene_name);
            write_scalar(writer, group_path + "/sum_of_interactions", sum_of_interactions);
            write_values(writer, group_path + "/relative_position_list", relative_distance_list, 9);
            write_values(writer, group_path + "/interaction_data_list", relative_interactions_list, 9);
            write_values(writer, group_path + "/pvalue", pvalue_list, 9);
            write_values(writer, group_path + "/xfold", xfold_list, 9);
            write_values(writer, group_path + "/raw", raw_target_list, 9);
            const py::object reference_point = getitem(pReferencePointsList, i);
            write_scalar(writer, group_path + "/reference_point_start", int_of(getitem(reference_point, 0)));
            write_scalar(writer, group_path + "/reference_point_end", int_of(getitem(reference_point, 1)));
            if (!writer.hard_link(group_path, genes_path + "/" + gene_key)) {
                log(kLogger, "debug", "exception Unable to create link (name already exists)");
                log(kLogger, "debug", "Gene group given: " + utf8(format(getitem(key, 2), "")));
                log(kLogger, "debug", "Gene group return: " + gene_key);
            }
            ++i;
        }
        writer.close();
    }, py::arg("pOutFileName"), py::arg("pSignificantDataList"), py::arg("pSignificantKeyList"),
       py::arg("pViewpointObj"), py::arg("pReferencePointsList"), py::arg("pArgs"));

    m.def("writeTargetHDF", [](py::handle pOutFileName, py::handle pTargetDataList, py::handle pTargetKeyList,
                               py::handle pViewpointObj, py::handle pResolution, py::handle pReferencePointsList,
                               py::handle pArgs) {
        (void)pViewpointObj;
        hc::Hdf5Writer writer = create_h5(pOutFileName);
        write_root_attributes(writer, "target", pArgs);
        std::map<std::string, std::set<std::string>> keys_seen;
        py::ssize_t i = 0;
        for (const py::handle pair : builtin("zip")(pTargetKeyList, pTargetDataList)) {
            const std::vector<py::object> key_data = unpack(pair, 2);
            const py::object& key = key_data[0];
            const py::object& data = key_data[1];
            if (length(data) == 0) {
                ++i;
                continue;
            }
            // pybedtools.BedTool(data).sort().merge(d=pResolution)
            std::vector<hc::TargetRegion> regions;
            for (const py::handle line : data) {
                const std::vector<py::object> fields =
                    unpack(getitem(line, py::slice(py::none(), py::int_(3), py::none())), 3);
                regions.push_back(hc::TargetRegion{utf8(str_object(fields[0])),
                                                   to_int64(parse_int(str_object(fields[1]))),
                                                   to_int64(parse_int(str_object(fields[2])))});
            }
            const std::vector<hc::TargetRegion> merged = hc::sort_merge_regions(regions, to_int64(pResolution));
            std::string chromosome;
            py::list start_list;
            py::list end_list;
            for (const hc::TargetRegion& region : merged) {
                chromosome = region.chromosome;
                start_list.append(py::str(std::to_string(region.start)));
                end_list.append(py::str(std::to_string(region.end)));
            }
            const std::string outer = utf8(str_object(getitem(key, 0)));
            std::string group = h5_path({outer});
            if (!writer.exists(group)) {
                writer.create_group(group);
                keys_seen[outer] = std::set<std::string>();
            }
            const py::ssize_t n = length(key);
            for (py::ssize_t k = 1; k < n; ++k) {
                const py::object name = getitem(key, k);
                if (eq(name, py::str("::"))) {
                    break;
                }
                const std::string text = utf8(str_object(name));
                const std::string path = group + "/" + text;
                if (!writer.exists(path)) {
                    writer.create_group(path);
                    keys_seen[text] = std::set<std::string>();
                }
                group = path;
            }
            const std::string genes_path = group + "/genes";
            ensure_group(writer, genes_path);
            const std::string chromosome_path = group + "/" + chromosome;
            ensure_group(writer, chromosome_path);
            const auto seen = keys_seen.find(outer);
            if (seen == keys_seen.end()) {
                raise_error("KeyError", repr_of(getitem(key, 0)));
            }
            const py::object last = getitem(key, -1);
            const std::string gene_key = unique_name(seen->second, utf8(str_object(last)));
            const std::string group_path = chromosome_path + "/" + gene_key;
            create_unique_group(writer, group_path);
            write_scalar(writer, group_path + "/chromosome", py::str(chromosome));
            write_values(writer, group_path + "/start_list", start_list, 0);
            write_values(writer, group_path + "/end_list", end_list, 0);
            const py::object reference_point = getitem(pReferencePointsList, i);
            write_scalar(writer, group_path + "/reference_point_start", int_of(getitem(reference_point, 0)));
            write_scalar(writer, group_path + "/reference_point_end", int_of(getitem(reference_point, 1)));
            if (!writer.hard_link(group_path, genes_path + "/" + gene_key)) {
                log(kLogger, "debug", "exception Unable to create link (name already exists)");
                log(kLogger, "debug", "Gene group: " + utf8(format(last, "")));
                log(kLogger, "debug", "Adjusted name " + gene_key);
            }
            ++i;
        }
        writer.close();
    }, py::arg("pOutFileName"), py::arg("pTargetDataList"), py::arg("pTargetKeyList"), py::arg("pViewpointObj"),
       py::arg("pResolution"), py::arg("pReferencePointsList"), py::arg("pArgs"));
}

}  // namespace hicx::python
