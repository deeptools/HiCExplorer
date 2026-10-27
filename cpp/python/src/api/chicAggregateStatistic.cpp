// hicexplorer._hicx.chicAggregateStatistic: the C++ bodies of the functions
// of hicexplorer.chicAggregateStatistic besides main() and parse_arguments()
// (which run the executable cpp/tools/chicAggregateStatistic.cpp).
//
// filter_scores_target_list keeps, per target region, the scores of the
// positions that overlap it, summed. The target regions of an HDF5 target
// file, of a bed4 target dict or of a bed3 file (run_target_list_compilation)
// are put into a C++ interval index with the semantics of hicmatrix's
// intervalListToIntervalTree and intervaltree 3: a position [start, end)
// overlaps the regions with begin < end and end > start, and the first of
// them in (begin, end, index) order is the target. An interval tree the
// caller passes (pTargetIntervalTree) is queried through its own interface.
// The aggregate file is written with hicx::chic::Hdf5Writer in h5py's layout.
//
// call_multi_core ran pFunctionName in one process per thread and collected
// the results through queues, waiting one second per round. Here the
// function is called once per thread in this process with a queue.Queue,
// with the same arguments, and the results are collected in thread order;
// there is no waiting. A worker that raised instead of putting its result
// made 3.7.6 wait forever; here its exception propagates.

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

#include "api.hpp"
#include "chic_common.hpp"
#include "hicx/chic_hdf5.hpp"
#include "hicx/utilities.hpp"

namespace hicx::python {

using namespace chic_api;
namespace hc = hicx::chic;

namespace {

constexpr const char* kLogger = "hicexplorer.chicAggregateStatistic";

// hicmatrix's intervalListToIntervalTree over (chromosome, start, end), as an
// index: per chromosome object the intervals with their position in the list.
struct Interval {
    std::int64_t begin = 0;
    std::int64_t end = 0;
    std::int64_t index = 0;
};

class IntervalIndex {
  public:
    // intervalListToIntervalTree(regions)[0]
    explicit IntervalIndex(py::handle regions) {
        std::int64_t index = 0;
        for (const py::handle region : regions) {
            const std::vector<py::object> fields = unpack(getitem(region, py::slice(py::int_(0), py::int_(3), py::none())), 3);
            const std::int64_t begin = to_int64(int_of(fields[1]));
            const std::int64_t end = to_int64(int_of(fields[2]));
            if (begin >= end) {
                raise_error("ValueError", "IntervalTree: Null Interval objects not allowed in IntervalTree: Interval(" +
                                              std::to_string(begin) + ", " + std::to_string(end) + ", " +
                                              std::to_string(index) + ")");
            }
            py::object slot;
            if (contains(chromosomes_, fields[0])) {
                slot = getitem(chromosomes_, fields[0]);
            } else {
                slot = py::int_(static_cast<py::ssize_t>(intervals_.size()));
                chromosomes_[fields[0]] = slot;
                intervals_.emplace_back();
            }
            intervals_[static_cast<std::size_t>(to_int64(slot))].push_back(Interval{begin, end, index});
            ++index;
        }
    }

    [[nodiscard]] bool has(py::handle chromosome) const { return contains(chromosomes_, chromosome); }

    // sorted(tree[chromosome][start:end])[0], or nullopt for no overlap
    [[nodiscard]] std::optional<Interval> first_overlap(py::handle chromosome, std::int64_t start,
                                                        std::int64_t end) const {
        if (start >= end) {
            return std::nullopt;
        }
        const auto& list = intervals_[static_cast<std::size_t>(to_int64(getitem(chromosomes_, chromosome)))];
        std::optional<Interval> best;
        for (const Interval& interval : list) {
            if (interval.begin < end && interval.end > start) {
                if (!best || std::tie(interval.begin, interval.end, interval.index) <
                                 std::tie(best->begin, best->end, best->index)) {
                    best = interval;
                }
            }
        }
        return best;
    }

  private:
    py::dict chromosomes_;
    std::vector<std::vector<Interval>> intervals_;
};

// sorted(keys) for the keys of the scores dict
py::list sorted_list(const py::list& keys) {
    py::dict ordered;
    for (const py::handle key : keys) {
        ordered[key] = py::none();
    }
    if (length(ordered) == length(keys)) {
        return sorted_keys(ordered);
    }
    return builtin("sorted")(keys);
}

py::dict filter_scores(py::handle pScoresDictionary, py::handle pTargetFType, py::handle pTargetPosDict,
                       py::handle pTargetList, py::handle pTargetIntervalTree, py::handle pTargetFile,
                       const IntervalIndex* given_index) {
    py::dict accepted_scores;
    std::optional<py::object> chromosome;
    std::optional<py::object> start_list;
    std::optional<py::object> end_list;
    py::object python_tree = py::none();
    const IntervalIndex* index = nullptr;
    if (eq(pTargetFType, py::str("hdf5"))) {
        hc::H5Reader reader = open_h5(pTargetFile);
        const std::string path = utf8(join("/", pTargetList));
        if (!reader.contains(path)) {
            raise_error("KeyError", reader.missing_message(path));
        }
        if (!reader.contains(path + "/chromosome")) {
            raise_error("TypeError", "'NoneType' object is not subscriptable");
        }
        chromosome = decode_utf8(h5_scalar(reader.read(path + "/chromosome")));
        start_list = as_list(h5_array(reader.read(path + "/start_list")));
        end_list = as_list(h5_array(reader.read(path + "/end_list")));
        reader.close();
    } else if (eq(pTargetFType, py::str("bed4"))) {
        const py::object entry = getitem(pTargetPosDict, getitem(pTargetList, -1));
        chromosome = getitem(entry, py::str("chromosome"));
        start_list = getitem(entry, py::str("start_list"));
        end_list = getitem(entry, py::str("end_list"));
    } else if (eq(pTargetFType, py::str("bed3"))) {
        if (given_index != nullptr) {
            index = given_index;
        } else {
            python_tree = py::reinterpret_borrow<py::object>(pTargetIntervalTree);
        }
    } else {
        log(kLogger, "error", "No target list given.");
        raise_error("Exception", "No target list given.");
    }
    std::optional<IntervalIndex> own_index;
    if (!pTargetList.is_none()) {
        if (!chromosome) {
            raise_error("UnboundLocalError",
                        "cannot access local variable 'chromosome' where it is not associated with a value");
        }
        const py::ssize_t n = length(*start_list);
        py::list chromosomes;
        for (py::ssize_t k = 0; k < n; ++k) {
            chromosomes.append(*chromosome);
        }
        const py::list target_regions = as_list(builtin("zip")(chromosomes, *start_list, *end_list));
        if (target_regions.empty()) {
            return accepted_scores;
        }
        own_index.emplace(target_regions);
        index = &*own_index;
        python_tree = py::none();
    }
    // per target the keys of its positions, in the order the targets are met
    py::dict same_target;
    for (const py::handle key : pScoresDictionary) {
        const py::object line = getitem(pScoresDictionary, key);
        const py::object chrom = getitem(line, 0);
        const std::int64_t start = to_int64(int_of(getitem(line, 1)));
        const std::int64_t end = to_int64(int_of(getitem(line, 2)));
        py::object target;
        if (index != nullptr) {
            if (!index->has(chrom)) {
                continue;
            }
            const std::optional<Interval> hit = index->first_overlap(chrom, start, end);
            if (!hit) {
                continue;
            }
            target = py::int_(hit->index);
        } else {
            if (python_tree.is_none()) {
                raise_error("TypeError", "argument of type 'NoneType' is not iterable");
            }
            if (!contains(python_tree, chrom)) {
                continue;
            }
            const py::object hits = getitem(getitem(python_tree, chrom), py::slice(py::int_(start), py::int_(end), py::none()));
            if (!truthy(hits)) {
                continue;
            }
            target = getitem(builtin("sorted")(hits), 0);
        }
        if (contains(same_target, target)) {
            getitem(same_target, target).attr("append")(key);
        } else {
            py::list keys;
            keys.append(key);
            same_target[target] = keys;
        }
    }
    for (const py::handle target : as_list(same_target)) {
        double values[3] = {0.0, 0.0, 0.0};
        const py::list keys = sorted_list(py::reinterpret_borrow<py::list>(getitem(same_target, target)));
        same_target[target] = keys;
        for (const py::handle key : keys) {
            const py::object last = getitem(getitem(pScoresDictionary, key), py::slice(py::int_(-3), py::none(), py::none()));
            std::vector<double> numbers;
            for (const py::handle x : last) {
                numbers.push_back(to_double(float_of(x)));
            }
            if (numbers.size() != 3 && numbers.size() != 1) {
                raise_error("ValueError", "operands could not be broadcast together with shapes (3,) (" +
                                              std::to_string(numbers.size()) + ",) (3,) ");
            }
            for (std::size_t k = 0; k < 3; ++k) {
                values[k] += numbers.size() == 1 ? numbers[0] : numbers[k];
            }
        }
        const py::object first_key = keys[0];
        const py::object last_key = keys[keys.size() - 1];
        const py::object new_data_line = getitem(pScoresDictionary, first_key);
        setitem(new_data_line, 2, getitem(getitem(pScoresDictionary, last_key), 2));
        setitem(new_data_line, -5, getitem(getitem(pScoresDictionary, last_key), -5));
        setitem(new_data_line, -3, np_float64(values[0]));
        setitem(new_data_line, -2, np_float64(values[1]));
        setitem(new_data_line, -1, np_float64(values[2]));
        accepted_scores[first_key] = new_data_line;
    }
    return accepted_scores;
}

void put_failure(py::handle queue, const std::exception& error) {
    std::string text;
    if (const auto* set = dynamic_cast<const py::error_already_set*>(&error)) {
        text = str_of(set->value());
    } else {
        text = fail_text(error);
    }
    queue.attr("put")(to_str("Fail: " + text + traceback_text(error)));
}

py::object run_target_list(py::handle pInteractionFilesList, py::handle pTargetList, py::handle pTargetFType,
                           py::handle pTargetPosDict, py::handle pArgs, py::handle pViewpointObj, py::handle pQueue,
                           py::handle pOneTarget) {
    py::list outfile_names_list;
    py::list accepted_scores_list;
    std::optional<py::list> accepted_scores_list_intern;
    std::optional<IntervalIndex> index;
    try {
        if (eq(pTargetFType, py::str("bed3"))) {
            py::list target_regions;
            try {
                for (const auto& region : hicx::utilities::read_bed(read_text_lines(pTargetList))) {
                    target_regions.append(py::make_tuple(to_str(region[0]), to_str(region[1]), to_str(region[2])));
                }
            } catch (py::error_already_set& error) {
                if (!is_exception(error)) {
                    throw;
                }
                put_failure(pQueue, error);
                return py::none();
            } catch (const std::exception& error) {
                put_failure(pQueue, error);
                return py::none();
            }
            index.emplace(target_regions);
        }
        py::ssize_t i = 0;
        for (const py::handle interaction_file : pInteractionFilesList) {
            py::list names;
            py::list scores;
            for (const py::handle sample : interaction_file) {
                const std::vector<py::object> data =
                    unpack(pViewpointObj.attr("readInteractionFile")(pArgs.attr("interactionFile"), sample), 3);
                const py::object target_file = eq(pOneTarget, py::bool_(true))
                                                   ? py::object(py::none())
                                                   : getitem(pTargetList, i);
                const py::dict accepted = filter_scores(data[1], pTargetFType, pTargetPosDict, target_file,
                                                        py::none(), pArgs.attr("targetFile"), index ? &*index : nullptr);
                names.append(sample);
                scores.append(accepted);
            }
            outfile_names_list.append(names);
            accepted_scores_list.append(scores);
            accepted_scores_list_intern = scores;
            ++i;
        }
    } catch (py::error_already_set& error) {
        if (!is_exception(error)) {
            throw;
        }
        put_failure(pQueue, error);
        return py::none();
    } catch (const std::exception& error) {
        put_failure(pQueue, error);
        return py::none();
    }
    if (pQueue.is_none()) {
        return py::none();
    }
    if (!accepted_scores_list_intern) {
        raise_error("UnboundLocalError",
                    "cannot access local variable 'accepted_scores_list_intern' where it is not associated with a "
                    "value");
    }
    py::list result;
    result.append(outfile_names_list);
    result.append(accepted_scores_list);
    pQueue.attr("put")(result);
    return py::none();
}

}  // namespace

void bind_chicAggregateStatistic(py::module_& m) {
    m.doc() = "C++ bodies of hicexplorer.chicAggregateStatistic.";

    m.def("filter_scores_target_list", [](py::handle pScoresDictionary, py::handle pTargetFType, py::handle pTargetPosDict,
                                          py::handle pTargetList, py::handle pTargetIntervalTree, py::handle pTargetFile) {
        return filter_scores(pScoresDictionary, pTargetFType, pTargetPosDict, pTargetList, pTargetIntervalTree,
                             pTargetFile, nullptr);
    }, py::arg("pScoresDictionary"), py::arg("pTargetFType"), py::arg("pTargetPosDict"),
       py::arg("pTargetList") = py::none(), py::arg("pTargetIntervalTree") = py::none(),
       py::arg("pTargetFile") = py::none());

    m.def("writeAggregateHDF", [](py::handle pOutFileName, py::handle pOutfileNamesList, py::handle pAcceptedScoresList,
                                  py::handle pArgs) {
        (void)pArgs;
        hc::Hdf5Writer writer = create_h5(pOutFileName);
        write_attribute(writer, "/", "type", py::str("aggregate"));
        write_attribute(writer, "/", "version", package_version());
        const py::object pairs = builtin("zip")(pOutfileNamesList, pAcceptedScoresList);
        for (const py::handle pair : pairs) {
            const std::vector<py::object> outer = unpack(pair, 2);
            const py::object combination =
                add(add(getitem(getitem(outer[0], 0), 0), py::str("_")), getitem(getitem(outer[0], 1), 0));
            const std::string combination_path = h5_path({utf8(combination)});
            ensure_group(writer, combination_path);
            for (const py::handle inner : builtin("zip")(outer[0], outer[1])) {
                const std::vector<py::object> key_data = unpack(inner, 2);
                const py::object& key = key_data[0];
                const py::object& data = key_data[1];
                if (length(data) == 0) {
                    continue;
                }
                py::object chromosome = py::none();
                py::list start_list;
                py::list end_list;
                py::object gene_name = py::none();
                py::object sum_of_interactions = py::none();
                py::list relative_distance_list;
                py::list raw_target_list;
                for (const py::handle accepted : data) {
                    const py::object line = getitem(data, accepted);
                    chromosome = getitem(line, 0);
                    start_list.append(getitem(line, 1));
                    end_list.append(getitem(line, 2));
                    gene_name = getitem(line, 3);
                    sum_of_interactions = getitem(line, 4);
                    relative_distance_list.append(getitem(line, 5));
                    raw_target_list.append(getitem(line, -1));
                }
                const std::string matrix_path = combination_path + "/" + utf8(getitem(key, 0));
                ensure_group(writer, matrix_path);
                std::string chromosome_path;
                const std::string key_chromosome = utf8(getitem(key, 1));
                if (!writer.exists(matrix_path + "/" + key_chromosome)) {
                    chromosome_path = matrix_path + "/" + key_chromosome;
                    writer.create_group(chromosome_path);
                } else {
                    chromosome_path = matrix_path + "/" + utf8(chromosome);
                    if (!writer.exists(chromosome_path)) {
                        raise_error("KeyError", "Unable to synchronously open object (object '" + utf8(chromosome) +
                                                    "' doesn't exist)");
                    }
                }
                const std::string genes_path = matrix_path + "/genes";
                ensure_group(writer, genes_path);
                const std::string gene_key = utf8(getitem(key, -1));
                const std::string group_path = chromosome_path + "/" + gene_key;
                create_group(writer, group_path);
                write_scalar(writer, group_path + "/chromosome", chromosome);
                write_values(writer, group_path + "/start_list", start_list, 9);
                write_values(writer, group_path + "/end_list", end_list, 9);
                write_scalar(writer, group_path + "/gene_name", gene_name);
                write_scalar(writer, group_path + "/sum_of_interactions", sum_of_interactions);
                write_values(writer, group_path + "/relative_distance_list", relative_distance_list, 9);
                write_values(writer, group_path + "/raw_target_list", raw_target_list, 9);
                hard_link(writer, group_path, genes_path + "/" + gene_key);
            }
        }
        writer.close();
    }, py::arg("pOutFileName"), py::arg("pOutfileNamesList"), py::arg("pAcceptedScoresList"), py::arg("pArgs"));

    m.def("run_target_list_compilation", &run_target_list, py::arg("pInteractionFilesList"), py::arg("pTargetList"),
          py::arg("pTargetFType"), py::arg("pTargetPosDict"), py::arg("pArgs"), py::arg("pViewpointObj"),
          py::arg("pQueue") = py::none(), py::arg("pOneTarget") = false);

    m.def("call_multi_core", [](py::handle pInteractionFilesList, py::handle pTargetFileList, py::handle pTargetFType,
                                py::handle pTargetPosDict, py::handle pFunctionName, py::handle pArgs,
                                py::handle pViewpointObj) {
        const py::ssize_t files = length(pInteractionFilesList);
        if (lt(py::int_(files), pArgs.attr("threads"))) {
            pArgs.attr("threads") = py::int_(files);
        }
        const py::object threads_object = pArgs.attr("threads");
        const std::int64_t threads = to_int64(threads_object);
        py::list outfile_names_list;
        py::list accepted_scores_list;
        for (std::int64_t i = 0; i < std::max<std::int64_t>(threads, 0); ++i) {
            outfile_names_list.append(py::none());
            accepted_scores_list.append(py::none());
        }
        // len(pInteractionFilesList) // pArgs.threads
        const py::object per_thread = floordiv(py::int_(files), threads_object);
        const bool one_target = length(pTargetFileList) == 1;
        bool fail_flag = false;
        py::object fail_message = py::str("");
        const py::object queue_type = py::module_::import("queue").attr("Queue");
        for (std::int64_t i = 0; i < threads; ++i) {
            py::object files_thread;
            py::object targets_thread;
            const py::object begin = steal(PyNumber_Multiply(py::int_(i).ptr(), per_thread.ptr()));
            if (i < threads - 1) {
                const py::object stop = steal(PyNumber_Multiply(py::int_(i + 1).ptr(), per_thread.ptr()));
                files_thread = getitem(pInteractionFilesList, py::slice(begin, stop, py::none()));
                targets_thread = one_target ? getitem(pTargetFileList, 0)
                                            : getitem(pTargetFileList, py::slice(begin, stop, py::none()));
            } else {
                files_thread = getitem(pInteractionFilesList, py::slice(begin, py::none(), py::none()));
                targets_thread = one_target ? getitem(pTargetFileList, 0)
                                            : getitem(pTargetFileList, py::slice(begin, py::none(), py::none()));
            }
            const py::object queue = queue_type();
            pFunctionName(py::arg("pInteractionFilesList") = files_thread, py::arg("pTargetList") = targets_thread,
                          py::arg("pTargetFType") = pTargetFType, py::arg("pTargetPosDict") = pTargetPosDict,
                          py::arg("pArgs") = pArgs, py::arg("pViewpointObj") = pViewpointObj,
                          py::arg("pQueue") = queue, py::arg("pOneTarget") = py::bool_(one_target));
            if (truthy(queue.attr("empty")())) {
                raise_error("RuntimeError", "call_multi_core: the worker function put no result into its queue");
            }
            const py::object result = queue.attr("get")();
            if (contains(result, py::str("Fail:"))) {
                fail_flag = true;
                fail_message = getitem(result, py::slice(py::int_(6), py::none(), py::none()));
            } else {
                const std::vector<py::object> parts = unpack(result, 2);
                setitem(outfile_names_list, static_cast<py::ssize_t>(i), parts[0]);
                setitem(accepted_scores_list, static_cast<py::ssize_t>(i), parts[1]);
            }
        }
        if (fail_flag) {
            log(kLogger, "error", utf8(str_object(fail_message)));
            python_exit(1);
        }
        py::list names;
        py::list scores;
        for (const py::handle sublist : outfile_names_list) {
            for (const py::handle item : sublist) {
                names.append(item);
            }
        }
        for (const py::handle sublist : accepted_scores_list) {
            for (const py::handle item : sublist) {
                scores.append(item);
            }
        }
        return py::make_tuple(names, scores);
    }, py::arg("pInteractionFilesList"), py::arg("pTargetFileList"), py::arg("pTargetFType"), py::arg("pTargetPosDict"),
       py::arg("pFunctionName"), py::arg("pArgs"), py::arg("pViewpointObj"));
}

}  // namespace hicx::python
