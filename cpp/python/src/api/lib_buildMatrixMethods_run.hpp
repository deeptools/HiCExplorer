// process_data and createMatrix of hicexplorer._hicx.lib_buildMatrixMethods.
// Included once, by lib_buildMatrixMethods.cpp, after its helpers.
//
// process_data classifies the pairs of its two buffers with
// hicx::build_matrix::classify_pairs, the loop the build matrix tools run,
// through a policy that reads the pysam objects and the ctypes arrays it is
// given. createMatrix is hicx::build_matrix::run with Arguments::python_api,
// its hooks doing what 3.7.6 did through Python objects: logging, closing and
// removing the file objects, reading pChromosomeSizes and the restriction
// site files, resolving pRegion with getUserRegion, upper-casing the sequence
// lists in place and calling hicPrepareQCreport.main.

#ifndef HICX_PYTHON_API_LIB_BUILD_MATRIX_METHODS_RUN_HPP
#define HICX_PYTHON_API_LIB_BUILD_MATRIX_METHODS_RUN_HPP

#include <unistd.h>

#include <cerrno>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicx/build_matrix.hpp"
#include "hicx/build_matrix_run.hpp"
#include "hicx/py_text.hpp"

namespace hicx::python {

namespace {

// --- ctypes arrays ---------------------------------------------------------------

// A ctypes array (RawArray, or the object of a synchronized Array) read and
// written through its buffer, with the indexing rules of ctypes arrays
// (negative indices from the end, IndexError "invalid index"); any other
// sequence is indexed with Python's getitem and setitem.
class SharedArray {
  public:
    // `element` is the ctypes type the fast path requires, `itemsize` its size.
    SharedArray(py::handle object, py::handle element, std::size_t itemsize) {
        py::object target = py::reinterpret_borrow<py::object>(object);
        if (py::hasattr(target, "get_obj") && !py::isinstance(target, ctypes_array())) {
            target = target.attr("get_obj")();
        }
        object_ = target;
        if (py::isinstance(target, ctypes_array()) && py::hasattr(target, "_type_") &&
            type_matches(target.attr("_type_"), element)) {
            if (PyObject_GetBuffer(target.ptr(), &view_, PyBUF_RECORDS) == 0) {
                // ctypes exports C-contiguous buffers without strides
                if (view_.ndim == 1 && view_.shape != nullptr &&
                    static_cast<std::size_t>(view_.itemsize) == itemsize &&
                    (view_.strides == nullptr || view_.strides[0] == view_.itemsize)) {
                    fast_ = true;
                } else {
                    PyBuffer_Release(&view_);
                }
            } else {
                PyErr_Clear();
            }
        }
    }
    SharedArray(const SharedArray&) = delete;
    SharedArray& operator=(const SharedArray&) = delete;
    ~SharedArray() {
        if (fast_) {
            PyBuffer_Release(&view_);
        }
    }

    [[nodiscard]] bool fast() const noexcept { return fast_; }
    [[nodiscard]] const py::object& object() const noexcept { return object_; }

    // The address of element `index` (fast path only).
    [[nodiscard]] char* at(std::int64_t index) const {
        const std::int64_t n = static_cast<std::int64_t>(view_.shape[0]);
        if (index < 0) {
            index += n;
        }
        if (index < 0 || index >= n) {
            raise_error("IndexError", "invalid index");
        }
        return static_cast<char*>(view_.buf) + index * view_.itemsize;
    }

  private:
    static py::object ctypes_array() { return py::module_::import("ctypes").attr("Array"); }
    static bool type_matches(py::handle type, py::handle element) {
        if (type.ptr() == element.ptr()) {
            return true;
        }
        const int sub = PyObject_IsSubclass(type.ptr(), element.ptr());
        if (sub < 0) {
            PyErr_Clear();
            return false;
        }
        return sub == 1;
    }

    py::object object_;
    Py_buffer view_{};
    bool fast_ = false;
};

std::uint32_t read_uint32(const char* address) {
    std::uint32_t value = 0;
    std::memcpy(&value, address, sizeof(value));
    return value;
}

// --- process_data -------------------------------------------------------------------

class ProcessDataPolicy {
  public:
    struct Inputs {
        py::object buffer1;
        py::object buffer2;
        py::object restriction_sequence;
        std::vector<py::object> sequences;  // the items of pRestrictionSequence
        py::object rf_positions;
        py::object ref_id2name;
        py::object dangling_sequences;
        py::object bin_index;
        py::object output_bam_set;
        py::object out_bam_index_buffer;
        hicx::QcCounters* counters = nullptr;
    };

    ProcessDataPolicy(Inputs inputs, py::handle shared_tree, py::handle coverage, py::handle coverage_index,
                      py::handle row, py::handle col, py::handle data, py::handle module)
        : in_(std::move(inputs)),
          tree_(shared_tree, module.attr("C_Interval"), 12),
          coverage_(coverage, py::module_::import("ctypes").attr("c_uint"), 4),
          coverage_index_(coverage_index, module.attr("C_Coverage"), 8),
          row_(row, py::module_::import("ctypes").attr("c_uint"), 4),
          col_(col, py::module_::import("ctypes").attr("c_uint"), 4),
          data_(data, py::module_::import("ctypes").attr("c_ushort"), 2) {
        rf_enabled_ = truthy(in_.rf_positions) && truthy(in_.restriction_sequence);
        output_bam_ = truthy(in_.output_bam_set);
    }

    void set_quick_qc(bool quick) { quick_qc_ = quick; }

    const py::object& mate(std::size_t i, int w) {
        if (!current_ || *current_ != i) {
            mate1_ = get_index(in_.buffer1, static_cast<std::int64_t>(i));
            mate2_ = get_index(in_.buffer2, static_cast<std::int64_t>(i));
            current_ = i;
        }
        return w == 0 ? mate1_ : mate2_;
    }
    std::int64_t tid(std::size_t i, int w) {
        return integer_argument(mate(i, w).attr("reference_id"), "process_data", "reference_id");
    }
    std::int64_t pos(std::size_t i, int w) {
        return integer_argument(mate(i, w).attr("pos"), "process_data", "pos");
    }
    std::int64_t qlen(std::size_t i, int w) {
        return integer_argument(mate(i, w).attr("qlen"), "process_data", "qlen");
    }
    std::int64_t seq_len(std::size_t i, int w) { return static_cast<std::int64_t>(length(mate(i, w).attr("seq"))); }
    bool is_reverse(std::size_t i, int w) { return truthy(mate(i, w).attr("is_reverse")); }

    std::optional<hicx::SearchInterval> bin(std::size_t i, int w, std::int64_t read_middle) {
        const py::object mate_ref = get_item(in_.ref_id2name, mate(i, w).attr("rname"));
        PyObject* range = PyObject_GetItem(in_.bin_index.ptr(), mate_ref.ptr());
        if (range == nullptr) {
            if (PyErr_ExceptionMatches(PyExc_KeyError)) {
                // small contigs without a restriction site have no bins
                PyErr_Clear();
                return std::nullopt;
            }
            throw_python_error();
        }
        const std::vector<py::object> bounds = unpack(py::reinterpret_steal<py::object>(range), 2);
        const std::int64_t start = integer_argument(bounds[0], "process_data", "a bin index");
        const std::int64_t end = integer_argument(bounds[1], "process_data", "a bin index");
        const auto found = hicx::search_bins(start, end, static_cast<std::uint64_t>(read_middle),
                                             [this](std::int64_t index) { return interval(index); });
        if (!found) {
            return std::nullopt;
        }
        return interval(*found);
    }

    [[nodiscard]] bool rf_enabled() const { return rf_enabled_; }

    bool rf_overlaps(std::size_t i, std::int64_t begin, std::int64_t end) {
        const py::object mate_ref = get_item(in_.ref_id2name, mate(i, 0).attr("rname"));
        const int contained = PySequence_Contains(in_.rf_positions.ptr(), mate_ref.ptr());
        if (contained < 0) {
            throw_python_error();
        }
        if (contained == 0) {
            return false;
        }
        // pRfPositions[mate_ref][frag_start:frag_end], an IntervalTree: read
        // once per chromosome into the index the tools use.
        std::size_t slot = 0;
        if (rf_slots_.contains(mate_ref)) {
            slot = rf_slots_[mate_ref].cast<std::size_t>();
        } else {
            std::vector<hicx::GenomeInterval> sites;
            for (py::handle interval : py::iter(get_item(in_.rf_positions, mate_ref))) {
                sites.push_back({0, integer_argument(interval.attr("begin"), "process_data", "a site begin"),
                                 integer_argument(interval.attr("end"), "process_data", "a site end")});
            }
            slot = rf_indices_.size();
            rf_indices_.push_back(std::make_unique<hicx::RestrictionSiteIndex>(sites, 1));
            rf_slots_[mate_ref] = py::int_(slot);
        }
        return rf_indices_[slot]->overlaps(0, begin, end);
    }

    bool dangling(std::size_t i, int w, std::size_t s) {
        return check_dangling_end(mate(i, w), get_item(in_.dangling_sequences, in_.sequences[s]));
    }

    std::pair<std::int64_t, std::int64_t> coverage_range(std::uint32_t bin_id) {
        if (coverage_index_.fast()) {
            const char* address = coverage_index_.at(bin_id);
            return {read_uint32(address), read_uint32(address + 4)};
        }
        const py::object entry = get_index(coverage_index_.object(), bin_id);
        return {integer_argument(entry.attr("begin"), "process_data", "a coverage begin"),
                integer_argument(entry.attr("end"), "process_data", "a coverage end")};
    }

    void coverage_add(std::int64_t cell) {
        if (coverage_.fast()) {
            char* address = coverage_.at(cell);
            const std::uint32_t value = read_uint32(address) + 1;
            std::memcpy(address, &value, sizeof(value));
            return;
        }
        const py::object value = get_index(coverage_.object(), cell);
        PyObject* sum = PyNumber_Add(value.ptr(), py::int_(1).ptr());
        if (sum == nullptr) {
            throw_python_error();
        }
        set_index(coverage_.object(), cell, py::reinterpret_steal<py::object>(sum));
    }

    void set_isize(std::size_t i) {
        const py::object pos1 = mate(i, 0).attr("pos");
        const py::object pos2 = mate(i, 1).attr("pos");
        const auto difference = [](const py::object& a, const py::object& b) {
            PyObject* value = PyNumber_Subtract(a.ptr(), b.ptr());
            if (value == nullptr) {
                throw_python_error();
            }
            return py::reinterpret_steal<py::object>(value);
        };
        mate(i, 0).attr("isize") = difference(pos2, pos1);
        mate(i, 1).attr("isize") = difference(pos1, pos2);
    }

    void add_pair(std::size_t i, std::uint32_t bin1, std::uint32_t bin2) {
        const std::int64_t index = in_.counters->pair_added;
        if (!quick_qc_) {
            write(row_, index, bin1, 4);
            write(col_, index, bin2, 4);
            write(data_, index, 1, 2);
        }
        if (output_bam_) {
            in_.out_bam_index_buffer.attr("append")(py::int_(i));
        }
    }

  private:
    hicx::SearchInterval interval(std::int64_t index) {
        if (tree_.fast()) {
            const char* address = tree_.at(index);
            return {read_uint32(address), read_uint32(address + 4), read_uint32(address + 8)};
        }
        const py::object entry = get_index(tree_.object(), index);
        return {static_cast<std::uint32_t>(integer_argument(entry.attr("begin"), "process_data", "a bin begin")),
                static_cast<std::uint32_t>(integer_argument(entry.attr("end"), "process_data", "a bin end")),
                static_cast<std::uint32_t>(integer_argument(entry.attr("data"), "process_data", "a bin id"))};
    }

    static void write(const SharedArray& array, std::int64_t index, std::uint32_t value, std::size_t width) {
        if (array.fast()) {
            char* address = array.at(index);
            if (width == 4) {
                std::memcpy(address, &value, 4);
            } else {
                const auto narrow = static_cast<std::uint16_t>(value);
                std::memcpy(address, &narrow, 2);
            }
            return;
        }
        set_index(array.object(), index, py::int_(value));
    }

    Inputs in_;
    SharedArray tree_;
    SharedArray coverage_;
    SharedArray coverage_index_;
    SharedArray row_;
    SharedArray col_;
    SharedArray data_;
    bool rf_enabled_ = false;
    bool output_bam_ = false;
    bool quick_qc_ = false;
    std::optional<std::size_t> current_;
    py::object mate1_;
    py::object mate2_;
    py::dict rf_slots_;
    std::vector<std::unique_ptr<hicx::RestrictionSiteIndex>> rf_indices_;
};

// 'Fail: ' + str(exp) + traceback.format_exc() of the except clause.
std::string failure_text(py::error_already_set& error) {
    const std::string message = str_of(error.value());
    std::string trace;
    try {
        const py::object lines = py::module_::import("traceback")
                                     .attr("format_exception")(error.type(), error.value(), error.trace());
        trace = utf8(py::str("").attr("join")(lines));
    } catch (py::error_already_set&) {
        trace = "Traceback (most recent call last):\n";
    }
    return "Fail: " + message + trace;
}

std::string failure_text(const std::string& type, const std::string& message) {
    return "Fail: " + message +
           "Traceback (most recent call last):\n  File \"hicexplorer._hicx.lib_buildMatrixMethods\", in "
           "process_data\n" + type + ": " + message + "\n";
}

void process_data(py::module_ module, py::handle pMateBuffer1, py::handle pMateBuffer2, py::handle /*pMinMappingQuality*/,
                  py::handle pKeepSelfCircles, py::handle pRestrictionSequence, py::handle pKeepSelfLigation,
                  py::handle /*pMatrixSize*/, py::handle pRfPositions, py::handle pRefId2name,
                  py::handle pDanglingSequences, py::handle pBinsize, py::handle pResultIndex, py::handle pQueueOut,
                  py::handle /*pTemplate*/, py::handle pOutputBamSet, py::handle pCounter,
                  py::handle pSharedBinIntvalTree, py::handle pDictBinIntervalTreeIndex, py::handle pCoverage,
                  py::handle pCoverageIndex, py::handle /*pOutputFileBufferDir*/, py::handle pRow, py::handle pCol,
                  py::handle pData, py::handle pMaxInsertSize, py::handle pQuickQCMode) {
    std::string failure;
    try {
        // dangling_end = {sequence: 0 for sequence in pRestrictionSequence}
        py::dict dangling_end;
        std::vector<py::object> sequences;
        if (!pRestrictionSequence.is_none()) {
            for (py::handle sequence : py::iter(pRestrictionSequence)) {
                dangling_end[sequence] = py::int_(0);
                sequences.push_back(py::reinterpret_borrow<py::object>(sequence));
            }
        }
        py::list out_bam_index_buffer;
        if (pMateBuffer1.is_none() || pMateBuffer2.is_none()) {
            py::list counts;
            for (int k = 0; k < 3; ++k) {
                counts.append(py::int_(0));
            }
            counts.append(dangling_end);
            for (int k = 0; k < 13; ++k) {
                counts.append(py::int_(0));
            }
            counts.append(pResultIndex);
            counts.append(out_bam_index_buffer);
            py::list message;
            message.append(py::none());
            message.append(counts);
            pQueueOut.attr("put")(message);
            return;
        }

        hicx::build_matrix::ClassifyOptions options;
        if (truthy(pRestrictionSequence)) {
            for (const py::object& sequence : sequences) {
                options.restriction_sequences.push_back(str_argument(sequence, "process_data", "the restriction sequences"));
            }
        }
        options.has_dangling = truthy(pDanglingSequences);
        options.keep_self_circles = truthy(pKeepSelfCircles);
        options.keep_self_ligation = truthy(pKeepSelfLigation);
        options.max_insert_size = integer_argument(pMaxInsertSize, "process_data", "pMaxInsertSize");
        options.quick_qc_mode = truthy(pQuickQCMode);
        options.coverage_bin_size = integer_argument(pBinsize, "process_data", "pBinsize");

        hicx::QcCounters counters;
        ProcessDataPolicy::Inputs inputs;
        inputs.buffer1 = py::reinterpret_borrow<py::object>(pMateBuffer1);
        inputs.buffer2 = py::reinterpret_borrow<py::object>(pMateBuffer2);
        inputs.restriction_sequence = py::reinterpret_borrow<py::object>(pRestrictionSequence);
        inputs.sequences = sequences;
        inputs.rf_positions = py::reinterpret_borrow<py::object>(pRfPositions);
        inputs.ref_id2name = py::reinterpret_borrow<py::object>(pRefId2name);
        inputs.dangling_sequences = py::reinterpret_borrow<py::object>(pDanglingSequences);
        inputs.bin_index = py::reinterpret_borrow<py::object>(pDictBinIntervalTreeIndex);
        inputs.output_bam_set = py::reinterpret_borrow<py::object>(pOutputBamSet);
        inputs.out_bam_index_buffer = out_bam_index_buffer;
        inputs.counters = &counters;
        ProcessDataPolicy policy(std::move(inputs), pSharedBinIntvalTree, pCoverage, pCoverageIndex, pRow, pCol,
                                 pData, module);
        policy.set_quick_qc(options.quick_qc_mode);
        const std::size_t pairs = std::min(length(pMateBuffer1), length(pMateBuffer2));
        hicx::build_matrix::classify_pairs(options, 0, pairs, policy, counters);

        for (std::size_t s = 0; s < options.restriction_sequences.size(); ++s) {
            if (counters.dangling_end[s] != 0) {
                const py::object key = sequences[s];
                PyObject* sum = PyNumber_Add(get_item(dangling_end, key).ptr(), py::int_(counters.dangling_end[s]).ptr());
                if (sum == nullptr) {
                    throw_python_error();
                }
                dangling_end[key] = py::reinterpret_steal<py::object>(sum);
            }
        }
        py::list counts;
        counts.append(py::int_(counters.one_mate_unmapped));
        counts.append(py::int_(counters.one_mate_low_quality));
        counts.append(py::int_(counters.one_mate_not_unique));
        counts.append(dangling_end);
        for (const std::int64_t value :
             {counters.self_circle, counters.self_ligation, counters.same_fragment, counters.mate_not_close_to_rf,
              counters.count_inward, counters.count_outward, counters.count_left, counters.count_right,
              counters.inter_chromosomal, counters.short_range, counters.long_range, counters.pair_added}) {
            counts.append(py::int_(value));
        }
        counts.append(py::int_(length(pMateBuffer1)));
        counts.append(pResultIndex);
        counts.append(pCounter);
        counts.append(out_bam_index_buffer);
        py::list message;
        message.append(counts);
        pQueueOut.attr("put")(message);
        return;
    } catch (py::error_already_set& error) {
        if (!error.matches(PyExc_Exception)) {
            throw;
        }
        failure = failure_text(error);
    } catch (const hicx::PythonError& error) {
        failure = failure_text(error.type(), error.what());
    } catch (const std::exception& error) {
        failure = failure_text("RuntimeError", error.what());
    }
    pQueueOut.attr("put")(to_str(failure));
}

// --- createMatrix ---------------------------------------------------------------------

std::vector<std::string> str_list(py::handle values, const char* what) {
    std::vector<std::string> out;
    for (py::handle value : py::iter(values)) {
        out.push_back(str_argument(value, "createMatrix", what));
    }
    return out;
}

hicx::ChromSizes chrom_sizes_of_list(py::handle sizes) {
    hicx::ChromSizes out;
    for (py::handle entry : py::iter(sizes)) {
        const std::vector<py::object> fields = unpack(entry, 2);
        out.emplace_back(str_argument(fields[0], "createMatrix", "chromosome names"),
                         integer_argument(fields[1], "createMatrix", "chromosome sizes"));
    }
    return out;
}

py::list chrom_sizes_list(const hicx::ChromSizes& sizes) {
    py::list out;
    for (const auto& [name, size] : sizes) {
        out.append(py::make_tuple(to_str(name), py::int_(size)));
    }
    return out;
}

// The chromosome sizes file as createMatrix reads it: readline() until the
// first line that is empty after strip(), each line split at '\t' into the
// name and int() of the size, kept in an OrderedDict.
hicx::ChromSizes read_chromosome_sizes_file(py::handle name) {
    hicx::ChromSizes sizes;
    for (const std::string& raw : read_text_lines(name)) {
        const std::string_view line = hicx::pytext::strip(raw);
        if (line.empty()) {
            break;
        }
        const std::vector<std::string_view> parts = hicx::pytext::split(line, "\t");
        if (parts.size() < 2) {
            raise_error("IndexError", "list index out of range");
        }
        const std::string chrom(parts[0]);
        const std::int64_t size = to_int64(python_int(to_str(parts[1])));
        const auto found = std::find_if(sizes.begin(), sizes.end(), [&](const auto& e) { return e.first == chrom; });
        if (found != sizes.end()) {
            found->second = size;
        } else {
            sizes.emplace_back(chrom, size);
        }
    }
    return sizes;
}

void create_matrix(py::handle pOutFileName, py::handle pMaxDistance, py::handle pMaxLibraryInsertSize,
                   py::handle pQCfolder, py::handle pThreads, py::handle pDanglingSequence,
                   py::handle pRestrictionSequence, py::handle pSamFiles, py::handle pDoTestRun, py::handle pOutBam,
                   py::handle pChromosomeSizes, py::handle pRestrictionCutFile, py::handle pRegion,
                   py::handle pBinSize, py::handle pInputBufferSize, py::handle pDoTestRunLines,
                   py::handle pSkipDuplicationCheck, py::handle pMinMappingQuality, py::handle pKeepSelfCircles,
                   py::handle pKeepSelfLigation, py::handle pMinDistance, py::handle pGenomeAssembly) {
    namespace bm = hicx::build_matrix;
    bm::Arguments args;
    args.python_api = true;
    args.out_file_name = str_argument(pOutFileName.attr("name"), "createMatrix", "pOutFileName.name");
    if (!pMaxDistance.is_none()) {
        args.max_distance = integer_argument(pMaxDistance, "createMatrix", "pMaxDistance");
    }
    args.max_library_insert_size = integer_argument(pMaxLibraryInsertSize, "createMatrix", "pMaxLibraryInsertSize");
    args.qc_folder = str_argument(pQCfolder, "createMatrix", "pQCfolder");
    args.threads = integer_argument(pThreads, "createMatrix", "pThreads");
    args.dangling_given = truthy(pDanglingSequence);
    if (args.dangling_given) {
        args.dangling_sequences = str_list(pDanglingSequence, "the dangling sequences");
    }
    if (!pRestrictionSequence.is_none()) {
        args.restriction_sequences = str_list(pRestrictionSequence, "the restriction sequences");
    }
    const py::object sam1 = get_index(pSamFiles, 0);
    const py::object sam2 = get_index(pSamFiles, 1);
    args.sam_files = {str_argument(sam1.attr("name"), "createMatrix", "the BAM file names"),
                      str_argument(sam2.attr("name"), "createMatrix", "the BAM file names")};
    args.do_test_run = truthy(pDoTestRun);
    const bool out_bam_given = truthy(pOutBam);
    if (out_bam_given) {
        args.out_bam = str_argument(pOutBam.attr("name"), "createMatrix", "pOutBam.name");
    }
    args.restriction_cut_file_given = !pRestrictionCutFile.is_none();
    if (truthy(pBinSize)) {
        for (py::handle size : py::iter(pBinSize)) {
            args.bin_size.push_back(integer_argument(size, "createMatrix", "the bin sizes"));
        }
    }
    args.input_buffer_size = integer_argument(pInputBufferSize, "createMatrix", "pInputBufferSize");
    args.do_test_run_lines = integer_argument(pDoTestRunLines, "createMatrix", "pDoTestRunLines");
    // readBamFiles tests `pSkipDuplicationCheck is False`
    args.skip_duplication_check = pSkipDuplicationCheck.ptr() != Py_False;
    args.min_mapping_quality = integer_argument(pMinMappingQuality, "createMatrix", "pMinMappingQuality");
    args.keep_self_circles = truthy(pKeepSelfCircles);
    args.keep_self_ligation = truthy(pKeepSelfLigation);
    if (pMinDistance.is_none()) {
        args.min_distance = 0;
        args.min_distance_none = true;
    } else {
        args.min_distance = integer_argument(pMinDistance, "createMatrix", "pMinDistance");
    }
    if (truthy(pGenomeAssembly)) {
        args.genome_assembly = str_argument(pGenomeAssembly, "createMatrix", "pGenomeAssembly");
    }

    bm::Hooks hooks;
    hooks.log = [](const char* level, const std::string& message) { log(kLogger, level, message); };
    hooks.warn = [](const std::string& message) { warn(PyExc_UserWarning, message, 1); };
    hooks.bams_opened = [&] {
        sam1.attr("close")();
        sam2.attr("close")();
        if (!args.do_test_run && out_bam_given) {
            pOutBam.attr("close")();
        }
    };
    if (!pChromosomeSizes.is_none()) {
        hooks.chromosome_sizes = [&] { return read_chromosome_sizes_file(pChromosomeSizes.attr("name")); };
    }
    if (truthy(pRestrictionCutFile)) {
        hooks.restriction_sites = [&](const hicx::ChromSizes& sizes) {
            const py::list sizes_object = chrom_sizes_list(sizes);
            std::vector<bm::RfInterval> sites;
            for (py::handle handle : py::iter(pRestrictionCutFile)) {
                for (py::handle site : bed2interval_list(handle, sizes_object, pRegion)) {
                    const std::vector<py::object> fields = unpack(site, 3);
                    sites.emplace_back(str_argument(fields[0], "createMatrix", "the restriction site chromosomes"),
                                       coordinate(fields[1], "createMatrix"), coordinate(fields[2], "createMatrix"));
                }
            }
            return sites;
        };
    }
    if (truthy(pRegion)) {
        hooks.user_region = [&](const hicx::ChromSizes& sizes) {
            const std::vector<py::object> user = unpack(user_region_function()(chrom_sizes_list(sizes), pRegion), 4);
            return std::make_pair(chrom_sizes_of_list(user[0]),
                                  integer_argument(user[1], "createMatrix", "the region start"));
        };
    }
    hooks.uppercase = [&](std::size_t index, bool dangling_too) {
        const auto upper = [](py::handle list, std::size_t i) {
            set_index(list, static_cast<std::int64_t>(i), get_index(list, static_cast<std::int64_t>(i)).attr("upper")());
        };
        upper(pRestrictionSequence, index);
        if (dangling_too) {
            upper(pDanglingSequence, index);
        }
    };
    hooks.after_qc_log = [&](const std::string& qc_log_file) {
        // QC.main("-l {} -o {}".format(log_file_name, pQCfolder).split())
        py::list argv;
        const std::string command = "-l " + qc_log_file + " -o " + formatted(pQCfolder);
        for (const std::string_view part : hicx::pytext::split_whitespace(command)) {
            argv.append(to_str(part));
        }
        py::module_::import("hicexplorer.hicPrepareQCreport").attr("main")(argv);
    };
    hooks.before_matrix = [&] {
        pOutFileName.attr("close")();
        const py::object name = pOutFileName.attr("name");
        if (::unlink(fs_path(name).c_str()) != 0) {
            raise_os_error(errno, name);
        }
    };

    try {
        bm::run(args, hooks);
    } catch (const bm::Exit& exit) {
        if (exit.logged()) {
            log(kLogger, "error", exit.what());
            raise_system_exit(py::int_(1));
        }
        raise_system_exit(to_str(exit.what()));
    }
}

}  // namespace

void bind_lib_buildMatrixMethods_process(py::module_& m) {
    m.def("process_data", [m](py::handle pMateBuffer1, py::handle pMateBuffer2, py::handle pMinMappingQuality,
                              py::handle pKeepSelfCircles, py::handle pRestrictionSequence,
                              py::handle pKeepSelfLigation, py::handle pMatrixSize, py::handle pRfPositions,
                              py::handle pRefId2name, py::handle pDanglingSequences, py::handle pBinsize,
                              py::handle pResultIndex, py::handle pQueueOut, py::handle pTemplate,
                              py::handle pOutputBamSet, py::handle pCounter, py::handle pSharedBinIntvalTree,
                              py::handle pDictBinIntervalTreeIndex, py::handle pCoverage, py::handle pCoverageIndex,
                              py::handle pOutputFileBufferDir, py::handle pRow, py::handle pCol, py::handle pData,
                              py::handle pMaxInsertSize, py::handle pQuickQCMode) {
        process_data(m, pMateBuffer1, pMateBuffer2, pMinMappingQuality, pKeepSelfCircles, pRestrictionSequence,
                     pKeepSelfLigation, pMatrixSize, pRfPositions, pRefId2name, pDanglingSequences, pBinsize,
                     pResultIndex, pQueueOut, pTemplate, pOutputBamSet, pCounter, pSharedBinIntvalTree,
                     pDictBinIntervalTreeIndex, pCoverage, pCoverageIndex, pOutputFileBufferDir, pRow, pCol, pData,
                     pMaxInsertSize, pQuickQCMode);
    }, py::arg("pMateBuffer1"), py::arg("pMateBuffer2"), py::arg("pMinMappingQuality"), py::arg("pKeepSelfCircles"),
       py::arg("pRestrictionSequence"), py::arg("pKeepSelfLigation"), py::arg("pMatrixSize"), py::arg("pRfPositions"),
       py::arg("pRefId2name"), py::arg("pDanglingSequences"), py::arg("pBinsize"), py::arg("pResultIndex"),
       py::arg("pQueueOut"), py::arg("pTemplate"), py::arg("pOutputBamSet"), py::arg("pCounter"),
       py::arg("pSharedBinIntvalTree"), py::arg("pDictBinIntervalTreeIndex"), py::arg("pCoverage"),
       py::arg("pCoverageIndex"), py::arg("pOutputFileBufferDir"), py::arg("pRow"), py::arg("pCol"), py::arg("pData"),
       py::arg("pMaxInsertSize"), py::arg("pQuickQCMode"));

    m.def("createMatrix", &create_matrix, py::arg("pOutFileName"), py::arg("pMaxDistance"),
          py::arg("pMaxLibraryInsertSize"), py::arg("pQCfolder"), py::arg("pThreads"), py::arg("pDanglingSequence"),
          py::arg("pRestrictionSequence"), py::arg("pSamFiles"), py::arg("pDoTestRun"), py::arg("pOutBam"),
          py::arg("pChromosomeSizes"), py::arg("pRestrictionCutFile"), py::arg("pRegion"), py::arg("pBinSize"),
          py::arg("pInputBufferSize"), py::arg("pDoTestRunLines"), py::arg("pSkipDuplicationCheck"),
          py::arg("pMinMappingQuality"), py::arg("pKeepSelfCircles"), py::arg("pKeepSelfLigation"),
          py::arg("pMinDistance"), py::arg("pGenomeAssembly"));
}

}  // namespace hicx::python

#endif  // HICX_PYTHON_API_LIB_BUILD_MATRIX_METHODS_RUN_HPP
