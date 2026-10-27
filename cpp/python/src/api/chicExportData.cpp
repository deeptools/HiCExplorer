// hicexplorer._hicx.chicExportData: the C++ body of exportData (main() and
// parse_arguments() run the executable cpp/tools/chicExportData.cpp).
//
// exportData turns the records of an interaction, significant, target,
// aggregate or differential file into the text the tool writes (or the
// columns of a bigwig file) and puts the file names and contents into the
// queue. The records are read with the Viewpoint object's methods as in
// 3.7.6, the target file with hicx::chic::H5Reader; the numbers are written
// with Python's format(), as '{:.{decimal_places}f}'.format(x) wrote them.
//
// A failure puts 'Fail: ' + str(exp) + traceback.format_exc() into the queue
// as in 3.7.6; the traceback names the frames the error passed in this
// implementation, which are not 3.7.6's.

#include <cstdint>
#include <string>
#include <vector>

#include "api.hpp"
#include "chic_common.hpp"
#include "hicx/chic_hdf5.hpp"

namespace hicx::python {

using namespace chic_api;
namespace hc = hicx::chic;

namespace {

// '\t'.join('{:.{decimal_places}f}'.format(x, ...) if isinstance(x, float)
// else str(x) for x in line) + '\n'
std::string format_line(py::handle line, const py::object& spec) {
    std::string out;
    bool first = true;
    for (const py::handle x : line) {
        if (!first) {
            out += "\t";
        }
        first = false;
        out += PyFloat_Check(x.ptr()) ? utf8(format(x, spec)) : utf8(str_object(x));
    }
    return out + "\n";
}

const char* const kInteractionHeader =
    "# Chromosome\tStart\tEnd\tGene\tSum of interactions\tRelative position\tRelative Interactions\tp-value\t"
    "x-fold\tRaw\n";

[[noreturn]] void unbound(const char* name) {
    raise_error("UnboundLocalError", std::string("cannot access local variable '") + name +
                                         "' where it is not associated with a value");
}

}  // namespace

void bind_chicExportData(py::module_& m) {
    m.doc() = "C++ bodies of hicexplorer.chicExportData.";

    m.def("exportData", [](py::object pFileList, py::object pArgs, py::object pViewpointObject,
                           py::object pDecimalPlace, py::object pChromosomeSizes, py::object pBackgroundData,
                           py::object pFileType, py::object pQueue) {
        py::list file_list;
        py::list file_content_list;
        try {
            const bool txt = eq(pArgs.attr("outputFileType"), py::str("txt"));
            const std::string file_ending = txt ? ".txt" : ".bigwig";
            const py::object spec = py::str("." + utf8(format(pDecimalPlace, "")) + "f");
            const py::object file_type_text = str_object(pFileType);
            if (eq(pFileType, py::str("interactions")) || eq(pFileType, py::str("significant"))) {
                // the function's locals: they keep their values from one
                // file to the next unless a file sets them again
                const bool bigwig = eq(pArgs.attr("outputFileType"), py::str("bigwig"));
                py::list chromosome_name;
                py::list start;
                py::list end;
                py::list values;
                py::dict relative_distance;
                bool bigwig_lists = false;
                py::object header;
                py::object chromosome_name_background;
                py::object start_background;
                py::object end_background;
                py::object value_background;
                std::optional<std::string> file_content_string;
                py::object sample;
                for (const py::handle file : pFileList) {
                    if (bigwig) {
                        chromosome_name = py::list();
                        start = py::list();
                        end = py::list();
                        values = py::list();
                        relative_distance = py::dict();
                        bigwig_lists = true;
                    }
                    for (const py::handle sample_handle : file) {
                        sample = py::reinterpret_borrow<py::object>(sample_handle);
                        const py::object data = pViewpointObject.attr("readInteractionFile")(pArgs.attr("file"), sample);
                        const py::object records = getitem(data, 1);
                        const py::list key_list = sorted_keys(records);
                        if (txt) {
                            std::string text = kInteractionHeader;
                            for (const py::handle key : key_list) {
                                text += format_line(getitem(records, key), spec);
                            }
                            file_content_string = text;
                        } else {
                            if (!bigwig_lists) {
                                unbound("chromosome_name");
                            }
                            const py::object value_kind = pArgs.attr("outputValueBigwig");
                            for (const py::handle key : key_list) {
                                const py::object record = getitem(records, key);
                                chromosome_name.append(str_object(getitem(record, 0)));
                                start.append(int_of(getitem(record, 1)));
                                end.append(int_of(getitem(record, 2)));
                                if (eq(value_kind, py::str("relative-interactions"))) {
                                    values.append(float_of(getitem(record, 6)));
                                } else if (eq(value_kind, py::str("p-value"))) {
                                    values.append(float_of(getitem(record, 7)));
                                } else if (eq(value_kind, py::str("x-fold"))) {
                                    values.append(float_of(getitem(record, 8)));
                                } else if (eq(value_kind, py::str("raw"))) {
                                    values.append(float_of(getitem(record, 9)));
                                }
                                py::list coordinates;
                                coordinates.append(str_object(getitem(record, 0)));
                                coordinates.append(int_of(getitem(record, 1)));
                                coordinates.append(int_of(getitem(record, 2)));
                                relative_distance[getitem(record, 5)] = coordinates;
                            }
                            const py::object first_chromosome = getitem(chromosome_name, 0);
                            py::list header_list;
                            header_list.append(py::make_tuple(first_chromosome, getitem(pChromosomeSizes, first_chromosome)));
                            header = header_list;
                            if (!pArgs.attr("backgroundModelFile").is_none()) {
                                py::list cb;
                                py::list sb;
                                py::list eb;
                                py::list vb;
                                for (const py::handle key : sorted_keys(pBackgroundData)) {
                                    if (contains(relative_distance, key)) {
                                        const py::object coordinates = getitem(relative_distance, key);
                                        cb.append(getitem(coordinates, 0));
                                        sb.append(getitem(coordinates, 1));
                                        eb.append(getitem(coordinates, 2));
                                        vb.append(float_of(getitem(getitem(pBackgroundData, key), 0)));
                                    }
                                }
                                chromosome_name_background = cb;
                                start_background = sb;
                                end_background = eb;
                                value_background = vb;
                            }
                        }
                    }
                    if (txt) {
                        if (!file_content_string) {
                            unbound("file_content_string");
                        }
                        file_content_list.append(to_str(*file_content_string));
                        if (!sample) {
                            unbound("sample");
                        }
                        file_list.append(add(add(join("_", sample), py::str("_")), add(file_type_text, py::str(file_ending))));
                    } else {
                        if (!sample) {
                            unbound("sample");
                        }
                        if (!header) {
                            unbound("header");
                        }
                        py::list columns;
                        columns.append(header);
                        columns.append(chromosome_name);
                        columns.append(start);
                        columns.append(end);
                        columns.append(values);
                        const py::object joined = join("_", sample);
                        py::list content;
                        content.append(columns);
                        py::list names;
                        if (!pArgs.attr("backgroundModelFile").is_none()) {
                            if (!chromosome_name_background) {
                                unbound("chromosome_name_background");
                            }
                            py::list background_columns;
                            background_columns.append(header);
                            background_columns.append(chromosome_name_background);
                            background_columns.append(start_background);
                            background_columns.append(end_background);
                            background_columns.append(value_background);
                            content.append(background_columns);
                            names.append(add(joined, py::str(file_ending)));
                            names.append(add(add(add(py::str("background_"), joined), py::str("_")),
                                             add(file_type_text, py::str(file_ending))));
                        } else {
                            names.append(add(add(joined, py::str("_")), add(file_type_text, py::str(file_ending))));
                        }
                        file_content_list.append(content);
                        file_list.append(names);
                    }
                }
            } else if (eq(pFileType, py::str("target"))) {
                for (const py::handle target_file : pFileList) {
                    hc::H5Reader reader = open_h5(pArgs.attr("file"));
                    const std::string path = utf8(join("/", target_file));
                    if (!reader.contains(path)) {
                        raise_error("KeyError", reader.missing_message(path));
                    }
                    if (!reader.contains(path + "/chromosome")) {
                        raise_error("TypeError", "'NoneType' object is not subscriptable");
                    }
                    const py::object chromosome = h5_scalar(reader.read(path + "/chromosome"));
                    const hc::H5Value starts = reader.read(path + "/start_list");
                    const hc::H5Value ends = reader.read(path + "/end_list");
                    reader.close();
                    std::string text;
                    const std::size_t n = std::min(starts.length, ends.length);
                    for (std::size_t i = 0; i < n; ++i) {
                        text += utf8(decode_utf8(chromosome)) + "\t" + utf8(decode_utf8(h5_element(starts, i))) + "\t" +
                                utf8(decode_utf8(h5_element(ends, i))) + "\n";
                    }
                    file_content_list.append(to_str(text));
                    file_list.append(add(add(join("_", target_file), py::str("_")), add(file_type_text, py::str(".txt"))));
                }
            } else if (eq(pFileType, py::str("aggregate"))) {
                const std::string header =
                    "# Chromosome\tStart\tEnd\tGene\tSum of interactions\tRelative position\tRaw\n";
                for (const py::handle file : pFileList) {
                    for (const py::handle sample : file) {
                        const std::vector<py::object> read =
                            unpack(pViewpointObject.attr("readAggregatedFileHDF")(pArgs.attr("file"), sample), 2);
                        std::string text = header;
                        for (const py::handle line : read[0]) {
                            text += format_line(line, spec);
                        }
                        file_content_list.append(to_str(text));
                        file_list.append(add(add(join("_", sample), py::str("_")), add(file_type_text, py::str(".txt"))));
                    }
                }
            } else if (eq(pFileType, py::str("differential"))) {
                const std::string header =
                    "# Chromosome\tStart\tEnd\tGene\tRelative distance\tsum of interactions 1\ttarget_1 raw\tsum of "
                    "interactions 2\ttarget_2 raw\tp-value\n";
                const char* const classification[] = {"accepted", "all", "rejected"};
                for (const py::handle file : pFileList) {
                    const py::object line_content = pViewpointObject.attr("readDifferentialFile")(pArgs.attr("file"), file);
                    std::size_t i = 0;
                    for (const py::handle item : line_content) {
                        std::string text = header;
                        for (const py::handle line : item) {
                            text += format_line(line, spec);
                        }
                        file_content_list.append(to_str(text));
                        if (i >= 3) {
                            raise_error("IndexError", "list index out of range");
                        }
                        file_list.append(add(add(add(join("_", file), py::str("_")), py::str(classification[i])),
                                             add(py::str("_"), add(file_type_text, py::str(".txt")))));
                        ++i;
                    }
                }
            }
        } catch (py::error_already_set& error) {
            if (!is_exception(error)) {
                throw;
            }
            pQueue.attr("put")(to_str("Fail: " + str_of(error.value()) + traceback_text(error)));
            return;
        } catch (const std::exception& error) {
            pQueue.attr("put")(to_str("Fail: " + fail_text(error) + traceback_text(error)));
            return;
        }
        py::list result;
        result.append(file_list);
        result.append(file_content_list);
        pQueue.attr("put")(result);
    }, py::arg("pFileList"), py::arg("pArgs"), py::arg("pViewpointObject"), py::arg("pDecimalPlace"),
       py::arg("pChromosomeSizes"), py::arg("pBackgroundData"), py::arg("pFileType"), py::arg("pQueue"));
}

}  // namespace hicx::python
