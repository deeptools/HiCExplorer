// hicexplorer._hicx.hicValidateLocations: the C++ bodies of the public
// functions of hicexplorer.hicValidateLocations.
//
// The functions take and return pandas objects, as in 3.7.6. A DataFrame
// argument is read column by column (its numpy values) into a
// hicx::TextTable, the computation is hicx::validate_locations (the code the
// hicValidateLocations tool runs: read_csv's dtype inference, pybedtools'
// sort, merge and intersect -c, the binning, the cool lookup), and the result
// is built as the DataFrame or Series 3.7.6 returned. The debug records show
// the str() of those pandas objects, as 3.7.6 logged them.
//
// Deviations from 3.7.6, recorded in hicexplorer/test/api/cases_hicValidateLocations.py:
//  * a DataFrame column that is neither int, float, bool nor object of str,
//    coordinates that are not numbers, a bin size that is not an int, and a
//    frame whose columns 1 and 2 are not at positions 1 and 2 raise TypeError
//    naming HiCExplorer 4;
//  * the errors pybedtools raises when bedtools rejects its input (unsorted
//    input to merge, non-integer coordinates) and the errors of a cool file
//    cooler cannot open are RuntimeError with the C++ message.

#include <optional>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicx/text_table.hpp"
#include "hicx/validate_locations.hpp"

namespace hicx::python {

namespace vl = hicx::validate_locations;

namespace {

constexpr const char* kLogger = "hicexplorer.hicValidateLocations";

py::module_ pandas() { return py::module_::import("pandas"); }

bool truthy(py::handle value) {
    const int result = PyObject_IsTrue(value.ptr());
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

bool equals_text(py::handle value, const char* text) {
    return is_str(value) && utf8(value) == text;
}

// pd.read_csv(file, sep='\t', header=None) of a path or a file object.
hicx::TextTable read_table(py::handle file) {
    std::string text;
    if (py::hasattr(file, "read")) {
        const py::object content = file.attr("read")();
        if (PyBytes_Check(content.ptr())) {
            text.assign(PyBytes_AS_STRING(content.ptr()), static_cast<std::size_t>(PyBytes_GET_SIZE(content.ptr())));
        } else {
            text = utf8(content);
        }
    } else {
        text = read_file(file);
    }
    hicx::TextTable table;
    try {
        table = hicx::TextTable::parse_tsv(text);
    } catch (const std::runtime_error& error) {
        raise_error("pandas.errors.ParserError", std::string(error.what()) + "\n");
    }
    if (table.cols() == 0) {
        raise_error("pandas.errors.EmptyDataError", "No columns to parse from file");
    }
    return table;
}

py::object column_object(const hicx::TableColumn& column) {
    switch (column.type) {
        case hicx::ColumnType::Int64: return from_vector(column.ints);
        case hicx::ColumnType::Float64: return from_vector(column.floats);
        default: {
            py::list values;
            for (const std::string& value : column.strings) {
                values.append(to_str(value));
            }
            return py::module_::import("numpy").attr("array")(values, py::arg("dtype") = py::module_::import("builtins").attr("object"));
        }
    }
}

// The DataFrame of a table with the given column labels (0, 1, ... when
// empty) and the default RangeIndex.
py::object frame_of(const hicx::TextTable& table, const std::vector<py::object>& labels = {}) {
    py::dict data;
    for (std::size_t c = 0; c < table.cols(); ++c) {
        const py::object label = labels.empty() ? py::object(py::int_(c)) : labels[c];
        data[label] = column_object(table.column(c));
    }
    return pandas().attr("DataFrame")(data);
}

std::vector<py::object> bed_names(std::size_t count) {
    // pybedtools' to_dataframe names the columns of a BED file
    static const char* const names[] = {"chrom", "start", "end", "name", "score", "strand",
                                        "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"};
    if (count > 12) {
        raise_error("TypeError", "hicValidateLocations: HiCExplorer 4 reads BED records of at most 12 columns");
    }
    std::vector<py::object> out;
    for (std::size_t c = 0; c < count; ++c) {
        out.push_back(py::str(names[c]));
    }
    return out;
}

py::list column_labels(py::handle frame) {
    return py::list(frame.attr("columns"));
}

py::object get_column(py::handle frame, py::handle label) {
    PyObject* column = PyObject_GetItem(frame.ptr(), label.ptr());
    if (column == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(column);
}

// df[[labels...]]: pandas' KeyError when some are missing.
void require_labels(py::handle frame, const std::vector<int>& wanted) {
    const py::list labels = column_labels(frame);
    py::list missing;
    for (const int label : wanted) {
        bool found = false;
        for (py::handle have : labels) {
            if (PyObject_RichCompareBool(have.ptr(), py::int_(label).ptr(), Py_EQ) == 1) {
                found = true;
                break;
            }
        }
        if (!found) {
            missing.append(py::int_(label));
        }
    }
    if (py::len(missing) > 0) {
        raise_error("KeyError", str_of(missing) + " not in index");
    }
}

// One column of a DataFrame as a TableColumn, from its numpy values.
hicx::TableColumn table_column(py::handle series, const std::string& function) {
    const py::array values = as_ndarray(series.attr("to_numpy")());
    const std::string kind = py::str(values.dtype().attr("kind"));
    hicx::TableColumn column;
    if (kind == "i" || kind == "u") {
        column.type = hicx::ColumnType::Int64;
        const py::array cast = values.attr("astype")("int64");
        column.ints = to_vector<std::int64_t>(cast);
    } else if (kind == "f") {
        column.type = hicx::ColumnType::Float64;
        const py::array cast = values.attr("astype")("float64");
        column.floats = to_vector<double>(cast);
    } else if (kind == "b") {
        column.type = hicx::ColumnType::String;
        for (py::handle value : values.attr("tolist")()) {
            column.strings.push_back(truthy(value) ? "True" : "False");
        }
    } else if (kind == "O") {
        column.type = hicx::ColumnType::String;
        for (py::handle value : values.attr("tolist")()) {
            if (!is_str(value)) {
                raise_error("TypeError", function + ": HiCExplorer 4 takes object columns of str, not " +
                                             type_name(value));
            }
            column.strings.push_back(utf8(value));
        }
    } else {
        raise_error("TypeError", function + ": HiCExplorer 4 takes int, float, bool and str columns, not " +
                                     dtype_name(values));
    }
    return column;
}

// The columns `labels` (all when empty) of a DataFrame as a table.
hicx::TextTable table_of_frame(py::handle frame, const std::string& function, const std::vector<int>& labels = {}) {
    std::vector<hicx::TableColumn> columns;
    if (labels.empty()) {
        for (py::handle label : column_labels(frame)) {
            columns.push_back(table_column(get_column(frame, label), function));
        }
    } else {
        require_labels(frame, labels);
        for (const int label : labels) {
            columns.push_back(table_column(get_column(frame, py::int_(label)), function));
        }
    }
    return hicx::TextTable::with_columns(std::move(columns));
}

void require_numeric_coordinates(const hicx::TextTable& table, const std::vector<std::size_t>& columns,
                                 const std::string& function) {
    for (const std::size_t c : columns) {
        if (c < table.cols() && table.column(c).type != hicx::ColumnType::Int64) {
            raise_error("TypeError", function + ": HiCExplorer 4 takes integer coordinates for bedtools");
        }
    }
}

py::object bool_series(const std::vector<bool>& values, py::handle name) {
    py::list items;
    for (const bool value : values) {
        items.append(py::bool_(value));
    }
    py::object array = py::module_::import("numpy").attr("array")(items, py::arg("dtype") = "bool");
    return pandas().attr("Series")(array, py::arg("name") = name);
}

// x['name'] >= 1 of an intersect -c result (column 3 of the read-back table)
std::vector<bool> name_at_least_one(const hicx::TextTable& x) {
    if (x.cols() < 4) {
        raise_error("KeyError", "'name'");
    }
    const hicx::TableColumn& column = x.column(3);
    std::vector<bool> mask(x.rows());
    for (std::size_t r = 0; r < x.rows(); ++r) {
        switch (column.type) {
            case hicx::ColumnType::Int64: mask[r] = column.ints[r] >= 1; break;
            case hicx::ColumnType::Float64: mask[r] = column.floats[r] >= 1.0; break;
            default: raise_error("TypeError", "'>=' not supported between instances of 'str' and 'int'");
        }
    }
    return mask;
}

void apply_prefix(hicx::TextTable& table, std::size_t column, py::handle mode) {
    if (equals_text(mode, "add")) {
        vl::apply_chr_prefix(table, column, "add");
    } else if (equals_text(mode, "remove")) {
        if (table.column(column).type != hicx::ColumnType::String) {
            raise_error("AttributeError", "Can only use .str accessor with string values!");
        }
        vl::apply_chr_prefix(table, column, "remove");
    }
}

hicx::TextTable empty_bed_check(const hicx::TextTable& table) {
    if (table.rows() == 0) {
        // pybedtools' to_dataframe of an empty result file
        raise_error("pandas.errors.EmptyDataError", "No columns to parse from file");
    }
    return table;
}

}  // namespace

void bind_hicValidateLocations(py::module_& m) {
    m.doc() = "C++ bodies of hicexplorer.hicValidateLocations.";

    m.def("readThreeColumn", [](py::handle pFile, py::handle pChrPrefixProtein) {
        const hicx::TextTable full = read_table(pFile);
        if (full.cols() < 3) {
            py::list missing;
            for (std::size_t c = full.cols(); c < 3; ++c) {
                missing.append(py::int_(c));
            }
            raise_error("KeyError", str_of(missing) + " not in index");
        }
        hicx::TextTable table = full.select_columns({0, 1, 2});
        log(kLogger, "debug", "protein_df " + str_of(frame_of(table)));
        apply_prefix(table, 0, pChrPrefixProtein);
        return frame_of(table);
    }, py::arg("pFile"), py::arg("pChrPrefixProtein"));

    m.def("readLoopFile", [](py::handle pInputFile, py::handle pChrPrefixLoops) {
        hicx::TextTable table = read_table(pInputFile);
        if (equals_text(pChrPrefixLoops, "add") || equals_text(pChrPrefixLoops, "remove")) {
            apply_prefix(table, 0, pChrPrefixLoops);
            if (table.cols() < 4) {
                raise_error("KeyError", "3");
            }
            apply_prefix(table, 3, pChrPrefixLoops);
        }
        return frame_of(table);
    }, py::arg("pInputFile"), py::arg("pChrPrefixLoops"));

    m.def("overlapLoop", [](py::handle pDataFrameLoop, py::handle pDataFrameProtein) {
        const hicx::TextTable x_side = table_of_frame(pDataFrameLoop, "overlapLoop", {0, 1, 2});
        const hicx::TextTable y_side = table_of_frame(pDataFrameLoop, "overlapLoop", {3, 4, 5});
        const hicx::TextTable protein = table_of_frame(pDataFrameProtein, "overlapLoop");
        require_numeric_coordinates(x_side, {1, 2}, "overlapLoop");
        require_numeric_coordinates(y_side, {1, 2}, "overlapLoop");
        require_numeric_coordinates(protein, {1, 2}, "overlapLoop");
        const std::vector<bool> mask_x = name_at_least_one(empty_bed_check(vl::intersect_count_table(x_side, protein)));
        const std::vector<bool> mask_y = name_at_least_one(empty_bed_check(vl::intersect_count_table(y_side, protein)));
        std::vector<bool> selection(mask_x.size());
        for (std::size_t r = 0; r < selection.size(); ++r) {
            selection[r] = mask_x[r] && mask_y[r];
        }
        return bool_series(selection, py::str("name"));
    }, py::arg("pDataFrameLoop"), py::arg("pDataFrameProtein"));

    m.def("overlapTAD", [](py::handle pDataFrameTAD, py::handle pDataFrameProtein) {
        const hicx::TextTable tads = table_of_frame(pDataFrameTAD, "overlapTAD");
        require_numeric_coordinates(tads, {1, 2}, "overlapTAD");
        log(kLogger, "debug", "tad_bedtool_x " + tads.to_tsv("."));
        const hicx::TextTable protein = table_of_frame(pDataFrameProtein, "overlapTAD");
        require_numeric_coordinates(protein, {1, 2}, "overlapTAD");
        const std::vector<bool> mask = name_at_least_one(empty_bed_check(vl::intersect_count_table(tads, protein)));
        return bool_series(mask, py::str("name"));
    }, py::arg("pDataFrameTAD"), py::arg("pDataFrameProtein"));

    m.def("correlateCool", [](py::handle pCoolFile, py::handle pDataFrameLoop) {
        const std::string path = fs_path(pCoolFile);
        const hicx::TextTable loops = table_of_frame(pDataFrameLoop, "correlateCool");
        vl::CoolPeaks peaks;
        std::vector<bool> matched(loops.rows(), false);
        if (loops.cols() >= 6) {
            for (const std::size_t c : {1, 2, 4, 5}) {
                if (loops.column(c).type == hicx::ColumnType::String) {
                    raise_error("TypeError", "correlateCool: HiCExplorer 4 takes numeric loop coordinates");
                }
            }
            peaks = vl::correlate_cool(path, loops);
            // str() of a float coordinate has a '.', which cooler's region
            // parser rejects, so those loops do not match
            for (std::size_t r = 0; r < loops.rows(); ++r) {
                bool numeric = true;
                for (const std::size_t c : {1, 2, 4, 5}) {
                    numeric = numeric && loops.column(c).type == hicx::ColumnType::Int64;
                }
                matched[r] = numeric && peaks.matched[r];
            }
        } else {
            // line[4] .. line[6] of itertuples raise inside the try block
            const hicx::TextTable no_loops = hicx::TextTable::with_columns({});
            peaks = vl::correlate_cool(path, no_loops);
        }
        const py::object selection = bool_series(matched, py::none());
        log(kLogger, "debug", str_of(selection));
        return py::make_tuple(numpy_scalar<std::int64_t>(peaks.number_of_peaks), selection);
    }, py::arg("pCoolFile"), py::arg("pDataFrameLoop"));

    m.def("applyBinning", [](py::handle pDataFrame, py::handle pBinSize, py::handle pMerge) {
        const py::list labels = column_labels(pDataFrame);
        if (py::len(labels) < 3 || PyObject_RichCompareBool(labels[1].ptr(), py::int_(1).ptr(), Py_EQ) != 1 ||
            PyObject_RichCompareBool(labels[2].ptr(), py::int_(2).ptr(), Py_EQ) != 1) {
            require_labels(pDataFrame, {1, 2});
            raise_error("TypeError", "applyBinning: HiCExplorer 4 takes the columns 1 and 2 at the positions 1 and 2");
        }
        if (!PyLong_Check(pBinSize.ptr())) {
            raise_error("TypeError", "applyBinning: HiCExplorer 4 takes pBinSize as an int, not " + type_name(pBinSize));
        }
        const hicx::TextTable table = table_of_frame(pDataFrame, "applyBinning");
        for (const std::size_t c : {1, 2}) {
            if (table.column(c).type == hicx::ColumnType::String) {
                raise_error("TypeError", "unsupported operand type(s) for /: 'str' and 'int'");
            }
            if (table.column(c).type == hicx::ColumnType::Float64) {
                for (const double value : table.column(c).floats) {
                    if (!std::isfinite(value)) {
                        raise_error("pandas.errors.IntCastingNaNError",
                                    "Cannot convert non-finite values (NA or inf) to integer");
                    }
                }
            }
        }
        const bool merge = truthy(pMerge);
        const vl::BinningSteps steps = vl::binning_steps(table, to_int64(pBinSize), merge);

        // out = pDataFrame.copy() with the columns 1 and 2 replaced, then
        // out.drop_duplicates(): the other columns keep their values and
        // dtypes, the index keeps the labels of the rows kept
        const py::object copy = pDataFrame.attr("copy")();
        copy.attr("__setitem__")(py::int_(1), from_vector(steps.binned.column(1).ints));
        copy.attr("__setitem__")(py::int_(2), from_vector(steps.binned.column(2).ints));
        log(kLogger, "debug", "pDataFrame_out " + str_of(copy));
        py::object out = copy;
        if (steps.kept.size() != table.rows()) {
            std::vector<std::int64_t> kept(steps.kept.begin(), steps.kept.end());
            out = copy.attr("take")(from_vector(kept));
        }
        log(kLogger, "debug", "pDataFrame_out " + str_of(out));
        if (!merge) {
            return out;
        }
        const hicx::TextTable& merged = empty_bed_check(steps.merged);
        log(kLogger, "debug", "bedtools_data " + str_of(frame_of(merged, bed_names(merged.cols()))));
        log(kLogger, "debug", "bedtools_data " + str_of(frame_of(steps.sorted, bed_names(steps.sorted.cols()))));
        return frame_of(steps.sorted, bed_names(steps.sorted.cols()));
    }, py::arg("pDataFrame"), py::arg("pBinSize"), py::arg("pMerge") = true);

    m.def("writeLoopFile", [](py::handle pOutFileName, py::handle pLoopDataFrame) {
        const std::string text = table_of_frame(pLoopDataFrame, "writeLoopFile").to_tsv("");
        if (py::hasattr(pOutFileName, "write")) {
            pOutFileName.attr("write")(to_str(text));
            return;
        }
        const py::object handle = py::module_::import("builtins").attr("open")(pOutFileName, "w");
        handle.attr("write")(to_str(text));
        handle.attr("close")();
    }, py::arg("pOutFileName"), py::arg("pLoopDataFrame"));
}

}  // namespace hicx::python
