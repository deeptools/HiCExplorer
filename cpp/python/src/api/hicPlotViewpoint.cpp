// hicexplorer._hicx.hicPlotViewpoint: the C++ bodies of relabelTicks and
// getViewpointValues (main() and parse_arguments() run the executable
// cpp/tools/hicPlotViewpoint.cpp).
//
// getViewpointValues loads the matrix as hm.hiCMatrix(pMatrix) does
// (hicx::ToolMatrix, the matrix layer the tool uses), keeps the chromosomes
// keepOnlyTheseChr keeps, finds the bins of the reference point and of the
// region as getRegionBinRange does (the bin containing a position, first in
// (start, end, index) order) and sums the reference point's rows in float64.
// The bins it returns have hiCMatrix's types: str and int for a cool file,
// numpy.str_ and numpy.int64 for an h5 file, and the matrix values the
// matrix's numpy dtype.

#include <cmath>
#include <cstdint>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include "api.hpp"
#include "chic_common.hpp"
#include "hicx/tool_matrix.hpp"

namespace hicx::python {

using namespace chic_api;

namespace {

constexpr const char* kLogger = "hicexplorer.hicPlotViewpoint";

// The bins of the matrix after keepOnlyTheseChr, by their original index.
struct View {
    std::vector<std::int64_t> original;
    std::vector<hicx::CutInterval> intervals;
};

View keep_chromosomes(const hicx::ToolMatrix& matrix, py::handle chromosomes) {
    View view;
    const std::vector<hicx::CutInterval>& all = matrix.cut_intervals();
    if (chromosomes.is_none()) {
        view.intervals = all;
        for (std::size_t i = 0; i < all.size(); ++i) {
            view.original.push_back(static_cast<std::int64_t>(i));
        }
        return view;
    }
    const auto& boundaries = matrix.boundaries();
    // [self.chrBinBoundaries[x] for x in chromosome_list]
    for (const py::handle name : chromosomes) {
        const std::string text = utf8(str_object(name));
        bool known = false;
        for (const auto& entry : boundaries) {
            known = known || entry.first == text;
        }
        if (!known) {
            raise_error("ValueError", "Chromosome name not in matrix. " + repr_of(name));
        }
    }
    // `chrName in chromosome_list` (a substring test for a str)
    std::vector<char> selected(all.size(), 0);
    for (const auto& entry : boundaries) {
        if (!contains(chromosomes, to_str(entry.first))) {
            continue;
        }
        for (std::int64_t i = entry.second.first; i < entry.second.last; ++i) {
            selected[static_cast<std::size_t>(i)] = 1;
        }
    }
    for (std::size_t i = 0; i < all.size(); ++i) {
        if (selected[i]) {
            view.original.push_back(static_cast<std::int64_t>(i));
            view.intervals.push_back(all[i]);
        }
    }
    return view;
}

// sorted(interval_trees[chrom][position:position + 1])[0].data
std::optional<std::int64_t> bin_containing(const View& view, const std::string& chrom, std::int64_t position) {
    std::optional<std::tuple<std::int64_t, std::int64_t, std::int64_t>> best;
    for (std::size_t i = 0; i < view.intervals.size(); ++i) {
        const hicx::CutInterval& bin = view.intervals[i];
        if (bin.chrom == chrom && bin.start <= position && position < bin.end) {
            const auto key = std::make_tuple(bin.start, bin.end, static_cast<std::int64_t>(i));
            if (!best || key < *best) {
                best = key;
            }
        }
    }
    if (!best) {
        return std::nullopt;
    }
    return std::get<2>(*best);
}

// hic.getRegionBinRange(chrname, startpos, endpos): (start bin, end bin) or None
py::object region_bin_range(const View& view, py::handle chromosome, py::handle start, py::handle end) {
    const std::string chrom = utf8(str_object(chromosome));
    bool known = false;
    std::vector<std::string> names;
    for (const hicx::CutInterval& bin : view.intervals) {
        if (names.empty() || names.back() != bin.chrom) {
            bool seen = false;
            for (const std::string& name : names) {
                seen = seen || name == bin.chrom;
            }
            if (!seen) {
                names.push_back(bin.chrom);
            }
        }
        known = known || bin.chrom == chrom;
    }
    if (!known) {
        std::string keys;
        for (std::size_t i = 0; i < names.size(); ++i) {
            keys += (i > 0 ? ", " : "") + repr_of(to_str(names[i]));
        }
        raise_error("ValueError", "chromosome: " + chrom + " name not found in matrixvalid names are: dict_keys([" +
                                      keys + "])");
    }
    const std::optional<std::int64_t> first = bin_containing(view, chrom, to_int64(int_of(start)));
    const std::optional<std::int64_t> last = bin_containing(view, chrom, to_int64(int_of(end)));
    if (!first || !last) {
        return py::none();
    }
    return py::make_tuple(py::int_(*first), py::int_(*last));
}

// hic.getBinPos(index) with hiCMatrix's field types
py::tuple bin_position(const View& view, std::int64_t index, bool h5) {
    if (index < 0 || index >= static_cast<std::int64_t>(view.intervals.size())) {
        raise_error("ValueError", "binIndex: " + std::to_string(index) + " not found");
    }
    const hicx::CutInterval& bin = view.intervals[static_cast<std::size_t>(index)];
    if (h5) {
        return py::make_tuple(numpy().attr("str_")(to_str(bin.chrom)), np_int64(bin.start), np_int64(bin.end),
                              bin.extra_text.empty() ? np_float64(bin.extra) : py::object(numpy().attr("str_")(to_str(bin.extra_text))));
    }
    return py::make_tuple(to_str(bin.chrom), py::int_(bin.start), py::int_(bin.end), py::float_(bin.extra));
}

// matrix[row, column] as the numpy scalar of the matrix dtype
py::object matrix_value(double value, const std::string& dtype) {
    if (dtype == "int64") {
        return numpy_scalar<std::int64_t>(static_cast<std::int64_t>(value));
    }
    if (dtype == "int32") {
        return numpy_scalar<std::int32_t>(static_cast<std::int32_t>(value));
    }
    if (dtype == "float32") {
        return numpy_scalar<float>(static_cast<float>(value));
    }
    return numpy_scalar<double>(value);
}

}  // namespace

void bind_hicPlotViewpoint(py::module_& m) {
    m.doc() = "C++ bodies of hicexplorer.hicPlotViewpoint.";

    m.def("relabelTicks", [](py::handle pTick) {
        // "{:.2f} Kb".format(int(pTick) / 1e3), Mb from 1e6 on
        if (lt(pTick, py::float_(1e6))) {
            return add(format(truediv(int_of(pTick), py::float_(1e3)), ".2f"), py::str(" Kb"));
        }
        return add(format(truediv(int_of(pTick), py::float_(1e6)), ".2f"), py::str(" Mb"));
    }, py::arg("pTick"));

    m.def("getViewpointValues", [](py::handle pMatrix, py::handle pReferencePoint, py::handle pChromViewpoint,
                                   py::handle pRegion_start, py::handle pRegion_end, py::handle pInteractionList,
                                   py::handle pChromosome) {
        const hicx::ToolMatrix matrix = hicx::ToolMatrix::load(fs_path(pMatrix));
        const bool h5 = matrix.input_is_h5();
        const View view = keep_chromosomes(matrix, pChromosome);
        const py::ssize_t n = length(pReferencePoint);
        py::object reference_range;
        if (n == 2) {
            reference_range = region_bin_range(view, getitem(pReferencePoint, 0), int_of(getitem(pReferencePoint, 1)),
                                               int_of(getitem(pReferencePoint, 1)));
        } else if (n == 3) {
            reference_range = region_bin_range(view, getitem(pReferencePoint, 0), int_of(getitem(pReferencePoint, 1)),
                                               int_of(getitem(pReferencePoint, 2)));
        } else {
            log(kLogger, "error", "No valid reference point given. " + utf8(format(pReferencePoint, "")));
            python_exit(1);
        }
        const std::vector<py::object> reference = unpack(reference_range, 2);
        const py::object view_point_range = region_bin_range(view, pChromViewpoint, pRegion_start, pRegion_end);
        if (view_point_range.is_none()) {
            raise_error("TypeError", "'NoneType' object is not subscriptable");
        }
        const std::int64_t range_start = to_int64(getitem(view_point_range, 0));
        const std::int64_t elements = to_int64(getitem(view_point_range, 1)) - range_start;
        if (elements < 0) {
            raise_error("ValueError", "negative dimensions are not allowed");
        }
        std::vector<double> data(static_cast<std::size_t>(elements), 0.0);
        py::object interactions = py::none();
        py::list interactions_list;
        const bool with_interactions = !pInteractionList.is_none();
        if (with_interactions) {
            interactions = interactions_list;
        }
        const hicx::CsrMatrix& csr = matrix.matrix();
        const std::string dtype = csr.dtype();
        for (std::int64_t vp = to_int64(reference[0]); vp <= to_int64(reference[1]); ++vp) {
            const py::tuple first = bin_position(view, vp, h5);
            const std::int64_t row = view.original[static_cast<std::size_t>(vp)];
            for (std::int64_t j = 0; j < elements; ++j) {
                const std::int64_t idx = range_start + j;
                const double value = csr.at(row, view.original[static_cast<std::size_t>(idx)]);
                data[static_cast<std::size_t>(j)] += value;
                if (with_interactions) {
                    const py::tuple second = bin_position(view, idx, h5);
                    interactions_list.append(py::make_tuple(first[0], first[1], first[2], second[0], second[1],
                                                            second[2], matrix_value(value, dtype)));
                }
            }
        }
        py::list out;
        out.append(reference[0]);
        out.append(reference[1]);
        out.append(view_point_range);
        out.append(from_vector(data));
        out.append(interactions);
        return out;
    }, py::arg("pMatrix"), py::arg("pReferencePoint"), py::arg("pChromViewpoint"), py::arg("pRegion_start"),
       py::arg("pRegion_end"), py::arg("pInteractionList") = py::none(), py::arg("pChromosome") = py::none());
}

}  // namespace hicx::python
