// hicexplorer._hicx.hicAdjustMatrix: the C++ body of adjustMatrix.
//
// The matrix is loaded by the C++ matrix layer into a hiCMatrix
// (hic_object.hpp); the argument checks, the region file parsing, the log
// records and the exit(1) calls are those of 3.7.6, and the operations
// 3.7.6 asks of the hiCMatrix (reorderChromosomes, maskChromosomes,
// reorderBins, maskBins, getRegionBinRange, getChrBinRange) are asked of it
// through its own methods. The inter and intra chromosomal handling assigns
// zero to matrix blocks as scipy's csr __setitem__ does: stored entries of
// the block become explicit zeros and the missing positions are inserted as
// explicit zeros before the old entries of their row, the rows sorted again
// when they were sorted (_set_many, _insert_many).
//
// Deviations: see hic_object.hpp for the loaded matrix (cool metadata as
// strings); a .hic input, which the hiCMatrix constructor could not read,
// loads through the C++ layer.

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "api.hpp"
#include "hic_object.hpp"
#include "hicx/cool_adapter.hpp"
#include "loops_correct_common.hpp"

namespace hicx::python {

namespace {

constexpr const char* kLogger = "hicexplorer.hicAdjustMatrix";

py::object load_matrix(py::handle path, py::handle chrname_list) {
    std::optional<std::string> chromosome;
    if (!chrname_list.is_none()) {
        chromosome = utf8(lc::item(chrname_list, 0));
    }
    const hic::Loaded loaded = hic::load(fs_path(path), chromosome, std::nullopt);
    return hic::object(loaded, path, chrname_list, true);
}

bool contains(py::handle sequence, py::handle value) {
    const int result = PySequence_Contains(sequence.ptr(), value.ptr());
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

// matrix[row_first:row_last, col_first:col_last] = 0 on a csr matrix, in
// place (hicx::scipy_coo::assign_zero_block).
void zero_block(py::handle matrix_object, std::int64_t row_first, std::int64_t row_last, std::int64_t col_first,
                std::int64_t col_last) {
    Sparse matrix = sparse(matrix_object, "adjustMatrix", "r");
    const auto clamp = [](std::int64_t v, std::int64_t n) { return std::max<std::int64_t>(0, std::min(v, n)); };
    row_first = clamp(row_first, matrix.rows);
    row_last = clamp(row_last, matrix.rows);
    col_first = clamp(col_first, matrix.cols);
    col_last = clamp(col_last, matrix.cols);
    if (row_last <= row_first || col_last <= col_first) {
        return;
    }
    const bool sorted = lc::truthy(matrix_object.attr("has_sorted_indices"));
    dispatch_arithmetic(matrix.data, "adjustMatrix", [&]<class T>() {
        const hicx::scipy_coo::BlockAssignment<T> result = hicx::scipy_coo::assign_zero_block(
            lc::csr_of<T>(matrix), sorted, row_first, row_last, col_first, col_last);
        if (!result.inserted) {
            // only existing cells are affected
            write_back(matrix.data, result.matrix.data);
            return;
        }
        py::module_::import("warnings").attr("warn")(
            "Changing the sparsity structure of a csr_matrix is expensive. lil_matrix is more efficient.",
            py::module_::import("scipy.sparse").attr("SparseEfficiencyWarning"), py::arg("stacklevel") = 3);
        std::vector<std::int32_t> indices32(result.matrix.indices.begin(), result.matrix.indices.end());
        std::vector<std::int32_t> indptr32(result.matrix.indptr.begin(), result.matrix.indptr.end());
        matrix_object.attr("indices") = from_vector(indices32);
        matrix_object.attr("indptr") = from_vector(indptr32);
        matrix_object.attr("data") = from_vector(result.matrix.data);
    });
}

py::object adjust_matrix(py::handle pArgs) {
    const py::object chromosomes = pArgs.attr("chromosomes");
    const py::object regions = pArgs.attr("regions");
    if (!chromosomes.is_none() && !regions.is_none()) {
        log(kLogger, "error", "Please specify either --chromosomes or --regions.");
        lc::python_exit(1);
    }
    const py::module_ module = py::module_::import("hicexplorer.hicAdjustMatrix");
    py::object hic_matrix = py::none();
    const py::object matrix_path = pArgs.attr("matrix");
    auto action_is = [&](const char* name) { return lc::equal(pArgs.attr("action"), py::str(name)); };
    if (lc::truthy(chromosomes)) {
        if (lc::truthy(module.attr("check_cooler")(matrix_path)) && lc::length(chromosomes) == 1 &&
            action_is("keep")) {
            const hicx::CoolFile cool(fs_path(matrix_path));
            py::list names;
            for (const std::string& name : cool.chrom_names()) {
                names.append(to_str(name));
            }
            const py::object first = lc::item(chromosomes, 0);
            if (contains(names, first)) {
                hic_matrix = load_matrix(matrix_path, chromosomes);
            } else {
                log(kLogger, "error",
                    "Chromosome not available in matrix: " + str_of(matrix_path) + " " + str_of(first));
                lc::python_exit(1);
            }
        } else {
            hic_matrix = load_matrix(matrix_path, py::none());
        }
        py::list chromosomes_list = py::list(hic_matrix.attr("chrBinBoundaries"));
        py::list to_operate;
        log(kLogger, "debug", "pArgs.chromosomes " + str_of(chromosomes));
        for (const py::handle chromosome : chromosomes) {
            if (contains(chromosomes_list, chromosome)) {
                to_operate.append(chromosome);
            } else {
                log(kLogger, "warning",
                    "Chromosome not available in matrix: " + str_of(matrix_path) + " " + str_of(chromosome));
            }
        }
        log(kLogger, "debug", "85: chromosomes_list_to_operate_on " + str_of(to_operate));
        if (to_operate.empty()) {
            log(kLogger, "error",
                "No valid chromosome given: " + str_of(chromosomes) + ". Available: " + str_of(chromosomes_list));
            lc::python_exit(1);
        }
        if (action_is("keep")) {
            hic_matrix.attr("reorderChromosomes")(to_operate);
        } else if (action_is("remove")) {
            log(kLogger, "debug", "chromosomes_list " + str_of(chromosomes_list));
            for (const py::handle chromosome : to_operate) {
                if (contains(chromosomes_list, chromosome)) {
                    chromosomes_list.attr("remove")(chromosome);
                    log(kLogger, "debug", "chromosome " + str_of(chromosome));
                }
            }
            log(kLogger, "debug", "chromosomes_list " + str_of(chromosomes_list));
            hic_matrix.attr("reorderChromosomes")(chromosomes_list);
        } else if (action_is("mask")) {
            hic_matrix.attr("maskChromosomes")(to_operate);
        }
    } else if (lc::truthy(regions)) {
        hic_matrix = load_matrix(matrix_path, py::none());
        const py::list chromosomes_list = py::list(hic_matrix.attr("chrBinBoundaries"));
        py::list genomic_regions;
        py::object chrom = py::none();
        bool chrom_bound = false;
        for (const std::string& text : read_text_lines(regions)) {
            const py::str line = to_str(text);
            const py::list fields = line.attr("strip")().attr("split")("\t");
            log(kLogger, "debug", "_line " + str_of(fields));
            if (lc::length(line) < 3) {
                log(kLogger, "warning", "An entry shorter than 3 columns has been found!");
                continue;
            }
            if (fields.size() >= 3) {
                chrom = fields[0];
                chrom_bound = true;
                const py::object start = parse_int(fields[1]);
                const py::object end = parse_int(fields[2]);
                log(kLogger, "debug", "chrom " + str_of(chrom));
                if (contains(chromosomes_list, chrom)) {
                    genomic_regions.append(py::make_tuple(chrom, start, end));
                } else {
                    log(kLogger, "warning",
                        "Chromosome not available in matrix, ignoring regions: " + str_of(matrix_path) + " " +
                            str_of(chrom));
                }
            }
        }
        if (genomic_regions.empty()) {
            log(kLogger, "error", "No valid chromosome given. Available: " + str_of(chromosomes_list));
            lc::python_exit(1);
        }
        py::list matrix_indices_regions;
        for (const py::handle region : genomic_regions) {
            log(kLogger, "debug", "region " + str_of(region));
            const py::object range = hic_matrix.attr("getRegionBinRange")(
                lc::item(region, 0), lc::item(region, 1), lc::item(region, 2));
            if (!range.is_none()) {
                const std::int64_t start = to_int64(lc::item(range, 0));
                const std::int64_t end = to_int64(lc::item(range, 1));
                for (std::int64_t i = start; i < end + 1; ++i) {
                    matrix_indices_regions.append(py::int_(i));
                }
                if (action_is("remove")) {
                    if (!chrom_bound) {
                        raise_error("UnboundLocalError",
                                    "cannot access local variable 'chrom' where it is not associated with a value");
                    }
                    const py::object chr_range = hic_matrix.attr("getChrBinRange")(chrom);
                    const std::int64_t chr_start = to_int64(lc::item(chr_range, 0));
                    const std::int64_t chr_end = to_int64(lc::item(chr_range, 1));
                    if (start > chr_start && end < chr_end - 1) {
                        log(kLogger, "warning",
                            str_of(chrom) + ":" + std::to_string(start) + "-" + std::to_string(end) +
                                " entry may generate discounted regions on a chromosome.Please consider using "
                                "`mask` action to deal with that.");
                    }
                }
            }
        }
        if (action_is("keep")) {
            hic_matrix.attr("reorderBins")(matrix_indices_regions);
        } else if (action_is("mask")) {
            hic_matrix.attr("maskBins")(matrix_indices_regions);
        } else if (action_is("remove")) {
            hic_matrix.attr("maskBins")(matrix_indices_regions);
            hic_matrix.attr("orig_bin_ids") = py::list();
            hic_matrix.attr("orig_cut_intervals") = py::list();
            hic_matrix.attr("nan_bins") = py::list();
        }
    } else if (lc::truthy(pArgs.attr("maskBadRegions"))) {
        if (lc::truthy(module.attr("check_cooler")(matrix_path)) && lc::length(chromosomes) == 1 &&
            action_is("keep")) {
            hic_matrix = load_matrix(matrix_path, chromosomes);
        } else {
            hic_matrix = load_matrix(matrix_path, py::none());
        }
    } else {
        log(kLogger, "info", "No data to adjust given. Please specify either --chromosomes or --region parameter.");
    }
    const py::object handling = pArgs.attr("interIntraHandling");
    if (!handling.is_none()) {
        const bool inter = lc::equal(handling, py::str("inter"));
        const bool intra = !inter && lc::equal(handling, py::str("intra"));
        if (inter || intra) {
            const py::list chromosomes_list = py::list(hic_matrix.attr("chrBinBoundaries"));
            for (const py::handle chromosome : chromosomes_list) {
                const py::object range = hic_matrix.attr("getChrBinRange")(chromosome);
                const std::int64_t start = to_int64(lc::item(range, 0));
                const std::int64_t end = to_int64(lc::item(range, 1));
                const py::object matrix = hic_matrix.attr("matrix");
                if (inter) {
                    zero_block(matrix, start, end, end, std::numeric_limits<std::int64_t>::max());
                } else {
                    zero_block(matrix, start, end, start, end);
                }
            }
        } else {
            log(kLogger, "warning", "Option not valid: " + str_of(handling));
        }
    }
    return hic_matrix;
}

}  // namespace

void bind_hicAdjustMatrix(py::module_& m) {
    m.def("adjustMatrix", &adjust_matrix, py::arg("pArgs"));
}

}  // namespace hicx::python
