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
// Fixed in HiCExplorer 4, as in the command line tool (cpp/tools/
// hicAdjustMatrix.cpp): --maskBadRegions masks the bins of its BED file; a
// chromosome name that is not in the matrix, and a call with none of
// chromosomes, regions and maskBadRegions, exit(1) with an error; the
// regions action remove removes exactly the bins of the regions (3.7.6 also
// dropped the NaN bins); the "discounted regions" warning names the
// chromosome of its region; the correction factors follow their bins when
// the bins are reordered (hiCMatrix.reorderBins of 3.7.6 left them).
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
// place (hicx::coo_ops::assign_zero_block).
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
        const hicx::coo_ops::BlockAssignment<T> result = hicx::coo_ops::assign_zero_block(
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

// hiCMatrix.reorderBins of 3.7.6 left the correction factors in their old
// order. Fixed in HiCExplorer 4: after the reordering the factors, taken from
// before it, follow their bins. `order` lists the old bin of every new bin.
void carry_factors(py::handle hic_matrix, py::handle factors_before, const std::vector<std::int64_t>& order) {
    if (factors_before.is_none()) {
        return;
    }
    py::list indices;
    for (const std::int64_t bin : order) {
        indices.append(py::int_(bin));
    }
    const py::object index = py::module_::import("numpy").attr("array")(indices, py::arg("dtype") = "int64");
    hic_matrix.attr("correction_factors") = factors_before.attr("__getitem__")(index);
}

std::vector<std::int64_t> bins_of_chromosomes(py::handle hic_matrix, py::handle names) {
    std::vector<std::int64_t> order;
    for (const py::handle name : names) {
        const py::object range = hic_matrix.attr("getChrBinRange")(name);
        const std::int64_t first = to_int64(lc::item(range, 0));
        const std::int64_t last = to_int64(lc::item(range, 1));
        for (std::int64_t bin = first; bin < last; ++bin) {
            order.push_back(bin);
        }
    }
    return order;
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
        py::list unknown;
        log(kLogger, "debug", "pArgs.chromosomes " + str_of(chromosomes));
        for (const py::handle chromosome : chromosomes) {
            if (contains(chromosomes_list, chromosome)) {
                to_operate.append(chromosome);
            } else {
                unknown.append(chromosome);
            }
        }
        if (!unknown.empty()) {
            // Fixed in HiCExplorer 4: 3.7.6 only warned and went on.
            log(kLogger, "error",
                "Chromosome not available in matrix: " + str_of(matrix_path) + " " + str_of(unknown) +
                    ". Available: " + str_of(chromosomes_list));
            lc::python_exit(1);
        }
        log(kLogger, "debug", "85: chromosomes_list_to_operate_on " + str_of(to_operate));
        if (to_operate.empty()) {
            log(kLogger, "error",
                "No valid chromosome given: " + str_of(chromosomes) + ". Available: " + str_of(chromosomes_list));
            lc::python_exit(1);
        }
        const py::object factors_before = hic_matrix.attr("correction_factors");
        if (action_is("keep")) {
            const std::vector<std::int64_t> order = bins_of_chromosomes(hic_matrix, to_operate);
            hic_matrix.attr("reorderChromosomes")(to_operate);
            carry_factors(hic_matrix, factors_before, order);
        } else if (action_is("remove")) {
            log(kLogger, "debug", "chromosomes_list " + str_of(chromosomes_list));
            for (const py::handle chromosome : to_operate) {
                if (contains(chromosomes_list, chromosome)) {
                    chromosomes_list.attr("remove")(chromosome);
                    log(kLogger, "debug", "chromosome " + str_of(chromosome));
                }
            }
            log(kLogger, "debug", "chromosomes_list " + str_of(chromosomes_list));
            const std::vector<std::int64_t> order = bins_of_chromosomes(hic_matrix, chromosomes_list);
            hic_matrix.attr("reorderChromosomes")(chromosomes_list);
            carry_factors(hic_matrix, factors_before, order);
        } else if (action_is("mask")) {
            hic_matrix.attr("maskChromosomes")(to_operate);
        }
    } else if (lc::truthy(regions) || lc::truthy(pArgs.attr("maskBadRegions"))) {
        // --maskBadRegions masks the bins of its BED file whatever the action
        // says (fixed in HiCExplorer 4; 3.7.6 never read the file).
        const bool bad_regions = !lc::truthy(regions);
        const py::object bed_file = bad_regions ? pArgs.attr("maskBadRegions") : regions;
        const bool mask_action = bad_regions;
        auto region_action_is = [&](const char* name) {
            return mask_action ? std::string(name) == "mask" : action_is(name);
        };
        hic_matrix = load_matrix(matrix_path, py::none());
        const py::list chromosomes_list = py::list(hic_matrix.attr("chrBinBoundaries"));
        py::list genomic_regions;
        py::object chrom = py::none();
        for (const std::string& text : read_text_lines(bed_file)) {
            const py::str line = to_str(text);
            const py::list fields = line.attr("strip")().attr("split")("\t");
            log(kLogger, "debug", "_line " + str_of(fields));
            if (lc::length(line) < 3) {
                log(kLogger, "warning", "An entry shorter than 3 columns has been found!");
                continue;
            }
            if (fields.size() >= 3) {
                chrom = fields[0];
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
                if (region_action_is("remove")) {
                    // Fixed in HiCExplorer 4: the chromosome of the region,
                    // not the leaked loop variable of the last BED line.
                    const py::object region_chrom = lc::item(region, 0);
                    const py::object chr_range = hic_matrix.attr("getChrBinRange")(region_chrom);
                    const std::int64_t chr_start = to_int64(lc::item(chr_range, 0));
                    const std::int64_t chr_end = to_int64(lc::item(chr_range, 1));
                    if (start > chr_start && end < chr_end - 1) {
                        log(kLogger, "warning",
                            str_of(region_chrom) + ":" + std::to_string(start) + "-" + std::to_string(end) +
                                " entry may generate discounted regions on a chromosome.Please consider using "
                                "`mask` action to deal with that.");
                    }
                }
            }
        }
        if (region_action_is("keep")) {
            std::vector<std::int64_t> order;
            for (const py::handle bin : matrix_indices_regions) {
                order.push_back(to_int64(bin));
            }
            const py::object factors_before = hic_matrix.attr("correction_factors");
            hic_matrix.attr("reorderBins")(matrix_indices_regions);
            carry_factors(hic_matrix, factors_before, order);
        } else if (region_action_is("mask")) {
            hic_matrix.attr("maskBins")(matrix_indices_regions);
        } else if (region_action_is("remove")) {
            // Fixed in HiCExplorer 4: exactly the bins of the regions are
            // removed; 3.7.6 masked them and cleared orig_bin_ids,
            // orig_cut_intervals and nan_bins, which dropped the NaN bins too.
            const std::int64_t rows = to_int64(lc::item(hic_matrix.attr("matrix").attr("shape"), 0));
            std::vector<char> drop(static_cast<std::size_t>(rows), 0);
            for (const py::handle bin : matrix_indices_regions) {
                drop[static_cast<std::size_t>(to_int64(bin))] = 1;
            }
            std::vector<std::int64_t> survivors;
            py::list survivor_list;
            for (std::int64_t bin = 0; bin < rows; ++bin) {
                if (drop[static_cast<std::size_t>(bin)] == 0) {
                    survivors.push_back(bin);
                    survivor_list.append(py::int_(bin));
                }
            }
            const py::object factors_before = hic_matrix.attr("correction_factors");
            hic_matrix.attr("reorderBins")(survivor_list);
            carry_factors(hic_matrix, factors_before, survivors);
        }
    } else {
        log(kLogger, "error",
            "No data to adjust given. Please specify one of --chromosomes, --regions or --maskBadRegions.");
        lc::python_exit(1);
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
