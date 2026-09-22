#include "hicx/tool_matrix.hpp"

#include <utility>

#include "hicx/h5_file.hpp"
#include "hicx/hic_adapter.hpp"
#include "hicx/matrix_ops.hpp"

namespace hicx {

namespace {

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

}  // namespace

ToolMatrix ToolMatrix::load(const std::string& path,
                            const std::optional<std::string>& chromosome) {
    ToolMatrix result;

    // hicmatrix.HiCMatrix.__init__:47-51 picks the format from the file name
    // suffix alone: '.h5' is the HiCExplorer format, everything else is cool.
    // v4 adds a third branch Python never had: a .hic source, detected by
    // content signature first and by extension as a fallback (hic_adapter.hpp
    // is_hic_path), read through hic_adapter's read_hic with the resolution
    // and normalisation selector documented there. `chromosome`, when given,
    // reaches read_hic exactly as it reaches CoolLoadOptions.chrom_name below:
    // a bare chromosome name or a "chrom:start-end" region loads only that
    // block, through hicfilecpp's block index, not the whole file (see
    // hic_adapter.hpp read_hic and hicPlotMatrix.cpp's fast path).
    if (is_hic_path(path)) {
        HicLoadResult loaded = read_hic(path, chromosome);
        result.data_ = std::move(loaded.data);
        result.cool_options_.correction_operator = loaded.correction_operator;
        result.cool_options_.hic_metadata = std::move(loaded.metadata);
        result.cool_options_.has_hic_metadata = true;
    } else if (ends_with(path, ".h5")) {
        result.data_ = read_hicexplorer_h5(path);
        result.input_is_h5_ = true;
    } else {
        CoolLoadOptions load_options;
        load_options.chrom_name = chromosome;
        CoolLoadResult loaded = read_cool(path, load_options);
        result.data_ = std::move(loaded.data);
        result.cool_options_.correction_operator = loaded.correction_operator;
        result.cool_options_.hic2cool_version = loaded.hic2cool_version;
        result.cool_options_.hic_metadata = std::move(loaded.metadata);
        result.cool_options_.has_hic_metadata = true;
    }

    // The load time field swap of cpp/PLAN.md 2.7 quirk 1. From here on the
    // two members mean what hiCMatrix means by them.
    std::swap(result.data_.correction_factors, result.data_.distance_counts);

    // fillLowerTriangle. The addition also drops entries stored as an exact
    // zero, so the C++ representation is normalised the same way.
    result.data_.matrix.symmetrize_in_place();
    if (result.data_.matrix.symmetry() == Symmetry::UpperTriangle) {
        result.data_.matrix.eliminate_zeros();
    }

    result.boundaries_ = chrom_bin_boundaries(result.data_.cut_intervals);
    return result;
}

void ToolMatrix::refresh_boundaries() {
    boundaries_ = chrom_bin_boundaries(data_.cut_intervals);
}

bool ToolMatrix::save(const std::string& path) {
    // hiCMatrix.save:99 guards the write with this test and does nothing when
    // it fails, so a name without a recognised suffix produces no output.
    if (!ends_with(path, "cool") && !ends_with(path, "h5")) {
        return false;
    }
    if (input_is_h5_) {
        // The handler built during the load is reused, so an h5 input is
        // written as h5 whatever the output name is; the writer appends '.h5'
        // when the name does not already end in it.
        H5SaveOptions options;
        options.symmetric = true;
        write_hicexplorer_h5(path, data_, options);
        return true;
    }
    CoolSaveOptions options = cool_options_;
    options.symmetric = true;
    // hiCMatrix.save's pApplyCorrection defaults to False, unlike the file
    // handler's own default.
    options.apply_correction = false;
    if (ends_with(path, "cool")) {
        // hiCMatrix.save:95 overwrites hic_metadata with pHiCInfo, which
        // defaults to None, so the source cooler's genome-assembly and
        // matrix-generated-by attributes are dropped. Only for a 'cool' name.
        options.has_hic_metadata = false;
        options.hic_metadata.clear();
    }
    write_cool(path, data_, options);
    return true;
}

}  // namespace hicx
