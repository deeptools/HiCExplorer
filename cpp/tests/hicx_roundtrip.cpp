// Reads a contact matrix and writes it out again, in either format.
//
// This is not a HiCExplorer tool. It is the smallest program that exercises
// the tier 0 file layer end to end, and it is what the equivalence harness
// runs against its Python counterpart, cpp/scripts/py_roundtrip.py. Both take
// the same path through their respective libraries, the one hicConvertFormat
// takes: load through the file handler for the input format, hand the five
// loader values to the handler for the output format, save. Nothing above the
// file layer, in particular nothing of hiCMatrix, is involved.
//
//   hicx_roundtrip <input> <output> [--enforce-integer] [--no-symmetric]
//                                   [--no-apply-correction] [--verbose]
//
// The format is chosen by the file name suffix, exactly as
// hicmatrix.HiCMatrix does it: '.h5' means the HiCExplorer h5 format and
// everything else means cool.

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "hicx/cool_file.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/matrix_data.hpp"
#include "hicx/resource_usage.hpp"

namespace {

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

int usage() {
    std::cerr << "usage: hicx_roundtrip <input> <output> [--enforce-integer] "
                 "[--no-symmetric] [--no-apply-correction] [--verbose]\n";
    return 2;
}

}  // namespace

int main(int argc, char** argv) {
    std::vector<std::string> positional;
    bool enforce_integer = false;
    bool symmetric = true;
    bool apply_correction = true;
    bool verbose = false;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "--enforce-integer") {
            enforce_integer = true;
        } else if (argument == "--no-symmetric") {
            symmetric = false;
        } else if (argument == "--no-apply-correction") {
            apply_correction = false;
        } else if (argument == "--verbose") {
            verbose = true;
        } else if (!argument.empty() && argument[0] == '-') {
            return usage();
        } else {
            positional.push_back(argument);
        }
    }
    if (positional.size() != 2) {
        return usage();
    }
    const std::string& input = positional[0];
    const std::string& output = positional[1];

    try {
        hicx::MatrixData data;
        hicx::CoolSaveOptions save_options;
        save_options.symmetric = symmetric;
        save_options.apply_correction = apply_correction;
        save_options.enforce_integer = enforce_integer;

        if (ends_with(input, ".h5")) {
            data = hicx::read_hicexplorer_h5(input);
            save_options.file_was_h5 = true;
        } else {
            hicx::CoolLoadResult loaded = hicx::read_cool(input);
            data = std::move(loaded.data);
            save_options.correction_operator = loaded.correction_operator;
            save_options.hic2cool_version = loaded.hic2cool_version;
            save_options.hic_metadata = std::move(loaded.metadata);
            save_options.has_hic_metadata = true;
        }

        if (ends_with(output, ".h5")) {
            hicx::H5SaveOptions h5_options;
            h5_options.symmetric = symmetric;
            hicx::write_hicexplorer_h5(output, data, h5_options);
        } else {
            hicx::write_cool(output, data, save_options);
        }
    } catch (const std::exception& error) {
        std::cerr << "hicx_roundtrip: " << error.what() << "\n";
        return 1;
    }

    if (verbose) {
        std::cerr << "peak resident set size: " << hicx::peak_rss_kb() << " kB\n";
    }
    return 0;
}
