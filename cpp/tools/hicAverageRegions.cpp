// Port of hicexplorer/hicAverageRegions.py.
//
// Averages the contact submatrix around each region of a BED file and writes
// the result as a scipy sparse .npz.
//
// The arithmetic is a sum and a division, but the details of *which* windows
// are summed and *in what precision* are not obvious from the description, so
// they are spelled out here. All of them are pinned by the characterization
// tests added to hicexplorer/test/general/test_hicAverageRegions.py.
//
//  1. **A window is used only when it is exactly the output size.**
//     hicAverageRegions.py:190 compares the shape of the submatrix against the
//     shape of the accumulator and skips the region otherwise, with a warning.
//     Since the accumulator is D by D, that means the window must be D bins
//     wide and must lie inside the matrix. A region whose window was clipped
//     at a chromosome start or end is therefore dropped, not padded, and the
//     `start_out`/`end_out` bookkeeping at ":185-188" turns out to be dead:
//     when the shape test passes, `_start` is always 0 and `_end` is always at
//     least D, so both slices cover the whole accumulator.
//
//  2. **The count matrix is uniform.** It is incremented over the same full
//     slices, so every cell ends at the number of accepted regions and the
//     final division is by one scalar. It is written as a D by D array
//     regardless.
//
//  3. **The accumulator is float32 and the sum is rounded to float32 after
//     every region.** `summed_matrix` is a lil_matrix of dtype float32
//     (":171"), and `lil[a:b, a:b] += csr` reads the block, adds in the
//     promoted dtype and writes the result back into float32 storage. With an
//     int32 cool matrix the addition happens in float64 and is then rounded,
//     which is what `add_window` below does.
//
//  4. **The output dtype is float64, not float32.** The division at ":200"
//     used to produce a float32 result and the checked-in masters are float32,
//     but scipy 1.14 returns a float64 coo_matrix for `lil(float32) /=
//     ndarray(float64)`. The reference the harness compares against is the
//     current one, so the port writes float64. The masters in
//     test_data/hicAverageRegions are float32 and the Python test compares
//     them at decimal=0, which is why nobody noticed.
//
//  5. **Two reading quirks are reproduced.** `if len(line) == 0` at ":142"
//     tests the raw line, not the split fields, so it never fires for a
//     non-empty line; and a line with exactly one field leaves `viewpoint`
//     bound to the previous region, or raises UnboundLocalError on the first
//     line. The port reports the second case instead of crashing.
//
// Threading and SIMD: none. The window is 40 by 40 or 200 by 200 and there are
// six regions in the corpus; the C++ run costs 0.01 s of CPU against the
// Python's 3.4 s, all of it in the file layer and the interpreter start-up.
// cpp/OPTIMIZATION.md 6 wants a measurement before an optimisation and the
// measurement says there is nothing here.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/npz_file.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicAverageRegions --matrix MATRIX --regions REGIONS\n"
    "                         (--range RANGE RANGE | --rangeInBins RANGEINBINS "
    "RANGEINBINS)\n"
    "                         --outFileName OUTFILENAME [--help]\n"
    "                         [--coordinatesToBinMapping {start,center,end}]\n"
    "                         [--considerStrandDirection] [--version]\n";

const char* const kHelp =
    "\n"
    "       Sums Hi-C contacts around given reference points and computes their "
    "average. This tool is useful to detect differences at certain reference "
    "points as for example TAD boundaries between samples.\n"
    "\n"
    "WARNING: This tool can only be used with fixed bin size Hi-C matrices. No "
    "guarantees how and if it works on restriction site interaction matrices.\n"
    "\n"
    "options:\n"
    "  --range RANGE RANGE, -ra RANGE RANGE\n"
    "                        Range of region up- and downstream of each region to\n"
    "                        include in genomic units.\n"
    "  --rangeInBins RANGEINBINS RANGEINBINS, -rib RANGEINBINS RANGEINBINS\n"
    "                        Range of region up- and downstream of each region to\n"
    "                        include in bin units.\n"
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        The matrix to use for the average of TAD regions.\n"
    "  --regions REGIONS, -r REGIONS\n"
    "                        BED file which stores a list of regions that are\n"
    "                        summed and averaged\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        File name to save the average regions TADs matrix.\n"
    "\n"
    "Optional arguments:\n"
    "  --help, -h            show this help message and exit\n"
    "  --coordinatesToBinMapping {start,center,end}, -cb {start,center,end}\n"
    "                        If the region contains start and end coordinates,\n"
    "                        define if the start, center (start + (end-start) / 2)\n"
    "                        or end bin should be used as start for range.This\n"
    "                        parameter is only important to set if the given start\n"
    "                        and end coordinates are not in the same bin (Default:\n"
    "                        start).\n"
    "  --considerStrandDirection\n"
    "                        This parameter specifies if the strand information is\n"
    "                        taken into account for the aggregation. It has the\n"
    "                        effect that the contacts of a reverse strand region\n"
    "                        are inverted e.g. [1,2,3] becomes [3,2,1].\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::string matrix;
    std::string regions;
    std::string out_file_name;
    std::int64_t range_up = 0;
    std::int64_t range_down = 0;
    bool has_range = false;
    std::int64_t bins_up = 0;
    std::int64_t bins_down = 0;
    bool has_range_in_bins = false;
    std::string coordinates_to_bin_mapping = "start";
    bool consider_strand_direction = false;
};

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicAverageRegions: error: %s\n", message.c_str());
    std::exit(2);
}

std::int64_t parse_int(const std::string& text, const std::string& option) {
    try {
        std::size_t used = 0;
        const long long value = std::stoll(text, &used);
        if (used != text.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return value;
    } catch (const std::exception&) {
        fail("argument " + option + ": invalid int value: '" + text + "'");
    }
}

Arguments parse_arguments(int argc, char** argv) {
    Arguments args;
    bool matrix_seen = false;
    bool regions_seen = false;
    bool out_seen = false;

    std::vector<std::string> tokens;
    for (int i = 1; i < argc; ++i) {
        const std::string token(argv[i]);
        const std::size_t equals = token.find('=');
        if (equals != std::string::npos && token.rfind("--", 0) == 0) {
            tokens.push_back(token.substr(0, equals));
            tokens.push_back(token.substr(equals + 1));
        } else {
            tokens.push_back(token);
        }
    }

    for (std::size_t i = 0; i < tokens.size(); ++i) {
        const std::string& token = tokens[i];
        const auto next = [&]() -> std::string {
            if (i + 1 >= tokens.size()) {
                fail("argument " + token + ": expected one argument");
            }
            return tokens[++i];
        };
        if (token == "-h" || token == "--help") {
            std::fputs(kUsage, stdout);
            std::fputs(kHelp, stdout);
            std::exit(0);
        }
        if (token == "--version") {
            std::printf("hicAverageRegions %s\n", hicx::kVersion);
            std::exit(0);
        }
        if (token == "-m" || token == "--matrix") {
            args.matrix = next();
            matrix_seen = true;
        } else if (token == "-r" || token == "--regions") {
            args.regions = next();
            regions_seen = true;
        } else if (token == "-o" || token == "--outFileName") {
            args.out_file_name = next();
            out_seen = true;
        } else if (token == "-ra" || token == "--range") {
            args.range_up = parse_int(next(), "--range/-ra");
            args.range_down = parse_int(next(), "--range/-ra");
            args.has_range = true;
        } else if (token == "-rib" || token == "--rangeInBins") {
            args.bins_up = parse_int(next(), "--rangeInBins/-rib");
            args.bins_down = parse_int(next(), "--rangeInBins/-rib");
            args.has_range_in_bins = true;
        } else if (token == "-cb" || token == "--coordinatesToBinMapping") {
            args.coordinates_to_bin_mapping = next();
        } else if (token == "--considerStrandDirection") {
            args.consider_strand_direction = true;
        } else {
            fail("unrecognized arguments: " + token);
        }
    }

    std::string missing;
    const auto add_missing = [&missing](const char* name) {
        missing += missing.empty() ? name : std::string(", ") + name;
    };
    if (!matrix_seen) {
        add_missing("--matrix/-m");
    }
    if (!regions_seen) {
        add_missing("--regions/-r");
    }
    if (!out_seen) {
        add_missing("--outFileName/-o");
    }
    if (!missing.empty()) {
        fail("the following arguments are required: " + missing);
    }
    if (args.has_range && args.has_range_in_bins) {
        fail("argument --rangeInBins/-rib: not allowed with argument --range/-ra");
    }
    if (!args.has_range && !args.has_range_in_bins) {
        fail("one of the arguments --range/-ra --rangeInBins/-rib is required");
    }
    if (args.coordinates_to_bin_mapping != "start" &&
        args.coordinates_to_bin_mapping != "center" &&
        args.coordinates_to_bin_mapping != "end") {
        fail("argument --coordinatesToBinMapping/-cb: invalid choice: '" +
             args.coordinates_to_bin_mapping +
             "' (choose from 'start', 'center', 'end')");
    }
    return args;
}

struct Region {
    std::int64_t start_bin = 0;
    std::int64_t end_bin = 0;
    char orientation = 0;  // 0 for "no strand column was read"
};

std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (true) {
        const std::size_t marker = line.find('\t', start);
        if (marker == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, marker - start));
        start = marker + 1;
    }
    return fields;
}

std::string strip(const std::string& text) {
    std::size_t begin = 0;
    std::size_t end = text.size();
    while (begin < end && std::isspace(static_cast<unsigned char>(text[begin])) != 0) {
        ++begin;
    }
    while (end > begin && std::isspace(static_cast<unsigned char>(text[end - 1])) != 0) {
        --end;
    }
    return text.substr(begin, end - begin);
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    try {
        hicx::ToolMatrix hic = hicx::ToolMatrix::load(args.matrix);
        const hicx::BinTable bins(hic.data().cut_intervals);
        const std::int64_t bin_size = bins.bin_size();

        std::ifstream regions_file(args.regions);
        if (!regions_file) {
            std::fprintf(stderr, "hicAverageRegions: cannot open '%s'\n",
                         args.regions.c_str());
            return 1;
        }

        std::vector<Region> regions;
        std::string line;
        bool have_viewpoint = false;
        std::string chrom;
        std::string viewpoint_start;
        std::string viewpoint_end;
        int line_number = 0;

        while (std::getline(regions_file, line)) {
            ++line_number;
            const std::vector<std::string> fields = split_tabs(strip(line));
            // hicAverageRegions.py:142 tests len(line), the raw line, so a
            // blank line is not skipped here; readlines() on a file ending in a
            // newline does not produce one either.
            if (fields.size() == 1 && fields[0].empty()) {
                continue;
            }
            std::string strand;
            if (fields.size() == 2) {
                chrom = fields[0];
                viewpoint_start = fields[1];
                viewpoint_end = fields[1];
                have_viewpoint = true;
            } else if (fields.size() >= 3) {
                chrom = fields[0];
                viewpoint_start = fields[1];
                viewpoint_end = fields[2];
                if (args.consider_strand_direction && fields.size() < 6) {
                    std::fputs(
                        "ERROR:hicexplorer.hicAverageRegions:Strand orientation should "
                        "be considered but file does not contain the 6th column of the "
                        "bed file containing this information. Exiting!\n",
                        stderr);
                    return 1;
                }
                have_viewpoint = true;
                if (args.consider_strand_direction) {
                    strand = fields[5];
                }
            } else if (!have_viewpoint) {
                std::fprintf(stderr,
                             "hicAverageRegions: line %d has one field and no region "
                             "has been read yet; hicAverageRegions.py:155 then reads an "
                             "unbound `viewpoint` and raises UnboundLocalError\n",
                             line_number);
                return 1;
            }

            const std::optional<hicx::BinRange> chrom_range = bins.chrom_bin_range(chrom);
            if (!chrom_range.has_value()) {
                std::fprintf(stderr,
                             "hicAverageRegions: chromosome '%s' of line %d is not in "
                             "the matrix; the Python raises KeyError here\n",
                             chrom.c_str(), line_number);
                return 1;
            }

            const std::int64_t start_value = std::stoll(viewpoint_start);
            const std::int64_t end_value = std::stoll(viewpoint_end);
            // int(float(start) + (float(end) - float(start)) / 2), truncated
            // towards zero.
            const std::int64_t center_value = static_cast<std::int64_t>(
                static_cast<double>(start_value) +
                (static_cast<double>(end_value) - static_cast<double>(start_value)) / 2.0);

            std::int64_t start_bin = 0;
            std::int64_t end_bin = 0;

            if (args.has_range) {
                // calculateViewpointRange
                const std::int64_t max_length =
                    bins.bin_pos(static_cast<std::size_t>(chrom_range->last - 1)).end;
                std::int64_t region_start = 0;
                std::int64_t region_end = 0;
                if (args.coordinates_to_bin_mapping == "start") {
                    region_start = start_value - args.range_up;
                    region_end = start_value + args.range_down;
                } else if (args.coordinates_to_bin_mapping == "end") {
                    region_start = end_value - args.range_up;
                    region_end = end_value + args.range_down;
                } else {
                    region_start = center_value - args.range_up;
                    region_end = center_value + args.range_down;
                }
                if (region_start < 0) {
                    region_start = 0;
                }
                if (region_end > max_length) {
                    region_end = max_length - 1;
                }
                const std::optional<std::pair<std::int64_t, std::int64_t>> range =
                    bins.region_bin_range(chrom, region_start, region_end);
                if (!range.has_value()) {
                    std::fprintf(stderr,
                                 "hicAverageRegions: region %s:%lld-%lld of line %d is "
                                 "not inside the matrix; hicAverageRegions.py:157 "
                                 "unpacks None here and raises TypeError\n",
                                 chrom.c_str(), static_cast<long long>(region_start),
                                 static_cast<long long>(region_end), line_number);
                    return 1;
                }
                start_bin = range->first;
                end_bin = range->second;
            } else {
                // calculateViewpointRangeBins
                std::int64_t viewpoint_index = 0;
                if (args.coordinates_to_bin_mapping == "start") {
                    const auto range =
                        bins.region_bin_range(chrom, start_value, end_value);
                    if (!range.has_value()) {
                        std::fprintf(stderr,
                                     "hicAverageRegions: region of line %d is not inside "
                                     "the matrix; the Python raises TypeError here\n",
                                     line_number);
                        return 1;
                    }
                    viewpoint_index = range->first;
                } else if (args.coordinates_to_bin_mapping == "end") {
                    const auto range =
                        bins.region_bin_range(chrom, start_value, end_value);
                    if (!range.has_value()) {
                        std::fprintf(stderr,
                                     "hicAverageRegions: region of line %d is not inside "
                                     "the matrix; the Python raises TypeError here\n",
                                     line_number);
                        return 1;
                    }
                    viewpoint_index = range->second;
                } else {
                    const auto range =
                        bins.region_bin_range(chrom, center_value, center_value);
                    if (!range.has_value()) {
                        std::fprintf(stderr,
                                     "hicAverageRegions: centre of line %d is not inside "
                                     "the matrix; the Python raises TypeError here\n",
                                     line_number);
                        return 1;
                    }
                    viewpoint_index = range->second;
                }
                start_bin = viewpoint_index - args.bins_up;
                end_bin = viewpoint_index + args.bins_down;
                if (start_bin < chrom_range->first) {
                    start_bin = chrom_range->first;
                }
                if (end_bin > chrom_range->last) {
                    end_bin = chrom_range->last;
                }
            }

            Region region;
            region.start_bin = start_bin;
            region.end_bin = end_bin;
            region.orientation =
                args.consider_strand_direction && !strand.empty() ? strand[0] : 0;
            regions.push_back(region);
        }

        std::int64_t dimension = 0;
        if (args.has_range) {
            if (bin_size <= 0) {
                std::fprintf(stderr,
                             "hicAverageRegions: the matrix has no usable bin size\n");
                return 1;
            }
            // Python's floor division on non negative values.
            dimension = args.range_up / bin_size + args.range_down / bin_size;
        } else {
            dimension = args.bins_up + args.bins_down;
        }
        if (dimension <= 0) {
            std::fprintf(stderr,
                         "hicAverageRegions: the requested range is smaller than one "
                         "bin, so the output matrix would be empty\n");
            return 1;
        }

        const std::size_t side = static_cast<std::size_t>(dimension);
        // The accumulator, float32 as the lil_matrix is.
        std::vector<float> accumulator(side * side, 0.0F);
        std::vector<double> window(side * side, 0.0);
        std::int64_t accepted = 0;
        const hicx::CsrMatrix& matrix = hic.matrix();
        const std::int64_t nbins = matrix.rows();

        for (const Region& region : regions) {
            const std::int64_t length = region.end_bin - region.start_bin;
            // summed_matrix.shape != submatrix.shape: the window has to be
            // exactly the output size and has to lie inside the matrix.
            if (length != dimension || region.start_bin < 0 || region.end_bin > nbins) {
                std::fputs("WARNING:hicexplorer.hicAverageRegions:Shape of a submatrix "
                           "does not match. It is ignored.\n",
                           stderr);
                continue;
            }
            std::fill(window.begin(), window.end(), 0.0);
            // The symmetric window out of the upper triangle storage: an entry
            // (r, c) of the triangle lands at (r, c) and, when r != c, at
            // (c, r), and both are inside the window exactly when r and c are.
            for (std::int64_t row = region.start_bin; row < region.end_bin; ++row) {
                const std::size_t begin =
                    static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row)]);
                const std::size_t end = static_cast<std::size_t>(
                    matrix.indptr()[static_cast<std::size_t>(row) + 1]);
                for (std::size_t k = begin; k < end; ++k) {
                    const std::int64_t column = matrix.indices()[k];
                    if (column < region.start_bin || column >= region.end_bin) {
                        continue;
                    }
                    const std::size_t local_row =
                        static_cast<std::size_t>(row - region.start_bin);
                    const std::size_t local_column =
                        static_cast<std::size_t>(column - region.start_bin);
                    window[local_row * side + local_column] = matrix.data()[k];
                    if (matrix.symmetry() == hicx::Symmetry::UpperTriangle &&
                        local_row != local_column) {
                        window[local_column * side + local_row] = matrix.data()[k];
                    }
                }
            }
            ++accepted;
            for (std::size_t row = 0; row < side; ++row) {
                for (std::size_t column = 0; column < side; ++column) {
                    // '-' transposes the window (the reverse strand case).
                    const double value =
                        region.orientation == '-'
                            ? window[column * side + row]
                            : window[row * side + column];
                    float& target = accumulator[row * side + column];
                    target = static_cast<float>(static_cast<double>(target) + value);
                }
            }
        }

        // summed_matrix /= count_matrix, where every cell of the count matrix
        // holds the number of accepted regions.
        //
        // scipy does not divide. _spbase._divide with a dense operand computes
        // `recip = np.true_divide(1., other)` and then `self.multiply(recip)`
        // (scipy/sparse/_base.py), so the result is `value * (1/N)` and not
        // `value / N`, and the two differ in the last bit: with N = 6 and a
        // value of 5 the division gives 0.8333333333333334 and the
        // multiplication 0.8333333333333333. Three of the 252 values of the
        // designated case sit on that difference, so it is reproduced rather
        // than absorbed by the tolerance. Only stored entries are touched, and
        // the result is float64 whatever the accumulator's dtype was.
        const double reciprocal = 1.0 / static_cast<double>(accepted);
        std::vector<std::int32_t> indptr(side + 1, 0);
        std::vector<std::int32_t> indices;
        std::vector<double> data;
        for (std::size_t row = 0; row < side; ++row) {
            for (std::size_t column = 0; column < side; ++column) {
                const float value = accumulator[row * side + column];
                if (value == 0.0F) {
                    continue;  // never stored, so never divided
                }
                indices.push_back(static_cast<std::int32_t>(column));
                data.push_back(static_cast<double>(value) * reciprocal);
            }
            indptr[row + 1] = static_cast<std::int32_t>(indices.size());
        }

        hicx::npz::save_csr_npz(args.out_file_name, dimension, dimension, indptr,
                                indices, data);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicAverageRegions: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicAverageRegions");
    return 0;
}
