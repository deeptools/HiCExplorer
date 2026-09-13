// Port of hicexplorer/hicValidateLocations.py.
//
// Overlaps a loop or TAD boundary file with a set of protein peaks and reports
// how many calls are supported. No matrix arithmetic at all: the whole tool is
// pandas plus pybedtools, and the port is therefore a port of pandas and
// pybedtools, not of numerics. Three pieces of that carry real behaviour and
// live in the core rather than here:
//
//   * hicx::TextTable, which is pandas.read_csv(sep='\t', header=None) with
//     its dtype inference and its precise_xstrtod float parser, plus
//     DataFrame.to_csv. See core/include/hicx/text_table.hpp for why the float
//     parser is not strtod.
//   * hicx::bedtools::sort_order, which reproduces the row order bedtools sort
//     produces, including the permutation of records that share a start
//     coordinate. loops_1.bedgraph has 3,730 such records and the output file
//     is written in that order, so it is not a detail. See
//     core/include/hicx/bedtools_ops.hpp.
//   * hicx::bedtools::merge and intersect_count.
//
// Every BedTool.from_dataframe(...).sort().to_dataframe(...) in the Python
// writes a temporary file and reads it back, so a column's dtype is re-inferred
// from the text at each step. `roundtrip` below does the same thing in memory
// rather than assuming the round trip is a no operation: it is not, for a float
// column holding NaN, because from_dataframe writes na_rep='.' and read_csv
// does not recognise '.' as a missing value, so the column comes back as text.
//
// Two defects of the reference are reproduced rather than fixed, both pinned by
// the characterization tests added to
// hicexplorer/test/general/test_hicValidateLocations.py:
//
//   * `--method tad` writes the matched TADs to the bare --outFileName and not
//     to <name>_matched_locations as the help text for the option says
//     (hicValidateLocations.py:284 against :232). The statistics file of the
//     TAD branch also calls them Loops (":279-281").
//   * `--method tad` indexes the unbinned TAD frame with a mask computed on the
//     binned and de-duplicated one (":264"). When binning collapses two
//     boundaries into one row the two frames have different lengths and pandas
//     raises IndexingError: Unalignable boolean Series. The port reports that
//     instead of crashing, and exits 1 as the traceback does.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "hicx/bedtools_ops.hpp"
#include "hicx/bins.hpp"
#include "hicx/cool_file.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/text_table.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicValidateLocations --data DATA --validationData VALIDATIONDATA\n"
    "                            [--validationType {bed,cool}]\n"
    "                            [--method {loops,tad}] --resolution RESOLUTION\n"
    "                            [--outFileName OUTFILENAME]\n"
    "                            [--chrPrefixLoops {None,add,remove}]\n"
    "                            [--chrPrefixProtein {None,add,remove}] [--help]\n"
    "                            [--version]\n";

const char* const kHelp =
    "\n"
    "This script overlaps the loop locations with protein locations to determine "
    "the accuracy of the loop detection.\n"
    "Loops need to have format as follows:\n"
    "\n"
    "`chr start end chr start end`\n"
    "\n"
    "The protein peaks need to be in narrowPeaks or broadPeak format.\n"
    "\n"
    "A protein match is successfull if at the bin of the x and y location a "
    "protein peak is overlapped.\n"
    "A bin is assumed to have a protein if one or more protein peaks falling "
    "within the bin region.\n"
    "The value of the protein is not considered, only match or non-match.\n"
    "\n"
    "Required arguments:\n"
    "  --data DATA, -d DATA  The loop file from hicDetectLoops. To use files from\n"
    "                        other sources, please follow 'chr start end chr start\n"
    "                        end' format. For TAD data use the boundaries.bed file\n"
    "                        and not the domains file!\n"
    "  --validationData VALIDATIONDATA, -vd VALIDATIONDATA\n"
    "                        The data file to validate the given locations. Can be\n"
    "                        narrowPeak, broadPeak (both in bed), or cool\n"
    "  --validationType {bed,cool}, -vt {bed,cool}\n"
    "                        The type of the validation data. Can be bed, or cool\n"
    "                        format\n"
    "  --method {loops,tad}, -m {loops,tad}\n"
    "                        The method used (for the moment only loop is possible)\n"
    "                        (Default: loops).\n"
    "  --resolution RESOLUTION, -r RESOLUTION\n"
    "                        The used resolution of the Hi-C interaction matrix.\n"
    "\n"
    "Optional arguments:\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        The prefix name of the output files. Two file are\n"
    "                        written: output_matched_locations and\n"
    "                        output_statistics.First file contains all loop\n"
    "                        locations with protein location matches, second file\n"
    "                        contains statistics about this matching.\n"
    "  --chrPrefixLoops {None,add,remove}, -cl {None,add,remove}\n"
    "                        Adding / removing / do nothing a 'chr'-prefix to\n"
    "                        chromosome name of the loops.\n"
    "  --chrPrefixProtein {None,add,remove}, -cp {None,add,remove}\n"
    "                        Adding / removing / do nothing a 'chr'-prefix to\n"
    "                        chromosome name of the protein.\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::string data;
    std::string validation_data;
    std::string validation_type = "bed";
    std::string method = "loops";
    std::int64_t resolution = 0;
    std::string out_file_name;
    bool has_out_file_name = false;
    std::string chr_prefix_loops = "None";
    std::string chr_prefix_protein = "None";
};

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicValidateLocations: error: %s\n", message.c_str());
    std::exit(2);
}

// argparse's option name prefix matching, which the Python test suite relies on:
// test_hicValidateLocations.py:78 passes --chrPrefixLoop for --chrPrefixLoops.
bool matches_option(const std::string& token, const std::string& full,
                    const std::string& shrt) {
    if (token == shrt) {
        return true;
    }
    return token.rfind("--", 0) == 0 && token.size() > 2 &&
           full.rfind(token, 0) == 0;
}

Arguments parse_arguments(int argc, char** argv) {
    Arguments args;
    bool data_seen = false;
    bool validation_seen = false;
    bool resolution_seen = false;

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
        const auto value = [&]() -> std::string {
            if (i + 1 >= tokens.size()) {
                fail("argument " + token + ": expected one argument");
            }
            return tokens[++i];
        };
        if (token == "-h" || matches_option(token, "--help", "-h")) {
            std::fputs(kUsage, stdout);
            std::fputs(kHelp, stdout);
            std::exit(0);
        }
        if (matches_option(token, "--version", "")) {
            std::printf("hicValidateLocations %s\n", hicx::kVersion);
            std::exit(0);
        }
        if (matches_option(token, "--data", "-d")) {
            args.data = value();
            data_seen = true;
        } else if (matches_option(token, "--validationData", "-vd")) {
            args.validation_data = value();
            validation_seen = true;
        } else if (matches_option(token, "--validationType", "-vt")) {
            args.validation_type = value();
        } else if (matches_option(token, "--method", "-m")) {
            args.method = value();
        } else if (matches_option(token, "--resolution", "-r")) {
            const std::string text = value();
            try {
                args.resolution = std::stoll(text);
            } catch (const std::exception&) {
                fail("argument --resolution/-r: invalid int value: '" + text + "'");
            }
            resolution_seen = true;
        } else if (matches_option(token, "--outFileName", "-o")) {
            args.out_file_name = value();
            args.has_out_file_name = true;
        } else if (matches_option(token, "--chrPrefixLoops", "-cl")) {
            args.chr_prefix_loops = value();
        } else if (matches_option(token, "--chrPrefixProtein", "-cp")) {
            args.chr_prefix_protein = value();
        } else {
            fail("unrecognized arguments: " + token);
        }
    }

    std::string missing;
    const auto add_missing = [&missing](const char* name) {
        missing += missing.empty() ? name : std::string(", ") + name;
    };
    if (!data_seen) {
        add_missing("--data/-d");
    }
    if (!validation_seen) {
        add_missing("--validationData/-vd");
    }
    if (!resolution_seen) {
        add_missing("--resolution/-r");
    }
    if (!missing.empty()) {
        fail("the following arguments are required: " + missing);
    }
    if (args.validation_type != "bed" && args.validation_type != "cool") {
        fail("argument --validationType/-vt: invalid choice: '" + args.validation_type +
             "' (choose from 'bed', 'cool')");
    }
    if (args.method != "loops" && args.method != "tad") {
        fail("argument --method/-m: invalid choice: '" + args.method +
             "' (choose from 'loops', 'tad')");
    }
    for (const std::string* choice : {&args.chr_prefix_loops, &args.chr_prefix_protein}) {
        if (*choice != "None" && *choice != "add" && *choice != "remove") {
            fail("argument --chrPrefixLoops/-cl: invalid choice: '" + *choice +
                 "' (choose from None, 'add', 'remove')");
        }
    }
    return args;
}

// BedTool.from_dataframe(df) writes the frame with na_rep='.' and every later
// to_dataframe() reads the text back with read_csv, which re-infers the dtype
// of every column. TextTable::reparse does the same thing in memory rather
// than assuming the round trip is the identity: it is not, for a float column
// holding NaN, because from_dataframe writes na_rep='.' and read_csv does not
// recognise '.' as a missing value, so the column comes back as text.

// BedTool.from_dataframe(df).sort().to_dataframe(...)
hicx::TextTable sorted_table(const hicx::TextTable& table) {
    std::vector<std::string> chrom(table.rows());
    std::vector<std::int64_t> start(table.rows());
    for (std::size_t r = 0; r < table.rows(); ++r) {
        chrom[r] = table.column(0).text(r, ".");
        start[r] = table.column(1).as_int(r);
    }
    return table.select_rows(hicx::bedtools::sort_order(chrom, start)).reparse();
}

hicx::bedtools::Intervals intervals_of(const hicx::TextTable& table, std::size_t chrom_col,
                                       std::size_t start_col, std::size_t end_col) {
    hicx::bedtools::Intervals out;
    out.chrom.reserve(table.rows());
    out.start.reserve(table.rows());
    out.end.reserve(table.rows());
    for (std::size_t r = 0; r < table.rows(); ++r) {
        out.push_back(table.column(chrom_col).text(r, "."),
                      table.column(start_col).as_int(r),
                      table.column(end_col).as_int(r));
    }
    return out;
}

hicx::TextTable table_of(const hicx::bedtools::Intervals& intervals) {
    hicx::TableColumn chrom;
    chrom.type = hicx::ColumnType::String;
    chrom.strings = intervals.chrom;
    hicx::TableColumn start;
    start.type = hicx::ColumnType::Int64;
    start.ints = intervals.start;
    hicx::TableColumn end;
    end.type = hicx::ColumnType::Int64;
    end.ints = intervals.end;
    return hicx::TextTable::with_columns({chrom, start, end});
}

// applyBinning: floor the start onto the bin grid, take the bin after the end,
// drop exact duplicates, and optionally merge.
hicx::TextTable apply_binning(const hicx::TextTable& table, std::int64_t bin_size,
                              bool do_merge) {
    hicx::TextTable binned = table;
    hicx::TableColumn& start = binned.column(1);
    hicx::TableColumn& end = binned.column(2);
    std::vector<std::int64_t> new_start(binned.rows());
    std::vector<std::int64_t> new_end(binned.rows());
    const double divisor = static_cast<double>(bin_size);
    for (std::size_t r = 0; r < binned.rows(); ++r) {
        // (df[1] / binsize).astype(int) is a float64 division truncated
        // towards zero, not an integer division.
        const double scaled_start = static_cast<double>(table.column(1).as_int(r)) / divisor;
        const double scaled_end = static_cast<double>(table.column(2).as_int(r)) / divisor;
        new_start[r] = static_cast<std::int64_t>(scaled_start) * bin_size;
        new_end[r] = (static_cast<std::int64_t>(scaled_end) + 1) * bin_size;
    }
    start.type = hicx::ColumnType::Int64;
    start.ints = std::move(new_start);
    start.floats.clear();
    start.strings.clear();
    end.type = hicx::ColumnType::Int64;
    end.ints = std::move(new_end);
    end.floats.clear();
    end.strings.clear();

    binned.drop_duplicates();
    if (!do_merge) {
        return binned;
    }
    const hicx::bedtools::Intervals merged =
        hicx::bedtools::merge(intervals_of(binned, 0, 1, 2));
    hicx::TextTable merged_table = table_of(merged);
    return sorted_table(merged_table);
}

void apply_chr_prefix(hicx::TextTable& table, std::size_t column,
                      const std::string& mode) {
    if (mode == "add") {
        table.add_chr_prefix(column);
    } else if (mode == "remove") {
        table.remove_chr_prefix(column);
    }
}

// pd.read_csv(file, sep='\t', header=None)[[0, 1, 2]] with the prefix applied.
hicx::TextTable read_three_column(const std::string& path, const std::string& prefix) {
    hicx::TextTable table = hicx::TextTable::read_tsv(path).select_columns({0, 1, 2});
    apply_chr_prefix(table, 0, prefix);
    return table;
}

// correlateCool: cooler.Cooler(...).matrix(balance=False, sparse=True).fetch()
// for every loop, reduced to the only question the caller asks, whether the
// fetched block holds a stored pixel. Any exception inside the Python loop is
// swallowed and recorded as no match (hicValidateLocations.py:145-146), which
// covers a chromosome the cooler does not have and a coordinate past the end
// of a chromosome; both occur on the designated input.
struct CoolPeaks {
    std::int64_t number_of_peaks = 0;
    std::vector<bool> matched;
};

CoolPeaks correlate_cool(const std::string& path, const hicx::TextTable& loops) {
    hicx::CoolFile cool(path);
    CoolPeaks result;
    result.number_of_peaks = cool.nnz();

    const hicx::json::Value* bin_size_value = cool.info_value("bin-size");
    std::int64_t bin_size = 0;
    if (bin_size_value != nullptr && bin_size_value->is_number()) {
        bin_size = static_cast<std::int64_t>(bin_size_value->as_double());
    }
    if (bin_size <= 0) {
        throw std::runtime_error(
            "the validation cooler has no fixed bin size; cooler's variable bin "
            "path of region_to_extent is not ported");
    }

    std::map<std::string, std::int64_t> chrom_length;
    for (std::size_t i = 0; i < cool.chrom_names().size(); ++i) {
        chrom_length[cool.chrom_names()[i]] = cool.chrom_lengths()[i];
    }
    std::map<std::string, std::int64_t> chrom_offset;
    {
        const hicx::BinTable bins(cool.read_bins());
        for (const auto& entry : bins.chrom_bin_boundaries()) {
            chrom_offset[entry.first] = entry.second.first;
        }
    }

    const hicx::CsrMatrix matrix = cool.read_matrix();
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();

    // Any stored pixel with row in [r0, r1) and column in [c0, c1)?
    const auto block_has_pixel = [&](std::int64_t r0, std::int64_t r1, std::int64_t c0,
                                     std::int64_t c1) {
        if (r0 >= r1 || c0 >= c1) {
            return false;
        }
        const std::int64_t rows = matrix.rows();
        r0 = std::max<std::int64_t>(r0, 0);
        r1 = std::min<std::int64_t>(r1, rows);
        for (std::int64_t row = r0; row < r1; ++row) {
            const auto begin = indices.begin() + indptr[static_cast<std::size_t>(row)];
            const auto end = indices.begin() + indptr[static_cast<std::size_t>(row) + 1];
            const auto lower =
                std::lower_bound(begin, end, static_cast<std::int32_t>(c0));
            if (lower != end && *lower < static_cast<std::int32_t>(c1)) {
                return true;
            }
        }
        return false;
    };

    const auto extent = [&](const std::string& chrom, std::int64_t start,
                            std::int64_t end,
                            std::pair<std::int64_t, std::int64_t>* out) {
        const auto length = chrom_length.find(chrom);
        if (length == chrom_length.end()) {
            return false;  // KeyError inside cooler.parse_region
        }
        if (end < start || end > length->second || start < 0) {
            return false;  // ValueError inside cooler.parse_region
        }
        const std::int64_t offset = chrom_offset.at(chrom);
        const std::int64_t first = offset + start / bin_size;
        const std::int64_t last =
            offset + (end + bin_size - 1) / bin_size;  // ceil
        *out = {first, last};
        return true;
    };

    result.matched.assign(loops.rows(), false);
    for (std::size_t r = 0; r < loops.rows(); ++r) {
        std::pair<std::int64_t, std::int64_t> x{};
        std::pair<std::int64_t, std::int64_t> y{};
        if (!extent(loops.column(0).text(r, "."), loops.column(1).as_int(r),
                    loops.column(2).as_int(r), &x)) {
            continue;
        }
        if (!extent(loops.column(3).text(r, "."), loops.column(4).as_int(r),
                    loops.column(5).as_int(r), &y)) {
            continue;
        }
        // The cooler is stored symmetric-upper, so the block the Python gets
        // back is the union of the upper-triangle hits in either orientation.
        result.matched[r] =
            block_has_pixel(x.first, x.second, y.first, y.second) ||
            block_has_pixel(y.first, y.second, x.first, x.second);
    }
    return result;
}

void write_statistics(const std::string& path, const std::string& kind,
                      const std::string& data, const std::string& validation_data,
                      std::int64_t number_of_proteins, std::size_t matched,
                      std::size_t total) {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("cannot write '" + path + "'");
    }
    out << "# HiCExplorer hicValidateLocations " << hicx::kVersion << "\n";
    out << "# Overlap of " << kind << " file " << data << " with protein file "
        << validation_data << "\n#\n";
    out << "Protein peaks: " << number_of_proteins << "\n";
    out << "Matched Loops: " << matched << "\n";
    out << "Total Loops: " << total << "\n";
    out << "Loops match protein: "
        << hicx::npy::float_repr(static_cast<double>(matched) /
                                 static_cast<double>(total))
        << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    try {
        if (args.method == "loops") {
            hicx::TextTable loops = hicx::TextTable::read_tsv(args.data);
            if (loops.cols() < 6) {
                std::fprintf(stderr,
                             "hicValidateLocations: the loop file needs six columns "
                             "'chr start end chr start end', found %zu\n",
                             loops.cols());
                return 1;
            }
            apply_chr_prefix(loops, 0, args.chr_prefix_loops);
            apply_chr_prefix(loops, 3, args.chr_prefix_loops);
            loops = sorted_table(loops);

            std::int64_t number_of_proteins = 0;
            std::vector<bool> selection(loops.rows(), false);

            if (args.validation_type == "bed") {
                hicx::TextTable protein =
                    read_three_column(args.validation_data, args.chr_prefix_protein);
                protein = sorted_table(protein);
                const hicx::TextTable protein_resolution =
                    apply_binning(protein, args.resolution, true);
                number_of_proteins = static_cast<std::int64_t>(protein_resolution.rows());

                const hicx::bedtools::Intervals peaks =
                    intervals_of(protein_resolution, 0, 1, 2);
                const std::vector<std::int64_t> x =
                    hicx::bedtools::intersect_count(intervals_of(loops, 0, 1, 2), peaks);
                const std::vector<std::int64_t> y =
                    hicx::bedtools::intersect_count(intervals_of(loops, 3, 4, 5), peaks);
                for (std::size_t r = 0; r < loops.rows(); ++r) {
                    selection[r] = x[r] >= 1 && y[r] >= 1;
                }
            } else {
                const CoolPeaks peaks = correlate_cool(args.validation_data, loops);
                number_of_proteins = peaks.number_of_peaks;
                for (std::size_t r = 0; r < loops.rows(); ++r) {
                    selection[r] = peaks.matched[r];
                }
            }

            std::vector<std::size_t> matched_rows;
            for (std::size_t r = 0; r < loops.rows(); ++r) {
                if (selection[r]) {
                    matched_rows.push_back(r);
                }
            }

            std::printf("Protein peaks: %lld\n", static_cast<long long>(number_of_proteins));
            std::printf("Matched Loops: %zu\n", matched_rows.size());
            std::printf("Total Loops: %zu\n", loops.rows());
            std::printf("Loops match protein: %s\n",
                        hicx::npy::float_repr(static_cast<double>(matched_rows.size()) /
                                              static_cast<double>(loops.rows()))
                            .c_str());

            if (args.has_out_file_name) {
                loops.select_rows(matched_rows)
                    .write_tsv(args.out_file_name + "_matched_locations");
                write_statistics(args.out_file_name + "_statistics", "loop", args.data,
                                 args.validation_data, number_of_proteins,
                                 matched_rows.size(), loops.rows());
            }
        } else {
            hicx::TextTable tads = read_three_column(args.data, args.chr_prefix_loops);
            hicx::TextTable protein =
                read_three_column(args.validation_data, args.chr_prefix_protein);
            tads = sorted_table(tads);
            protein = sorted_table(protein);

            const hicx::TextTable tads_resolution =
                apply_binning(tads, args.resolution, false);
            const hicx::TextTable protein_resolution =
                apply_binning(protein, args.resolution, true);

            const std::vector<std::int64_t> counts = hicx::bedtools::intersect_count(
                intervals_of(tads_resolution, 0, 1, 2),
                intervals_of(protein_resolution, 0, 1, 2));

            if (counts.size() != tads.rows()) {
                // hicValidateLocations.py:264 indexes the unbinned frame with a
                // mask of a different length, which pandas rejects.
                std::fprintf(stderr,
                             "hicValidateLocations: binning collapsed %zu TAD rows into "
                             "%zu, and hicValidateLocations.py:264 indexes the unbinned "
                             "frame with the binned mask, which raises "
                             "pandas.errors.IndexingError: Unalignable boolean Series\n",
                             tads.rows(), counts.size());
                return 1;
            }

            std::vector<std::size_t> matched_rows;
            for (std::size_t r = 0; r < counts.size(); ++r) {
                if (counts[r] >= 1) {
                    matched_rows.push_back(r);
                }
            }

            std::printf("Protein peaks: %zu\n", protein_resolution.rows());
            std::printf("Matched TADs: %zu\n", matched_rows.size());
            std::printf("Total TADs: %zu\n", tads.rows());
            std::printf("TADs match protein: %s\n",
                        hicx::npy::float_repr(static_cast<double>(matched_rows.size()) /
                                              static_cast<double>(tads.rows()))
                            .c_str());

            if (args.has_out_file_name) {
                write_statistics(args.out_file_name + "_statistics", "TAD", args.data,
                                 args.validation_data,
                                 static_cast<std::int64_t>(protein_resolution.rows()),
                                 matched_rows.size(), tads.rows());
                // The bare name, not <name>_matched_locations: reproduced, see
                // the file comment.
                tads.select_rows(matched_rows).write_tsv(args.out_file_name);
            }
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicValidateLocations: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicValidateLocations");
    return 0;
}
