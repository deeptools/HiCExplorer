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
// These live in core/include/hicx/validate_locations.hpp, which the public
// Python functions of hicexplorer.hicValidateLocations call too.
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

#include "hicx/argparse.hpp"
#include "hicx/bedtools_ops.hpp"
#include "hicx/bins.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/text_table.hpp"
#include "hicx/validate_locations.hpp"
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

// hicValidateLocations.py parse_arguments. --chrPrefixLoops and
// --chrPrefixProtein have None among their choices, which only the default can
// be: the string 'None' on the command line is an invalid choice, as in Python.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicValidateLocations",
                       "This script overlaps the loop locations with protein locations to determine "
                       "the accuracy of the loop detection.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    const std::vector<hicx::json::Value> prefix_choices{
        hicx::json::Value(), hicx::json::Value::string("add"), hicx::json::Value::string("remove")};
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--data", "-d"})
        .required()
        .input({"bedgraph", "bed"})
        .help("The loop file from hicDetectLoops, 'chr start end chr start end'; for TAD data "
              "the boundaries.bed file.");
    required.add({"--validationData", "-vd"})
        .required()
        .input({"bed", "cool"})
        .help("The data file to validate the given locations. Can be narrowPeak, broadPeak (both "
              "in bed), or cool");
    required.add({"--validationType", "-vt"})
        .choices({"bed", "cool"})
        .default_value("bed")
        .help("The type of the validation data. Can be bed, or cool format");
    required.add({"--method", "-m"})
        .choices({"loops", "tad"})
        .default_value("loops")
        .help("The method used.");
    required.add({"--resolution", "-r"})
        .type("int")
        .required()
        .help("The used resolution of the Hi-C interaction matrix.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--outFileName", "-o"})
        .output({"txt"}, "prefix")
        .help("The prefix name of the output files output_matched_locations and "
              "output_statistics.");
    optional.add({"--chrPrefixLoops", "-cl"})
        .choices(prefix_choices)
        .help("Adding / removing / do nothing a 'chr'-prefix to chromosome name of the loops.");
    optional.add({"--chrPrefixProtein", "-cp"})
        .choices(prefix_choices)
        .help("Adding / removing / do nothing a 'chr'-prefix to chromosome name of the protein.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.data = ns.str("data");
    args.validation_data = ns.str("validationData");
    args.validation_type = ns.str("validationType");
    args.method = ns.str("method");
    args.resolution = ns.integer("resolution");
    const std::optional<std::string> out_file_name = ns.opt_str("outFileName");
    args.has_out_file_name = out_file_name.has_value();
    args.out_file_name = out_file_name.value_or("");
    args.chr_prefix_loops = ns.opt_str("chrPrefixLoops").value_or("None");
    args.chr_prefix_protein = ns.opt_str("chrPrefixProtein").value_or("None");
    return args;
}

using hicx::validate_locations::apply_binning;
using hicx::validate_locations::apply_chr_prefix;
using hicx::validate_locations::correlate_cool;
using hicx::validate_locations::CoolPeaks;
using hicx::validate_locations::intervals_of;
using hicx::validate_locations::read_three_column;
using hicx::validate_locations::sorted_table;

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
