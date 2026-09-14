// Port of hicexplorer/hicCompartmentalization.py, numeric outputs only.
//
// The tool sorts the bins into quantiles of a first principal component,
// averages the obs/exp matrix over every pair of quantiles, and reports how
// much stronger the interaction is within the A and B ends of that ordering
// than between them (Schwarzer et al. 2017). It writes three things:
//
//   <outputFileName>        a matplotlib line plot of the ratios
//   <outputFileName>_dat    np.savetxt of the ratios, one row per matrix
//   --outputMatrix          np.savez of the normalised sum per quantile
//
// ---------------------------------------------------------------------------
// The plot
// ---------------------------------------------------------------------------
// The figure belongs to tier 7 of cpp/PLAN.md, whose strategy the project owner
// has not decided yet, so the port draws nothing. Because --outputFileName is
// required and the Python always draws it, a plain invocation is a request for
// the figure, and the port refuses it with exit status 1 before reading any
// input, rather than writing the two numeric files and quietly leaving the
// figure out. --noPlot, which the Python does not have, is the explicit way to
// ask for the numeric outputs alone; --outputFileName then only names the
// _dat file. Both files are produced by exactly the computation the figure
// would plot.
//
// ---------------------------------------------------------------------------
// Reproduced quirks, pinned by hicexplorer/test/general/test_hicCompartmentalization.py
// ---------------------------------------------------------------------------
//  1. The first ratio divides by an empty block. Quantile 0 is always empty
//     (quirk 5), so for q = 1 the between-compartment block [0, Q-1] is zero and
//     np.savetxt writes `inf`, or `nan` when the within block is zero too.
//  2. Every quantile pair is counted twice on the diagonal: :125-128 add each
//     block to [qi, qj] and to [qj, qi], so a block with qi == qj lands in the
//     same cell twice, bin count included. The normalised value is unchanged.
//  3. The loop `for chrom in chromosomes` (:106) never uses `chrom`: the per
//     chromosome slicing is commented out, so the whole genome-wide count is
//     repeated once per chromosome of the pca file and both accumulators scale
//     with the chromosome count, which leaves the normalised sums unchanged to
//     rounding but not bit for bit: s + s + s is not always 3 * s, so the last
//     bits depend on the count. Its comment says "It is only handeling cis
//     contacts", but the blocks mix the bins of all chromosomes, so trans
//     contacts are counted as well. The port computes each block once and
//     replays the additions once per chromosome, which is the same floating
//     point sequence at a fraction of the cost.
//  4. :199 `pc1.loc[pc1["pc1"] == np.nan]["quantile"] = args.quantile + 1` does
//     nothing twice over: `== np.nan` is never true, and the chained indexing
//     assigns into a copy. A bin whose pc1 is NaN gets np.searchsorted's
//     answer for NaN, Q, and is silently left out of every block.
//  5. The quantile boundaries run from the minimum to the maximum pc1 value, and
//     np.searchsorted(side='right') puts a value equal to a boundary *after*
//     it. So no bin ever lands in quantile 0 (the first row and column of
//     --outputMatrix are always zero), and the bin or bins holding the maximum
//     pc1 value land in quantile Q, outside the range, and are dropped like the
//     NaN bins of quirk 4. On hicCompartmentalization/pca1.bedgraph that is one
//     50 kb bin of 449. With --outliers the boundaries are a linspace between
//     two trimmed quantiles and the trimmed tail above the upper one is
//     dropped the same way.
//  6. np.savez appends '.npz' to an --outputMatrix name that lacks it, and is
//     handed a list, so the file holds one array `arr_0` of shape
//     (matrices, Q, Q).
//  7. --quantile 1 without --outliers raises ZeroDivisionError at :192;
//     --quantile 0 succeeds and writes one empty line per matrix.
//  8. --offset writes NaN into the named diagonals of the whole genome matrix,
//     so it also blanks the inter-chromosomal cells that sit on those diagonals
//     next to a chromosome border.
//
// ---------------------------------------------------------------------------
// Arithmetic
// ---------------------------------------------------------------------------
// Every sum is numpy's: the dense block of :121-125 is streamed through the
// 8192 element buffer and the pairwise kernel of hicx::npy, the slices of
// within_vs_between are summed as the contiguous copies numpy reduces, the
// quantiles are np.nanquantile's linear interpolation and the boundaries under
// --outliers np.linspace's. Nothing is threaded: see compartmentalization_impl.

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "compartmentalization_impl.hpp"
#include "hicx/argparse.hpp"
#include "hicx/bins.hpp"
#include "hicx/npz_file.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

namespace cm = hicx::compartments;

const char* const kUsage =
    "usage: hicCompartmentalization --obsexp_matrices OBSEXP_MATRICES\n"
    "                               [OBSEXP_MATRICES ...] --pca PCA\n"
    "                               --outputFileName OUTPUTFILENAME\n"
    "                               [--quantile QUANTILE] [--outliers OUTLIERS]\n"
    "                               [--outputMatrix OUTPUTMATRIX]\n"
    "                               [--offset OFFSET [OFFSET ...]] [--noPlot] [-h]\n"
    "                               [--version]\n";

const char* const kHelp =
    "\n"
    "Rearrange the average interaction frequencies using the first PC values to\n"
    "represent the global compartmentalization signal. To our knowledge this has been\n"
    "first introduced and implemented by Wibke Schwarzer et al. 2017 (Nature. 2017 Nov\n"
    "2; 551(7678): 51-56)\n"
    "\n"
    "$ hicCompartmentalization --obsexp_matrices obsExpMatrix.h5 --pca pc1.bedgraph -o\n"
    "global_signal.png --noPlot\n"
    "\n"
    "Required arguments:\n"
    "  --obsexp_matrices OBSEXP_MATRICES [OBSEXP_MATRICES ...], -m OBSEXP_MATRICES ...\n"
    "                        HiCExplorer matrices in h5/cool format.\n"
    "  --pca PCA             a PCA vector as a bedgraph file with no header.\n"
    "  --outputFileName OUTPUTFILENAME, -o OUTPUTFILENAME\n"
    "                        Plot to represent the polarization of A/B compartments.\n"
    "                        The ratios are written to OUTPUTFILENAME_dat.\n"
    "\n"
    "Optional arguments:\n"
    "  --quantile QUANTILE, -q QUANTILE\n"
    "                        number of quantiles. (Default: 30).\n"
    "  --outliers OUTLIERS   precentage of outlier to remove. (Default: 0).\n"
    "  --outputMatrix OUTPUTMATRIX\n"
    "                        output .npz file includes all the generated matrices\n"
    "  --offset OFFSET [OFFSET ...]\n"
    "                        set nan for the distances mentioned as offset from main\n"
    "                        diagonal, only positive values are accepted!\n"
    "  --noPlot              C++ port only: write the numeric outputs\n"
    "                        (OUTPUTFILENAME_dat and --outputMatrix) without the plot.\n"
    "                        Plotting is not yet available in the C++ port, so a run\n"
    "                        without this flag is refused.\n"
    "  -h                    show the help message and exit.\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::vector<std::string> matrices;
    std::string pca;
    std::string output_file_name;
    std::int64_t quantile = 30;
    double outliers = 0.0;
    std::optional<std::string> output_matrix;
    std::vector<std::int64_t> offset;
    bool no_plot = false;
};

// hicCompartmentalization.py parse_arguments, plus the C++-only --noPlot.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicCompartmentalization",
                       "Rearrange the average interaction frequencies using the first PC values to "
                       "represent the global compartmentalization signal.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--obsexp_matrices", "-m"})
        .nargs("+")
        .required()
        .input({"h5", "cool", "mcool"})
        .help("HiCExplorer matrices in h5/cool format.");
    required.add({"--pca"})
        .required()
        .input({"bedgraph"})
        .help("a PCA vector as a bedgraph file with no header.");
    required.add({"--outputFileName", "-o"})
        .required()
        .output({"png", "pdf", "svg"})
        .help("Plot to represent the polarization of A/B compartments.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--quantile", "-q"}).type("int").default_value(30).help("number of quantiles.");
    optional.add({"--outliers"})
        .type("float")
        .default_value(0)
        .help("precentage of outlier to remove.");
    optional.add({"--outputMatrix"})
        .output({"npz"})
        .help("output .npz file includes all the generated matrices");
    optional.add({"--offset"})
        .nargs("+")
        .type("int")
        .help("set nan for the distances mentioned as offset from main diagonal.");
    optional.add({"--noPlot"})
        .action(cli::Action::StoreTrue)
        .cpp_only("Plotting is not yet available in the C++ port; the flag writes the numeric "
                  "outputs without the required figure.")
        .help("write the numeric outputs without the plot.");
    optional.add({"-h"}).action(cli::Action::Help).help("show the help message and exit.");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.matrices = ns.strs("obsexp_matrices");
    args.pca = ns.str("pca");
    args.output_file_name = ns.str("outputFileName");
    args.quantile = ns.integer("quantile");
    args.outliers = ns.real("outliers");
    args.output_matrix = ns.opt_str("outputMatrix");
    args.offset = ns.integers("offset");
    args.no_plot = ns.flag("noPlot");
    return args;
}

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

void write_text_file(const std::string& path, const std::string& content) {
    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    if (!file) {
        throw std::runtime_error("[Errno " + std::to_string(errno) + "] " +
                                 std::strerror(errno) + ": '" + path + "'");
    }
    file << content;
    file.close();
    if (!file) {
        throw std::runtime_error("could not write '" + path + "'");
    }
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);

    if (!args.no_plot) {
        std::fprintf(
            stderr,
            "hicCompartmentalization: the polarization plot '%s' was requested, but "
            "plotting is not yet available in the C++ port (tier 7 of cpp/PLAN.md). "
            "Nothing was written. Pass --noPlot to write only the numeric outputs, "
            "'%s_dat' and --outputMatrix, or use the Python hicCompartmentalization "
            "for the figure.\n",
            args.output_file_name.c_str(), args.output_file_name.c_str());
        return 1;
    }

    try {
        const std::vector<cm::PcaRow> rows = cm::read_pca_bedgraph(args.pca);
        const std::vector<double> boundaries =
            cm::quantile_boundaries(rows, args.quantile, args.outliers);
        std::vector<double> keys;
        keys.reserve(rows.size());
        for (const auto& row : rows) {
            keys.push_back(static_cast<double>(row.pc1));
        }
        const std::vector<std::int64_t> quantile_of_row =
            cm::searchsorted_right(boundaries, keys);
        std::set<std::string> chromosomes;
        for (const auto& row : rows) {
            chromosomes.insert(row.chrom);
        }
        const auto chromosome_count = static_cast<std::int64_t>(chromosomes.size());

        std::vector<std::vector<double>> output_matrices;
        std::vector<std::vector<double>> polarization_ratio;
        for (const auto& path : args.matrices) {
            hicx::ToolMatrix obs_exp = hicx::ToolMatrix::load(path);
            const hicx::BinTable bins(obs_exp.cut_intervals());

            // get_indices (:87-90), evaluated for every row by pc1.apply.
            std::vector<std::vector<std::int64_t>> bin_ids(rows.size());
            for (std::size_t k = 0; k < rows.size(); ++k) {
                const cm::PcaRow& row = rows[k];
                if (!bins.chrom_bin_range(row.chrom).has_value()) {
                    throw std::runtime_error("ValueError: chromosome: " + row.chrom +
                                             " name not found in matrix " + path);
                }
                const auto range = bins.region_bin_range(row.chrom, row.start, row.end - 1);
                if (!range.has_value()) {
                    throw std::runtime_error(
                        "TypeError: getRegionBinRange found no bin for " + row.chrom + ":" +
                        std::to_string(row.start) + "-" + std::to_string(row.end) + " in " +
                        path + " and returned None");
                }
                for (std::int64_t bin = range->first; bin <= range->second; ++bin) {
                    bin_ids[k].push_back(bin);
                }
            }

            if (args.quantile < 0) {
                // np.zeros((quantiles_number, quantiles_number)) at :98
                throw std::runtime_error("ValueError: negative dimensions are not allowed");
            }
            for (const std::int64_t dist : args.offset) {
                if (dist < 0) {
                    // assert (dist >= 0) at :102
                    throw std::runtime_error("AssertionError: --offset " +
                                             std::to_string(dist) + " is negative");
                }
            }

            std::vector<double> normalised = cm::normalised_sum_per_quantile(
                obs_exp.matrix(), bin_ids, quantile_of_row, args.quantile, args.offset,
                chromosome_count);
            polarization_ratio.push_back(cm::within_vs_between(normalised, args.quantile));
            if (args.output_matrix.has_value()) {
                output_matrices.push_back(std::move(normalised));
            }
        }

        if (args.output_matrix.has_value()) {
            std::string npz_path = *args.output_matrix;
            if (!ends_with(npz_path, ".npz")) {
                npz_path += ".npz";  // quirk 6
            }
            hicx::npz::Array array;
            array.name = "arr_0";
            array.dtype = "<f8";
            array.shape = {static_cast<std::int64_t>(output_matrices.size()), args.quantile,
                           args.quantile};
            for (const auto& matrix : output_matrices) {
                array.data.append(reinterpret_cast<const char*>(matrix.data()),
                                  matrix.size() * sizeof(double));
            }
            // np.savez stores its entries uncompressed.
            hicx::npz::write_npz(npz_path, {array}, false);
        }

        std::string dat;
        for (const auto& ratios : polarization_ratio) {
            dat += cm::savetxt_line(ratios);
        }
        write_text_file(args.output_file_name + "_dat", dat);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicCompartmentalization: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicCompartmentalization");
    return 0;
}
