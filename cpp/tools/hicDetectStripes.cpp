// hicDetectStripes: architectural stripe detection, PLAN.md tier 9 section
// 9.3. HiCExplorer 3.7 has no stripe caller; there is no Python reference and
// the tool is validated as class EX (planted-stripe recovery and agreement
// with Stripenn, cpp/scripts/stripe_calibration.py).
//
// The statistical method lives in stripes_impl.{hpp,cpp}; this file is the
// argument parser, the per-chromosome loader (mirroring hicDetectLoops' band
// loader) and the writer.
//
// Output format (new, no Python precedent to match): tab separated, one call
// per line,
//
//   chrom  anchor_start  anchor_end  orientation  extent_start  extent_end
//   enrichment  pvalue  qvalue
//
// `orientation` is "vertical" or "horizontal". `anchor_start`/`anchor_end`
// are the anchor bin's genomic interval; `extent_start`/`extent_end` the
// genomic interval the stripe covers, anchor included, always with
// extent_start < extent_end regardless of orientation.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

#include "hicx/adjust_ops.hpp"
#include "hicx/argparse.hpp"
#include "hicx/bins.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/matrix_ops.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/stats_ops.hpp"
#include "hicx/version.hpp"
#include "stripes_impl.hpp"

namespace {

const char* const kUsage =
    "usage: hicDetectStripes --matrix MATRIX --outFileName OUTFILENAME\n"
    "                        [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]\n"
    "                        [--minStripeLength MINSTRIPELENGTH]\n"
    "                        [--maxStripeLength MAXSTRIPELENGTH]\n"
    "                        [--stripeLengthStep STRIPELENGTHSTEP]\n"
    "                        [--backgroundWindow BACKGROUNDWINDOW]\n"
    "                        [--backgroundGap BACKGROUNDGAP]\n"
    "                        [--obsExpThreshold OBSEXPTHRESHOLD]\n"
    "                        [--zScoreThreshold ZSCORETHRESHOLD]\n"
    "                        [--minRawCount MINRAWCOUNT]\n"
    "                        [--mergeWindow MERGEWINDOW] [--fdr FDR]\n"
    "                        [--threads THREADS] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Detects architectural stripes (vertical and horizontal) on the given\n"
    "contact matrix. New in HiCExplorer v4; PLAN.md tier 9 section 9.3.\n"
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        The matrix to compute the stripe detection on.\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        Outfile name to store the detected stripes\n"
    "                        (tab separated).\n"
    "\n"
    "Optional arguments:\n"
    "  --chromosomes CHROMOSOMES [CHROMOSOMES ...]\n"
    "                        Chromosomes to include in the analysis. If not set,\n"
    "                        all chromosomes are included.\n"
    "  --minStripeLength MINSTRIPELENGTH\n"
    "                        Shortest candidate stripe length, in base pairs.\n"
    "                        (Default: 100000).\n"
    "  --maxStripeLength MAXSTRIPELENGTH\n"
    "                        Longest candidate stripe length, in base pairs.\n"
    "                        (Default: 3000000).\n"
    "  --stripeLengthStep STRIPELENGTHSTEP\n"
    "                        Spacing of the candidate length grid, in base\n"
    "                        pairs. (Default: 100000).\n"
    "  --backgroundWindow BACKGROUNDWINDOW\n"
    "                        Flanking anchors on each side used as the local\n"
    "                        background, in bins. (Default: 15).\n"
    "  --backgroundGap BACKGROUNDGAP\n"
    "                        Anchors next to the stripe body excluded from the\n"
    "                        background, in bins. (Default: 2).\n"
    "  --obsExpThreshold OBSEXPTHRESHOLD\n"
    "                        Minimum enrichment of the stripe body's mean\n"
    "                        obs/exp over the local background's mean.\n"
    "                        (Default: 1.5).\n"
    "  --zScoreThreshold ZSCORETHRESHOLD\n"
    "                        Minimum preselection z-score against the local\n"
    "                        background. (Default: 2.0).\n"
    "  --minRawCount MINRAWCOUNT\n"
    "                        Minimum mean raw count over the stripe body.\n"
    "                        (Default: 1.0).\n"
    "  --mergeWindow MERGEWINDOW\n"
    "                        Non-maximum suppression window, in bins.\n"
    "                        (Default: 5).\n"
    "  --fdr FDR             Benjamini-Hochberg q-value threshold for the final\n"
    "                        call set. (Default: 0.05).\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads to use. (Default: 4).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::string matrix;
    std::string out_file_name;
    std::optional<std::vector<std::string>> chromosomes;
    std::int64_t min_stripe_length = 100000;
    std::int64_t max_stripe_length = 3000000;
    std::int64_t stripe_length_step = 100000;
    std::int64_t background_window = 15;
    std::int64_t background_gap = 2;
    double obs_exp_threshold = 1.5;
    double z_score_threshold = 2.0;
    double min_raw_count = 1.0;
    std::int64_t merge_window = 5;
    double fdr = 0.05;
    int threads = 4;
};

Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicDetectStripes",
                       "Detects architectural stripes (vertical and horizontal) on the given "
                       "contact matrix.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrix", "-m"})
        .required()
        .input({"h5", "cool", "mcool"})
        .help("The matrix to compute the stripe detection on.");
    required.add({"--outFileName", "-o"})
        .required()
        .output({"tsv"})
        .help("Outfile name to store the detected stripes.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--chromosomes"})
        .nargs("+")
        .help("Chromosomes to include in the analysis. If not set, all chromosomes are included.");
    optional.add({"--minStripeLength"})
        .type("int")
        .default_value(100000)
        .help("Shortest candidate stripe length, in base pairs.");
    optional.add({"--maxStripeLength"})
        .type("int")
        .default_value(3000000)
        .help("Longest candidate stripe length, in base pairs.");
    optional.add({"--stripeLengthStep"})
        .type("int")
        .default_value(100000)
        .help("Spacing of the candidate length grid, in base pairs.");
    optional.add({"--backgroundWindow"})
        .type("int")
        .default_value(15)
        .help("Flanking anchors on each side used as the local background, in bins.");
    optional.add({"--backgroundGap"})
        .type("int")
        .default_value(2)
        .help("Anchors next to the stripe body excluded from the background, in bins.");
    optional.add({"--obsExpThreshold"})
        .type("float")
        .default_value(1.5)
        .help("Minimum enrichment of the stripe body's mean obs/exp over the local "
              "background's mean.");
    optional.add({"--zScoreThreshold"})
        .type("float")
        .default_value(2.0)
        .help("Minimum preselection z-score against the local background.");
    optional.add({"--minRawCount"})
        .type("float")
        .default_value(1.0)
        .help("Minimum mean raw count over the stripe body.");
    optional.add({"--mergeWindow"})
        .type("int")
        .default_value(5)
        .help("Non-maximum suppression window, in bins.");
    optional.add({"--fdr"})
        .type("float")
        .default_value(0.05)
        .help("Benjamini-Hochberg q-value threshold for the final call set.");
    optional.add({"--threads", "-t"})
        .type("int")
        .default_value(4)
        .help("Number of threads to use.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.matrix = ns.str("matrix");
    args.out_file_name = ns.str("outFileName");
    if (ns.given("chromosomes")) {
        args.chromosomes = ns.strs("chromosomes");
    }
    args.min_stripe_length = ns.integer("minStripeLength");
    args.max_stripe_length = ns.integer("maxStripeLength");
    args.stripe_length_step = ns.integer("stripeLengthStep");
    args.background_window = ns.integer("backgroundWindow");
    args.background_gap = ns.integer("backgroundGap");
    args.obs_exp_threshold = ns.real("obsExpThreshold");
    args.z_score_threshold = ns.real("zScoreThreshold");
    args.min_raw_count = ns.real("minRawCount");
    args.merge_window = ns.integer("mergeWindow");
    args.fdr = ns.real("fdr");
    args.threads = static_cast<int>(ns.integer("threads"));
    return args;
}

struct ChromosomeMatrix {
    hicx::CsrMatrix matrix;  // triu, diagonal removed, Symmetry::Full
    std::vector<hicx::CutInterval> cut_intervals;
    std::int64_t bin_size = 0;
};

// The same triu-and-band cut hicDetectLoops.cpp uses, kept local because the
// two tools are independent (PLAN.md 5.0 item 3: implementation is free).
hicx::CsrMatrix upper_band(const hicx::CsrMatrix& matrix, std::int64_t distance_limit_bins) {
    hicx::CsrMatrix::Arrays arrays;
    arrays.rows = matrix.rows();
    arrays.cols = matrix.cols();
    arrays.dtype = matrix.dtype();
    arrays.symmetry = hicx::Symmetry::Full;
    arrays.indptr.assign(static_cast<std::size_t>(matrix.rows()) + 1, 0);

    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& values = matrix.data();
    for (std::int64_t row = 0; row < matrix.rows(); ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t column = static_cast<std::int64_t>(indices[k]);
            if (column <= row) {
                continue;
            }
            if (column - row > distance_limit_bins) {
                continue;
            }
            if (values[k] == 0.0) {
                continue;
            }
            arrays.indices.push_back(indices[k]);
            arrays.data.push_back(values[k]);
            ++arrays.indptr[static_cast<std::size_t>(row) + 1];
        }
    }
    for (std::size_t i = 1; i < arrays.indptr.size(); ++i) {
        arrays.indptr[i] += arrays.indptr[i - 1];
    }
    return hicx::CsrMatrix::adopt(std::move(arrays));
}

ChromosomeMatrix load_cool_chromosome(const hicx::CoolFile& cool,
                                      const std::vector<hicx::CutInterval>& bins,
                                      const std::string& chromosome,
                                      std::int64_t max_distance_bins) {
    std::int64_t first = -1;
    std::int64_t last = 0;
    for (std::size_t bin = 0; bin < bins.size(); ++bin) {
        if (bins[bin].chrom != chromosome) {
            continue;
        }
        if (first < 0) {
            first = static_cast<std::int64_t>(bin);
        }
        last = static_cast<std::int64_t>(bin) + 1;
    }
    if (first < 0) {
        throw std::runtime_error("Chromosome name not in matrix. '" + chromosome + "'");
    }

    hicx::CsrMatrix::Arrays arrays;
    arrays.rows = last - first;
    arrays.cols = last - first;
    arrays.symmetry = hicx::Symmetry::Full;
    arrays.indptr.assign(static_cast<std::size_t>(arrays.rows) + 1, 0);
    arrays.dtype = cool.count_dtype();

    cool.for_each_pixel_chunk(first, last, first, last, [&](const hicx::PixelChunk& chunk) {
        for (std::size_t k = 0; k < chunk.bin1.size(); ++k) {
            const std::int64_t row = chunk.bin1[k];
            const std::int64_t column = chunk.bin2[k];
            if (column - row > max_distance_bins) {
                continue;
            }
            const double value = chunk.count[k];
            if (value == 0.0) {
                continue;
            }
            arrays.indices.push_back(static_cast<std::int32_t>(column - first));
            arrays.data.push_back(value);
            ++arrays.indptr[static_cast<std::size_t>(row - first) + 1];
        }
    });
    for (std::size_t i = 1; i < arrays.indptr.size(); ++i) {
        arrays.indptr[i] += arrays.indptr[i - 1];
    }

    ChromosomeMatrix result;
    result.matrix = hicx::CsrMatrix::adopt(std::move(arrays));
    result.cut_intervals.assign(bins.begin() + static_cast<std::ptrdiff_t>(first),
                                bins.begin() + static_cast<std::ptrdiff_t>(last));
    result.bin_size = hicx::BinTable(result.cut_intervals).bin_size();
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);

    if (args.max_stripe_length < args.min_stripe_length) {
        std::fprintf(stderr,
                     "ERROR:hicexplorer.hicDetectStripes:--maxStripeLength must be at least "
                     "--minStripeLength\n");
        return 1;
    }
    if (args.stripe_length_step <= 0) {
        std::fprintf(stderr,
                     "ERROR:hicexplorer.hicDetectStripes:--stripeLengthStep must be positive\n");
        return 1;
    }

    const unsigned int workers = std::max(1U, static_cast<unsigned int>(std::max(1, args.threads)));
    const bool is_cooler = hicx::check_cooler(args.matrix);

    std::vector<std::string> chromosomes;
    hicx::MatrixData whole;
    std::optional<hicx::CoolFile> cool;
    std::vector<hicx::CutInterval> cool_bins;
    if (!is_cooler) {
        whole = hicx::read_hicexplorer_h5(args.matrix);
        whole.matrix.symmetrize_in_place();
        if (args.chromosomes.has_value()) {
            chromosomes = *args.chromosomes;
        } else {
            for (const auto& [name, range] : hicx::chrom_bin_boundaries(whole.cut_intervals)) {
                (void)range;
                chromosomes.push_back(name);
            }
        }
    } else {
        cool.emplace(args.matrix);
        cool_bins = cool->read_bins();
        if (args.chromosomes.has_value()) {
            chromosomes = *args.chromosomes;
        } else {
            chromosomes = cool->chrom_names();
        }
    }

    std::vector<hicx::stripes::Candidate> all_candidates;
    // The genomic mapping for every kept candidate, indexed the same way, so
    // the writer does not need to search the bin table again.
    std::vector<std::pair<const std::vector<hicx::CutInterval>*, std::int64_t>> origin;

    // The full per-chromosome context is kept only long enough to map the
    // surviving candidates; a chromosome's bands are freed before the next one
    // is loaded so the peak stays at one chromosome's dense arrays rather than
    // the whole genome's.
    struct ChromosomeContext {
        std::vector<hicx::CutInterval> cut_intervals;
        std::int64_t bin_size = 0;
    };
    std::vector<ChromosomeContext> contexts;
    contexts.reserve(chromosomes.size());

    for (const std::string& chromosome : chromosomes) {
        ChromosomeMatrix block;
        // A generous band: the longest candidate length plus the background
        // window and gap, so every candidate's background has room.
        const std::int64_t approx_bin_size_guess = 10000;  // refined once bin_size is known
        (void)approx_bin_size_guess;
        try {
            if (is_cooler) {
                // First pass to learn bin_size cheaply: read_bins already
                // gives cut_intervals, so bin_size is known before the band
                // width is decided.
                const std::int64_t bin_size_hint =
                    hicx::BinTable(cool_bins).bin_size() > 0 ? hicx::BinTable(cool_bins).bin_size()
                                                             : 10000;
                const std::int64_t max_distance_bins =
                    args.max_stripe_length / bin_size_hint +
                    args.background_gap + args.background_window + 1;
                block = load_cool_chromosome(*cool, cool_bins, chromosome, max_distance_bins);
            } else {
                std::vector<std::int64_t> selection;
                for (std::size_t bin = 0; bin < whole.cut_intervals.size(); ++bin) {
                    if (whole.cut_intervals[bin].chrom == chromosome) {
                        selection.push_back(static_cast<std::int64_t>(bin));
                    }
                }
                if (selection.empty()) {
                    throw std::runtime_error("Chromosome name not in matrix. '" + chromosome + "'");
                }
                block.matrix = hicx::select_bins(whole.matrix, selection);
                block.cut_intervals.reserve(selection.size());
                for (const std::int64_t bin : selection) {
                    block.cut_intervals.push_back(
                        whole.cut_intervals[static_cast<std::size_t>(bin)]);
                }
                block.bin_size = hicx::BinTable(block.cut_intervals).bin_size();
            }
        } catch (const std::exception& error) {
            std::fprintf(stderr, "ERROR:hicexplorer.hicDetectStripes:%s\n", error.what());
            return 1;
        }

        if (block.bin_size <= 0 || block.matrix.rows() < 5) {
            continue;
        }

        const std::int64_t bin_size = block.bin_size;
        const std::int64_t max_distance_bins =
            args.max_stripe_length / bin_size + args.background_gap + args.background_window + 1;

        hicx::CsrMatrix matrix;
        if (is_cooler) {
            matrix = upper_band(block.matrix, max_distance_bins);
        } else {
            matrix = upper_band(block.matrix, max_distance_bins);
        }
        block.matrix = hicx::CsrMatrix();

        if (matrix.stored_nnz() == 0) {
            continue;
        }

        hicx::stripes::DetectOptions options;
        options.background_window_bins = args.background_window;
        options.background_gap_bins = args.background_gap;
        options.min_obs_exp = args.obs_exp_threshold;
        options.preselect_z = args.z_score_threshold;
        options.min_raw_count = args.min_raw_count;
        options.merge_window_bins = args.merge_window;
        options.fdr_q = args.fdr;
        for (std::int64_t length_bp = args.min_stripe_length; length_bp <= args.max_stripe_length;
            length_bp += args.stripe_length_step) {
            const std::int64_t length_bins = std::max<std::int64_t>(1, length_bp / bin_size);
            if (options.length_grid_bins.empty() || options.length_grid_bins.back() != length_bins) {
                options.length_grid_bins.push_back(length_bins);
            }
        }
        if (options.length_grid_bins.empty()) {
            continue;
        }

        const hicx::stripes::Band horizontal =
            hicx::stripes::build_horizontal_band(matrix, matrix.rows(), max_distance_bins);
        const hicx::stripes::Band vertical = hicx::stripes::build_vertical_band(horizontal);
        const hicx::stripes::RunningSums horizontal_sums = hicx::stripes::build_running_sums(horizontal);
        const hicx::stripes::RunningSums vertical_sums = hicx::stripes::build_running_sums(vertical);

        std::vector<hicx::stripes::Candidate> candidates = hicx::stripes::preselect(
            horizontal, vertical, horizontal_sums, vertical_sums, options, workers);
        candidates = hicx::stripes::suppress_non_maximal(std::move(candidates), args.merge_window);
        hicx::stripes::compute_pvalues(candidates, horizontal, vertical, options, workers);

        contexts.push_back(ChromosomeContext{block.cut_intervals, bin_size});
        const std::size_t context_index = contexts.size() - 1;
        for (auto& candidate : candidates) {
            all_candidates.push_back(candidate);
            origin.emplace_back(&contexts[context_index].cut_intervals, context_index);
        }
    }

    // Genome-wide Benjamini-Hochberg, done here (not through
    // hicx::stripes::apply_fdr) so that the surviving candidates keep their
    // index into `origin`, which apply_fdr's plain Candidate vector would
    // lose.
    std::vector<double> pvalues;
    pvalues.reserve(all_candidates.size());
    for (const hicx::stripes::Candidate& candidate : all_candidates) {
        pvalues.push_back(candidate.pvalue);
    }
    const std::vector<double> adjusted = hicx::stats::benjamini_hochberg_adjusted(pvalues);

    std::ofstream output(args.out_file_name);
    if (!output) {
        std::fprintf(stderr, "hicDetectStripes: error: cannot write %s\n",
                     args.out_file_name.c_str());
        return 1;
    }
    std::size_t kept_count = 0;
    for (std::size_t found = 0; found < all_candidates.size(); ++found) {
        if (adjusted[found] > args.fdr) {
            continue;
        }
        hicx::stripes::Candidate candidate = all_candidates[found];
        candidate.qvalue = adjusted[found];
        ++kept_count;
        const auto& [intervals_ptr, context_index] = origin[found];
        (void)context_index;
        const std::vector<hicx::CutInterval>& intervals = *intervals_ptr;
        const hicx::CutInterval& anchor_interval =
            intervals[static_cast<std::size_t>(candidate.anchor)];
        std::int64_t extent_start = 0;
        std::int64_t extent_end = 0;
        if (!candidate.vertical) {
            const std::int64_t far_bin = candidate.anchor + candidate.length_bins;
            const std::int64_t clamped =
                std::min<std::int64_t>(far_bin, static_cast<std::int64_t>(intervals.size()) - 1);
            extent_start = anchor_interval.start;
            extent_end = intervals[static_cast<std::size_t>(clamped)].end;
        } else {
            const std::int64_t far_bin = candidate.anchor - candidate.length_bins;
            const std::int64_t clamped = std::max<std::int64_t>(far_bin, 0);
            extent_start = intervals[static_cast<std::size_t>(clamped)].start;
            extent_end = anchor_interval.end;
        }
        output << anchor_interval.chrom << '\t' << anchor_interval.start << '\t'
              << anchor_interval.end << '\t' << (candidate.vertical ? "vertical" : "horizontal")
              << '\t' << extent_start << '\t' << extent_end << '\t' << candidate.enrichment << '\t'
              << candidate.pvalue << '\t' << candidate.qvalue << '\n';
    }
    std::fprintf(stderr, "INFO:hicexplorer.hicDetectStripes:Number of detected stripes: %zu\n",
                kept_count);
    hicx::report_resource_usage("hicDetectStripes");
    return 0;
}
