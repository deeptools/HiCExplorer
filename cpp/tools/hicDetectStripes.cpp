// hicDetectStripes: architectural stripe detection, PLAN.md tier 9 section
// 9.3. HiCExplorer 3.7 has no stripe caller, and this port is a faithful
// reimplementation of Stripenn 1.1.65.22's own detection method (Yoon et al.
// 2022), not an invented stand-in -- see cpp/tools/stripes_impl.hpp for the
// full rationale and the two classes of documented simplification (the
// Canny threshold heuristic, and the background/null sampling scheme).
//
// This file is the argument parser, the per-chromosome loader (the same
// band loader hicDetectLoops.cpp uses), the frame scan that drives
// stripes_impl's StripeSearch port, the background model, and the writer.
//
// Output format (new, no Python precedent): tab separated,
//
//   chrom  anchor_start  anchor_end  orientation  extent_start  extent_end
//   enrichment  pvalue  qvalue
//
// `orientation` is "horizontal" when the candidate's narrow axis sits near
// the *start* of its long axis (pos1 close to pos3, one bin apart at most --
// stripenn.getStripe.pvalue's own x1==y1 "downward" case, generalised from
// an exact match to a distance comparison because a real box's narrow axis
// is never literally the same bin as its long axis's start: the closest
// pixel to the diagonal is one bin away, not zero) and "vertical" when it
// sits near the *end* (pos2 close to pos4, Stripenn's x2==y2 "upward" case).
// This mapping, and which of Stripenn's two cases is "horizontal" versus
// "vertical" here, is chosen to match cpp/scripts/stripe_calibration.py's
// own planting convention exactly (plant_one() builds a "horizontal" plant
// as (anchor, anchor + d) for d = 1..length -- narrow axis at the start --
// and a "vertical" one as (anchor - d, anchor) -- narrow axis at the end).
// `enrichment` is the box's mean raw count divided by the chromosome's
// expected count at its offset from the diagonal.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include "hicx/adjust_ops.hpp"
#include "hicx/argparse.hpp"
#include "hicx/bins.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/matrix_ops.hpp"
#include "hicx/parallel.hpp"
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
    "                        [--maxWidth MAXWIDTH] [--canny CANNY]\n"
    "                        [--blurFilter BLURFILTER]\n"
    "                        [--maxPixelPercentiles P [P ...]]\n"
    "                        [--backgroundSamples BACKGROUNDSAMPLES]\n"
    "                        [--pValue PVALUE] [--fdr FDR] [--seed SEED]\n"
    "                        [--threads THREADS] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Detects architectural stripes on the given contact matrix. New in\n"
    "HiCExplorer v4 (PLAN.md tier 9 section 9.3): a faithful C++\n"
    "reimplementation of Stripenn 1.1.65.22's detection method (Yoon et al.\n"
    "2022), not an invented stand-in.\n"
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
    "                        Shortest candidate stripe length, in base pairs\n"
    "                        (Stripenn's minL, converted to bins).\n"
    "                        (Default: 100000).\n"
    "  --maxWidth MAXWIDTH   Maximum stripe width, in bins (Stripenn's maxW).\n"
    "                        (Default: 8).\n"
    "  --canny CANNY         Canny edge detection sigma (Stripenn's canny).\n"
    "                        (Default: 2.0).\n"
    "  --blurFilter BLURFILTER\n"
    "                        Mean filter size, an odd number (Stripenn's\n"
    "                        bfilter). (Default: 3).\n"
    "  --maxPixelPercentiles P [P ...]\n"
    "                        Percentiles of the contact frequency data to\n"
    "                        saturate the image (Stripenn's maxpixel).\n"
    "                        (Default: 0.95 0.96 0.97 0.98 0.99).\n"
    "  --backgroundSamples BACKGROUNDSAMPLES\n"
    "                        Random samples per chromosome for the background\n"
    "                        model the p-value is ranked against.\n"
    "                        (Default: 200000).\n"
    "  --pValue PVALUE       Raw p-value cutoff applied before FDR\n"
    "                        (Stripenn's pvalue). (Default: 0.1).\n"
    "  --fdr FDR             Benjamini-Hochberg q-value threshold for the final\n"
    "                        call set. (Default: 0.05).\n"
    "  --seed SEED           Seed for the background sample. (Default: 20260915).\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads to use. (Default: 4).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::string matrix;
    std::string out_file_name;
    std::optional<std::vector<std::string>> chromosomes;
    std::int64_t min_stripe_length = 100000;
    std::int64_t max_width = 8;
    double canny_sigma = 2.0;
    std::int64_t blur_filter = 3;
    std::vector<double> maxpixel_percentiles{0.95, 0.96, 0.97, 0.98, 0.99};
    std::int64_t background_samples = 200000;
    double p_value = 0.1;
    double fdr = 0.05;
    std::int64_t seed = 20260915;
    int threads = 4;
};

Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicDetectStripes",
                       "Detects architectural stripes on the given contact matrix.");
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
    optional.add({"--maxWidth"})
        .type("int")
        .default_value(8)
        .help("Maximum stripe width, in bins.");
    optional.add({"--canny"})
        .type("float")
        .default_value(2.0)
        .help("Canny edge detection sigma.");
    optional.add({"--blurFilter"})
        .type("int")
        .default_value(3)
        .help("Mean filter size, an odd number.");
    optional.add({"--maxPixelPercentiles"})
        .nargs("+")
        .type("float")
        .help("Percentiles of the contact frequency data to saturate the image. "
              "(Default: 0.95 0.96 0.97 0.98 0.99).");
    optional.add({"--backgroundSamples"})
        .type("int")
        .default_value(200000)
        .help("Random samples per chromosome for the background model.");
    optional.add({"--pValue"})
        .type("float")
        .default_value(0.1)
        .help("Raw p-value cutoff applied before FDR.");
    optional.add({"--fdr"})
        .type("float")
        .default_value(0.05)
        .help("Benjamini-Hochberg q-value threshold for the final call set.");
    optional.add({"--seed"})
        .type("int")
        .default_value(20260915)
        .help("Seed for the background sample.");
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
    args.max_width = ns.integer("maxWidth");
    args.canny_sigma = ns.real("canny");
    args.blur_filter = ns.integer("blurFilter");
    if (ns.given("maxPixelPercentiles")) {
        args.maxpixel_percentiles = ns.reals("maxPixelPercentiles");
    }
    args.background_samples = ns.integer("backgroundSamples");
    args.p_value = ns.real("pValue");
    args.fdr = ns.real("fdr");
    args.seed = ns.integer("seed");
    args.threads = static_cast<int>(ns.integer("threads"));
    return args;
}

struct ChromosomeMatrix {
    hicx::CsrMatrix matrix;
    std::vector<hicx::CutInterval> cut_intervals;
    std::int64_t bin_size = 0;
};

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
        const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t column = static_cast<std::int64_t>(indices[k]);
            if (column <= row || column - row > distance_limit_bins || values[k] == 0.0) {
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
                                      const std::string& chromosome, std::int64_t max_distance_bins) {
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

// A simple (anchor, distance) -> raw count lookup over the chromosome's band,
// used only for the background sample and the enrichment denominator, kept
// local to this file (stripes_impl.hpp's own API is the detection engine).
struct AnchorDistanceBand {
    std::int64_t n_bins = 0;
    std::int64_t max_distance = 0;
    std::vector<double> raw;  // [a * max_distance + (d - 1)]
    std::vector<double> expected;  // per distance, mean raw count

    [[nodiscard]] double at(std::int64_t anchor, std::int64_t distance) const {
        if (anchor < 0 || anchor >= n_bins || distance < 1 || distance > max_distance) {
            return 0.0;
        }
        return raw[static_cast<std::size_t>(anchor * max_distance + (distance - 1))];
    }

    // A single matrix cell by absolute (row, column), symmetric: only the
    // upper triangle (row < col) is actually stored, so a query with
    // row > col reads the mirrored entry, and the diagonal (row == col,
    // which this band never stores self-counts for) reads 0.
    [[nodiscard]] double cell(std::int64_t row, std::int64_t col) const {
        if (row == col) {
            return 0.0;
        }
        return row < col ? at(row, col - row) : at(col, row - col);
    }

    // The mean over a genuine rectangular row x column block, [r0, r1) x
    // [c0, c1), out-of-range rows/columns contributing 0 (matching Stripenn
    // treating anything past the chromosome's own bins, via nantozero, as
    // zero rather than skipping it, since hicDetectStripes.cpp's own
    // stripenn.getStripe.nulldist port below clamps its own row and column
    // ranges to stay in bounds in the same places nulldist does). This is
    // stripenn.getStripe.nulldist's own tableft_up/tabcenter_up/tabright_up
    // (and pvalue's mat_center/mat_left/mat_right): a real 2-D area, not a
    // window in (anchor, distance) space, which an earlier version of this
    // background model used and which does not correspond to any rectangle
    // Stripenn actually averages over.
    [[nodiscard]] double block_mean(std::int64_t r0, std::int64_t r1, std::int64_t c0,
                                    std::int64_t c1) const {
        double sum = 0.0;
        std::int64_t count = 0;
        for (std::int64_t r = r0; r < r1; ++r) {
            for (std::int64_t c = c0; c < c1; ++c) {
                if (r >= 0 && r < n_bins && c >= 0 && c < n_bins) {
                    sum += cell(r, c);
                }
                ++count;
            }
        }
        return count > 0 ? sum / static_cast<double>(count) : 0.0;
    }
};

AnchorDistanceBand build_anchor_distance_band(const hicx::CsrMatrix& upper, std::int64_t n_bins,
                                              std::int64_t max_distance) {
    AnchorDistanceBand band;
    band.n_bins = n_bins;
    band.max_distance = max_distance;
    band.raw.assign(static_cast<std::size_t>(n_bins * max_distance), 0.0);
    const std::vector<std::int64_t>& indptr = upper.indptr();
    const std::vector<std::int32_t>& indices = upper.indices();
    const std::vector<double>& values = upper.data();
    for (std::int64_t row = 0; row < upper.rows() && row < n_bins; ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t col = static_cast<std::int64_t>(indices[k]);
            const std::int64_t d = col - row;
            if (d >= 1 && d <= max_distance) {
                band.raw[static_cast<std::size_t>(row * max_distance + (d - 1))] = values[k];
            }
        }
    }
    band.expected.assign(static_cast<std::size_t>(max_distance), 0.0);
    for (std::int64_t d = 1; d <= max_distance; ++d) {
        const std::int64_t positions = n_bins - d;
        if (positions <= 0) {
            continue;
        }
        double sum = 0.0;
        for (std::int64_t a = 0; a < positions; ++a) {
            sum += band.raw[static_cast<std::size_t>(a * max_distance + (d - 1))];
        }
        band.expected[static_cast<std::size_t>(d - 1)] = sum / static_cast<double>(positions);
    }
    return band;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
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
        chromosomes = args.chromosomes.has_value() ? *args.chromosomes : cool->chrom_names();
    }

    // The background/null model needs distance coverage to at least 399
    // bins (Stripenn's own table size) plus its background_size margin; the
    // frame construction (200-bin steps, ~100-bin overlap) needs full
    // pairwise coverage up to about 400 bins too, so one band width serves
    // both.
    constexpr std::int64_t kBackgroundBins = 400;

    std::vector<hicx::stripes::Candidate> all_candidates;
    struct ChromosomeContext {
        std::vector<hicx::CutInterval> cut_intervals;
        std::int64_t bin_size = 0;
    };
    std::vector<ChromosomeContext> contexts;
    contexts.reserve(chromosomes.size());
    std::mt19937_64 rng(static_cast<std::uint64_t>(args.seed));

    for (const std::string& chromosome : chromosomes) {
        ChromosomeMatrix block;
        try {
            if (is_cooler) {
                const std::int64_t bin_size_hint =
                    hicx::BinTable(cool_bins).bin_size() > 0 ? hicx::BinTable(cool_bins).bin_size() : 10000;
                const std::int64_t background_margin = 50000 / bin_size_hint + 1;
                const std::int64_t max_distance_bins = kBackgroundBins + background_margin + 10;
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
                    block.cut_intervals.push_back(whole.cut_intervals[static_cast<std::size_t>(bin)]);
                }
                block.bin_size = hicx::BinTable(block.cut_intervals).bin_size();
            }
        } catch (const std::exception& error) {
            std::fprintf(stderr, "ERROR:hicexplorer.hicDetectStripes:%s\n", error.what());
            return 1;
        }
        if (block.bin_size <= 0 || block.matrix.rows() < 10) {
            continue;
        }

        const std::int64_t bin_size = block.bin_size;
        const std::int64_t background_margin = 50000 / bin_size + 1;
        const std::int64_t max_distance_bins = kBackgroundBins + background_margin + 10;
        const std::int64_t n_bins = block.matrix.rows();

        const hicx::CsrMatrix upper = upper_band(block.matrix, max_distance_bins);
        block.matrix = hicx::CsrMatrix();
        if (upper.stored_nnz() == 0) {
            continue;
        }
        const AnchorDistanceBand band = build_anchor_distance_band(upper, n_bins, max_distance_bins);

        // The maxpixel quantiles, over the band's nonzero raw counts (a
        // documented simplification of getQuantile_original, which uses the
        // whole chromosome's nonzero counts; see stripes_impl.hpp).
        std::vector<double> nonzero;
        nonzero.reserve(band.raw.size());
        for (const double v : band.raw) {
            if (v > 0.0) {
                nonzero.push_back(v);
            }
        }
        std::vector<double> maxpixel_values;
        for (const double p : args.maxpixel_percentiles) {
            maxpixel_values.push_back(hicx::stripes::quantile(nonzero, p));
        }

        // The frame scan: 200-bin steps, +-100-bin overlap (stripenn.getStripe.extract).
        const std::int64_t nframes = (n_bins + 199) / 200;
        const std::int64_t min_length_bins = std::max<std::int64_t>(1, args.min_stripe_length / bin_size);

        std::vector<std::vector<hicx::stripes::Candidate>> per_frame(static_cast<std::size_t>(nframes));
        hicx::parallel_for(static_cast<std::size_t>(nframes), workers, [&](std::size_t frame_index) {
            const auto idx = static_cast<std::int64_t>(frame_index);
            std::int64_t start = idx * 200 - 100;
            std::int64_t end = (idx + 1) * 200 + 99;
            if (idx == 0) {
                start = 0;
            }
            if (end >= n_bins) {
                end = n_bins - 1;
            }
            if (end <= start) {
                return;
            }
            const int S = static_cast<int>(end - start + 1);
            hicx::stripes::Image submat(S, S, 0.0);
            for (int i = 0; i < S; ++i) {
                for (int j = i + 1; j < S; ++j) {
                    const std::int64_t a = start + i;
                    const std::int64_t b = start + j;
                    const double v = band.at(a, b - a);
                    submat.at(i, j) = v;
                    submat.at(j, i) = v;
                }
            }
            std::vector<hicx::stripes::Candidate> frame_candidates;
            for (const double M : maxpixel_values) {
                const std::vector<hicx::stripes::FrameCandidate> found = hicx::stripes::stripe_search_frame(
                    submat, M, args.canny_sigma, static_cast<int>(min_length_bins),
                    static_cast<int>(args.max_width), static_cast<int>(args.blur_filter));
                for (const hicx::stripes::FrameCandidate& fc : found) {
                    const std::int64_t abs_x0 = start + fc.x;
                    const std::int64_t abs_x1 = start + fc.x + fc.w - 1;
                    const std::int64_t abs_y0 = start + fc.y;
                    const std::int64_t abs_y1 = start + fc.y + fc.h - 1;
                    if (abs_x1 >= n_bins || abs_y1 >= n_bins) {
                        continue;
                    }
                    hicx::stripes::Candidate candidate;
                    candidate.chrom = chromosome;
                    candidate.pos1 = abs_x0;  // kept as bin indices until genomic mapping below
                    candidate.pos2 = abs_x1;
                    candidate.pos3 = abs_y0;
                    candidate.pos4 = abs_y1;
                    candidate.frame_index = static_cast<int>(idx);
                    const auto box_mean = [&](int col0) {
                        if (col0 < 0 || col0 + fc.w > S) {
                            return -1.0;  // out of frame: no usable shifted box
                        }
                        double sum = 0.0;
                        for (int r = fc.y; r < fc.y + fc.h; ++r) {
                            for (int c = col0; c < col0 + fc.w; ++c) {
                                sum += submat.at(r, c);
                            }
                        }
                        return sum / static_cast<double>(fc.h) / static_cast<double>(fc.w);
                    };
                    candidate.mean = box_mean(fc.x);
                    const int background_bins = static_cast<int>(std::max<std::int64_t>(1, 50000 / bin_size));
                    candidate.left_mean = box_mean(fc.x - background_bins);
                    candidate.right_mean = box_mean(fc.x + background_bins);
                    frame_candidates.push_back(candidate);
                }
            }
            per_frame[frame_index] = std::move(frame_candidates);
        });

        std::vector<hicx::stripes::Candidate> chrom_candidates;
        for (auto& frame_candidates : per_frame) {
            for (auto& candidate : frame_candidates) {
                chrom_candidates.push_back(std::move(candidate));
            }
        }
        chrom_candidates = hicx::stripes::remove_redundant(std::move(chrom_candidates), false, workers);

        // Background model: two direction-specific tables, matching
        // stripenn.getStripe.nulldist's own bgleft_down/bgright_down (x1==y1,
        // "horizontal" here) and bgleft_up/bgright_up (x2==y2, "vertical"),
        // each built from genuine background_size x background_size
        // rectangular windows -- nulldist's own tableft_up/tabcenter_up/
        // tabright_up shape exactly, read from stripenn/getStripe.py again
        // to get this right, not from memory of the earlier read. What
        // still differs from nulldist(): the anchors sampled here are drawn
        // uniformly from this one chromosome, not from nulldist's
        // genome-wide sample stratified by chromosome size and by each
        // chromosome's own non-empty-column density (still the simplified,
        // documented part of this port -- see stripes_impl.hpp).
        const std::int64_t background_size = std::max<std::int64_t>(1, 50000 / bin_size);
        const std::int64_t background_up = background_size / 2;
        const std::int64_t background_down = background_size - background_up;
        hicx::stripes::BackgroundModel model_down;  // x1==y1, "horizontal"
        hicx::stripes::BackgroundModel model_up;    // x2==y2, "vertical"
        model_down.left.assign(kBackgroundBins, {});
        model_down.right.assign(kBackgroundBins, {});
        model_up.left.assign(kBackgroundBins, {});
        model_up.right.assign(kBackgroundBins, {});
        const std::int64_t margin = background_up + background_size + kBackgroundBins;
        std::uniform_int_distribution<std::int64_t> anchor_dist(
            margin, std::max<std::int64_t>(margin, n_bins - margin - 1));
        std::uniform_int_distribution<std::int64_t> distance_dist(1, kBackgroundBins);
        if (n_bins > 2 * margin) {
            for (std::int64_t sample = 0; sample < args.background_samples; ++sample) {
                const std::int64_t x = anchor_dist(rng);
                const std::int64_t j = distance_dist(rng);
                const std::size_t d = static_cast<std::size_t>(j - 1);
                for (const bool down : {true, false}) {
                    const std::int64_t y = down ? x + j : x - j;
                    const double center = band.block_mean(x - background_up, x + background_down,
                                                          y - background_up, y + background_down);
                    const double left =
                        band.block_mean(x - background_up - background_size, x - background_up,
                                        y - background_up, y + background_down);
                    const double right =
                        band.block_mean(x + background_down, x + background_down + background_size,
                                        y - background_up, y + background_down);
                    hicx::stripes::BackgroundModel& target = down ? model_down : model_up;
                    target.left[d].push_back(center - left);
                    target.right[d].push_back(center - right);
                }
            }
        }

        // Genomic mapping and the p-value: stripenn.getStripe.pvalue's own
        // per-row test (one center/left/right comparison per row of the
        // candidate's long axis, ranked against the matching direction's
        // background model at that row's own distance from the anchor end)
        // and its median aggregation across rows, not a single whole-box
        // comparison against a single distance bucket.
        for (hicx::stripes::Candidate& candidate : chrom_candidates) {
            const std::int64_t x0 = candidate.pos1;
            const std::int64_t x1 = candidate.pos2;
            const std::int64_t y0 = candidate.pos3;
            const std::int64_t y1 = candidate.pos4;
            const std::int64_t d_start = std::llabs(x0 - y0);
            const std::int64_t d_end = std::llabs(x1 - y1);
            const bool vertical = d_end < d_start;
            const hicx::stripes::BackgroundModel& model = vertical ? model_up : model_down;
            std::vector<double> row_pvalues;
            row_pvalues.reserve(static_cast<std::size_t>(y1 - y0 + 1));
            for (std::int64_t row = y0; row <= y1; ++row) {
                const double center = band.block_mean(row, row + 1, x0, x1 + 1);
                const double left = band.block_mean(row, row + 1, x0 - background_size, x0);
                const double right = band.block_mean(row, row + 1, x1 + 1, x1 + 1 + background_size);
                const std::int64_t distance_from_anchor = vertical ? (y1 - row) : (row - y0);
                const int clamped = static_cast<int>(
                    std::clamp<std::int64_t>(distance_from_anchor, 0, kBackgroundBins - 1));
                row_pvalues.push_back(
                    hicx::stripes::candidate_pvalue(model, clamped, center - left, center - right));
            }
            std::sort(row_pvalues.begin(), row_pvalues.end());
            candidate.pvalue = row_pvalues.empty() ? 1.0 : row_pvalues[row_pvalues.size() / 2];

            const std::int64_t distance_bin = std::max<std::int64_t>(1, y0 - x0);
            const hicx::CutInterval& xc0 = block.cut_intervals[static_cast<std::size_t>(x0)];
            const hicx::CutInterval& xc1 = block.cut_intervals[static_cast<std::size_t>(x1)];
            const hicx::CutInterval& yc0 = block.cut_intervals[static_cast<std::size_t>(y0)];
            const hicx::CutInterval& yc1 = block.cut_intervals[static_cast<std::size_t>(y1)];
            const double expected = band.expected[static_cast<std::size_t>(
                std::clamp<std::int64_t>(distance_bin, 1, max_distance_bins) - 1)];
            candidate.mean = expected > 0.0 ? candidate.mean / expected : 0.0;  // now holds enrichment
            candidate.pos1 = xc0.start;
            candidate.pos2 = xc1.end;
            candidate.pos3 = yc0.start;
            candidate.pos4 = yc1.end;
        }
        chrom_candidates.erase(std::remove_if(chrom_candidates.begin(), chrom_candidates.end(),
                                              [&](const hicx::stripes::Candidate& c) {
                                                  return c.pvalue >= args.p_value;
                                              }),
                               chrom_candidates.end());

        contexts.push_back(ChromosomeContext{block.cut_intervals, bin_size});
        for (auto& candidate : chrom_candidates) {
            all_candidates.push_back(std::move(candidate));
        }
    }

    all_candidates =
        hicx::stripes::remove_redundant(std::move(all_candidates), true, workers);

    std::vector<double> pvalues;
    pvalues.reserve(all_candidates.size());
    for (const auto& c : all_candidates) {
        pvalues.push_back(c.pvalue);
    }
    const std::vector<double> adjusted = hicx::stats::benjamini_hochberg_adjusted(pvalues);

    std::ofstream output(args.out_file_name);
    if (!output) {
        std::fprintf(stderr, "hicDetectStripes: error: cannot write %s\n", args.out_file_name.c_str());
        return 1;
    }
    std::size_t kept_count = 0;
    for (std::size_t i = 0; i < all_candidates.size(); ++i) {
        if (adjusted[i] > args.fdr) {
            continue;
        }
        const hicx::stripes::Candidate& c = all_candidates[i];
        // Which end of the box sits nearer the diagonal (stripenn.getStripe.
        // pvalue's own "downward" x1==y1 / "upward" x2==y2 distinction,
        // generalised to a distance rather than an exact bin match: a
        // symmetric Hi-C matrix means a box never literally touches the
        // diagonal at bin resolution, since the pixel closest to it is one
        // bin away, not zero, so pos1==pos3 essentially never holds at
        // real genomic coordinates -- verified this was firing on almost
        // nothing on real data, which is what caught the bug). "vertical"
        // is the end tag cpp/scripts/stripe_calibration.py's plant_one()
        // uses for a box built as (anchor - d, anchor), d = 1..length: its
        // narrow axis (pos1/pos2) sits at the *end* of its long axis
        // (pos4 is one bin from pos1/pos2, pos3 is far), so d_end is the
        // small one. "horizontal" is plant_one()'s (anchor, anchor + d):
        // its narrow axis sits at the *start* of the long axis, so d_start
        // is the small one.
        const std::int64_t d_start = std::llabs(c.pos1 - c.pos3);
        const std::int64_t d_end = std::llabs(c.pos2 - c.pos4);
        const bool vertical = d_end < d_start;
        ++kept_count;
        output << c.chrom << '\t' << c.pos1 << '\t' << c.pos2 << '\t' << (vertical ? "vertical" : "horizontal")
              << '\t' << std::min(c.pos1, c.pos3) << '\t' << std::max(c.pos2, c.pos4) << '\t' << c.mean << '\t'
              << c.pvalue << '\t' << adjusted[i] << '\n';
    }
    std::fprintf(stderr, "INFO:hicexplorer.hicDetectStripes:Number of detected stripes: %zu\n", kept_count);
    hicx::report_resource_usage("hicDetectStripes");
    return 0;
}
