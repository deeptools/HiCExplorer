// hicDetectStripes: architectural stripe detection by the Stripenn method
// (Yoon S, Vo C, Chandra V, et al. Stripenn detects architectural stripes from
// chromatin conformation data using computer vision. Nature Communications
// 2022; Stripenn 1.1.65.22). HiCExplorer 3.7 has no stripe caller. This tool
// implements Stripenn's own algorithm stage by stage (cpp/tools/stripes_impl.*
// and stripes_numpy.*), and its acceptance gate is identical output to
// Stripenn on the same matrix: the same candidate stripes, p-values and
// stripiness scores, with Stripenn's defaults (cpp/scripts/stripenn_oracle
// compares every stage with Stripenn itself).
//
// Reading the matrix: balanced values come through the C++ matrix layer
// (hicx::for_each_corrected_pixel_chunk for a cooler), exactly as
// cooler.matrix(balance=...) applies the weights. Stripenn's own selection of
// chromosomes is kept: the matrix's chromosomes without the unplaced contigs,
// chrM and chrY, and its null distribution is estimated from all of them
// whichever chromosomes were requested.
//
// Output (tab separated, one row per called stripe, Stripenn's
// result_filtered order: descending stripiness):
//
//   chrom  anchor_start  anchor_end  orientation  extent_start  extent_end
//   enrichment  pvalue  qvalue  stripiness  mean  maxpixel
//
// anchor_* is the stripe's narrow axis (Stripenn pos1 - 1, pos2), extent_*
// the span of its box along the diagonal axis, `orientation` is "horizontal"
// when the box starts on the diagonal (Stripenn's downward stripe, pos1 ==
// pos3) and "vertical" otherwise (the upward stripe). `enrichment` is the
// mean observed over expected value of the box, `pvalue` and `stripiness`
// are Stripenn's, `mean` and `maxpixel` its Mean and maxpixel columns.
// --unfilteredFileName writes Stripenn's result_unfiltered.tsv (all
// candidates, its columns). False-positive control follows Stripenn: a raw
// p-value cutoff (--pValue, default 0.1) and no multiple-testing correction;
// --fdr is an option, off by default.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <optional>
#include <sstream>
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
#include "stripes_npy.hpp"

namespace {

namespace st = hicx::stripes;

const char* const kUsage =
    "usage: hicDetectStripes --matrix MATRIX --outFileName OUTFILENAME\n"
    "                        [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]\n"
    "                        [--norm NORM] [--canny CANNY] [--minL MINL]\n"
    "                        [--maxWidth MAXWIDTH]\n"
    "                        [--maxPixelPercentiles P [P ...]]\n"
    "                        [--pValue PVALUE] [--mask MASK]\n"
    "                        [--blurFilter BLURFILTER] [--seed SEED]\n"
    "                        [--fdr FDR] [--unfilteredFileName UNFILTERED]\n"
    "                        [--threads THREADS] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Detects architectural stripes on the given contact matrix with the\n"
    "method of Stripenn (Yoon S, Vo C, Chandra V, et al., Stripenn detects\n"
    "architectural stripes from chromatin conformation data using computer\n"
    "vision, Nature Communications 2022; Stripenn 1.1.65.22). HiCExplorer 3.7\n"
    "has no stripe caller; this is Stripenn's algorithm in C++, with its\n"
    "defaults, and gives the same stripes, p-values and stripiness scores as\n"
    "Stripenn on the same matrix.\n"
    "\n"
    "The contact image of each chromosome is saturated at several percentiles\n"
    "of its contact values (--maxPixelPercentiles) and scanned in overlapping\n"
    "frames of 200 bins. In every frame Canny edge detection (--canny, the\n"
    "sigma of its Gaussian), a vertical line filter and a run-length scan of\n"
    "the edge columns (--minL) find pairs of edges at most --maxWidth bins\n"
    "apart, which are the sides of a candidate stripe. Overlapping candidates\n"
    "are reduced to the more elongated one. Each candidate gets an empirical\n"
    "p-value against a background distribution that is sampled, with the seeded\n"
    "random generator (--seed), from all chromosomes of the matrix, and a\n"
    "stripiness score from its observed over expected contacts. A stripe is\n"
    "called when its p-value is below --pValue.\n"
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        The matrix to compute the stripe detection on\n"
    "                        (cool, mcool or h5).\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        Outfile name to store the called stripes\n"
    "                        (tab separated, Stripenn's result_filtered\n"
    "                        order).\n"
    "\n"
    "Optional arguments:\n"
    "  --chromosomes CHROMOSOMES [CHROMOSOMES ...]\n"
    "                        Chromosomes to search for stripes (Stripenn's\n"
    "                        chrom). If not set, all chromosomes are used.\n"
    "                        The background distribution always uses all\n"
    "                        chromosomes of the matrix.\n"
    "  --norm NORM           Normalisation: 'weight' (the balancing weights\n"
    "                        of the cool file), the name of another column of\n"
    "                        its bin table, or 'None' for the raw counts\n"
    "                        (Stripenn's norm). (Default: weight).\n"
    "  --canny CANNY         Canny edge detection parameter, the sigma of its\n"
    "                        Gaussian filter (Stripenn's canny).\n"
    "                        (Default: 2.0).\n"
    "  --minL MINL           Minimum length of a stripe in bins (Stripenn's\n"
    "                        minL). (Default: 10).\n"
    "  --minStripeLength MINSTRIPELENGTH\n"
    "                        The same in base pairs; when given it replaces\n"
    "                        --minL by MINSTRIPELENGTH divided by the bin\n"
    "                        size.\n"
    "  --maxWidth MAXWIDTH   Maximum width of a stripe in bins (Stripenn's\n"
    "                        maxW). (Default: 8).\n"
    "  --maxPixelPercentiles P [P ...]\n"
    "                        Percentiles of the chromosome's contact values\n"
    "                        at which the image is saturated (Stripenn's\n"
    "                        maxpixel). (Default: 0.95 0.96 0.97 0.98 0.99).\n"
    "  --pValue PVALUE       P-value cutoff for a called stripe (Stripenn's\n"
    "                        pvalue). No multiple-testing correction is\n"
    "                        applied, as in Stripenn. (Default: 0.1).\n"
    "  --mask MASK           Column coordinates to be masked in the\n"
    "                        stripiness score, e.g. chr9:12345678-12345789\n"
    "                        (Stripenn's mask). (Default: 0, none).\n"
    "  --blurFilter BLURFILTER\n"
    "                        Mean filter size, an odd number (Stripenn's\n"
    "                        bfilter). (Default: 3).\n"
    "  --seed SEED           Seed of the pseudo random number generator that\n"
    "                        samples the background (Stripenn's seed, Python's\n"
    "                        random module). (Default: 123456789).\n"
    "  --sequentialRandom    Draw the background samples of all chromosomes\n"
    "                        from one generator in sequence. Stripenn does so\n"
    "                        when it runs on a single core; with several cores\n"
    "                        (and here by default) every chromosome starts\n"
    "                        from the seeded initial state.\n"
    "  --fdr FDR             Benjamini-Hochberg q-value threshold, applied\n"
    "                        instead of --pValue when given. Off by default\n"
    "                        (Stripenn applies none).\n"
    "  --unfilteredFileName UNFILTEREDFILENAME\n"
    "                        Also write all candidate stripes with Stripenn's\n"
    "                        result_unfiltered.tsv columns.\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads to use. (Default: 4).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::string matrix;
    std::string out_file_name;
    std::optional<std::vector<std::string>> chromosomes;
    std::string norm = "weight";
    st::Parameters params;
    std::optional<std::int64_t> min_stripe_length;
    std::optional<double> fdr;
    std::string unfiltered_file_name;
    std::string dump_stages;
    int threads = 4;
};

Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicDetectStripes",
                       "Detects architectural stripes on the given contact matrix with the Stripenn method.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrix", "-m"})
        .required()
        .input({"h5", "cool", "mcool"})
        .help("The matrix to compute the stripe detection on.");
    required.add({"--outFileName", "-o"})
        .required()
        .output({"tsv"})
        .help("Outfile name to store the called stripes.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--chromosomes"})
        .nargs("+")
        .help("Chromosomes to search for stripes. If not set, all chromosomes are used.");
    optional.add({"--norm"})
        .default_value("weight")
        .help("Normalisation: 'weight', another bin table column, or 'None'.");
    optional.add({"--canny"})
        .type("float")
        .default_value(2.0)
        .help("Canny edge detection parameter (sigma of its Gaussian filter).");
    optional.add({"--minL"})
        .type("int")
        .default_value(10)
        .help("Minimum length of a stripe, in bins.");
    optional.add({"--minStripeLength"})
        .type("int")
        .help("Minimum length of a stripe, in base pairs; replaces --minL when given.");
    optional.add({"--maxWidth"})
        .type("int")
        .default_value(8)
        .help("Maximum width of a stripe, in bins.");
    optional.add({"--maxPixelPercentiles"})
        .nargs("+")
        .type("float")
        .help("Percentiles of the contact values to saturate the image. "
              "(Default: 0.95 0.96 0.97 0.98 0.99).");
    optional.add({"--pValue"})
        .type("float")
        .default_value(0.1)
        .help("P-value cutoff for a called stripe.");
    optional.add({"--mask"})
        .default_value("0")
        .help("Column coordinates to be masked in the stripiness score, e.g. chr9:12345678-12345789.");
    optional.add({"--blurFilter"})
        .type("int")
        .default_value(3)
        .help("Mean filter size, an odd number.");
    optional.add({"--seed"})
        .type("int")
        .default_value(123456789)
        .help("Seed of the pseudo random number generator that samples the background.");
    optional.add({"--sequentialRandom"})
        .action(cli::Action::StoreTrue)
        .help("Draw the background samples of all chromosomes from one generator in sequence.");
    optional.add({"--fdr"})
        .type("float")
        .help("Benjamini-Hochberg q-value threshold, replacing --pValue when given. Off by default.");
    optional.add({"--unfilteredFileName"})
        .output({"tsv"})
        .help("Also write all candidate stripes with Stripenn's result_unfiltered.tsv columns.");
    optional.add({"--dumpStages"})
        .help("Directory to write the intermediate arrays of every stage to (for the comparison with Stripenn).");
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
    args.norm = ns.str("norm");
    args.params.canny = ns.real("canny");
    args.params.min_length = static_cast<int>(ns.integer("minL"));
    if (ns.given("minStripeLength")) {
        args.min_stripe_length = ns.integer("minStripeLength");
    }
    args.params.max_width = static_cast<int>(ns.integer("maxWidth"));
    if (ns.given("maxPixelPercentiles")) {
        args.params.maxpixel = ns.reals("maxPixelPercentiles");
    }
    args.params.pvalue = ns.real("pValue");
    args.params.mask = ns.str("mask");
    args.params.bfilter = static_cast<int>(ns.integer("blurFilter"));
    args.params.seed = static_cast<std::uint64_t>(ns.integer("seed"));
    args.params.chained_rng = ns.flag("sequentialRandom");
    if (ns.given("fdr")) {
        args.fdr = ns.real("fdr");
    }
    if (ns.given("unfilteredFileName")) {
        args.unfiltered_file_name = ns.str("unfilteredFileName");
    }
    if (ns.given("dumpStages")) {
        args.dump_stages = ns.str("dumpStages");
    }
    args.threads = static_cast<int>(ns.integer("threads"));
    return args;
}

// ---------------------------------------------------------------------------
// Where the balanced pixels come from.

class PixelSource {
  public:
    virtual ~PixelSource() = default;
    [[nodiscard]] std::int64_t resolution() const { return resolution_; }
    // Every chromosome of the matrix in matrix order.
    [[nodiscard]] const std::vector<st::ChromInfo>& chromosomes() const { return chromosomes_; }
    // The balanced upper triangle pixels of a chromosome with their chromosome
    // relative bins; zero and NaN values are not passed on.
    virtual void for_each_pixel(const st::ChromInfo& chrom,
                                const std::function<void(std::int64_t, std::int64_t, double)>& visit) const = 0;
    // Bins whose balancing weight is NaN (NaN rows and columns in cooler's matrix).
    [[nodiscard]] virtual std::vector<std::uint8_t> bad_bins(const st::ChromInfo& chrom) const = 0;

  protected:
    std::int64_t resolution_ = 0;
    std::vector<st::ChromInfo> chromosomes_;
};

class CoolSource final : public PixelSource {
  public:
    CoolSource(const std::string& uri, const std::string& norm) : cool_(uri), norm_(norm) {
        const hicx::json::Value* bin_size = cool_.info_value("bin-size");
        if (bin_size == nullptr || !bin_size->is_number() || bin_size->as_double() <= 0) {
            throw std::runtime_error("the cooler has no fixed bin size");
        }
        resolution_ = static_cast<std::int64_t>(bin_size->as_double());
        const hicx::json::Value* mode = cool_.info_value("storage-mode");
        square_ = mode != nullptr && mode->is_string() && mode->as_string() == "square";
        std::int64_t first = 0;
        const auto& names = cool_.chrom_names();
        const auto& lengths = cool_.chrom_lengths();
        for (std::size_t i = 0; i < names.size(); ++i) {
            st::ChromInfo c;
            c.name = names[i];
            c.size = lengths[i];
            c.first_bin = first;
            c.nbins = st::ceil_div(lengths[i], resolution_);
            first += c.nbins;
            chromosomes_.push_back(c);
        }
        options_.apply_correction = norm != "None";
        if (options_.apply_correction) {
            if (!cool_.has_column(norm)) {
                if (norm == "weight") {
                    throw std::runtime_error("No column 'bins/weight' found in the cooler. Balance it (cooler balance) or use --norm None");
                }
                std::fprintf(stderr,
                             "WARNING:hicexplorer.hicDetectStripes:Invalid normalization method. "
                             "Normalization method is forced to None\n");
                options_.apply_correction = false;
            } else {
                options_.correction_factor_table = norm;
                // cooler's matrix(balance=<column>) multiplies by the weights.
                options_.correction_operator = '*';
            }
        }
        if (options_.apply_correction) {
            factors_ = cool_.read_column(options_.correction_factor_table);
        }
    }

    void for_each_pixel(const st::ChromInfo& chrom,
                        const std::function<void(std::int64_t, std::int64_t, double)>& visit) const override {
        const std::int64_t first = chrom.first_bin;
        const std::int64_t last = chrom.first_bin + chrom.nbins;
        hicx::for_each_corrected_pixel_chunk(cool_, first, last, first, last, options_,
                                             [&](const hicx::PixelChunk& chunk) {
            for (std::size_t k = 0; k < chunk.bin1.size(); ++k) {
                if (chunk.bin1[k] > chunk.bin2[k]) {
                    if (square_) {
                        continue;
                    }
                    visit(chunk.bin2[k] - first, chunk.bin1[k] - first, chunk.count[k]);
                    continue;
                }
                const double v = chunk.count[k];
                if (v == 0.0 || v != v) {
                    continue;
                }
                visit(chunk.bin1[k] - first, chunk.bin2[k] - first, v);
            }
        });
    }

    [[nodiscard]] std::vector<std::uint8_t> bad_bins(const st::ChromInfo& chrom) const override {
        std::vector<std::uint8_t> bad;
        if (!options_.apply_correction || factors_.empty()) {
            return bad;
        }
        bad.assign(static_cast<std::size_t>(chrom.nbins), 0);
        for (std::int64_t i = 0; i < chrom.nbins; ++i) {
            if (std::isnan(factors_[static_cast<std::size_t>(chrom.first_bin + i)])) {
                bad[static_cast<std::size_t>(i)] = 1;
            }
        }
        return bad;
    }

  private:
    hicx::CoolFile cool_;
    std::string norm_;
    hicx::CoolLoadOptions options_;
    std::vector<double> factors_;
    bool square_ = false;
};

class H5Source final : public PixelSource {
  public:
    explicit H5Source(const std::string& path) {
        whole_ = hicx::read_hicexplorer_h5(path);
        whole_.matrix.symmetrize_in_place();
        resolution_ = hicx::BinTable(whole_.cut_intervals).bin_size();
        if (resolution_ <= 0) {
            throw std::runtime_error("the matrix has no fixed bin size");
        }
        for (const auto& [name, range] : hicx::chrom_bin_boundaries(whole_.cut_intervals)) {
            st::ChromInfo c;
            c.name = name;
            c.first_bin = range.first;
            c.nbins = range.last - range.first;
            c.size = whole_.cut_intervals[static_cast<std::size_t>(range.last - 1)].end;
            chromosomes_.push_back(c);
        }
    }

    void for_each_pixel(const st::ChromInfo& chrom,
                        const std::function<void(std::int64_t, std::int64_t, double)>& visit) const override {
        std::vector<std::int64_t> selection(static_cast<std::size_t>(chrom.nbins));
        for (std::int64_t i = 0; i < chrom.nbins; ++i) {
            selection[static_cast<std::size_t>(i)] = chrom.first_bin + i;
        }
        const hicx::CsrMatrix block = hicx::select_bins(whole_.matrix, selection);
        const std::vector<std::int64_t>& indptr = block.indptr();
        const std::vector<std::int32_t>& indices = block.indices();
        const std::vector<double>& values = block.data();
        for (std::int64_t row = 0; row < block.rows(); ++row) {
            for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
                 k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
                const std::int64_t col = indices[k];
                if (col < row || values[k] == 0.0 || values[k] != values[k]) {
                    continue;
                }
                visit(row, col, values[k]);
            }
        }
    }

    [[nodiscard]] std::vector<std::uint8_t> bad_bins(const st::ChromInfo&) const override { return {}; }

  private:
    hicx::MatrixData whole_;
};

// Stripenn's list of chromosomes (compute(): the matrix's chromosomes without
// unplaced contigs, chrM and chrY).
std::vector<st::ChromInfo> stripenn_chromosomes(const std::vector<st::ChromInfo>& all) {
    std::vector<st::ChromInfo> kept;
    for (const st::ChromInfo& c : all) {
        const std::string& n = c.name;
        if (n.find("JH5") != std::string::npos || n.find("GL4") != std::string::npos ||
            n.find("RANDOM") != std::string::npos || n == "M" || n == "chrM" || n == "Y" || n == "chrY") {
            continue;
        }
        kept.push_back(c);
    }
    return kept;
}

// ---------------------------------------------------------------------------
// Dumps for the oracle comparison.

class Dumper {
  public:
    explicit Dumper(std::string dir) : dir_(std::move(dir)) {
        if (!dir_.empty()) {
            std::filesystem::create_directories(dir_);
        }
    }
    [[nodiscard]] bool enabled() const { return !dir_.empty(); }
    [[nodiscard]] const std::string& dir() const { return dir_; }

  private:
    std::string dir_;
};

void write_table(const std::string& path, const std::vector<st::Stripe>& rows,
                 const std::vector<st::ChromInfo>& requested, bool with_internal) {
    std::ofstream out(path);
    out << "chr\tpos1\tpos2\tchr2\tpos3\tpos4\tlength\twidth\tMean\tmaxpixel\tpvalue\tStripiness";
    if (with_internal) {
        out << "\ttotal\tnum\tstart\tend\tx\ty\th\tw";
    }
    out << '\n';
    for (const st::Stripe& s : rows) {
        const std::string& name = requested[static_cast<std::size_t>(s.chrom)].name;
        out << name << '\t' << s.pos1 << '\t' << s.pos2 << '\t' << name << '\t' << s.pos3 << '\t' << s.pos4 << '\t'
            << s.length << '\t' << s.width << '\t' << st::py_repr(s.mean) << '\t'
            << st::py_repr(s.maxpixel * 100.0) << "%\t"
            << (std::isnan(s.pvalue) ? std::string() : st::py_repr(s.pvalue)) << '\t'
            << (std::isnan(s.stripiness) ? std::string() : st::py_repr(s.stripiness));
        if (with_internal) {
            out << '\t' << st::py_repr(s.total) << '\t' << s.num << '\t' << s.start << '\t' << s.end << '\t' << s.x
                << '\t' << s.y << '\t' << s.h << '\t' << s.w;
        }
        out << '\n';
    }
}

void dump_array(const std::string& path, const std::vector<double>& data, const std::vector<std::size_t>& shape) {
    st::npy::save(path, data, shape);
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    const auto threads = static_cast<unsigned int>(std::max(1, args.threads));
    st::Parameters params = args.params;
    const Dumper dump(args.dump_stages);

    try {
        const bool is_cooler = hicx::check_cooler(args.matrix);
        std::unique_ptr<PixelSource> source;
        if (is_cooler) {
            source = std::make_unique<CoolSource>(args.matrix, args.norm);
        } else {
            source = std::make_unique<H5Source>(args.matrix);
        }
        const std::int64_t res = source->resolution();
        if (args.min_stripe_length.has_value()) {
            params.min_length = static_cast<int>(*args.min_stripe_length / res);
        }
        const std::vector<st::ChromInfo> all = stripenn_chromosomes(source->chromosomes());
        if (all.empty()) {
            throw std::runtime_error("no chromosomes left in the matrix");
        }
        std::vector<st::ChromInfo> requested;
        if (args.chromosomes.has_value()) {
            for (const std::string& name : *args.chromosomes) {
                const auto found = std::find_if(all.begin(), all.end(), [&](const st::ChromInfo& c) { return c.name == name; });
                if (found == all.end()) {
                    throw std::runtime_error("Chromosome name not in matrix. '" + name + "'");
                }
                requested.push_back(*found);
            }
        } else {
            requested = all;
        }
        const std::int64_t depth = st::required_depth(res, params.max_width);

        // ---- the null distribution, from all chromosomes (nulldist)
        std::vector<std::int64_t> available(all.size(), 0);
        hicx::parallel_for(all.size(), threads, [&](std::size_t i) {
            st::BandMatrix flags(all[i].nbins, -1, res, all[i].size);
            source->for_each_pixel(all[i], [&](std::int64_t a, std::int64_t b, double v) { flags.add(a, b, v); });
            available[i] = st::available_columns(flags);
        });
        const std::vector<std::int64_t> samples = st::allocate_samples(all, available);
        std::vector<std::size_t> null_chroms;
        for (std::size_t i = 0; i < all.size(); ++i) {
            if (samples[i] >= 0) {
                null_chroms.push_back(i);
            }
        }
        std::vector<st::NullColumns> parts(null_chroms.size());
        auto null_for = [&](std::size_t k, hicx::stripes::np::PyRandom& rng) {
            const st::ChromInfo& c = all[null_chroms[k]];
            st::BandMatrix band(c.nbins, depth, res, c.size);
            source->for_each_pixel(c, [&](std::int64_t a, std::int64_t b, double v) { band.add(a, b, v); });
            band.set_bad(source->bad_bins(c));
            parts[k] = st::compute_null(band, samples[null_chroms[k]], rng);
        };
        if (params.chained_rng) {
            st::np::PyRandom rng(params.seed);
            for (std::size_t k = 0; k < null_chroms.size(); ++k) {
                null_for(k, rng);
            }
        } else {
            hicx::parallel_for(null_chroms.size(), threads, [&](std::size_t k) {
                st::np::PyRandom rng(params.seed);
                null_for(k, rng);
            });
        }
        const st::NullTables null = st::join_null(parts);
        parts.clear();
        if (dump.enabled()) {
            const std::vector<std::size_t> shape{400, static_cast<std::size_t>(null.columns)};
            dump_array(dump.dir() + "/null_leftup.npy", null.leftup, shape);
            dump_array(dump.dir() + "/null_rightup.npy", null.rightup, shape);
            dump_array(dump.dir() + "/null_leftdown.npy", null.leftdown, shape);
            dump_array(dump.dir() + "/null_rightdown.npy", null.rightdown, shape);
        }

        // ---- per requested chromosome
        const std::size_t n_perc = params.maxpixel.size();
        std::vector<std::vector<st::Stripe>> extract(n_perc);  // per percentile, chromosomes in order
        std::vector<st::Stripe> final_rows;
        std::vector<std::string> maxpixel_json;
        for (std::size_t ci = 0; ci < requested.size(); ++ci) {
            const st::ChromInfo& chrom = requested[ci];
            st::BandMatrix band(chrom.nbins, depth, res, chrom.size);
            st::QuantileSelector selector(params.maxpixel);
            source->for_each_pixel(chrom, [&](std::int64_t a, std::int64_t b, double v) {
                band.add(a, b, v);
                selector.count(v, a == b);
            });
            band.set_bad(source->bad_bins(chrom));
            selector.prepare();
            source->for_each_pixel(chrom, [&](std::int64_t a, std::int64_t b, double v) { selector.gather(v, a == b); });
            const std::vector<double> mp = selector.result();
            const std::vector<double> expected = st::expected_values(band);
            if (dump.enabled()) {
                dump_array(dump.dir() + "/expected_" + chrom.name + ".npy", expected, {expected.size()});
                std::string line = "\"" + chrom.name + "\": [";
                for (std::size_t p = 0; p < mp.size(); ++p) {
                    line += (p ? ", " : "") + st::py_repr(mp[p]);
                }
                maxpixel_json.push_back(line + "]");
            }

            const std::int64_t nframes = st::ceil_div(chrom.nbins, 200);
            std::vector<st::FrameResult> frames(n_perc * static_cast<std::size_t>(nframes));
            hicx::parallel_for(frames.size(), threads, [&](std::size_t task) {
                const std::size_t p = task / static_cast<std::size_t>(nframes);
                const std::int64_t f = static_cast<std::int64_t>(task % static_cast<std::size_t>(nframes));
                st::FrameTrace trace;
                frames[task] = st::search_frame(band, chrom, f, mp[p], params, dump.enabled() ? &trace : nullptr);
                if (dump.enabled() && frames[task].searched) {
                    const std::string dir = dump.dir() + "/frames/p" + st::py_repr(params.maxpixel[p]) + "_f" + std::to_string(f);
                    std::filesystem::create_directories(dir);
                    const std::size_t S = static_cast<std::size_t>(trace.size);
                    const std::string prefix = dir + "/" + chrom.name + "_";
                    st::npy::save(prefix + "submat.npy", frames[task].submat, {S, S});
                    st::npy::save(prefix + "start_array.npy", frames[task].start_array, {S});
                    st::npy::save(prefix + "end_array.npy", frames[task].end_array, {S});
                    st::npy::save(prefix + "gray.npy", trace.gray, {6, S, S});
                    st::npy::save(prefix + "edges.npy", trace.edges, {6, S, S});
                    st::npy::save(prefix + "vert.npy", trace.vert, {6, S, S});
                    std::vector<st::Stripe> rows = frames[task].stripes;
                    for (st::Stripe& s : rows) {
                        s.chrom = static_cast<std::int32_t>(ci);
                        s.maxpixel = params.maxpixel[p];
                    }
                    write_table(prefix + "result.tsv", rows, requested, true);
                }
            });

            std::vector<st::Stripe> all_percs;
            for (std::size_t p = 0; p < n_perc; ++p) {
                std::vector<st::Stripe> merged;
                for (std::int64_t f = 0; f < nframes; ++f) {
                    for (const st::Stripe& s : frames[p * static_cast<std::size_t>(nframes) + static_cast<std::size_t>(f)].stripes) {
                        merged.push_back(s);
                    }
                }
                std::vector<st::Stripe> kept = st::remove_redundant(merged, st::RedundantBy::Size);
                for (st::Stripe& s : kept) {
                    s.chrom = static_cast<std::int32_t>(ci);
                    s.maxpixel = params.maxpixel[p];
                    s.perc_index = static_cast<int>(p);
                }
                st::compute_pvalues(band, chrom, null, kept, threads);
                for (const st::Stripe& s : kept) {
                    extract[p].push_back(s);
                    all_percs.push_back(s);
                }
            }
            std::vector<st::Stripe> chrom_final = st::remove_redundant(all_percs, st::RedundantBy::Pvalue);
            st::compute_scores(band, chrom, expected, params, chrom_final, threads);
            for (const st::Stripe& s : chrom_final) {
                final_rows.push_back(s);
            }
            std::fprintf(stderr, "INFO:hicexplorer.hicDetectStripes:%s: %zu candidates\n", chrom.name.c_str(),
                         chrom_final.size());
        }
        if (dump.enabled()) {
            std::ofstream json(dump.dir() + "/maxpixel.json");
            json << "{";
            for (std::size_t i = 0; i < maxpixel_json.size(); ++i) {
                json << (i ? ", " : "") << maxpixel_json[i];
            }
            json << "}\n";
            for (std::size_t p = 0; p < n_perc; ++p) {
                write_table(dump.dir() + "/extract_" + std::to_string(p) + ".tsv", extract[p], requested, true);
            }
        }

        // Stripenn's table order: percentile, then chromosome, then position.
        std::stable_sort(final_rows.begin(), final_rows.end(), [](const st::Stripe& a, const st::Stripe& b) {
            if (a.perc_index != b.perc_index) {
                return a.perc_index < b.perc_index;
            }
            return a.chrom < b.chrom;
        });
        if (dump.enabled()) {
            write_table(dump.dir() + "/result_unfiltered.tsv", final_rows, requested, true);
        }
        if (!args.unfiltered_file_name.empty()) {
            write_table(args.unfiltered_file_name, final_rows, requested, false);
        }

        std::vector<double> pvalues;
        pvalues.reserve(final_rows.size());
        for (const st::Stripe& s : final_rows) {
            pvalues.push_back(s.pvalue);
        }
        const std::vector<double> qvalues = hicx::stats::benjamini_hochberg_adjusted(pvalues);

        std::vector<std::size_t> called;
        for (std::size_t i = 0; i < final_rows.size(); ++i) {
            const bool pass = args.fdr.has_value() ? qvalues[i] <= *args.fdr : final_rows[i].pvalue < params.pvalue;
            if (pass) {
                called.push_back(i);
            }
        }
        // result_filtered.sort_values(by='Stripiness', ascending=False): NaN last.
        std::stable_sort(called.begin(), called.end(), [&](std::size_t a, std::size_t b) {
            const double sa = final_rows[a].stripiness;
            const double sb = final_rows[b].stripiness;
            if (std::isnan(sa) || std::isnan(sb)) {
                return !std::isnan(sa) && std::isnan(sb);
            }
            return sa > sb;
        });

        std::ofstream output(args.out_file_name);
        if (!output) {
            std::fprintf(stderr, "hicDetectStripes: error: cannot write %s\n", args.out_file_name.c_str());
            return 1;
        }
        output << "chrom\tanchor_start\tanchor_end\torientation\textent_start\textent_end\tenrichment\tpvalue\tqvalue\t"
                  "stripiness\tmean\tmaxpixel\n";
        for (const std::size_t i : called) {
            const st::Stripe& s = final_rows[i];
            const bool horizontal = s.pos1 == s.pos3;
            output << requested[static_cast<std::size_t>(s.chrom)].name << '\t' << s.pos1 - 1 << '\t' << s.pos2 << '\t'
                   << (horizontal ? "horizontal" : "vertical") << '\t' << std::min(s.pos1, s.pos3) - 1 << '\t'
                   << std::max(s.pos2, s.pos4) << '\t' << st::py_repr(s.oe_mean) << '\t' << st::py_repr(s.pvalue)
                   << '\t' << st::py_repr(qvalues[i]) << '\t' << st::py_repr(s.stripiness) << '\t'
                   << st::py_repr(s.mean) << '\t' << st::py_repr(s.maxpixel * 100.0) << "%\n";
        }
        std::fprintf(stderr, "INFO:hicexplorer.hicDetectStripes:Number of detected stripes: %zu\n", called.size());
    } catch (const std::exception& error) {
        std::fprintf(stderr, "ERROR:hicexplorer.hicDetectStripes:%s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicDetectStripes");
    return 0;
}
