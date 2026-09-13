// Port of hicexplorer/hicTransform.py.
//
// Turns a contact matrix into one of five transforms: three flavours of
// observed over expected, the Pearson correlation matrix, or the covariance
// matrix.
//
// Nothing here is new machinery. hicPCA already needed the obs/exp kernels,
// the covariance identity and the numpy-corrcoef scaling, and they live in
// core/transform_ops. What hicTransform adds is the one thing hicPCA does not
// need: the dense result has to be *written out*, and on the designated large
// input it does not fit in memory as a CsrMatrix.
//
// The memory problem, which is the point of this tool
// ---------------------------------------------------
// `hicTransform --method pearson Li_et_al_2015.h5` peaks at 5,150 MB in the
// Python against a 19.9 MB working set, a factor of 259 and the second worst
// ratio in the corpus (cpp/PLAN.md 4.3). The accounting: the branch densifies
// the whole 11,104-bin matrix (987 MB), np.corrcoef allocates its own copy and
// two intermediates of the same size, the result is copied into a lil_matrix
// and then into a CSR. Five to six live copies of a 987 MB block.
//
// Removing four of them is not enough. Two remain that are not obviously
// avoidable: the dense input block the correlation is computed from, and the
// dense output, which as a CsrMatrix is 123 M stored entries over 11,104 bins,
// 1.5 GB, against a budget of 1,221 MB (cpp/PLAN.md 4.5). This port has
// neither:
//
//  * the input block never exists, because the covariance comes straight out
//    of the CSR through cov(A) = (A A^T - n m m^T) / (n - 1), which for a
//    symmetric A is what transform_ops::DenseCorrelationRows evaluates row by
//    row;
//  * the output block never exists, because the rows are handed to the h5
//    writer as they are produced (cpp/PLAN.md 4.4 rule 7).
//
// The resident set is then the sparse matrix and O(n) vectors, whatever the
// density of the result. Measured on Li_et_al_2015.h5, whose Pearson matrix
// has 52,526,125 stored entries and is 373 MB on disk: **96.7 MB peak RSS
// against the Python's 5,534 MB**, a factor of 57, and 4.06 s of CPU against
// 66.3 s, a factor of 16. Every one of the 52.5 M values agrees at ED.
//
// The cool writer has no streaming entry point, so a .cool output does hold
// the result as a CsrMatrix. That is fine for every cool case in the corpus
// and it is recorded rather than hidden: a whole-genome pearson to cool on a
// matrix with many bins would be as large as the Python's, and the h5 path is
// the one the budget table describes.
//
// One flag the Python does not have, as in hicPCA:
//   --threads N   how many workers the dense rows are split over. The rows are
//                 taken in fixed contiguous ranges and emitted in row order,
//                 so the output is byte-identical for any N.
//
// What the profile says, on Li_et_al_2015.h5 --method pearson, 4.08 s of CPU
// in total. Reading and symmetrising the matrix is 1.5 s (measured as the whole
// of --method obs_exp, 1.82 s, minus its small write), the three row passes are
// 1.3 s (measured as the 0.44 s gap between --method covariance at two passes
// and --method pearson at three), and the remaining 1.2 s is blosc compressing
// the 373 MB result. The arithmetic is therefore a third of the run and the
// file layer is the rest, which is where the next work is if this tool ever
// needs to be faster.
//
// Two optimisations measured, one kept:
//
//  * **Flattening the inner accumulation when the block is the whole matrix**,
//    which removes two bounds tests per stored entry and lets the loop
//    auto-vectorise (cpp/OPTIMIZATION.md section 5). Interleaved A/B, eight
//    runs each on a loaded machine: 4.317 s against 4.076 s of CPU, 5.6 % of
//    the whole run and about 19 % of the accumulation, output byte-identical.
//    Kept.
//  * **More threads.** 1, 4 and 16 workers cost 4.09, 4.06 and 4.74 s of CPU
//    for 2.56, 1.05 and 0.78 s of wall time. Sixteen workers buy 26 % of wall
//    clock for 16 % more CPU on a machine that is already oversubscribed, and
//    the harness gates on CPU. The default stays 4.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "hicx/adjust_ops.hpp"
#include "hicx/bins.hpp"
#include "hicx/cool_file.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/matrix_ops.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/transform_ops.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicTransform --matrix MATRIX --outFileName OUTFILENAME\n"
    "                    [--method {obs_exp,obs_exp_lieberman,obs_exp_non_zero,"
    "pearson,covariance}]\n"
    "                    [--ligation_factor]\n"
    "                    [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]\n"
    "                    [--perChromosome] [--threads THREADS]\n"
    "                    [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Converts the (interaction) matrix to different types of obs/exp, pearson\n"
    "or covariance matrix.\n"
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        input file. The computation is done per chromosome.\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        File name to save the exported matrix.\n"
    "\n"
    "Optional arguments:\n"
    "  --method {obs_exp,obs_exp_lieberman,obs_exp_non_zero,pearson,covariance},"
    " -me ...\n"
    "                        Transformation method to use for the input matrix.\n"
    "                        (Default: obs_exp).\n"
    "  --ligation_factor     Multiply a scaling factor to each entry of the\n"
    "                        expected matrix, as the Homer software does. Only\n"
    "                        effective with obs_exp_non_zero.\n"
    "  --chromosomes CHROMOSOMES [CHROMOSOMES ...]\n"
    "                        List of chromosomes to be included in the\n"
    "                        computation.\n"
    "  --perChromosome, -pc  Each chromosome is processed individually,\n"
    "                        inter-chromosomal interactions are ignored. Option\n"
    "                        not valid for obs_exp_lieberman.\n"
    "  --threads THREADS     Workers the dense pearson and covariance rows are\n"
    "                        split over. The output is byte-identical for any\n"
    "                        value. Not a Python option; see\n"
    "                        cpp/OPTIMIZATION.md. (Default: 4).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::string matrix;
    std::string out_file_name;
    std::string method = "obs_exp";
    bool ligation_factor = false;
    std::vector<std::string> chromosomes;
    bool per_chromosome = false;
    int threads = 4;
};

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicTransform: error: %s\n", message.c_str());
    std::exit(2);
}

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

Arguments parse_arguments(int argc, char** argv) {
    Arguments args;
    bool matrix_seen = false;
    bool output_seen = false;

    std::vector<std::string>* collecting = nullptr;
    std::string* pending_string = nullptr;
    int* pending_int = nullptr;

    // argparse nargs='+' keeps consuming tokens until the next option.
    const auto looks_like_option = [](const std::string& token) {
        return token.size() > 1 && token[0] == '-' &&
               (std::isdigit(static_cast<unsigned char>(token[1])) == 0);
    };

    for (int i = 1; i < argc; ++i) {
        const std::string token(argv[i]);
        if (pending_string != nullptr) {
            *pending_string = token;
            pending_string = nullptr;
            continue;
        }
        if (pending_int != nullptr) {
            try {
                *pending_int = std::stoi(token);
            } catch (const std::exception&) {
                fail("argument --threads: invalid int value: '" + token + "'");
            }
            pending_int = nullptr;
            continue;
        }
        if (collecting != nullptr && !looks_like_option(token)) {
            collecting->push_back(token);
            continue;
        }
        collecting = nullptr;

        std::string name = token;
        std::optional<std::string> inline_value;
        const std::size_t equals = token.find('=');
        if (equals != std::string::npos && token.rfind("--", 0) == 0) {
            name = token.substr(0, equals);
            inline_value = token.substr(equals + 1);
        }

        if (name == "-h" || name == "--help") {
            std::fputs(kUsage, stdout);
            std::fputs(kHelp, stdout);
            std::exit(0);
        }
        if (name == "--version") {
            std::printf("hicTransform %s\n", hicx::kVersion);
            std::exit(0);
        }
        if (name == "--ligation_factor") {
            args.ligation_factor = true;
            continue;
        }
        if (name == "--perChromosome" || name == "-pc") {
            args.per_chromosome = true;
            continue;
        }
        if (name == "-m" || name == "--matrix") {
            matrix_seen = true;
            if (inline_value) {
                args.matrix = *inline_value;
            } else {
                pending_string = &args.matrix;
            }
            continue;
        }
        if (name == "-o" || name == "--outFileName") {
            output_seen = true;
            if (inline_value) {
                args.out_file_name = *inline_value;
            } else {
                pending_string = &args.out_file_name;
            }
            continue;
        }
        if (name == "-me" || name == "--method") {
            if (inline_value) {
                args.method = *inline_value;
            } else {
                pending_string = &args.method;
            }
            continue;
        }
        if (name == "--chromosomes") {
            if (inline_value) {
                args.chromosomes.push_back(*inline_value);
            } else {
                collecting = &args.chromosomes;
            }
            continue;
        }
        if (name == "--threads" || name == "-t") {
            if (inline_value) {
                try {
                    args.threads = std::stoi(*inline_value);
                } catch (const std::exception&) {
                    fail("argument --threads: invalid int value: '" + *inline_value + "'");
                }
            } else {
                pending_int = &args.threads;
            }
            continue;
        }
        fail("unrecognized arguments: " + token);
    }
    if (pending_string != nullptr || pending_int != nullptr) {
        fail("expected one argument");
    }

    std::string missing;
    if (!matrix_seen) {
        missing = "--matrix/-m";
    }
    if (!output_seen) {
        missing += missing.empty() ? "--outFileName/-o" : ", --outFileName/-o";
    }
    if (!missing.empty()) {
        fail("the following arguments are required: " + missing);
    }
    static const char* const kMethods[] = {"obs_exp", "obs_exp_lieberman",
                                           "obs_exp_non_zero", "pearson", "covariance"};
    if (std::find(std::begin(kMethods), std::end(kMethods), args.method) ==
        std::end(kMethods)) {
        fail("argument --method/-me: invalid choice: '" + args.method +
             "' (choose from 'obs_exp', 'obs_exp_lieberman', 'obs_exp_non_zero', "
             "'pearson', 'covariance')");
    }
    if (args.threads < 1) {
        fail("argument --threads: must be at least 1");
    }
    return args;
}

// A CSR matrix assembled one complete row at a time, in increasing row order.
// This is what replaces hicTransform.py:158's lil_matrix accumulator: no per
// element Python object, no second copy at .tocsr() time (cpp/PLAN.md 4.4
// rules 4 and 7).
class RowwiseCsrBuilder {
  public:
    explicit RowwiseCsrBuilder(std::int64_t n) : n_(n), indptr_(1, 0) {}

    void skip_to(std::int64_t row) {
        while (next_row_ < row) {
            indptr_.push_back(static_cast<std::int64_t>(data_.size()));
            ++next_row_;
        }
    }
    void push(std::int64_t column, double value) {
        indices_.push_back(static_cast<std::int32_t>(column));
        data_.push_back(value);
    }
    void end_row() {
        indptr_.push_back(static_cast<std::int64_t>(data_.size()));
        ++next_row_;
    }
    [[nodiscard]] hicx::CsrMatrix finish(std::string dtype) {
        skip_to(n_);
        return hicx::CsrMatrix(n_, n_, std::move(indptr_), std::move(indices_),
                               std::move(data_), std::move(dtype));
    }

  private:
    std::int64_t n_;
    std::int64_t next_row_ = 0;
    std::vector<std::int64_t> indptr_;
    std::vector<std::int32_t> indices_;
    std::vector<double> data_;
};

// hicTransform's four private helpers, which all begin with
// `if len(pSubmatrix.data) == 0: return pSubmatrix` and end with
// convertNansToZeros followed by convertInfsToZeros.
void apply_obs_exp(hicx::CsrMatrix& block, const std::string& method,
                   bool ligation_factor, std::int64_t length_chromosome,
                   std::int64_t chromosome_count) {
    if (block.stored_nnz() == 0) {
        return;
    }
    if (method == "obs_exp") {
        hicx::obs_exp_in_place(block);
    } else if (method == "obs_exp_non_zero") {
        hicx::obs_exp_non_zero_in_place(block, ligation_factor);
    } else {
        hicx::obs_exp_lieberman_in_place(block, length_chromosome, chromosome_count);
    }
    hicx::convert_nans_and_infs_to_zeros(block);
}

// transform_ops::DenseCorrelationRows behind the writer's interface.
class CorrelationRowSource : public hicx::DenseRowSource {
  public:
    CorrelationRowSource(const hicx::CsrMatrix& matrix, std::vector<hicx::BinRange> blocks,
                         hicx::DenseCorrelationRows::Kind kind, int threads)
        : rows_(matrix, std::move(blocks), kind, threads), dtype_("float64") {}

    [[nodiscard]] std::int64_t rows() const override { return rows_.size(); }
    [[nodiscard]] const std::string& dtype() const override { return dtype_; }
    void fill_row(std::int64_t i, double* out) const override { rows_.fill_row(i, out); }

  private:
    hicx::DenseCorrelationRows rows_;
    std::string dtype_;
};

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);

    // hicTransform.py:145-149. The check is on the *name* only; which writer
    // actually runs is decided by the input format, see below.
    if (!ends_with(args.out_file_name, ".h5") && !ends_with(args.out_file_name, ".cool")) {
        std::fprintf(stderr, "ERROR:hicexplorer.hicTransform:Output filetype not known.\n");
        std::fprintf(stderr, "ERROR:hicexplorer.hicTransform:It is: %s\n",
                     args.out_file_name.c_str());
        std::fprintf(stderr, "ERROR:hicexplorer.hicTransform:Accepted is .h5 or .cool\n");
        return 1;
    }

    try {
        // hicTransform.py:151-156. A cool input with exactly one chromosome is
        // loaded through the cool reader's pChrnameList; anything else is
        // loaded whole and then reduced with keepOnlyTheseChr. The two are
        // different code paths in hicmatrix and are kept different here.
        std::optional<std::string> load_chromosome;
        const bool cool_input = ends_with(args.matrix, "cool");
        if (cool_input && args.chromosomes.size() == 1) {
            load_chromosome = args.chromosomes.front();
        }
        hicx::ToolMatrix hic = hicx::ToolMatrix::load(args.matrix, load_chromosome);
        if (load_chromosome.has_value()) {
            // The cool reader returns an empty matrix for a chromosome the
            // file does not have, and the run then fails much later with an
            // unrelated message from the writer. cooler raises
            // "Unknown sequence label" at load time, so the check belongs
            // here. Reported as a defect of the shared cool reader, which this
            // workstream does not own; the exit status was already 1 either
            // way, only the message was wrong.
            bool found = false;
            for (const std::pair<std::string, hicx::BinRange>& entry : hic.boundaries()) {
                if (entry.first == *load_chromosome) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw std::runtime_error("Unknown sequence label: " + *load_chromosome);
            }
        }
        if (!load_chromosome.has_value() && !args.chromosomes.empty()) {
            hicx::keep_only_chromosomes(hic.data(), args.chromosomes);
            hic.refresh_boundaries();
        }

        // Every kernel below reads whole rows, so the matrix needs its
        // explicit symmetric form: this is one of the tools cpp/PLAN.md 4.4
        // rule 2 exempts from upper-triangle-only storage. It is the same
        // matrix the Python holds after fillLowerTriangle, and it is still
        // half of what the Python's smallest live set is, because the Python
        // then densifies it.
        hic.matrix().materialize_full();

        const std::vector<std::pair<std::string, hicx::BinRange>>& boundaries =
            hic.boundaries();
        const std::int64_t nbins = hic.matrix().rows();

        // The row-wise builders and the streaming writer both require the
        // chromosome blocks to arrive in increasing bin order, which is the
        // order a bin table taken from a file has.
        std::int64_t previous_end = 0;
        for (const std::pair<std::string, hicx::BinRange>& entry : boundaries) {
            if (entry.second.first < previous_end) {
                throw std::runtime_error(
                    "the chromosomes of this matrix are not in bin order, which the "
                    "transform writers require");
            }
            previous_end = entry.second.last;
        }

        const bool is_correlation = args.method == "pearson" || args.method == "covariance";
        // obs_exp_lieberman ignores --perChromosome: hicTransform.py:194-211
        // always runs the per chromosome loop for it.
        const bool per_chromosome =
            args.per_chromosome || args.method == "obs_exp_lieberman";

        if (is_correlation) {
            std::vector<hicx::BinRange> blocks;
            if (per_chromosome) {
                for (const std::pair<std::string, hicx::BinRange>& entry : boundaries) {
                    blocks.push_back(entry.second);
                }
            } else {
                blocks.push_back(hicx::BinRange{0, nbins});
            }
            const hicx::DenseCorrelationRows::Kind kind =
                args.method == "pearson" ? hicx::DenseCorrelationRows::Kind::Pearson
                                         : hicx::DenseCorrelationRows::Kind::Covariance;
            const CorrelationRowSource source(hic.matrix(), std::move(blocks), kind,
                                              args.threads);

            if (hic.input_is_h5()) {
                // The streaming path: nothing dense is ever resident.
                hicx::H5SaveOptions options;
                options.symmetric = true;
                hicx::write_hicexplorer_h5(args.out_file_name, hic.data(), source,
                                           args.threads, options);
            } else {
                // cooler has no streaming writer here, so the result is
                // materialised. Recorded in the header comment as the one
                // path that still holds a dense result.
                hicx::CsrMatrix result =
                    hicx::materialize_dense_row_source(source, args.threads);
                hic.matrix() = std::move(result);
                hic.save(args.out_file_name);
            }
        } else if (!per_chromosome) {
            // The whole-matrix branch. The transform runs on the genome-wide
            // matrix in place, and the dtype the kernel leaves behind is the
            // dtype that is written: obs_exp casts back to the input dtype,
            // obs_exp_non_zero does not and stays float32.
            apply_obs_exp(hic.matrix(), args.method, args.ligation_factor, 0, 0);
            hic.save(args.out_file_name);
        } else {
            // The per chromosome branch. Two things the accumulator decides
            // and the kernels do not: only the diagonal blocks survive, and
            // the result is float64 whatever the blocks were, because
            // lil_matrix(shape) is float64 and assigning an integer block into
            // it does not change that (hicTransform.py:158).
            const std::int64_t chromosome_count = static_cast<std::int64_t>(boundaries.size());
            std::int64_t length_chromosome = 0;
            for (const std::pair<std::string, hicx::BinRange>& entry : boundaries) {
                length_chromosome += entry.second.last - entry.second.first;
            }

            RowwiseCsrBuilder builder(nbins);
            for (const std::pair<std::string, hicx::BinRange>& entry : boundaries) {
                const std::int64_t first = entry.second.first;
                const std::int64_t last = entry.second.last;
                std::vector<std::int64_t> order(static_cast<std::size_t>(last - first));
                for (std::int64_t i = first; i < last; ++i) {
                    order[static_cast<std::size_t>(i - first)] = i;
                }
                hicx::CsrMatrix block = hicx::select_bins(hic.matrix(), order);
                block.materialize_full();
                apply_obs_exp(block, args.method, args.ligation_factor, length_chromosome,
                              chromosome_count);

                builder.skip_to(first);
                for (std::int64_t i = 0; i < last - first; ++i) {
                    const std::size_t begin =
                        static_cast<std::size_t>(block.indptr()[static_cast<std::size_t>(i)]);
                    const std::size_t end = static_cast<std::size_t>(
                        block.indptr()[static_cast<std::size_t>(i) + 1]);
                    for (std::size_t k = begin; k < end; ++k) {
                        const double value = block.data()[k];
                        // Assigning a block into a lil_matrix drops every
                        // explicitly stored zero, unlike lil_matrix(block).
                        if (value != 0.0) {
                            builder.push(first + block.indices()[k], value);
                        }
                    }
                    builder.end_row();
                }
            }
            hic.matrix() = builder.finish("float64");
            hic.save(args.out_file_name);
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicTransform: %s\n", error.what());
        return 1;
    }

    hicx::report_resource_usage("hicTransform");
    return 0;
}
