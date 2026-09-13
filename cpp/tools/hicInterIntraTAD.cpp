// Port of hicexplorer/hicInterIntraTAD.py.
//
// For every TAD of a hicFindTADs domain file, writes the contact sum, the
// number of cells, the number of stored cells and the density of the TAD
// itself and of its left and right inter-TAD blocks, and the ratios of the
// inter sums to the intra sum.
//
// The geometry, and every quirk of it, is shared with hicDifferentialTAD and
// documented in tad_contacts_impl.hpp. What is specific to this tool is the
// text, which is Python's str() of whatever object the Python happens to hold:
//
//  * A column that was not computed because the neighbour does not exist is
//    the literal int 0 (hicInterIntraTAD.py:179-186) and prints as "0", while
//    a computed density is a Python float and prints as "1.0". The same
//    column therefore switches format between rows.
//  * A sum is the numpy scalar scipy's `.sum()` returns, so its type follows
//    the matrix dtype: int64 for an integer matrix, float32 or float64
//    otherwise, and a ratio follows numpy's promotion of the two operands.
//    tad_contacts_impl's PyNumber carries that.
//  * A block with no cells makes the density a ZeroDivisionError. The worker
//    reports it and the tool exits with status 1 before writing anything. On
//    an h5 matrix this happens for the last TAD of every chromosome with at
//    least three TADs (point 3 of tad_contacts_impl.hpp), so the Python tool
//    cannot process an h5 matrix at all; pinned by
//    test_hicInterIntraTAD.py::test_h5_input_raises_zero_division and
//    reproduced here.
//
// Not ported: the scatter plot. It is a matplotlib figure and belongs to the
// tier 7 plotting shell (cpp/PLAN.md, tier 7, option (a)). The Python writes
// it on every run, to ratio.png when --outFileNameRatioPlot is not given
// (hicInterIntraTAD.py:39-42, :513). The interim policy until plotting is
// decided, set by the project owner on 2026-09-13: a tool never exits 0
// without writing every file the user explicitly asked for, so
//
//  * --outFileNameRatioPlot / -op given explicitly: the tool exits 1 before it
//    reads or writes anything;
//  * -op not given: the default ratio.png is a side effect nobody asked for,
//    so it is skipped with a note on stderr and the table is written;
//  * --fontsize and --dpi request no file and are accepted and ignored.
//
// Threading: one independent problem per TAD, hicx::parallel_for into
// preallocated slots, written in file order, so the output does not depend on
// --threads. The Python's does in the degenerate cases listed in
// tad_contacts_impl.hpp; this port reproduces --threads 1.

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <string>
#include <vector>

#include "hicx/cool_adapter.hpp"
#include "hicx/parallel.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/version.hpp"
#include "tad_contacts_impl.hpp"

namespace {

const char* const kUsage =
    "usage: hicInterIntraTAD [--matrix MATRIX] [--tadDomains TADDOMAINS]\n"
    "                        [--outFileName OUTFILENAME]\n"
    "                        [--outFileNameRatioPlot OUTFILENAMERATIOPLOT]\n"
    "                        [--fontsize FONTSIZE] [--dpi DPI] [--threads THREADS]\n"
    "                        [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Extracts and computes different inter and intra TAD values and ratios.\n"
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        The matrix which was used to compute the TADs\n"
    "  --tadDomains TADDOMAINS, -td TADDOMAINS\n"
    "                        The TADs domain file computed by hicFindTADs.\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        Outfile name\n"
    "\n"
    "Optional arguments:\n"
    "  --outFileNameRatioPlot OUTFILENAMERATIOPLOT, -op OUTFILENAMERATIOPLOT\n"
    "                        Outfile name for the inter-left/intra vs inter-\n"
    "                        right/intra ratio plot. The plot is not yet\n"
    "                        available in the C++ port: giving this option makes\n"
    "                        the tool exit with status 1 before writing anything.\n"
    "                        Without it no plot is written.\n"
    "  --fontsize FONTSIZE   Fontsize in the plot for x and y axis. Accepted and\n"
    "                        ignored by the C++ port.\n"
    "  --dpi DPI             The dpi of the scatter plot. Accepted and ignored by\n"
    "                        the C++ port.\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads to use, the parallelization is\n"
    "                        implemented per chromosome (Default: 4).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::optional<std::string> matrix;
    std::optional<std::string> domains;
    std::string out_file = "output_interintra_tad.tzt";
    std::string plot_file = "ratio.png";
    // Whether --outFileNameRatioPlot / -op appeared on the command line, as
    // opposed to plot_file holding the argparse default.
    bool plot_requested = false;
    long long threads = 4;
};

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicInterIntraTAD: error: %s\n", message.c_str());
    std::exit(2);
}

std::string trim(const std::string& text) {
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

bool python_float(const std::string& text) {
    const std::string trimmed = trim(text);
    char* end = nullptr;
    (void)std::strtod(trimmed.c_str(), &end);
    return !trimmed.empty() && end == trimmed.c_str() + trimmed.size();
}

bool python_int(const std::string& text, long long* value) {
    const std::string trimmed = trim(text);
    char* end = nullptr;
    *value = std::strtoll(trimmed.c_str(), &end, 10);
    return !trimmed.empty() && end == trimmed.c_str() + trimmed.size();
}

bool looks_like_option(const std::string& token) {
    return token.size() > 1 && token[0] == '-' &&
           std::isdigit(static_cast<unsigned char>(token[1])) == 0 && token[1] != '.';
}

Arguments parse_arguments(int argc, char** argv) {
    Arguments args;
    const std::vector<std::string> tokens(argv + 1, argv + argc);
    for (std::size_t i = 0; i < tokens.size(); ++i) {
        const std::string& token = tokens[i];
        std::string name = token;
        std::optional<std::string> inline_value;
        if (token.rfind("--", 0) == 0) {
            const std::size_t equals = token.find('=');
            if (equals != std::string::npos) {
                name = token.substr(0, equals);
                inline_value = token.substr(equals + 1);
            }
        }
        if (name == "-h" || name == "--help") {
            std::fputs(kUsage, stdout);
            std::fputs(kHelp, stdout);
            std::exit(0);
        }
        if (name == "--version") {
            std::printf("hicInterIntraTAD %s\n", hicx::kVersion);
            std::exit(0);
        }
        const auto value = [&](const char* label) {
            if (inline_value.has_value()) {
                return *inline_value;
            }
            if (i + 1 >= tokens.size() || looks_like_option(tokens[i + 1])) {
                fail(std::string("argument ") + label + ": expected one argument");
            }
            return tokens[++i];
        };
        if (name == "--matrix" || name == "-m") {
            args.matrix = value("--matrix/-m");
        } else if (name == "--tadDomains" || name == "-td") {
            args.domains = value("--tadDomains/-td");
        } else if (name == "--outFileName" || name == "-o") {
            args.out_file = value("--outFileName/-o");
        } else if (name == "--outFileNameRatioPlot" || name == "-op") {
            args.plot_file = value("--outFileNameRatioPlot/-op");
            args.plot_requested = true;
        } else if (name == "--fontsize") {
            const std::string text = value("--fontsize");
            if (!python_float(text)) {
                fail("argument --fontsize: invalid float value: '" + text + "'");
            }
        } else if (name == "--dpi") {
            const std::string text = value("--dpi");
            long long ignored = 0;
            if (!python_int(text, &ignored)) {
                fail("argument --dpi: invalid int value: '" + text + "'");
            }
        } else if (name == "--threads" || name == "-t") {
            const std::string text = value("--threads/-t");
            if (!python_int(text, &args.threads)) {
                fail("argument --threads/-t: invalid int value: '" + text + "'");
            }
        } else {
            fail("unrecognized arguments: " + token);
        }
    }
    return args;
}

struct Side {
    hicx::tads::PyNumber sum = hicx::tads::PyNumber::py_int(0);
    hicx::tads::PyNumber density = hicx::tads::PyNumber::py_int(0);
    std::int64_t contacts = 0;
    std::int64_t nnz = 0;
};

struct TadResult {
    Side left;
    Side right;
    Side intra;
    hicx::tads::PyNumber left_ratio;
    hicx::tads::PyNumber right_ratio;
    hicx::tads::PyNumber both_ratio;
    std::optional<std::string> error;
};

// sum, shape[0] * shape[1], nnz and nnz / (shape[0] * shape[1]) of one block.
Side measure(const hicx::tads::ContactMatrix& matrix, const hicx::tads::Block& block) {
    namespace tads = hicx::tads;
    const tads::DenseBlock values = tads::extract_block(matrix, block);
    Side side;
    side.sum = tads::block_sum(values, tads::block_dtype(matrix, block));
    side.contacts = values.rows * values.cols;
    side.nnz = values.nnz;
    if (side.contacts == 0) {
        throw std::runtime_error("division by zero");
    }
    side.density = tads::PyNumber::float64(static_cast<double>(side.nnz) /
                                           static_cast<double>(side.contacts));
    return side;
}

bool write_file(const std::string& path, const std::string& content) {
    std::FILE* handle = std::fopen(path.c_str(), "wb");
    if (handle == nullptr) {
        return false;
    }
    const bool written = std::fwrite(content.data(), 1, content.size(), handle) ==
                         content.size();
    return std::fclose(handle) == 0 && written;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    namespace tads = hicx::tads;

    if (args.plot_requested) {
        // Before anything is read or written: the user named a file this port
        // cannot produce, so it must not report success, and it must not leave
        // a table behind that looks like a completed run.
        std::fprintf(stderr,
                     "hicInterIntraTAD: error: --outFileNameRatioPlot '%s' was requested, "
                     "but the ratio plot is not yet available in the C++ port (it is a "
                     "matplotlib figure, see cpp/PLAN.md tier 7). Nothing was written. "
                     "Run without --outFileNameRatioPlot to get the table, or use the "
                     "Python hicInterIntraTAD for the plot.\n",
                     args.plot_file.c_str());
        return 1;
    }

    try {
        if (!args.domains.has_value() || !args.matrix.has_value()) {
            throw std::runtime_error("--matrix and --tadDomains must both be given");
        }
        const std::vector<tads::Domain> domains = tads::read_domains(*args.domains);
        const std::vector<std::vector<tads::Domain>> chromosomes =
            tads::group_by_chromosome(domains);
        const bool is_cooler = hicx::check_cooler(*args.matrix);
        if (args.threads == 0) {
            throw std::runtime_error("integer division or modulo by zero");
        }
        const tads::ContactMatrix matrix = tads::load_contact_matrix(*args.matrix, is_cooler);

        struct Item {
            std::size_t chromosome = 0;
            std::size_t index = 0;
        };
        std::vector<Item> items;
        if (args.threads > 0) {
            for (std::size_t c = 0; c < chromosomes.size(); ++c) {
                for (std::size_t i = 0; i < chromosomes[c].size(); ++i) {
                    items.push_back(Item{c, i});
                }
            }
        }

        std::vector<TadResult> results(items.size());
        const unsigned int workers = static_cast<unsigned int>(
            std::clamp<long long>(args.threads, 1, hicx::hardware_threads()));
        hicx::parallel_for(items.size(), workers, [&](std::size_t k) {
            TadResult& result = results[k];
            try {
                const tads::TadGeometry geometry =
                    tads::tad_geometry(matrix, chromosomes[items[k].chromosome], items[k].index);
                result.intra = measure(matrix, geometry.intra);
                if (geometry.left.has_value()) {
                    result.left = measure(matrix, *geometry.left);
                }
                if (geometry.right.has_value()) {
                    result.right = measure(matrix, *geometry.right);
                }
                result.left_ratio = tads::py_divide(result.left.sum, result.intra.sum);
                result.right_ratio = tads::py_divide(result.right.sum, result.intra.sum);
                result.both_ratio = tads::py_divide(
                    tads::py_add(result.left.sum, result.right.sum), result.intra.sum);
            } catch (const std::exception& error) {
                result.error = error.what();
            }
        });
        for (const TadResult& result : results) {
            if (result.error.has_value()) {
                std::fprintf(stderr, "ERROR:hicexplorer.hicInterIntraTAD:%s\n",
                             result.error->c_str());
                return 1;
            }
        }

        std::string text = "# Created with HiCExplorer's hicInterIntraTAD version ";
        text += hicx::kVersion;
        text += "\n";
        text += "# Chromosome\tstart\tend\tname\tscore\tstrand\tinter_left_sum\t"
                "inter_right_sum\tinter_left_density\tinter_right_density\t"
                "inter_left_number_of_contacts\tinter_right_number_of_contacts\t"
                "inter_left_number_of_contacts_nnz\tinter_right_number_of_contacts_nnz\t"
                "intra_sum\tintra_number_of_contacts\tintra_number_of_contacts_nnz\t"
                "intra_density\tinter_left_intra_ratio\tinter_right_intra_ratio\t"
                "inter_left_inter_right_intra_ratio\n";
        for (std::size_t k = 0; k < items.size(); ++k) {
            const tads::Domain& domain = chromosomes[items[k].chromosome][items[k].index];
            const TadResult& r = results[k];
            std::string line;
            for (std::size_t column = 0; column < domain.text.size(); ++column) {
                if (column > 0) {
                    line += '\t';
                }
                line += domain.text[column];
            }
            const auto field = [&line](const std::string& value) {
                line += '\t';
                line += value;
            };
            field(r.left.sum.str());
            field(r.right.sum.str());
            field(r.left.density.str());
            field(r.right.density.str());
            field(std::to_string(r.left.contacts));
            field(std::to_string(r.right.contacts));
            field(std::to_string(r.left.nnz));
            field(std::to_string(r.right.nnz));
            field(r.intra.sum.str());
            field(std::to_string(r.intra.contacts));
            field(std::to_string(r.intra.nnz));
            field(r.intra.density.str());
            field(r.left_ratio.str());
            field(r.right_ratio.str());
            field(r.both_ratio.str());
            line += '\n';
            text += line;
        }
        if (!write_file(args.out_file, text)) {
            throw std::runtime_error("cannot write '" + args.out_file + "'");
        }
        std::fprintf(stderr,
                     "hicInterIntraTAD: no ratio plot was written; the plot is a matplotlib "
                     "figure and is not yet available in the C++ port. The table in '%s' "
                     "holds its data.\n",
                     args.out_file.c_str());
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicInterIntraTAD: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicInterIntraTAD");
    return 0;
}
