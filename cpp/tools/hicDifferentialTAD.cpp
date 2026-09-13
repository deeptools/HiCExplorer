// Port of hicexplorer/hicDifferentialTAD.py.
//
// For every TAD of a hicFindTADs domain file, tests whether the contacts of a
// target and a control matrix differ, with three Wilcoxon rank sum tests: the
// TAD itself (intra-TAD), the TAD against its left neighbour and against its
// right neighbour (inter-TAD). Regions whose selected tests reach --pValue are
// written to <prefix>_rejected.diff_tad, the rest to <prefix>_accepted.diff_tad.
//
// The only scipy.stats function the Python calls is ranksums
// (hicDifferentialTAD.py:3), reused from hicx::stats, with the NaN propagation
// of scipy's _axis_nan_policy decorator added in tad_contacts_impl. Each test
// ranks the *dense* blocks, zeros included, which is what `.toarray().flatten()`
// hands to ranksums, so most ranks are one large tie group of zeros.
//
// Which blocks are tested, and the quirks of how they are cut, is shared with
// hicInterIntraTAD and documented once in tad_contacts_impl.hpp. The two
// points that are specific to this tool:
//
//  * A test that did not run (no neighbour) or whose p-value is NaN prints
//    'nan' for both its p-value and its statistic and never rejects, whatever
//    --modeReject says (hicDifferentialTAD.py:228-262). With `-mr all` a TAD
//    with a missing neighbour can therefore never be rejected, which is why
//    the first TAD of every chromosome sits in the accepted file in mode all.
//  * `p <= --pValue` rejects, so a p-value exactly on the threshold rejects.
//    No p-value of the corpus runs sits exactly on a threshold. The closest to
//    0.05 is 0.04934747494654432, 1.3 % below it. With --pValue 0.01 one value
//    lies within 1 % of the threshold: 0.010074561112534944, the left
//    inter-TAD test of chr1:148600000, 0.75 % above it, so that test does not
//    reject. Measured over every output in test_data/hicDifferentialTAD on
//    2026-09-13; that is where two implementations could legitimately
//    disagree, and on this corpus they do not.
//
// Threading: each TAD is an independent problem. The TADs of all chromosomes
// are processed by hicx::parallel_for into preallocated result slots and
// written in file order afterwards, so the output is byte-identical at any
// --threads. The Python's own output depends on --threads in three degenerate
// situations; this port always reproduces --threads 1 (tad_contacts_impl.hpp).

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "hicx/cool_file.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/parallel.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/version.hpp"
#include "tad_contacts_impl.hpp"

namespace {

const char* const kUsage =
    "usage: hicDifferentialTAD [--targetMatrix TARGETMATRIX]\n"
    "                          [--controlMatrix CONTROLMATRIX]\n"
    "                          [--tadDomains TADDOMAINS]\n"
    "                          [--outFileNamePrefix OUTFILENAMEPREFIX]\n"
    "                          [--pValue PVALUE]\n"
    "                          [--mode {intra-TAD,left-inter-TAD,right-inter-TAD,all}]\n"
    "                          [--modeReject {all,one}] [--threads THREADS] [--help]\n"
    "                          [--version]\n";

const char* const kHelp =
    "\n"
    "Computes differential TADs by comparing the precomputed TAD regions of the target "
    "matrix with the same regions of the control matrix.\n"
    "Please notice that the matrices need to have the same read coverage, this can be "
    "achieved with hicNormalize and the 'smallest'-mode.\n"
    "H0 is the assumption that two regions are identical, the rejected files contain "
    "therefore the as differential considered regions.\n"
    "\n"
    "Required arguments:\n"
    "  --targetMatrix TARGETMATRIX, -tm TARGETMATRIX\n"
    "                        The matrix which was used to compute the TADs\n"
    "  --controlMatrix CONTROLMATRIX, -cm CONTROLMATRIX\n"
    "                        The control matrix to test the TADs for a differential\n"
    "                        interaction pattern.\n"
    "  --tadDomains TADDOMAINS, -td TADDOMAINS\n"
    "                        The TADs domain file computed by hicFindTADs.\n"
    "  --outFileNamePrefix OUTFILENAMEPREFIX, -o OUTFILENAMEPREFIX\n"
    "                        Outfile name prefix to store the accepted / rejected H0\n"
    "                        TADs.\n"
    "\n"
    "Optional arguments:\n"
    "  --pValue PVALUE, -p PVALUE\n"
    "                        H0 is considered as 'two regions are identical.' i.e.\n"
    "                        all regions with a test result of <= p-value are\n"
    "                        rejected and considered as differential (Default:\n"
    "                        0.05).\n"
    "  --mode {intra-TAD,left-inter-TAD,right-inter-TAD,all}, -m "
    "{intra-TAD,left-inter-TAD,right-inter-TAD,all}\n"
    "                        Consider only intra-TAD interactions, or additional\n"
    "                        left inter-TAD, right inter-TAD or all (Default: all).\n"
    "  --modeReject {all,one}, -mr {all,one}\n"
    "                        All test of a mode must be rejected (all) or reject\n"
    "                        region (and accept it is differential) as soon as at\n"
    "                        least one region is having a p-value <= --pValue\n"
    "                        (Default: one).\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads to use, the parallelization is\n"
    "                        implemented per chromosome (Default: 4).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::optional<std::string> target;
    std::optional<std::string> control;
    std::optional<std::string> domains;
    std::string prefix = "output_differential_tad";
    double p_value = 0.05;
    std::string mode = "all";
    std::string mode_reject = "one";
    long long threads = 4;
};

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicDifferentialTAD: error: %s\n", message.c_str());
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

// float(text) and int(text), which is what argparse's type= calls.
bool python_float(const std::string& text, double* value) {
    const std::string trimmed = trim(text);
    if (trimmed.empty()) {
        return false;
    }
    char* end = nullptr;
    *value = std::strtod(trimmed.c_str(), &end);
    return end == trimmed.c_str() + trimmed.size();
}

bool python_int(const std::string& text, long long* value) {
    const std::string trimmed = trim(text);
    if (trimmed.empty()) {
        return false;
    }
    char* end = nullptr;
    *value = std::strtoll(trimmed.c_str(), &end, 10);
    return end == trimmed.c_str() + trimmed.size();
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
            std::printf("hicDifferentialTAD %s\n", hicx::kVersion);
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
        if (name == "--targetMatrix" || name == "-tm") {
            args.target = value("--targetMatrix/-tm");
        } else if (name == "--controlMatrix" || name == "-cm") {
            args.control = value("--controlMatrix/-cm");
        } else if (name == "--tadDomains" || name == "-td") {
            args.domains = value("--tadDomains/-td");
        } else if (name == "--outFileNamePrefix" || name == "-o") {
            args.prefix = value("--outFileNamePrefix/-o");
        } else if (name == "--pValue" || name == "-p") {
            const std::string text = value("--pValue/-p");
            if (!python_float(text, &args.p_value)) {
                fail("argument --pValue/-p: invalid float value: '" + text + "'");
            }
        } else if (name == "--mode" || name == "-m") {
            args.mode = value("--mode/-m");
            if (args.mode != "intra-TAD" && args.mode != "left-inter-TAD" &&
                args.mode != "right-inter-TAD" && args.mode != "all") {
                fail("argument --mode/-m: invalid choice: '" + args.mode +
                     "' (choose from 'intra-TAD', 'left-inter-TAD', 'right-inter-TAD', "
                     "'all')");
            }
        } else if (name == "--modeReject" || name == "-mr") {
            args.mode_reject = value("--modeReject/-mr");
            if (args.mode_reject != "all" && args.mode_reject != "one") {
                fail("argument --modeReject/-mr: invalid choice: '" + args.mode_reject +
                     "' (choose from 'all', 'one')");
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

struct TadResult {
    // left inter-TAD, right inter-TAD, intra-TAD, the order of the columns.
    std::array<double, 3> pvalue{};
    std::array<double, 3> statistic{};
    std::array<bool, 3> rejected{};
    std::optional<std::string> error;
};

std::string header(bool accepted, const Arguments& args) {
    std::string text = "# Created with HiCExplorer's hicDifferentialTAD version ";
    text += hicx::kVersion;
    text += "\n";
    if (accepted) {
        text += "# H0 'regions are equal' H0 is accepted for all p-value greater the user "
                "given p-value threshold; i.e. regions in this file are not considered as "
                "differential.\n";
        text += "# Accepted regions with Wilcoxon rank-sum test to p-value: ";
    } else {
        text += "# H0 'regions are equal' H0 is rejected for all p-value smaller or equal "
                "the user given p-value threshold; i.e. regions in this file are "
                "considered as differential.\n";
        text += "# Rejected regions with Wilcoxon rank-sum test to p-value: ";
    }
    text += hicx::npy::float_repr(args.p_value) + "  with used mode: " + args.mode +
            " and modeReject: " + args.mode_reject + " \n";
    text += "# Chromosome\tstart\tend\tname\tscore\tstrand\tp-value left-inter-TAD\t"
            "p-value right-inter-TAD\tp-value intra-TAD\tW left-inter-TAD\t"
            "W right-inter-TAD\tW intra-TAD\n";
    return text;
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

    try {
        if (!args.domains.has_value() || !args.target.has_value() ||
            !args.control.has_value()) {
            // All three are declared required=False; the Python then fails on
            // the first use of the missing None.
            throw std::runtime_error(
                "--targetMatrix, --controlMatrix and --tadDomains must all be given");
        }
        const std::vector<tads::Domain> domains = tads::read_domains(*args.domains);
        const std::vector<std::vector<tads::Domain>> chromosomes =
            tads::group_by_chromosome(domains);

        const bool target_is_cooler = hicx::check_cooler(*args.target);
        const bool control_is_cooler = hicx::check_cooler(*args.control);
        if (target_is_cooler != control_is_cooler) {
            std::fputs("ERROR:hicexplorer.hicDifferentialTAD:Matrices are not given in the "
                       "same format!\n",
                       stderr);
            return 1;
        }
        if (args.threads == 0) {
            // len(chromosome) // args.threads
            throw std::runtime_error("integer division or modulo by zero");
        }

        const tads::ContactMatrix target =
            tads::load_contact_matrix(*args.target, target_is_cooler);
        const tads::ContactMatrix control =
            tads::load_contact_matrix(*args.control, control_is_cooler);

        struct Item {
            std::size_t chromosome = 0;
            std::size_t index = 0;
        };
        std::vector<Item> items;
        // A negative --threads starts no process at all and the Python writes
        // two files holding only their headers.
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
                const std::vector<tads::Domain>& list = chromosomes[items[k].chromosome];
                const tads::TadGeometry geometry_target =
                    tads::tad_geometry(target, list, items[k].index);
                const tads::TadGeometry geometry_control =
                    tads::tad_geometry(control, list, items[k].index);
                const auto test = [&](const tads::Block& block_target,
                                      const tads::Block& block_control) {
                    const tads::DenseBlock values_target =
                        tads::extract_block(target, block_target);
                    const tads::DenseBlock values_control =
                        tads::extract_block(control, block_control);
                    return tads::rank_sum_test(values_target.values, values_control.values);
                };
                std::array<std::optional<tads::RankTest>, 3> tests;
                if (geometry_target.left.has_value()) {
                    tests[0] = test(*geometry_target.left, *geometry_control.left);
                }
                if (geometry_target.right.has_value()) {
                    tests[1] = test(*geometry_target.right, *geometry_control.right);
                }
                tests[2] = test(geometry_target.intra, geometry_control.intra);

                const double nan = std::numeric_limits<double>::quiet_NaN();
                for (std::size_t t = 0; t < 3; ++t) {
                    if (!tests[t].has_value() || std::isnan(tests[t]->pvalue)) {
                        result.pvalue[t] = nan;
                        result.statistic[t] = nan;
                        result.rejected[t] = false;
                    } else {
                        result.pvalue[t] = tests[t]->pvalue;
                        result.statistic[t] = tests[t]->statistic;
                        result.rejected[t] = tests[t]->pvalue <= args.p_value;
                    }
                }
            } catch (const std::exception& error) {
                result.error = error.what();
            }
        });
        for (const TadResult& result : results) {
            if (result.error.has_value()) {
                std::fprintf(stderr, "ERROR:hicexplorer.hicDifferentialTAD:%s\n",
                             result.error->c_str());
                return 1;
            }
        }

        const bool reject_all = args.mode_reject == "all";
        std::string accepted = header(true, args);
        std::string rejected = header(false, args);
        for (std::size_t k = 0; k < items.size(); ++k) {
            const TadResult& result = results[k];
            const bool left = result.rejected[0];
            const bool right = result.rejected[1];
            const bool intra = result.rejected[2];
            bool mask = false;
            if (args.mode == "intra-TAD") {
                mask = intra;
            } else if (args.mode == "left-inter-TAD") {
                mask = reject_all ? (left && intra) : (left || intra);
            } else if (args.mode == "right-inter-TAD") {
                mask = reject_all ? (intra && right) : (intra || right);
            } else {
                mask = reject_all ? (left && right && intra) : (left || right || intra);
            }

            const tads::Domain& domain = chromosomes[items[k].chromosome][items[k].index];
            std::string line;
            for (std::size_t column = 0; column < domain.text.size(); ++column) {
                line += domain.text[column];
                line += '\t';
            }
            for (std::size_t t = 0; t < 3; ++t) {
                line += hicx::npy::float_repr(result.pvalue[t]);
                line += '\t';
            }
            for (std::size_t t = 0; t < 3; ++t) {
                line += hicx::npy::float_repr(result.statistic[t]);
                line += t + 1 < 3 ? '\t' : '\n';
            }
            (mask ? rejected : accepted) += line;
        }

        const std::string accepted_path = args.prefix + "_accepted.diff_tad";
        const std::string rejected_path = args.prefix + "_rejected.diff_tad";
        if (!write_file(accepted_path, accepted)) {
            throw std::runtime_error("cannot write '" + accepted_path + "'");
        }
        if (!write_file(rejected_path, rejected)) {
            throw std::runtime_error("cannot write '" + rejected_path + "'");
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicDifferentialTAD: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicDifferentialTAD");
    return 0;
}
