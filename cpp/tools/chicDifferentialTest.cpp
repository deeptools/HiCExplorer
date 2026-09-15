// Port of hicexplorer/chicDifferentialTest.py.
//
// For every pair of samples and reference point of a chicAggregateStatistic
// file, each aggregated location is tested for a different interaction count
// with Fisher's exact test on the rounded up 2x2 table [[sum 1, raw 1], [sum 2,
// raw 2]], or with the chi-squared contingency test without Yates' correction
// on the same table unrounded. The results go to one HDF5 file with the
// accepted, rejected and all locations per reference point.
//
// The statistics are scipy 1.14's, routine for routine (hicx/scipy_stats.hpp):
// fisher_exact on Boost's hypergeometric distribution under scipy's policy,
// chi2_contingency with the Cephes chdtrc, and the chi-squared test's
// critical value chi2.ppf(1 - alpha, 1) with the Cephes igami.
//
// Behaviour reproduced as the Python has it:
//
//   * The chi-squared test rejects on the statistic against the critical
//     value, Fisher's test on the p-value against --alpha.
//   * A table scipy refuses (a negative value, or an expected frequency of 0
//     in the chi-squared test) is accepted with p-value 1.0 and listed with
//     NaN among all results.
//   * A reference point whose aggregated data is empty in either sample is
//     skipped, but writeResultHDF indexes the results by position in the list
//     of all reference points, so every later reference point is written with
//     the results of the next one, and the last raises IndexError after its
//     groups have been created (exit 1).
//   * --threads 0 raises ZeroDivisionError; a negative --threads starts no
//     worker, and the first reference point then raises IndexError when it is
//     written.
//   * The gene dataset holds the gene group name, not the aggregate file's
//     gene_name.
//
// --correctForMultipleTesting {none,fdr,bonferroni} (C++ only, PLAN.md 9.7
// work item 3): every p-value the tool computes, over all reference points of
// all sample pairs in computation order, is adjusted together (NaN left NaN
// and not counted), and a location is rejected when its adjusted p-value is
// at most --alpha, for both tests. A refused table stays accepted with 1.0.
// Every result group gains pvalue_adjusted_list and the file the attribute
// correctForMultipleTesting; `none` writes the Python output. The reference is
// cpp/scripts/py_chicDifferentialTest_calibrated.py.
//
// Threading: as for chicViewpoint (cpp/OPTIMIZATION.md 6) the loop is
// sequential. The Python collects its workers in order, so its output does not
// depend on --threads either.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/chic_hdf5.hpp"
#include "hicx/hdf5_util.hpp"
#include "hicx/scipy_stats.hpp"
#include "hicx/stats_ops.hpp"
#include "hicx/version.hpp"

namespace {

namespace chic = hicx::chic;
namespace h5 = hicx::h5;
namespace cli = hicx::cli;

const char* const kUsage =
    "usage: chicDifferentialTest --aggregatedFile AGGREGATEDFILE --alpha ALPHA\n"
    "                            [--outFileName OUTFILENAME]\n"
    "                            [--statisticTest {fisher,chi2}]\n"
    "                            [--threads THREADS]\n"
    "                            [--correctForMultipleTesting {none,fdr,bonferroni}]\n"
    "                            [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "chicDifferentialTest tests if two locations under consideration of the reference point "
    "have a different interaction count. For this either Fisher's test or the chi2 contingency "
    "test can be used.\n"
    "The file that is accepted for this test can be created with `chicAggregateStatistic`. H0 "
    "assumes the interactions are not different. Therefore the differential interaction counts "
    "are all where H0 was rejected.\n"
    "\n"
    "Required arguments:\n"
    "  --aggregatedFile AGGREGATEDFILE, -af AGGREGATEDFILE\n"
    "                        path to the aggregated files which should be used for\n"
    "                        the differential test.\n"
    "  --alpha ALPHA, -a ALPHA\n"
    "                        define a significance level (alpha) for accepting\n"
    "                        samples\n"
    "\n"
    "Optional arguments:\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        Output file for the differential results (Default:\n"
    "                        differentialResults.hdf5).\n"
    "  --statisticTest {fisher,chi2}\n"
    "                        Type of test used: fisher's exact test or chi2\n"
    "                        contingency (Default: fisher).\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads (uses the python multiprocessing\n"
    "                        module) (Default: 4).\n"
    "  --correctForMultipleTesting {none,fdr,bonferroni}\n"
    "                        Adjust the p-values of all tested locations of all\n"
    "                        reference points, Benjamini-Hochberg (fdr) or\n"
    "                        Bonferroni, and reject where the adjusted value is at\n"
    "                        most --alpha; the result groups gain\n"
    "                        pvalue_adjusted_list. Not in the Python tool; none\n"
    "                        gives its output (Default: none).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

class PythonError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

using Quadruple = std::vector<std::string>;

// Viewpoint.readAggregatedFileHDF: per location [chromosome, start, end,
// gene_name, sum_of_interactions, relative distance, raw], and [sum, raw].
struct Aggregated {
    std::string chromosome;
    std::vector<std::int64_t> starts;
    std::vector<std::int64_t> ends;
    std::vector<std::int64_t> relative;
    std::vector<double> raw;
    double sum = 0.0;
    [[nodiscard]] std::size_t size() const { return starts.size(); }
};

Aggregated read_aggregated(const h5::File& file, const Quadruple& path) {
    std::string internal;
    for (std::size_t i = 0; i < path.size(); ++i) {
        internal += (i > 0 ? "/" : "") + path[i];
    }
    const std::string base = "/" + internal + "/";
    const auto require = [&](const char* name) {
        if (!chic::contains(file, internal + "/" + name)) {
            throw PythonError(std::string("'NoneType' object is not subscriptable (") + name + ")");
        }
    };
    Aggregated out;
    require("chromosome");
    out.chromosome = file.read_strings(base + "chromosome").at(0);
    require("gene_name");
    for (const char* name : {"start_list", "end_list", "relative_distance_list", "raw_target_list",
                             "sum_of_interactions"}) {
        require(name);
    }
    const auto as_int = [](const std::vector<double>& values) {
        std::vector<std::int64_t> out_values;
        out_values.reserve(values.size());
        for (const double value : values) {
            out_values.push_back(static_cast<std::int64_t>(value));
        }
        return out_values;
    };
    out.starts = as_int(file.read_doubles(base + "start_list"));
    const std::vector<std::int64_t> ends = as_int(file.read_doubles(base + "end_list"));
    const std::vector<std::int64_t> relative =
        as_int(file.read_doubles(base + "relative_distance_list"));
    const std::vector<double> raw = file.read_doubles(base + "raw_target_list");
    out.sum = file.read_doubles(base + "sum_of_interactions").at(0);
    if (ends.size() < out.size() || relative.size() < out.size() || raw.size() < out.size()) {
        throw PythonError("index out of bounds");
    }
    out.ends = ends;
    out.relative = relative;
    out.raw = raw;
    return out;
}

// np.asarray(np.ceil(value), dtype=np.int64) on x86-64: out of range and
// non-finite values become INT64_MIN.
std::int64_t numpy_int64(double value) {
    const double ceiled = std::ceil(value);
    if (!std::isfinite(ceiled) || ceiled >= 9223372036854775808.0 ||
        ceiled < -9223372036854775808.0) {
        return std::numeric_limits<std::int64_t>::min();
    }
    return static_cast<std::int64_t>(ceiled);
}

// One tested location.
struct TestResult {
    double pvalue = 0.0;       // test_result[i]: NaN for a refused table
    bool refused = false;      // accepted with 1.0
    bool rejected = false;
};

struct PairResult {
    std::size_t length = 0;  // min of the two sizes
    std::vector<TestResult> tests;
};

PairResult run_tests(const Aggregated& first, const Aggregated& second, const std::string& test,
                     double alpha) {
    PairResult result;
    result.length = std::min(first.size(), second.size());
    const double critical = hicx::scipy::chi2_ppf(1.0 - alpha, 1.0);
    for (std::size_t i = 0; i < result.length; ++i) {
        TestResult entry;
        if (test == "chi2") {
            const auto chi2 = hicx::scipy::chi2_contingency_2x2(first.sum, first.raw[i],
                                                                second.sum, second.raw[i]);
            if (!chi2.has_value()) {
                entry.pvalue = std::numeric_limits<double>::quiet_NaN();
                entry.refused = true;
            } else {
                entry.pvalue = chi2->pvalue;
                entry.rejected = chi2->statistic >= critical;
            }
        } else {
            const auto pvalue = hicx::scipy::fisher_exact_pvalue(
                numpy_int64(first.sum), numpy_int64(first.raw[i]), numpy_int64(second.sum),
                numpy_int64(second.raw[i]));
            if (!pvalue.has_value()) {
                entry.pvalue = std::numeric_limits<double>::quiet_NaN();
                entry.refused = true;
            } else {
                entry.pvalue = *pvalue;
                entry.rejected = *pvalue <= alpha;
            }
        }
        result.tests.push_back(entry);
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    cli::Parser parser("chicDifferentialTest",
                       "chicDifferentialTest tests if two locations under consideration of the "
                       "reference point have a different interaction count.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--aggregatedFile", "-af"})
        .required()
        .input({"hdf5"})
        .help("path to the aggregated files which should be used for the differential test.");
    required.add({"--alpha", "-a"})
        .required()
        .type("float")
        .help("define a significance level (alpha) for accepting samples");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--outFileName", "-o"})
        .default_value("differentialResults.hdf5")
        .output({"hdf5"})
        .help("Output file for the differential results.");
    optional.add({"--statisticTest"})
        .choices({"fisher", "chi2"})
        .default_value("fisher")
        .help("Type of test used: fisher's exact test or chi2 contingency.");
    optional.add({"--threads", "-t"})
        .type("int")
        .default_value(4)
        .help("Number of threads.");
    optional.add({"--correctForMultipleTesting"})
        .choices({"none", "fdr", "bonferroni"})
        .default_value("none")
        .cpp_only("Multiple testing correction of the p-values across all reference points "
                  "(PLAN.md 9.7); none gives the Python output.")
        .help("Adjust the p-values of all tested locations, Benjamini-Hochberg (fdr) or "
              "Bonferroni, and reject where the adjusted value is at most --alpha.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    const cli::Namespace args = parser.parse(argc, argv);

    try {
        const double alpha = args.real("alpha");
        const std::string test = args.str("statisticTest");
        const std::string correction = args.str("correctForMultipleTesting");
        const h5::File file(args.str("aggregatedFile"));
        const auto root = file.attributes("/");
        const auto type = root.find("type");
        if (type == root.end()) {
            throw PythonError("KeyError: \"Can't open attribute (can't locate attribute: 'type')\"");
        }
        const auto* type_text = std::get_if<std::string>(&type->second);
        if (type_text == nullptr || *type_text != "aggregate") {
            std::fprintf(stderr, "Please provide a file created by chicAggregateStatistic for the "
                                 "parameter --aggregatedFile.\n");
            return 1;
        }

        const auto sorted_children = [&](const std::string& path) {
            std::vector<std::string> names = file.children(path);
            std::sort(names.begin(), names.end());
            return names;
        };
        std::vector<std::vector<Quadruple>> aggregated_list;
        for (const std::string& combination : file.children("/")) {
            const std::vector<std::string> matrices = file.children("/" + combination);
            if (matrices.empty()) {
                continue;
            }
            if (matrices.size() < 2) {
                throw PythonError("IndexError: list index out of range");
            }
            const std::string& matrix1 = matrices[0];
            const std::string& matrix2 = matrices[1];
            std::vector<std::string> chromosomes1 = sorted_children("/" + combination + "/" + matrix1);
            std::vector<std::string> chromosomes2 = sorted_children("/" + combination + "/" + matrix2);
            for (auto* list : {&chromosomes1, &chromosomes2}) {
                const auto genes = std::find(list->begin(), list->end(), "genes");
                if (genes == list->end()) {
                    throw PythonError("ValueError: list.remove(x): x not in list");
                }
                list->erase(genes);
            }
            for (std::size_t c = 0; c < std::min(chromosomes1.size(), chromosomes2.size()); ++c) {
                const std::vector<std::string> genes1 =
                    sorted_children("/" + combination + "/" + matrix1 + "/" + chromosomes1[c]);
                const std::vector<std::string> genes2 =
                    sorted_children("/" + combination + "/" + matrix2 + "/" + chromosomes2[c]);
                for (std::size_t g = 0; g < std::min(genes1.size(), genes2.size()); ++g) {
                    aggregated_list.push_back(
                        {{combination, matrix1, chromosomes1[c], genes1[g]},
                         {combination, matrix2, chromosomes2[c], genes2[g]}});
                }
            }
        }

        const std::int64_t threads = args.integer("threads");
        if (threads == 0) {
            throw PythonError("ZeroDivisionError: integer division or modulo by zero");
        }

        // run_statistical_tests, in order, skipping empty reference points.
        struct Computed {
            Aggregated first;
            Aggregated second;
            PairResult result;
        };
        std::vector<Computed> computed;
        if (threads > 0) {
            for (const auto& pair : aggregated_list) {
                Aggregated first = read_aggregated(file, pair[0]);
                Aggregated second = read_aggregated(file, pair[1]);
                if (first.size() == 0 || second.size() == 0) {
                    continue;
                }
                PairResult result = run_tests(first, second, test, alpha);
                computed.push_back({std::move(first), std::move(second), std::move(result)});
            }
        }

        // --correctForMultipleTesting
        std::vector<std::vector<double>> adjusted(computed.size());
        if (correction != "none") {
            std::vector<double> family;
            for (const Computed& entry : computed) {
                for (const TestResult& t : entry.result.tests) {
                    family.push_back(t.pvalue);
                }
            }
            const std::vector<double> corrected =
                correction == "fdr" ? hicx::stats::benjamini_hochberg_adjusted(family)
                                    : hicx::stats::bonferroni_adjusted(family);
            std::size_t k = 0;
            for (std::size_t n = 0; n < computed.size(); ++n) {
                for (TestResult& t : computed[n].result.tests) {
                    adjusted[n].push_back(corrected[k++]);
                    if (!t.refused) {
                        t.rejected = !std::isnan(adjusted[n].back()) && adjusted[n].back() <= alpha;
                    }
                }
            }
        }

        // writeResultHDF
        chic::Hdf5Writer writer(args.str("outFileName"));
        writer.set_attribute("/", "type", std::string("differential"));
        writer.set_attribute("/", "version", std::string(hicx::kVersion));
        writer.set_attribute("/", "alpha", alpha);
        writer.set_attribute("/", "test", test);
        if (correction != "none") {
            writer.set_attribute("/", "correctForMultipleTesting", correction);
        }
        for (std::size_t i = 0; i < aggregated_list.size(); ++i) {
            const Quadruple& input = aggregated_list[i][0];
            const std::string matrix1 = input.at(1);
            const std::string matrix2 = aggregated_list[i][1].at(1);
            const std::string& chromosome = input.at(2);
            const std::string& gene = input.at(3);
            if (!writer.exists(matrix1)) {
                writer.create_group(matrix1);
            }
            const std::string pair_group = matrix1 + "/" + matrix2;
            if (!writer.exists(pair_group)) {
                writer.create_group(pair_group);
            }
            if (!writer.exists(pair_group + "/genes")) {
                writer.create_group(pair_group + "/genes");
            }
            const std::string chromosome_group = pair_group + "/" + chromosome;
            if (!writer.exists(chromosome_group)) {
                writer.create_group(chromosome_group);
            }
            const std::string gene_group = chromosome_group + "/" + gene;
            if (writer.exists(gene_group)) {
                throw PythonError("ValueError: Unable to synchronously create group (name already "
                                  "exists)");
            }
            writer.create_group(gene_group);
            for (const char* category : {"accepted", "rejected", "all"}) {
                writer.create_group(gene_group + "/" + category);
            }
            for (const std::string category : {"accepted", "rejected", "all"}) {
                if (i >= computed.size()) {
                    throw PythonError("IndexError: list index out of range");
                }
                const Computed& entry = computed[i];
                std::vector<std::size_t> rows;
                for (std::size_t t = 0; t < entry.result.tests.size(); ++t) {
                    const TestResult& r = entry.result.tests[t];
                    const bool accepted = r.refused || !r.rejected;
                    if (category == "all" || (category == "accepted" && accepted) ||
                        (category == "rejected" && !accepted)) {
                        rows.push_back(t);
                    }
                }
                if (rows.empty()) {
                    continue;
                }
                std::vector<std::int64_t> starts;
                std::vector<std::int64_t> ends;
                std::vector<std::int64_t> relative;
                std::vector<double> raw1;
                std::vector<double> raw2;
                std::vector<double> pvalues;
                std::vector<double> adjusted_values;
                for (const std::size_t t : rows) {
                    starts.push_back(entry.first.starts[t]);
                    ends.push_back(entry.first.ends[t]);
                    relative.push_back(entry.first.relative[t]);
                    raw1.push_back(entry.first.raw[t]);
                    raw2.push_back(entry.second.raw[t]);
                    const TestResult& r = entry.result.tests[t];
                    pvalues.push_back(category != "all" && r.refused ? 1.0 : r.pvalue);
                    if (correction != "none") {
                        adjusted_values.push_back(adjusted[i][t]);
                    }
                }
                const std::string group = gene_group + "/" + category;
                writer.write_string(group + "/chromosome", entry.first.chromosome);
                writer.write_array(group + "/start_list", std::span<const std::int64_t>(starts), 9);
                writer.write_array(group + "/end_list", std::span<const std::int64_t>(ends), 9);
                writer.write_string(group + "/gene", gene);
                writer.write_array(group + "/relative_distance_list",
                                   std::span<const std::int64_t>(relative), 9);
                writer.write_scalar(group + "/sum_of_interactions_1", entry.first.sum);
                writer.write_scalar(group + "/sum_of_interactions_2", entry.second.sum);
                writer.write_array(group + "/raw_target_list_1", std::span<const double>(raw1), 9);
                writer.write_array(group + "/raw_target_list_2", std::span<const double>(raw2), 9);
                writer.write_array(group + "/pvalue_list", std::span<const double>(pvalues), 9);
                if (correction != "none") {
                    writer.write_array(group + "/pvalue_adjusted_list",
                                       std::span<const double>(adjusted_values), 9);
                }
            }
            (void)writer.hard_link(gene_group, pair_group + "/genes/" + gene);
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chicDifferentialTest: %s\n", error.what());
        return 1;
    }
    return 0;
}
