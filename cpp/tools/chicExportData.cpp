// Port of hicexplorer/chicExportData.py.
//
// Exports the cHi-C HDF5 files (interactions, significant, target, aggregate,
// differential) to text files per reference point, packed into a tar.gz or,
// with --outputMode geneName, written next to --outFileName; an interaction
// file can also go to bigWig, a file per viewpoint and one of the background
// model, through libBigWig 0.4.6 as pyBigWig 0.3.22 vendors it, with the calls
// pyBigWig makes (bwCreateHdr with 10 zoom levels, bwCreateChromList,
// bwWriteHdr, bwAddIntervals, bwClose).
//
// Behaviour reproduced as the Python has it:
//
//   * Text values follow the Python types: floats with --decimalPlaces, ints
//     and strings as they print; in aggregate and differential files the
//     stored dtype decides (numpy int64 prints as an integer).
//   * A bigWig header lists only the viewpoint's chromosome, with its size
//     from --chromosomeSizes, read up to the first blank line. The background
//     model is read with --range as the fixate range and holds the positions
//     the viewpoint has.
//   * --oneTargetFile concatenates the target files into targets.tsv, and only
//     for target files.
//   * --threads 0 raises ZeroDivisionError; a negative count, a gene name
//     that is not present, or a file type without an export all end in
//     'Contains not the requested data!' and exit 1.
//   * With --outputMode geneName the files go to the directory of
//     --outFileName, and the archive is not written.
//
// Deviation: the bigWig archive of --outputMode all lists its members in
// sorted order. The Python adds them in the os.walk order of a temporary
// directory, which is not reproducible; the comparison is by name.
//
// Threading: the export is sequential (cpp/OPTIMIZATION.md 6). The Python
// collects its workers in order, so the text archive does not depend on
// --threads either.

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

extern "C" {
#include "bigWig.h"
}

#include "hicx/argparse.hpp"
#include "hicx/chic_hdf5.hpp"
#include "hicx/chic_viewpoint.hpp"
#include "hicx/hdf5_util.hpp"
#include "hicx/tar_gz.hpp"
#include "hicx/version.hpp"

namespace {

namespace chic = hicx::chic;
namespace h5 = hicx::h5;
namespace cli = hicx::cli;
namespace fs = std::filesystem;

const char* const kUsage =
    "usage: chicExportData --file FILE [--outFileName OUTFILENAME]\n"
    "                      [--outputFileType {txt,bigwig}]\n"
    "                      [--outputMode {all,geneName}]\n"
    "                      [--outputModeName OUTPUTMODENAME]\n"
    "                      [--decimalPlaces DECIMALPLACES]\n"
    "                      [--chromosomeSizes txt file]\n"
    "                      [--backgroundModelFile BACKGROUNDMODELFILE]\n"
    "                      [--oneTargetFile] [--range RANGE RANGE]\n"
    "                      [--outputValueBigwig {relative-interactions,p-value,x-fold,raw}]\n"
    "                      [--threads THREADS] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "chicExportData exports the data stored in the intermediate hdf5 files to text files per "
    "reference point.\n"
    "\n"
    "Required arguments:\n"
    "  --file FILE, -f FILE  path to the file which should be used for data export\n"
    "\n"
    "Optional arguments:\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        Output tar.gz of the files. In case of --outputMode ==\n"
    "                        geneName it is ignored. (Default: data.tar.gz).\n"
    "  --outputFileType {txt,bigwig}, -oft {txt,bigwig}\n"
    "                        Output file type can be set for all file types to txt;\n"
    "                        except 'interaction' supports also bigwig (Default:\n"
    "                        txt).\n"
    "  --outputMode {all,geneName}, -om {all,geneName}\n"
    "                        Output mode: Either all date is written or a gene name\n"
    "                        must be specified. (Default: all).\n"
    "  --outputModeName OUTPUTMODENAME, -omn OUTPUTMODENAME\n"
    "                        ONLY valid if --outputMode geneName! Define the name\n"
    "                        of the gene\n"
    "  --decimalPlaces DECIMALPLACES\n"
    "                        Decimal places for all output floating numbers in the\n"
    "                        viewpoint files (Default: 12).\n"
    "  --chromosomeSizes txt file, -cs txt file\n"
    "                        File with the chromosome sizes for your genome. A tab-\n"
    "                        delimited two column layout \"chr_name size\" is\n"
    "                        expectedUsually the sizes can be determined from the\n"
    "                        SAM/BAM input files, however, for cHi-C or scHi-C it\n"
    "                        can be that at the start or end no data is present.\n"
    "                        Please consider that this option causes that only\n"
    "                        reads are considered which are on the listed\n"
    "                        chromosomes.Use this option to guarantee fixed sizes.\n"
    "                        An example file is available via UCSC: http://hgdownlo\n"
    "                        ad.soe.ucsc.edu/goldenPath/dm3/bigZips/dm3.chrom.sizes\n"
    "  --backgroundModelFile BACKGROUNDMODELFILE, -bmf BACKGROUNDMODELFILE\n"
    "                        Path to the background model file. Required only for\n"
    "                        fileType=interactions and outputFileTypeBigwig.\n"
    "  --oneTargetFile, -otf\n"
    "                        Compile all target files to one. Applies only if\n"
    "                        --fileType is target\n"
    "  --range RANGE RANGE   Defines the region upstream and downstream of a\n"
    "                        reference point which should be included. Format is\n"
    "                        --range upstream downstream, e.g.: --range 500000\n"
    "                        500000 plots 500kb up- and 500kb downstream. This\n"
    "                        value should not exceed the range used in the other\n"
    "                        chic-tools. Applies only for interaction files in the\n"
    "                        combination with bigwig and a background model file!\n"
    "  --outputValueBigwig {relative-interactions,p-value,x-fold,raw}, -ovb {relative-interactions,p-value,x-fold,raw}\n"
    "                        Select which value the bigwig file should contain:\n"
    "                        'relative-interactions', 'p-value', 'x-fold', 'raw'\n"
    "                        (Default: relative-interactions).\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads (uses the python multiprocessing\n"
    "                        module) (Default: 4).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

class PythonError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

// Exit 1 with the Python's log message.
class LoggedExit : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

using Path = std::vector<std::string>;

std::string join(const Path& parts, const std::string& separator) {
    std::string out;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        out += (i > 0 ? separator : "") + parts[i];
    }
    return out;
}

std::string fixed(double value, std::int64_t decimals) {
    if (decimals < 0) {
        throw PythonError("ValueError: Format specifier missing precision");
    }
    return chic::format_fixed(value, static_cast<int>(decimals));
}

// A numeric dataset read with its dtype: an integer dtype prints as an
// integer, a float dtype with the decimal places.
struct Column {
    std::vector<double> values;
    bool integral = false;
};

Column read_column(const h5::File& file, const std::string& path) {
    Column column;
    column.values = file.read_doubles(path);
    const std::string dtype = file.dataset_dtype(path);
    column.integral = dtype.find("int") != std::string::npos;
    return column;
}

std::string cell(const Column& column, std::size_t i, std::int64_t decimals) {
    const double value = column.values.at(i);
    if (column.integral) {
        return std::to_string(static_cast<std::int64_t>(value));
    }
    return fixed(value, decimals);
}

std::vector<std::string> sorted_children(const h5::File& file, const std::string& path) {
    std::vector<std::string> names = file.children(path);
    std::sort(names.begin(), names.end());
    return names;
}

std::vector<std::string> without_genes(std::vector<std::string> names) {
    const auto genes = std::find(names.begin(), names.end(), "genes");
    if (genes == names.end()) {
        throw PythonError("ValueError: list.remove(x): x not in list");
    }
    names.erase(genes);
    return names;
}

std::vector<std::string> genes_of(const h5::File& file, const std::string& group) {
    if (!chic::contains(file, group.substr(1) + "/genes")) {
        throw PythonError("KeyError: \"Unable to synchronously open object (object 'genes' doesn't "
                          "exist)\"");
    }
    return sorted_children(file, group + "/genes");
}

// The while loop over gene_name, gene_name_1, gene_name_2, ...
std::vector<std::string> numbered(const std::vector<std::string>& present,
                                  const std::string& name,
                                  const std::vector<std::string>* also = nullptr) {
    std::vector<std::string> out;
    std::string gene = name;
    for (int counter = 1;; ++counter) {
        const bool here = std::find(present.begin(), present.end(), gene) != present.end();
        const bool there =
            also == nullptr || std::find(also->begin(), also->end(), gene) != also->end();
        if (!here || !there) {
            return out;
        }
        out.push_back(gene);
        gene = name + "_" + std::to_string(counter);
    }
}

struct Entries {
    std::vector<std::string> chromosomes;
    std::vector<std::int64_t> starts;
    std::vector<std::int64_t> ends;
    std::vector<double> values;
};

struct BigWigFile {
    std::string name;
    std::string chromosome;
    std::int64_t size = 0;
    Entries entries;
};

// pyBigWig: addHeader, addEntries (validated as addEntriesInputOK does), close.
void write_bigwig(const std::string& path, const BigWigFile& file) {
    bigWigFile_t* bw = bwOpen(const_cast<char*>(path.c_str()), nullptr, "w");
    if (bw == nullptr) {
        throw PythonError("RuntimeError: Received an error during file opening!");
    }
    struct Closer {
        bigWigFile_t* bw;
        ~Closer() { bwClose(bw); }
    } closer{bw};
    if (file.size > 0xFFFFFFFFLL) {
        throw PythonError("RuntimeError: Length out of bounds for a bigWig file!");
    }
    if (bwCreateHdr(bw, 10)) {
        throw PythonError("RuntimeError: Received an error in bwCreateHdr");
    }
    char* chrom = const_cast<char*>(file.chromosome.c_str());
    uint32_t length = static_cast<uint32_t>(file.size);
    bw->cl = bwCreateChromList(&chrom, &length, 1);
    if (bw->cl == nullptr) {
        throw PythonError("RuntimeError: Received an error in bwCreateChromList");
    }
    if (bwWriteHdr(bw)) {
        throw PythonError("RuntimeError: Received an error while writing the bigWig header");
    }
    const Entries& e = file.entries;
    const std::size_t n = e.starts.size();
    const auto invalid = [] {
        return PythonError("RuntimeError: The entries you tried to add are out of order, precede "
                           "already added entries, or otherwise use illegal values.");
    };
    if (n == 0) {
        throw invalid();
    }
    std::vector<char*> chroms(n);
    std::vector<uint32_t> starts(n);
    std::vector<uint32_t> ends(n);
    std::vector<float> values(n);
    uint32_t last_tid = static_cast<uint32_t>(-1);
    uint32_t last_end = 0;
    for (std::size_t i = 0; i < n; ++i) {
        chroms[i] = const_cast<char*>(e.chromosomes[i].c_str());
        const uint32_t tid = bwGetTid(bw, chroms[i]);
        if (tid == static_cast<uint32_t>(-1)) {
            throw invalid();
        }
        if (e.starts[i] > 0xFFFFFFFFLL || e.ends[i] > 0xFFFFFFFFLL) {
            throw PythonError("RuntimeError: Length out of bounds for a bigWig file!");
        }
        starts[i] = static_cast<uint32_t>(e.starts[i]);
        ends[i] = static_cast<uint32_t>(e.ends[i]);
        if (starts[i] >= ends[i]) {
            throw invalid();
        }
        if (last_tid != static_cast<uint32_t>(-1)) {
            if (last_tid > tid || (last_tid == tid && starts[i] < last_end)) {
                throw invalid();
            }
        }
        last_tid = tid;
        last_end = ends[i];
        values[i] = static_cast<float>(e.values[i]);
    }
    if (bwAddIntervals(bw, chroms.data(), starts.data(), ends.data(), values.data(),
                       static_cast<uint32_t>(n))) {
        throw PythonError("RuntimeError: Received an error while adding the intervals.");
    }
}

}  // namespace

int main(int argc, char** argv) {
    cli::Parser parser("chicExportData",
                       "chicExportData exports the data stored in the intermediate hdf5 files to "
                       "text files per reference point.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--file", "-f"})
        .required()
        .input({"hdf5"})
        .help("path to the file which should be used for data export");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--outFileName", "-o"})
        .default_value("data.tar.gz")
        .output({"tar.gz"})
        .help("Output tar.gz of the files. In case of --outputMode == geneName it is ignored.");
    optional.add({"--outputFileType", "-oft"})
        .default_value("txt")
        .choices({"txt", "bigwig"})
        .help("Output file type can be set for all file types to txt; except 'interaction' "
              "supports also bigwig.");
    optional.add({"--outputMode", "-om"})
        .default_value("all")
        .choices({"all", "geneName"})
        .help("Output mode: Either all date is written or a gene name must be specified.");
    optional.add({"--outputModeName", "-omn"})
        .help("ONLY valid if --outputMode geneName! Define the name of the gene");
    optional.add({"--decimalPlaces"})
        .type("int")
        .default_value(12)
        .help("Decimal places for all output floating numbers in the viewpoint files.");
    optional.add({"--chromosomeSizes", "-cs"})
        .file_type("r")
        .help("File with the chromosome sizes for your genome.");
    optional.add({"--backgroundModelFile", "-bmf"})
        .input({"txt"})
        .help("Path to the background model file.");
    optional.add({"--oneTargetFile", "-otf"})
        .action(cli::Action::StoreTrue)
        .help("Compile all target files to one. Applies only if --fileType is target");
    optional.add({"--range"})
        .type("int")
        .nargs(2)
        .help("Defines the region upstream and downstream of a reference point which should be "
              "included.");
    optional.add({"--outputValueBigwig", "-ovb"})
        .default_value("relative-interactions")
        .choices({"relative-interactions", "p-value", "x-fold", "raw"})
        .help("Select which value the bigwig file should contain.");
    optional.add({"--threads", "-t"})
        .type("int")
        .default_value(4)
        .help("Number of threads.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    const cli::Namespace args = parser.parse(argc, argv);

    try {
        const std::string input = args.str("file");
        std::string out_file_name = args.str("outFileName");
        const std::string output_type = args.str("outputFileType");
        const std::string output_mode = args.str("outputMode");
        const std::optional<std::string> mode_name = args.opt_str("outputModeName");
        const std::int64_t decimals = args.integer("decimalPlaces");
        const std::optional<std::string> background_file = args.opt_str("backgroundModelFile");
        const std::string value_kind = args.str("outputValueBigwig");

        const h5::File file(input);
        const auto root = file.attributes("/");
        const auto type_attribute = root.find("type");
        if (type_attribute == root.end()) {
            throw PythonError("KeyError: \"Can't open attribute (can't locate attribute: 'type')\"");
        }
        const std::string file_type = std::holds_alternative<std::string>(type_attribute->second)
                                          ? std::get<std::string>(type_attribute->second)
                                          : std::string();

        if (output_mode == "geneName" && !mode_name.has_value()) {
            throw LoggedExit("Output mode is 'geneName'. Please specify a gene name via "
                             "--outputModeName too!");
        }

        std::optional<chic::BackgroundModel> background;
        std::map<std::string, std::int64_t> chromosome_sizes;
        std::vector<std::int64_t> range;
        if (output_type == "bigwig") {
            if (file_type != "interactions") {
                throw LoggedExit("Only file type 'interactions' supports bigwig. Exiting.");
            }
            if (!args.given("range")) {
                throw LoggedExit("Bigwig files require the argument '--range upstream downstream'. "
                                 "Exiting.");
            }
            range = args.integers("range");
            if (background_file.has_value()) {
                if (!background_file->empty()) {
                    background = chic::read_background_model(*background_file, range.at(0),
                                                             range.at(1), range.at(1), true);
                }
            } else {
                throw LoggedExit("Please define a background file via --backgroundModelFile.");
            }
            const std::optional<std::string> sizes = args.opt_str("chromosomeSizes");
            if (!sizes.has_value()) {
                throw LoggedExit("Bigwig files require the argument '--chromosomeSizes'. Exiting.");
            }
            for (const std::string& raw : chic::read_lines(*sizes)) {
                const std::string line(chic::strip(raw));
                if (line.empty()) {
                    break;
                }
                const std::size_t tab = line.find('\t');
                if (tab == std::string::npos) {
                    throw PythonError("IndexError: list index out of range");
                }
                const std::size_t next = line.find('\t', tab + 1);
                chromosome_sizes[line.substr(0, tab)] = chic::python_int(
                    line.substr(tab + 1, next == std::string::npos ? std::string::npos
                                                                   : next - tab - 1));
            }
        }

        // The file list per file type.
        const std::vector<std::string> keys = file.children("/");
        std::vector<std::vector<Path>> file_list;  // each entry: one or two paths
        if (file_type == "interactions" || file_type == "significant") {
            for (const std::string& sample : keys) {
                if (output_mode == "all") {
                    for (const std::string& chromosome :
                         without_genes(sorted_children(file, "/" + sample))) {
                        for (const std::string& gene :
                             sorted_children(file, "/" + sample + "/" + chromosome)) {
                            file_list.push_back({{sample, chromosome, gene}});
                        }
                    }
                } else {
                    for (const std::string& gene :
                         numbered(genes_of(file, "/" + sample), *mode_name)) {
                        file_list.push_back({{sample, "genes", gene}});
                    }
                }
            }
        } else if (file_type == "target") {
            const auto mode = root.find("combinationMode");
            if (mode == root.end()) {
                throw PythonError(
                    "KeyError: \"Can't open attribute (can't locate attribute: 'combinationMode')\"");
            }
            const std::string combination = std::holds_alternative<std::string>(mode->second)
                                                ? std::get<std::string>(mode->second)
                                                : std::string();
            if (combination == "dual") {
                for (const std::string& outer : keys) {
                    for (const std::string& inner : file.children("/" + outer)) {
                        const std::vector<std::string> genes = file.children(
                            "/" + outer + "/" + inner + "/genes");
                        if (output_mode == "all") {
                            for (const std::string& gene : genes) {
                                file_list.push_back({{outer, inner, "genes", gene}});
                            }
                        } else {
                            for (const std::string& gene : numbered(genes, *mode_name)) {
                                file_list.push_back({{outer, inner, "genes", gene}});
                            }
                        }
                    }
                }
            } else if (combination == "single") {
                for (const std::string& outer : keys) {
                    const std::vector<std::string> genes = file.children("/" + outer + "/genes");
                    const std::vector<std::string> selected =
                        output_mode == "all" ? genes : numbered(genes, *mode_name);
                    for (const std::string& gene : selected) {
                        file_list.push_back({{outer, "genes", gene}});
                    }
                }
            }
        } else if (file_type == "aggregate") {
            for (const std::string& combination : keys) {
                const std::vector<std::string> matrices = file.children("/" + combination);
                if (matrices.empty()) {
                    continue;
                }
                if (matrices.size() < 2) {
                    throw PythonError("IndexError: list index out of range");
                }
                const std::string base1 = "/" + combination + "/" + matrices[0];
                const std::string base2 = "/" + combination + "/" + matrices[1];
                if (output_mode == "all") {
                    const std::vector<std::string> chromosomes1 =
                        without_genes(sorted_children(file, base1));
                    const std::vector<std::string> chromosomes2 =
                        without_genes(sorted_children(file, base2));
                    for (std::size_t c = 0; c < std::min(chromosomes1.size(), chromosomes2.size());
                         ++c) {
                        const auto genes1 = sorted_children(file, base1 + "/" + chromosomes1[c]);
                        const auto genes2 = sorted_children(file, base2 + "/" + chromosomes2[c]);
                        for (std::size_t g = 0; g < std::min(genes1.size(), genes2.size()); ++g) {
                            file_list.push_back(
                                {{combination, matrices[0], chromosomes1[c], genes1[g]},
                                 {combination, matrices[1], chromosomes2[c], genes2[g]}});
                        }
                    }
                } else {
                    const std::vector<std::string> genes1 = genes_of(file, base1);
                    const std::vector<std::string> genes2 = genes_of(file, base2);
                    for (const std::string& gene : numbered(genes1, *mode_name, &genes2)) {
                        file_list.push_back({{combination, matrices[0], "genes", gene},
                                             {combination, matrices[1], "genes", gene}});
                    }
                }
            }
        } else if (file_type == "differential") {
            for (const std::string& outer : keys) {
                for (const std::string& inner : file.children("/" + outer)) {
                    const std::string base = "/" + outer + "/" + inner;
                    if (output_mode == "all") {
                        for (const std::string& chromosome :
                             without_genes(sorted_children(file, base))) {
                            for (const std::string& gene :
                                 sorted_children(file, base + "/" + chromosome)) {
                                file_list.push_back({{outer, inner, chromosome, gene}});
                            }
                        }
                    } else {
                        for (const std::string& gene : numbered(genes_of(file, base), *mode_name)) {
                            file_list.push_back({{outer, inner, "genes", gene}});
                        }
                    }
                }
            }
        }

        const std::int64_t threads = args.integer("threads");
        if (threads == 0) {
            throw PythonError("ZeroDivisionError: integer division or modulo by zero");
        }

        // exportData, in order.
        std::vector<std::string> names;
        std::vector<std::string> contents;
        std::vector<std::vector<BigWigFile>> bigwigs;
        if (threads > 0) {
            if (file_type == "interactions" || file_type == "significant") {
                const std::string header =
                    "# Chromosome\tStart\tEnd\tGene\tSum of interactions\tRelative position\t"
                    "Relative Interactions\tp-value\tx-fold\tRaw\n";
                for (const auto& entry : file_list) {
                    const Path& sample = entry.front();
                    const chic::InteractionTable table = chic::read_interaction_table(file, sample);
                    std::vector<double> sorted_keys = table.keys;
                    std::sort(sorted_keys.begin(), sorted_keys.end());
                    if (output_type == "txt") {
                        std::string content = header;
                        for (const double key : sorted_keys) {
                            const chic::InteractionRecord& r = table.records.at(key);
                            content += r.chromosome + "\t" + std::to_string(r.start) + "\t" +
                                       std::to_string(r.end) + "\t" + r.gene + "\t" +
                                       fixed(r.sum_of_interactions, decimals) + "\t" +
                                       std::to_string(r.relative_position) + "\t" +
                                       fixed(r.interaction, decimals) + "\t" +
                                       fixed(r.pvalue, decimals) + "\t" + fixed(r.xfold, decimals) +
                                       "\t" + fixed(r.raw, decimals) + "\n";
                        }
                        contents.push_back(std::move(content));
                        names.push_back(join(sample, "_") + "_" + file_type + ".txt");
                    } else {
                        BigWigFile viewpoint;
                        std::map<std::int64_t, std::pair<std::int64_t, std::int64_t>> relative;
                        std::map<std::int64_t, std::string> relative_chromosome;
                        for (const double key : sorted_keys) {
                            const chic::InteractionRecord& r = table.records.at(key);
                            viewpoint.entries.chromosomes.push_back(r.chromosome);
                            viewpoint.entries.starts.push_back(r.start);
                            viewpoint.entries.ends.push_back(r.end);
                            viewpoint.entries.values.push_back(
                                value_kind == "relative-interactions" ? r.interaction
                                : value_kind == "p-value"             ? r.pvalue
                                : value_kind == "x-fold"              ? r.xfold
                                                                      : r.raw);
                            relative[r.relative_position] = {r.start, r.end};
                            relative_chromosome[r.relative_position] = r.chromosome;
                        }
                        if (viewpoint.entries.chromosomes.empty()) {
                            throw PythonError("IndexError: list index out of range");
                        }
                        viewpoint.chromosome = viewpoint.entries.chromosomes.front();
                        const auto size = chromosome_sizes.find(viewpoint.chromosome);
                        if (size == chromosome_sizes.end()) {
                            throw PythonError("KeyError: '" + viewpoint.chromosome + "'");
                        }
                        viewpoint.size = size->second;
                        viewpoint.name = join(sample, "_") + ".bigwig";
                        if (!background.has_value()) {
                            throw PythonError(
                                "AttributeError: 'NoneType' object has no attribute 'keys'");
                        }
                        BigWigFile model;
                        model.chromosome = viewpoint.chromosome;
                        model.size = viewpoint.size;
                        model.name = "background_" + join(sample, "_") + "_" + file_type + ".bigwig";
                        for (const std::int64_t key : background->sorted_keys()) {
                            const auto it = relative.find(key);
                            if (it != relative.end()) {
                                model.entries.chromosomes.push_back(relative_chromosome.at(key));
                                model.entries.starts.push_back(it->second.first);
                                model.entries.ends.push_back(it->second.second);
                                model.entries.values.push_back(background->at(key).at(0));
                            }
                        }
                        bigwigs.push_back({std::move(viewpoint), std::move(model)});
                    }
                }
            } else if (file_type == "target") {
                for (const auto& entry : file_list) {
                    const Path& target = entry.front();
                    const std::string joined = join(target, "/");
                    if (!chic::contains(file, joined)) {
                        throw PythonError("KeyError: \"Unable to synchronously open object\"");
                    }
                    const std::string base = "/" + joined + "/";
                    const std::string chromosome = file.read_strings(base + "chromosome").at(0);
                    const std::vector<std::string> starts = file.read_strings(base + "start_list");
                    const std::vector<std::string> ends = file.read_strings(base + "end_list");
                    std::string content;
                    for (std::size_t i = 0; i < std::min(starts.size(), ends.size()); ++i) {
                        content += chromosome + "\t" + starts[i] + "\t" + ends[i] + "\n";
                    }
                    contents.push_back(std::move(content));
                    names.push_back(join(target, "_") + "_target.txt");
                }
            } else if (file_type == "aggregate") {
                const std::string header =
                    "# Chromosome\tStart\tEnd\tGene\tSum of interactions\tRelative position\tRaw\n";
                for (const auto& entry : file_list) {
                    for (const Path& sample : entry) {
                        const std::string base = "/" + join(sample, "/") + "/";
                        for (const char* name : {"chromosome", "gene_name", "start_list", "end_list",
                                                 "relative_distance_list", "raw_target_list",
                                                 "sum_of_interactions"}) {
                            if (!file.exists(base + name)) {
                                throw PythonError(std::string("'NoneType' object (") + name + ")");
                            }
                        }
                        const std::string chromosome = file.read_strings(base + "chromosome").at(0);
                        const std::string gene = file.read_strings(base + "gene_name").at(0);
                        const Column starts = read_column(file, base + "start_list");
                        const Column ends = read_column(file, base + "end_list");
                        const Column relative = read_column(file, base + "relative_distance_list");
                        const Column raw = read_column(file, base + "raw_target_list");
                        const double sum = file.read_doubles(base + "sum_of_interactions").at(0);
                        std::string content = header;
                        for (std::size_t i = 0; i < starts.values.size(); ++i) {
                            content += chromosome + "\t" + cell(starts, i, decimals) + "\t" +
                                       cell(ends, i, decimals) + "\t" + gene + "\t" +
                                       fixed(sum, decimals) + "\t" + cell(relative, i, decimals) +
                                       "\t" + cell(raw, i, decimals) + "\n";
                        }
                        contents.push_back(std::move(content));
                        names.push_back(join(sample, "_") + "_aggregate.txt");
                    }
                }
            } else if (file_type == "differential") {
                const std::string header =
                    "# Chromosome\tStart\tEnd\tGene\tRelative distance\tsum of interactions 1\t"
                    "target_1 raw\tsum of interactions 2\ttarget_2 raw\tp-value\n";
                for (const auto& entry : file_list) {
                    const Path& quadruple = entry.front();
                    const char* const items[3] = {"accepted", "all", "rejected"};
                    for (const char* item : items) {
                        const std::string base = "/" + join(quadruple, "/") + "/" + item + "/";
                        std::string content = header;
                        bool complete = chic::contains(file, join(quadruple, "/") + "/" + item);
                        for (const char* name : {"chromosome", "start_list", "end_list",
                                                 "relative_distance_list", "gene", "pvalue_list",
                                                 "raw_target_list_1", "raw_target_list_2",
                                                 "sum_of_interactions_1", "sum_of_interactions_2"}) {
                            complete = complete && file.exists(base + name);
                        }
                        if (complete) {
                            const std::string chromosome =
                                file.read_strings(base + "chromosome").at(0);
                            const std::string gene = file.read_strings(base + "gene").at(0);
                            const Column starts = read_column(file, base + "start_list");
                            const Column ends = read_column(file, base + "end_list");
                            const Column relative = read_column(file, base + "relative_distance_list");
                            const Column pvalues = read_column(file, base + "pvalue_list");
                            const Column raw1 = read_column(file, base + "raw_target_list_1");
                            const Column raw2 = read_column(file, base + "raw_target_list_2");
                            const double sum1 = file.read_doubles(base + "sum_of_interactions_1").at(0);
                            const double sum2 = file.read_doubles(base + "sum_of_interactions_2").at(0);
                            const std::size_t length = std::min(
                                {starts.values.size(), ends.values.size(), relative.values.size(),
                                 raw1.values.size(), raw2.values.size(), pvalues.values.size()});
                            for (std::size_t i = 0; i < length; ++i) {
                                content += chromosome + "\t" + cell(starts, i, decimals) + "\t" +
                                           cell(ends, i, decimals) + "\t" + gene + "\t" +
                                           cell(relative, i, decimals) + "\t" +
                                           fixed(sum1, decimals) + "\t" + cell(raw1, i, decimals) +
                                           "\t" + fixed(sum2, decimals) + "\t" +
                                           cell(raw2, i, decimals) + "\t" +
                                           cell(pvalues, i, decimals) + "\n";
                            }
                        }
                        contents.push_back(std::move(content));
                        names.push_back(join(quadruple, "_") + "_" + item + "_differential.txt");
                    }
                }
            }
        }

        if (contents.empty() && bigwigs.empty()) {
            throw LoggedExit("Contains not the requested data!");
        }
        // tarfile stamps every member with time.time(); the port stamps 0, so
        // that two runs write the same archive. Nothing reads the time.
        const std::int64_t now = 0;
        if (output_type == "txt") {
            if (output_mode == "geneName") {
                const std::string basepath = fs::path(out_file_name).parent_path().string();
                for (std::size_t i = 0; i < contents.size(); ++i) {
                    const std::string target = basepath + "/" + names[i];
                    std::ofstream out(target, std::ios::binary);
                    if (!out) {
                        throw PythonError("FileNotFoundError: [Errno 2] No such file or directory: '" +
                                          target + "'");
                    }
                    out << contents[i];
                }
            } else {
                hicx::TarGzWriter tar(out_file_name);
                if (args.flag("oneTargetFile") && file_type == "target") {
                    std::string all;
                    for (const std::string& content : contents) {
                        all += content;
                    }
                    tar.add("targets.tsv", all, now);
                } else {
                    for (std::size_t i = 0; i < contents.size(); ++i) {
                        tar.add(names[i], contents[i], now);
                    }
                }
                tar.close();
            }
        } else {
            if (bwInit(128000)) {
                throw PythonError("RuntimeError: bwInit failed");
            }
            fs::path folder;
            if (output_mode == "geneName") {
                folder = fs::path(out_file_name).parent_path();
            } else {
                std::string pattern = (fs::temp_directory_path() / "bigwig_folderXXXXXX").string();
                if (::mkdtemp(pattern.data()) == nullptr) {
                    throw PythonError("OSError: cannot create a temporary directory");
                }
                folder = pattern;
            }
            std::vector<std::string> written;
            try {
                for (const auto& group : bigwigs) {
                    for (const BigWigFile& bigwig : group) {
                        const std::string target = (folder.string()) + "/" + bigwig.name;
                        write_bigwig(target, bigwig);
                        if (std::find(written.begin(), written.end(), bigwig.name) == written.end()) {
                            written.push_back(bigwig.name);
                        }
                    }
                }
                if (output_mode == "all") {
                    if (out_file_name.size() < 7 ||
                        out_file_name.compare(out_file_name.size() - 7, 7, ".tar.gz") != 0) {
                        out_file_name += ".tar.gz";
                    }
                    std::sort(written.begin(), written.end());
                    hicx::TarGzWriter tar(out_file_name);
                    for (const std::string& name : written) {
                        tar.add_file(name, (folder / name).string(), now);
                    }
                    tar.close();
                }
            } catch (...) {
                if (output_mode == "all") {
                    std::error_code ignored;
                    fs::remove_all(folder, ignored);
                }
                bwCleanup();
                throw;
            }
            if (output_mode == "all") {
                std::error_code ignored;
                fs::remove_all(folder, ignored);
            }
            bwCleanup();
        }
    } catch (const LoggedExit& error) {
        std::fprintf(stderr, "ERROR:hicexplorer.chicExportData:%s\n", error.what());
        return 1;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chicExportData: %s\n", error.what());
        return 1;
    }
    return 0;
}
