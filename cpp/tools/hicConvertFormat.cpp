// Port of hicexplorer/hicConvertFormat.py.
//
// Converts a Hi-C contact matrix between the formats HiCExplorer understands.
// The tool is a thin driver over the file layer: load through the reader for
// the input format, hand the five loader values to the writer for the output
// format, save. Only two paths do anything else, and both are reproduced here
// rather than in the file layer because that is where the Python has them:
//
//   * homer and ginteractions output first replaces the matrix with
//     triu(m).maximum(triu(m).T) (hicConvertFormat.py:258-260)
//   * mcool output rebuilds the matrix through hiCMatrix.setMatrix and merges
//     its bins once per requested resolution (:294-329)
//
// .hic files (cpp/PLAN.md tier 9, item 9.1) go through hicfilecpp and
// core/src/hic_adapter.cpp:
//
//   * --inputFormat hic --outputFormat cool is hic2cool_convert, as in the
//     Python (hicConvertFormat.py:124-138): every resolution into one mcool
//     file without --resolutions, one cool file per resolution with it.
//   * --inputFormat hic into mcool, h5, homer, ginteractions and hicpro goes
//     beyond the Python, which refuses them. mcool is hic2cool's multi
//     resolution layout of every resolution or of those in --resolutions. The
//     other formats convert the single resolution in --resolutions into a
//     temporary cool file with hic2cool and continue as --inputFormat cool,
//     which is the Python route in two steps.
//   * --outputFormat hic, also beyond the Python, writes the loaded matrix as
//     a Juicer .hic file of version --hicVersion with the normalizations of
//     --hicNormalizations; coarser --resolutions are binned from the matrix.
//     A cool file with a /resolutions group given to --inputFormat cool is
//     written with every resolution, or with those in --resolutions.
//
// Not supported, and refused rather than half done:
//
//   * --chromosome. Loading a single chromosome out of a cooler is a distinct
//     cooler code path (cool.py:119-151) that the C++ cool reader does not
//     have yet.
//
// Defects reproduced on purpose, each pinned by
// hicexplorer/test/general/test_hicConvertFormat.py:
//
//   * --enforce_integer rounds an already corrected matrix to all zeros and
//     writes a cool file of zero valued pixels, exit status 0, no warning.
//     See the comment at apply_enforce_integer_note below.
//   * --outputFormat ginteractions writes <outFileName>.tsv and leaves
//     <outFileName> untouched.
//   * hicpro output writes all stored entries, not the upper triangle.
//   * the hic_metadata dictionary is consumed in place, so of several mcool
//     resolutions only the first carries genome-assembly in its metadata
//     attribute.

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/hdf5_util.hpp"
#include "hicx/hic_adapter.hpp"
#include "hicx/matrix_data.hpp"
#include "hicx/reduce_matrix.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/text_formats.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicConvertFormat --matrices MATRICES [MATRICES ...] --outFileName\n"
    "                        OUTFILENAME [OUTFILENAME ...] --inputFormat\n"
    "                        {h5,cool,hic,homer,hicpro,2D-text} --outputFormat\n"
    "                        {cool,h5,homer,ginteractions,mcool,hicpro,hic}\n"
    "                        [--correction_name CORRECTION_NAME]\n"
    "                        [--correction_division] [--store_applied_correction]\n"
    "                        [--chromosome CHROMOSOME] [--enforce_integer]\n"
    "                        [--load_raw_values] [--resolutions RESOLUTIONS "
    "[RESOLUTIONS ...]]\n"
    "                        [--help] [--chromosomeSizes txt file] [--version]\n"
    "                        [--bedFileHicpro BEDFILEHICPRO [BEDFILEHICPRO ...]]\n"
    "                        [--hicVersion {8,9}]\n"
    "                        [--hicNormalizations {VC,VC_SQRT,KR,SCALE,none} "
    "[{VC,VC_SQRT,KR,SCALE,none} ...]]\n"
    "                        [--threads THREADS]\n";

const char* const kHelp =
    "\n"
    "Conversion of Hi-C matrices of different file formats. We support the\n"
    "conversion of hic to cool format via hic2cool, and homer, HicPro, h5 and cool\n"
    "format to h5, cool, homer or ginteractions format. Moreover, hicConvertFormat\n"
    "accepts multiple input files from one format with different resolutions and\n"
    "creates a mcool file. Each original file is stored under the path e.g.\n"
    "::/resolutions/10000. A batch computation is possible, the number of input\n"
    "files and output files needs to match, all input files need to be of the same\n"
    "format type and all output files too. For input and output of cooler files\n"
    "special options are available, for all other formats they will be ignored.\n"
    "HiCPro file format needs an additional bed file as input.\n"
    "\n"
    "Required arguments:\n"
    "  --matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]\n"
    "                        input file(s). Could be one or many files.\n"
    "  --outFileName OUTFILENAME [OUTFILENAME ...], -o OUTFILENAME [OUTFILENAME ...]\n"
    "                        File name to save the exported matrix.\n"
    "  --inputFormat {h5,cool,hic,homer,hicpro,2D-text}\n"
    "                        File format of the input matrix file.\n"
    "  --outputFormat {cool,h5,homer,ginteractions,mcool,hicpro,hic}\n"
    "                        Output format. (Default: cool).\n"
    "\n"
    "Optional arguments:\n"
    "  --correction_name CORRECTION_NAME\n"
    "                        Name of the column which stores the correction\n"
    "                        factors. Option only for cool input files.\n"
    "                        (Default: weight).\n"
    "  --correction_division\n"
    "                        If set, division is applied for correction.\n"
    "  --store_applied_correction\n"
    "                        Store the applied correction and do not set\n"
    "                        correction factors.\n"
    "  --chromosome CHROMOSOME\n"
    "                        Load only one chromosome.\n"
    "  --enforce_integer     Enforce datatype of counts to integer.\n"
    "  --load_raw_values     Load only 'count' data and do not apply a correction.\n"
    "  --resolutions RESOLUTIONS [RESOLUTIONS ...], -r RESOLUTIONS [RESOLUTIONS ...]\n"
    "                        List of resolutions that should be added.\n"
    "  --help, -h            show this help message and exit.\n"
    "  --chromosomeSizes txt file, -cs txt file\n"
    "                        For the input format `2D-text` only.\n"
    "  --version             show program's version number and exit\n"
    "  --bedFileHicpro BEDFILEHICPRO [BEDFILEHICPRO ...], -bf BEDFILEHICPRO "
    "[BEDFILEHICPRO ...]\n"
    "                        Bed file(s) of hicpro file format.\n"
    "  --hicVersion {8,9}    Version of a .hic output file. (Default: 8).\n"
    "  --hicNormalizations {VC,VC_SQRT,KR,SCALE,none} [{VC,VC_SQRT,KR,SCALE,none} ...]\n"
    "                        Normalizations a .hic output file stores, computed\n"
    "                        as Juicer tools addNorm does. (Default: VC VC_SQRT KR\n"
    "                        SCALE).\n"
    "  --threads THREADS     Threads for compressing a .hic output file; the file\n"
    "                        does not depend on the number. (Default: 1).\n";

struct Arguments {
    std::vector<std::string> matrices;
    std::vector<std::string> out_file_names;
    std::vector<std::string> resolutions;
    std::vector<std::string> bed_file_hicpro;
    std::string input_format;
    std::string output_format = "cool";
    std::string correction_name = "weight";
    std::string chromosome;
    std::string chromosome_sizes;
    bool correction_division = false;
    bool store_applied_correction = false;
    bool enforce_integer = false;
    bool load_raw_values = false;
    std::string hic_version = "8";
    std::vector<std::string> hic_normalizations{"VC", "VC_SQRT", "KR", "SCALE"};
    bool hic_normalizations_seen = false;
    std::string threads = "1";
};

[[noreturn]] void argument_error(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicConvertFormat: error: %s\n", message.c_str());
    std::exit(2);
}

// log.error followed by exit(1), which is how the tool reports a bad
// combination of otherwise valid arguments.
int reject(const std::string& message) {
    std::fprintf(stderr, "ERROR:hicexplorer.hicConvertFormat:%s\n", message.c_str());
    return 1;
}

void report(const std::string& level, const std::string& message) {
    std::fprintf(stderr, "%s:hicexplorer.hicConvertFormat:%s\n", level.c_str(),
                 message.c_str());
}

bool one_of(const std::string& value, const std::vector<std::string>& choices) {
    return std::find(choices.begin(), choices.end(), value) != choices.end();
}

std::string choice_error(const std::string& option, const std::string& value,
                         const std::vector<std::string>& choices) {
    std::string text = "argument " + option + ": invalid choice: '" + value +
                       "' (choose from ";
    for (std::size_t i = 0; i < choices.size(); ++i) {
        text += "'" + choices[i] + "'";
        if (i + 1 < choices.size()) {
            text += ", ";
        }
    }
    return text + ")";
}

const std::vector<std::string> kInputFormats{"h5",    "cool",   "hic",
                                             "homer", "hicpro", "2D-text"};
const std::vector<std::string> kOutputFormats{"cool",          "h5",    "homer",
                                              "ginteractions", "mcool", "hicpro",
                                              "hic"};
const std::vector<std::string> kHicVersions{"8", "9"};
const std::vector<std::string> kHicNormalizations{"VC", "VC_SQRT", "KR", "SCALE", "none"};

Arguments parse_arguments(int argc, char** argv) {
    Arguments args;
    bool matrices_seen = false;
    bool out_files_seen = false;
    bool input_format_seen = false;
    bool output_format_seen = false;

    // The option that is currently collecting values, for nargs='+'.
    std::vector<std::string>* collecting = nullptr;
    // The option that still needs exactly one value.
    std::string* pending = nullptr;
    std::string pending_name;

    for (int i = 1; i < argc; ++i) {
        std::string token(argv[i]);
        const bool looks_like_option =
            token.size() > 1 && token[0] == '-' &&
            std::isdigit(static_cast<unsigned char>(token[1])) == 0;

        if (pending != nullptr && !looks_like_option) {
            *pending = token;
            pending = nullptr;
            continue;
        }
        if (!looks_like_option) {
            if (collecting != nullptr) {
                collecting->push_back(token);
                continue;
            }
            argument_error("unrecognized arguments: " + token);
        }
        if (pending != nullptr) {
            argument_error("argument " + pending_name + ": expected one argument");
        }
        collecting = nullptr;

        std::string name = token;
        std::optional<std::string> inline_value;
        if (const std::size_t equals = token.find('=');
            equals != std::string::npos && token.rfind("--", 0) == 0) {
            name = token.substr(0, equals);
            inline_value = token.substr(equals + 1);
        }

        const auto take_list = [&](std::vector<std::string>& target, bool& seen) {
            seen = true;
            if (inline_value.has_value()) {
                target.push_back(*inline_value);
            } else {
                collecting = &target;
            }
        };
        const auto take_value = [&](std::string& target, const std::string& option) {
            if (inline_value.has_value()) {
                target = *inline_value;
            } else {
                pending = &target;
                pending_name = option;
            }
        };

        if (name == "-h" || name == "--help") {
            std::fputs(kUsage, stdout);
            std::fputs(kHelp, stdout);
            std::exit(0);
        }
        if (name == "--version") {
            std::printf("hicConvertFormat %s\n", hicx::kVersion);
            std::exit(0);
        }
        if (name == "-m" || name == "--matrices") {
            take_list(args.matrices, matrices_seen);
            continue;
        }
        if (name == "-o" || name == "--outFileName") {
            take_list(args.out_file_names, out_files_seen);
            continue;
        }
        if (name == "-r" || name == "--resolutions") {
            bool seen = false;
            take_list(args.resolutions, seen);
            continue;
        }
        if (name == "-bf" || name == "--bedFileHicpro") {
            bool seen = false;
            take_list(args.bed_file_hicpro, seen);
            continue;
        }
        if (name == "--inputFormat") {
            input_format_seen = true;
            take_value(args.input_format, "--inputFormat");
            continue;
        }
        if (name == "--outputFormat") {
            output_format_seen = true;
            take_value(args.output_format, "--outputFormat");
            continue;
        }
        if (name == "--correction_name") {
            take_value(args.correction_name, "--correction_name");
            continue;
        }
        if (name == "--chromosome") {
            take_value(args.chromosome, "--chromosome");
            continue;
        }
        if (name == "-cs" || name == "--chromosomeSizes") {
            take_value(args.chromosome_sizes, "--chromosomeSizes");
            continue;
        }
        if (name == "--hicVersion") {
            take_value(args.hic_version, "--hicVersion");
            continue;
        }
        if (name == "--hicNormalizations") {
            if (!args.hic_normalizations_seen) {
                args.hic_normalizations.clear();
                args.hic_normalizations_seen = true;
            }
            bool seen = false;
            take_list(args.hic_normalizations, seen);
            continue;
        }
        if (name == "--threads") {
            take_value(args.threads, "--threads");
            continue;
        }
        if (name == "--correction_division") {
            args.correction_division = true;
            continue;
        }
        if (name == "--store_applied_correction") {
            args.store_applied_correction = true;
            continue;
        }
        if (name == "--enforce_integer") {
            args.enforce_integer = true;
            continue;
        }
        if (name == "--load_raw_values") {
            args.load_raw_values = true;
            continue;
        }
        argument_error("unrecognized arguments: " + token);
    }
    if (pending != nullptr) {
        argument_error("argument " + pending_name + ": expected one argument");
    }

    std::vector<std::string> missing;
    if (!matrices_seen) {
        missing.push_back("--matrices/-m");
    }
    if (!out_files_seen) {
        missing.push_back("--outFileName/-o");
    }
    if (!input_format_seen) {
        missing.push_back("--inputFormat");
    }
    if (!output_format_seen) {
        missing.push_back("--outputFormat");
    }
    if (!missing.empty()) {
        std::string text = "the following arguments are required: ";
        for (std::size_t i = 0; i < missing.size(); ++i) {
            text += missing[i];
            if (i + 1 < missing.size()) {
                text += ", ";
            }
        }
        argument_error(text);
    }
    if (args.matrices.empty()) {
        argument_error("argument --matrices/-m: expected at least one argument");
    }
    if (args.out_file_names.empty()) {
        argument_error("argument --outFileName/-o: expected at least one argument");
    }
    if (!one_of(args.input_format, kInputFormats)) {
        argument_error(choice_error("--inputFormat", args.input_format, kInputFormats));
    }
    if (!one_of(args.output_format, kOutputFormats)) {
        argument_error(
            choice_error("--outputFormat", args.output_format, kOutputFormats));
    }
    if (args.hic_version == "6" || args.hic_version == "7") {
        // .hic input of versions 6 to 9 is read; writing 6 or 7 is refused
        // because no Juicer tools release that writes them can be obtained to
        // validate against (hicfilecpp docs/PROVENANCE.md).
        argument_error("argument --hicVersion: writing .hic version " + args.hic_version +
                       " is not supported: no Juicer tools release that writes it can be "
                       "obtained to validate against (choose from '8', '9')");
    }
    if (!one_of(args.hic_version, kHicVersions)) {
        argument_error(choice_error("--hicVersion", args.hic_version, kHicVersions));
    }
    if (args.hic_normalizations_seen && args.hic_normalizations.empty()) {
        argument_error("argument --hicNormalizations: expected at least one argument");
    }
    for (const auto& norm : args.hic_normalizations) {
        if (!one_of(norm, kHicNormalizations)) {
            argument_error(choice_error("--hicNormalizations", norm, kHicNormalizations));
        }
    }
    {
        const bool digits = !args.threads.empty() &&
                            std::all_of(args.threads.begin(), args.threads.end(),
                                        [](unsigned char c) { return std::isdigit(c) != 0; });
        if (!digits || std::stoll(args.threads) < 1 || std::stoll(args.threads) > 1024) {
            argument_error("argument --threads: invalid value: '" + args.threads + "'");
        }
    }
    if (!args.chromosome_sizes.empty()) {
        // argparse.FileType('r') opens the file while parsing.
        std::FILE* handle = std::fopen(args.chromosome_sizes.c_str(), "r");
        if (handle == nullptr) {
            argument_error("argument --chromosomeSizes/-cs: can't open '" +
                           args.chromosome_sizes + "'");
        }
        std::fclose(handle);
    }
    return args;
}

// What a load produced, plus the two pieces of cool provenance that the output
// handler is given (hicConvertFormat.py:263-268).
struct Loaded {
    hicx::MatrixData data;
    std::optional<std::string> hic2cool_version;
    std::map<std::string, std::string> metadata;
    bool has_metadata = false;
};

Loaded load_input(const Arguments& args, std::size_t index) {
    Loaded loaded;
    const std::string& matrix = args.matrices[index];
    if (args.input_format == "hicpro") {
        loaded.data = hicx::read_hicpro(matrix, args.bed_file_hicpro[index]);
    } else if (args.input_format == "homer") {
        loaded.data = hicx::read_homer(matrix);
    } else if (args.input_format == "2D-text") {
        const std::vector<std::pair<std::string, std::int64_t>> sizes =
            hicx::read_chromosome_sizes(args.chromosome_sizes);
        const std::int64_t resolution = std::stoll(args.resolutions[0]);
        loaded.data = hicx::read_two_dimensional_text(matrix, sizes, resolution);
    } else if (args.input_format == "h5") {
        loaded.data = hicx::read_hicexplorer_h5(matrix);
    } else {
        hicx::CoolLoadOptions options;
        options.correction_factor_table = args.correction_name;
        options.apply_correction = !args.load_raw_values;
        if (args.correction_division) {
            options.correction_operator = '/';
        }
        hicx::CoolLoadResult result = hicx::read_cool(matrix, options);
        loaded.data = std::move(result.data);
        loaded.hic2cool_version = result.hic2cool_version;
        loaded.metadata = std::move(result.metadata);
        loaded.has_metadata = true;
    }
    return loaded;
}

// The save options the output Cool object carries. Note what is *not* set: the
// output handler is never given a correction operator (hicConvertFormat.py:268
// passes no pCorrectionOperator), so the inversion test at cool.py:307 sees
// None and only fileWasH5 or a hic2cool version can trigger it.
hicx::CoolSaveOptions cool_options(const Arguments& args, const Loaded& loaded,
                                   bool format_was_h5, bool enforce_integer,
                                   bool append,
                                   const std::map<std::string, std::string>& metadata,
                                   bool has_metadata) {
    hicx::CoolSaveOptions options;
    options.symmetric = true;
    options.apply_correction = !args.store_applied_correction;
    options.enforce_integer = enforce_integer;
    options.file_was_h5 = format_was_h5;
    options.hic2cool_version = loaded.hic2cool_version;
    options.hic_metadata = metadata;
    options.has_hic_metadata = has_metadata;
    options.append = append;
    return options;
}

// create_cooler_input pops these three out of the dictionary it was handed,
// and the dictionary is shared between the resolutions of one mcool, so the
// second resolution no longer sees them.
void consume_metadata(std::map<std::string, std::string>& metadata) {
    metadata.erase("matrix-generated-by");
    metadata.erase("matrix-generated-by-url");
    metadata.erase("genome-assembly");
}

int run(const Arguments& args);

// Removes a file when it goes out of scope.
struct TemporaryFile {
    std::string path;
    ~TemporaryFile() {
        std::error_code ignored;
        std::filesystem::remove(path, ignored);
    }
};

std::vector<std::string> split_on_dots(const std::string& text) {
    std::vector<std::string> parts{""};
    for (const char c : text) {
        if (c == '.') {
            parts.emplace_back();
        } else {
            parts.back() += c;
        }
    }
    return parts;
}

// --inputFormat hic.
int run_hic_input(const Arguments& args) {
    if (args.output_format == "hic") {
        return reject("hic to hic conversion is not supported.");
    }
    std::vector<std::int64_t> resolutions;
    for (const auto& text : args.resolutions) {
        resolutions.push_back(std::stoll(text));
    }
    try {
        if (args.output_format == "cool") {
            // hicConvertFormat.py:124-138. The Python indexes outFileName by
            // the matrix and splits it on dots, so a missing name or a name
            // without an extension raises IndexError there.
            if (args.out_file_names.size() < args.matrices.size()) {
                return reject("Number of input matrices is larger than the number of output "
                              "file names.");
            }
            for (std::size_t i = 0; i < args.matrices.size(); ++i) {
                if (resolutions.empty()) {
                    hicx::hic2cool_convert(args.matrices[i], args.out_file_names[i], 0);
                    continue;
                }
                for (std::size_t j = 0; j < resolutions.size(); ++j) {
                    std::vector<std::string> parts = split_on_dots(args.out_file_names[i]);
                    if (parts.size() < 2) {
                        return reject("The output file name " + args.out_file_names[i] +
                                      " has no extension to insert the resolution before.");
                    }
                    parts[parts.size() - 2] += "_" + args.resolutions[j];
                    std::string name = parts[0];
                    for (std::size_t k = 1; k < parts.size(); ++k) {
                        name += "." + parts[k];
                    }
                    hicx::hic2cool_convert(args.matrices[i], name, resolutions[j]);
                }
            }
            return 0;
        }
        if (args.matrices.size() != args.out_file_names.size()) {
            return reject("Number of input matrices does not match number output "
                          "matrices!: Input matrices " +
                          std::to_string(args.matrices.size()) + "; output matrices " +
                          std::to_string(args.out_file_names.size()));
        }
        if (args.output_format == "mcool") {
            for (std::size_t i = 0; i < args.matrices.size(); ++i) {
                hicx::hic2cool_convert_mcool(args.matrices[i], args.out_file_names[i],
                                             resolutions);
            }
            return 0;
        }
        if (resolutions.size() != 1) {
            return reject("--inputFormat hic with --outputFormat " + args.output_format +
                          " needs exactly one resolution in --resolutions.");
        }
        for (std::size_t i = 0; i < args.matrices.size(); ++i) {
            const std::filesystem::path out(args.out_file_names[i]);
            TemporaryFile temporary{
                (out.parent_path() / ("." + out.filename().string() + ".hic2cool.cool")).string()};
            hicx::hic2cool_convert(args.matrices[i], temporary.path, resolutions[0]);
            Arguments cool = args;
            cool.input_format = "cool";
            cool.matrices = {temporary.path};
            cool.out_file_names = {args.out_file_names[i]};
            cool.resolutions.clear();
            if (args.output_format == "hicpro" && i < args.bed_file_hicpro.size()) {
                cool.bed_file_hicpro = {args.bed_file_hicpro[i]};
            }
            const int status = run(cool);
            if (status != 0) {
                return status;
            }
        }
        return 0;
    } catch (const hicx::Hic2coolExit& exit) {
        // hic2cool's force_exit: the message on stderr, then sys.exit(1).
        std::fprintf(stderr, "%s\n", exit.what());
        return 1;
    }
}

// The resolution groups of a multi resolution cool file, numerically sorted;
// empty for anything else.
std::vector<std::string> mcool_resolutions(const std::string& path) {
    if (path.find("::") != std::string::npos || !hicx::h5::is_hdf5(path)) {
        return {};
    }
    const hicx::h5::File file(path);
    if (!file.exists("/resolutions")) {
        return {};
    }
    std::vector<std::string> groups = file.children("/resolutions");
    std::sort(groups.begin(), groups.end(), [](const std::string& a, const std::string& b) {
        return std::stoll(a) < std::stoll(b);
    });
    return groups;
}

// --outputFormat hic.
int run_hic_output(const Arguments& args) {
    hicx::HicWriteOptions options;
    options.version = std::stoi(args.hic_version);
    options.threads = static_cast<int>(std::stoll(args.threads));
    options.normalizations.clear();
    if (!one_of("none", args.hic_normalizations)) {
        for (const auto& norm : {"VC", "VC_SQRT", "KR", "SCALE"}) {
            if (one_of(norm, args.hic_normalizations)) {
                options.normalizations.emplace_back(norm);
            }
        }
    }
    std::vector<std::int64_t> resolutions;
    for (const auto& text : args.resolutions) {
        resolutions.push_back(std::stoll(text));
    }
    for (std::size_t i = 0; i < args.matrices.size(); ++i) {
        std::vector<Loaded> loaded;
        std::vector<std::int64_t> extra = resolutions;
        const std::vector<std::string> groups =
            args.input_format == "cool" ? mcool_resolutions(args.matrices[i])
                                        : std::vector<std::string>{};
        if (!groups.empty()) {
            for (const auto& group : groups) {
                if (!resolutions.empty() &&
                    std::find(resolutions.begin(), resolutions.end(), std::stoll(group)) ==
                        resolutions.end()) {
                    continue;
                }
                Arguments single = args;
                single.matrices = {args.matrices[i] + "::/resolutions/" + group};
                loaded.push_back(load_input(single, 0));
            }
            if (loaded.empty()) {
                return reject("None of --resolutions is a resolution of " + args.matrices[i]);
            }
            extra.clear();
        } else {
            loaded.push_back(load_input(args, i));
        }
        const auto assembly = loaded.front().metadata.find("genome-assembly");
        if (assembly != loaded.front().metadata.end() && !assembly->second.empty()) {
            options.genome = assembly->second;
        }
        std::vector<const hicx::MatrixData*> matrices;
        for (const auto& entry : loaded) {
            matrices.push_back(&entry.data);
        }
        hicx::write_hic(args.out_file_names[i], matrices, extra, options);
    }
    return 0;
}

int run(const Arguments& args) {
    if (args.input_format != "hic" && args.output_format != "mcool") {
        if (args.matrices.size() != args.out_file_names.size()) {
            return reject("Number of input matrices does not match number output "
                          "matrices!: Input matrices " +
                          std::to_string(args.matrices.size()) + "; output matrices " +
                          std::to_string(args.out_file_names.size()));
        }
    }
    if (args.input_format == "hic") {
        return run_hic_input(args);
    }
    if (!args.chromosome.empty()) {
        return reject("--chromosome is not supported yet: loading a single "
                      "chromosome out of a cooler is a distinct cooler code path "
                      "that the C++ cool reader does not have.");
    }

    const bool format_was_h5 = args.input_format == "h5";
    const bool apply_correction = !args.store_applied_correction;

    if (args.input_format == "hicpro") {
        if (args.matrices.size() != args.bed_file_hicpro.size()) {
            reject("Number of matrices and associated bed files need to be the same.");
            return reject("Matrices: " + std::to_string(args.matrices.size()) +
                          "; Bed files: " + std::to_string(args.bed_file_hicpro.size()));
        }
    }
    if (args.input_format == "2D-text") {
        if (args.resolutions.empty()) {
            return reject("The resolution must be defined via --resolutions");
        }
        if (args.chromosome_sizes.empty()) {
            return reject("The sizes of the chromosomes must be defined via "
                          "--chromosomeSizes.");
        }
    }

    if (args.output_format == "hic") {
        return run_hic_output(args);
    }

    for (std::size_t i = 0; i < args.matrices.size(); ++i) {
        Loaded loaded = load_input(args, i);

        if (args.output_format == "cool" || args.output_format == "h5" ||
            args.output_format == "homer" || args.output_format == "ginteractions") {
            if (args.output_format == "homer") {
                // triu then maximum with the transpose, not a sum: a sum would
                // double the diagonal.
                loaded.data.matrix = hicx::maximum_with_transpose(loaded.data.matrix);
            } else if (args.output_format == "ginteractions") {
                // Ginteractions.save takes the upper triangle of that same
                // symmetric matrix straight away, so the mirror is never
                // needed. Only the triangle is built.
                loaded.data.matrix =
                    hicx::upper_triangle_after_maximum(loaded.data.matrix);
            }
            const std::string& out = args.out_file_names[i];
            if (args.output_format == "cool") {
                std::map<std::string, std::string> metadata = loaded.metadata;
                hicx::CoolSaveOptions options =
                    cool_options(args, loaded, format_was_h5, args.enforce_integer,
                                 false, metadata, loaded.has_metadata);
                hicx::write_cool(out, loaded.data, options);
            } else if (args.output_format == "h5") {
                hicx::H5SaveOptions options;
                options.symmetric = true;
                hicx::write_hicexplorer_h5(out, loaded.data, options);
            } else if (args.output_format == "homer") {
                hicx::write_homer(out, loaded.data);
            } else {
                hicx::write_ginteractions(out, loaded.data);
            }
        }

        if (args.output_format == "hicpro") {
            if (args.matrices.size() == args.out_file_names.size() &&
                args.out_file_names.size() == args.bed_file_hicpro.size()) {
                hicx::write_hicpro(args.out_file_names[i], args.bed_file_hicpro[i],
                                   loaded.data);
            } else {
                // The Python builds this message from args.matrix, an attribute
                // that does not exist, so it raises AttributeError instead of
                // exiting cleanly. The port reports the mismatch and exits 1.
                return reject("The number of input matrices, output files and bed "
                              "files does not match: Input: " +
                              std::to_string(args.matrices.size()) + "; Output: " +
                              std::to_string(args.out_file_names.size()) + "; Bed: " +
                              std::to_string(args.bed_file_hicpro.size()));
            }
        } else if (args.output_format == "mcool") {
            if (!args.resolutions.empty() && args.matrices.size() > 1) {
                // Logged and not acted on, as in the Python: the run continues
                // and every input matrix rewrites the same output file.
                report("ERROR", "Please define one matrix and many resolutions which "
                                "should be created or multiple matrices.");
            }
            if (!args.resolutions.empty()) {
                report("INFO", "Correction factors are removed. They are not valid "
                               "for any new created resolution.");
                // hiCMatrix.setMatrix keeps only the matrix and the bin table,
                // so the nan bins, correction factors and distance counts of
                // the input do not reach the output.
                hicx::MatrixData base;
                base.matrix = loaded.data.matrix;
                base.cut_intervals = loaded.data.cut_intervals;
                const std::int64_t bin_size =
                    hicx::BinTable(base.cut_intervals).bin_size();
                std::map<std::string, std::string> metadata = loaded.metadata;

                for (std::size_t j = 0; j < args.resolutions.size(); ++j) {
                    const std::int64_t resolution = std::stoll(args.resolutions[j]);
                    const std::int64_t merge_factor =
                        bin_size == 0 ? 0 : resolution / bin_size;
                    hicx::MatrixData merged =
                        resolution != bin_size ? hicx::merge_bins(base, merge_factor)
                                               : base;
                    hicx::CoolSaveOptions options = cool_options(
                        args, loaded, format_was_h5, args.enforce_integer, j > 0,
                        metadata, loaded.has_metadata);
                    hicx::write_cool(args.out_file_names[0] + "::/resolutions/" +
                                         args.resolutions[j],
                                     merged, options);
                    consume_metadata(metadata);
                }
            } else {
                const std::int64_t bin_size =
                    hicx::BinTable(loaded.data.cut_intervals).bin_size();
                // This branch builds no output handler metadata at all: it
                // passes neither pHiCInfo nor pHic2CoolVersion nor
                // pEnforceInteger (hicConvertFormat.py:338-339).
                hicx::CoolSaveOptions options;
                options.symmetric = true;
                options.apply_correction = apply_correction;
                options.file_was_h5 = format_was_h5;
                options.append = i > 0;
                hicx::write_cool(args.out_file_names[0] + "::/resolutions/" +
                                     std::to_string(bin_size),
                                 loaded.data, options);
            }
        }
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    int status = 0;
    try {
        status = run(args);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicConvertFormat: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicConvertFormat");
    return status;
}
