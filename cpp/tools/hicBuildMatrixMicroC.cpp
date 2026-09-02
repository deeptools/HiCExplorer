// Port of hicexplorer/hicBuildMatrixMicroC.py.
//
// Micro-C libraries are digested with MNase rather than a restriction enzyme,
// so there are no cut sites, no restriction sequence and no dangling end. The
// Python file is hicBuildMatrix.py with those arguments deleted, --binSize
// promoted to required, and a main() that calls the same
// lib/buildMatrixMethods.createMatrix with pRestrictionCutFile,
// pRestrictionSequence, pDanglingSequence, pMinDistance, pMaxDistance and
// pKeepSelfLigation all None. This port is the same delta: the argument parser
// below, and build_matrix_impl.hpp for everything else.
//
// What the missing arguments change inside createMatrix:
//
//   * rf_positions is None, so the self circle branch and the restriction site
//     lookup of the inward branch are never entered. Every close inward pair
//     that is not a dangling end becomes "same fragment", and self circles are
//     never counted at all.
//   * pMinDistance is None, so the QC log header omits the
//     "Min rest. site distance" line.
//   * pRestrictionCutFile is None, so the QC log omits "self ligation",
//     "One mate not close to rest site" and "self circle".
//   * the dangling end dictionary is empty, so no dangling end line is
//     written either.
//   * --keepSelfCircles is still accepted and, exactly as in hicBuildMatrix,
//     still does nothing.
//
// Bins can only come from --binSize here, which is why the Python makes it
// required: get_rf_bins on an empty cut site list would fail on the unpacking.

#include "build_matrix_impl.hpp"

namespace {

const char* const kUsage =
    "usage: hicBuildMatrixMicroC --samFiles two sam files two sam files\n"
    "                            --outFileName FILENAME --QCfolder FOLDER\n"
    "                            [--outBam bam file] --binSize BINSIZE\n"
    "                            [BINSIZE ...]\n"
    "                            [--maxLibraryInsertSize MAXLIBRARYINSERTSIZE]\n"
    "                            [--genomeAssembly GENOMEASSEMBLY]\n"
    "                            [--region CHR:START-END] [--keepSelfCircles]\n"
    "                            [--minMappingQuality MINMAPPINGQUALITY]\n"
    "                            [--threads THREADS]\n"
    "                            [--inputBufferSize INPUTBUFFERSIZE] [--doTestRun]\n"
    "                            [--doTestRunLines DOTESTRUNLINES]\n"
    "                            [--skipDuplicationCheck]\n"
    "                            [--chromosomeSizes txt file] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Using an alignment from a program that supports local alignment (eg. Bowtie2)\n"
    "where both PE reads are mapped using the --local option, this program reads such\n"
    "file and creates a matrix of interactions.\n"
    "\n"
    "Required arguments:\n"
    "  --samFiles two sam files two sam files, -s two sam files two sam files\n"
    "                        The two PE alignment sam files to process\n"
    "  --outFileName FILENAME, -o FILENAME\n"
    "                        Output file name for the Hi-C matrix.\n"
    "  --QCfolder FOLDER     Path of folder to save the quality control data for the\n"
    "                        matrix.\n"
    "  --binSize BINSIZE [BINSIZE ...], -bs BINSIZE [BINSIZE ...]\n"
    "                        Size in bp for the bins. (default: None)\n"
    "\n"
    "Optional arguments:\n"
    "  --outBam bam file, -b bam file\n"
    "                        Output bam file with all valid Hi-C reads.\n"
    "  --maxLibraryInsertSize MAXLIBRARYINSERTSIZE\n"
    "                        The maximum library insert size. (default: 1000)\n"
    "  --genomeAssembly GENOMEASSEMBLY, -ga GENOMEASSEMBLY\n"
    "                        The genome the reads were mapped to.\n"
    "  --region CHR:START-END, -r CHR:START-END\n"
    "                        Region of the genome to limit the operation to.\n"
    "  --keepSelfCircles     Keep self circles. (default: False)\n"
    "  --minMappingQuality MINMAPPINGQUALITY\n"
    "                        Minimum mapping quality. (default: 15)\n"
    "  --threads THREADS     Number of threads. (default: 4)\n"
    "  --inputBufferSize INPUTBUFFERSIZE\n"
    "                        Size of the input buffer of each thread.\n"
    "                        (default: 400000)\n"
    "  --doTestRun           Test only --doTestRunLines reads. (default: False)\n"
    "  --doTestRunLines DOTESTRUNLINES\n"
    "                        Number of lines for the qc test run. (default: 1000000)\n"
    "  --skipDuplicationCheck\n"
    "                        Skip the identification of duplicated read pairs.\n"
    "  --chromosomeSizes txt file, -cs txt file\n"
    "                        File with the chromosome sizes for your genome.\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

Arguments parse_micro_c_arguments(int argc, char** argv) {
    Arguments args;
    // pMinDistance is None in hicBuildMatrixMicroC.main, which is falsy, so
    // the QC log takes the branch without the "Min rest. site distance" line.
    args.min_distance = 0;
    bool sam_seen = false;
    bool out_seen = false;
    bool qc_seen = false;
    bool bin_size_seen = false;

    for (int i = 1; i < argc; ++i) {
        std::string name(argv[i]);
        std::optional<std::string> inline_value;
        const std::size_t equals = name.find('=');
        if (equals != std::string::npos && name.rfind("--", 0) == 0) {
            inline_value = name.substr(equals + 1);
            name = name.substr(0, equals);
        }
        auto take_one = [&](const std::string& option) -> std::string {
            if (inline_value.has_value()) {
                return *inline_value;
            }
            if (i + 1 >= argc) {
                fail("argument " + option + ": expected one argument");
            }
            return std::string(argv[++i]);
        };
        auto take_many = [&](const std::string& option,
                             int minimum) -> std::vector<std::string> {
            std::vector<std::string> values;
            if (inline_value.has_value()) {
                values.push_back(*inline_value);
            }
            while (i + 1 < argc && !is_option(argv[i + 1])) {
                values.emplace_back(argv[++i]);
            }
            if (static_cast<int>(values.size()) < minimum) {
                fail("argument " + option + ": expected at least " +
                     std::to_string(minimum) + " arguments");
            }
            return values;
        };

        if (name == "-h" || name == "--help") {
            std::fputs(kUsage, stdout);
            std::fputs(kHelp, stdout);
            std::exit(0);
        }
        if (name == "--version") {
            std::printf("hicBuildMatrixMicroC %s\n", hicx::kVersion);
            std::exit(0);
        }
        if (name == "-s" || name == "--samFiles") {
            args.sam_files = take_many("--samFiles/-s", 2);
            if (args.sam_files.size() != 2) {
                fail("argument --samFiles/-s: expected 2 arguments");
            }
            sam_seen = true;
        } else if (name == "-o" || name == "--outFileName") {
            args.out_file_name = take_one("--outFileName/-o");
            out_seen = true;
        } else if (name == "--QCfolder") {
            args.qc_folder = take_one("--QCfolder");
            qc_seen = true;
        } else if (name == "-b" || name == "--outBam") {
            args.out_bam = take_one("--outBam/-b");
        } else if (name == "-bs" || name == "--binSize") {
            for (const auto& value : take_many("--binSize/-bs", 1)) {
                args.bin_size.push_back(parse_int("--binSize/-bs", value));
            }
            bin_size_seen = true;
        } else if (name == "--maxLibraryInsertSize") {
            args.max_library_insert_size = parse_int(name, take_one(name));
        } else if (name == "-ga" || name == "--genomeAssembly") {
            args.genome_assembly = take_one("--genomeAssembly/-ga");
        } else if (name == "-r" || name == "--region") {
            args.region = hicx::normalise_region(take_one("--region/-r"));
        } else if (name == "--keepSelfCircles") {
            args.keep_self_circles = true;
        } else if (name == "--minMappingQuality") {
            args.min_mapping_quality = parse_int(name, take_one(name));
        } else if (name == "--threads") {
            args.threads = parse_int(name, take_one(name));
        } else if (name == "--inputBufferSize") {
            args.input_buffer_size = parse_int(name, take_one(name));
        } else if (name == "--doTestRun") {
            args.do_test_run = true;
        } else if (name == "--doTestRunLines") {
            args.do_test_run_lines = parse_int(name, take_one(name));
        } else if (name == "--skipDuplicationCheck") {
            args.skip_duplication_check = true;
        } else if (name == "-cs" || name == "--chromosomeSizes") {
            args.chromosome_sizes = take_one("--chromosomeSizes/-cs");
        } else {
            fail("unrecognized arguments: " + name);
        }
    }

    std::string missing;
    auto require = [&missing](bool seen, const char* text) {
        if (!seen) {
            missing += missing.empty() ? text : (std::string(", ") + text);
        }
    };
    require(sam_seen, "--samFiles/-s");
    require(out_seen, "--outFileName/-o");
    require(qc_seen, "--QCfolder");
    require(bin_size_seen, "--binSize/-bs");
    if (!missing.empty()) {
        fail("the following arguments are required: " + missing);
    }
    return args;
}

}  // namespace

int main(int argc, char** argv) {
    g_tool = "hicBuildMatrixMicroC";
    g_usage = kUsage;
    return run_build_matrix(parse_micro_c_arguments(argc, argv));
}
