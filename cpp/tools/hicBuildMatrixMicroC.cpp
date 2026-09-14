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
#include "hicx/argparse.hpp"

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

// hicBuildMatrixMicroC.py parse_arguments.
Arguments parse_micro_c_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicBuildMatrixMicroC",
                       "Using an alignment from a program that supports local alignment (eg. "
                       "Bowtie2) where both PE reads are mapped using the --local option, this "
                       "program reads such file and creates a matrix of interactions.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--samFiles", "-s"})
        .metavar("two sam files")
        .nargs(2)
        .file_type("r")
        .required()
        .input({"bam", "sam"})
        .help("The two PE alignment sam files to process");
    required.add({"--outFileName", "-o"})
        .metavar("FILENAME")
        .file_type("w")
        .required()
        .output({"h5", "cool", "mcool"})
        .help("Output file name for the Hi-C matrix.");
    required.add({"--QCfolder"})
        .metavar("FOLDER")
        .required()
        .output({}, "directory")
        .help("Path of folder to save the quality control data for the matrix.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--outBam", "-b"})
        .metavar("bam file")
        .file_type("w")
        .output({"bam"})
        .help("Output bam file with all valid Hi-C reads.");
    optional.add({"--binSize", "-bs"})
        .type("int")
        .nargs("+")
        .required()
        .help("Size in bp for the bins.");
    optional.add({"--maxLibraryInsertSize"})
        .type("int")
        .default_value(1000)
        .help("The maximum library insert size.");
    optional.add({"--genomeAssembly", "-ga"}).help("The genome the reads were mapped to.");
    optional.add({"--region", "-r"})
        .metavar("CHR:START-END")
        .type("genomicRegion")
        .check(genomic_region_check)
        .help("Region of the genome to limit the operation to.");
    optional.add({"--keepSelfCircles"}).action(cli::Action::StoreTrue).help("Keep self circles.");
    optional.add({"--minMappingQuality"})
        .default_value(15)
        .type("int")
        .help("Minimum mapping quality.");
    optional.add({"--threads"}).default_value(4).type("int").help("Number of threads.");
    optional.add({"--inputBufferSize"})
        .default_value(400000)
        .type("int")
        .help("Size of the input buffer of each thread.");
    optional.add({"--doTestRun"})
        .action(cli::Action::StoreTrue)
        .help("Test only --doTestRunLines reads.");
    optional.add({"--doTestRunLines"})
        .default_value(1000000)
        .type("int")
        .help("Number of lines for the qc test run.");
    optional.add({"--skipDuplicationCheck"})
        .action(cli::Action::StoreTrue)
        .help("Skip the identification of duplicated read pairs.");
    optional.add({"--chromosomeSizes", "-cs"})
        .file_type("r")
        .metavar("txt file")
        .input({"txt"})
        .help("File with the chromosome sizes for your genome.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    // pMinDistance is None in hicBuildMatrixMicroC.main, which is falsy, so
    // the QC log takes the branch without the "Min rest. site distance" line.
    args.min_distance = 0;
    args.sam_files = ns.strs("samFiles");
    args.out_file_name = ns.str("outFileName");
    args.qc_folder = ns.str("QCfolder");
    args.out_bam = ns.opt_str("outBam").value_or("");
    args.bin_size = ns.integers("binSize");
    args.max_library_insert_size = ns.integer("maxLibraryInsertSize");
    args.genome_assembly = ns.opt_str("genomeAssembly").value_or("");
    args.region = hicx::normalise_region(ns.opt_str("region").value_or(""));
    args.keep_self_circles = ns.flag("keepSelfCircles");
    args.min_mapping_quality = ns.integer("minMappingQuality");
    args.threads = ns.integer("threads");
    args.input_buffer_size = ns.integer("inputBufferSize");
    args.do_test_run = ns.flag("doTestRun");
    args.do_test_run_lines = ns.integer("doTestRunLines");
    args.skip_duplication_check = ns.flag("skipDuplicationCheck");
    args.chromosome_sizes = ns.opt_str("chromosomeSizes").value_or("");
    return args;
}

}  // namespace

int main(int argc, char** argv) {
    g_tool = "hicBuildMatrixMicroC";
    const Arguments args = parse_micro_c_arguments(argc, argv);
    // The QC report is drawn after every successful run.
    if (const int refused = hicx::plot::preflight("hicPrepareQCreport", true); refused != 0) {
        return refused;
    }
    const int status = run_build_matrix(args);
    return status != 0 ? status : draw_qc_report(args.qc_folder);
}
