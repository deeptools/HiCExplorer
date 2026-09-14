// Port of hicexplorer/hicBuildMatrix.py and hicexplorer/lib/buildMatrixMethods.py.
//
// The entry point of every Hi-C workflow: two name ordered BAM files holding
// the two ends of each ligation product go in, a contact matrix and a QC
// folder come out.
//
// ---------------------------------------------------------------------------
// Structure
// ---------------------------------------------------------------------------
// The Python spawns --threads minus one multiprocessing workers, each of which
// receives a fork of the parent holding a share of the read buffer, and
// consumes their results through a queue in completion order. The port keeps
// the shape of that pipeline and drops the process boundary
// (cpp/PLAN.md 4.4 rule 8):
//
//   * one reader, serial, because the PCR duplicate check is stateful and
//     because htslib decompresses one BGZF stream in order anyway;
//   * a chunk of up to --inputBufferSize accepted pairs is then split into
//     `threads` *contiguous* index ranges, each classified by one thread into
//     its own counters, its own pixel list and its own out-BAM index list;
//   * the ranges are combined in index order, so the result is exactly the
//     sequential result for any thread count (cpp/OPTIMIZATION.md 3).
//
// The only shared mutable state is the coverage vector, and it is touched with
// relaxed integer fetch_add: a sum of ones is the same number whatever order
// it is accumulated in, so this is deterministic by construction rather than
// by luck. The Python's equivalent, a multiprocessing Array incremented with a
// non atomic read-modify-write from several processes, is not.
//
// ---------------------------------------------------------------------------
// Reproduced defects (never silently fixed)
// ---------------------------------------------------------------------------
//  Q1 buildMatrixMethods.py:523. `if mate1.mapq == 0 & mate2.mapq == 0` parses
//     as `mate1.mapq == (0 & mate2.mapq) == 0`, a chained comparison that is
//     simply `mate1.mapq == 0`. So a pair is charged to "One mate not unique"
//     on the first mate's quality alone and the second mate is never looked
//     at. Reproduced.
//  Q2 buildMatrixMethods.py:759-762. The self circle test's `continue` sits
//     inside the `for restrictionSequence` loop, not inside the read loop, so
//     a self circle is *counted* and then kept regardless of
//     --keepSelfCircles. The option changes the QC log and nothing else.
//     With two restriction sequences the same pair can be counted twice.
//     Reproduced.
//  Q3 buildMatrixMethods.py:841-851. The coverage vector is filled using
//     `mate_bin`, which after the bin lookup loop holds the *second* mate's
//     bin for both mates. So the first mate's offset is measured against the
//     wrong bin start and its coverage lands in the second mate's bin.
//     Reproduced.
//  Q4 buildMatrixMethods.py:1231. The per bin maximum is taken over
//     range(begin, end) with end already the last valid index, so the last
//     coverage cell of every bin is never read. Reproduced.
//  Q5 buildMatrixMethods.py:814-815 against :1149-1170. The insert size is
//     patched in the worker process and the record is written from the
//     master's own copy, so --outBam never carries the patched isize.
//     Reproduced: the port does not patch it either.
//  Q6 buildMatrixMethods.py:180. The --region filter on the restriction cut
//     file demands `region_end <= bed_end`, and region_end defaults to the
//     chromosome length, so --region chrX leaves the cut site list empty and
//     every close inward pair is classified as "same fragment". Reproduced.
//  Q7 buildMatrixMethods.py:1379. The mcool branch is taken only for
//     `len(--binSize) > 2`. Two resolutions therefore write a single cool at
//     the file root of the .mcool file, because ".mcool".endswith("cool").
//     Reproduced.
//  Q10 hicmatrix/lib/cool.py:394-402 against buildMatrixMethods.py:1381.
//     The cool writer deletes the provenance keys out of the metadata dict it
//     is given, and the mcool branch reuses that one dict for every
//     resolution, so only the first resolution of an mcool carries
//     matrix-generated-by. Reproduced.
//  Q9 buildMatrixMethods.py:1371-1376. The matrix provenance strings are
//     built with np.string_, so they are bytes, and hicmatrix's cool writer
//     stores str(bytes), which is the repr. The cool attribute really reads
//     b'HiCExplorer-3.7.6'. Reproduced.
//  Q8 buildMatrixMethods.py:1116-1170. The master consumes worker queues in
//     completion order, so the --outBam record order depends on --threads.
//     Measured: `--threads 8 --inputBufferSize 20000` on
//     small_test_R1/R2_unsorted.bam produces the same 74,642 records as
//     `--threads 2` in a different order. NOT reproduced, because it cannot
//     be: the port writes input order, which is what the Python produces with
//     a single worker.
//
// ---------------------------------------------------------------------------
// Deliberate deviations, all of them from a Python crash
// ---------------------------------------------------------------------------
//  D1 buildMatrixMethods.py:553-557 against :1083. readBamFiles returns
//     (None, None, True, ...) when a call finds no accepted pair left, which
//     happens whenever the number of accepted pairs is an exact multiple of
//     --inputBufferSize. createMatrix then calls len() on that None and dies
//     with a TypeError. Measured: R1_1000/R2_1000 accept 300 pairs, so
//     --inputBufferSize 50 reproduces it and 40 does not. The port completes
//     the run instead. Pinned in hicexplorer/test/general/test_hicBuildMatrix.py.
//  D2 buildMatrixMethods.py:907-910 against :1169. With --doTestRun the output
//     BAM is never opened, but the write loop is still entered when --outBam
//     was given, so `hicBuildMatrix --doTestRun --outBam x.bam` dies with an
//     UnboundLocalError. The port ignores --outBam in test run mode, which is
//     what the guard at :907 intends.
//  Neither is fixed silently: both are reported and both are pinned.
//
// ---------------------------------------------------------------------------
// Not produced
// ---------------------------------------------------------------------------
// The QC folder holds six text artifacts and six rendered ones. The port
// writes QC.log and the five *_table.txt files that hicPrepareQCreport derives
// from it; the five PNGs and hicQC.html are matplotlib figures and a pandas
// Styler HTML table and belong to the tier 7 plot work on hicPrepareQCreport
// itself. See the report.

// The body of the tool lives in build_matrix_impl.hpp, which hicBuildMatrixMicroC
// shares: the two Python tools differ in their argument parser and call the
// same lib/buildMatrixMethods.createMatrix.

#include "build_matrix_impl.hpp"
#include "hicx/argparse.hpp"

namespace {

const char* const kUsage =
    "usage: hicBuildMatrix --samFiles two sam files two sam files --outFileName\n"
    "                      FILENAME --QCfolder FOLDER --restrictionCutFile BED file\n"
    "                      [BED file ...] --restrictionSequence RESTRICTIONSEQUENCE\n"
    "                      [RESTRICTIONSEQUENCE ...] --danglingSequence\n"
    "                      DANGLINGSEQUENCE [DANGLINGSEQUENCE ...] [--outBam bam file]\n"
    "                      [--binSize BINSIZE [BINSIZE ...]] [--minDistance MINDISTANCE]\n"
    "                      [--maxDistance MAXDISTANCE]\n"
    "                      [--maxLibraryInsertSize MAXLIBRARYINSERTSIZE]\n"
    "                      [--genomeAssembly GENOMEASSEMBLY] [--region CHR:START-END]\n"
    "                      [--keepSelfLigation] [--keepSelfCircles]\n"
    "                      [--minMappingQuality MINMAPPINGQUALITY] [--threads THREADS]\n"
    "                      [--inputBufferSize INPUTBUFFERSIZE] [--doTestRun]\n"
    "                      [--doTestRunLines DOTESTRUNLINES] [--skipDuplicationCheck]\n"
    "                      [--chromosomeSizes txt file] [--help] [--version]\n";

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
    "  --restrictionCutFile BED file [BED file ...], -rs BED file [BED file ...]\n"
    "                        BED file(s) with all restriction cut sites.\n"
    "  --restrictionSequence RESTRICTIONSEQUENCE [RESTRICTIONSEQUENCE ...], -seq ...\n"
    "                        Sequence of the restriction site.\n"
    "  --danglingSequence DANGLINGSEQUENCE [DANGLINGSEQUENCE ...]\n"
    "                        Sequence left by the restriction enzyme after cutting.\n"
    "\n"
    "Optional arguments:\n"
    "  --outBam bam file, -b bam file\n"
    "                        Output bam file with all valid Hi-C reads.\n"
    "  --binSize BINSIZE [BINSIZE ...], -bs BINSIZE [BINSIZE ...]\n"
    "                        Size in bp for the bins. (default: None)\n"
    "  --minDistance MINDISTANCE\n"
    "                        Minimum distance between restriction sites.\n"
    "                        (default: 300)\n"
    "  --maxDistance MAXDISTANCE\n"
    "                        Obsolete. Use --maxLibraryInsertSize instead.\n"
    "  --maxLibraryInsertSize MAXLIBRARYINSERTSIZE\n"
    "                        The maximum library insert size. (default: 1000)\n"
    "  --genomeAssembly GENOMEASSEMBLY, -ga GENOMEASSEMBLY\n"
    "                        The genome the reads were mapped to.\n"
    "  --region CHR:START-END, -r CHR:START-END\n"
    "                        Region of the genome to limit the operation to.\n"
    "  --keepSelfLigation    Keep self ligations. (default: False)\n"
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

// hicBuildMatrix.py parse_arguments.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicBuildMatrix",
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
    required.add({"--restrictionCutFile", "-rs"})
        .file_type("r")
        .metavar("BED file")
        .nargs("+")
        .required()
        .input({"bed"})
        .help("BED file(s) with all restriction cut sites.");
    required.add({"--restrictionSequence", "-seq"})
        .type("str")
        .nargs("+")
        .required()
        .help("Sequence of the restriction site.");
    required.add({"--danglingSequence"})
        .type("str")
        .nargs("+")
        .required()
        .help("Sequence left by the restriction enzyme after cutting.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--outBam", "-b"})
        .metavar("bam file")
        .file_type("w")
        .output({"bam"})
        .help("Output bam file with all valid Hi-C reads.");
    optional.add({"--binSize", "-bs"}).type("int").nargs("+").help("Size in bp for the bins.");
    optional.add({"--minDistance"})
        .type("int")
        .default_value(300)
        .help("Minimum distance between restriction sites.");
    optional.add({"--maxDistance"})
        .type("int")
        .help("Obsolete. Use --maxLibraryInsertSize instead.");
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
    optional.add({"--keepSelfLigation"})
        .action(cli::Action::StoreTrue)
        .help("Keep self ligations.");
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
    args.sam_files = ns.strs("samFiles");
    args.out_file_name = ns.str("outFileName");
    args.qc_folder = ns.str("QCfolder");
    args.restriction_cut_files = ns.strs("restrictionCutFile");
    args.restriction_sequences = ns.strs("restrictionSequence");
    args.dangling_sequences = ns.strs("danglingSequence");
    args.out_bam = ns.opt_str("outBam").value_or("");
    args.bin_size = ns.integers("binSize");
    args.min_distance = ns.integer("minDistance");
    args.max_distance = ns.opt_integer("maxDistance");
    args.max_library_insert_size = ns.integer("maxLibraryInsertSize");
    args.genome_assembly = ns.opt_str("genomeAssembly").value_or("");
    args.region = hicx::normalise_region(ns.opt_str("region").value_or(""));
    args.keep_self_ligation = ns.flag("keepSelfLigation");
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
    g_tool = "hicBuildMatrix";
    return run_build_matrix(parse_arguments(argc, argv));
}
