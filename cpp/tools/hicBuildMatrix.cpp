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

Arguments parse_arguments(int argc, char** argv) {
    Arguments args;
    bool sam_seen = false;
    bool out_seen = false;
    bool qc_seen = false;
    bool cut_seen = false;
    bool seq_seen = false;
    bool dangling_seen = false;

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
            std::printf("hicBuildMatrix %s\n", hicx::kVersion);
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
        } else if (name == "-rs" || name == "--restrictionCutFile") {
            args.restriction_cut_files = take_many("--restrictionCutFile/-rs", 1);
            cut_seen = true;
        } else if (name == "-seq" || name == "--restrictionSequence") {
            args.restriction_sequences = take_many("--restrictionSequence/-seq", 1);
            seq_seen = true;
        } else if (name == "--danglingSequence") {
            args.dangling_sequences = take_many("--danglingSequence", 1);
            dangling_seen = true;
        } else if (name == "-b" || name == "--outBam") {
            args.out_bam = take_one("--outBam/-b");
        } else if (name == "-bs" || name == "--binSize") {
            for (const auto& value : take_many("--binSize/-bs", 1)) {
                args.bin_size.push_back(parse_int("--binSize/-bs", value));
            }
        } else if (name == "--minDistance") {
            args.min_distance = parse_int(name, take_one(name));
        } else if (name == "--maxDistance") {
            args.max_distance = parse_int(name, take_one(name));
        } else if (name == "--maxLibraryInsertSize") {
            args.max_library_insert_size = parse_int(name, take_one(name));
        } else if (name == "-ga" || name == "--genomeAssembly") {
            args.genome_assembly = take_one("--genomeAssembly/-ga");
        } else if (name == "-r" || name == "--region") {
            args.region = hicx::normalise_region(take_one("--region/-r"));
        } else if (name == "--keepSelfLigation") {
            args.keep_self_ligation = true;
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
    require(cut_seen, "--restrictionCutFile/-rs");
    require(seq_seen, "--restrictionSequence/-seq");
    require(dangling_seen, "--danglingSequence");
    if (!missing.empty()) {
        fail("the following arguments are required: " + missing);
    }
    return args;
}

}  // namespace

int main(int argc, char** argv) {
    g_tool = "hicBuildMatrix";
    g_usage = kUsage;
    return run_build_matrix(parse_arguments(argc, argv));
}
