// Port of hicexplorer/hicFindRestSite.py.
//
// Scans a FASTA for the occurrences of one or more restriction site patterns,
// on both strands, and writes a sorted BED file.
//
// The reference shells the sorting out to GNU sort
// (hicFindRestSite.py:129-132):
//
//     sort -k1,1 -k2,2n -u <tmpfile>          with LC_ALL=' C'
//
// Three properties of that command decide the output and are reproduced here:
//
//   1. The key is the chromosome name compared byte by byte, then the start
//      compared numerically. Sorting the whole line as text would put 10000
//      before 402.
//   2. `-u` de-duplicates on the **key**, not on the line. Two sites that share
//      a chromosome and a start collapse into one even when their end and
//      their strand differ. On dm3 chrM with the patterns AAGCTT and AAGC that
//      is four collisions out of 74 sites.
//   3. GNU sort's merge is stable and `-u` disables the whole-line last-resort
//      comparison, so the survivor of a collision is the first of them in the
//      temporary file. The temporary file is written pattern by pattern and,
//      within a pattern, forward strand before reverse strand, so the survivor
//      is the site of whichever pattern came first on the command line.
//      Verified against GNU coreutils 9.4.
//
//  Note LC_ALL=' C' at hicFindRestSite.py:131: the value has a leading space
//  and is not a valid locale name, so setlocale fails and the process stays in
//  the C locale. The intended byte ordering is what happens, by accident.
//
// The port sorts in memory with std::stable_sort and then keeps the first
// record of every equal-key run, which is exactly the three properties above
// and needs no external binary. The scan and the sort are
// hicx::fasta::find_restriction_sites (core/include/hicx/fasta_reader.hpp),
// which the public Python function hicFindRestSite.find_pattern calls too.
//
// Pattern matching: the reference uses Python's `re` with re.IGNORECASE and
// finds non-overlapping matches left to right. Every pattern in the corpus, and
// the one in the tool's own doctest, consists of letters and '.', so that case
// is matched by a direct scanner. Anything else falls back to std::regex in
// ECMAScript mode, which is not Python's dialect; the difference cannot show up
// on a pattern built from letters, '.', character classes and the usual
// quantifiers, but a pattern using a Python-specific construct would behave
// differently and that is recorded here rather than hidden.
//
// Threading and SIMD: none. The tool reads a FASTA once per pattern and the
// scan is a byte comparison loop the compiler already vectorises; on the
// designated input the C++ run costs 0.00 s of CPU against the Python's 1.5 s,
// all of which is interpreter start-up. cpp/OPTIMIZATION.md 6 asks for a
// measurement before an optimisation and the measurement says there is nothing
// to optimise.

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <optional>
#include <regex>
#include <string>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/fasta_reader.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicFindRestSite --fasta mm10.fa --searchPattern AAGCTT -o "
    "rest_site_positions.bed\n";

const char* const kHelp =
    "\n"
    "Identifies the genomic locations of restriction sites.\n"
    "\n"
    "Required arguments:\n"
    "  --fasta FASTA, -f FASTA\n"
    "                        Path to fasta file for the organism genome.\n"
    "  --searchPattern SEARCHPATTERN [SEARCHPATTERN ...], -p SEARCHPATTERN "
    "[SEARCHPATTERN ...]\n"
    "                        Search pattern. For example, for HindIII this pattern\n"
    "                        is \"AAGCTT\". Both, forward and reverse strand are\n"
    "                        searched for a match. The pattern is a regexp and can\n"
    "                        contain regexp specif syntax (see\n"
    "                        https://docs.python.org/2/library/re.html). For\n"
    "                        example the patternCG..GC will find all occurrence of\n"
    "                        CG followed by any two bases and then GC.\n"
    "  --outFile OUTFILE, -o OUTFILE\n"
    "                        Name for the resulting bed file.\n"
    "\n"
    "Optional arguments:\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::string fasta;
    std::vector<std::string> patterns;
    std::string out_file;
};

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kUsage, stderr);
    std::fprintf(stderr, "hicFindRestSite: error: %s\n", message.c_str());
    std::exit(2);
}

// hicFindRestSite.py parse_arguments.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicFindRestSite", "Identifies the genomic locations of restriction sites. ");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    // argparse.FileType('r') opens the fasta while parsing, so a missing file
    // is an exit 2 from the parser and not a later error.
    required.add({"--fasta", "-f"})
        .file_type("r")
        .required()
        .input({"fasta"})
        .help("Path to fasta file for the organism genome.");
    required.add({"--searchPattern", "-p"})
        .required()
        .nargs("+")
        .help("Search pattern, a regular expression; both strands are searched.");
    required.add({"--outFile", "-o"})
        .file_type("w")
        .required()
        .output({"bed"})
        .help("Name for the resulting bed file.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.fasta = ns.str("fasta");
    args.patterns = ns.strs("searchPattern");
    args.out_file = ns.str("outFile");
    return args;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);

    // argparse.FileType('w') creates and truncates the output before anything
    // is computed. Reproduced so that a failing run leaves the same trace.
    {
        std::ofstream create(args.out_file, std::ios::binary | std::ios::trunc);
        if (!create) {
            fail("argument --outFile/-o: can't open '" + args.out_file + "'");
        }
    }

    try {
        const std::string bed = hicx::fasta::restriction_sites_bed(
            hicx::fasta::find_restriction_sites(args.patterns, args.fasta));
        std::ofstream out(args.out_file, std::ios::binary | std::ios::trunc);
        if (!out) {
            std::fprintf(stderr, "hicFindRestSite: cannot write '%s'\n",
                         args.out_file.c_str());
            return 1;
        }
        out << bed;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicFindRestSite: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicFindRestSite");
    return 0;
}
