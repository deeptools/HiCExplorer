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
// and needs no external binary.
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

Arguments parse_arguments(int argc, char** argv) {
    Arguments args;
    bool fasta_seen = false;
    bool pattern_seen = false;
    bool out_seen = false;
    bool collecting_patterns = false;
    std::string* pending = nullptr;

    for (int i = 1; i < argc; ++i) {
        const std::string token(argv[i]);
        if (pending != nullptr) {
            *pending = token;
            pending = nullptr;
            continue;
        }
        const bool is_option = token.size() > 1 && token[0] == '-' &&
                               std::isdigit(static_cast<unsigned char>(token[1])) == 0;
        if (!is_option) {
            if (collecting_patterns) {
                args.patterns.push_back(token);
                continue;
            }
            fail("unrecognized arguments: " + token);
        }
        collecting_patterns = false;
        std::string name = token;
        std::optional<std::string> inline_value;
        const std::size_t equals = token.find('=');
        if (equals != std::string::npos && token.rfind("--", 0) == 0) {
            name = token.substr(0, equals);
            inline_value = token.substr(equals + 1);
        }
        if (name == "-h" || name == "--help") {
            std::fputs(kUsage, stdout);
            std::fputs(kHelp, stdout);
            std::exit(0);
        }
        if (name == "--version") {
            std::printf("hicFindRestSite %s\n", hicx::kVersion);
            std::exit(0);
        }
        if (name == "-f" || name == "--fasta") {
            fasta_seen = true;
            if (inline_value.has_value()) {
                args.fasta = *inline_value;
            } else {
                pending = &args.fasta;
            }
            continue;
        }
        if (name == "-o" || name == "--outFile") {
            out_seen = true;
            if (inline_value.has_value()) {
                args.out_file = *inline_value;
            } else {
                pending = &args.out_file;
            }
            continue;
        }
        if (name == "-p" || name == "--searchPattern") {
            pattern_seen = true;
            if (inline_value.has_value()) {
                args.patterns.push_back(*inline_value);
            } else {
                collecting_patterns = true;
            }
            continue;
        }
        fail("unrecognized arguments: " + token);
    }
    if (pending != nullptr) {
        fail("expected one argument");
    }
    std::string missing;
    const auto add_missing = [&missing](const char* name) {
        missing += missing.empty() ? name : std::string(", ") + name;
    };
    if (!fasta_seen) {
        add_missing("--fasta/-f");
    }
    if (!pattern_seen || args.patterns.empty()) {
        add_missing("--searchPattern/-p");
    }
    if (!out_seen) {
        add_missing("--outFile/-o");
    }
    if (!missing.empty()) {
        fail("the following arguments are required: " + missing);
    }
    return args;
}

struct Site {
    std::string chrom;
    std::int64_t start = 0;
    std::int64_t end = 0;
    char strand = '+';
    std::size_t order = 0;  // position in the temporary file
};

bool is_simple_pattern(const std::string& pattern) {
    for (char c : pattern) {
        if (std::isalpha(static_cast<unsigned char>(c)) == 0 && c != '.') {
            return false;
        }
    }
    return !pattern.empty();
}

// re.finditer(pattern, sequence, re.IGNORECASE) for a pattern of letters and
// '.': non-overlapping, left to right, resuming at the end of a match.
void find_simple(const std::string& pattern, const std::string& sequence,
                 const std::function<void(std::size_t, std::size_t)>& emit) {
    const std::size_t width = pattern.size();
    if (sequence.size() < width) {
        return;
    }
    std::string upper = pattern;
    for (char& c : upper) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    const std::size_t last = sequence.size() - width;
    for (std::size_t i = 0; i <= last;) {
        std::size_t k = 0;
        for (; k < width; ++k) {
            if (upper[k] == '.') {
                continue;
            }
            if (static_cast<char>(std::toupper(
                    static_cast<unsigned char>(sequence[i + k]))) != upper[k]) {
                break;
            }
        }
        if (k == width) {
            emit(i, i + width);
            i += width;
        } else {
            ++i;
        }
    }
}

void find_regex(const std::string& pattern, const std::string& sequence,
                const std::function<void(std::size_t, std::size_t)>& emit) {
    const std::regex expression(pattern, std::regex::ECMAScript | std::regex::icase);
    auto it = std::sregex_iterator(sequence.begin(), sequence.end(), expression);
    const auto end = std::sregex_iterator();
    for (; it != end; ++it) {
        const std::smatch& match = *it;
        emit(static_cast<std::size_t>(match.position(0)),
             static_cast<std::size_t>(match.position(0) + match.length(0)));
    }
}

void find_all(const std::string& pattern, const std::string& sequence,
              const std::function<void(std::size_t, std::size_t)>& emit) {
    if (is_simple_pattern(pattern)) {
        find_simple(pattern, sequence, emit);
        return;
    }
    find_regex(pattern, sequence, emit);
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);

    // argparse.FileType('r') opens the fasta while parsing, so a missing file
    // is an exit 2 from the parser and not a later error.
    {
        std::ifstream probe(args.fasta, std::ios::binary);
        if (!probe) {
            fail("argument --fasta/-f: can't open '" + args.fasta +
                 "': [Errno 2] No such file or directory: '" + args.fasta + "'");
        }
    }
    // argparse.FileType('w') creates and truncates the output before anything
    // is computed. Reproduced so that a failing run leaves the same trace.
    {
        std::ofstream create(args.out_file, std::ios::binary | std::ios::trunc);
        if (!create) {
            fail("argument --outFile/-o: can't open '" + args.out_file + "'");
        }
    }

    try {
        std::vector<Site> sites;
        const bool gzipped = hicx::fasta::looks_gzipped(args.fasta);

        for (const std::string& pattern : args.patterns) {
            const std::string rev_compl = hicx::fasta::reverse_complement(pattern);
            // The fasta is reopened once per pattern, exactly as the Python
            // does (hicFindRestSite.py:111): the file order, and with it the
            // order that decides which site survives -u, depends on it.
            hicx::fasta::read_fasta(args.fasta, gzipped,
                             [&](const std::string& name, const std::string& sequence) {
                                 find_all(pattern, sequence,
                                          [&](std::size_t begin, std::size_t end) {
                                              sites.push_back(
                                                  Site{name,
                                                       static_cast<std::int64_t>(begin),
                                                       static_cast<std::int64_t>(end), '+',
                                                       sites.size()});
                                          });
                                 if (rev_compl != pattern) {
                                     find_all(rev_compl, sequence,
                                              [&](std::size_t begin, std::size_t end) {
                                                  sites.push_back(Site{
                                                      name,
                                                      static_cast<std::int64_t>(begin),
                                                      static_cast<std::int64_t>(end), '-',
                                                      sites.size()});
                                              });
                                 }
                             });
        }

        // sort -k1,1 -k2,2n, stable, so that -u keeps the first of an equal run.
        std::stable_sort(sites.begin(), sites.end(), [](const Site& a, const Site& b) {
            if (a.chrom != b.chrom) {
                return a.chrom < b.chrom;
            }
            return a.start < b.start;
        });

        std::ofstream out(args.out_file, std::ios::binary | std::ios::trunc);
        if (!out) {
            std::fprintf(stderr, "hicFindRestSite: cannot write '%s'\n",
                         args.out_file.c_str());
            return 1;
        }
        for (std::size_t i = 0; i < sites.size(); ++i) {
            if (i > 0 && sites[i].chrom == sites[i - 1].chrom &&
                sites[i].start == sites[i - 1].start) {
                continue;  // -u de-duplicates on the key only
            }
            out << sites[i].chrom << '\t' << sites[i].start << '\t' << sites[i].end
                << "\t.\t0\t" << sites[i].strand << '\n';
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicFindRestSite: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicFindRestSite");
    return 0;
}
