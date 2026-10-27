#include "hicx/fasta_reader.hpp"

#include <zlib.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <functional>
#include <memory>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace hicx::fasta {

namespace {

// Bio.Seq._dna_complement_table. Every character outside it maps to itself.
const std::array<char, 256>& complement_table() {
    static const std::array<char, 256> table = [] {
        std::array<char, 256> values{};
        for (int i = 0; i < 256; ++i) {
            values[static_cast<std::size_t>(i)] = static_cast<char>(i);
        }
        const std::string source = "ACGTacgtRYWSMKHBVDNryswmkhbvdn";
        const std::string target = "TGCAtgcaYRWSKMDVBHNyrswkmdvbhn";
        for (std::size_t i = 0; i < source.size(); ++i) {
            values[static_cast<unsigned char>(source[i])] = target[i];
        }
        return values;
    }();
    return table;
}

// A line source that reads either a plain file or a gzip stream.
class LineReader {
  public:
    LineReader(const std::string& path, bool gzipped) : gzipped_(gzipped) {
        if (gzipped_) {
            gz_ = gzopen(path.c_str(), "rb");
            if (gz_ == nullptr) {
                throw std::runtime_error("cannot open '" + path + "'");
            }
            return;
        }
        plain_.open(path, std::ios::binary);
        if (!plain_) {
            throw std::runtime_error("cannot open '" + path + "'");
        }
    }
    LineReader(const LineReader&) = delete;
    LineReader& operator=(const LineReader&) = delete;
    ~LineReader() {
        if (gz_ != nullptr) {
            gzclose(gz_);
        }
    }

    bool next(std::string* line) {
        if (!gzipped_) {
            return static_cast<bool>(std::getline(plain_, *line));
        }
        line->clear();
        char buffer[65536];
        while (true) {
            if (gzgets(gz_, buffer, sizeof(buffer)) == nullptr) {
                return !line->empty();
            }
            const std::size_t length = std::char_traits<char>::length(buffer);
            if (length > 0 && buffer[length - 1] == '\n') {
                line->append(buffer, length - 1);
                return true;
            }
            line->append(buffer, length);
        }
    }

  private:
    bool gzipped_ = false;
    std::ifstream plain_;
    gzFile gz_ = nullptr;
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

std::string clean_sequence_line(const std::string& line) {
    // SimpleFastaParser rstrips each line and then removes spaces and carriage
    // returns from the joined sequence.
    std::string cleaned;
    cleaned.reserve(line.size());
    for (char c : line) {
        if (c == ' ' || c == '\r') {
            continue;
        }
        cleaned.push_back(c);
    }
    while (!cleaned.empty() &&
           std::isspace(static_cast<unsigned char>(cleaned.back())) != 0) {
        cleaned.pop_back();
    }
    return cleaned;
}

}  // namespace

bool looks_gzipped(const std::string& path) {
    return path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
}

std::string reverse_complement(const std::string& text) {
    const std::array<char, 256>& table = complement_table();
    std::string result;
    result.reserve(text.size());
    for (auto it = text.rbegin(); it != text.rend(); ++it) {
        result.push_back(table[static_cast<unsigned char>(*it)]);
    }
    return result;
}

void read_fasta(const std::string& path, bool gzipped,
                const std::function<void(const std::string&, const std::string&)>& visit) {
    LineReader reader(path, gzipped);
    std::string line;
    std::string name;
    std::string sequence;
    bool in_record = false;

    while (reader.next(&line)) {
        if (!line.empty() && line[0] == '>') {
            if (in_record) {
                visit(name, sequence);
            }
            // record.name is the header up to the first whitespace.
            std::string header = line.substr(1);
            const std::size_t space = header.find_first_of(" \t\r\n");
            name = space == std::string::npos ? header : header.substr(0, space);
            while (!name.empty() &&
                   std::isspace(static_cast<unsigned char>(name.back())) != 0) {
                name.pop_back();
            }
            sequence.clear();
            in_record = true;
            continue;
        }
        if (!in_record) {
            continue;  // text before the first header, which SeqIO ignores
        }
        sequence += clean_sequence_line(line);
    }
    if (in_record) {
        visit(name, sequence);
    }
}

std::vector<RestrictionSite> find_restriction_sites(const std::vector<std::string>& patterns,
                                                   const std::string& fasta_path) {
    std::vector<RestrictionSite> sites;
    const bool gzipped = looks_gzipped(fasta_path);
    for (const std::string& pattern : patterns) {
        const std::string rev_compl = reverse_complement(pattern);
        read_fasta(fasta_path, gzipped, [&](const std::string& name, const std::string& sequence) {
            find_all(pattern, sequence, [&](std::size_t begin, std::size_t end) {
                sites.push_back(RestrictionSite{name, static_cast<std::int64_t>(begin),
                                                static_cast<std::int64_t>(end), '+'});
            });
            if (rev_compl != pattern) {
                find_all(rev_compl, sequence, [&](std::size_t begin, std::size_t end) {
                    sites.push_back(RestrictionSite{name, static_cast<std::int64_t>(begin),
                                                    static_cast<std::int64_t>(end), '-'});
                });
            }
        });
    }
    // sort -k1,1 -k2,2n, stable, so that -u keeps the first of an equal run.
    std::stable_sort(sites.begin(), sites.end(), [](const RestrictionSite& a, const RestrictionSite& b) {
        if (a.chrom != b.chrom) {
            return a.chrom < b.chrom;
        }
        return a.start < b.start;
    });
    std::vector<RestrictionSite> unique;
    unique.reserve(sites.size());
    for (RestrictionSite& site : sites) {
        // -u de-duplicates on the key only; equal keys are adjacent
        if (!unique.empty() && site.chrom == unique.back().chrom && site.start == unique.back().start) {
            continue;
        }
        unique.push_back(std::move(site));
    }
    return unique;
}

std::string restriction_sites_bed(const std::vector<RestrictionSite>& sites) {
    std::ostringstream out;
    for (const RestrictionSite& site : sites) {
        out << site.chrom << '\t' << site.start << '\t' << site.end << "\t.\t0\t" << site.strand << '\n';
    }
    return out.str();
}

}  // namespace hicx::fasta
