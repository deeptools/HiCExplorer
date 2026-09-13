#include "hicx/fasta_reader.hpp"

#include <zlib.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <memory>
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

}  // namespace hicx::fasta
