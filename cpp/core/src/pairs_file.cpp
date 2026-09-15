#include "hicx/pairs_file.hpp"

#include <htslib/hts.h>
#include <htslib/kstring.h>

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <limits>
#include <set>
#include <stdexcept>

namespace hicx {
namespace {

std::string_view trim(std::string_view text) {
    while (!text.empty() && std::isspace(static_cast<unsigned char>(text.front())) != 0) {
        text.remove_prefix(1);
    }
    while (!text.empty() && std::isspace(static_cast<unsigned char>(text.back())) != 0) {
        text.remove_suffix(1);
    }
    return text;
}

std::vector<std::string> split_whitespace(std::string_view text) {
    std::vector<std::string> out;
    std::size_t i = 0;
    while (i < text.size()) {
        while (i < text.size() && std::isspace(static_cast<unsigned char>(text[i])) != 0) {
            ++i;
        }
        std::size_t j = i;
        while (j < text.size() && std::isspace(static_cast<unsigned char>(text[j])) == 0) {
            ++j;
        }
        if (j > i) {
            out.emplace_back(text.substr(i, j - i));
        }
        i = j;
    }
    return out;
}

std::string canonical_column(const std::string& name) {
    if (name == "chr1") {
        return "chrom1";
    }
    if (name == "chr2") {
        return "chrom2";
    }
    return name;
}

// The value of a "#key: value" header line.
bool header_field(std::string_view line, std::string_view key, std::string_view& value) {
    if (line.size() < key.size() + 2 || line.substr(1, key.size()) != key ||
        line[1 + key.size()] != ':') {
        return false;
    }
    value = trim(line.substr(key.size() + 2));
    return true;
}

bool parse_non_negative(std::string_view text, std::int64_t& value) {
    if (text.empty()) {
        return false;
    }
    const auto [ptr, ec] = std::from_chars(text.data(), text.data() + text.size(), value);
    return ec == std::errc() && ptr == text.data() + text.size() && value >= 0;
}

std::int8_t strand_of(std::string_view text) {
    if (text == "+") {
        return 1;
    }
    if (text == "-") {
        return -1;
    }
    return 0;
}

}  // namespace

bool PairsHeader::sorted_by_chroms_then_positions() const {
    return sorted_keys == std::vector<std::string>{"chrom1", "chrom2", "pos1", "pos2"};
}

bool PairsHeader::upper_triangle() const { return shape == "upper triangle"; }

PairsHeader parse_pairs_header(const std::vector<std::string>& lines) {
    PairsHeader header;
    const std::string_view magic = "## pairs format v";
    if (lines.empty() || std::string_view(lines.front()).substr(0, magic.size()) != magic) {
        throw std::runtime_error(
            "not a .pairs file: its first line must be '## pairs format v1.0' "
            "(4DN pairs specification)");
    }
    header.format_version = std::string(trim(std::string_view(lines.front()).substr(magic.size())));

    std::set<std::string> chromosomes;
    for (std::size_t i = 1; i < lines.size(); ++i) {
        const std::string_view line = lines[i];
        std::string_view value;
        if (header_field(line, "sorted", value)) {
            header.sorted_keys.clear();
            std::size_t start = 0;
            while (start <= value.size()) {
                const std::size_t dash = value.find('-', start);
                const std::string_view key = value.substr(
                    start, dash == std::string_view::npos ? std::string_view::npos : dash - start);
                header.sorted_keys.push_back(canonical_column(std::string(trim(key))));
                if (dash == std::string_view::npos) {
                    break;
                }
                start = dash + 1;
            }
        } else if (header_field(line, "shape", value)) {
            header.shape = std::string(value);
        } else if (header_field(line, "chromsize", value)) {
            const std::vector<std::string> parts = split_whitespace(value);
            std::int64_t length = 0;
            if (parts.size() != 2 || !parse_non_negative(parts[1], length) || length == 0) {
                throw std::runtime_error("malformed #chromsize header line " +
                                         std::to_string(i + 1) + ": " + std::string(line));
            }
            if (!chromosomes.insert(parts[0]).second) {
                throw std::runtime_error("chromosome " + parts[0] +
                                         " has more than one #chromsize header line");
            }
            header.chrom_sizes.emplace_back(parts[0], length);
        } else if (header_field(line, "columns", value)) {
            header.columns.clear();
            for (const std::string& name : split_whitespace(value)) {
                header.columns.push_back(canonical_column(name));
            }
            header.columns_declared = true;
        } else if (header_field(line, "genome_assembly", value)) {
            header.genome_assembly = std::string(value);
        }
    }
    if (!header.columns_declared) {
        header.columns = {"readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"};
    }

    std::set<std::string> seen;
    int largest = -1;
    for (std::size_t k = 0; k < header.columns.size(); ++k) {
        const std::string& name = header.columns[k];
        if (!seen.insert(name).second) {
            throw std::runtime_error("the #columns header line names " + name + " twice");
        }
        const int index = static_cast<int>(k);
        int* slot = nullptr;
        if (name == "chrom1") {
            slot = &header.chrom1;
        } else if (name == "pos1") {
            slot = &header.pos1;
        } else if (name == "chrom2") {
            slot = &header.chrom2;
        } else if (name == "pos2") {
            slot = &header.pos2;
        } else if (name == "strand1") {
            slot = &header.strand1;
        } else if (name == "strand2") {
            slot = &header.strand2;
        } else if (name == "mapq1") {
            slot = &header.mapq1;
        } else if (name == "mapq2") {
            slot = &header.mapq2;
        } else if (name == "pair_type") {
            slot = &header.pair_type;
        }
        if (slot != nullptr) {
            *slot = index;
            largest = std::max(largest, index);
        }
    }
    if (header.chrom1 < 0 || header.pos1 < 0 || header.chrom2 < 0 || header.pos2 < 0) {
        throw std::runtime_error(
            "the #columns header line must name chrom1 (or chr1), pos1, chrom2 (or chr2) "
            "and pos2");
    }
    header.fields_needed = static_cast<std::size_t>(largest + 1);
    return header;
}

std::optional<std::string> parse_pairs_record(std::string_view line, const PairsHeader& header,
                                              std::vector<std::string_view>& fields,
                                              PairsRecord& out) {
    fields.clear();
    std::size_t start = 0;
    while (fields.size() < header.fields_needed) {
        const std::size_t tab = line.find('\t', start);
        if (tab == std::string_view::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }
    if (fields.size() < header.fields_needed) {
        return "has " + std::to_string(fields.size()) +
               " tab separated fields, its column layout needs at least " +
               std::to_string(header.fields_needed);
    }

    out.chrom1 = fields[static_cast<std::size_t>(header.chrom1)];
    out.chrom2 = fields[static_cast<std::size_t>(header.chrom2)];
    const std::string_view pos1 = fields[static_cast<std::size_t>(header.pos1)];
    const std::string_view pos2 = fields[static_cast<std::size_t>(header.pos2)];
    if (!parse_non_negative(pos1, out.pos1)) {
        return "pos1 '" + std::string(pos1) + "' is not a non-negative integer";
    }
    if (!parse_non_negative(pos2, out.pos2)) {
        return "pos2 '" + std::string(pos2) + "' is not a non-negative integer";
    }

    out.mapq1 = -1;
    out.mapq2 = -1;
    for (int side = 0; side < 2; ++side) {
        const int column = side == 0 ? header.mapq1 : header.mapq2;
        if (column < 0) {
            continue;
        }
        const std::string_view text = fields[static_cast<std::size_t>(column)];
        std::int64_t value = 0;
        if (!parse_non_negative(text, value)) {
            return std::string(side == 0 ? "mapq1" : "mapq2") + " '" + std::string(text) +
                   "' is not a non-negative integer";
        }
        const auto clipped = static_cast<std::int32_t>(
            std::min<std::int64_t>(value, std::numeric_limits<std::int32_t>::max()));
        (side == 0 ? out.mapq1 : out.mapq2) = clipped;
    }
    out.strand1 = header.strand1 >= 0 ? strand_of(fields[static_cast<std::size_t>(header.strand1)]) : 0;
    out.strand2 = header.strand2 >= 0 ? strand_of(fields[static_cast<std::size_t>(header.strand2)]) : 0;

    out.kind = PairKind::Mapped;
    if (header.pair_type >= 0) {
        const std::string_view type = fields[static_cast<std::size_t>(header.pair_type)];
        if (type == "DD") {
            out.kind = PairKind::MarkedDuplicate;
        } else if (type != ".") {
            bool multi = false;
            bool unmapped = false;
            for (const char letter : type) {
                if (letter == 'M') {
                    multi = true;
                } else if (letter != 'U' && letter != 'R') {
                    unmapped = true;
                }
            }
            if (unmapped) {
                out.kind = PairKind::Unmapped;
            } else if (multi) {
                out.kind = PairKind::NotUnique;
            }
        }
    }
    if (out.kind == PairKind::Mapped && (out.chrom1 == "!" || out.chrom2 == "!")) {
        out.kind = PairKind::Unmapped;
    }
    return std::nullopt;
}

struct PairsReader::Impl {
    htsFile* file = nullptr;
    kstring_t buffer = {0, 0, nullptr};
    bool pending = false;

    ~Impl() {
        if (file != nullptr) {
            hts_close(file);
        }
        std::free(buffer.s);
    }

    // hts_getline: >= 0 a line, -1 the end, below -1 an error.
    int read_line() {
        const int status = hts_getline(file, '\n', &buffer);
        if (status >= 0 && buffer.l > 0 && buffer.s[buffer.l - 1] == '\r') {
            buffer.s[--buffer.l] = '\0';
        }
        return status;
    }
};

PairsReader::PairsReader(const std::string& path, int threads)
    : impl_(std::make_unique<Impl>()), path_(path) {
    impl_->file = hts_open(path.c_str(), "r");
    if (impl_->file == nullptr) {
        throw std::runtime_error("could not open " + path);
    }
    if (threads > 1) {
        // Applies to bgzip input only; plain and gzip input read on one thread.
        hts_set_threads(impl_->file, threads);
    }
    std::vector<std::string> header_lines;
    while (true) {
        const int status = impl_->read_line();
        if (status == -1) {
            break;
        }
        if (status < -1) {
            throw std::runtime_error("error reading " + path + " at line " +
                                     std::to_string(line_number_ + 1));
        }
        ++line_number_;
        if (impl_->buffer.l > 0 && impl_->buffer.s[0] == '#') {
            header_lines.emplace_back(impl_->buffer.s, impl_->buffer.l);
            continue;
        }
        impl_->pending = true;
        break;
    }
    try {
        header_ = parse_pairs_header(header_lines);
    } catch (const std::runtime_error& error) {
        throw std::runtime_error(path + ": " + error.what());
    }
}

PairsReader::~PairsReader() = default;

bool PairsReader::next(std::string_view& line) {
    while (true) {
        if (impl_->pending) {
            impl_->pending = false;
        } else {
            const int status = impl_->read_line();
            if (status == -1) {
                return false;
            }
            if (status < -1) {
                throw std::runtime_error("error reading " + path_ + " at line " +
                                         std::to_string(line_number_ + 1) +
                                         " (truncated or corrupt compressed data?)");
            }
            ++line_number_;
        }
        if (impl_->buffer.l == 0 || impl_->buffer.s[0] == '#') {
            continue;
        }
        line = std::string_view(impl_->buffer.s, impl_->buffer.l);
        return true;
    }
}

}  // namespace hicx
