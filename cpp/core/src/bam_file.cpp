#include "hicx/bam_file.hpp"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <stdexcept>
#include <utility>

namespace hicx {
namespace {

samFile* as_file(void* handle) { return static_cast<samFile*>(handle); }
sam_hdr_t* as_header(void* handle) { return static_cast<sam_hdr_t*>(handle); }

// pysam AlignedSegment.query_alignment_start / query_alignment_end. htslib
// keeps hard clips in the CIGAR but not in SEQ, so a leading hard clip is
// skipped and only soft clips move the boundary, which is exactly what
// pysam's getQueryStart and getQueryEnd do.
std::pair<std::int32_t, std::int32_t> query_alignment_bounds(const bam1_t* record) {
    const std::int32_t seq_len = record->core.l_qseq;
    const std::uint32_t n_cigar = record->core.n_cigar;
    if (n_cigar == 0) {
        return {0, seq_len};
    }
    const std::uint32_t* cigar = bam_get_cigar(record);
    std::int32_t start = 0;
    for (std::uint32_t i = 0; i < n_cigar; ++i) {
        const int op = static_cast<int>(bam_cigar_op(cigar[i]));
        if (op == BAM_CHARD_CLIP) {
            continue;
        }
        if (op == BAM_CSOFT_CLIP) {
            start += static_cast<std::int32_t>(bam_cigar_oplen(cigar[i]));
            continue;
        }
        break;
    }
    std::int32_t end = seq_len;
    for (std::uint32_t i = n_cigar; i-- > 0;) {
        const int op = static_cast<int>(bam_cigar_op(cigar[i]));
        if (op == BAM_CHARD_CLIP) {
            continue;
        }
        if (op == BAM_CSOFT_CLIP) {
            end -= static_cast<std::int32_t>(bam_cigar_oplen(cigar[i]));
            continue;
        }
        break;
    }
    return {start, end};
}

void copy_probe(const bam1_t* record, std::int32_t from, std::int32_t count,
                char* out) {
    const std::uint8_t* seq = bam_get_seq(record);
    for (std::int32_t i = 0; i < count; ++i) {
        const char base = seq_nt16_str[bam_seqi(seq, from + i)];
        // read.seq.upper(): htslib already produces uppercase IUPAC letters,
        // but the Python calls upper() unconditionally, so do the same.
        out[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    }
}

}  // namespace

BamRecord::BamRecord() : record_(bam_init1()) {
    if (record_ == nullptr) {
        throw std::runtime_error("bam_init1 failed");
    }
}

BamRecord::~BamRecord() {
    if (record_ != nullptr) {
        bam_destroy1(record_);
    }
}

BamRecord::BamRecord(BamRecord&& other) noexcept
    : record_(std::exchange(other.record_, nullptr)) {}

BamRecord& BamRecord::operator=(BamRecord&& other) noexcept {
    if (this != &other) {
        if (record_ != nullptr) {
            bam_destroy1(record_);
        }
        record_ = std::exchange(other.record_, nullptr);
    }
    return *this;
}

BamRecord::BamRecord(const BamRecord& other)
    : record_(other.record_ != nullptr ? bam_dup1(other.record_) : bam_init1()) {
    if (record_ == nullptr) {
        throw std::runtime_error("bam_dup1 failed");
    }
}

BamRecord& BamRecord::operator=(const BamRecord& other) {
    if (this != &other) {
        BamRecord copy(other);
        std::swap(record_, copy.record_);
    }
    return *this;
}

std::string BamRecord::qname() const { return std::string(bam_get_qname(record_)); }

std::uint16_t BamRecord::flag() const { return record_->core.flag; }

std::int32_t BamRecord::tid() const { return record_->core.tid; }

std::int32_t BamRecord::pos() const {
    return static_cast<std::int32_t>(record_->core.pos);
}

std::uint8_t BamRecord::mapq() const { return record_->core.qual; }

bool BamRecord::has_tag(const char* tag) const {
    return bam_aux_get(record_, tag) != nullptr;
}

std::string BamRecord::tag_string(const char* tag) const {
    const std::uint8_t* value = bam_aux_get(record_, tag);
    if (value == nullptr) {
        return {};
    }
    const char* text = bam_aux2Z(value);
    return text != nullptr ? std::string(text) : std::string();
}

ReadFields BamRecord::fields() const {
    ReadFields out;
    out.tid = record_->core.tid;
    out.pos = static_cast<std::int32_t>(record_->core.pos);
    out.flag = record_->core.flag;
    out.mapq = record_->core.qual;
    out.seq_len = record_->core.l_qseq;
    const auto [start, end] = query_alignment_bounds(record_);
    out.qlen = end - start;
    const std::int32_t head = std::min<std::int32_t>(kDanglingProbe, out.seq_len);
    copy_probe(record_, 0, head, out.head);
    const std::int32_t tail = std::min<std::int32_t>(kDanglingProbe, out.seq_len);
    copy_probe(record_, out.seq_len - tail, tail, out.tail);
    return out;
}

std::int64_t BamRecord::leading_bases_before_match() const {
    const std::uint32_t n_cigar = record_->core.n_cigar;
    if (n_cigar == 0) {
        return 0;
    }
    const std::uint32_t* cigar = bam_get_cigar(record_);
    const bool reverse = (record_->core.flag & 0x10) != 0;
    std::int64_t sum = 0;
    for (std::uint32_t k = 0; k < n_cigar; ++k) {
        const std::uint32_t entry = reverse ? cigar[n_cigar - 1 - k] : cigar[k];
        if (static_cast<int>(bam_cigar_op(entry)) == BAM_CMATCH) {
            break;
        }
        sum += bam_cigar_oplen(entry);
    }
    return sum;
}

void BamRecord::set_paired_flag_first() { record_->core.flag |= 0x1 | 0x40; }

void BamRecord::set_paired_flag_second() { record_->core.flag |= 0x1 | 0x80; }

void BamRecord::set_mate(std::int32_t mate_tid, std::int32_t mate_pos) {
    record_->core.mtid = mate_tid;
    record_->core.mpos = mate_pos;
}

BamReader::BamReader(const std::string& path) : path_(path) {
    file_ = sam_open(path.c_str(), "rb");
    if (file_ == nullptr) {
        throw std::runtime_error("could not open " + path);
    }
    header_ = sam_hdr_read(as_file(file_));
    if (header_ == nullptr) {
        sam_close(as_file(file_));
        file_ = nullptr;
        throw std::runtime_error("could not read the header of " + path);
    }
    sam_hdr_t* header = as_header(header_);
    const int count = sam_hdr_nref(header);
    references_.reserve(static_cast<std::size_t>(count));
    lengths_.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
        references_.emplace_back(sam_hdr_tid2name(header, i));
        lengths_.push_back(static_cast<std::int64_t>(sam_hdr_tid2len(header, i)));
    }
}

BamReader::~BamReader() {
    if (header_ != nullptr) {
        sam_hdr_destroy(as_header(header_));
    }
    if (file_ != nullptr) {
        sam_close(as_file(file_));
    }
}

bool BamReader::read(BamRecord& record) {
    const int status = sam_read1(as_file(file_), as_header(header_), record.handle());
    if (status >= 0) {
        return true;
    }
    if (status == -1) {
        return false;
    }
    throw std::runtime_error("truncated or malformed BAM: " + path_);
}

BamWriter::BamWriter(const std::string& path, const BamReader& template_reader)
    : path_(path) {
    file_ = sam_open(path.c_str(), "wb");
    if (file_ == nullptr) {
        throw std::runtime_error("could not open " + path + " for writing");
    }
    header_ = sam_hdr_dup(
        const_cast<sam_hdr_t*>(static_cast<const sam_hdr_t*>(template_reader.header())));
    if (header_ == nullptr) {
        throw std::runtime_error("could not duplicate the template header");
    }
    if (sam_hdr_write(as_file(file_), as_header(header_)) < 0) {
        throw std::runtime_error("could not write the header of " + path);
    }
}

BamWriter::~BamWriter() {
    close();
    if (header_ != nullptr) {
        sam_hdr_destroy(as_header(header_));
        header_ = nullptr;
    }
}

void BamWriter::write(const BamRecord& record) {
    if (sam_write1(as_file(file_), as_header(header_), record.handle()) < 0) {
        throw std::runtime_error("could not write a record to " + path_);
    }
}

void BamWriter::close() {
    if (file_ != nullptr) {
        sam_close(as_file(file_));
        file_ = nullptr;
    }
}

std::vector<std::pair<std::string, std::int64_t>> chrom_sizes_of(
    const BamReader& reader) {
    std::vector<std::pair<std::string, std::int64_t>> out;
    for (std::size_t i = 0; i < reader.references().size(); ++i) {
        const std::string& name = reader.references()[i];
        const auto found = std::find_if(
            out.begin(), out.end(),
            [&name](const auto& entry) { return entry.first == name; });
        if (found != out.end()) {
            found->second = reader.lengths()[i];  // OrderedDict: last value wins
        } else {
            out.emplace_back(name, reader.lengths()[i]);
        }
    }
    return out;
}

}  // namespace hicx
