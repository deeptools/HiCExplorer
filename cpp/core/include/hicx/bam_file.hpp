// BAM input and output, the htslib replacement for pysam.
//
// hicBuildMatrix and hicBuildMatrixMicroC are the only tools in the suite that
// read alignments. They use pysam for exactly five things, and nothing else:
//
//   1. pysam.Samfile(path, 'rb') and iteration over it, which is sam_open,
//      sam_hdr_read and sam_read1.
//   2. the reference name and length lists, which are the @SQ lines.
//   3. six fields per record: flag, reference id, position, mapping quality,
//      the CIGAR, and the query sequence.
//   4. two derived quantities: pysam's `qlen` (the aligned length of the
//      query, soft clips excluded) and `len(read.seq)` (the full SEQ length).
//   5. the SA auxiliary tag, to find supplementary alignments.
//
// plus, on the --outBam path, writing records back out with the template
// header and four patched fields.
//
// htslib is not visible in this header. bam1_t is forward declared so that a
// BamRecord can hold one; every accessor is defined in bam_file.cpp. That
// keeps the include of htslib/sam.h in one translation unit and lets the tool
// and the tests be compiled without the htslib include path.
//
// Memory note (cpp/PLAN.md 4.4 rule 8). The Python buffers up to
// --inputBufferSize whole pysam AlignedSegment objects per worker process,
// then forks, so the buffer is resident once per worker. The port buffers a
// PairedReadFields per read, which is 40 bytes and holds only what the
// classification actually reads, and it keeps the full BamRecord only when
// --outBam asks for the records back.

#ifndef HICX_BAM_FILE_HPP
#define HICX_BAM_FILE_HPP

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

struct bam1_t;  // htslib

namespace hicx {

// The longest dangling sequence the classification ever compares against. Only
// the first and last kDanglingProbe bases of a read are needed, so the buffered
// form of a read never stores the sequence itself. HindIII's AGCT is 4 long,
// the longest restriction site in the corpus is AAGCTT at 6.
inline constexpr int kDanglingProbe = 16;

// Everything the read pair classification of buildMatrixMethods.process_data
// reads off one alignment, in a form that costs no allocation.
struct ReadFields {
    std::int32_t tid = -1;
    std::int32_t pos = -1;
    std::uint16_t flag = 0;
    std::uint8_t mapq = 0;
    // pysam AlignedSegment.qlen, that is query_alignment_length: the SEQ
    // length minus the leading and trailing soft clips. A record without a
    // CIGAR reports the full SEQ length, which is what pysam does as well.
    std::int32_t qlen = 0;
    // len(read.seq).
    std::int32_t seq_len = 0;
    // The first and last kDanglingProbe bases of SEQ, uppercase, unpadded.
    // check_dangling_end only ever asks for a prefix or a suffix.
    char head[kDanglingProbe] = {};
    char tail[kDanglingProbe] = {};

    [[nodiscard]] bool is_reverse() const noexcept { return (flag & 0x10) != 0; }
    [[nodiscard]] bool is_unmapped() const noexcept { return (flag & 0x4) != 0; }
    [[nodiscard]] bool is_secondary() const noexcept { return (flag & 0x100) != 0; }
};

// One alignment record, owning its htslib bam1_t.
class BamRecord {
  public:
    BamRecord();
    ~BamRecord();
    BamRecord(BamRecord&& other) noexcept;
    BamRecord& operator=(BamRecord&& other) noexcept;
    BamRecord(const BamRecord& other);
    BamRecord& operator=(const BamRecord& other);

    [[nodiscard]] bam1_t* handle() const noexcept { return record_; }

    [[nodiscard]] std::string qname() const;
    [[nodiscard]] std::uint16_t flag() const;
    [[nodiscard]] std::int32_t tid() const;
    [[nodiscard]] std::int32_t pos() const;
    [[nodiscard]] std::uint8_t mapq() const;
    [[nodiscard]] bool has_tag(const char* tag) const;
    [[nodiscard]] std::string tag_string(const char* tag) const;

    // The plain fields the classification needs, filled in one pass.
    [[nodiscard]] ReadFields fields() const;

    // get_correct_map's per-read key: the number of query bases before the
    // first CIGAR M, with the CIGAR read backwards for a reverse read. A read
    // without a CIGAR contributes 0, as an empty Python cigartuples loop does.
    [[nodiscard]] std::int64_t leading_bases_before_match() const;

    // The four patches buildMatrixMethods.createMatrix applies in the master
    // process before writing a pair to --outBam. The fifth patch, the one to
    // isize, is applied in the worker process and is therefore lost; see
    // hicBuildMatrix.cpp.
    void set_paired_flag_first();
    void set_paired_flag_second();
    void set_mate(std::int32_t mate_tid, std::int32_t mate_pos);

  private:
    bam1_t* record_ = nullptr;
};

// sam_open plus sam_hdr_read, iterated with sam_read1.
class BamReader {
  public:
    explicit BamReader(const std::string& path);
    ~BamReader();
    BamReader(const BamReader&) = delete;
    BamReader& operator=(const BamReader&) = delete;

    // get_chrom_sizes: the @SQ names and lengths in file order.
    [[nodiscard]] const std::vector<std::string>& references() const noexcept {
        return references_;
    }
    [[nodiscard]] const std::vector<std::int64_t>& lengths() const noexcept {
        return lengths_;
    }

    // Reads the next record. False at end of file, which is Python's
    // StopIteration.
    bool read(BamRecord& record);

    // Opaque handle to the header, for BamWriter's template argument.
    [[nodiscard]] const void* header() const noexcept { return header_; }

  private:
    void* file_ = nullptr;
    void* header_ = nullptr;
    std::vector<std::string> references_;
    std::vector<std::int64_t> lengths_;
    std::string path_;
};

// pysam.Samfile(name, 'wb', template=other): a BGZF compressed BAM whose
// header is a verbatim copy of the template's.
class BamWriter {
  public:
    BamWriter(const std::string& path, const BamReader& template_reader);
    ~BamWriter();
    BamWriter(const BamWriter&) = delete;
    BamWriter& operator=(const BamWriter&) = delete;

    void write(const BamRecord& record);
    void close();

  private:
    void* file_ = nullptr;
    void* header_ = nullptr;
    std::string path_;
};

// get_chrom_sizes: OrderedDict(zip(references, lengths)) then list(items()).
// A reference name that appears twice keeps its first position and its last
// length, which is what the OrderedDict does.
[[nodiscard]] std::vector<std::pair<std::string, std::int64_t>> chrom_sizes_of(
    const BamReader& reader);

}  // namespace hicx

#endif  // HICX_BAM_FILE_HPP
