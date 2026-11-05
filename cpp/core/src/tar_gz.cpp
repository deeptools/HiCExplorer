// See hicx/tar_gz.hpp.

#include "hicx/tar_gz.hpp"

#include <zlib.h>

#include <array>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <vector>

namespace hicx {

namespace {

constexpr std::size_t kBlock = 512;
constexpr std::size_t kRecord = 20 * kBlock;

// An octal field of `width` bytes: width - 1 zero-padded digits and a
// terminating NUL. A value that needs more digits is refused, because a
// truncated field (a size of 8 GiB or more in the 12-byte size field) would
// silently corrupt the archive.
void octal(char* field, std::size_t width, std::uint64_t value) {
    const std::size_t digits = width - 1;
    if (digits < 22 && (value >> (3 * digits)) != 0) {
        throw std::runtime_error("tar header field of " + std::to_string(width) +
                                 " bytes cannot hold " + std::to_string(value));
    }
    for (std::size_t i = digits; i-- > 0;) {
        field[i] = static_cast<char>('0' + (value & 7u));
        value >>= 3;
    }
    field[digits] = '\0';
}

}  // namespace

TarGzWriter::TarGzWriter(const std::string& path) : path_(path) {
    gz_ = gzopen(path.c_str(), "wb9");
    if (gz_ == nullptr) {
        throw std::runtime_error("cannot create " + path);
    }
}

TarGzWriter::~TarGzWriter() {
    if (gz_ != nullptr) {
        gzclose(static_cast<gzFile>(gz_));
    }
}

void TarGzWriter::write(const char* data, std::size_t length) {
    std::size_t offset = 0;
    while (offset < length) {
        const unsigned chunk = static_cast<unsigned>(std::min<std::size_t>(length - offset, 1u << 30));
        if (gzwrite(static_cast<gzFile>(gz_), data + offset, chunk) != static_cast<int>(chunk)) {
            throw std::runtime_error("cannot write " + path_);
        }
        offset += chunk;
    }
    written_ += length;
}

void TarGzWriter::header(const std::string& name, std::uint64_t size, std::int64_t mtime,
                         unsigned mode, char type) {
    std::array<char, kBlock> block{};
    std::memcpy(block.data(), name.data(), std::min<std::size_t>(name.size(), 100));
    octal(block.data() + 100, 8, mode);
    octal(block.data() + 108, 8, 0);
    octal(block.data() + 116, 8, 0);
    octal(block.data() + 124, 12, size);
    octal(block.data() + 136, 12, static_cast<std::uint64_t>(mtime < 0 ? 0 : mtime));
    std::memset(block.data() + 148, ' ', 8);
    block[156] = type;
    std::memcpy(block.data() + 257, "ustar", 6);
    std::memcpy(block.data() + 263, "00", 2);
    unsigned checksum = 0;
    for (const char c : block) {
        checksum += static_cast<unsigned char>(c);
    }
    std::snprintf(block.data() + 148, 7, "%06o", checksum);
    block[154] = '\0';
    block[155] = ' ';
    write(block.data(), block.size());
}

void TarGzWriter::add(const std::string& name, std::string_view content, std::int64_t mtime,
                      unsigned mode) {
    if (name.size() > 100) {
        // pax extended header: "<length> path=<name>\n", where the length
        // counts the whole record including its own digits.
        const std::string body = " path=" + name + "\n";
        std::size_t length = body.size() + 1;
        while (std::to_string(length).size() + body.size() != length) {
            ++length;
        }
        const std::string record = std::to_string(length) + body;
        header("././@PaxHeader", record.size(), mtime, mode, 'x');
        write(record.data(), record.size());
        const std::size_t padding = (kBlock - record.size() % kBlock) % kBlock;
        const std::array<char, kBlock> zeros{};
        write(zeros.data(), padding);
    }
    header(name, content.size(), mtime, mode, '0');
    write(content.data(), content.size());
    const std::size_t padding = (kBlock - content.size() % kBlock) % kBlock;
    const std::array<char, kBlock> zeros{};
    write(zeros.data(), padding);
}

void TarGzWriter::add_file(const std::string& name, const std::string& source, std::int64_t mtime,
                           unsigned mode) {
    std::ifstream in(source, std::ios::binary);
    if (!in) {
        throw std::runtime_error("cannot read " + source);
    }
    const std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    add(name, content, mtime, mode);
}

void TarGzWriter::close() {
    if (gz_ == nullptr) {
        return;
    }
    const std::array<char, kBlock> zeros{};
    write(zeros.data(), kBlock);
    write(zeros.data(), kBlock);
    while (written_ % kRecord != 0) {
        write(zeros.data(), kBlock);
    }
    const int status = gzclose(static_cast<gzFile>(gz_));
    gz_ = nullptr;
    if (status != Z_OK) {
        throw std::runtime_error("cannot close " + path_);
    }
}

}  // namespace hicx
