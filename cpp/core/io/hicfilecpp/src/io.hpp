#ifndef HICFILECPP_DETAIL_IO_HPP
#define HICFILECPP_DETAIL_IO_HPP

// Little-endian binary input: positional reads on a file descriptor (safe to
// share between threads), a buffered sequential reader over it, and a bounds
// checked reader over memory.

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "hicfilecpp/errors.hpp"

#if !(defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
#error "hicfilecpp reads and writes .hic files on little-endian hosts only"
#endif

namespace hicfilecpp::detail {

class File {
public:
    explicit File(const std::string& path);
    ~File();
    File(const File&) = delete;
    File& operator=(const File&) = delete;

    int64_t size() const { return size_; }
    const std::string& path() const { return path_; }
    // Reads exactly n bytes or throws.
    void readAt(int64_t offset, char* out, size_t n) const;
    // Reads up to n bytes; fewer only at the end of the file.
    size_t readSomeAt(int64_t offset, char* out, size_t n) const;

private:
    std::string path_;
    int fd_ = -1;
    int64_t size_ = 0;
};

class SequentialReader {
public:
    SequentialReader(const File& file, int64_t offset, size_t bufferSize = size_t{1} << 20);

    int64_t tell() const { return base_ + static_cast<int64_t>(pos_); }
    void seek(int64_t offset);
    void skip(int64_t n) { seek(tell() + n); }
    bool eof();
    void read(char* out, size_t n);
    // A NUL-terminated string; at most maxLength bytes before the NUL.
    std::string cstr(size_t maxLength = size_t{1} << 30);

    template <class T>
    T get() {
        T value;
        read(reinterpret_cast<char*>(&value), sizeof(T));
        return value;
    }

private:
    bool fill(size_t need);

    const File& file_;
    std::vector<char> buffer_;
    int64_t base_ = 0;
    size_t pos_ = 0;
    size_t len_ = 0;
};

class MemReader {
public:
    MemReader(const char* data, size_t size) : data_(data), size_(size) {}

    size_t remaining() const { return size_ - pos_; }
    const char* cursor() const { return data_ + pos_; }
    void advance(size_t n) { pos_ += n; }
    void require(size_t n) const {
        if (size_ - pos_ < n) {
            throw HicError("truncated block: needed " + std::to_string(n) + " more bytes");
        }
    }

    template <class T>
    T get() {
        require(sizeof(T));
        T value;
        std::memcpy(&value, data_ + pos_, sizeof(T));
        pos_ += sizeof(T);
        return value;
    }

    // Unchecked; the caller has called require().
    template <class T>
    T getUnchecked() {
        T value;
        std::memcpy(&value, data_ + pos_, sizeof(T));
        pos_ += sizeof(T);
        return value;
    }

private:
    const char* data_;
    size_t size_;
    size_t pos_ = 0;
};

}  // namespace hicfilecpp::detail

#endif  // HICFILECPP_DETAIL_IO_HPP
