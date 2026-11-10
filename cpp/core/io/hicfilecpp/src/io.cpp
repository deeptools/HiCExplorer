#include "io.hpp"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>

namespace hicfilecpp::detail {

File::File(const std::string& path) : path_(path) {
    fd_ = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
    if (fd_ < 0) {
        throw HicError("File " + path + " cannot be opened for reading");
    }
    struct stat st {};
    if (::fstat(fd_, &st) != 0) {
        ::close(fd_);
        throw HicError("File " + path + " cannot be opened for reading");
    }
    size_ = static_cast<int64_t>(st.st_size);
}

File::~File() {
    if (fd_ >= 0) {
        ::close(fd_);
    }
}

size_t File::readSomeAt(int64_t offset, char* out, size_t n) const {
    size_t done = 0;
    while (done < n) {
        const ssize_t got = ::pread(fd_, out + done, n - done, static_cast<off_t>(offset) + static_cast<off_t>(done));
        if (got < 0) {
            if (errno == EINTR) {
                continue;
            }
            throw HicError("read error in " + path_ + " at offset " + std::to_string(offset));
        }
        if (got == 0) {
            break;
        }
        done += static_cast<size_t>(got);
    }
    return done;
}

void File::readAt(int64_t offset, char* out, size_t n) const {
    if (offset < 0 || readSomeAt(offset, out, n) != n) {
        throw HicError("unexpected end of file " + path_ + " reading " + std::to_string(n) + " bytes at offset " +
                       std::to_string(offset));
    }
}

SequentialReader::SequentialReader(const File& file, int64_t offset, size_t bufferSize)
    : file_(file), buffer_(bufferSize), base_(offset) {}

void SequentialReader::seek(int64_t offset) {
    if (offset >= base_ && offset <= base_ + static_cast<int64_t>(len_)) {
        pos_ = static_cast<size_t>(offset - base_);
    } else {
        base_ = offset;
        pos_ = 0;
        len_ = 0;
    }
}

bool SequentialReader::fill(size_t need) {
    if (len_ - pos_ >= need) {
        return true;
    }
    const size_t keep = len_ - pos_;
    if (keep > 0) {
        std::memmove(buffer_.data(), buffer_.data() + pos_, keep);
    }
    base_ += static_cast<int64_t>(pos_);
    pos_ = 0;
    len_ = keep;
    if (buffer_.size() < need) {
        buffer_.resize(need);
    }
    if (base_ < 0) {
        return false;
    }
    len_ += file_.readSomeAt(base_ + static_cast<int64_t>(keep), buffer_.data() + keep, buffer_.size() - keep);
    return len_ >= need;
}

bool SequentialReader::eof() {
    return !fill(1);
}

void SequentialReader::read(char* out, size_t n) {
    if (n > buffer_.size() && len_ == pos_) {
        file_.readAt(tell(), out, n);
        base_ += static_cast<int64_t>(pos_ + n);
        pos_ = 0;
        len_ = 0;
        return;
    }
    if (!fill(n)) {
        throw HicError("unexpected end of file " + file_.path() + " at offset " + std::to_string(tell()));
    }
    std::memcpy(out, buffer_.data() + pos_, n);
    pos_ += n;
}

std::string SequentialReader::cstr(size_t maxLength) {
    std::string text;
    while (true) {
        if (!fill(1)) {
            throw HicError("unexpected end of file " + file_.path() + " inside a string");
        }
        const char* begin = buffer_.data() + pos_;
        const size_t available = len_ - pos_;
        const auto* nul = static_cast<const char*>(std::memchr(begin, '\0', available));
        if (nul != nullptr) {
            text.append(begin, static_cast<size_t>(nul - begin));
            pos_ += static_cast<size_t>(nul - begin) + 1;
            break;
        }
        text.append(begin, available);
        pos_ = len_;
        if (text.size() > maxLength) {
            throw HicError("unterminated string in " + file_.path());
        }
    }
    if (text.size() > maxLength) {
        throw HicError("unterminated string in " + file_.path());
    }
    return text;
}

}  // namespace hicfilecpp::detail
