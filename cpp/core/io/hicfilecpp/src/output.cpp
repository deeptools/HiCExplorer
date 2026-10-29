#include "output.hpp"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstring>

#include "hicfilecpp/errors.hpp"

namespace hicfilecpp::detail {

namespace {

constexpr size_t kBufferSize = size_t{4} << 20;

}  // namespace

OutputFile::OutputFile(const std::string& path, std::optional<int64_t> truncateAt) : path_(path) {
    if (truncateAt) {
        fd_ = ::open(path.c_str(), O_RDWR | O_CLOEXEC);
        if (fd_ < 0 || ::ftruncate(fd_, static_cast<off_t>(*truncateAt)) != 0) {
            if (fd_ >= 0) {
                ::close(fd_);
            }
            fd_ = -1;
            throw HicError("File " + path + " cannot be opened for writing");
        }
        flushed_ = *truncateAt;
    } else {
        fd_ = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0644);
        if (fd_ < 0) {
            throw HicError("File " + path + " cannot be opened for writing");
        }
    }
    buffer_.reserve(kBufferSize);
}

OutputFile::~OutputFile() {
    if (fd_ >= 0) {
        ::close(fd_);
    }
}

void OutputFile::writeAt(int64_t offset, const char* data, size_t size) {
    size_t done = 0;
    while (done < size) {
        const ssize_t wrote =
            ::pwrite(fd_, data + done, size - done, static_cast<off_t>(offset) + static_cast<off_t>(done));
        if (wrote < 0) {
            if (errno == EINTR) {
                continue;
            }
            throw HicError("write error in " + path_ + ": " + std::strerror(errno));
        }
        done += static_cast<size_t>(wrote);
    }
}

void OutputFile::flush() {
    if (!buffer_.empty()) {
        writeAt(flushed_, buffer_.data(), buffer_.size());
        flushed_ += static_cast<int64_t>(buffer_.size());
        buffer_.clear();
    }
}

void OutputFile::write(const char* data, size_t size) {
    if (buffer_.size() + size > kBufferSize) {
        flush();
    }
    if (size >= kBufferSize) {
        writeAt(flushed_, data, size);
        flushed_ += static_cast<int64_t>(size);
        return;
    }
    buffer_.insert(buffer_.end(), data, data + size);
}

void OutputFile::patch(int64_t offset, const char* data, size_t size) {
    const int64_t end = offset + static_cast<int64_t>(size);
    if (offset < 0 || end > position()) {
        throw HicError("patch outside the written part of " + path_);
    }
    if (offset < flushed_) {
        const auto inFile = static_cast<size_t>(std::min(end, flushed_) - offset);
        writeAt(offset, data, inFile);
        offset += static_cast<int64_t>(inFile);
        data += inFile;
        size -= inFile;
    }
    if (size > 0) {
        std::memcpy(buffer_.data() + (offset - flushed_), data, size);
    }
}

void OutputFile::close() {
    flush();
    if (::close(fd_) != 0) {
        fd_ = -1;
        throw HicError("close error in " + path_);
    }
    fd_ = -1;
}

ThreadPool::ThreadPool(int threads) {
    for (int i = 1; i < threads; ++i) {
        workers_.emplace_back([this] { loop(); });
    }
}

ThreadPool::~ThreadPool() {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        stop_ = true;
    }
    wake_.notify_all();
    for (auto& worker : workers_) {
        worker.join();
    }
}

void ThreadPool::work(Job& job) {
    while (true) {
        const size_t i = job.next.fetch_add(1);
        if (i >= job.n) {
            return;
        }
        try {
            (*job.fn)(i);
        } catch (...) {
            std::lock_guard<std::mutex> lock(mutex_);
            if (!job.error) {
                job.error = std::current_exception();
            }
        }
    }
}

void ThreadPool::loop() {
    uint64_t seen = 0;
    std::unique_lock<std::mutex> lock(mutex_);
    while (true) {
        wake_.wait(lock, [&] { return stop_ || (job_ != nullptr && generation_ != seen); });
        if (stop_) {
            return;
        }
        Job* job = job_;
        seen = generation_;
        ++job->active;
        lock.unlock();
        work(*job);
        lock.lock();
        if (--job->active == 0) {
            done_.notify_all();
        }
    }
}

void ThreadPool::parallelFor(size_t n, const std::function<void(size_t)>& fn) {
    if (n == 0) {
        return;
    }
    if (workers_.empty() || n == 1) {
        for (size_t i = 0; i < n; ++i) {
            fn(i);
        }
        return;
    }
    Job job;
    job.fn = &fn;
    job.n = n;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        job_ = &job;
        ++generation_;
        ++job.active;
    }
    wake_.notify_all();
    work(job);
    {
        std::unique_lock<std::mutex> lock(mutex_);
        --job.active;
        done_.wait(lock, [&] { return job.active == 0; });
        job_ = nullptr;
    }
    if (job.error) {
        std::rethrow_exception(job.error);
    }
}

}  // namespace hicfilecpp::detail
