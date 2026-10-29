#ifndef HICFILECPP_DETAIL_OUTPUT_HPP
#define HICFILECPP_DETAIL_OUTPUT_HPP

// Little-endian output: a byte buffer, a buffered file that can patch bytes
// it has already written, and a small thread pool whose parallelFor runs one
// function over an index range and waits.

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <functional>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <vector>

namespace hicfilecpp::detail {

class ByteWriter {
public:
    template <class T>
    void put(T value) {
        const auto* p = reinterpret_cast<const char*>(&value);
        bytes.insert(bytes.end(), p, p + sizeof(T));
    }
    void cstr(const std::string& text) {
        bytes.insert(bytes.end(), text.begin(), text.end());
        bytes.push_back('\0');
    }
    size_t size() const { return bytes.size(); }

    std::vector<char> bytes;
};

class OutputFile {
public:
    // Creates or truncates path; with truncateAt, opens the existing file and
    // cuts it at that offset, positioned there.
    OutputFile(const std::string& path, std::optional<int64_t> truncateAt);
    ~OutputFile();
    OutputFile(const OutputFile&) = delete;
    OutputFile& operator=(const OutputFile&) = delete;

    int64_t position() const { return flushed_ + static_cast<int64_t>(buffer_.size()); }
    void write(const char* data, size_t size);
    void write(const ByteWriter& writer) { write(writer.bytes.data(), writer.bytes.size()); }
    void write(const std::vector<char>& data) { write(data.data(), data.size()); }
    void patch(int64_t offset, const char* data, size_t size);
    template <class T>
    void patchValue(int64_t offset, T value) {
        patch(offset, reinterpret_cast<const char*>(&value), sizeof(T));
    }
    void close();

private:
    void flush();
    void writeAt(int64_t offset, const char* data, size_t size);

    std::string path_;
    int fd_ = -1;
    std::vector<char> buffer_;
    int64_t flushed_ = 0;
};

class ThreadPool {
public:
    explicit ThreadPool(int threads);
    ~ThreadPool();
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    int size() const { return static_cast<int>(workers_.size()) + 1; }
    // Runs fn(i) for every i in [0, n) on the pool and the calling thread,
    // and rethrows the first exception after all have finished.
    void parallelFor(size_t n, const std::function<void(size_t)>& fn);

private:
    struct Job {
        const std::function<void(size_t)>* fn = nullptr;
        size_t n = 0;
        std::atomic<size_t> next{0};
        int active = 0;
        std::exception_ptr error;
    };
    void work(Job& job);
    void loop();

    std::vector<std::thread> workers_;
    std::mutex mutex_;
    std::condition_variable wake_;
    std::condition_variable done_;
    Job* job_ = nullptr;
    uint64_t generation_ = 0;
    bool stop_ = false;
};

}  // namespace hicfilecpp::detail

#endif  // HICFILECPP_DETAIL_OUTPUT_HPP
