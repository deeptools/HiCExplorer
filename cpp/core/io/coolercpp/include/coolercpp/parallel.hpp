// A fixed pool of worker threads for the parallel parts of the threaded parts of the library:
// decoding .hic blocks and compressing and decompressing HDF5 chunks. The
// calling thread takes part in every job, so a pool of n threads starts n - 1
// workers, and a pool of one thread runs everything inline.

#ifndef COOLERCPP_PARALLEL_HPP
#define COOLERCPP_PARALLEL_HPP

#include <condition_variable>
#include <cstddef>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace coolercpp {

// The number of CPUs the process may use.
[[nodiscard]] int available_threads();

class ThreadPool {
  public:
    explicit ThreadPool(int threads);
    ~ThreadPool();
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    [[nodiscard]] int threads() const noexcept { return threads_; }

    // Runs task(i) for i in [0, count) and returns when all have finished.
    // Tasks are taken in ascending order, one at a time, by the workers and
    // the calling thread. The first exception a task throws is rethrown here,
    // after the remaining tasks have been skipped.
    void run(std::size_t count, const std::function<void(std::size_t)>& task);

  private:
    struct Job {
        const std::function<void(std::size_t)>* task = nullptr;
        std::size_t count = 0;
        std::size_t next = 0;
        std::size_t active = 0;
        std::exception_ptr error;
    };

    void work();
    // Takes and runs tasks of the current job until none is left. Called with
    // the lock held; returns with it held.
    void drain(std::unique_lock<std::mutex>& lock);

    int threads_;
    std::vector<std::thread> workers_;
    std::mutex mutex_;
    std::condition_variable wake_;
    std::condition_variable done_;
    Job* job_ = nullptr;
    bool stop_ = false;
};

}  // namespace coolercpp

#endif  // COOLERCPP_PARALLEL_HPP
