// A persistent pool of worker threads for kernels that are called many times
// in a row, such as the sparse products of a Krylov solver or the per column
// matrix-vector product of a Hessenberg reduction.
//
// hicx::parallel_for (parallel.hpp) starts and joins its threads on every
// call, which is right for one pass over a matrix and too expensive for a
// kernel called ten thousand times. The contract is the same: run(count, body)
// calls body(index) for every index in [0, count) and returns when all are
// done; the bodies write disjoint memory, so which thread runs which index
// changes no result (cpp/OPTIMIZATION.md section 3).

#ifndef HICX_WORKER_POOL_HPP
#define HICX_WORKER_POOL_HPP

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace hicx {

class WorkerPool {
  public:
    // `workers` threads in total, the calling thread included.
    explicit WorkerPool(int workers) : workers_(std::max(1, workers)) {
        threads_.reserve(static_cast<std::size_t>(workers_ - 1));
        for (int w = 1; w < workers_; ++w) {
            threads_.emplace_back([this] { loop(); });
        }
    }
    WorkerPool(const WorkerPool&) = delete;
    WorkerPool& operator=(const WorkerPool&) = delete;
    ~WorkerPool() {
        {
            const std::lock_guard<std::mutex> lock(mutex_);
            stop_ = true;
        }
        wake_.notify_all();
        for (std::thread& thread : threads_) {
            thread.join();
        }
    }

    [[nodiscard]] int workers() const noexcept { return workers_; }

    void run(int count, const std::function<void(int)>& body) {
        if (workers_ == 1 || count <= 1) {
            for (int index = 0; index < count; ++index) {
                body(index);
            }
            return;
        }
        {
            const std::lock_guard<std::mutex> lock(mutex_);
            body_ = &body;
            count_ = count;
            next_.store(0);
            active_ = static_cast<int>(threads_.size());
            error_ = nullptr;
            ++generation_;
        }
        wake_.notify_all();
        work();
        std::unique_lock<std::mutex> lock(mutex_);
        done_.wait(lock, [this] { return active_ == 0; });
        body_ = nullptr;
        if (error_) {
            std::rethrow_exception(error_);
        }
    }

  private:
    void work() {
        while (true) {
            const int index = next_.fetch_add(1);
            if (index >= count_) {
                return;
            }
            try {
                (*body_)(index);
            } catch (...) {
                const std::lock_guard<std::mutex> lock(mutex_);
                if (!error_) {
                    error_ = std::current_exception();
                }
            }
        }
    }

    void loop() {
        std::uint64_t seen = 0;
        while (true) {
            {
                std::unique_lock<std::mutex> lock(mutex_);
                wake_.wait(lock, [&] { return stop_ || generation_ != seen; });
                if (stop_) {
                    return;
                }
                seen = generation_;
            }
            work();
            const std::lock_guard<std::mutex> lock(mutex_);
            if (--active_ == 0) {
                done_.notify_all();
            }
        }
    }

    int workers_;
    std::vector<std::thread> threads_;
    std::mutex mutex_;
    std::condition_variable wake_;
    std::condition_variable done_;
    bool stop_ = false;
    std::uint64_t generation_ = 0;
    const std::function<void(int)>* body_ = nullptr;
    int count_ = 0;
    std::atomic<int> next_{0};
    int active_ = 0;
    std::exception_ptr error_;
};

}  // namespace hicx

#endif  // HICX_WORKER_POOL_HPP
