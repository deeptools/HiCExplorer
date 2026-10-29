#include "coolercpp/parallel.hpp"

#include <algorithm>

#include <sched.h>

namespace coolercpp {

int available_threads() {
    cpu_set_t set;
    CPU_ZERO(&set);
    if (sched_getaffinity(0, sizeof(set), &set) == 0) {
        const int count = CPU_COUNT(&set);
        if (count > 0) {
            return count;
        }
    }
    return static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
}

ThreadPool::ThreadPool(int threads) : threads_(std::max(1, threads)) {
    for (int i = 1; i < threads_; ++i) {
        workers_.emplace_back([this] { work(); });
    }
}

ThreadPool::~ThreadPool() {
    {
        std::lock_guard<std::mutex> guard(mutex_);
        stop_ = true;
    }
    wake_.notify_all();
    for (auto& worker : workers_) {
        worker.join();
    }
}

void ThreadPool::drain(std::unique_lock<std::mutex>& lock) {
    Job* job = job_;
    while (job->next < job->count) {
        const std::size_t item = job->next++;
        if (job->error) {
            continue;
        }
        job->active++;
        lock.unlock();
        std::exception_ptr error;
        try {
            (*job->task)(item);
        } catch (...) {
            error = std::current_exception();
        }
        lock.lock();
        job->active--;
        if (error && !job->error) {
            job->error = error;
        }
    }
}

void ThreadPool::work() {
    std::unique_lock<std::mutex> lock(mutex_);
    while (true) {
        wake_.wait(lock, [this] { return stop_ || (job_ != nullptr && job_->next < job_->count); });
        if (stop_) {
            return;
        }
        drain(lock);
        done_.notify_all();
    }
}

void ThreadPool::run(std::size_t count, const std::function<void(std::size_t)>& task) {
    if (count == 0) {
        return;
    }
    if (workers_.empty() || count == 1) {
        for (std::size_t i = 0; i < count; ++i) {
            task(i);
        }
        return;
    }
    Job job;
    job.task = &task;
    job.count = count;
    std::unique_lock<std::mutex> lock(mutex_);
    job_ = &job;
    wake_.notify_all();
    drain(lock);
    done_.wait(lock, [&] { return job.next >= job.count && job.active == 0; });
    job_ = nullptr;
    if (job.error) {
        std::rethrow_exception(job.error);
    }
}

}  // namespace coolercpp
