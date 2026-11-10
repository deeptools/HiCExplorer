// A thread pool whose results do not depend on the number of threads.
//
// cpp/OPTIMIZATION.md section 3 is the specification: two runs of the same
// binary on the same input must be byte-identical, and --threads 1 must be
// byte-identical to --threads 16. The way to get that is not to be careful, it
// is to make the arithmetic independent of scheduling in the first place:
//
//   * work is cut into a *fixed* number of index ranges that depends only on
//     the problem size and the requested thread count is used only to decide
//     how many ranges are in flight at once,
//   * every range is reduced sequentially inside one thread,
//   * range results are written into a preallocated slot indexed by the range,
//     never accumulated into a shared value, and are combined afterwards in
//     index order by the caller.
//
// Nothing here uses an atomic float, a work stealing queue, or a completion
// order. The only shared mutable state is the next-range counter, and which
// thread takes which range cannot change any result.
//
// The Python this replaces uses multiprocessing.Pool (hicFindTADs.py:1107),
// which forks a copy of the matrix per worker and serialises every result
// through a pickle. Threads share the one matrix, which is cpp/PLAN.md 4.4
// rule 8.

#ifndef HICX_PARALLEL_HPP
#define HICX_PARALLEL_HPP

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace hicx {

// The number of hardware threads, or 1 when that cannot be determined.
[[nodiscard]] inline unsigned int hardware_threads() {
    const unsigned int count = std::thread::hardware_concurrency();
    return count == 0 ? 1U : count;
}

// Runs body(index) for index in [0, count), on at most `threads` threads.
//
// The bodies must be independent: each writes only to its own slot of a
// result container the caller allocated up front. An exception thrown by a
// body is rethrown from parallel_for after every started body has finished, so
// a failure cannot leave a half-written result behind.
template <class Body>
void parallel_for(std::size_t count, unsigned int threads, Body&& body) {
    if (count == 0) {
        return;
    }
    const unsigned int workers = std::max(
        1U, static_cast<unsigned int>(
                std::min<std::size_t>(count, threads == 0 ? 1U : threads)));
    if (workers == 1) {
        for (std::size_t index = 0; index < count; ++index) {
            body(index);
        }
        return;
    }

    std::atomic<std::size_t> next{0};
    std::mutex error_mutex;
    std::exception_ptr error;

    const auto run = [&]() {
        while (true) {
            const std::size_t index = next.fetch_add(1, std::memory_order_relaxed);
            if (index >= count) {
                return;
            }
            try {
                body(index);
            } catch (...) {
                const std::lock_guard<std::mutex> guard(error_mutex);
                if (!error) {
                    error = std::current_exception();
                }
                // Drain the remaining indices so the other threads stop soon,
                // without leaving any half-finished slot behind.
                next.store(count, std::memory_order_relaxed);
                return;
            }
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(workers - 1);
    for (unsigned int i = 0; i + 1 < workers; ++i) {
        pool.emplace_back(run);
    }
    run();
    for (std::thread& worker : pool) {
        worker.join();
    }
    if (error) {
        std::rethrow_exception(error);
    }
}

// Splits [0, count) into `parts` contiguous ranges exactly as
// numpy.array_split does: the first count % parts ranges are one element
// longer than the rest. The split depends only on count and parts, never on
// the thread count, which is what makes a partitioned reduction thread count
// invariant.
struct IndexRange {
    std::size_t begin = 0;
    std::size_t end = 0;
    [[nodiscard]] std::size_t size() const noexcept { return end - begin; }
};

[[nodiscard]] inline std::vector<IndexRange> array_split(std::size_t count,
                                                         std::size_t parts) {
    std::vector<IndexRange> ranges;
    if (parts == 0) {
        return ranges;
    }
    ranges.reserve(parts);
    const std::size_t base = count / parts;
    const std::size_t remainder = count % parts;
    std::size_t begin = 0;
    for (std::size_t i = 0; i < parts; ++i) {
        const std::size_t length = base + (i < remainder ? 1 : 0);
        ranges.push_back(IndexRange{begin, begin + length});
        begin += length;
    }
    return ranges;
}

}  // namespace hicx

#endif  // HICX_PARALLEL_HPP
