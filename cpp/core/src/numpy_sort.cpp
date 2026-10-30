#include "hicx/numpy_sort.hpp"

#include <cstddef>
#include <numeric>
#include <utility>

// Index sorting and index selection with a fixed, documented tie order.
//
// argsort is an introspective quicksort (D. R. Musser, "Introspective sorting
// and selection algorithms", Software: Practice and Experience 27(8), 1997):
// quicksort with a median-of-three pivot (R. Sedgewick, "Implementing
// Quicksort programs", CACM 21(10), 1978), insertion sort for short ranges and
// heapsort (J. W. J. Williams, 1964) for ranges that were reached after too
// many partitioning steps. The choices that decide where equal keys end up
// are the ones numpy documents for kind='quicksort' (pivot rule, range length
// below which insertion sort takes over, depth limit of twice the binary
// logarithm of the length). They are reproduced here so that the permutation
// is the one numpy.argsort(kind='quicksort') returns on a CPU without the
// AVX-512 sorting path.
//
// argpartition is an introspective selection: Hoare's quickselect (1961) with
// a median-of-three pivot, a fall back to the median-of-medians pivot of groups
// of five (Blum, Floyd, Pratt, Rivest, Tarjan 1973) when partitioning does not
// shrink the range quickly enough, and selection sort for the first three
// ranks.
//
// Identity with numpy 1.26.4 was checked on randomised inputs, see
// cpp/tests/test_numpy_oracle.cpp and cpp/scripts/gen_numpy_oracle_fixture.py.

namespace hicx::npy {

namespace detail {

namespace {

using Index = std::int64_t;

// Order of float64 keys with NaN sorted after every number.
struct DoubleOrder {
    static bool before(double a, double b) { return a < b || (b != b && a == a); }
};

struct IntOrder {
    static bool before(std::int64_t a, std::int64_t b) { return a < b; }
};

int floor_log2(std::uint64_t x) {
    int result = 0;
    while (x >>= 1) {
        ++result;
    }
    return result;
}

// Ranges of at most this many entries beyond the first are insertion sorted.
constexpr Index kInsertionSpan = 15;

// Keys addressed through an index array.
template <class Key, class Order>
struct Indexed {
    const Key* key;
    Index* idx;
    bool before(Index i, Index j) const { return Order::before(key[idx[i]], key[idx[j]]); }
    bool before_value(const Key& k, Index j) const { return Order::before(k, key[idx[j]]); }
    bool value_before(Index i, const Key& k) const { return Order::before(key[idx[i]], k); }
    void swap(Index i, Index j) const { std::swap(idx[i], idx[j]); }
};

// Heapsort of idx[first .. first+count-1] as a binary max-heap: build with
// Floyd's bottom-up construction, then repeatedly move the maximum to the end.
// Of two children the right one is taken only when it is strictly larger.
template <class Key, class Order>
void heap_sort(const Key* key, Index* first, Index count) {
    auto sift_down = [&](Index root, Index size, Index moving) {
        Index hole = root;
        for (Index child = 2 * hole + 1; child < size; child = 2 * hole + 1) {
            if (child + 1 < size && Order::before(key[first[child]], key[first[child + 1]])) {
                ++child;
            }
            if (!Order::before(key[moving], key[first[child]])) {
                break;
            }
            first[hole] = first[child];
            hole = child;
        }
        first[hole] = moving;
    };
    for (Index root = count / 2 - 1; root >= 0; --root) {
        sift_down(root, count, first[root]);
    }
    for (Index size = count; size > 1;) {
        const Index last = first[size - 1];
        first[size - 1] = first[0];
        --size;
        sift_down(0, size, last);
    }
}

// Sorts idx[lo..hi]. `budget` is the number of partitioning steps still
// allowed on the path that led to this range. It is tested only when a range
// is entered, never while its own partitioning loop continues on the smaller
// part, which is what decides which ranges fall back to heapsort.
template <class Key, class Order>
void introsort_range(const Key* key, Index* idx, Index lo, Index hi, int budget) {
    if (budget < 0) {
        heap_sort<Key, Order>(key, idx + lo, hi - lo + 1);
        return;
    }
    const Indexed<Key, Order> a{key, idx};
    while (hi - lo > kInsertionSpan) {
        // Median of first, middle and last entry, left in order lo <= mid <= hi.
        const Index mid = lo + ((hi - lo) >> 1);
        if (a.before(mid, lo)) {
            a.swap(mid, lo);
        }
        if (a.before(hi, mid)) {
            a.swap(hi, mid);
        }
        if (a.before(mid, lo)) {
            a.swap(mid, lo);
        }
        const Key pivot = key[idx[mid]];
        // Park the pivot next to the last entry and partition lo+1 .. hi-2.
        a.swap(mid, hi - 1);
        Index i = lo;
        Index j = hi - 1;
        for (;;) {
            do {
                ++i;
            } while (a.value_before(i, pivot));
            do {
                --j;
            } while (a.before_value(pivot, j));
            if (i >= j) {
                break;
            }
            a.swap(i, j);
        }
        a.swap(i, hi - 1);
        --budget;
        // The larger part is sorted by a nested call, the smaller one here.
        if (i - lo < hi - i) {
            introsort_range<Key, Order>(key, idx, i + 1, hi, budget);
            hi = i - 1;
        } else {
            introsort_range<Key, Order>(key, idx, lo, i - 1, budget);
            lo = i + 1;
        }
    }
    for (Index i = lo + 1; i <= hi; ++i) {
        const Index moving = idx[i];
        Index j = i;
        while (j > lo && Order::before(key[moving], key[idx[j - 1]])) {
            idx[j] = idx[j - 1];
            --j;
        }
        idx[j] = moving;
    }
}

template <class Key, class Order>
void introsort_all(const Key* key, Index* idx, Index n) {
    if (n <= 1) {
        return;
    }
    introsort_range<Key, Order>(key, idx, 0, n - 1, floor_log2(static_cast<std::uint64_t>(n)) * 2);
}

// ---------------------------------------------------------------------------
// Selection.

// Median of three for the selection: the median of the entries at lo, mid and
// hi is left at lo.
template <class Key, class Order>
void select_order3(const Indexed<Key, Order>& v, Index lo, Index mid, Index hi) {
    if (v.before(hi, mid)) {
        v.swap(hi, mid);
    }
    if (v.before(hi, lo)) {
        v.swap(hi, lo);
    }
    if (v.before(lo, mid)) {
        v.swap(lo, mid);
    }
    // The median is now at lo and serves as the pivot. The entry at mid goes
    // to lo + 1, so that the partition below starts behind the pivot.
    v.swap(mid, lo + 1);
}

// Position (0..4) of the median of the five entries starting at p, found with
// a fixed comparison network; the entries are partially reordered.
template <class Key, class Order>
Index median_of_five(const Indexed<Key, Order>& v, Index p) {
    auto cx = [&](Index x, Index y) {
        if (v.before(p + y, p + x)) {
            v.swap(p + y, p + x);
        }
    };
    cx(0, 1);
    cx(3, 4);
    cx(0, 3);
    cx(1, 4);
    cx(1, 2);
    if (v.before(p + 3, p + 2)) {
        return v.before(p + 3, p + 1) ? 1 : 3;
    }
    return 2;
}

// After the call, positions 0..k have the k+1 smallest entries of idx[0..n-1]
// in order, by repeated minimum search.
template <class Key, class Order>
void select_by_minima(const Indexed<Key, Order>& v, Index base, Index n, Index k) {
    for (Index i = 0; i <= k; ++i) {
        Index best = i;
        Key best_key = v.key[v.idx[base + i]];
        for (Index j = i + 1; j < n; ++j) {
            if (Order::before(v.key[v.idx[base + j]], best_key)) {
                best = j;
                best_key = v.key[v.idx[base + j]];
            }
        }
        v.swap(base + i, base + best);
    }
}

template <class Key, class Order>
void select_kth(const Key* key, Index* idx, Index n, Index k);

// Pivot by median of medians: the medians of consecutive groups of five are
// gathered at the front and their median is selected there. Returns the
// position of that pivot relative to p.
template <class Key, class Order>
Index pivot_median_of_medians(const Key* key, Index* idx, Index p, Index count) {
    const Indexed<Key, Order> v{key, idx};
    const Index groups = count / 5;
    for (Index g = 0; g < groups; ++g) {
        const Index m = median_of_five<Key, Order>(v, p + 5 * g);
        v.swap(p + 5 * g + m, p + g);
    }
    if (groups > 2) {
        select_kth<Key, Order>(key, idx + p, groups, groups / 2);
    }
    return groups / 2;
}

template <class Key, class Order>
void select_kth(const Key* key, Index* idx, Index n, Index k) {
    const Indexed<Key, Order> v{key, idx};
    Index low = 0;
    Index high = n - 1;

    if (k < 3) {
        select_by_minima<Key, Order>(v, 0, n, k);
        return;
    }
    if (k == n - 1) {
        // Largest key; among equal ones the last is taken.
        Index best = 0;
        Key best_key = key[idx[0]];
        for (Index j = 1; j < n; ++j) {
            if (!Order::before(key[idx[j]], best_key)) {
                best = j;
                best_key = key[idx[j]];
            }
        }
        v.swap(k, best);
        return;
    }

    int budget = floor_log2(static_cast<std::uint64_t>(n)) * 2;
    while (low + 1 < high) {
        Index i = low + 1;
        Index j = high;
        if (budget > 0 || j - i < 5) {
            select_order3<Key, Order>(v, low, low + (high - low) / 2, high);
        } else {
            const Index pivot_pos = i + pivot_median_of_medians<Key, Order>(key, idx, i, j - i);
            v.swap(pivot_pos, low);
            --i;
            ++j;
        }
        --budget;
        const Key pivot = key[idx[low]];
        for (;;) {
            do {
                ++i;
            } while (v.value_before(i, pivot));
            do {
                --j;
            } while (v.before_value(pivot, j));
            if (j < i) {
                break;
            }
            v.swap(i, j);
        }
        v.swap(low, j);
        if (j >= k) {
            high = j - 1;
        }
        if (j <= k) {
            low = i;
        }
    }
    if (high == low + 1 && v.before(high, low)) {
        v.swap(high, low);
    }
}

}  // namespace

void argsort_introsort(const double* values, std::int64_t* index, std::int64_t n) {
    introsort_all<double, DoubleOrder>(values, index, n);
}

void argsort_introsort(const std::int64_t* values, std::int64_t* index, std::int64_t n) {
    introsort_all<std::int64_t, IntOrder>(values, index, n);
}

void argpartition_select(const double* values, std::int64_t* index, std::int64_t n,
                         std::int64_t kth) {
    select_kth<double, DoubleOrder>(values, index, n, kth);
}

}  // namespace detail

namespace {

std::vector<std::int64_t> identity(std::size_t n) {
    std::vector<std::int64_t> index(n);
    std::iota(index.begin(), index.end(), std::int64_t{0});
    return index;
}

}  // namespace

std::vector<std::int64_t> argsort(const std::vector<double>& values) {
    std::vector<std::int64_t> index = identity(values.size());
    detail::argsort_introsort(values.data(), index.data(), static_cast<std::int64_t>(values.size()));
    return index;
}

std::vector<std::int64_t> argsort(const std::vector<std::int64_t>& values) {
    std::vector<std::int64_t> index = identity(values.size());
    detail::argsort_introsort(values.data(), index.data(), static_cast<std::int64_t>(values.size()));
    return index;
}

std::vector<std::int64_t> argpartition(const std::vector<double>& values, std::int64_t kth) {
    const auto n = static_cast<std::int64_t>(values.size());
    if (kth < 0) {
        kth += n;
    }
    std::vector<std::int64_t> index = identity(values.size());
    if (n > 0) {
        detail::argpartition_select(values.data(), index.data(), n, kth);
    }
    return index;
}

}  // namespace hicx::npy
