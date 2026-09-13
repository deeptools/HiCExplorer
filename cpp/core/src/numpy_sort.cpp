#include "hicx/numpy_sort.hpp"

#include <array>
#include <cstddef>
#include <numeric>
#include <utility>

namespace hicx::npy {

namespace detail {

namespace {

// npysort_common.h
inline bool double_less(double a, double b) { return a < b || (b != b && a == a); }
inline bool longlong_less(std::int64_t a, std::int64_t b) { return a < b; }

// npy_sort.h.src
int get_msb(std::uint64_t unum) {
    int depth_limit = 0;
    while (unum >>= 1) {
        ++depth_limit;
    }
    return depth_limit;
}

constexpr int kSmallQuicksort = 15;       // SMALL_QUICKSORT
constexpr int kStack = 64 * 2;            // PYA_QS_STACK = NPY_BITSOF_INTP * 2

// aheapsort_<Tag> (npysort/heapsort.hpp): a[] is tosort offset by one.
template <class T, class Less>
void aheapsort(const T* v, std::int64_t* tosort, std::int64_t n, Less less) {
    std::int64_t* a = tosort - 1;
    std::int64_t i = 0;
    std::int64_t j = 0;
    std::int64_t l = 0;
    std::int64_t tmp = 0;
    for (l = n >> 1; l > 0; --l) {
        tmp = a[l];
        for (i = l, j = l << 1; j <= n;) {
            if (j < n && less(v[a[j]], v[a[j + 1]])) {
                ++j;
            }
            if (less(v[tmp], v[a[j]])) {
                a[i] = a[j];
                i = j;
                j += j;
            } else {
                break;
            }
        }
        a[i] = tmp;
    }
    for (; n > 1;) {
        tmp = a[n];
        a[n] = a[1];
        n -= 1;
        for (i = 1, j = 2; j <= n;) {
            if (j < n && less(v[a[j]], v[a[j + 1]])) {
                ++j;
            }
            if (less(v[tmp], v[a[j]])) {
                a[i] = a[j];
                i = j;
                j += j;
            } else {
                break;
            }
        }
        a[i] = tmp;
    }
}

// aquicksort_<Tag> (npysort/quicksort.cpp), statement for statement.
template <class T, class Less>
void aquicksort(const T* v, std::int64_t* tosort, std::int64_t num, Less less) {
    if (num <= 1) {
        return;
    }
    T vp;
    std::int64_t* pl = tosort;
    std::int64_t* pr = tosort + num - 1;
    std::array<std::int64_t*, kStack> stack{};
    std::int64_t** sptr = stack.data();
    std::int64_t* pm = nullptr;
    std::int64_t* pi = nullptr;
    std::int64_t* pj = nullptr;
    std::int64_t* pk = nullptr;
    std::int64_t vi = 0;
    std::array<int, kStack> depth{};
    int* psdepth = depth.data();
    int cdepth = get_msb(static_cast<std::uint64_t>(num)) * 2;

    for (;;) {
        if (cdepth < 0) {
            aheapsort(v, pl, pr - pl + 1, less);
            goto stack_pop;
        }
        while ((pr - pl) > kSmallQuicksort) {
            pm = pl + ((pr - pl) >> 1);
            if (less(v[*pm], v[*pl])) {
                std::swap(*pm, *pl);
            }
            if (less(v[*pr], v[*pm])) {
                std::swap(*pr, *pm);
            }
            if (less(v[*pm], v[*pl])) {
                std::swap(*pm, *pl);
            }
            vp = v[*pm];
            pi = pl;
            pj = pr - 1;
            std::swap(*pm, *pj);
            for (;;) {
                do {
                    ++pi;
                } while (less(v[*pi], vp));
                do {
                    --pj;
                } while (less(vp, v[*pj]));
                if (pi >= pj) {
                    break;
                }
                std::swap(*pi, *pj);
            }
            pk = pr - 1;
            std::swap(*pi, *pk);
            if (pi - pl < pr - pi) {
                *sptr++ = pi + 1;
                *sptr++ = pr;
                pr = pi - 1;
            } else {
                *sptr++ = pl;
                *sptr++ = pi - 1;
                pl = pi + 1;
            }
            *psdepth++ = --cdepth;
        }

        for (pi = pl + 1; pi <= pr; ++pi) {
            vi = *pi;
            vp = v[vi];
            pj = pi;
            pk = pi - 1;
            while (pj > pl && less(vp, v[*pk])) {
                *pj-- = *pk--;
            }
            *pj = vi;
        }
    stack_pop:
        if (sptr == stack.data()) {
            break;
        }
        pr = *(--sptr);
        pl = *(--sptr);
        cdepth = *(--psdepth);
    }
}

// ---------------------------------------------------------------------------
// introselect_<Tag, arg=true> (npysort/selection.cpp), without the pivot
// stack, which np.argpartition with a single kth starts empty and never reads
// back.

template <class T, class Less>
void median3_swap(const T* v, std::int64_t* tosort, std::int64_t low, std::int64_t mid,
                  std::int64_t high, Less less) {
    if (less(v[tosort[high]], v[tosort[mid]])) {
        std::swap(tosort[high], tosort[mid]);
    }
    if (less(v[tosort[high]], v[tosort[low]])) {
        std::swap(tosort[high], tosort[low]);
    }
    if (less(v[tosort[low]], v[tosort[mid]])) {
        std::swap(tosort[low], tosort[mid]);
    }
    std::swap(tosort[mid], tosort[low + 1]);
}

template <class T, class Less>
std::int64_t median5(const T* v, std::int64_t* tosort, Less less) {
    if (less(v[tosort[1]], v[tosort[0]])) {
        std::swap(tosort[1], tosort[0]);
    }
    if (less(v[tosort[4]], v[tosort[3]])) {
        std::swap(tosort[4], tosort[3]);
    }
    if (less(v[tosort[3]], v[tosort[0]])) {
        std::swap(tosort[3], tosort[0]);
    }
    if (less(v[tosort[4]], v[tosort[1]])) {
        std::swap(tosort[4], tosort[1]);
    }
    if (less(v[tosort[2]], v[tosort[1]])) {
        std::swap(tosort[2], tosort[1]);
    }
    if (less(v[tosort[3]], v[tosort[2]])) {
        if (less(v[tosort[3]], v[tosort[1]])) {
            return 1;
        }
        return 3;
    }
    return 2;
}

template <class T, class Less>
void unguarded_partition(const T* v, std::int64_t* tosort, const T pivot, std::int64_t* ll,
                         std::int64_t* hh, Less less) {
    for (;;) {
        do {
            (*ll)++;
        } while (less(v[tosort[*ll]], pivot));
        do {
            (*hh)--;
        } while (less(pivot, v[tosort[*hh]]));
        if (*hh < *ll) {
            break;
        }
        std::swap(tosort[*ll], tosort[*hh]);
    }
}

template <class T, class Less>
void dumb_select(const T* v, std::int64_t* tosort, std::int64_t num, std::int64_t kth,
                 Less less) {
    for (std::int64_t i = 0; i <= kth; i++) {
        std::int64_t minidx = i;
        T minval = v[tosort[i]];
        for (std::int64_t k = i + 1; k < num; k++) {
            if (less(v[tosort[k]], minval)) {
                minidx = k;
                minval = v[tosort[k]];
            }
        }
        std::swap(tosort[i], tosort[minidx]);
    }
}

template <class T, class Less>
void introselect(const T* v, std::int64_t* tosort, std::int64_t num, std::int64_t kth,
                 bool inexact, Less less);

template <class T, class Less>
std::int64_t median_of_median5(const T* v, std::int64_t* tosort, const std::int64_t num,
                               bool inexact, Less less) {
    std::int64_t i = 0;
    std::int64_t subleft = 0;
    std::int64_t right = num - 1;
    std::int64_t nmed = (right + 1) / 5;
    for (i = 0, subleft = 0; i < nmed; i++, subleft += 5) {
        std::int64_t m = median5(v, tosort + subleft, less);
        std::swap(tosort[subleft + m], tosort[i]);
    }
    if (nmed > 2) {
        introselect(v, tosort, nmed, nmed / 2, inexact, less);
    }
    return nmed / 2;
}

template <class T, class Less>
void introselect(const T* v, std::int64_t* tosort, std::int64_t num, std::int64_t kth,
                 bool inexact, Less less) {
    std::int64_t low = 0;
    std::int64_t high = num - 1;

    if (kth - low < 3) {
        dumb_select(v, tosort + low, high - low + 1, kth - low, less);
        return;
    }
    if (inexact && kth == num - 1) {
        std::int64_t maxidx = low;
        T maxval = v[tosort[low]];
        for (std::int64_t k = low + 1; k < num; k++) {
            if (!less(v[tosort[k]], maxval)) {
                maxidx = k;
                maxval = v[tosort[k]];
            }
        }
        std::swap(tosort[kth], tosort[maxidx]);
        return;
    }

    int depth_limit = get_msb(static_cast<std::uint64_t>(num)) * 2;
    for (; low + 1 < high;) {
        std::int64_t ll = low + 1;
        std::int64_t hh = high;
        if (depth_limit > 0 || hh - ll < 5) {
            const std::int64_t mid = low + (high - low) / 2;
            median3_swap(v, tosort, low, mid, high, less);
        } else {
            std::int64_t mid = ll + median_of_median5(v, tosort + ll, hh - ll, inexact, less);
            std::swap(tosort[mid], tosort[low]);
            ll--;
            hh++;
        }
        depth_limit--;
        unguarded_partition(v, tosort, v[tosort[low]], &ll, &hh, less);
        std::swap(tosort[low], tosort[hh]);
        if (hh >= kth) {
            high = hh - 1;
        }
        if (hh <= kth) {
            low = ll;
        }
    }
    if (high == low + 1) {
        if (less(v[tosort[high]], v[tosort[low]])) {
            std::swap(tosort[high], tosort[low]);
        }
    }
}

bool detect_avx512_skx() {
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512cd") &&
           __builtin_cpu_supports("avx512bw") && __builtin_cpu_supports("avx512dq") &&
           __builtin_cpu_supports("avx512vl");
}

}  // namespace

bool avx512_skx_available() {
    static const bool available = detect_avx512_skx();
    return available;
}

void argsort_introsort(const double* values, std::int64_t* index, std::int64_t n) {
    aquicksort(values, index, n, double_less);
}

void argsort_introsort(const std::int64_t* values, std::int64_t* index, std::int64_t n) {
    aquicksort(values, index, n, longlong_less);
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
    const auto n = static_cast<std::int64_t>(values.size());
    if (detail::avx512_skx_available()) {
        detail::argsort_avx512_skx(values.data(), index.data(), n);
    } else {
        detail::argsort_introsort(values.data(), index.data(), n);
    }
    return index;
}

std::vector<std::int64_t> argsort(const std::vector<std::int64_t>& values) {
    std::vector<std::int64_t> index = identity(values.size());
    const auto n = static_cast<std::int64_t>(values.size());
    if (detail::avx512_skx_available()) {
        detail::argsort_avx512_skx(values.data(), index.data(), n);
    } else {
        detail::argsort_introsort(values.data(), index.data(), n);
    }
    return index;
}

std::vector<std::int64_t> argpartition(const std::vector<double>& values, std::int64_t kth) {
    const auto n = static_cast<std::int64_t>(values.size());
    if (kth < 0) {
        kth += n;
    }
    std::vector<std::int64_t> index = identity(values.size());
    if (n > 0) {
        detail::introselect(values.data(), index.data(), n, kth, true, detail::double_less);
    }
    return index;
}

const char* argsort_dispatch_name() {
    return detail::avx512_skx_available() ? "avx512_skx" : "introsort";
}

}  // namespace hicx::npy
