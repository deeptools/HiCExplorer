// The iteration order of a Python 3 set, reproduced.
//
// Several HiCExplorer outputs are written by iterating a Python set (the
// filtered bins of hicCorrectMatrix, the candidate intervals that
// hicMergeLoops takes from intervaltree), so their line order is the slot
// order of the interpreter's hash table. The order is observable, hence
// reproduced; it is not part of any language guarantee.
//
// Behaviour implemented, from the interpreter's documented hashing rules and
// from black-box measurement:
//   * hash(int): the absolute value modulo the Mersenne prime 2^61 - 1 with
//     the sign restored, and -1 replaced by -2 (sys.hash_info documents the
//     modulus);
//   * hash(tuple): an xxHash64 style accumulation over the hashes of the
//     items, finished with the length (the scheme of CPython 3.8 and later);
//   * the table: a power of two number of slots, probe sequence of ten
//     consecutive slots (when they lie inside the table) followed by a jump
//     i -> 5 i + 1 + perturb with perturb = hash >> 5, 10, 15, ...;
//   * removal leaves a tombstone; an insertion that finds none of its key
//     reuses the tombstone it passed last;
//   * the table grows when live + tombstone slots reach three fifths of its
//     size, to the power of two above 4 x live (2 x live beyond 50000), and a
//     difference update that leaves more than a quarter of the slots as
//     tombstones rebuilds the table the same way;
//   * set.update(other_set) and set(other_set) take shortcuts that place the
//     keys differently from adding them one by one.
//
// Verified against Python 3.12.7 (the interpreter of the HiCExplorer 3.7.6
// reference environment): cpp/tests/test_cpython_oracles.cpp replays
// cpp/tests/fixtures/cpython/set_order.txt, written by
// cpp/scripts/gen_pyset_order_fixture.py: random add / discard / update /
// difference_update sequences on sets of ints (negative, up to the full 64 bit
// range) and of pairs and triples, plus bulk cases of up to 10^6 keys.
//
// A key is an opaque 64 bit identifier; two keys are the same element when
// the identifiers are equal. The caller supplies hash(key) (hash_of_int or
// hash_of_tuple below), because the table position depends on the hash only.

#ifndef HICX_PY_SET_ORDER_HPP
#define HICX_PY_SET_ORDER_HPP

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <utility>
#include <vector>

namespace hicx::pyorder {

// hash(n) for a Python int that fits 64 bits.
[[nodiscard]] inline std::int64_t hash_of_int(std::int64_t value) {
    constexpr std::uint64_t kModulus = (std::uint64_t{1} << 61) - 1;
    const bool negative = value < 0;
    // Negating in unsigned arithmetic is exact for INT64_MIN as well.
    const std::uint64_t magnitude = negative ? 0 - static_cast<std::uint64_t>(value)
                                             : static_cast<std::uint64_t>(value);
    const auto residue = static_cast<std::int64_t>(magnitude % kModulus);
    const std::int64_t signed_residue = negative ? -residue : residue;
    return signed_residue == -1 ? -2 : signed_residue;
}

// hash(t) for a tuple whose items hash to item_hashes.
[[nodiscard]] inline std::int64_t hash_of_tuple(std::initializer_list<std::int64_t> item_hashes) {
    // The multipliers are the primes of xxHash64 (public specification).
    constexpr std::uint64_t kPrimeA = 11400714785074694791ULL;
    constexpr std::uint64_t kPrimeB = 14029467366897019727ULL;
    constexpr std::uint64_t kPrimeC = 2870177450012600261ULL;
    std::uint64_t acc = kPrimeC;
    for (const std::int64_t item : item_hashes) {
        acc += static_cast<std::uint64_t>(item) * kPrimeB;
        acc = (acc << 31) | (acc >> 33);
        acc *= kPrimeA;
    }
    acc += static_cast<std::uint64_t>(item_hashes.size()) ^ (kPrimeC ^ 3527539ULL);
    // -1 is reserved by the interpreter as the error value.
    return acc == ~std::uint64_t{0} ? 1546275796 : static_cast<std::int64_t>(acc);
}

[[nodiscard]] inline std::int64_t hash_of_int_pair(std::int64_t first, std::int64_t second) {
    return hash_of_tuple({hash_of_int(first), hash_of_int(second)});
}

class OrderedSet {
public:
    OrderedSet() : slots_(kMinSlots) {}

    [[nodiscard]] std::size_t size() const { return live_; }
    [[nodiscard]] bool empty() const { return live_ == 0; }

    // Adds key (a no-op when present).
    void add(std::int64_t key, std::int64_t hash) {
        std::size_t reusable = kNone;
        std::size_t found = kNone;
        const std::size_t empty_slot = probe(key, hash, &reusable, &found);
        if (found != kNone) {
            return;
        }
        if (reusable != kNone) {
            fill_slot(reusable, key, hash);
            ++live_;
            return;
        }
        fill_slot(empty_slot, key, hash);
        ++live_;
        ++used_;
        if (used_ * 5 >= mask() * 3) {
            rebuild(grown_size());
        }
    }

    [[nodiscard]] bool contains(std::int64_t key, std::int64_t hash) const {
        std::size_t found = kNone;
        probe(key, hash, nullptr, &found);
        return found != kNone;
    }

    // Removes key; returns whether it was present.
    bool discard(std::int64_t key, std::int64_t hash) {
        std::size_t found = kNone;
        probe(key, hash, nullptr, &found);
        if (found == kNone) {
            return false;
        }
        slots_[found].state = State::Tombstone;
        --live_;
        return true;
    }

    // The keys in iteration order.
    [[nodiscard]] std::vector<std::int64_t> items() const {
        std::vector<std::int64_t> out;
        out.reserve(live_);
        for (const Slot& slot : slots_) {
            if (slot.state == State::Live) {
                out.push_back(slot.key);
            }
        }
        return out;
    }

    // self.update(other) for a set argument, and set(other) when self is
    // empty.
    void update(const OrderedSet& other) {
        if (&other == this || other.live_ == 0) {
            return;
        }
        if ((used_ + other.live_) * 5 >= mask() * 3) {
            rebuild(next_pow2_above((live_ + other.live_) * 2));
        }
        if (used_ == 0) {
            if (mask() == other.mask() && other.used_ == other.live_) {
                // Same size and no tombstones: the table is copied slot for slot.
                slots_ = other.slots_;
                used_ = other.used_;
                live_ = other.live_;
                return;
            }
            // Distinct keys into an empty table: no equality probing.
            used_ = live_ = other.live_;
            for (const Slot& slot : other.slots_) {
                if (slot.state == State::Live) {
                    place_new(slots_, slot.key, slot.hash);
                }
            }
            return;
        }
        for (const Slot& slot : other.slots_) {
            if (slot.state == State::Live) {
                add(slot.key, slot.hash);
            }
        }
    }

    // self.difference_update(other) for a set argument.
    void difference_update(const OrderedSet& other) {
        for (const Slot& slot : other.slots_) {
            if (slot.state == State::Live) {
                discard(slot.key, slot.hash);
            }
        }
        if (used_ - live_ > mask() / 4) {
            rebuild(grown_size());
        }
    }

private:
    enum class State : std::uint8_t { Free, Tombstone, Live };
    struct Slot {
        std::int64_t key = 0;
        std::int64_t hash = 0;
        State state = State::Free;
    };

    static constexpr std::size_t kMinSlots = 8;
    static constexpr std::size_t kRun = 10;  // consecutive slots tried per probe step
    static constexpr int kPerturbShift = 5;
    static constexpr std::size_t kNone = static_cast<std::size_t>(-1);

    [[nodiscard]] std::size_t mask() const { return slots_.size() - 1; }

    [[nodiscard]] static std::size_t next_pow2_above(std::size_t n) {
        std::size_t size = kMinSlots;
        while (size <= n) {
            size <<= 1;
        }
        return size;
    }
    [[nodiscard]] std::size_t grown_size() const {
        return next_pow2_above(live_ > 50000 ? live_ * 2 : live_ * 4);
    }

    // Visits the probe sequence of hash. Returns the first free slot reached;
    // *found is the slot holding key when it comes first, *reusable the last
    // tombstone passed before that.
    std::size_t probe(std::int64_t key, std::int64_t hash, std::size_t* reusable,
                      std::size_t* found) const {
        const std::size_t m = mask();
        std::size_t perturb = static_cast<std::size_t>(hash);
        std::size_t start = perturb & m;
        while (true) {
            // The run of ten is used only when it fits inside the table.
            const std::size_t run = (start + kRun - 1 <= m) ? kRun : 1;
            for (std::size_t slot = start; slot < start + run; ++slot) {
                const Slot& s = slots_[slot];
                if (s.state == State::Free) {
                    return slot;
                }
                if (s.state == State::Live && s.key == key) {
                    *found = slot;
                    return slot;
                }
                if (s.state == State::Tombstone && reusable != nullptr) {
                    *reusable = slot;
                }
            }
            perturb >>= kPerturbShift;
            start = (start * 5 + 1 + perturb) & m;
        }
    }

    void fill_slot(std::size_t slot, std::int64_t key, std::int64_t hash) {
        slots_[slot] = Slot{key, hash, State::Live};
    }

    // Insertion of a key known to be absent into a table without tombstones.
    static void place_new(std::vector<Slot>& table, std::int64_t key, std::int64_t hash) {
        const std::size_t m = table.size() - 1;
        std::size_t perturb = static_cast<std::size_t>(hash);
        std::size_t start = perturb & m;
        while (true) {
            const std::size_t run = (start + kRun - 1 <= m) ? kRun : 1;
            for (std::size_t slot = start; slot < start + run; ++slot) {
                if (table[slot].state == State::Free) {
                    table[slot] = Slot{key, hash, State::Live};
                    return;
                }
            }
            perturb >>= kPerturbShift;
            start = (start * 5 + 1 + perturb) & m;
        }
    }

    // Moves the live keys, in slot order, into a fresh table of new_size.
    void rebuild(std::size_t new_size) {
        std::vector<Slot> old = std::move(slots_);
        slots_.assign(new_size, Slot{});
        for (const Slot& slot : old) {
            if (slot.state == State::Live) {
                place_new(slots_, slot.key, slot.hash);
            }
        }
        used_ = live_;
    }

    std::vector<Slot> slots_;
    std::size_t used_ = 0;  // live slots plus tombstones
    std::size_t live_ = 0;
};

}  // namespace hicx::pyorder

#endif  // HICX_PY_SET_ORDER_HPP
