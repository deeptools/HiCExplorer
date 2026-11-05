// See hicx/reduce_matrix_entries.hpp.

#include "hicx/reduce_matrix_entries.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>
#include <utility>

namespace hicx::reduce_entries {

Entries upper_triangle(const Entries& entries) {
    Entries out;
    out.rows = entries.rows;
    out.cols = entries.cols;
    for (std::size_t k = 0; k < entries.data.size(); ++k) {
        if (entries.row[k] <= entries.col[k]) {
            out.row.push_back(entries.row[k]);
            out.col.push_back(entries.col[k]);
            out.data.push_back(entries.data[k]);
        }
    }
    return out;
}

std::size_t count_nan(const Entries& entries) {
    return static_cast<std::size_t>(
        std::count_if(entries.data.begin(), entries.data.end(), [](double v) { return std::isnan(v); }));
}

namespace {

using Key = std::pair<std::int64_t, std::int64_t>;

struct Operands {
    bool has_a = false;
    bool has_b = false;
    double a = 0.0;
    double b = 0.0;
};

// scipy's csr_binop_csr: op(a, b), op(a, 0) or op(0, b) per cell, cells whose
// result is zero dropped.
std::map<Key, double> binop(const std::map<Key, Operands>& cells,
                            const std::function<double(double, double)>& op) {
    std::map<Key, double> out;
    for (const auto& [key, operands] : cells) {
        const double value = op(operands.has_a ? operands.a : 0.0, operands.has_b ? operands.b : 0.0);
        if (value != 0.0 || std::isnan(value)) {
            out.emplace(key, value);
        }
    }
    return out;
}

}  // namespace

Cells reduce(const Entries& entries, const std::vector<std::int64_t>& group_of, std::int64_t groups,
             bool symmetrise, bool diagonal, const Arithmetic& arithmetic) {
    // np.bincount over the unique (row, col) of the mapped entries, summed in
    // entry order, then the cells sorted as np.unique sorts them
    std::unordered_map<std::uint64_t, std::size_t> slot_of;
    std::vector<Key> keys;
    std::vector<double> sums;
    for (std::size_t k = 0; k < entries.data.size(); ++k) {
        const std::int64_t r = group_of[static_cast<std::size_t>(entries.row[k])];
        const std::int64_t c = group_of[static_cast<std::size_t>(entries.col[k])];
        if (r < 0 || c < 0) {
            continue;
        }
        const std::uint64_t code = static_cast<std::uint64_t>(r) * static_cast<std::uint64_t>(groups) +
                                   static_cast<std::uint64_t>(c);
        auto [it, inserted] = slot_of.emplace(code, sums.size());
        if (inserted) {
            keys.emplace_back(r, c);
            sums.push_back(0.0);
        }
        sums[it->second] += entries.data[k];
    }
    std::vector<std::size_t> order(keys.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return keys[a] < keys[b]; });

    Cells out;
    out.size = groups;
    if (!symmetrise) {
        for (const std::size_t i : order) {
            const double value = arithmetic.cast(sums[i]);
            if (value != 0.0 || std::isnan(value)) {
                out.row.push_back(keys[i].first);
                out.col.push_back(keys[i].second);
                out.data.push_back(value);
            }
        }
        return out;
    }

    std::map<Key, Operands> sum_cells;
    std::map<Key, double> diag;
    for (const std::size_t i : order) {
        const double value = arithmetic.cast(sums[i]);
        const Key key = keys[i];
        Operands& a = sum_cells[key];
        a.has_a = true;
        a.a = value;
        Operands& b = sum_cells[{key.second, key.first}];
        b.has_b = true;
        b.b = value;
        if (key.first == key.second) {
            const double d = diagonal ? value : arithmetic.twice(value);
            if (d != 0.0 || std::isnan(d)) {
                diag[key] = d;
            }
        }
    }
    // result + result.T keeps the entries that are zero in the CSR of
    // `result` out of the sum only when the sum is zero as well
    const std::map<Key, double> symmetric = binop(sum_cells, arithmetic.add);
    std::map<Key, Operands> difference;
    for (const auto& [key, value] : symmetric) {
        difference[key].has_a = true;
        difference[key].a = value;
    }
    for (const auto& [key, value] : diag) {
        difference[key].has_b = true;
        difference[key].b = value;
    }
    for (const auto& [key, value] : binop(difference, arithmetic.subtract)) {
        out.row.push_back(key.first);
        out.col.push_back(key.second);
        out.data.push_back(value);
    }
    return out;
}

}  // namespace hicx::reduce_entries
