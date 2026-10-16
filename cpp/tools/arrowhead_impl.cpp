#include "arrowhead_impl.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>

namespace hicx::arrowhead {

namespace {

double sign_of(double v) { return v > 0.0 ? 1.0 : (v < 0.0 ? -1.0 : 0.0); }

// DynamicProgrammingUtils.right(matrix, maxSize=n): right[i][j] = sum of
// matrix[k][j] for k = i..j (the column's running sum from the diagonal
// upward), defined for i <= j.
Matrix dp_right(const Matrix& m) {
    const int n = m.n;
    Matrix out(n);
    for (int k = 0; k < n; ++k) out.at(k, k) = m.at(k, k);
    for (int j = 1; j < n; ++j) {
        for (int i = j - 1; i >= 0; --i) {
            out.at(i, j) = m.at(i, j) + out.at(i + 1, j);
        }
    }
    return out;
}

// DynamicProgrammingUtils.upper(matrix, maxSize=n): upper[i][j] = sum of
// matrix[i][k] for k = i..j (the row's running sum away from the
// diagonal), defined for i <= j.
Matrix dp_upper(const Matrix& m) {
    const int n = m.n;
    Matrix out(n);
    for (int k = 0; k < n; ++k) out.at(k, k) = m.at(k, k);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            out.at(i, j) = m.at(i, j) + out.at(i, j - 1);
        }
    }
    return out;
}

// MatrixTriangles's up/upSign/upSquared/lo/loSign/loSquared, upVar, loVar
// and the final thresholded blockScore, all bundled together since nothing
// downstream of block_results needs the intermediate matrices individually.
struct Triangles {
    Matrix up_var;
    Matrix lo_var;
    Matrix up_sign;
    Matrix lo_sign;
    Matrix block_score;
};

Triangles compute_triangles(Matrix matrix) {
    const int n = matrix.n;
    for (double& v : matrix.v) {
        if (std::isnan(v)) v = 0.0;  // MatrixTools.setNaNs(matrix, 0)
    }

    Matrix squared(n);
    Matrix signed_m(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double v = matrix.at(i, j);
            squared.at(i, j) = v * v;
            signed_m.at(i, j) = sign_of(v);
        }
    }
    Matrix ones(n);
    std::fill(ones.v.begin(), ones.v.end(), 1.0);

    const Matrix r_sum = dp_right(matrix);
    const Matrix r_sign = dp_right(signed_m);
    const Matrix r_squared = dp_right(squared);
    const Matrix r_count = dp_right(ones);

    const Matrix u_sum = dp_upper(matrix);
    const Matrix u_sign = dp_upper(signed_m);
    const Matrix u_squared = dp_upper(squared);
    const Matrix u_count = dp_upper(ones);

    Matrix up(n), up_sign(n), up_squared(n), up_count(n);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const int bottom = (j - i + 1) / 2;
            up.at(i, j) = up.at(i, j - 1) + r_sum.at(i, j) - r_sum.at(i + bottom, j);
            up_sign.at(i, j) = up_sign.at(i, j - 1) + r_sign.at(i, j) - r_sign.at(i + bottom, j);
            up_squared.at(i, j) = up_squared.at(i, j - 1) + r_squared.at(i, j) - r_squared.at(i + bottom, j);
            up_count.at(i, j) = up_count.at(i, j - 1) + r_count.at(i, j) - r_count.at(i + bottom, j);
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double c = up_count.at(i, j) == 0.0 ? 1.0 : up_count.at(i, j);
            up.at(i, j) /= c;
            up_sign.at(i, j) /= c;
            up_squared.at(i, j) /= c;
        }
    }

    Matrix lo(n), lo_sign(n), lo_squared(n), lo_count(n);
    for (int a = 0; a < n; ++a) {
        for (int b = a + 1; b < n; ++b) {
            const int val = (b - a + 1) / 2;
            const int endpt = std::min(2 * b - a, n - 1);
            lo_count.at(a, b) = lo_count.at(a, b - 1) + u_count.at(b, endpt) - r_count.at(a + val, b);
            lo.at(a, b) = lo.at(a, b - 1) + u_sum.at(b, endpt) - r_sum.at(a + val, b);
            lo_sign.at(a, b) = lo_sign.at(a, b - 1) + u_sign.at(b, endpt) - r_sign.at(a + val, b);
            lo_squared.at(a, b) = lo_squared.at(a, b - 1) + u_squared.at(b, endpt) - r_squared.at(a + val, b);
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double c = lo_count.at(i, j) == 0.0 ? 1.0 : lo_count.at(i, j);
            lo.at(i, j) /= c;
            lo_sign.at(i, j) /= c;
            lo_squared.at(i, j) /= c;
        }
    }

    Triangles t;
    t.up_var = Matrix(n);
    t.lo_var = Matrix(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            t.up_var.at(i, j) = up_squared.at(i, j) - up.at(i, j) * up.at(i, j);
            t.lo_var.at(i, j) = lo_squared.at(i, j) - lo.at(i, j) * lo.at(i, j);
        }
    }
    t.up_sign = std::move(up_sign);
    t.lo_sign = std::move(lo_sign);

    auto normalize_by_max = [&](const Matrix& m) {
        double max_v = m.v.empty() ? 0.0 : m.v[0];
        for (const double v : m.v) max_v = std::max(max_v, v);
        Matrix out(n);
        for (int i = 0; i < n * n; ++i) out.v[static_cast<std::size_t>(i)] = m.v[static_cast<std::size_t>(i)] / max_v;
        return out;
    };
    Matrix diff(n), diff_sign(n), sum_var(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            diff.at(i, j) = lo.at(i, j) - up.at(i, j);
            diff_sign.at(i, j) = t.lo_sign.at(i, j) - t.up_sign.at(i, j);
            sum_var.at(i, j) = t.up_var.at(i, j) + t.lo_var.at(i, j);
        }
    }
    diff = normalize_by_max(diff);
    diff_sign = normalize_by_max(diff_sign);
    sum_var = normalize_by_max(sum_var);

    t.block_score = Matrix(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            t.block_score.at(i, j) = (diff.at(i, j) + diff_sign.at(i, j)) - sum_var.at(i, j);
        }
    }
    return t;
}

void apply_thresholds(Triangles& t, double var_threshold, double sign_threshold) {
    const int n = t.block_score.n;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if ((-t.up_sign.at(i, j)) < sign_threshold || t.lo_sign.at(i, j) < sign_threshold) {
                t.block_score.at(i, j) = 0.0;
            }
        }
    }
    if (!std::isnan(var_threshold)) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (t.up_var.at(i, j) + t.lo_var.at(i, j) > var_threshold) {
                    t.block_score.at(i, j) = 0.0;
                }
            }
        }
    }
}

// BinaryConnectedComponents.detection(image, threshold=0): 8-connectivity
// components of pixels with value > threshold. Any correct
// connected-components algorithm partitions the same pixels into the same
// sets, so this is a plain BFS flood fill rather than the real tool's
// two-pass union-find mechanics (see arrowhead_impl.hpp deviation 2).
std::vector<std::vector<std::pair<int, int>>> connected_components(const Matrix& score) {
    const int n = score.n;
    std::vector<char> visited(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0);
    std::vector<std::vector<std::pair<int, int>>> components;
    std::vector<std::pair<int, int>> stack;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const std::size_t idx = static_cast<std::size_t>(i) * static_cast<std::size_t>(n) + j;
            // BinaryConnectedComponents.detection's own condition is
            // `image[i][j] > threshold`, which Java's IEEE 754 comparison
            // makes false for NaN (so NaN pixels never seed or join a
            // component); `score <= 0.0` is also false for NaN, so the
            // skip test is the negation, not the mirror, of that condition.
            if (visited[idx] || !(score.at(i, j) > 0.0)) continue;
            std::vector<std::pair<int, int>> component;
            stack.clear();
            stack.emplace_back(i, j);
            visited[idx] = 1;
            while (!stack.empty()) {
                const auto [ci, cj] = stack.back();
                stack.pop_back();
                component.emplace_back(ci, cj);
                for (int di = -1; di <= 1; ++di) {
                    for (int dj = -1; dj <= 1; ++dj) {
                        if (di == 0 && dj == 0) continue;
                        const int ni = ci + di;
                        const int nj = cj + dj;
                        if (ni < 0 || ni >= n || nj < 0 || nj >= n) continue;
                        const std::size_t nidx = static_cast<std::size_t>(ni) * static_cast<std::size_t>(n) + nj;
                        if (visited[nidx] || !(score.at(ni, nj) > 0.0)) continue;
                        visited[nidx] = 1;
                        stack.emplace_back(ni, nj);
                    }
                }
            }
            components.push_back(std::move(component));
        }
    }
    return components;
}

}  // namespace

Matrix directionality_index_upstream(const Matrix& observed, int gap) {
    const int n = observed.n;
    Matrix d(n);
    for (int i = 0; i < n; ++i) {
        int window = std::min(n - (i + gap), i - gap);
        window = std::min(window, n);
        if (window < gap) continue;
        for (int off = 0; off <= window - gap; ++off) {
            const int left = i - gap - off;
            const int right = i + gap + off;
            const double a = observed.at(i, left);
            const double b = observed.at(i, right);
            d.at(i, right) = (a - b) / (a + b);
        }
    }
    return d;
}

std::vector<HighScore> block_results(const Matrix& observed, double var_threshold, double sign_threshold) {
    const Matrix d_upstream = directionality_index_upstream(observed, 7);
    Triangles t = compute_triangles(d_upstream);
    apply_thresholds(t, var_threshold, sign_threshold);

    const auto components = connected_components(t.block_score);
    std::vector<HighScore> results;
    results.reserve(components.size());
    for (const auto& component : components) {
        auto best = component.front();
        double best_score = t.block_score.at(best.first, best.second);
        for (const auto& [i, j] : component) {
            const double s = t.block_score.at(i, j);
            if (s > best_score) {
                best_score = s;
                best = {i, j};
            }
        }
        HighScore hs;
        hs.i = best.first;
        hs.j = best.second;
        hs.score = best_score;
        hs.u_var_score = t.up_var.at(best.first, best.second);
        hs.l_var_score = t.lo_var.at(best.first, best.second);
        hs.up_sign = -t.up_sign.at(best.first, best.second);
        hs.lo_sign = t.lo_sign.at(best.first, best.second);
        results.push_back(hs);
    }
    return results;
}

namespace {

// Exact-field equality, matching HighScore.equals/hashCode (every field,
// not just i/j), used by orderedSetDifference's HashSet-based set
// difference.
struct HighScoreKey {
    int i, j;
    double score, u_var, l_var, up_sign, lo_sign;

    friend bool operator<(const HighScoreKey& a, const HighScoreKey& b) {
        if (a.i != b.i) return a.i < b.i;
        if (a.j != b.j) return a.j < b.j;
        if (a.score != b.score) return a.score < b.score;
        if (a.u_var != b.u_var) return a.u_var < b.u_var;
        if (a.l_var != b.l_var) return a.l_var < b.l_var;
        if (a.up_sign != b.up_sign) return a.up_sign < b.up_sign;
        return a.lo_sign < b.lo_sign;
    }
};

HighScoreKey key_of(const HighScore& h) {
    return {h.i, h.j, h.score, h.u_var_score, h.l_var_score, h.up_sign, h.lo_sign};
}

// BlockBuster.orderedSetDifference: the (deduplicated) elements of listA
// that are not present (by full field equality) anywhere in listB.
std::vector<HighScore> ordered_set_difference(const std::vector<HighScore>& a, const std::vector<HighScore>& b) {
    std::set<HighScoreKey> set_a;
    for (const HighScore& h : a) set_a.insert(key_of(h));
    std::set<HighScoreKey> set_b;
    for (const HighScore& h : b) set_b.insert(key_of(h));

    std::vector<HighScore> diff;
    std::map<HighScoreKey, HighScore> by_key;
    for (const HighScore& h : a) by_key.emplace(key_of(h), h);
    for (const HighScoreKey& k : set_a) {
        if (set_b.find(k) == set_b.end()) {
            diff.push_back(by_key.at(k));
        }
    }
    return diff;
}

// BlockBuster.filterBlocksBySize: drop blocks whose width is <= minWidth.
std::vector<HighScore> filter_by_size(std::vector<HighScore> blocks, int min_width) {
    std::vector<HighScore> out;
    for (HighScore& h : blocks) {
        if (h.width() > min_width) out.push_back(h);
    }
    return out;
}

// BlockBuster.appendNonConflictingBlocks: a candidate is rejected if any
// bin in its [i, j] span is already an edge (i or j) of a block already in
// mainList; accepted candidates are appended and become edges themselves,
// so later candidates in `additions` see them too.
void append_non_conflicting(std::vector<HighScore>& main_list, const std::vector<HighScore>& additions) {
    std::map<int, HighScore> block_edges;
    for (const HighScore& h : main_list) {
        block_edges[h.i] = h;
        block_edges[h.j] = h;
    }
    for (const HighScore& h : additions) {
        const int lo = std::min(h.i, h.j);
        const int hi = std::max(h.i, h.j);
        bool conflicts = false;
        for (int k = lo; k <= hi; ++k) {
            if (block_edges.find(k) != block_edges.end()) {
                conflicts = true;
                break;
            }
        }
        if (!conflicts) {
            main_list.push_back(h);
            block_edges[h.i] = h;
            block_edges[h.j] = h;
        }
    }
}

// BinnedScore + BlockBuster.binScoresByDistance: merge scores whose (i, j)
// corners are within `dist` bins of an existing bin's corner range, each
// merged bin collapsing to the mean of its members' scores at the bin's own
// (max i, max j) corner.
std::vector<HighScore> bin_by_distance(const std::vector<HighScore>& scores, int dist) {
    struct Bin {
        int min_x, max_x, min_y, max_y;
        double sum_score = 0.0, sum_u = 0.0, sum_l = 0.0, sum_up = 0.0, sum_lo = 0.0;
        int count = 0;

        void add(const HighScore& h) {
            if (h.i < min_x) {
                min_x = h.i;
            } else if (h.i > max_x) {
                max_x = h.i;
            }
            if (h.j < min_y) {
                min_y = h.j;
            } else if (h.j > max_y) {
                max_y = h.j;
            }
            sum_score += h.score;
            sum_u += h.u_var_score;
            sum_l += h.l_var_score;
            sum_up += h.up_sign;
            sum_lo += h.lo_sign;
            ++count;
        }

        [[nodiscard]] bool is_near(const HighScore& h, int threshold) const {
            const bool near_x = std::abs(min_x - h.i) < threshold || std::abs(max_x - h.i) < threshold;
            const bool near_y = std::abs(min_y - h.j) < threshold || std::abs(max_y - h.j) < threshold;
            return near_x && near_y;
        }

        [[nodiscard]] HighScore to_high_score() const {
            HighScore h;
            h.i = max_x;
            h.j = max_y;
            h.score = sum_score / count;
            h.u_var_score = sum_u / count;
            h.l_var_score = sum_l / count;
            h.up_sign = sum_up / count;
            h.lo_sign = sum_lo / count;
            return h;
        }
    };

    std::vector<Bin> bins;
    for (const HighScore& h : scores) {
        bool placed = false;
        for (Bin& bin : bins) {
            if (bin.is_near(h, dist)) {
                bin.add(h);
                placed = true;
                break;
            }
        }
        if (!placed) {
            Bin bin;
            bin.min_x = bin.max_x = h.i;
            bin.min_y = bin.max_y = h.j;
            bin.add(h);
            bins.push_back(bin);
        }
    }
    std::vector<HighScore> out;
    out.reserve(bins.size());
    for (const Bin& bin : bins) out.push_back(bin.to_high_score());
    return out;
}

}  // namespace

std::vector<HighScore> call_domains(int chrom_length_bins, const Options& options,
                                    const std::function<Matrix(int, int)>& observed_window) {
    const int matrix_width = options.matrix_width;
    const int increment = matrix_width / 2;

    auto call_sub_blockbuster = [&](double var_threshold, double sign_threshold) {
        std::vector<HighScore> cumulative;
        for (int lim_start = 0; lim_start < chrom_length_bins; lim_start += increment) {
            int window_start = lim_start;
            int lim_end = std::min(lim_start + matrix_width, chrom_length_bins);
            if (lim_end == chrom_length_bins && chrom_length_bins > increment) {
                window_start = std::max(0, lim_end - matrix_width);
            }
            const Matrix observed = observed_window(window_start, lim_end);
            std::vector<HighScore> window_results = block_results(observed, var_threshold, sign_threshold);
            for (HighScore& h : window_results) {
                h.i += window_start;
                h.j += window_start;
            }
            cumulative.insert(cumulative.end(), window_results.begin(), window_results.end());
            // The real BlockBuster keeps incrementing limStart by `increment`
            // even after a window has already been clamped to the
            // chromosome's end, reprocessing the same clamped window until
            // limStart passes chrLength. That only produces duplicate
            // candidates (deduplicated downstream by the set-based merge
            // steps), so this stops as soon as a window reaches the end
            // instead.
            if (lim_end == chrom_length_bins) break;
        }
        return cumulative;
    };

    // Low-confidence sweep: BlockBuster.run's maxLowSignThreshold down to
    // minLowSignThreshold, stepping by decrementLowSignThreshold, stopping
    // at the first non-empty result (var_threshold NaN: no variance gate).
    std::vector<HighScore> results;
    for (double sign_threshold = options.max_low_sign_threshold;
         sign_threshold >= options.min_low_sign_threshold - 1e-12;
         sign_threshold -= options.decrement_low_sign_threshold) {
        results = call_sub_blockbuster(std::numeric_limits<double>::quiet_NaN(), sign_threshold);
        if (!results.empty()) break;
    }

    // High-confidence pass: fixed thresholds, fewer but more reliable blocks.
    std::vector<HighScore> high_confidence =
        call_sub_blockbuster(options.var_threshold, options.high_sign_threshold);

    std::vector<HighScore> unique_blocks = ordered_set_difference(results, high_confidence);
    std::vector<HighScore> filtered_unique = filter_by_size(unique_blocks, options.min_block_size);
    append_non_conflicting(high_confidence, filtered_unique);

    if (high_confidence.empty()) {
        return {};
    }

    std::vector<HighScore> binned = bin_by_distance(high_confidence, 5);
    binned = bin_by_distance(binned, 10);
    std::sort(binned.begin(), binned.end(),
             [](const HighScore& a, const HighScore& b) { return a.sort_value() > b.sort_value(); });
    return binned;
}

}  // namespace hicx::arrowhead
