#include "stripes_impl.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace hicx::stripes {

namespace {
constexpr double kPi = 3.14159265358979323846;
}  // namespace

double quantile(std::vector<double> values, double q) {
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    if (values.size() == 1) {
        return values[0];
    }
    const double idx = q * static_cast<double>(values.size() - 1);
    const auto lo = static_cast<std::size_t>(std::floor(idx));
    const auto hi = static_cast<std::size_t>(std::ceil(idx));
    const double frac = idx - static_cast<double>(lo);
    return values[lo] + (values[hi] - values[lo]) * frac;
}

Image box_blur(const Image& img, int k) {
    Image out(img.rows, img.cols);
    const int half = k / 2;
    for (int r = 0; r < img.rows; ++r) {
        for (int c = 0; c < img.cols; ++c) {
            double sum = 0.0;
            for (int dr = -half; dr <= half; ++dr) {
                for (int dc = -half; dc <= half; ++dc) {
                    const int rr = std::clamp(r + dr, 0, img.rows - 1);
                    const int cc = std::clamp(c + dc, 0, img.cols - 1);
                    sum += img.at(rr, cc);
                }
            }
            out.at(r, c) = sum / static_cast<double>(k * k);
        }
    }
    return out;
}

namespace {

Image gaussian_blur(const Image& img, double sigma) {
    const int radius = std::max(1, static_cast<int>(std::lround(3.0 * sigma)));
    std::vector<double> kernel(static_cast<std::size_t>(2 * radius + 1));
    double sum = 0.0;
    for (int i = -radius; i <= radius; ++i) {
        const double value = std::exp(-(static_cast<double>(i) * i) / (2.0 * sigma * sigma));
        kernel[static_cast<std::size_t>(i + radius)] = value;
        sum += value;
    }
    for (double& v : kernel) {
        v /= sum;
    }
    // Separable: horizontal pass then vertical pass, edge-replicated border.
    Image temp(img.rows, img.cols);
    for (int r = 0; r < img.rows; ++r) {
        for (int c = 0; c < img.cols; ++c) {
            double acc = 0.0;
            for (int i = -radius; i <= radius; ++i) {
                const int cc = std::clamp(c + i, 0, img.cols - 1);
                acc += img.at(r, cc) * kernel[static_cast<std::size_t>(i + radius)];
            }
            temp.at(r, c) = acc;
        }
    }
    Image out(img.rows, img.cols);
    for (int r = 0; r < img.rows; ++r) {
        for (int c = 0; c < img.cols; ++c) {
            double acc = 0.0;
            for (int i = -radius; i <= radius; ++i) {
                const int rr = std::clamp(r + i, 0, img.rows - 1);
                acc += temp.at(rr, c) * kernel[static_cast<std::size_t>(i + radius)];
            }
            out.at(r, c) = acc;
        }
    }
    return out;
}

void sobel(const Image& img, Image& gx, Image& gy) {
    gx = Image(img.rows, img.cols);
    gy = Image(img.rows, img.cols);
    static constexpr int kSx[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    static constexpr int kSy[3][3] = {{1, 2, 1}, {0, 0, 0}, {-1, -2, -1}};
    for (int r = 0; r < img.rows; ++r) {
        for (int c = 0; c < img.cols; ++c) {
            double sx = 0.0;
            double sy = 0.0;
            for (int dr = -1; dr <= 1; ++dr) {
                for (int dc = -1; dc <= 1; ++dc) {
                    const int rr = std::clamp(r + dr, 0, img.rows - 1);
                    const int cc = std::clamp(c + dc, 0, img.cols - 1);
                    const double v = img.at(rr, cc);
                    sx += v * kSx[dr + 1][dc + 1];
                    sy += v * kSy[dr + 1][dc + 1];
                }
            }
            gx.at(r, c) = sx;
            gy.at(r, c) = sy;
        }
    }
}

}  // namespace

Image canny(const Image& gray01, double sigma) {
    const Image smoothed = gaussian_blur(gray01, sigma);
    Image gx;
    Image gy;
    sobel(smoothed, gx, gy);

    const int rows = gray01.rows;
    const int cols = gray01.cols;
    Image mag(rows, cols);
    Image dir(rows, cols);  // quantised to 0, 45, 90, 135 degrees
    double max_mag = 0.0;
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            const double x = gx.at(r, c);
            const double y = gy.at(r, c);
            const double m = std::sqrt(x * x + y * y);
            mag.at(r, c) = m;
            max_mag = std::max(max_mag, m);
            double angle = std::atan2(y, x) * 180.0 / kPi;
            if (angle < 0) {
                angle += 180.0;
            }
            double quantised = 0.0;
            if ((angle >= 0 && angle < 22.5) || (angle >= 157.5 && angle <= 180.0)) {
                quantised = 0.0;
            } else if (angle >= 22.5 && angle < 67.5) {
                quantised = 45.0;
            } else if (angle >= 67.5 && angle < 112.5) {
                quantised = 90.0;
            } else {
                quantised = 135.0;
            }
            dir.at(r, c) = quantised;
        }
    }
    if (max_mag <= 0.0) {
        return Image(rows, cols, 0.0);
    }

    // Non-maximum suppression along the gradient direction.
    Image nms(rows, cols);
    for (int r = 1; r < rows - 1; ++r) {
        for (int c = 1; c < cols - 1; ++c) {
            const double m = mag.at(r, c);
            double n1 = 0.0;
            double n2 = 0.0;
            const double d = dir.at(r, c);
            if (d == 0.0) {
                n1 = mag.at(r, c - 1);
                n2 = mag.at(r, c + 1);
            } else if (d == 45.0) {
                n1 = mag.at(r - 1, c + 1);
                n2 = mag.at(r + 1, c - 1);
            } else if (d == 90.0) {
                n1 = mag.at(r - 1, c);
                n2 = mag.at(r + 1, c);
            } else {
                n1 = mag.at(r - 1, c - 1);
                n2 = mag.at(r + 1, c + 1);
            }
            nms.at(r, c) = (m >= n1 && m >= n2) ? m : 0.0;
        }
    }

    // Automatic hysteresis thresholds on the gradient magnitude (the file
    // header explains why, and why not skimage's own default): the strong
    // threshold is the 80th percentile of the nonzero, non-maximum-suppressed
    // magnitude and the weak threshold half of that, the classic
    // percentile-based automatic-Canny heuristic. A median-based heuristic
    // (Stripenn's own auto_canny, applied to the image rather than the
    // gradient) was tried first and rejected: on a single sharp step edge
    // the surviving gradient magnitudes are all close to the same value, so
    // 1.5 times their median exceeds every one of them and no pixel ever
    // reaches the strong threshold.
    std::vector<double> nonzero;
    nonzero.reserve(static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols));
    for (const double v : nms.v) {
        if (v > 0.0) {
            nonzero.push_back(v);
        }
    }
    if (nonzero.empty()) {
        return Image(rows, cols, 0.0);
    }
    const double high = quantile(nonzero, 0.8);
    const double low = std::max(0.0, 0.4 * high);

    // Hysteresis: strong pixels seed a BFS through weak-but-connected ones,
    // the real multi-hop hysteresis a Canny detector needs (the ImageProcessing.Canny
    // helper the Python codebase carries but does not use only checks one
    // hop, which is not what skimage.feature.canny -- the function actually
    // called -- does).
    Image edges(rows, cols, 0.0);
    std::vector<std::pair<int, int>> stack;
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            if (nms.at(r, c) >= high) {
                edges.at(r, c) = 1.0;
                stack.emplace_back(r, c);
            }
        }
    }
    while (!stack.empty()) {
        const auto [r, c] = stack.back();
        stack.pop_back();
        for (int dr = -1; dr <= 1; ++dr) {
            for (int dc = -1; dc <= 1; ++dc) {
                if (dr == 0 && dc == 0) {
                    continue;
                }
                const int rr = r + dr;
                const int cc = c + dc;
                if (rr < 0 || rr >= rows || cc < 0 || cc >= cols) {
                    continue;
                }
                if (edges.at(rr, cc) != 0.0) {
                    continue;
                }
                if (nms.at(rr, cc) >= low) {
                    edges.at(rr, cc) = 1.0;
                    stack.emplace_back(rr, cc);
                }
            }
        }
    }
    return edges;
}

Image vertical_line(const Image& edges, double low_degrees, double high_degrees) {
    const int rows = edges.rows;
    const int cols = edges.cols;
    static constexpr int kGx[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    static constexpr int kGy[3][3] = {{1, 2, 1}, {0, 0, 0}, {-1, -2, -1}};
    Image out(rows, cols, 0.0);
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            double fx = 0.0;
            double fy = 0.0;
            for (int dr = -1; dr <= 1; ++dr) {
                for (int dc = -1; dc <= 1; ++dc) {
                    const int rr = r + dr;
                    const int cc = c + dc;
                    const double v = (rr >= 0 && rr < rows && cc >= 0 && cc < cols) ? edges.at(rr, cc) : 0.0;
                    fx += v * kGx[dr + 1][dc + 1];
                    fy += v * kGy[dr + 1][dc + 1];
                }
            }
            double orientation = std::atan2(fx, fy) * 180.0 / kPi;
            if (orientation < 0) {
                orientation += 360.0;
            }
            if (orientation > low_degrees && orientation < high_degrees) {
                const int shifted = c - 1;
                if (shifted >= 0 && shifted < cols) {
                    out.at(r, shifted) = 1.0;
                }
            }
        }
    }
    return out;
}

BlockResult block_scan(const Image& vert, int column) {
    const int rows = vert.rows;
    const int cols = vert.cols;
    std::vector<double> value(static_cast<std::size_t>(rows), 0.0);
    for (int r = 0; r < rows; ++r) {
        double m = 0.0;
        for (int c = std::max(0, column - 1); c < std::min(cols, column + 2); ++c) {
            m = std::max(m, vert.at(r, c));
        }
        value[static_cast<std::size_t>(r)] = m;
    }
    int count = 0;
    int max_count = 0;
    int end = column;
    int buffer = 0;
    int last_one = column;
    bool seen_one = false;
    for (int r = 0; r < rows; ++r) {
        if (value[static_cast<std::size_t>(r)] == 1.0) {
            ++count;
            last_one = r;
            seen_one = true;
        } else if (buffer < 5) {
            ++buffer;
        } else {
            if (count > max_count) {
                max_count = count;
                end = last_one;
            }
            count = 0;
            buffer = 0;
        }
    }
    if (count > max_count) {
        max_count = count;
        end = last_one;
    }
    if (!seen_one) {
        return {0, column};
    }
    if (end < column) {
        end = end - max_count + 1;
    }
    return {max_count, end};
}

namespace {

double adjust_brightness(double v, double b) {
    if (v <= 0.0) {
        return 0.0;
    }
    if (v > b) {
        return 1.0;
    }
    return v / b;
}

}  // namespace

std::vector<FrameCandidate> stripe_search_frame(const Image& submat, double maxpixel_value,
                                                 double canny_sigma, int min_length, int max_width,
                                                 int blur_filter) {
    std::vector<FrameCandidate> out;
    const int S = submat.rows;
    if (S < 5 || maxpixel_value <= 0.0) {
        return out;
    }

    // The blue = green channel intensity: bright where the count is low, dark
    // near or above M. Red is constant, so the RGB2GRAY combination reduces
    // to a fixed affine map of this single channel (see stripes_impl.hpp).
    Image intensity(S, S);
    for (int r = 0; r < S; ++r) {
        for (int c = 0; c < S; ++c) {
            const double value = 255.0 * (maxpixel_value - submat.at(r, c)) / maxpixel_value;
            intensity.at(r, c) = std::clamp(value / 255.0, 0.0, 1.0);
        }
    }

    for (double b = 0.5; b <= 1.0001; b += 0.1) {
        Image adjusted(S, S);
        for (std::size_t i = 0; i < adjusted.v.size(); ++i) {
            adjusted.v[i] = adjust_brightness(intensity.v[i], b);
        }
        // gray = 0.299 * red(=1) + 0.587 * green + 0.114 * blue, green == blue
        // here, so gray = 0.299 + 0.701 * adjusted. Blurring commutes with this
        // affine map, so blurring `adjusted` directly (rather than the three
        // channels separately, then combining) is exact, not an approximation.
        const Image blurred = box_blur(adjusted, std::max(1, blur_filter));
        Image gray(S, S);
        for (std::size_t i = 0; i < gray.v.size(); ++i) {
            gray.v[i] = std::clamp(0.299 + 0.701 * blurred.v[i], 0.0, 1.0);
        }

        const Image edges = canny(gray, canny_sigma);
        const Image vert = vertical_line(edges, 60.0, 120.0);

        std::vector<int> test_column;
        std::vector<int> end_points;
        std::vector<int> updown;
        for (int c = 0; c < S; ++c) {
            const BlockResult block = block_scan(vert, c);
            const int above = std::min(c, block.end);
            const int bottom = std::max(c, block.end);
            double column_sum = 0.0;
            for (int r = above; r <= bottom; ++r) {
                column_sum += vert.at(r, c);
            }
            if (block.length > min_length && column_sum != 0.0) {
                test_column.push_back(c);
                end_points.push_back(block.end);
                updown.push_back(block.end > c ? 2 : 1);
            }
        }

        for (const int ud : {1, 2}) {
            std::vector<std::uint8_t> testmat(static_cast<std::size_t>(S) * static_cast<std::size_t>(S), 0);
            const auto set = [&](int r, int c) { testmat[static_cast<std::size_t>(r) * S + c] = 1; };
            const auto get = [&](int r, int c) { return testmat[static_cast<std::size_t>(r) * S + c]; };
            for (std::size_t i = 0; i < test_column.size(); ++i) {
                if (updown[i] != ud) {
                    continue;
                }
                int st = test_column[i];
                int en = end_points[i];
                if (ud == 1) {
                    std::swap(st, en);
                }
                const int lo = std::min(st, en);
                const int hi = std::max(st, en);
                for (int r = lo; r <= hi; ++r) {
                    set(r, test_column[i]);
                }
            }

            // Per-column run lengths, columns with length >= 3, grouped into
            // contiguous runs of adjacent columns (stripenn.py:980-1028).
            std::map<int, int> column_length;
            for (int c = 0; c < S; ++c) {
                int len = 0;
                for (int r = 0; r < S; ++r) {
                    len += get(r, c);
                }
                if (len >= 3) {
                    column_length[c] = len;
                }
            }
            std::vector<double> boundary_positions;
            {
                std::vector<int> run;
                std::vector<int> run_len;
                const auto flush = [&]() {
                    if (run.empty()) {
                        return;
                    }
                    double weight_sum = 0.0;
                    for (const int w : run_len) {
                        weight_sum += w;
                    }
                    double acc = 0.0;
                    if (weight_sum > 0.0) {
                        for (std::size_t i = 0; i < run.size(); ++i) {
                            acc += static_cast<double>(run[i]) * (static_cast<double>(run_len[i]) / weight_sum);
                        }
                    } else {
                        acc = run.back();
                    }
                    boundary_positions.push_back(std::round(acc));
                    run.clear();
                    run_len.clear();
                };
                int previous = std::numeric_limits<int>::min();
                for (const auto& [col, len] : column_length) {
                    if (!run.empty() && col - previous != 1) {
                        flush();
                    }
                    run.push_back(col);
                    run_len.push_back(len);
                    previous = col;
                }
                flush();
            }
            std::sort(boundary_positions.begin(), boundary_positions.end());
            boundary_positions.erase(std::unique(boundary_positions.begin(), boundary_positions.end()),
                                     boundary_positions.end());

            for (std::size_t i = 0; i + 1 < boundary_positions.size(); ++i) {
                const int n = static_cast<int>(boundary_positions[i]);
                const int m = static_cast<int>(boundary_positions[i + 1]);
                const int width = std::abs(m - n);
                if (width <= 1 || width > max_width) {
                    continue;
                }
                int pair_start = n;
                int pair_end = m;
                if (width > 4) {
                    pair_end = m - 2;
                }
                const auto column_extent = [&](int center) {
                    int lo = std::numeric_limits<int>::max();
                    int hi = std::numeric_limits<int>::min();
                    for (int c = std::max(0, center - 1); c < std::min(S, center + 2); ++c) {
                        for (int r = 0; r < S; ++r) {
                            if (get(r, c) == 1) {
                                lo = std::min(lo, r);
                                hi = std::max(hi, r);
                            }
                        }
                    }
                    return std::pair<int, int>(lo, hi);
                };
                const auto [min1, max1] = column_extent(n);
                const auto [min2, max2] = column_extent(m);
                if (min1 > max1 || min2 > max2) {
                    continue;
                }
                int lo = std::min(min1, min2);
                int hi = std::max(max1, max2);
                if (ud == 1) {
                    hi = (width > 4) ? (m - 2) : m;
                } else {
                    lo = n;
                }
                if (hi < lo) {
                    continue;
                }
                FrameCandidate candidate;
                candidate.x = std::min(pair_start, pair_end);
                candidate.w = std::abs(pair_end - pair_start) + 1;
                candidate.y = lo;
                candidate.h = hi - lo + 1;
                if (candidate.x >= 0 && candidate.y >= 0 && candidate.x + candidate.w <= S &&
                    candidate.y + candidate.h <= S) {
                    out.push_back(candidate);
                }
            }
        }
    }
    return out;
}

std::vector<Candidate> remove_redundant(std::vector<Candidate> candidates, bool by_pvalue,
                                        unsigned int threads) {
    (void)threads;  // the per-chromosome group is small enough to run in one thread
    std::map<std::string, std::vector<std::size_t>> by_chrom;
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        by_chrom[candidates[i].chrom].push_back(i);
    }
    std::vector<bool> keep(candidates.size(), true);
    for (auto& [chrom, indices] : by_chrom) {
        (void)chrom;
        std::sort(indices.begin(), indices.end(), [&](std::size_t a, std::size_t b) {
            return candidates[a].frame_index < candidates[b].frame_index;
        });
        for (std::size_t ii = 0; ii < indices.size(); ++ii) {
            for (std::size_t jj = ii + 1; jj < indices.size(); ++jj) {
                const std::size_t i = indices[ii];
                const std::size_t j = indices[jj];
                if (std::abs(candidates[i].frame_index - candidates[j].frame_index) > 1) {
                    continue;
                }
                const Candidate& a = candidates[i];
                const Candidate& c = candidates[j];
                const std::int64_t ix_lo = std::max(a.pos1, c.pos1);
                const std::int64_t ix_hi = std::min(a.pos2, c.pos2);
                const std::int64_t iy_lo = std::max(a.pos3, c.pos3);
                const std::int64_t iy_hi = std::min(a.pos4, c.pos4);
                if (ix_hi < ix_lo || iy_hi < iy_lo) {
                    continue;
                }
                const double sx = static_cast<double>(ix_hi - ix_lo + 1) /
                                  static_cast<double>(std::min(a.pos2 - a.pos1, c.pos2 - c.pos1) + 1);
                const double sy = static_cast<double>(iy_hi - iy_lo + 1) /
                                  static_cast<double>(std::min(a.pos4 - a.pos3, c.pos4 - c.pos3) + 1);
                if (sx <= 0.2 || sy <= 0.2) {
                    continue;
                }
                if (!keep[i] || !keep[j]) {
                    continue;
                }
                if (by_pvalue) {
                    if (a.pvalue > c.pvalue) {
                        keep[i] = false;
                    } else {
                        keep[j] = false;
                    }
                } else {
                    const double ratio_a = static_cast<double>(a.pos4 - a.pos3) /
                                          static_cast<double>(std::max<std::int64_t>(1, a.pos2 - a.pos1));
                    const double ratio_c = static_cast<double>(c.pos4 - c.pos3) /
                                          static_cast<double>(std::max<std::int64_t>(1, c.pos2 - c.pos1));
                    if (ratio_a <= ratio_c) {
                        keep[i] = false;
                    } else {
                        keep[j] = false;
                    }
                }
            }
        }
    }
    std::vector<Candidate> result;
    result.reserve(candidates.size());
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        if (keep[i]) {
            result.push_back(std::move(candidates[i]));
        }
    }
    return result;
}

double candidate_pvalue(const BackgroundModel& model, int distance, double left_diff,
                        double right_diff) {
    const std::size_t d = static_cast<std::size_t>(std::clamp(distance, 0, 399));
    if (d >= model.left.size() || model.left[d].empty()) {
        return 1.0;
    }
    const auto rank = [](const std::vector<double>& sample, double observed) {
        std::size_t at_least = 0;
        for (const double v : sample) {
            if (v >= observed) {
                ++at_least;
            }
        }
        double p = static_cast<double>(at_least) / static_cast<double>(sample.size());
        if (p == 0.0) {
            p = 1.0 / static_cast<double>(sample.size());
        }
        return p;
    };
    const double p1 = rank(model.left[d], left_diff);
    const double p2 = rank(model.right[d], right_diff);
    return std::max(p1, p2);
}

}  // namespace hicx::stripes
