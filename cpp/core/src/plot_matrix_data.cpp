#include "hicx/plot_matrix_data.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
#include <stdexcept>

namespace hicx::plot_matrix {

namespace {

// Python's '{:.Nf}'.format(x)
std::string fixed(double x, int digits) {
    if (std::isnan(x)) {
        return "nan";
    }
    if (std::isinf(x)) {
        return x > 0 ? "inf" : "-inf";
    }
    char buffer[512];
    std::snprintf(buffer, sizeof(buffer), "%.*f", digits, x);
    return buffer;
}

// np.nanmin over the values selected by `keep`; nullopt for an empty
// selection (numpy raises ValueError).
std::optional<double> nanmin(const std::vector<double>& values, const std::vector<char>& excluded) {
    bool any = false;
    double best = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t k = 0; k < values.size(); ++k) {
        if (excluded[k]) {
            continue;
        }
        any = true;
        const double v = values[k];
        if (std::isnan(v)) {
            continue;
        }
        if (std::isnan(best) || v < best) {
            best = v;
        }
    }
    if (!any) {
        return std::nullopt;
    }
    return best;
}

}  // namespace

std::vector<std::string> relabel_ticks(const std::vector<double>& ticks) {
    if (ticks.empty()) {
        throw std::out_of_range("index -1 is out of bounds");
    }
    std::vector<std::string> labels;
    const double last = ticks.back();
    const char* unit = " bp";
    if (last > 1.5e6) {
        for (const double x : ticks) {
            labels.push_back(fixed(x / 1e6, 2) + " ");
        }
        unit = " Mbp";
    } else if (last > 1500) {
        for (const double x : ticks) {
            labels.push_back(fixed(x / 1e3, 0));
        }
        unit = " Kbp";
    } else {
        for (const double x : ticks) {
            labels.push_back(fixed(x, 2) + " ");
        }
    }
    if (labels.size() < 2) {
        throw std::out_of_range("list index out of range");
    }
    labels[labels.size() - 2] += unit;
    return labels;
}

std::vector<std::string> region_fields(std::string_view region) {
    std::string cleaned;
    for (const char c : region) {
        if (c == ',' || c == ';' || c == '!') {
            continue;
        }
        cleaned.push_back(c == '-' ? ':' : c);
    }
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (true) {
        const std::size_t colon = cleaned.find(':', start);
        if (colon == std::string::npos) {
            fields.push_back(cleaned.substr(start));
            break;
        }
        fields.push_back(cleaned.substr(start, colon - start));
        start = colon + 1;
    }
    return fields;
}

LogClearing clear_for_log(std::vector<double>& matrix) {
    LogClearing out;
    const std::size_t n = matrix.size();
    std::vector<char> zero(n);
    std::vector<char> nan(n);
    std::vector<char> inf(n);
    for (std::size_t k = 0; k < n; ++k) {
        zero[k] = matrix[k] == 0.0 ? 1 : 0;
        nan[k] = std::isnan(matrix[k]) ? 1 : 0;
        inf[k] = std::isinf(matrix[k]) ? 1 : 0;
        out.any_nan_before = out.any_nan_before || nan[k];
        out.any_inf_before = out.any_inf_before || inf[k];
    }
    for (const std::vector<char>* mask : {&zero, &nan, &inf}) {
        const std::optional<double> value = nanmin(matrix, *mask);
        if (!value) {
            out.failed = true;
            break;
        }
        for (std::size_t k = 0; k < n; ++k) {
            if ((*mask)[k]) {
                matrix[k] = *value;
            }
        }
    }
    for (const double v : matrix) {
        out.any_nan_after = out.any_nan_after || std::isnan(v);
        out.any_inf_after = out.any_inf_after || std::isinf(v);
    }
    return out;
}

void clip(std::vector<double>& values, const double* low, const double* high) {
    for (double& v : values) {
        if (low != nullptr) {
            // np.maximum propagates NaN from either side
            v = (std::isnan(v) || std::isnan(*low)) ? std::numeric_limits<double>::quiet_NaN()
                                                    : (v < *low ? *low : v);
        }
        if (high != nullptr) {
            v = (std::isnan(v) || std::isnan(*high)) ? std::numeric_limits<double>::quiet_NaN()
                                                     : (v > *high ? *high : v);
        }
    }
}

}  // namespace hicx::plot_matrix
