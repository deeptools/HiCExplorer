// Port of hicexplorer/lib/viewpoint.py. See hicx/chic_viewpoint.hpp.

#include "hicx/chic_viewpoint.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>

#include "hicx/numpy_compat.hpp"
#include "hicx/scipy_special.hpp"

namespace hicx::chic {

namespace {

bool is_python_space(char c) {
    // The ASCII characters str.strip() and str.split() treat as whitespace.
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\x0b' || c == '\x0c' ||
           c == '\x1c' || c == '\x1d' || c == '\x1e' || c == '\x1f';
}

std::vector<std::string> split_tab(std::string_view text) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true) {
        const std::size_t tab = text.find('\t', begin);
        if (tab == std::string_view::npos) {
            fields.emplace_back(text.substr(begin));
            return fields;
        }
        fields.emplace_back(text.substr(begin, tab - begin));
        begin = tab + 1;
    }
}

// numpy's normalisation of a slice bound for an array of length n.
std::int64_t slice_index(std::int64_t index, std::int64_t n) {
    if (index < 0) {
        index += n;
        if (index < 0) {
            index = 0;
        }
    } else if (index > n) {
        index = n;
    }
    return index;
}

std::int64_t slice_length(std::int64_t start, std::int64_t stop, std::int64_t n) {
    const std::int64_t begin = slice_index(start, n);
    const std::int64_t end = slice_index(stop, n);
    return std::max<std::int64_t>(0, end - begin);
}

// target[a:b] = source[c:d] with numpy's broadcasting: the lengths must be
// equal, or the source must have exactly one element.
void slice_assign(std::vector<double>& target, std::int64_t a, std::int64_t b,
                  const std::vector<double>& source, std::int64_t c, std::int64_t d) {
    const std::int64_t target_n = static_cast<std::int64_t>(target.size());
    const std::int64_t source_n = static_cast<std::int64_t>(source.size());
    const std::int64_t target_len = slice_length(a, b, target_n);
    const std::int64_t source_len = slice_length(c, d, source_n);
    const std::int64_t target_begin = slice_index(a, target_n);
    const std::int64_t source_begin = slice_index(c, source_n);
    if (source_len == target_len) {
        for (std::int64_t i = 0; i < target_len; ++i) {
            target[static_cast<std::size_t>(target_begin + i)] =
                source[static_cast<std::size_t>(source_begin + i)];
        }
        return;
    }
    if (source_len == 1) {
        for (std::int64_t i = 0; i < target_len; ++i) {
            target[static_cast<std::size_t>(target_begin + i)] =
                source[static_cast<std::size_t>(source_begin)];
        }
        return;
    }
    throw ValueError("could not broadcast input array from shape (" +
                     std::to_string(source_len) + ",) into shape (" +
                     std::to_string(target_len) + ",)");
}

template <class T>
std::vector<double> smooth_impl(std::span<const T> data, std::int64_t window) {
    // window_size = np.int32(np.floor(pWindowSize / 2))
    const std::int64_t half = static_cast<std::int64_t>(
        std::floor(static_cast<double>(window) / 2.0));
    std::int64_t upstream = half;
    if (window % 2 == 0) {
        upstream -= 1;
    }
    const std::int64_t n = static_cast<std::int64_t>(data.size());
    std::vector<double> average(static_cast<std::size_t>(n), 0.0);

    // np.mean(data[begin:end]) with Python slice clipping.
    const auto mean = [&](std::int64_t begin, std::int64_t end) -> double {
        const std::int64_t b = slice_index(begin, n);
        const std::int64_t e = slice_index(end, n);
        const std::size_t count = e > b ? static_cast<std::size_t>(e - b) : 0;
        if (count == 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        const T sum = npy::pairwise_sum(data.data() + b, count);
        return static_cast<double>(sum / static_cast<T>(count));
    };
    // average[i] for a Python index that may be negative.
    const auto slot = [&](std::int64_t index) -> double& {
        std::int64_t position = index < 0 ? index + n : index;
        if (position < 0 || position >= n) {
            throw IndexError("index " + std::to_string(index) +
                             " is out of bounds for axis 0 with size " + std::to_string(n));
        }
        return average[static_cast<std::size_t>(position)];
    };

    for (std::int64_t i = upstream; i < n - half; ++i) {
        slot(i) = mean(i - upstream, i + half + 1);
    }
    for (std::int64_t i = 0; i < half; ++i) {
        std::int64_t start = i - upstream;
        if (start < 0) {
            start = 0;
        }
        const std::int64_t end = i + half + 1;
        slot(i) = mean(start, end);
        // np.mean(pData[-end:])
        slot(-(i + 1)) = mean(-end, n);
    }
    return average;
}

}  // namespace

// ---------------------------------------------------------------------------
// Text helpers

std::vector<std::string> read_lines(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("[Errno 2] No such file or directory: '" + path + "'");
    }
    std::ostringstream buffer;
    buffer << in.rdbuf();
    const std::string text = buffer.str();
    std::vector<std::string> lines;
    std::string current;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];
        if (c == '\r') {
            if (i + 1 < text.size() && text[i + 1] == '\n') {
                ++i;
            }
            current.push_back('\n');
            lines.push_back(std::move(current));
            current.clear();
        } else if (c == '\n') {
            current.push_back('\n');
            lines.push_back(std::move(current));
            current.clear();
        } else {
            current.push_back(c);
        }
    }
    if (!current.empty()) {
        lines.push_back(std::move(current));
    }
    return lines;
}

std::string_view strip(std::string_view text) {
    std::size_t begin = 0;
    std::size_t end = text.size();
    while (begin < end && is_python_space(text[begin])) {
        ++begin;
    }
    while (end > begin && is_python_space(text[end - 1])) {
        --end;
    }
    return text.substr(begin, end - begin);
}

std::int64_t python_int(std::string_view text) {
    const std::string_view body = strip(text);
    std::string digits;
    std::size_t i = 0;
    bool negative = false;
    if (i < body.size() && (body[i] == '+' || body[i] == '-')) {
        negative = body[i] == '-';
        ++i;
    }
    bool previous_digit = false;
    for (; i < body.size(); ++i) {
        const char c = body[i];
        if (c >= '0' && c <= '9') {
            digits.push_back(c);
            previous_digit = true;
        } else if (c == '_' && previous_digit && i + 1 < body.size() && body[i + 1] >= '0' &&
                   body[i + 1] <= '9') {
            previous_digit = false;
        } else {
            digits.clear();
            break;
        }
    }
    std::int64_t value = 0;
    if (digits.empty() ||
        std::from_chars(digits.data(), digits.data() + digits.size(), value).ec != std::errc()) {
        throw ValueError("invalid literal for int() with base 10: '" + std::string(text) + "'");
    }
    return negative ? -value : value;
}

double python_float(std::string_view text) {
    std::string body(strip(text));
    std::string lowered;
    for (const char c : body) {
        lowered.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    std::string unsigned_part = lowered;
    double sign = 1.0;
    if (!unsigned_part.empty() && (unsigned_part[0] == '+' || unsigned_part[0] == '-')) {
        sign = unsigned_part[0] == '-' ? -1.0 : 1.0;
        unsigned_part = unsigned_part.substr(1);
    }
    if (unsigned_part == "nan") {
        return std::copysign(std::numeric_limits<double>::quiet_NaN(), sign);
    }
    if (unsigned_part == "inf" || unsigned_part == "infinity") {
        return sign * std::numeric_limits<double>::infinity();
    }
    double value = 0.0;
    const auto result = std::from_chars(body.data(), body.data() + body.size(), value);
    bool ok = result.ec == std::errc() && result.ptr == body.data() + body.size() && !body.empty();
    if (!ok && !body.empty() && body[0] == '+') {
        const auto retry = std::from_chars(body.data() + 1, body.data() + body.size(), value);
        ok = retry.ec == std::errc() && retry.ptr == body.data() + body.size() && body.size() > 1;
    }
    if (!ok) {
        throw ValueError("could not convert string to float: '" + std::string(text) + "'");
    }
    return value;
}

std::string format_fixed(double value, int decimals) {
    if (std::isnan(value)) {
        return "nan";
    }
    if (std::isinf(value)) {
        return value < 0 ? "-inf" : "inf";
    }
    // glibc's printf is exact, and so is Python's float formatting, so the
    // two agree on every correctly rounded fixed notation.
    // The largest finite double has 309 integer digits; 64 decimals are more
    // than any caller asks for.
    decimals = std::clamp(decimals, 0, 64);
    char buffer[400];
    const int size = std::snprintf(buffer, sizeof(buffer), "%.*f", decimals, value);
    return std::string(buffer, static_cast<std::size_t>(size));
}

// ---------------------------------------------------------------------------
// Reference points

ReferencePoints read_reference_points(const std::string& path) {
    ReferencePoints result;
    for (const std::string& line : read_lines(path)) {
        const std::vector<std::string> fields = split_tab(strip(line));
        if (fields.size() == 3) {
            result.points.push_back(ReferencePoint{fields[0], fields[1], fields[1]});
            result.genes.push_back(fields[2]);
        } else if (fields.size() > 3) {
            result.points.push_back(ReferencePoint{fields[0], fields[1], fields[2]});
            result.genes.push_back(fields[3]);
        }
    }
    return result;
}

std::string reference_point_string(const ReferencePoint& point) {
    return point.chromosome + "_" + point.start + "_" + point.end;
}

// ---------------------------------------------------------------------------
// ViewpointMatrix

ViewpointMatrix ViewpointMatrix::load(const std::string& path) {
    ViewpointMatrix result;
    result.matrix_ = HiCMatrix::load(path);
    return result;
}

std::vector<std::string> ViewpointMatrix::chromosome_names() const {
    return bins().chrom_names();
}

bool ViewpointMatrix::has_chromosome(const std::string& chromosome) const {
    return bins().chrom_bin_range(chromosome).has_value();
}

std::optional<std::pair<std::int64_t, std::int64_t>> ViewpointMatrix::region_bin_range(
    const std::string& chromosome, std::int64_t start, std::int64_t end) const {
    if (!has_chromosome(chromosome)) {
        throw ValueError("chromosome: " + chromosome + " name not found in matrix");
    }
    return bins().region_bin_range(chromosome, start, end);
}

std::pair<std::int64_t, std::int64_t> ViewpointMatrix::reference_point_indices(
    const ReferencePoint& point) const {
    const std::int64_t start = python_int(point.start);
    const std::int64_t end = python_int(point.end);
    const auto range = region_bin_range(point.chromosome, start, end);
    if (!range.has_value()) {
        throw TypeError("cannot unpack non-iterable NoneType object");
    }
    return *range;
}

const CutInterval& ViewpointMatrix::bin_position(std::int64_t index) const {
    if (index < 0) {
        // Python indexes the cut interval list from the end; no caller in the
        // cHi-C tools can produce a negative bin, so this is not reproduced.
        throw ValueError("binIndex: " + std::to_string(index) + " not found");
    }
    if (static_cast<std::size_t>(index) >= bins().size()) {
        throw ValueError("binIndex: " + std::to_string(index) + " not found");
    }
    return bins().bin_pos(static_cast<std::size_t>(index));
}

double ViewpointMatrix::value(std::int64_t row, std::int64_t column) const {
    return matrix().at(row, column);
}

// ---------------------------------------------------------------------------
// Ranges and viewpoints

ViewpointRange calculate_viewpoint_range(const ViewpointMatrix& matrix,
                                         const ReferencePoint& point, std::int64_t upstream,
                                         std::int64_t downstream) {
    const std::optional<BinRange> chromosome = matrix.bins().chrom_bin_range(point.chromosome);
    if (!chromosome.has_value()) {
        throw ValueError("chrName: " + point.chromosome + " not found in chrBinBoundaries");
    }
    const std::int64_t max_length = matrix.bin_position(chromosome->last - 1).end;
    const std::int64_t bin_size = matrix.bin_size();
    ViewpointRange result;
    result.upstream = upstream;
    result.downstream = downstream;
    result.region_start = python_int(point.start) - upstream;
    if (result.region_start < 0) {
        result.region_start = 0;
        result.upstream = python_int(point.start);
    }
    result.region_end = python_int(point.end) + downstream;
    if (result.region_end > max_length) {
        // "-1 is important, otherwise self.hicMatrix.getRegionBinRange will crash"
        result.region_end = max_length - 1;
        result.downstream = (max_length - python_int(point.end)) + bin_size;
    }
    return result;
}

ComputedViewpoint compute_viewpoint(const ViewpointMatrix& matrix, const ReferencePoint& point,
                                    const std::string& chromosome, std::int64_t region_start,
                                    std::int64_t region_end) {
    const auto [view_point_start, view_point_end] = matrix.reference_point_indices(point);
    const auto range = matrix.region_bin_range(chromosome, region_start, region_end);
    if (!range.has_value()) {
        throw TypeError("'NoneType' object is not iterable");
    }
    const std::int64_t range_start = range->first;
    const std::int64_t range_end = range->second + 1;

    const std::int64_t elements = range_end - range_start;
    if (elements < 0) {
        throw ValueError("negative dimensions are not allowed");
    }
    std::vector<double> data(static_cast<std::size_t>(elements), 0.0);
    for (std::int64_t row = view_point_start; row <= view_point_end; ++row) {
        (void)matrix.bin_position(row);
        for (std::int64_t j = 0; j < elements; ++j) {
            data[static_cast<std::size_t>(j)] += matrix.value(row, range_start + j);
        }
    }

    const std::int64_t elements_new = elements - (view_point_end - view_point_start);
    if (elements_new < 0) {
        throw ValueError("negative dimensions are not allowed");
    }
    std::vector<double> data_new(static_cast<std::size_t>(elements_new), 0.0);
    const std::int64_t index_before_viewpoint = view_point_start - range_start;
    const std::int64_t width = view_point_end - view_point_start;

    // data_list_new[0:index_before_viewpoint] = data_list[0:index_before_viewpoint]
    slice_assign(data_new, 0, index_before_viewpoint, data, 0, index_before_viewpoint);

    // data_list_new[index_before_viewpoint] = np.sum(data_list[i : i + width + 1])
    {
        const std::int64_t n = static_cast<std::int64_t>(data.size());
        const std::int64_t begin = slice_index(index_before_viewpoint, n);
        const std::int64_t end = slice_index(index_before_viewpoint + width + 1, n);
        const double sum = end > begin
                               ? npy::pairwise_sum(data.data() + begin,
                                                   static_cast<std::size_t>(end - begin))
                               : 0.0;
        std::int64_t position = index_before_viewpoint;
        if (position < 0) {
            position += elements_new;
        }
        if (position < 0 || position >= elements_new) {
            throw IndexError("index " + std::to_string(index_before_viewpoint) +
                             " is out of bounds for axis 0 with size " +
                             std::to_string(elements_new));
        }
        data_new[static_cast<std::size_t>(position)] = sum;
    }

    // data_list_new[i + 1:] = data_list[i + width + 1:]
    slice_assign(data_new, index_before_viewpoint + 1, elements_new, data,
                 index_before_viewpoint + width + 1, elements);

    return ComputedViewpoint{std::move(data_new), index_before_viewpoint};
}

std::vector<double> smooth_interaction_values(std::span<const double> data,
                                              std::int64_t window) {
    return smooth_impl(data, window);
}

std::vector<double> smooth_interaction_values(std::span<const float> data,
                                              std::int64_t window) {
    return smooth_impl(data, window);
}

std::vector<double> compute_relative_values(std::span<const double> data, double denominator) {
    std::vector<double> output(data.begin(), data.end());
    // `if pDenominator:` is False for 0.0 only; NaN is truthy.
    const double divisor = denominator != 0.0 || std::isnan(denominator)
                               ? denominator
                               : npy::pairwise_sum(output.data(), output.size());
    for (double& value : output) {
        value /= divisor;
    }
    return output;
}

// ---------------------------------------------------------------------------
// BackgroundModel

void BackgroundModel::set(std::int64_t key, std::vector<double> values) {
    const auto [it, inserted] = values_.insert_or_assign(key, std::move(values));
    (void)it;
    if (inserted) {
        order_.push_back(key);
    }
}

bool BackgroundModel::contains(std::int64_t key) const {
    return values_.find(key) != values_.end();
}

const std::vector<double>& BackgroundModel::at(std::int64_t key) const {
    const auto it = values_.find(key);
    if (it == values_.end()) {
        throw ValueError("KeyError: " + std::to_string(key));
    }
    return it->second;
}

std::int64_t BackgroundModel::min_key() const {
    if (values_.empty()) {
        throw ValueError("min() arg is an empty sequence");
    }
    return values_.begin()->first;
}

std::int64_t BackgroundModel::max_key() const {
    if (values_.empty()) {
        throw ValueError("max() arg is an empty sequence");
    }
    return values_.rbegin()->first;
}

std::vector<std::int64_t> BackgroundModel::sorted_keys() const {
    std::vector<std::int64_t> keys;
    keys.reserve(values_.size());
    for (const auto& entry : values_) {
        keys.push_back(entry.first);
    }
    return keys;
}

BackgroundModel read_background_model(const std::string& path, std::int64_t range_upstream,
                                      std::int64_t range_downstream, std::int64_t fixate_range,
                                      bool mean) {
    const std::vector<std::string> lines = read_lines(path);
    BackgroundModel model;
    for (std::size_t index = 1; index < lines.size(); ++index) {
        const std::vector<std::string> fields = split_tab(lines[index]);
        const std::int64_t key = python_int(fields[0]);
        if (!mean) {
            if (fields.size() < 4) {
                throw IndexError("list index out of range");
            }
            model.set(key, {python_float(fields[1]), python_float(fields[2]),
                            python_float(fields[3])});
        } else {
            model.set(key, {python_float(fields.back())});
        }
    }
    std::int64_t max_key = model.max_key();
    std::int64_t min_key = model.min_key();
    if (max_key > fixate_range) {
        max_key = fixate_range;
    }
    if (min_key < -fixate_range) {
        min_key = -fixate_range;
    }
    const std::vector<std::int64_t>& keys = model.insertion_order();
    if (keys.size() < 2) {
        throw IndexError("list index out of range");
    }
    const std::int64_t increment = std::llabs(std::llabs(keys[0]) - std::llabs(keys[1]));
    if (max_key < range_downstream) {
        if (increment == 0) {
            // The Python loops forever here: i never grows.
            throw ValueError("the background model's first two positions have the same "
                             "distance, so it cannot be extended to --range; the Python "
                             "reference does not terminate on this input");
        }
        const std::vector<double> value = model.at(max_key);
        std::int64_t i = max_key;
        while (i < range_downstream) {
            i += increment;
            model.set(i, value);
        }
    }
    if (min_key > -range_upstream) {
        if (increment == 0) {
            throw ValueError("the background model's first two positions have the same "
                             "distance, so it cannot be extended to --range; the Python "
                             "reference does not terminate on this input");
        }
        const std::vector<double> value = model.at(min_key);
        std::int64_t i = min_key;
        while (i > -range_upstream) {
            i -= increment;
            model.set(i, value);
        }
    }
    return model;
}

std::vector<double> interaction_background_data(const BackgroundModel& model,
                                                std::int64_t upstream,
                                                std::int64_t downstream) {
    std::vector<double> result;
    for (const std::int64_t key : model.sorted_keys()) {
        if (key >= -upstream && key <= downstream) {
            const std::vector<double>& values = model.at(key);
            result.insert(result.end(), values.begin(), values.end());
        }
    }
    return result;
}

std::vector<double> p_values(const BackgroundModel& model, std::span<const double> data,
                             std::int64_t index_reference_point, bool data_is_float32) {
    std::vector<double> result(data.size(), 0.0);
    for (std::size_t i = 0; i < data.size(); ++i) {
        std::int64_t relative_distance = static_cast<std::int64_t>(i) - index_reference_point;
        if (!model.contains(relative_distance)) {
            relative_distance = relative_distance < 0 ? model.min_key() : model.max_key();
        }
        const double element = data[i];
        double cdf = 0.0;
        if (element != 0.0) {
            const std::vector<double>& parameters = model.at(relative_distance);
            const double shifted = data_is_float32
                                       ? static_cast<double>(static_cast<float>(element) + 1.0F)
                                       : element + 1.0;
            cdf = scipy::betainc(parameters[0], shifted, parameters[1]);
        }
        double p = 1.0 - cdf;
        if (std::isnan(p) || std::isinf(p)) {
            p = 1.0;
        }
        result[i] = p;
    }
    return result;
}

// ---------------------------------------------------------------------------
// Interaction file content

InteractionFileData create_interaction_file_data(
    const ViewpointMatrix& matrix, const ReferencePoint& point, const std::string& chromosome,
    std::int64_t region_start, std::int64_t region_end,
    std::span<const double> interaction_data, std::span<const double> raw,
    const std::string& gene, double sum_of_interactions, std::vector<double> pvalues,
    std::vector<double> xfold) {
    const auto [view_point_start, view_point_end] = matrix.reference_point_indices(point);
    const auto range = matrix.region_bin_range(chromosome, region_start, region_end);
    if (!range.has_value()) {
        throw TypeError("'NoneType' object is not iterable");
    }
    const std::int64_t range_start = range->first;
    const std::int64_t range_end = range->second + 1;

    InteractionFileData result;
    result.gene = gene;
    result.sum_of_interactions = sum_of_interactions;
    result.pvalues = std::move(pvalues);
    result.xfold = std::move(xfold);

    const std::int64_t start = matrix.bin_position(view_point_start).start;
    const std::int64_t end = matrix.bin_position(view_point_end).end;

    std::vector<std::int64_t> positions;
    for (std::int64_t i = range_start; i < view_point_start; ++i) {
        positions.push_back(i);
    }
    positions.push_back(view_point_start);
    for (std::int64_t i = view_point_end + 1; i < range_end; ++i) {
        positions.push_back(i);
    }

    const std::size_t count = std::min(interaction_data.size(), positions.size());
    std::int64_t relative_position = -1;
    for (std::size_t j = 0; j < count; ++j) {
        const CutInterval& bin = matrix.bin_position(positions[j]);
        if (relative_position < 0) {
            relative_position = bin.start - start;
        } else {
            relative_position = bin.end - end;
        }
        result.chromosome = bin.chrom;
        result.starts.push_back(bin.start);
        result.ends.push_back(bin.end);
        result.relative_positions.push_back(relative_position);
        result.interaction_data.push_back(interaction_data[j]);
        result.raw.push_back(raw[j]);
    }
    return result;
}

}  // namespace hicx::chic
