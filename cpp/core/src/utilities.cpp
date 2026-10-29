#include "hicx/utilities.hpp"

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <limits>
#include <type_traits>
#include <unordered_map>

#include "hicx/numpy_compat.hpp"
#include "hicx/py_text.hpp"
#include "hicx/python_error.hpp"

namespace hicx::utilities {

namespace unidecode_table {
extern const char kBlob[];
extern const std::uint32_t kOffsets[];
extern const std::int32_t kBlockFirst[];
extern const std::uint16_t kBlockCount[];
extern const int kBlocks;
}  // namespace unidecode_table

namespace {

template <class T>
constexpr bool is_float_v = std::is_floating_point_v<T>;

// numpy's result dtype of np.sum over T: float stays float, integers sum in
// int64 (np.int_).
template <class T>
using SumType = std::conditional_t<is_float_v<T>, T, std::int64_t>;

// Wrapping integer addition and multiplication, as numpy's int64 does.
inline std::int64_t wrap_add(std::int64_t a, std::int64_t b) noexcept {
    return static_cast<std::int64_t>(static_cast<std::uint64_t>(a) + static_cast<std::uint64_t>(b));
}
inline std::int32_t wrap_add(std::int32_t a, std::int32_t b) noexcept {
    return static_cast<std::int32_t>(static_cast<std::uint32_t>(a) + static_cast<std::uint32_t>(b));
}
inline std::int64_t wrap_mul(std::int64_t a, std::int64_t b) noexcept {
    return static_cast<std::int64_t>(static_cast<std::uint64_t>(a) * static_cast<std::uint64_t>(b));
}

// np.sum of a contiguous array of T, converted to double.
template <class T>
double np_sum(const std::vector<T>& values) {
    if constexpr (is_float_v<T>) {
        return static_cast<double>(npy::pairwise_sum(values.data(), values.size()));
    } else {
        std::int64_t sum = 0;
        for (const T v : values) {
            sum = wrap_add(sum, static_cast<std::int64_t>(v));
        }
        return static_cast<double>(sum);
    }
}

inline double zero_non_finite(double value) noexcept {
    return std::isfinite(value) ? value : 0.0;
}

[[noreturn]] void raise_index(std::int64_t index, std::int64_t size) {
    throw PythonError("IndexError", index_out_of_bounds_message(index, size));
}

[[noreturn]] void raise_mismatch(std::size_t dimension, std::size_t boolean_dimension) {
    throw PythonError("IndexError", boolean_index_mismatch_message(dimension, boolean_dimension));
}

std::string format_fixed(double value, int decimals) {
    if (std::isnan(value)) {
        return "nan";  // Python drops the sign of a NaN
    }
    char buffer[512];
    std::snprintf(buffer, sizeof(buffer), "%.*f", decimals, value);
    return buffer;
}

// Decodes one code point of surrogatepass UTF-8 at text[pos].
char32_t decode_at(std::string_view text, std::size_t pos, std::size_t& length) noexcept {
    const auto b0 = static_cast<unsigned char>(text[pos]);
    if (b0 < 0x80U) {
        length = 1;
        return b0;
    }
    auto cont = [&](std::size_t k) {
        return pos + k < text.size() ? static_cast<unsigned char>(text[pos + k]) & 0x3FU : 0U;
    };
    if (b0 < 0xE0U) {
        length = 2;
        return (static_cast<char32_t>(b0 & 0x1FU) << 6) | cont(1);
    }
    if (b0 < 0xF0U) {
        length = 3;
        return (static_cast<char32_t>(b0 & 0x0FU) << 12) | (cont(1) << 6) | cont(2);
    }
    length = 4;
    return (static_cast<char32_t>(b0 & 0x07U) << 18) | (cont(1) << 12) | (cont(2) << 6) | cont(3);
}

}  // namespace

// --------------------------------------------------------------------------
// casts, Python arithmetic, messages

template <>
float cast_from_double<float>(double value) noexcept {
    return static_cast<float>(value);
}
template <>
double cast_from_double<double>(double value) noexcept {
    return value;
}
template <>
std::int64_t cast_from_double<std::int64_t>(double value) noexcept {
    if (!(value >= -9223372036854775808.0 && value < 9223372036854775808.0)) {
        return std::numeric_limits<std::int64_t>::min();
    }
    return static_cast<std::int64_t>(value);
}
template <>
std::int32_t cast_from_double<std::int32_t>(double value) noexcept {
    if (!(value > -2147483649.0 && value < 2147483648.0)) {
        return std::numeric_limits<std::int32_t>::min();
    }
    return static_cast<std::int32_t>(value);
}

std::int64_t python_mod(std::int64_t a, std::int64_t b) {
    if (b == 0) {
        throw PythonError("ZeroDivisionError", "integer modulo by zero");
    }
    if (b == -1) {
        return 0;
    }
    std::int64_t r = a % b;
    if (r != 0 && ((r < 0) != (b < 0))) {
        r += b;
    }
    return r;
}

double python_float_mod(double a, double b) noexcept {
    double mod = std::fmod(a, b);
    if (mod != 0.0) {
        if ((b < 0) != (mod < 0)) {
            mod += b;
        }
    } else {
        mod = std::copysign(0.0, b);
    }
    return mod;
}

// Python compares int and float exactly. The integer is not converted to
// long double (that type is a plain double on arm64 and on MSVC); the float is
// rounded to an integer on the side that keeps the comparison exact instead.
bool python_greater(std::int64_t a, double b) noexcept {
    if (std::isnan(b)) {
        return false;
    }
    if (b >= 9223372036854775808.0) {
        return false;
    }
    if (b < -9223372036854775808.0) {
        return true;
    }
    // a > b  <=>  a > floor(b) for an integer a
    return a > static_cast<std::int64_t>(std::floor(b));
}

// a < b  <=>  a < ceil(b) for an integer a
bool python_less(std::int64_t a, double b) noexcept {
    if (std::isnan(b)) {
        return false;
    }
    if (b >= 9223372036854775808.0) {
        return true;
    }
    if (b < -9223372036854775808.0) {
        return false;
    }
    return a < static_cast<std::int64_t>(std::ceil(b));
}

std::string index_out_of_bounds_message(std::int64_t index, std::int64_t size) {
    return "index " + std::to_string(index) + " is out of bounds for axis 0 with size " +
           std::to_string(size);
}

std::string boolean_index_mismatch_message(std::size_t dimension, std::size_t boolean_dimension) {
    return "boolean index did not match indexed array along dimension 0; dimension is " +
           std::to_string(dimension) + " but corresponding boolean dimension is " +
           std::to_string(boolean_dimension);
}

std::string broadcast_message(std::size_t a, std::size_t b) {
    return "operands could not be broadcast together with shapes (" + std::to_string(a) +
           ",) (" + std::to_string(b) + ",) ";
}

// --------------------------------------------------------------------------
// text and files

std::vector<std::vector<std::string>> read_bed(const std::vector<std::string>& lines) {
    std::vector<std::vector<std::string>> viewpoints;
    viewpoints.reserve(lines.size());
    for (const std::string& line : lines) {
        if (!line.empty() && line[0] == '#') {
            continue;
        }
        const std::string_view stripped = pytext::strip(line);
        std::vector<std::string_view> fields = pytext::split(stripped, "\t");
        if (fields.size() < 3) {
            // `chrom, start, end = _line[:3]` failed: split on whitespace
            fields = pytext::split_whitespace(stripped);
            if (fields.size() < 3) {
                throw PythonError("ValueError", "not enough values to unpack (expected 3, got " +
                                                    std::to_string(fields.size()) + ")");
            }
        }
        viewpoints.push_back({std::string(fields[0]), std::string(fields[1]), std::string(fields[2])});
    }
    return viewpoints;
}

TargetBed read_target_bed(const std::vector<std::string>& lines) {
    TargetBed result;
    std::unordered_map<std::string, std::size_t> index;
    std::optional<std::string> gene;
    for (const std::string& line : lines) {
        if (!line.empty() && line[0] == '#') {
            continue;
        }
        const std::string_view stripped = pytext::strip(line);
        std::vector<std::string_view> fields = pytext::split(stripped, "\t");
        std::string chrom;
        std::string start;
        std::string end;
        if (fields.size() == 3) {
            // `chrom, start, end = _line[:4]` succeeds and leaves gene as it was
            chrom = fields[0];
            start = fields[1];
            end = fields[2];
        } else {
            fields = pytext::split_whitespace(stripped);
            if (fields.size() < 4) {
                throw PythonError("ValueError", "not enough values to unpack (expected 4, got " +
                                                    std::to_string(fields.size()) + ")");
            }
            chrom = fields[0];
            start = fields[1];
            end = fields[2];
            gene = std::string(fields[3]);
        }
        if (!gene) {
            throw PythonError("UnboundLocalError",
                              "cannot access local variable 'gene' where it is not associated "
                              "with a value");
        }
        auto [it, inserted] = index.try_emplace(*gene, result.genes.size());
        if (inserted) {
            result.genes.push_back(*gene);
            result.positions.emplace_back();
        }
        TargetBed::Positions& positions = result.positions[it->second];
        positions.chromosome = std::move(chrom);
        positions.starts.push_back(std::move(start));
        positions.ends.push_back(std::move(end));
    }
    return result;
}

int truncate_for_writing(const std::string& path) noexcept {
    const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0666);
    if (fd < 0) {
        return errno;
    }
    ::close(fd);
    return 0;
}

std::optional<std::string> genomic_region(std::string_view text) {
    const std::string region = pytext::remove_whitespace(text);
    if (region.empty()) {
        return std::nullopt;
    }
    std::string cleaned;
    cleaned.reserve(region.size());
    for (const char c : region) {
        switch (c) {
            case ',': case ';': case '|': case '!': case '{': case '}': case '(': case ')':
                break;
            case '-':
                cleaned.push_back(':');
                break;
            default:
                cleaned.push_back(c);
        }
    }
    if (cleaned.empty()) {
        throw PythonError("argparse.ArgumentTypeError", std::string(text) + " is not a valid region");
    }
    return cleaned;
}

std::string change_chrom_name(std::string_view chrom) {
    if (chrom.substr(0, 3) == "chr") {
        return std::string(chrom.substr(3));
    }
    return "chr" + std::string(chrom);
}

std::string in_units(double base_position) {
    if (base_position > 1.5e6) {
        return format_fixed(base_position / 1e6, 2) + " " + " Mbp";
    }
    if (base_position > 1500) {
        return format_fixed(base_position / 1e3, 0) + " Kbp";
    }
    return format_fixed(base_position, 2) + " " + " bp";
}

bool gzip_magic(std::string_view first_bytes) noexcept {
    return first_bytes == std::string_view("\x1f\x8b", 2);
}

bool is_ascii(std::string_view text) noexcept {
    return std::all_of(text.begin(), text.end(),
                       [](char c) { return static_cast<unsigned char>(c) < 0x80U; });
}

std::string unidecode(std::string_view text, std::vector<char32_t>* surrogates) {
    std::string out;
    out.reserve(text.size());
    std::size_t pos = 0;
    while (pos < text.size()) {
        std::size_t length = 1;
        const char32_t cp = decode_at(text, pos, length);
        pos += length;
        if (cp < 0x80U) {
            out.push_back(static_cast<char>(cp));
            continue;
        }
        if (cp > 0xEFFFFU) {
            continue;
        }
        if (cp >= 0xD800U && cp <= 0xDFFFU && surrogates != nullptr) {
            surrogates->push_back(cp);
        }
        const int block = static_cast<int>(cp >> 8);
        const unsigned position = cp & 0xFFU;
        if (block >= unidecode_table::kBlocks) {
            continue;
        }
        const std::int32_t first = unidecode_table::kBlockFirst[block];
        if (first < 0 || position >= unidecode_table::kBlockCount[block]) {
            continue;
        }
        const std::uint32_t begin = unidecode_table::kOffsets[first + position];
        const std::uint32_t end = unidecode_table::kOffsets[first + position + 1];
        out.append(unidecode_table::kBlob + begin, end - begin);
    }
    return out;
}

UserRegion user_region_bounds(std::int64_t chrom_size, std::optional<std::int64_t> start,
                              std::optional<std::int64_t> end, std::string_view region_string) {
    UserRegion region;
    region.start = start.value_or(0);
    if (end && *end <= chrom_size) {
        region.end = *end;
    } else {
        region.end = chrom_size;
        region.end_is_chrom_size = true;
    }
    if (region.start > region.end || region.start < 0) {
        throw PythonError("NameError", std::string(region_string) +
                                           " not valid. The format is chrom:start:end. Without "
                                           "comas, dashes or dots. ");
    }
    return region;
}

void user_region_finish(UserRegion& region, std::optional<std::int64_t> tile_size,
                        const ChunkLimit& limit) {
    const std::int64_t tile = tile_size.value_or(0);
    if (tile != 0) {
        region.start -= python_mod(region.start, tile);
        std::int64_t moved = 0;
        if (__builtin_add_overflow(region.end, tile - python_mod(region.end, tile), &moved)) {
            throw PythonError("OverflowError", "region end out of the int64 range");
        }
        region.end = moved;
        region.end_tiled = true;
    }
    std::int64_t chunk = 0;
    if (__builtin_sub_overflow(region.end, region.start, &chunk)) {
        throw PythonError("OverflowError", "region size out of the int64 range");
    }
    if (limit.is_float) {
        if (python_greater(chunk, limit.floating)) {
            double chunk_float = limit.floating;
            if (tile != 0 && python_less(tile, chunk_float)) {
                chunk_float -= python_float_mod(chunk_float, static_cast<double>(tile));
            }
            if (std::isnan(chunk_float)) {
                throw PythonError("ValueError", "cannot convert float NaN to integer");
            }
            if (std::isinf(chunk_float)) {
                throw PythonError("OverflowError", "cannot convert float infinity to integer");
            }
            chunk = static_cast<std::int64_t>(chunk_float);
        }
    } else if (chunk > limit.integer) {
        chunk = limit.integer;
        if (tile != 0 && tile < chunk) {
            chunk -= python_mod(chunk, tile);
        }
    }
    region.chunk_size = chunk;
}

std::vector<EnlargedBin> enlarge_bins(const std::vector<std::int64_t>& starts,
                                      const std::vector<std::int64_t>& ends,
                                      const std::vector<bool>& same_chrom) {
    const std::size_t n = starts.size();
    std::vector<EnlargedBin> bins(n);
    for (std::size_t i = 0; i < n; ++i) {
        bins[i].start = starts[i];
        bins[i].end = ends[i];
        bins[i].chrom_from = i;
    }
    if (n == 0) {
        return bins;
    }
    bool chr_start = true;
    for (std::size_t idx = 0; idx + 1 < n; ++idx) {
        EnlargedBin& current = bins[idx];
        EnlargedBin& next = bins[idx + 1];
        if (chr_start) {
            current.start = 0;
            current.start_kind = Coordinate::Zero;
            current.rebuilt = true;
            chr_start = false;
        }
        const std::int64_t start_next = next.start;
        if (same_chrom[idx] && current.end != start_next) {
            const double half =
                static_cast<double>(start_next - current.end) / 2.0;  // true division
            const std::int64_t middle = start_next - static_cast<std::int64_t>(half);
            current.end = middle;
            current.end_kind = Coordinate::Middle;
            current.rebuilt = true;
            next.start = middle;
            next.start_kind = Coordinate::Middle;
            next.chrom_from = current.chrom_from;
            next.rebuilt = true;
        }
        if (!same_chrom[idx]) {
            chr_start = true;
        }
    }
    bins[n - 1].rebuilt = true;
    return bins;
}

// --------------------------------------------------------------------------
// statistics

template <class T>
double median(const std::vector<T>& values) {
    const std::size_t n = values.size();
    if (n == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if constexpr (is_float_v<T>) {
        for (const T v : values) {
            if (std::isnan(v)) {
                return std::numeric_limits<double>::quiet_NaN();
            }
        }
    }
    std::vector<T> work(values);
    const std::size_t half = n / 2;
    std::nth_element(work.begin(), work.begin() + static_cast<std::ptrdiff_t>(half), work.end());
    const T upper = work[half];
    using Acc = std::conditional_t<std::is_same_v<T, float>, float, double>;
    Acc sum = Acc(0);
    std::size_t count = 1;
    if (n % 2 == 0) {
        const T lower = *std::max_element(work.begin(), work.begin() + static_cast<std::ptrdiff_t>(half));
        sum = static_cast<Acc>(sum + static_cast<Acc>(lower));
        sum = static_cast<Acc>(sum + static_cast<Acc>(upper));
        count = 2;
    } else {
        sum = static_cast<Acc>(sum + static_cast<Acc>(upper));
    }
    const double result = static_cast<double>(sum) / static_cast<double>(count);
    if constexpr (std::is_same_v<T, float>) {
        return static_cast<double>(static_cast<float>(result));
    }
    return result;
}

template <class T>
double mean(const std::vector<T>& values) {
    const std::size_t n = values.size();
    if constexpr (std::is_same_v<T, float>) {
        const float sum = npy::pairwise_sum(values.data(), n);
        return static_cast<double>(static_cast<float>(static_cast<double>(sum) / static_cast<double>(n)));
    } else if constexpr (std::is_same_v<T, double>) {
        return npy::pairwise_sum(values.data(), n) / static_cast<double>(n);
    } else {
        // np.mean of integers or booleans accumulates in float64
        std::vector<double> converted(values.begin(), values.end());
        return npy::pairwise_sum(converted.data(), n) / static_cast<double>(n);
    }
}

template <class T>
std::optional<std::vector<std::size_t>> remove_outliers_keep(const std::vector<T>& values,
                                                             double max_deviation) {
    const std::size_t n = values.size();
    if constexpr (std::is_same_v<T, float>) {
        // A float32 array stays float32 against its float32 median.
        const float m = static_cast<float>(median(values));
        std::vector<float> deviation(n);
        for (std::size_t i = 0; i < n; ++i) {
            deviation[i] = std::fabs(values[i] - m);
        }
        const double mad = 1.4826 * median(deviation);  // np.float32 * float is float64
        if (!(mad > 0)) {
            return std::nullopt;
        }
        const float mad32 = static_cast<float>(mad);
        const float limit = static_cast<float>(max_deviation);
        std::vector<std::size_t> keep;
        for (std::size_t i = 0; i < n; ++i) {
            if (deviation[i] / mad32 <= limit) {
                keep.push_back(i);
            }
        }
        return keep;
    } else {
        const double m = median(values);
        std::vector<double> deviation(n);
        for (std::size_t i = 0; i < n; ++i) {
            deviation[i] = std::fabs(static_cast<double>(values[i]) - m);
        }
        const double mad = 1.4826 * median(deviation);
        if (!(mad > 0)) {
            return std::nullopt;
        }
        std::vector<std::size_t> keep;
        for (std::size_t i = 0; i < n; ++i) {
            if (deviation[i] / mad <= max_deviation) {
                keep.push_back(i);
            }
        }
        return keep;
    }
}

// --------------------------------------------------------------------------
// sparse matrices

std::vector<std::int64_t> Nonzero::distances() const {
    std::vector<std::int64_t> d(row.size());
    for (std::size_t i = 0; i < row.size(); ++i) {
        d[i] = std::llabs(row[i] - col[i]);
    }
    return d;
}

template <class T>
Nonzero nonzero(const SparseEntries& entries, const T* data) {
    Nonzero nz;
    for (std::size_t k = 0; k < entries.row.size(); ++k) {
        if (data[k] != T(0)) {
            nz.row.push_back(entries.row[k]);
            nz.col.push_back(entries.col[k]);
        }
    }
    return nz;
}

template <class T>
std::vector<double> expected_interactions_in_distance(std::int64_t n,
                                                      const std::vector<std::int64_t>& distance,
                                                      const T* data, std::int32_t length,
                                                      std::int32_t count) {
    std::vector<double> expected(static_cast<std::size_t>(n), 0.0);
    for (std::size_t i = 0; i < distance.size(); ++i) {
        if (distance[i] >= n) {
            raise_index(distance[i], n);
        }
        expected[static_cast<std::size_t>(distance[i])] += static_cast<double>(data[i]);
    }
    for (std::int64_t i = 0; i < n; ++i) {
        const double divisor =
            (static_cast<double>(i) * static_cast<double>(count) - static_cast<double>(length)) * -1.0;
        expected[static_cast<std::size_t>(i)] /= divisor;
    }
    return expected;
}

template <class T>
std::vector<double> expected_interactions_non_zero(std::int64_t n,
                                                   const std::vector<std::int64_t>& distance,
                                                   const T* data) {
    std::vector<double> expected(static_cast<std::size_t>(n), 0.0);
    std::vector<double> occurrences(static_cast<std::size_t>(n), 0.0);
    for (std::size_t i = 0; i < distance.size(); ++i) {
        if (distance[i] >= n) {
            raise_index(distance[i], n);
        }
        expected[static_cast<std::size_t>(distance[i])] += static_cast<double>(data[i]);
        occurrences[static_cast<std::size_t>(distance[i])] += 1.0;
    }
    for (std::size_t i = 0; i < expected.size(); ++i) {
        expected[i] = zero_non_finite(expected[i] / occurrences[i]);
    }
    return expected;
}

namespace {

// The values of data grouped by distance, in storage order, for the
// distances in [low, high].
template <class T>
std::vector<std::vector<T>> group_by_distance(const std::vector<std::int64_t>& distance,
                                              const T* data, std::int64_t low, std::int64_t high) {
    std::vector<std::vector<T>> groups(static_cast<std::size_t>(std::max<std::int64_t>(high - low + 1, 0)));
    for (std::size_t k = 0; k < distance.size(); ++k) {
        if (distance[k] >= low && distance[k] <= high) {
            groups[static_cast<std::size_t>(distance[k] - low)].push_back(data[k]);
        }
    }
    return groups;
}

}  // namespace

template <class T>
std::vector<double> distance_sums(std::int64_t size, const std::vector<std::int64_t>& distance,
                                  const T* data, std::size_t data_size, std::int64_t low,
                                  std::int64_t high) {
    std::vector<double> sums(static_cast<std::size_t>(size), 0.0);
    if (low >= high) {
        return sums;
    }
    if (distance.size() != data_size) {
        raise_mismatch(data_size, distance.size());
    }
    const std::vector<std::vector<T>> groups = group_by_distance(distance, data, low, high - 1);
    for (std::int64_t i = low; i < high; ++i) {
        const double sum = np_sum(groups[static_cast<std::size_t>(i - low)]);
        const std::int64_t index = i < 0 ? i + size : i;
        if (index < 0 || index >= size) {
            raise_index(i, size);
        }
        sums[static_cast<std::size_t>(index)] = sum;
    }
    return sums;
}

template <class T>
std::optional<std::vector<double>> expected_interactions(std::int64_t n,
                                                         const std::vector<std::int64_t>& distance,
                                                         const T* data, std::size_t data_size) {
    if (distance.empty()) {
        return std::nullopt;
    }
    const auto [lo, hi] = std::minmax_element(distance.begin(), distance.end());
    std::vector<double> expected = distance_sums(n, distance, data, data_size, *lo, *hi + 1);
    for (std::int64_t i = 0; i < n; ++i) {
        const double occurrences = static_cast<double>(n + 1 - i);
        expected[static_cast<std::size_t>(i)] =
            zero_non_finite(expected[static_cast<std::size_t>(i)] / occurrences);
    }
    return expected;
}

template <class T>
void compute_zscore(std::int64_t n, const std::vector<std::int64_t>& distance, T* data,
                    std::size_t data_size, std::optional<double> depth, bool depth_is_float,
                    const std::optional<std::vector<double>>& expected) {
    std::int64_t max_distance = n;
    if (depth) {
        if (*depth < static_cast<double>(n)) {
            if (depth_is_float) {
                // max_distance = min(n, pDepth) is the float, which range() refuses;
                // the mask below is applied first
                max_distance = -1;
            } else {
                max_distance = static_cast<std::int64_t>(*depth);
            }
        }
        if (distance.size() != data_size) {
            raise_mismatch(data_size, distance.size());
        }
        for (std::size_t k = 0; k < data_size; ++k) {
            if (static_cast<double>(distance[k]) > *depth) {
                data[k] = T(0);
            }
        }
        if (max_distance == -1) {
            throw PythonError("TypeError", "'float' object cannot be interpreted as an integer");
        }
    }
    if (max_distance <= 0) {
        return;
    }
    if (distance.size() != data_size) {
        raise_mismatch(data_size, distance.size());
    }
    if (!expected) {
        throw PythonError("TypeError", "'NoneType' object is not subscriptable");
    }
    std::vector<std::vector<std::size_t>> positions(static_cast<std::size_t>(max_distance));
    for (std::size_t k = 0; k < distance.size(); ++k) {
        if (distance[k] >= 0 && distance[k] < max_distance) {
            positions[static_cast<std::size_t>(distance[k])].push_back(k);
        }
    }
    for (std::int64_t d = 0; d < max_distance; ++d) {
        const std::vector<std::size_t>& mask = positions[static_cast<std::size_t>(d)];
        const double e = (*expected)[static_cast<std::size_t>(d)];
        const double x = 0.0 - e;  // x_minus_mu starts at zero
        std::vector<double> squares(mask.size(), x * x);
        double sum_for_sigma = npy::pairwise_sum(squares.data(), squares.size());
        const std::int64_t occurrences = n + 1 - d;
        const double missing =
            (0.0 - e) * static_cast<double>(occurrences - static_cast<std::int64_t>(mask.size()));
        sum_for_sigma += missing * missing;
        const double sigma = std::sqrt(sum_for_sigma / static_cast<double>(occurrences));
        for (const std::size_t k : mask) {
            data[k] = cast_from_double<T>(x / sigma);
        }
    }
}

template <class T>
std::vector<T> obs_exp_divide(const T* data, std::size_t data_size,
                              const std::vector<double>& expected,
                              const std::vector<std::int64_t>& distance, double non_finite) {
    const std::size_t L = distance.size();
    std::vector<double> gathered(L);
    for (std::size_t i = 0; i < L; ++i) {
        // np.ceil(np.absolute(row - col) / 2).astype(np.int32)
        const std::int64_t index = (distance[i] + 1) / 2;
        if (index >= static_cast<std::int64_t>(expected.size())) {
            raise_index(index, static_cast<std::int64_t>(expected.size()));
        }
        gathered[i] = expected[static_cast<std::size_t>(index)];
    }
    std::size_t size = 0;
    if (data_size == L || L == 1) {
        size = data_size;
    } else if (data_size == 1) {
        size = L;
    } else {
        throw PythonError("ValueError", broadcast_message(data_size, L));
    }
    std::vector<T> out(size);
    for (std::size_t i = 0; i < size; ++i) {
        const double observed = static_cast<double>(static_cast<float>(data[data_size == 1 ? 0 : i]));
        double q = observed / gathered[L == 1 ? 0 : i];
        if (!std::isfinite(q)) {
            q = non_finite;
        }
        out[i] = cast_from_double<T>(q);
    }
    return out;
}

template <class T>
void zero_from_distance(const std::vector<std::int64_t>& distance, T* data, std::size_t data_size,
                        double limit) {
    if (distance.size() != data_size) {
        raise_mismatch(data_size, distance.size());
    }
    for (std::size_t k = 0; k < data_size; ++k) {
        if (static_cast<double>(distance[k]) >= limit) {
            data[k] = T(0);
        }
    }
}

namespace {

// scipy's sum dtype (get_sum_dtype) for T.
template <class T>
using RowSumType = std::conditional_t<is_float_v<T>, T, std::int64_t>;

template <class R>
R add(R a, R b) noexcept {
    if constexpr (std::is_integral_v<R>) {
        return wrap_add(a, b);
    } else {
        return a + b;
    }
}

}  // namespace

template <class T>
std::vector<float> obs_exp_non_zero_values(const SparseEntries& entries,
                                           const std::vector<std::int64_t>& indptr, bool csr,
                                           const T* data, std::size_t data_size,
                                           const Nonzero& nonzero,
                                           const std::vector<double>& expected,
                                           bool ligation_factor, bool to_epsilon) {
    using R = RowSumType<T>;
    const std::size_t rows = static_cast<std::size_t>(entries.rows);
    std::vector<R> row_sums(rows, R(0));
    R total = R(0);
    if (ligation_factor) {
        if (csr) {
            // sum(axis=1): np.add.reduceat over every non-empty row, in T
            for (std::size_t r = 0; r < rows; ++r) {
                const std::size_t begin = static_cast<std::size_t>(indptr[r]);
                const std::size_t end = static_cast<std::size_t>(indptr[r + 1]);
                if (end <= begin) {
                    continue;
                }
                if constexpr (is_float_v<T>) {
                    // a segment of one value is copied, not added to zero
                    row_sums[r] = end - begin == 1
                                      ? data[begin]
                                      : data[begin] + npy::pairwise_sum_unbuffered(
                                                          data + begin + 1, end - begin - 1);
                } else {
                    T acc = data[begin];
                    for (std::size_t k = begin + 1; k < end; ++k) {
                        acc = wrap_add(acc, data[k]);
                    }
                    row_sums[r] = static_cast<R>(acc);
                }
            }
            // sum(): csr_matvec against ones in R, each row from zero in order
            std::vector<R> matvec(rows, R(0));
            for (std::size_t r = 0; r < rows; ++r) {
                R acc = R(0);
                for (std::size_t k = static_cast<std::size_t>(indptr[r]);
                     k < static_cast<std::size_t>(indptr[r + 1]); ++k) {
                    acc = add<R>(acc, static_cast<R>(data[k]));
                }
                matvec[r] = acc;
            }
            if constexpr (is_float_v<R>) {
                total = npy::pairwise_sum(matvec.data(), matvec.size());
            } else {
                for (const R v : matvec) {
                    total = wrap_add(total, v);
                }
            }
        } else {
            // csc: both sums are csc_matvec, rows accumulated in storage order
            for (std::size_t k = 0; k < entries.row.size(); ++k) {
                const std::size_t r = static_cast<std::size_t>(entries.row[k]);
                row_sums[r] = add<R>(row_sums[r], static_cast<R>(data[k]));
            }
            if constexpr (is_float_v<R>) {
                total = npy::pairwise_sum(row_sums.data(), row_sums.size());
            } else {
                for (const R v : row_sums) {
                    total = wrap_add(total, v);
                }
            }
        }
    }

    std::vector<float> values(data_size);
    for (std::size_t k = 0; k < data_size; ++k) {
        values[k] = static_cast<float>(data[k]);
    }
    for (std::size_t i = 0; i < nonzero.size(); ++i) {
        const std::int64_t r = nonzero.row[i];
        const std::int64_t c = nonzero.col[i];
        double divisor = expected[static_cast<std::size_t>(std::llabs(r - c))];
        if (ligation_factor) {
            if (c >= static_cast<std::int64_t>(rows)) {
                raise_index(c, static_cast<std::int64_t>(rows));
            }
            const R a = row_sums[static_cast<std::size_t>(r)];
            const R b = row_sums[static_cast<std::size_t>(c)];
            if constexpr (std::is_integral_v<R>) {
                divisor *= static_cast<double>(wrap_mul(a, b)) / static_cast<double>(total);
            } else {
                const R product = a * b;
                divisor *= static_cast<double>(static_cast<R>(product / total));
            }
        }
        values[i] = static_cast<float>(static_cast<double>(values[i]) / divisor);
    }
    const float epsilon = to_epsilon ? static_cast<float>(0.000000001) : 0.0F;
    for (float& v : values) {
        if (!std::isfinite(v)) {
            v = epsilon;
        }
    }
    return values;
}

template <class T>
std::vector<std::size_t> eliminate_zeros(std::vector<std::int64_t>& indptr, const T* data) {
    std::vector<std::size_t> keep;
    std::int64_t nnz = 0;
    std::int64_t row_end = 0;
    for (std::size_t i = 0; i + 1 < indptr.size(); ++i) {
        std::int64_t jj = row_end;
        row_end = indptr[i + 1];
        for (; jj < row_end; ++jj) {
            if (data[jj] != T(0)) {
                keep.push_back(static_cast<std::size_t>(jj));
                ++nnz;
            }
        }
        indptr[i + 1] = nnz;
    }
    return keep;
}

// --------------------------------------------------------------------------
// instantiations for the dtypes the Python API accepts

#define HICX_UTILITIES_INSTANTIATE(T)                                                          \
    template double median<T>(const std::vector<T>&);                                          \
    template double mean<T>(const std::vector<T>&);                                            \
    template std::optional<std::vector<std::size_t>> remove_outliers_keep<T>(                  \
        const std::vector<T>&, double);                                                        \
    template Nonzero nonzero<T>(const SparseEntries&, const T*);                               \
    template std::vector<double> expected_interactions_in_distance<T>(                         \
        std::int64_t, const std::vector<std::int64_t>&, const T*, std::int32_t, std::int32_t); \
    template std::vector<double> expected_interactions_non_zero<T>(                            \
        std::int64_t, const std::vector<std::int64_t>&, const T*);                             \
    template std::vector<double> distance_sums<T>(std::int64_t, const std::vector<std::int64_t>&, \
                                                  const T*, std::size_t, std::int64_t,         \
                                                  std::int64_t);                               \
    template std::optional<std::vector<double>> expected_interactions<T>(                      \
        std::int64_t, const std::vector<std::int64_t>&, const T*, std::size_t);                \
    template void compute_zscore<T>(std::int64_t, const std::vector<std::int64_t>&, T*,        \
                                    std::size_t, std::optional<double>, bool,                  \
                                    const std::optional<std::vector<double>>&);                \
    template std::vector<T> obs_exp_divide<T>(const T*, std::size_t, const std::vector<double>&, \
                                              const std::vector<std::int64_t>&, double);       \
    template void zero_from_distance<T>(const std::vector<std::int64_t>&, T*, std::size_t,      \
                                        double);                                               \
    template std::vector<float> obs_exp_non_zero_values<T>(                                    \
        const SparseEntries&, const std::vector<std::int64_t>&, bool, const T*, std::size_t,   \
        const Nonzero&, const std::vector<double>&, bool, bool);                               \
    template std::vector<std::size_t> eliminate_zeros<T>(std::vector<std::int64_t>&, const T*);

HICX_UTILITIES_INSTANTIATE(std::int32_t)
HICX_UTILITIES_INSTANTIATE(std::int64_t)
HICX_UTILITIES_INSTANTIATE(float)
HICX_UTILITIES_INSTANTIATE(double)

#undef HICX_UTILITIES_INSTANTIATE

}  // namespace hicx::utilities
