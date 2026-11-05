#include "hicx/validate_locations.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/json_lite.hpp"

namespace hicx::validate_locations {

// BedTool.from_dataframe(df).sort().to_dataframe(...)
TextTable sorted_table(const TextTable& table) {
    std::vector<std::string> chrom(table.rows());
    std::vector<std::int64_t> start(table.rows());
    for (std::size_t r = 0; r < table.rows(); ++r) {
        chrom[r] = table.column(0).text(r, ".");
        start[r] = table.column(1).as_int(r);
    }
    return table.select_rows(bedtools::sort_order(chrom, start)).reparse();
}

bedtools::Intervals intervals_of(const TextTable& table, std::size_t chrom_col,
                                       std::size_t start_col, std::size_t end_col) {
    bedtools::Intervals out;
    out.chrom.reserve(table.rows());
    out.start.reserve(table.rows());
    out.end.reserve(table.rows());
    for (std::size_t r = 0; r < table.rows(); ++r) {
        out.push_back(table.column(chrom_col).text(r, "."),
                      table.column(start_col).as_int(r),
                      table.column(end_col).as_int(r));
    }
    return out;
}

TextTable table_of(const bedtools::Intervals& intervals) {
    TableColumn chrom;
    chrom.type = ColumnType::String;
    chrom.strings = intervals.chrom;
    TableColumn start;
    start.type = ColumnType::Int64;
    start.ints = intervals.start;
    TableColumn end;
    end.type = ColumnType::Int64;
    end.ints = intervals.end;
    return TextTable::with_columns({chrom, start, end});
}

BinningSteps binning_steps(const TextTable& table, std::int64_t bin_size, bool merge) {
    BinningSteps steps;
    TextTable binned = table;
    TableColumn& start = binned.column(1);
    TableColumn& end = binned.column(2);
    std::vector<std::int64_t> new_start(binned.rows());
    std::vector<std::int64_t> new_end(binned.rows());
    const double divisor = static_cast<double>(bin_size);
    for (std::size_t r = 0; r < binned.rows(); ++r) {
        // (df[1] / binsize).astype(int) is a float64 division truncated
        // towards zero, not an integer division.
        const double scaled_start = static_cast<double>(table.column(1).as_int(r)) / divisor;
        const double scaled_end = static_cast<double>(table.column(2).as_int(r)) / divisor;
        new_start[r] = static_cast<std::int64_t>(scaled_start) * bin_size;
        new_end[r] = (static_cast<std::int64_t>(scaled_end) + 1) * bin_size;
    }
    start.type = ColumnType::Int64;
    start.ints = std::move(new_start);
    start.floats.clear();
    start.strings.clear();
    end.type = ColumnType::Int64;
    end.ints = std::move(new_end);
    end.floats.clear();
    end.strings.clear();
    steps.binned = binned;
    steps.kept = binned.first_occurrences();
    steps.deduped = binned.select_rows(steps.kept);
    if (merge) {
        const bedtools::Intervals merged = bedtools::merge(intervals_of(steps.deduped, 0, 1, 2));
        steps.merged = table_of(merged).reparse();
        steps.sorted = sorted_table(steps.merged);
    }
    return steps;
}

TextTable apply_binning(const TextTable& table, std::int64_t bin_size, bool merge) {
    BinningSteps steps = binning_steps(table, bin_size, merge);
    return merge ? std::move(steps.sorted) : std::move(steps.deduped);
}

void apply_chr_prefix(TextTable& table, std::size_t column,
                      const std::string& mode) {
    if (mode == "add") {
        table.add_chr_prefix(column);
    } else if (mode == "remove") {
        table.remove_chr_prefix(column);
    }
}

TextTable read_three_column(const std::string& path, const std::string& prefix) {
    TextTable table = TextTable::read_tsv(path).select_columns({0, 1, 2});
    apply_chr_prefix(table, 0, prefix);
    return table;
}

CoolPeaks correlate_cool(const std::string& path, const TextTable& loops) {
    CoolFile cool(path);
    CoolPeaks result;
    result.number_of_peaks = cool.nnz();

    const json::Value* bin_size_value = cool.info_value("bin-size");
    std::int64_t bin_size = 0;
    if (bin_size_value != nullptr && bin_size_value->is_number()) {
        bin_size = static_cast<std::int64_t>(bin_size_value->as_double());
    }
    if (bin_size <= 0) {
        throw std::runtime_error(
            "the validation cooler has no fixed bin size; cooler's variable bin "
            "path of region_to_extent is not ported");
    }

    std::map<std::string, std::int64_t> chrom_length;
    for (std::size_t i = 0; i < cool.chrom_names().size(); ++i) {
        chrom_length[cool.chrom_names()[i]] = cool.chrom_lengths()[i];
    }
    std::map<std::string, std::int64_t> chrom_offset;
    {
        const BinTable bins(cool.read_bins());
        for (const auto& entry : bins.chrom_bin_boundaries()) {
            chrom_offset[entry.first] = entry.second.first;
        }
    }

    const CsrMatrix matrix = cool.read_matrix();
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();

    // Any stored pixel with row in [r0, r1) and column in [c0, c1)?
    const auto block_has_pixel = [&](std::int64_t r0, std::int64_t r1, std::int64_t c0,
                                     std::int64_t c1) {
        if (r0 >= r1 || c0 >= c1) {
            return false;
        }
        const std::int64_t rows = matrix.rows();
        r0 = std::max<std::int64_t>(r0, 0);
        r1 = std::min<std::int64_t>(r1, rows);
        for (std::int64_t row = r0; row < r1; ++row) {
            const auto begin = indices.begin() + indptr[static_cast<std::size_t>(row)];
            const auto end = indices.begin() + indptr[static_cast<std::size_t>(row) + 1];
            const auto lower =
                std::lower_bound(begin, end, static_cast<std::int32_t>(c0));
            if (lower != end && *lower < static_cast<std::int32_t>(c1)) {
                return true;
            }
        }
        return false;
    };

    const auto extent = [&](const std::string& chrom, std::int64_t start,
                            std::int64_t end,
                            std::pair<std::int64_t, std::int64_t>* out) {
        const auto length = chrom_length.find(chrom);
        if (length == chrom_length.end()) {
            return false;  // KeyError inside cooler.parse_region
        }
        if (end < start || end > length->second || start < 0) {
            return false;  // ValueError inside cooler.parse_region
        }
        const std::int64_t offset = chrom_offset.at(chrom);
        const std::int64_t first = offset + start / bin_size;
        const std::int64_t last =
            offset + (end + bin_size - 1) / bin_size;  // ceil
        *out = {first, last};
        return true;
    };

    result.matched.assign(loops.rows(), false);
    for (std::size_t r = 0; r < loops.rows(); ++r) {
        std::pair<std::int64_t, std::int64_t> x{};
        std::pair<std::int64_t, std::int64_t> y{};
        if (!extent(loops.column(0).text(r, "."), loops.column(1).as_int(r),
                    loops.column(2).as_int(r), &x)) {
            continue;
        }
        if (!extent(loops.column(3).text(r, "."), loops.column(4).as_int(r),
                    loops.column(5).as_int(r), &y)) {
            continue;
        }
        // The cooler is stored symmetric-upper, so the block the Python gets
        // back is the union of the upper-triangle hits in either orientation.
        result.matched[r] =
            block_has_pixel(x.first, x.second, y.first, y.second) ||
            block_has_pixel(y.first, y.second, x.first, x.second);
    }
    return result;
}

TextTable intersect_count_table(const TextTable& a, const TextTable& b) {
    const std::vector<std::int64_t> counts = bedtools::intersect_count(intervals_of(a, 0, 1, 2), intervals_of(b, 0, 1, 2));
    std::vector<TableColumn> columns = a.columns();
    TableColumn count;
    count.type = ColumnType::Int64;
    count.ints = counts;
    columns.push_back(std::move(count));
    return TextTable::with_columns(std::move(columns)).reparse();
}

}  // namespace hicx::validate_locations
